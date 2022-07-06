import math
import numpy as np
import ete3

# tangleType = "V"
tangleType = "E"
if tangleType == "V":
    import VertexTangleSet as tset
elif tangleType == "E":
    import EdgeTangleSet as tset

import pandas as pd
import collections as coll
import colorsys
from matplotlib import cm
import matplotlib.pyplot as plt
import igraph as ig
import itertools as iter
import protChecker
import sklearn.metrics as skl
import os
import socket
import random

from py2cytoscape.data.cyrest_client import CyRestClient

# from igraph import arpack_options
# arpack_options.maxiter=300000

class graphCD():
    def __init__(self, job, log):

        # todo Chuck some error checking in here later ******
        self.job = job
        self.log = log

        self.log.log("{}".format(job['outName']))

        self.doPrint = False

        if "construct" in job and job["construct"]:
            graph = self.constructRandom(job)
            # print("Edges: ", graph.ecount())
            # ig.plot(graph)
            # exit()
        else:
            graph = ig.Graph.Read_Ncol(job['inFile'], names=True, directed=False)

        if not graph.is_weighted():
            graph.es["weight"] = 1

        graph.simplify(combine_edges="sum")
        self.giantComp = graph.clusters().giant()

        # todo remove after testing
        # self.igPrint()
        # exit()

        self.protChecker = None
        self.colmaps = None
        # initialising to None as easier to check for than not existing. So later only creates if needed.

        print("Number of nodes: {}".format(self.giantComp.vcount()))
        print("Number of edges: {}".format(self.giantComp.ecount()))

    def constructRandom(self, job):
        if job["M"] < job["N"]:
            exit("Error: asking for random graph with fewer edges than nodes")

        e = 0
        # create random graph until get one with at least enough edges,
        # then remove randomly (as long as they don't disconnect the graph)
        # until get exactly the right number of edges
        # need to include the fudge because otherwise it might keep running and running
        # and never get enough edges.
        fudge = 0
        while e < job["M"]:
            # the +1 is there so don't get p = 1

            # these numbers are based on linear regression done previously
            p = (job["M"] - 3.52663056*job["N"] +201.35998895)/(206.97060932) + fudge
            if p >= 1:
                p = 0.99999999
            if p <= 0:
                p = 0.00000001

            graph = ig.Graph.Forest_Fire(job["N"], p)
            e = graph.ecount()
            fudge+=((job["M"]-e)/(job["M"]*100)*(1-p))
            # fudge += ((e - e_actual)/(e*100) * (1 - p))

        while e > job["M"]:
            edgeDel = random.choice([e for e in graph.es.indices if e not in graph.bridges()])
            graph.delete_edges(edgeDel)
            e = graph.ecount()

        graph.vs["name"] = [str(v) for v in graph.vs.indices]
        return graph

    ## Creates a list of colours evenly distributed around the hue spectrum
    # with fixed saturation and value
    # def getColourList(self, format="rgb"):
    #     self.colList = []
    #
    #     if format=="hsv":
    #         for i in range(0,commRange):
    #             self.colList.append(
    #                 "{:0<4f} {:0<4f} {:0<4f}".format(i/commRange, 1, 0.7))
    #     elif format=="rgb":
    #         for i in range(0,commRange):
    #             rgblist = list(colorsys.hsv_to_rgb(i/commRange, 1, 0.7))
    #             self.colList.append(
    #                 "#{}".format("".join(map(lambda x: hex(int(x*255)).\
    #                     split("0x")[1].zfill(2),rgblist))))
    #     else:
    #         print("Incorrect colour format, press any key to continue")
    #         input()


    def getColour(self, nodedep, stream = 0):

        if self.colmaps is None:
            self.colStreams = 1  # todo fix later
            self.colDep = self.TangleSet.kmax - self.TangleSet.kmin + 2
            # note, +1 for tangle depth, +1 again for leaving white for not in tangle

            self.colmaps = [cm.get_cmap("Purples", self.colDep+1)]


        colCode = self.colmaps[stream](nodedep/self.colDep)

        # convert to #ffffff format
        return("#{}".format("".join(map(lambda x: hex(int(x*255)).\
                split("0x")[1].zfill(2),colCode[0:3]))))




    def findTangleComms(self, dep = 4, sepsOnly = False):

        #### dep is added to kmin to give *separation* order
        #### add 1 to get *tangle* order
        #### note, finds dep+1 total levels

        if tangleType == "V":
            self.groundset = {"{}_{}".format(self.giantComp.vs[edge.source]['name'], self.giantComp.vs[edge.target]['name']) for edge in self.giantComp.es}
            self.TangleSet = tset.VertexTangleSet(self.giantComp, self.job, self.log)
        elif tangleType == "E":
            self.groundset = set(self.giantComp.vs["name"])
            self.TangleSet = tset.EdgeTangleSet(self.giantComp, self.job, self.log)
            print("Finding Edge connectivity tangles")
        else:
            print("incorrect tangle type {} specified. Use V or E".format(tangleType))

        timings = self.TangleSet.findAllTangles(depth=dep, sepsOnly=sepsOnly)

        self.assignCommunities(thres = 0.95)

        # if "Yeast" in self.job["outName"]:
        #     quality = self.evaluateCommunities()
            # todo something with quality

        self.doPrint = False
        # if self.doPrint and socket.gethostname().find("ginger") > -1:
        #     self.cytoPrint()
        if self.doPrint:
            self.igPrint()

        # todo return timings for each k
        return(self.giantComp.vcount(), self.giantComp.ecount(), self.TangleSet.getTangleCounts(), timings)

    # todo this probably doesn't work for vertex connectivity tangles!
    # also note currently only does *once* at end for all k levels -
    # could change to do each level?
    def assignCommunities(self, thres):
        self.foundcover = pd.DataFrame(index = sorted(self.giantComp.vs["name"]), columns=range(self.TangleSet.currentTangle), dtype=int)

        tangOrders = []

        # Note that these orders are the order of the *separations*, not the tangles.
        tangNum = 0
        for order in range(self.TangleSet.kmin, self.TangleSet.kmax+1):
            for tang in self.TangleSet.TangleLists[order]:
                numSeps = len(tang.smallSides)
                tangOrders.append(order)

                # tang.smallSides is list of sets
                smallCounter = coll.Counter([x for s in tang.smallSides for x in s])

                for v in range(self.TangleSet.groundsetSize):
                    self.foundcover.loc[self.TangleSet.names[v], tangNum] = \
                        1 if (1-smallCounter[v]/numSeps >= thres) else 0

                tangNum+=1

        outfile = "{}/{}-TangNodes.csv". \
            format(self.job['outputFolder'], self.job['outName'])
        # making copy to add row for orders without messing with anything else.
        cover_copy = self.foundcover.copy(deep = True)
        cover_copy = cover_copy.append(pd.Series(tangOrders, index=self.foundcover.columns, name="order"), ignore_index=False)
        cover_copy.to_csv(outfile)


        # this makes sure only those communities with at least 3 modes are included.
        # the astype is necessary as df init as NaNs, which are stored as floats.
        self.foundcover = self.foundcover.loc[:, (self.foundcover.sum(axis=0) >= 3)].astype(np.int8)

        for node in self.giantComp.vs:
            node["color"] = "#ffffff"


        # This bit is working with the tangle tree to assign colours to the comms.
        # note that traverse default is BFS, which is what we want - deeper comms will
        # overwrite colours for parent comms.
        # todo if change to every level, need to do on *copy* of tangletree.
        for treenode in self.TangleSet.TangleTree.traverse():
            if "T" not in treenode.name:
                # keeps only complete tangles
                # print("deleting non-tangle treenodes")
                treenode.delete(prevent_nondicotomic=False)
            elif int(treenode.name.replace("T", "")) not in self.foundcover.columns:
                # keeps only tangles with >= 3 nodes
                # print("deleting trivial tangle treenodes")
                treenode.delete(prevent_nondicotomic=False)
            else:
                # these are the actual communities we want to colour
                commIndex = int(treenode.name.replace("T", ""))
                treedep = treenode.get_distance(self.TangleSet.TangleTree)
                for nodeName in self.foundcover.index[self.foundcover[commIndex]==1].tolist():
                    self.giantComp.vs.find(nodeName)["color"] = self.getColour(treedep)
        self.printCommTree()

    def printCommTree(self):
        outfile = "{}/{}-CommTree.png". \
            format(self.job['outputFolder'], self.job['outName'])

        style = ete3.NodeStyle()
        style["size"] = 0
        for n in self.TangleSet.TangleTree.traverse():
            n.set_style(style)
            # n.dist *= 4

        ts = ete3.TreeStyle()
        # ts.show_branch_length = False
        ts.show_scale = False
        ts.show_leaf_name = False

        def my_layout(node):
            F = ete3.TextFace(node.name, tight_text=True)
            # F.fsize = 16
            F.rotation = -90
            ete3.add_face_to_node(F, node, column=0, position="branch-right")

        ts.layout_fn = my_layout
        ts.rotation = (90)

        ts.branch_vertical_margin = 10
        ts.show_branch_length = False

        for node in self.TangleSet.TangleTree.iter_descendants():
            node.dist *= 4

        ts.branch_vertical_margin = 8
        ts.scale = 360

        try:
            self.TangleSet.TangleTree.render(outfile, tree_style=ts)
        except Exception as rendError:
            print(rendError)


    # def analyseOverlapComms(self):

        # quality = coll.defaultdict(float)

        # if "Yeast" in self.job["outName"]:
        # # if False:
        #     quality = self.evaluateCommunities(self.foundcover)
        # else:
        #     quality["enrichment"] = 0
        #     quality["mi"] = 0
        #     quality["NMI"] = 0

        # quality["commcover"] = len(graphData[graphData["commCount"] > 0]) / self.giantComp.vcount()
        # quality["overlapcover"] = graphData["commCount"].mean()

        # qualFile = "{}/{}quality.csv".format(self.job['outputFolder'], self.job['outName'])
        # if not os.path.exists(qualFile):
        #     text = "numComms\tcommQ\toverlapQ\tNMI\tcommCover\toverlapCover\n"
        #     with open(qualFile, 'w') as the_file:
        #         the_file.write(text)

        # text = "{:d}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.2f}\t{:.2f}\n".format(
        #     self.numOverlapComms,
        #     quality["enrichment"],
        #     quality["mi"],
        #     quality["NMI"],
        #     quality["commcover"],
        #     quality["overlapcover"]
        # )

    def igPrint(self):
        visual_style = {
            "vertex_label": self.giantComp.vs["name"],
            "vertex_shape": "rectangle",
            "vertex_width": max(map(len, self.giantComp.vs["name"])) * 10 + 10,
            "vertex_height": 20
        }

        ig.plot(self.giantComp, **visual_style)

        dummy = 1

    def cytoPrint(self, graphToPrint=None):
        if graphToPrint is None:
            graphToPrint = self.giantComp

        cy = CyRestClient()
        cy.session.delete() ### WHy do I need to do this?
        cyGraph = cy.network.create_from_igraph(graphToPrint)

        # layouts = cy.layout.get_all()
        # print("--------------")
        # print(json.dumps(layouts, indent=4))
        #
        # print(cy.layout.apply.__code__.co_varnames)
        # print(dir(cy.layout))
        #
        # print("--------------")

        # cy.layout.apply(name='kamada-kawai', network=cyGraph)
        cy.layout.apply(name='force-directed', network=cyGraph)
        cyStyle = cy.style.create('default')
        # cyStyle._Style__url = "http://red-ginger.sms.vuw.ac.nz/v1/styles/CommColours/"

        # print(cy.style.vps.get_all())
        # vps = pd.Series(cy.style.vps.get_all())
        # vps.to_csv("VisualProperties.csv")

        basic_settings = {
            # You can set default values as key-value pairs.
            # 'NODE_SIZE': 80,
            # 'NODE_SHAPE': 'Ellipse',
            # 'NODE_SHAPE': 'Rectangle',
            # 'NODE_HEIGHT': 20,
            # 'NODE_WIDTH': 80,
            # 'NODE_BORDER_WIDTH': 1,
            # 'NODE_PAINT': '#FFFFFF',
            # 'NODE_BORDER_PAINT': '#000000',
            'NODE_LABEL_COLOR': '#000000',
            'EDGE_WIDTH': 1,
            'EDGE_TRANSPARENCY': 255,
            # 'NODE_FILL_COLOR': '#DCDCDC',
        }

        cyStyle.update_defaults(basic_settings)

        # cyStyle.create_passthrough_mapping(column='format',
        #     vp='NODE_CUSTOMGRAPHICS_1', col_type='String')
        # cyStyle.create_passthrough_mapping(column='edgecol',
        #     vp='EDGE_STROKE_UNSELECTED_PAINT', col_type='String')
        # cyStyle.create_passthrough_mapping(column='id',
        #     vp='NODE_LABEL', col_type='String')
        cyStyle.create_passthrough_mapping(column='name',
            vp='NODE_LABEL', col_type='String')
        cyStyle.create_passthrough_mapping(column='color',
            vp='NODE_FILL_COLOR', col_type='String')
        # cyStyle.create_passthrough_mapping(column='flow',
        #     vp='EDGE_LABEL', col_type='String')
        cy.style.apply(cyStyle, cyGraph)

        cyPDF = cyGraph.get_pdf()

        PDFfilename = '{}/{}-CommLayout.pdf'.\
            format(self.job['outputFolder'], self.job['outName'])
        PDFoutfile = open(PDFfilename, 'wb')
        PDFoutfile.write(cyPDF)
        PDFoutfile.close()

        # os.system("pdfcrop {} {}".format(PDFfilename, PDFfilename))

    def evaluateCommunities(self):
        # see analyseOverlapComms

        quality = coll.defaultdict(float)

        if self.protChecker is None:
            self.protChecker = protChecker.protChecker(self.giantComp.vs['name'])

        NMI = self.protChecker.calculateNMI(self.foundcover)
        print("Found NMI: ", NMI)
        commQual = self.protChecker.getSimilarityRatio(self.foundcover)
        print("Found commQual: ", commQual)
        # overlapMI = self.protChecker.getOverlapMI(self.foundcover)
        # print("Found numComms MI: ", overlapMI)
        # note mi is Mutual inf for (per node) num of comms vs num of GO terms
        # NMI is Normalised mutual inf between assigned comms and comms by GO terms

        # todo - do somthing with the qual measures
        # todo also re-add the coverage ratio
        return(quality)


    def overLapCliquePercolation(self):
        # https://stackoverflow.com/questions/20063927/overlapping-community-detection-with-igraph-or-other-libaries

        minClique = 3
        maxClique = 6

        cliques = list(map(set, self.giantComp.maximal_cliques(min=minClique)))
        # so, each clique is assigned an index in the list.

        for k in range(minClique,maxClique+1):
            filteredCliques = [c for c in cliques if len(c) >= k]

            edgelist = []
            # iter.combinations(iterable, 2) returns every possible pair of values
            # in the iterable (irrespective of order), but returns them ordered
            for i, j in iter.combinations(range(len(filteredCliques)), 2):
                if len(filteredCliques[i].intersection(filteredCliques[j])) >= k-1:
                    edgelist.append((i, j))

            cliqueLinks = ig.Graph(edgelist, directed=False)
            cliqueComps = cliqueLinks.components()

            numComps = len(cliqueComps)
            if numComps == 0:
                break

            self.foundcover = pd.DataFrame(index = sorted(self.giantComp.vs["name"]), columns=range(numComps), dtype=int)
            for col in self.foundcover.columns:
                self.foundcover[col].values[:] = 0

            commIndex = 0
            for comp in cliqueComps:
                # each component is a list of clique indices

                for cliqueIndex in comp:
                    for node in filteredCliques[cliqueIndex]:
                        # todo fix this so not refer to tangleset
                        self.foundcover.loc[self.TangleSet.names[node], commIndex] = 1
                commIndex+=1

            # this makes sure only those communities with at least 3 modes are included.
            # the astype is necessary as df init as NaNs, which are stored as floats.
            self.foundcover = self.foundcover.loc[:, (self.foundcover.sum(axis=0) >= 3)].astype(np.int8)

            print("CPM: {}".format(k))
            self.evaluateCommunities()

