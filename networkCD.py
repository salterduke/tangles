import math
import numpy as np

# tangleType = "V"
tangleType = "E"
if tangleType == "V":
    import VertexTangleSet as tset
elif tangleType == "E":
    import EdgeTangleSet as tset

import pandas as pd
import collections as coll
import colorsys
import matplotlib.pyplot as plt
import igraph as ig
import itertools as iter
import protChecker
import sklearn.metrics as skl
import os
import socket

from py2cytoscape.data.cyrest_client import CyRestClient

from igraph import arpack_options
arpack_options.maxiter=300000

class graphCD():
    def __init__(self, job, log):

        # todo Chuck some error checking in here later ******
        self.job = job
        self.log = log

        self.log.log("{}".format(job['outName']))

        self.doPrint = False

        graph = ig.Graph.Read_Ncol(job['inFile'],names=True, directed=False)

        graph.simplify()
        self.giantComp = graph.clusters().giant()

        print("Number of nodes: {}".format(self.giantComp.vcount()))
        print("Number of edges: {}".format(self.giantComp.ecount()))


    ## Creates a list of colours evenly distributed around the hue spectrum
    # with fixed saturation and value
    def getColourList(self, format="rgb", overlap=False):
        self.colList = []
        if overlap:
            commRange = self.numOverlapComms
        else:
            commRange = self.CDnumComps

        if format=="hsv":
            for i in range(0,commRange):
                self.colList.append(
                    "{:0<4f} {:0<4f} {:0<4f}".format(i/commRange, 1, 0.7))
        elif format=="rgb":
            for i in range(0,commRange):
                rgblist = list(colorsys.hsv_to_rgb(i/commRange, 1, 0.7))
                self.colList.append(
                    "#{}".format("".join(map(lambda x: hex(int(x*255)).\
                        split("0x")[1].zfill(2),rgblist))))
        else:
            print("Incorrect colour format, press any key to continue")
            input()


    def findTangleComms(self):

        #### dep is added to kmin to give *separation* order
        #### add 1 to get *tangle* order
        dep = 4
        if tangleType == "V":
            self.groundset = {"{}_{}".format(self.giantComp.vs[edge.source]['name'], self.giantComp.vs[edge.target]['name']) for edge in self.giantComp.es}
            self.TangleSet = tset.VertexTangleSet(self.giantComp, self.job, self.log)
        elif tangleType == "E":
            self.groundset = set(self.giantComp.vs["name"])
            self.TangleSet = tset.EdgeTangleSet(self.giantComp, self.job, self.log)
            print("Finding Edge connectivity tangles")
        else:
            print("incorrect tangle type {} specified. Use V or E".format(tangleType))

        self.TangleSet.findAllTangles(depth=dep)

        self.assignCommunities(thres = 1)

        if "Yeast" in self.job["outName"]:
            quality = self.evaluateCommunities()
            # todo something with quality

        if self.doPrint and socket.gethostname().find("ginger") > -1:
            self.cytoPrint()

        return(self.giantComp.vcount(), self.giantComp.ecount(), self.TangleSet.getTangleCounts())

    # todo this probably doesn't work for vertex connectivity tangles!
    def assignCommunities(self, thres):
        self.foundcover = pd.DataFrame(index = sorted(self.giantComp.vs["name"]), columns=range(self.TangleSet.currentTangle), dtype=int)
        tangNum = 0
        for order in range(self.TangleSet.kmin, self.TangleSet.kmax+1):
            for tang in self.TangleSet.TangleLists[order]:
                numSeps = len(tang.smallSides)
                # tang.smallSides is list of sets
                smallCounter = coll.Counter([x for s in tang.smallSides for x in s])

                for v in range(self.TangleSet.groundsetSize):
                    self.foundcover.loc[self.TangleSet.names[v], tangNum] = \
                        1 if (1-smallCounter[v]/numSeps >= thres) else 0

                tangNum+=1

        # this makes sure only those communities with at least 3 modes are included.
        self.foundcover = self.foundcover.loc[:, (self.foundcover.sum(axis=0) >= 3)]

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

    def cytoPrint(self, graphToPrint=None):
        if graphToPrint is None:
            graphToPrint = self.giantComp


        cy = CyRestClient()

        cy.session.delete() ### WHy do I need to do this?

        for nodeObj in graphToPrint.vs:
            if "overlapComms" in nodeObj.attributes() and nodeObj['overlapComms'] is not None:
                nodeObj['overlapComms'] = "_".join(map(str, nodeObj['overlapComms']))

        cyGraph = cy.network.create_from_igraph(graphToPrint)

        # layouts = cy.layout.get_all()
        # print("--------------")
        # print(json.dumps(layouts, indent=4))
        #
        # print(cy.layout.apply.__code__.co_varnames)
        # print(dir(cy.layout))
        #
        # print("--------------")

        cy.layout.apply(name='kamada-kawai', network=cyGraph)
        # cy.layout.apply(name='force-directed', network=cyGraph)
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
            # 'NODE_LABEL_COLOR': '#FFFFFF',
            'NODE_LABEL_COLOR': '#000000',
            'EDGE_WIDTH': 1,
            'EDGE_TRANSPARENCY': 255,
            # 'NODE_FILL_COLOR': '#DCDCDC',
        }

        cyStyle.update_defaults(basic_settings)

        # cyStyle.create_passthrough_mapping(column='format',
        #     vp='NODE_CUSTOMGRAPHICS_1', col_type='String')
        cyStyle.create_passthrough_mapping(column='edgecol',
            vp='EDGE_STROKE_UNSELECTED_PAINT', col_type='String')
        # cyStyle.create_passthrough_mapping(column='id',
        #     vp='NODE_LABEL', col_type='String')
        cyStyle.create_passthrough_mapping(column='name',
            vp='NODE_LABEL', col_type='String')
        cyStyle.create_passthrough_mapping(column='flow',
            vp='EDGE_LABEL', col_type='String')
        cy.style.apply(cyStyle, cyGraph)

        cyPDF = cyGraph.get_pdf()

        PDFfilename = '{}/{}-GiantComponentLayout.pdf'.\
            format(self.job['outputFolder'], self.job['outName'])
        PDFoutfile = open(PDFfilename, 'wb')
        PDFoutfile.write(cyPDF)
        PDFoutfile.close()

        os.system("pdfcrop {} {}".format(PDFfilename, PDFfilename))

    def evaluateCommunities(self):
        # see analyseOverlapComms

        quality = coll.defaultdict(float)

        self.protChecker = protChecker.protChecker(self.giantComp.vs['name'])
        NMI = self.protChecker.calculateNMI(self.foundcover)
        print("Found NMI: ", NMI)
        # commQual = self.protChecker.getSimilarityRatio(self.foundcover)
        # overlapMI = self.protChecker.getOverlapMI(self.foundcover)
        # note mi is Mutual inf for (per node) num of comms vs num of GO terms
        # NMI is Normalised mutual inf between assigned comms and comms by GO terms

        # todo - do somthing with the qual measures
        # todo also re-add the coverage ratio
        return(quality)

