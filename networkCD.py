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

    def findOverLapCommunities(self):
        # todo unwrap this to get rid of extras
        self.commLists = coll.defaultdict(set)

        self.graphData = pd.DataFrame({"nodeID": [idx for idx, v in
            enumerate(self.giantComp.vs)], "degree": self.giantComp.degree()})

        self.overlapByTangles()

        if self.doPrint and socket.gethostname().find("ginger") > -1:
            self.cytoPrint()


    def overlapByTangles(self):

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
        self.TangleSet.getTangleCounts()
        # todo consider what return value we want?
        # todo need to actually assign to communities here. Maybe new function? Reads in orientations?
        # remember that in adding comms to nodes, we have to clear comm assignments if want to assign differently later
        # save to file somehow, then read that to evaluate communities

    def addCommToNode(self, nodeIndex, comm):

        nodeObj = self.giantComp.vs[nodeIndex]

        # initialise if necessary
        if 'overlapComms' not in nodeObj.attribute_names() or \
            nodeObj['overlapComms'] is None:
            nodeObj['overlapComms'] = set()

        nodeObj['overlapComms'].add(comm)

        commList = []
        if 'comm0' not in nodeObj.attribute_names():
            for c in range(self.numOverlapComms):
                nodeObj["comm{}".format(c)] = "0"
                commList.append("comm{}".format(c))
            nodeObj['format'] = \
            'piechart: colorlist="{}" showlabels=false attributelist="{}"'.\
                format(",".join(self.colList), ",".join(commList)) ######

        nodeObj["comm{}".format(comm)] = "1"

    def clearCommAssignments(self):
        for node in self.giantComp.vs:
            node['overlapComms'] = None
        self.numOverlapComms = 0
        self.CDnumComps = 0

    def calculateNMI(self, coverX, coverY):

        coverIntersection = coverX.transpose() @ coverY

        HXgivenY = self.conditionalEntropy(coverX, coverY, coverIntersection)
        HYgivenX = self.conditionalEntropy(coverY, coverX, coverIntersection.transpose())
        NMI = 1 - 0.5*(HXgivenY + HYgivenX)
        print("NMI: {}".format(NMI))
        # exit()  ##### ******
        return(NMI)

    def conditionalEntropy(self, coverX, coverY, coverIntersection):
        def h(p):
            # print("p: {}".format(p))
            if p == 0:
                return(0)
            else:
                return(-p*np.log2(p))   ###  ***** Check if correct base????

        N = self.giantComp.vcount()
        coverXsizes = coverX.sum(axis = 0).tolist()
        coverYsizes = coverY.sum(axis = 0).tolist()

        HXkbigYNormSum = 0

        for k in range(coverX.shape[1]):
            HXkYlList = []
            PX1 = coverXsizes[k]/N
            PX0 = 1 - PX1
            if (h(PX1) + h(PX0)) == 0:
                continue
            for l in range(coverY.shape[1]):
                PY1 = coverYsizes[l]/N
                PY0 = 1 - PY1
                PX1Y1 = coverIntersection.iloc[k, l] / N
                PX1Y0 = (coverXsizes[k] - coverIntersection.iloc[k, l]) / N
                PX0Y1 = (coverYsizes[l] - coverIntersection.iloc[k, l]) / N
                PX0Y0 = 1 - PX1Y1 - PX1Y0 - PX0Y1

                if h(PX1Y1) + h(PX0Y0) > h(PX1Y0) + h(PX0Y1):
                    HXkYl = h(PX1Y1) + h(PX0Y0) + h(PX1Y0) + h(PX0Y1) - \
                            h(PY1) - h(PY0)
                    HXkYlList.append(HXkYl)
            if len(HXkYlList) == 0:
                HXkbigY = h(PX1) + h(PX0) # H(Xk)
            else:
                HXkbigY = min(HXkYlList)

            HXkbigYNorm = HXkbigY / (h(PX1) + h(PX0)) # H(Xk)
            HXkbigYNormSum += HXkbigYNorm

        HbigXbigYNorm = HXkbigYNormSum / coverX.shape[1]
        return(HbigXbigYNorm)

    def analyseOverlapComms(self):
        # todo tidy this

        columns = ['node', 'degree', 'commCount']

        graphData = pd.DataFrame(columns = columns)
        # foundcover = pd.DataFrame(index = sorted(self.giantComp.vs["name"]), columns = range(self.numOverlapComms)).fillna(0)
        foundcover = pd.DataFrame(index = sorted(self.giantComp.vs["name"]), dtype=int)

        for node in self.giantComp.vs:
            nodeData = coll.defaultdict()
            if 'name' in node.attribute_names():
                nodeData['node'] = node['name']
            else:
                node['name'] = node.index
                nodeData['node'] = node.index

            nodeData['degree'] = node.degree()

            if 'overlapComms' in node.attribute_names() and \
                node['overlapComms'] is not None:
                nodeData['commCount'] = len(node['overlapComms'])
                for comm in node['overlapComms']:
                    foundcover.loc[node["name"], comm] = 1
            else:
                nodeData['commCount'] = 0
            ## ** Could put something in here adding the list of comms to the
            ## node, but unweildy for large graphs
            graphData = graphData.append(dict(nodeData), ignore_index = True)

        foundcover.fillna(0, inplace = True)
        # print(foundcover)

        # ffs, why should I have to do this? This seems stupid
        graphData['degree'] = pd.to_numeric(graphData['degree'])
        graphData['commCount'] = pd.to_numeric(graphData['commCount'])

        graphData.to_csv("{}/{}-CommDetailsPerNode.csv".\
            format(self.job['outputFolder'], self.job['outName']))

        graphData.set_index("node", inplace=True)
        self.graphData = graphData

        quality = coll.defaultdict(float)


        ###### todo fix problem with protchecker and cutfinder! ****** removed this shit for testing
        # if "Yeast" in self.job["outName"]:
        if False:
            quality = self.evaluateCommunities(foundcover)
        else:
            quality["enrichment"] = 0
            quality["mi"] = 0
            quality["NMI"] = 0

        quality["commcover"] = len(graphData[graphData["commCount"] > 0]) / self.giantComp.vcount()
        quality["overlapcover"] = graphData["commCount"].mean()


        qualFile = "{}/{}quality.csv".format(self.job['outputFolder'], self.job['outName'])
        if not os.path.exists(qualFile):
            text = "numComms\tcommQ\toverlapQ\tNMI\tcommCover\toverlapCover\n"
            with open(qualFile, 'w') as the_file:
                the_file.write(text)

        text = "{:d}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.2f}\t{:.2f}\n".format(
            self.numOverlapComms,
            quality["enrichment"],
            quality["mi"],
            quality["NMI"],
            quality["commcover"],
            quality["overlapcover"]
        )
        print("{}, writing to quality>>>>>>>>>>>>>>>>>>>>>")

        with open(qualFile, 'a') as the_file:
            the_file.write(text)

        return(quality)

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

    # todo sort this out. Lots of todo items!
    def evaluateCommunities(self, foundcover, doOverlap=True):


        quality = coll.defaultdict(float)

        self.protChecker = protChecker.protChecker(self.giantComp.vs['name'])

        realcover = self.protChecker.getRealCover()

        ###### teting this shit *****
        # shortcover = realcover.drop(realcover.columns[10:126], axis=1)
        # self.calculateNMI(shortcover, foundcover)

        NMI = self.calculateNMI(realcover, foundcover)

        # Note: totalSim is only i vs j, not j vs i, does not include i vs i
        similarity = \
            self.protChecker.makeSimMatrix(self.giantComp.vs['name'])
        aveSim = self.protChecker.getAveSim()

        totalCommSim = 0
        numPairs = 0
        for commIndex, commNodes in self.commLists.items():
            if len(commNodes) >= 3:
                for pair in iter.combinations(commNodes, 2):
                    totalCommSim += similarity[pair[0]][pair[1]]
                    numPairs+=1
        aveCommSim = totalCommSim / numPairs
        enrichment = aveCommSim / aveSim

        if doOverlap:
            overlapMeta = self.protChecker.getOverlapMetadata()
            self.graphData = self.graphData.join(overlapMeta, on="node")

            #### get rid of the NANS.
            self.graphData.dropna(inplace = True)

            ### ? Why?
            self.graphData['GOcount'] = pd.to_numeric(self.graphData['GOcount'])

            #### fucking no idea here on number of bins.
            # https://stats.stackexchange.com/questions/179674/number-of-bins-when-computing-mutual-information
            bins = math.floor(math.sqrt(len(self.graphData)/5))
            ## from https://stackoverflow.com/questions/20491028/optimal-way-to-compute-pairwise-mutual-information-using-numpy
            contingencyTable = np.histogram2d(self.graphData["commCount"], self.graphData["GOcount"], bins=bins)[0]
            mi = skl.mutual_info_score(None, None, contingency=contingencyTable)
            print("MI: {}".format(mi))

            ax = self.graphData.plot.scatter(x='GOcount',y='commCount')
            ax.set_xlabel("Num of GO annotations")
            ax.set_ylabel("Community membership count")
            plt.title("Comm count vs GO count: {}, E: {}, MI: {}"\
                .format(self.job['outName'], enrichment, mi))
            plt.savefig('{}/{}-ScatterPlotCommVsGO.pdf'.\
                format(self.job['outputFolder'], self.job['outName']),
                bbox_inches='tight')
            plt.close()
            ##### working out NMI (Lancichinetti 2009)

        quality["enrichment"] = enrichment
        quality["mi"] = mi
        quality["NMI"] = NMI


        return(quality)

