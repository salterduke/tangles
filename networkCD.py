import numpy as np
import ete3

import EdgeTangleSet_VY
import EdgeTangleSet_YWS


import pandas as pd
from matplotlib import cm
import igraph as ig
import Modules.dataviz as dataviz
from Modules import cdChecker
import random
import os

class graphCD():
    def __init__(self, job, log):

        # todo Chuck some error checking in here later ******
        self.job = job
        self.log = log

        self.log.log("{}".format(job['outName']))

        # parse bonus arguments
        for arg in job.get("args", "").split(","):
            if len(arg) > 0:
                if "=" in arg:
                    argname, argval = arg.split("=")
                    if argval.isnumeric():
                        job[argname] = int(argval)
                    else:
                        job[argname] = argval
                else:
                    print("Args not properly specified. Need to be arg1=val1,arg2=val2 etc")
                    print(arg)
                    exit()

        try:
            fileExt = job["inFile"].split(".")[-1]
            if fileExt in ("png", "ico", "jpg", "jpeg", "bmp"):
                job["doImage"] = True
                job["imType"] = "IMAGE"
            if "MNIST" in job["outName"].upper():
                job["doImage"] = True
                job["imType"] = "MNIST"
                job["MNISTid"] = job["outName"].split(".")[-1]
        except:
            print("moocow")
            exit("Error checking file extension")

        if "construct" in job and job["construct"]:
            graph = self.constructRandom(job)
            fname = "../NetworkData/Constructed/{}.ncol".format(job["outName"])
            graph.write_ncol(fname)
            with open("const.txt", 'a') as the_file:
                the_file.write("{};{}".format(fname, job["outName"]))
        elif "doImage" in job and job["doImage"]:
            print(job["outName"])
            graph = job["imParser"].fetchSingleImage(job)
            print(job["outName"])
            dummy = 1
        else:
            graph = ig.Graph.Read_Ncol(job['inFile'], names=True, directed=False)

        if not graph.is_weighted():
            graph.es["weight"] = 1

        graph.simplify(combine_edges="sum")
        self.giantComp = graph.clusters().giant()

        self.protChecker = None
        self.cdChecker = None
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


    def findTangleComms(self, dep = 4, maxEmptyOrders = 4, sepsOnly = False):
        # edited to find dep total levels
        # added maxEmptyOrders, so can still limit number of partial seps that must be added to heap
        # (looser limit though)
        # but need to keep eye out for if this is an issue.

        self.groundset = set(self.giantComp.vs["name"])
        if "YWS" in self.job["testName"]:
            self.TangleSet = EdgeTangleSet_YWS.EdgeTangleSet(self.giantComp, self.job, self.log)
        else:
            self.TangleSet = EdgeTangleSet_VY.EdgeTangleSet(self.giantComp, self.job, self.log)

        timings, sepCounts = self.TangleSet.findAllTangles(depth=dep, maxEmptyOrders=maxEmptyOrders, sepsOnly=sepsOnly)

        if not sepsOnly:
            self.assignCommunities(rule="dist")
            # if not self.job.get("doImage"):
            #     quality = self.evaluateCommunities()
            # todo add evaluation for images

        self.doPrint = True
        if self.doPrint:
            self.igPrint()

        return(self.giantComp.vcount(), self.giantComp.ecount(), self.TangleSet.getTangleCounts(), timings, sepCounts)

    # also note currently only does *once* at end for all k levels -
    # could change to do each level?
    def assignCommunities(self, rule="dist"):
        if rule == "dist":
            self.assignCommunitiesDistSeps()
        else:
            self.assignCommunitiesAllSeps() # not currently working, but ignore

        # self.printCommTree()


    def assignCommunitiesDistSeps(self):
        self.foundcover = pd.DataFrame(index = sorted(self.giantComp.vs["name"]), columns=range(self.TangleSet.currentTangle), dtype=int)
        # we are going to work out which seps are distinguishing
        # I don't think we can mark them off as they're added, as we might be adding prematurely
        self.distinguishingSeps = []

        sepOrders = []

        nonDistSeps = []
        distSepsOneSided = []
        # Note that these orders are the order of the *separations*, not the tangles.
        tangNum = 0

        try:
            for order in range(self.TangleSet.kmin, self.TangleSet.kmax+1):
                donothing = 1
        except:
            print("moocow")

        for order in range(self.TangleSet.kmin, self.TangleSet.kmax+1):
            # print("Order: -------------------- {}".format(order))
            for side in self.TangleSet.separations[order]:
                # print("Side: {}".format(side))
                # if side == {6,7,8}:
                #     print("Found trouble sep")
                #     print("moocow")
                sideIn = False
                compIn = False
                comp = set(self.giantComp.vs.indices) - side
                for tang in self.TangleSet.TangleLists[order]:
                    # print(tang.smallSides)
                    # if side in tang.smallSides:
                    if any(side.issubset(smallside) for smallside in tang.smallSides):
                        sideIn = True
                        if compIn:
                            break
                            # stop looking altogether
                        else:
                            continue
                            # stop looking in this tangle, keep looking for comp
                    # elif comp in tang.smallSides:
                    elif any(comp.issubset(smallside) for smallside in tang.smallSides):
                        compIn = True
                        if sideIn:
                            break
                        else:
                            continue
                if sideIn and compIn:
                    # print("Adding dist sep")
                    # add both for easier checking
                    self.distinguishingSeps.append(side)
                    self.distinguishingSeps.append(comp)
                    distSepsOneSided.append([self.giantComp.vs[v]["name"] for v in side])
                elif sideIn and not compIn:
                    nonDistSeps.append([self.giantComp.vs[v]["name"] for v in side])
                elif not sideIn and compIn:
                    nonDistSeps.append([self.giantComp.vs[v]["name"] for v in comp])
                else:
                    print("What the hell - neither side in there")


            for tang in self.TangleSet.TangleLists[order]:

                sepOrders.append(order)

                distSmallSides = []
                for sep in tang.smallSides:
                    for distSep in self.distinguishingSeps:
                        # if sep supset of distsep, then we know which way distsep has to be pointing
                        # note that if a dist sep is not in smallSides, it must be a subset of another sep -
                        if sep.issuperset(distSep):
                            distSmallSides.append(distSep)

                if len(distSmallSides) > 0:
                    onAllBig = set(self.giantComp.vs.indices) - set.union(*distSmallSides)
                else:
                    onAllBig = set(self.giantComp.vs.indices)

                for v in range(self.TangleSet.groundsetSize):
                    self.foundcover.loc[self.TangleSet.names[v], tangNum] = \
                        1 if (v in onAllBig) else 0

                tangNum+=1
        nonDistDF = pd.DataFrame({"nonDist": nonDistSeps})
        outfile = "{}/{}-NonDist.csv". \
            format(self.job['outputFolder'], self.job['outName'])
        nonDistDF.to_csv(outfile)

        distDF = pd.DataFrame({"Dist": distSepsOneSided})
        outfile = "{}/{}-Dist.csv". \
            format(self.job['outputFolder'], self.job['outName'])
        distDF.to_csv(outfile)


        outfile = "{}/{}-TangNodes.csv". \
            format(self.job['outputFolder'], self.job['outName'])
        # making copy to add row for orders without messing with anything else.
        # cover_copy = cover_copy.append(pd.Series(sepOrders, index=self.foundcover.columns, name="order"), ignore_index=False)
        # not sure why the above was removed, but eh.
        # the astype is necessary as results_wide init as NaNs, which are stored as floats.
        self.foundcover = self.foundcover.astype(np.int8)
        self.foundcover.loc[len(self.foundcover)] = [ord+1 for ord in sepOrders]
        # as tangle orders are +1 to the highest order seps in them

        self.foundcover.index.values[len(self.foundcover)-1] = "order"
        #
        print("Found {} tangles total".format(self.foundcover.shape[1]))
        # self.job["depth"]
        self.foundcover.to_csv(outfile)
        # I *think* we can keep the "order" row in and it won't bugger anything up,
        # except yeast?


    def assignCommunitiesAllSeps(self):
        self.foundcover = pd.DataFrame(index = sorted(self.giantComp.vs["name"]), columns=range(self.TangleSet.currentTangle), dtype=int)

        sepOrders = []

        # Note that these orders are the order of the *separations*, not the tangles.
        tangNum = 0
        for order in range(self.TangleSet.kmin, self.TangleSet.kmax+1):
            # print("Order: -------------------- {}".format(order))

            for tang in self.TangleSet.TangleLists[order]:

                # probably not necessary, but helps my brain in comparison to dist version
                allSmallSides = list(tang.smallSides)

                if len(allSmallSides) > 0:
                    onAllBig = set(self.giantComp.vs.indices) - set.union(*allSmallSides)
                else:
                    onAllBig = set(self.giantComp.vs.indices)


                for v in range(self.TangleSet.groundsetSize):
                    self.foundcover.loc[self.TangleSet.names[v], tangNum] = \
                        1 if (v in onAllBig) else 0

                tangNum+=1
                sepOrders.append(order)



        outfile = "{}/{}-TangNodesAllSeps.csv". \
            format(self.job['outputFolder'], self.job['outName'])
        # the astype is necessary as results_wide init as NaNs, which are stored as floats.
        self.foundcover = self.foundcover.astype(np.int8)


        self.foundcover.loc[len(self.foundcover)] = [ord+1 for ord in sepOrders]
        # as tangle orders are +1 to the highest order seps in them

        # For if we want to remove communities with zero nodes assigned.
        # not currently working correctly.
        # self.foundcover = self.foundcover.loc[:, (self.foundcover.sum(axis=0) > 0)]

        self.foundcover.index.values[len(self.foundcover)-1] = "order"
        #
        self.foundcover.to_csv(outfile)
        # I *think* we can keep the "order" row in and it won't bugger anything up,
        # except yeast?

    def printCommTree(self):
        # This bit is working with the tangle tree to assign colours to the comms.
        # note that traverse default is BFS, which is what we want - deeper comms will
        # overwrite colours for parent comms.
        # todo if change to every level, need to do on *copy* of tangletree.
        for node in self.giantComp.vs:
            node["color"] = "#ffffff"

        # for treenode in self.TangleSet.TangleTree.traverse():
        #     if "T" not in treenode.name:
        #         # keeps only complete tangles
        #         # print("deleting non-tangle treenodes")
        #         treenode.delete(prevent_nondicotomic=False)
        #     elif int(treenode.name.replace("T", "")) not in self.foundcover.columns:
        #         # keeps only tangles with >= 3 nodes
        #         # print("deleting trivial tangle treenodes")
        #         treenode.delete(prevent_nondicotomic=False)
        #     else:
        #         # these are the actual communities we want to colour
        #         commIndex = int(treenode.name.replace("T", ""))
        #         treedep = treenode.get_distance(self.TangleSet.TangleTree)
        #         for nodeName in self.foundcover.index[self.foundcover[commIndex]==1].tolist():
        #             if nodeName != "order":
        #                 self.giantComp.vs.find(nodeName)["color"] = self.getColour(treedep)


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
            node.dist /= 4

        ts.branch_vertical_margin = 8
        ts.scale = 360

        try:
            # Note that this shit doesn't work correctly in Python 3.10. It's a known issue.
            self.TangleSet.TangleTree.render(outfile, tree_style=ts)

        except Exception as rendError:
            print("Render doesn't work correctly in Python 3.10. Use 3.9 or lower.")
            # print(rendError)


    def igPrint(self, doGrid=True):

        if self.cdChecker is None:
            self.cdChecker = cdChecker.cdChecker(self.giantComp)

        plotCover = self.foundcover.loc[:, self.foundcover.loc["order"].duplicated(keep=False)]

        tangPlotFolder = "{}/tanglePlots/".format(self.job['outputFolder'])
        if not os.path.isdir(tangPlotFolder):
            print("About to create dir {}".format(tangPlotFolder))
            os.makedirs(tangPlotFolder)
        plotter = dataviz.Grapher(outputFolder=tangPlotFolder)

        for order in np.unique(plotCover.loc["order"]):
            orderCover = plotCover.loc[:, plotCover.loc["order"] == order]
            orderCover = orderCover.drop(index="order")
            mship, expMship = self.cdChecker.getMembershipFromCover(orderCover, expandVs=True, useNegative=True)

            outname = "{}_{}cols".format(self.job["outName"], self.job["numColours"])
            plotter.plotSingleClustering(dataName = outname,
                                         clusteringName = str(order),
                                         G = self.giantComp,
                                         clustering = expMship,
                                         inherit=0,
                                         doGrid=doGrid,
                                         gridWidth=self.job["crop2"])

    def evaluateCommunities(self):
        # see analyseOverlapComms (currently commented out)

        # quality = coll.defaultdict(float)

        if self.cdChecker is None:
            self.cdChecker = cdChecker.cdChecker(self.giantComp)

        # todo if re-adding Yeast check, do I need to remove "order" row from foundcover?
        # if "Yeast" in self.job["outName"]:
        #     if self.protChecker is None:
        #         self.protChecker = protChecker.protChecker(self.giantComp.vs['name'])
        #     quality["NMI"] = self.protChecker.calculateNMI(self.foundcover)
        #     print("Found NMI: ", quality["NMI"])
        #     quality["commQual"] = self.protChecker.getSimilarityRatio(self.foundcover)
        #     print("Found commQual: ", quality["commQual"])
        #     # NMI is Normalised mutual inf between assigned comms and comms by GO terms

        # quality = self.cdChecker.compareCDMethods(self.foundcover, methods=["CPM3", "CPM4"])


        # todo - do somthing with the qual measures
        # todo also re-add the coverage ratio?
        return(quality)



