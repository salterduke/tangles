import pandas as pd
from collections import defaultdict
import csv
import numpy as np
import baseChecker as bch
import itertools as iter
import igraph as ig
import Modules.CliquePercolationMethod as cpm
import Modules.tools as tools
import cdlib
import platform

class cdChecker(bch.commChecker):
    def __init__(self, G):
        before = pd.Series(G.vs["name"])

        self.oldVs = G.vcount()

        self.G = tools.pruneToStubs(G)
        after = pd.Series(self.G.vs["name"])
        bch.commChecker.__init__(self, self.G.vs["name"])

    def compareCovers(self, cov1, cov2):
        # print(otherFile)
        fvo = tools.matchCols(cov1, cov2)
        # print(fvo)
        if all(fvo):
            print("All comms from first cover found in second cover")
        else:
            unfound = np.where(np.logical_not(fvo))[0]
            print("Column ids {} from first cover missing from second cover".format(unfound))
            # print("Column names {} from first cover missing from second cover".format(list(cov1.columns[unfound])))
        ovf = tools.matchCols(cov2, cov1)
        # print(ovf)
        if all(ovf):
            print("All comms from second cover found in first cover")
        else:
            unfound = np.where(np.logical_not(ovf))[0]
            print("Column ids {} from second cover missing from first cover".format(unfound))
            # print("Column names {} from second cover missing from first cover".format(list(cov1.columns[unfound])))

        if all(fvo) and all(ovf):
            return True
        else:
            return False

    def compareCDMethods(self, foundcover, othercover = None, inheritParent = False,
                         methods = ["between", "fastgreedy", "infomap", "labelprop", "eigen",
                                    "leiden", "multilevel", "modularity", "spinglass", "walktrap"]):
        resList = [] # will be list of dicts, then convert to DF
        mshipList = []
        refValList = []

        modularityList = []
        objFunctionsList = []

        # optimal modularity doesn't run for larger graphs, so remove
        if self.G.vcount() > 100:
            methods = [m for m in methods if m != "modularity"]

        if othercover is not None:
            methods = methods + ["otherCover"]

        expandedVsDF = pd.DataFrame()

        # for method in [method in methods if "CPM" not in method]:
        for method in methods:
            print("Running method {}".format(method))
            if method == "otherCover":
                orderOther = othercover.loc[:, othercover.loc["order"] == order]
                orderOther = orderOther.drop(index="order")
                CD_mship = self.getMembershipFromCover(orderOther)
            elif method == "between":
                dendro = self.G.community_edge_betweenness(directed=False)
                CD_mship = self.getMembershipFromDendro(dendro)
            elif method == "fastgreedy":
                dendro = self.G.community_fastgreedy()
                CD_mship = self.getMembershipFromDendro(dendro)
                modularityList.append({
                    "method": method,
                    "modularity": self.G.modularity(CD_mship)
                })
            elif method == "infomap":
                cluster = self.G.community_infomap()
                infoCodelength = cluster.codelength
                objFunctionsList.append({
                    "objFunction": "infomap",
                    "value": infoCodelength
                })
                CD_mship = self.getMembershipFromClustering(cluster)
            elif method == "labelprop":
                cluster = self.G.community_label_propagation()
                CD_mship = self.getMembershipFromClustering(cluster)
            elif method == "eigen":
                cluster = self.G.community_leading_eigenvector()
                CD_mship = self.getMembershipFromClustering(cluster)
                modularityList.append({
                    "method": method,
                    "modularity": self.G.modularity(CD_mship)
                })
            elif method == "leiden":
                cluster = self.G.community_leiden()
                CD_mship = self.getMembershipFromClustering(cluster)
            elif method == "multilevel":
                cluster = self.G.community_multilevel()
                CD_mship = self.getMembershipFromClustering(cluster)
                modularityList.append({
                    "method": method,
                    "modularity": self.G.modularity(CD_mship)
                })
                # yes, I know it'd redundant; I don't give a shit.
                objFunctionsList.append({
                    "objFunction": method,
                    "value": self.G.modularity(CD_mship)
                })
            elif method == "modularity":
                cluster = self.G.community_optimal_modularity()
                CD_mship = self.getMembershipFromClustering(cluster)
                modularityList.append({
                    "method": method,
                    "modularity": self.G.modularity(CD_mship)
                })
            elif method == "spinglass":
                cluster = self.G.community_spinglass(update_rule="simple")
                CD_mship = self.getMembershipFromClustering(cluster)
            elif method == "walktrap":
                dendro = self.G.community_walktrap()
                CD_mship = self.getMembershipFromDendro(dendro)
            elif "CPM" in method:
                pass  # will handle separately, but don't want to crack the sads
            else:
                print("Unknown method: {}".format(method))

            if "CPM" not in method:
                mshipList.append({"method": method,
                                  "mship": CD_mship
                })
                expandedVsDF[method] = self.expandMshipList(CD_mship)

            ############################# new line here!
            self.noneID = self.G.vcount() - 1
            last_Mship = [self.noneID]*self.G.vcount()
            last_Expanded = [-1]*self.oldVs
            # if any vs unassigned even at lowest order, they should *stay* unassigned

            # print(dataName)
            # for order in range(min(foundcover.loc["order"]), max(foundcover.loc["order"]) + 1):
            for order in np.unique(foundcover.loc["order"]):
                orderCover = foundcover.loc[:, foundcover.loc["order"] == order]
                orderCover = orderCover.drop(index="order")
                Tangle_mship, expandedMship = self.getMembershipFromCover(orderCover, expandVs = True, useNegative = True)
                # print("Order {}----------------------------------------------".format(order))
                # print(Tangle_mship)
                Tangle_mship, expandedMship = self.dealwithUnassigned(Tangle_mship, last_Mship, expandedMship, last_Expanded, inheritParent)
                last_Mship = Tangle_mship
                last_Expanded = expandedMship
                # print(Tangle_mship)

                Tangle_commList = self.getCommListFromMship(Tangle_mship)
                expandedVsDF[order] = expandedMship

                if "CPM" in method:
                    # method should be CPM3, CPM4, etc, so just get last char
                    # todo add error checking
                    cliqueSize = int(method[3])
                    commList_his = cpm.clique_percolation_method(self.G, k=cliqueSize)

                    for metric in ("omega", "LFK", "MGH"):
                        value = self.compareOverlapping(commList_his, Tangle_commList, metric)
                        resList.append({"order": order,
                                        "method": method,
                                        "metric": metric,
                                        "value": value})
                else:  # CD method other than CPM
                    for metric in ("vi", "nmi", "split-join", "rand", "adjusted_rand"):
                        value = ig.compare_communities(Tangle_mship, CD_mship, method=metric, remove_none=False)
                        resList.append({"order": order,
                                        "method": method,
                                        "metric": metric,
                                        "value": value})

                        referenceVal = self.getReferenceVal(metric)
                        # at this stage *always doing these two algs - too complicated to do otherwise
                        refValList.append({"method1": "infomap",
                                           "method2": "multilevel",
                                           "metric": metric,
                                           "value": referenceVal})


        referenceDF = pd.DataFrame(refValList).drop_duplicates()
        resDF = pd.DataFrame(resList)
        modDF = pd.DataFrame(modularityList)
        objDF = pd.DataFrame(objFunctionsList)

        mshipDF = pd.DataFrame(mshipList)
        return resDF, modDF, mshipDF, expandedVsDF, referenceDF, objDF

    def dealwithUnassigned(self, thisMship, lastMship, thisExp, lastExp, inheritParent=True):

        if inheritParent:

                for i in range(len(thisMship)):
                    if thisMship[i] == -1:
                        thisMship[i] = lastMship[i]
                for i in range(len(thisExp)):
                    if thisExp[i] == -1:
                        thisExp[i] = lastExp[i]
        else:
            # have to use non negative for compare_communities, but want -1 in expanded for output
            for i in range(len(thisMship)):
                if thisMship[i] == -1:
                    thisMship[i] = self.noneID

        return thisMship, thisExp



    def getReferenceVal(self, metric):
        # at this stage at least *always* doing same two algs - too complicated to do otherwise
        # "between", "fastgreedy", "infomap", "labelprop", "eigen",
        # "leiden", "multilevel", "modularity", "spinglass", "walktrap"

        InfomapCluster = self.G.community_infomap()
        InfomapMship = self.getMembershipFromClustering(InfomapCluster)

        MultilevelCluster = self.G.community_multilevel()
        MultilevelMship = self.getMembershipFromClustering(MultilevelCluster)

        refVal = ig.compare_communities(InfomapMship, MultilevelMship, method=metric, remove_none=False)
        return refVal


    def compareOverlapping(self, commList1, commList2, metric):

        # need to create new community for unassigned nodes
        for commList in (commList1, commList2):
            assignedNodes = {j for i in commList for j in i}
            newComm = [vid for vid in self.G.vs.indices if vid not in assignedNodes]
            commList.append(newComm)

        commList1 = [comm for comm in commList1 if len(comm) > 0] # probably a neater way of doing this, but eh.
        commList2 = [comm for comm in commList2 if len(comm) > 0] # and probably a way of doing it in the loop, but...

        # convert to NodeClustering
        clustering1 = cdlib.NodeClustering(commList1, self.G, overlap=True)
        clustering2 = cdlib.NodeClustering(commList2, self.G, overlap=True)

        if metric == "omega":
            value = cdlib.evaluation.omega(clustering1, clustering2)
        elif metric == "LFK":
            value = cdlib.evaluation.overlapping_normalized_mutual_information_LFK(clustering1, clustering2)
        elif metric == "MGH":
            value = cdlib.evaluation.overlapping_normalized_mutual_information_MGH(clustering1, clustering2)
        else:
            input("invalid metric specified {}".format(metric))

        return value.score

    def getCommListFromMship(self, mship):
        nComms = max(mship)+1
        commList = [ [] for _ in range(nComms) ]    # create list of empty lists
        for vid, commID in enumerate(mship):
            commList[commID].append(vid)
        return commList

    def expandMshipList(self, membership):
        unmergedVs = [v for mergedv in self.G.vs["name"] for v in mergedv.split(";")]
        unmergedVs.sort()
        vSeries = pd.Series(index=unmergedVs, dtype=int)
        for vid, mergedv in enumerate(self.G.vs["name"]):
            # I *think* it's safe to enumerate rather than calling .indices explicitly....
            for v in mergedv.split(";"):
                vSeries[v] = membership[vid]
        vSeries = vSeries.astype(int)
        return(vSeries)
        # vSeries = vSeries.replace(noneID, -1)

    def getMembershipFromCover(self, cover, expandVs = False, useNegative = False):
        if any(cover.sum(axis="columns") > 1):
            print("Ooops, node assigned to too many communities")

        if useNegative:
            localNone = -1
        else:
            localNone = self.noneID

        # all vertices not assigned to any comm should still be included in the comparison,
        # as "unassigned" is still information. Therefore give them all the same id so None's not removed
        membership = [localNone]*self.G.vcount()

        for vid in self.G.vs.indices:
            if any(cover.loc[self.G.vs[vid]["name"], :]==1):
                membership[vid] = \
                    cover.loc[:, cover.loc[self.G.vs[vid]["name"], :] == 1].columns[0]

        if expandVs:
            # yes I know this is redundant...
            unmergedVs = [v for mergedv in self.G.vs["name"] for v in mergedv.split(";")]
            unmergedVs.sort()
            vSeries = pd.Series(index = unmergedVs, dtype=int)
            for vid, mergedv in enumerate(self.G.vs["name"]):
                # I *think* it's safe to enumerate rather than calling .indices explicitly....
                for v in mergedv.split(";"):
                    vSeries[v] = membership[vid]
            vSeries = vSeries.astype(int)
            vSeries = vSeries.replace(localNone, -1)
        else:
            vSeries = None

        return membership, vSeries

    def getMembershipFromCommList(self, commList):
        noneID = len(commList)
        membership = [noneID] * self.G.vcount()
        for commID, comm in enumerate(commList):
            for vid in comm:
                if membership[vid] == noneID:
                    membership[vid] = commID
                else:
                    print("What the hell, vid already assigned!")
        return membership

    def getMembershipFromModTuple(self, modTuple):
        membership = modTuple[0]
        noneID = max(membership) + 1
        membership = [noneID if v is None else v for v in membership]
        return membership

    def getMembershipFromClustering(self, cluster):
        membership = cluster.membership
        noneID = max(membership) + 1
        membership = [noneID if v is None else v for v in membership]
        return membership

    def getMembershipFromDendro(self, dendro):
        membership = np.array(dendro.as_clustering().membership)
        noneID = max(membership) + 1
        membership = [noneID if v is None else v for v in membership]
        return membership

    def getCoverFromCommList(self, commList):
        cover = pd.DataFrame(index=sorted(self.G.vs["name"]), columns=range(len(commList)), dtype='Int64').fillna(0)
        for id, comm in enumerate(commList):
            cover.loc[self.G.vs[list(comm)]["name"], id] = 1
        return cover


    def getCoverFromDendro(self, dendro):
        membership = np.array(dendro.as_clustering().membership)

        cover = pd.DataFrame(index = sorted(self.G.vs["name"]), columns=range(max(membership)), dtype=int)
        for col in cover.columns:
            cover[col].values[:] = 0

        for id in range(max(membership)):
            members = np.where(membership == id)[0]
            cover.loc[self.G.vs[members.tolist()]["name"], id] = 1

        # this makes sure only those communities with at least 3 modes are included.
        # the astype is necessary as results_wide init as NaNs, which are stored as floats.
        # cover = cover.loc[:, (cover.sum(axis=0) >= 3)].astype(np.int8)
        cover = cover.astype(np.int8)
        return cover

    # todo NOTE currently has error in commList - DO NOT USE!
    def overlapCliquePercolation(self, cliqueSize):
        # https://stackoverflow.com/questions/20063927/overlapping-community-detection-with-igraph-or-other-libaries

        cliques = list(map(set, self.G.maximal_cliques(min=cliqueSize)))
        # so, each clique is assigned an index in the list.

        edgelist = []

        # This bit is flat-out wrong. Was an attempt to fix error with comms with a single clique being ignored,
        # but I got confused between edges in G and edges in the clique graph
        # replaced with ***** below
        # adding all the edges in each clique
        # for clique in cliques:
        #     for i, j in iter.combinations(clique, 2):
        #         edgelist.append((i, j))


        # adding any edges between two cliques with sufficient overlap
        for i, j in iter.combinations(range(len(cliques)), 2):
            # iter.combinations(iterable, 2) returns every possible pair of values
            # in the iterable (irrespective of order), but returns them ordered
            if len(cliques[i].intersection(cliques[j])) >= cliqueSize - 1:
                edgelist.append((i, j))

        # *****
        cliqueLinks = ig.Graph()
        cliqueLinks.add_vertices(list(range(0, len(cliques))))
        cliqueLinks.add_edges(edgelist)
        # the above, instead of the below line
        # cliqueLinks = ig.Graph(edgelist, directed=False)

        cliqueComps = cliqueLinks.components()

        numComps = len(cliqueComps)
        if numComps == 0:
            return None

        cpmCover = pd.DataFrame(index=self.G.vs["name"], columns=range(numComps), dtype='Int64').fillna(0)

        commList = []

        for id, comp in enumerate(cliqueComps):
            # each component is a list of vertices in the new clique graph, NOT original graph
            singleComm = set()
            for cid in comp:
                singleComm.update(cliques[cid])
                for vid in cliques[cid]:
                    cpmCover.loc[self.nodeNames[vid], id] = 1
            commList.append(singleComm)

        return cpmCover, commList

    def getInterestingOrders(self, cover):
        return cover.loc[:, cover.loc["order"].duplicated(keep=False)]


if __name__ == '__main__':



    dataSets = {
        # "Karate": ("../NetworkData/SmallNWs/KarateEdges.csv","Karate-TangNodes.csv"),
        # "YeastB": ("../NetworkData/BioDBs/YeastPPI/YuEtAlGSCompB.csv","YeastGSCompB_core-TangNodes.csv"),
        # "YeastA": ("../NetworkData/BioDBs/YeastPPI/YuEtAlGSCompA.csv","YeastGSCompA-TangNodes.csv"),
        # "Celegans": ("../NetworkData/Celegans/NeuronConnect.csv","Celegans-TangNodes.csv"),
        # "Jazz": ("../NetworkData/MediumSize/Jazz.csv","Jazz-TangNodes.csv"),
        "Copperfield": ("../NetworkData/MediumSize/Copperfield.csv","Copperfield-TangNodes.csv"),
        # "Football": ("../NetworkData/MediumSize/Football.csv","Football-TangNodes.csv"),
        "Bsubtilis": ("../NetworkData/BioDBs/HINTformatted/BacillusSubtilisSubspSubtilisStr168-htb-hq.txt","BSubtilis-htb-TangNodes.csv"),
        "Iceland": ("../NetworkData/MediumSize/Iceland.csv", "Iceland-TangNodes.csv"),
        "Zebra": ("../NetworkData/MediumSize/Zebra.csv", "Zebra-TangNodes.csv"),
        "Dolphins": ("../NetworkData/MediumSize/Dolphins.csv", "Dolphins-TangNodes.csv")
    }
    coverFolder = "./outputdevResults_VY/"
    subfolder = "inheritComparisons/"
    # subfolder = ""
    outfileLabel = "Disj"

    if platform.system() != "Linux":
        for dname in dataSets.keys():
            dataSets[dname] = (
                dataSets[dname][0].replace("../", "C:/Users/mrousset/Documents/PhDThesisLaptop/Code/"),
                dataSets[dname][1]
            )

    # todo
    # for methodClass in ("CPM", "Disj"):
    for methodClass in ("Disj"):
        for inheritParent in (True, False):
            allComparisons = []
            allModularities = []
            allMships = []
            allRefVals = []
            allobjVals = []

            for dataName, dataFiles in dataSets.items():
                print("Running data {}".format(dataName))
                graphFile = dataFiles[0]
                G = ig.Graph.Read_Ncol(graphFile, names=True, directed=False)
                G = G.connected_components().giant()
                checker = cdChecker(G) # should condense nodes in __init__

                colTypes = defaultdict(lambda:int, {0:str})
                fullFilename = coverFolder + dataFiles[1]
                foundcover = pd.read_csv(fullFilename, index_col=0, dtype=colTypes)
                foundcover.columns = foundcover.columns.astype(int)

                # make sure vertices match
                graphVSeries = pd.Series(checker.G.vs["name"]).sort_values(ignore_index=True)
                coverVSeries = pd.Series(foundcover.index)
                order = foundcover.loc["order"]
                coverVSeries = coverVSeries[coverVSeries != "order"]
                coverVSeries = coverVSeries.sort_values(ignore_index=True)
                # tweak so sorted *after* removing order so doesn't get sorted into the midddle.
                if not graphVSeries.equals(coverVSeries):
                    print("vertices don't match up for datafile {} ----------------------------.".format(dataName))
                    continue
                else:
                    print("yay vertices match")
                    # stuff it, I couldn't be stuffed to avoid the redundancy here
                    sortOrder = pd.Series(checker.G.vs["name"])
                    foundcover = foundcover.loc[sortOrder]
                    foundcover.loc["order"] = order
                    # The order in foundcover now matches the order in the graph.

                foundcover = checker.getInterestingOrders(foundcover) # remove orders where only one tangle
                if len(foundcover.columns) != 0:
                    if inheritParent:
                        outfileLabel = "{}_inherit".format(methodClass)
                    else:
                        outfileLabel = "{}_negatives".format(methodClass)
                    if methodClass == "Disj":
                        checkresults, modularityVals, mships, expMships, refVals, objVals = checker.compareCDMethods(foundcover, inheritParent = inheritParent)
                    elif methodClass == "CPM":
                        methods = ["CPM3", "CPM4", "CPM5"]
                        checkresults, modularityVals, mships, expMships, refVals, objVals  = checker.compareCDMethods(foundcover, inheritParent = inheritParent, methods = methods)
                    else:
                        exit("crack the sads, invalid methodClass {}".format(methodClass))

                    expMships.to_csv("{}{}_expandedTangleComms{}.csv".format(coverFolder+subfolder, dataName, outfileLabel))
                    checkresults["dataName"] = dataName
                    modularityVals["dataName"] = dataName
                    refVals["dataName"] = dataName
                    objVals["dataName"] = dataName
                    mships["dataName"] = dataName
                    allComparisons.append(checkresults)
                    allModularities.append(modularityVals)
                    allMships.append(mships)
                    allRefVals.append(refVals)
                    allobjVals.append(objVals)

            comparisonsDF = pd.concat(allComparisons)
            comparisonsDF.to_csv("{}ComparisonValues{}_All.csv".format(coverFolder+subfolder, outfileLabel))
            modularityDF = pd.concat(allModularities)
            refDF = pd.concat(allRefVals)
            objDF = pd.concat(allobjVals)
            objDF.to_csv("{}objectiveFns.csv".format(coverFolder+subfolder))
            refDF.to_csv("{}referenceMetrics_{}.csv".format(coverFolder+subfolder, outfileLabel))
            if modularityDF.size > 0:
                modularityDF.to_csv("{}ModularitiesAll.csv".format(coverFolder+subfolder))
                modularityDF.groupby(["dataName"])["modularity"].max()
                modularityDF[modularityDF['modularity'] == modularityDF.groupby(['dataName'])['modularity'].transform(max)]
            mshipsDF = pd.concat(allMships)
            mshipsDF.to_csv("{}Memberships.csv".format(coverFolder+subfolder))



        # This stuff is for comparing two covers to ensure they're the same.
        # Not what we want to do to get CD comparison results
        # otherCoverFiles = [
        #     "outputDevYWS/YeastGSCompB_core-TangNodes.csv",
        #     "outputDevVY/YeastGSCompB_core-TangNodes.csv"
        # ]
        #
        # for otherFile in otherCoverFiles:
        #     othercover = pd.read_csv(otherFile, index_col=0, dtype=colTypes)
        #     othercover.columns = othercover.columns.astype(int)
        #
        #     # checkresults, mships = checker.compareCDMethods(foundcover, othercover, methods=["CPM3", "CPM4"])
        #     checkresults, mships = checker.compareCDMethods(foundcover, methods=["CPM3", "CPM4"])
