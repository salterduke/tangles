import pandas as pd
from collections import defaultdict
import csv
import numpy as np
import baseChecker as bch
import itertools as iter
import igraph as ig
import Modules.CliquePercolationMethod as cpm
import Modules.tools as tools

class cdChecker(bch.commChecker):
    def __init__(self, G):
        self.G = tools.pruneToStubs(G)
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


    def compareCDMethods(self, foundcover, othercover = None,
                         methods = ["between", "fastgreedy", "infomap", "labelprop", "eigen",
                                    "leiden", "multilevel", "modularity", "spinglass", "walktrap"]):
        resList = [] # will be list of dicts, then convert to DF
        mshipList = []

        # optimal modularity doesn't run for larger graphs, so remove
        if self.G.vcount() > 100:
            methods = [m for m in methods if m != "modularity"]

        if othercover is not None:
            methods = methods + ["otherCover"]

        for method in methods:
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
            elif method == "infomap":
                cluster = self.G.community_infomap()
                CD_mship = self.getMembershipFromClustering(cluster)
            elif method == "labelprop":
                cluster = self.G.community_label_propagation()
                CD_mship = self.getMembershipFromClustering(cluster)
            elif method == "eigen":
                cluster = self.G.community_leading_eigenvector()
                CD_mship = self.getMembershipFromClustering(cluster)
            elif method == "leiden":
                cluster = self.G.community_leiden()
                CD_mship = self.getMembershipFromClustering(cluster)
            elif method == "multilevel":
                cluster = self.G.community_multilevel()
                CD_mship = self.getMembershipFromClustering(cluster)
            elif method == "modularity":
                cluster = self.G.community_optimal_modularity()
                CD_mship = self.getMembershipFromClustering(cluster)
            elif method == "spinglass":
                cluster = self.G.community_spinglass()
                CD_mship = self.getMembershipFromClustering(cluster)
            elif method == "walktrap":
                dendro = self.G.community_walktrap()
                CD_mship = self.getMembershipFromDendro(dendro)
            elif "CPM" in method:
                # method should be CPM3, CPM4, etc, so just get last char
                # todo add error checking
                cliqueSize = int(method[3])
                CD_cover = self.overlapCliquePercolation(cliqueSize)
                CD_mship1 = self.getMembershipFromCover(CD_cover)
                # testing
                commList = cpm.clique_percolation_method(self.G, k=cliqueSize)
                CD_mship = self.getMembershipFromCommList(commList)

                CD_cover_his = self.getCoverFromCommList(commList)

                print("Mine is first, his is second, method {}".format(method))
                self.compareCovers(CD_cover, CD_cover_his)

                dummy = 1
            else:
                print("Unknown method: {}".format(method))

            mshipList.append({"method": method,
                              "mship": CD_mship
            })

            # note, double [[]] in .loc gives df, [] gives series
            for order in range(min(foundcover.loc["order"]), max(foundcover.loc["order"]) + 1):
                orderCover = foundcover.loc[:, foundcover.loc["order"] == order]
                orderCover = orderCover.drop(index="order")
                Tangle_mship = self.getMembershipFromCover(orderCover)

                for metric in ("vi", "nmi", "split-join", "rand", "adjusted_rand"):
                    try:
                        value = ig.compare_communities(Tangle_mship, CD_mship, method=metric, remove_none=False)
                    except:
                        print("moocow")

                    # print(order, method, metric, value)
                    resList.append({"order": order,
                                    "method": method,
                                    "metric": metric,
                                    "value": value})

        resDF = pd.DataFrame(resList)
        mshipDF = pd.DataFrame(mshipList)
        return resDF, mshipDF

    def getMembershipFromCover(self, cover):
        if any(cover.sum(axis="columns") > 1):
            print("Ooops, node assigned to too many communities")

        noneID = cover.shape[1]
        # all vertices not assigned to any comm should still be included in the comparison,
        # as "unassigned" is still information. Therefore give them all the same id so None's not removed
        membership = [noneID]*self.G.vcount()

        for vid in self.G.vs.indices:
            try:
                if any(cover.loc[self.G.vs[vid]["name"], :]==1):
                    membership[vid] = \
                        cover.loc[:, cover.loc[self.G.vs[vid]["name"], :] == 1].columns[0]
            except:
                print(vid)

        return membership

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
        # not sure about data types.  If it works, probably not important.
        membership = np.array(dendro.as_clustering().membership)
        noneID = max(membership) + 1
        membership = [noneID if v is None else v for v in membership]
        return membership

    def getCoverFromCommList(self, commList):
        cover = pd.DataFrame(index=sorted(self.G.vs["name"]), columns=range(len(commList)), dtype='Int64').fillna(0)

        # for col in cover.columns:
        #     cover[col].values[:] = 0

        try:
            for id, comm in enumerate(commList):
                cover.loc[self.G.vs[list(comm)]["name"], id] = 1
        except:
            dummy = 1

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

        for id, comp in enumerate(cliqueComps):
            # print(comp)
            # each component is a list of vertices
            for cid in comp:
                # print(cliques[cid])
                for vid in cliques[cid]:
                    # print(self.nodeNames[vid])
                    cpmCover.loc[self.nodeNames[vid], id] = 1

        # the astype is necessary as results_wide init as NaNs, which are stored as floats.
        # try:
        #     cpmCover = cpmCover.astype(np.int8)
        # except:
        #     dummy = 1
        return cpmCover


if __name__ == '__main__':

    dataSets = {
        # "Karate": ("../NetworkData/SmallNWs/KarateEdges.csv","Karate-TangNodes.csv"),
        "YeastB": ("../NetworkData/BioDBs/YeastPPI/YuEtAlGSCompB.csv","YeastGSCompB_core-TangNodes.csv"),
        "YeastA": ("../NetworkData/BioDBs/YeastPPI/YuEtAlGSCompA.csv","YeastGSCompA-TangNodes.csv"),
        # "Celegans": ("../NetworkData/Celegans/NeuronConnect.csv","Celegans-TangNodes.csv"),
        "Jazz": ("../NetworkData/MediumSize/Jazz.csv","Jazz-TangNodes.csv")
        # COMMA !!!!!!!!!!
        # "Copperfield": ("../NetworkData/MediumSize/Copperfield.csv","Copperfield-TangNodes.csv"),
        # "Football": ("../NetworkData/MediumSize/Football.csv","Football-TangNodes.csv"),
        # "Bsubtilis": ("../NetworkData/BioDBs/HINTformatted/BacillusSubtilisSubspSubtilisStr168-htb-hq.txt","BSubtilis-htb-TangNodes.csv")
    }
    coverFolder = "./outputdevResults_VY/"

    # "Karate-TangNodes.csv",
    # "YeastGSCompB_core-TangNodes.csv",
    # "YeastGSCompA-TangNodes.csv",
    # "Celegans-TangNodes.csv",
    # "Jazz-TangNodes.csv",
    # "Copperfield-TangNodes.csv",
    # "Football-TangNodes.csv",
    # "BSubtilis-htb-TangNodes.csv"

    allComparisons = []

    for dataName, dataFiles in dataSets.items():
        graphFile = dataFiles[0]
        G = ig.Graph.Read_Ncol(graphFile, names=True, directed=False)
        checker = cdChecker(G)

        colTypes = defaultdict(lambda:int, {0:str})
        fullFilename = coverFolder + dataFiles[1]
        foundcover = pd.read_csv(fullFilename, index_col=0, dtype=colTypes)
        foundcover.columns = foundcover.columns.astype(int)

        checkresults, mships = checker.compareCDMethods(foundcover)
        checkresults["dataName"] = dataName
        allComparisons.append(checkresults)
        dummy = 1

    comparisonsDF = pd.concat(allComparisons)
    comparisonsDF.to_csv("{}ComparisonValues.csv".format(coverFolder))


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
