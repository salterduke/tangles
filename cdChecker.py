import pandas as pd
from collections import defaultdict
import csv
import numpy as np
import baseChecker as bch
import itertools as iter
import igraph as ig

class cdChecker(bch.commChecker):
    def __init__(self, G):
        bch.commChecker.__init__(self, G.vs["name"])
        self.G = G

    def compareCDMethods(self, foundcover, methods = ["between", "fastgreedy", "infomap", "labelprop", "eigen", "leiden", "multilevel", "modularity", "spinglass", "walktrap"]):
        resList = [] # will be list of dicts, then convert to DF
        mshipList = []

        # optimal modularity doesn't run for larger graphs, so remove
        if self.G.vcount() > 100:
            methods = [m for m in methods if m != "modularity"]

        # note, double [[]] in .loc gives df, [] gives series
        for order in range(min(foundcover.loc["order"]), max(foundcover.loc["order"]) + 1):
            orderCover = foundcover.loc[:,foundcover.loc["order"]==order]
            orderCover = orderCover.drop(index="order")

            Tangle_mship = self.getMembershipFromCover(orderCover)
            # todo look at built in comparison methods

            for method in methods:
                if method == "between":
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
                elif method == "CPM":
                    pass
                else:
                    print("Unknown method: ", method)

                mshipList.append({"method": method,
                                  "mship": CD_mship
                })

                for metric in ("vi", "nmi", "split-join", "rand", "adjusted_rand"):
                    value = ig.compare_communities(Tangle_mship, CD_mship, method=metric, remove_none=False)

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
            if any(cover.loc[self.G.vs[vid]["name"], :]==1):
                membership[vid] = \
                    cover.loc[:, cover.loc[self.G.vs[vid]["name"], :] == 1].columns[0]

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

    # CPM method removed, as error, and take too long to fix and validate.
    # For code, see master branch (I think?) or commit previous to this