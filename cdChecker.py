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

    def compareCDMethods(self, foundcover, methods = ["CPM","between"]):
        resList = []


        # note, double [[]] in .loc gives df, [] gives series
        for order in range(min(foundcover.loc["order"]), max(foundcover.loc["order"]) + 1):
            orderCover = foundcover.loc[:,foundcover.loc["order"]==order]
            orderCover = orderCover.drop(index="order")

            mm = self.getMembershipFromCover(orderCover)
            # todo look at built in comparison methods

            for method in methods:
                if method == "CPM":
                    for cliqueSize in range(3,7):
                        cpmCover = self.overlapCliquePercolation(cliqueSize)
                        if cpmCover is not None:
                            print("CPM{}".format(cliqueSize))
                            NMI = self.calculateNMI(orderCover, cpmCover, remove_none=True)

                            cpm = self.getMembershipFromCover(cpmCover)
                            NMIcalc = ig.compare_communities(mm, cpm, method="nmi", remove_none=True)
                            NMI_alt = self.calculateNMI_alt(mm, cpm)
                            print(NMI, NMIcalc, NMI_alt)
                            resList.append({"method": "CPM{}".format(cliqueSize),
                                            "order": order,
                                            "NMI": NMI})
                        else:
                            break
                else:
                    if method == "between":
                        print("between")
                        dendro = self.G.community_edge_betweenness(directed=False)
                        cover = self.getCoverFromDendro(dendro)
                    else:
                        print("Unknown method: ", method)

                    NMI = self.calculateNMI(orderCover, cover)
                    resList.append({"method": method,
                                    "order": order,
                                    "NMI": NMI})



        resDF = pd.DataFrame(resList)
        dummy = 1
        # print to file here

    def getMembershipFromCover(self, cover):
        if any(cover.sum(axis="columns") > 1):
            print("Ooops, node assigned to too many communities")

        membership = [None]*self.G.vcount()

        for vid in self.G.vs.indices:
            if any(cover.loc[self.G.vs[vid]["name"], :]==1):
                membership[vid] = \
                    cover.loc[:, cover.loc[self.G.vs[vid]["name"], :] == 1].columns[0]

        return membership

    def getCoverFromDendro(self, dendro):
        membership = np.array(dendro.as_clustering().membership)

        cover = pd.DataFrame(index = sorted(self.G.vs["name"]), columns=range(max(membership)), dtype=int)
        for col in cover.columns:
            cover[col].values[:] = 0

        for id in range(max(membership)):
            members = np.where(membership == id)[0]
            cover.loc[self.G.vs[members.tolist()]["name"], id] = 1

        # cpmCover = cpmCover.astype(np.int8)
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

        # adding all the edges in each clique
        for clique in cliques:
            for i, j in iter.combinations(clique, 2):
                edgelist.append((i, j))

        # adding any edges between two cliques with sufficient overlap
        for i, j in iter.combinations(range(len(cliques)), 2):
            # iter.combinations(iterable, 2) returns every possible pair of values
            # in the iterable (irrespective of order), but returns them ordered
            if len(cliques[i].intersection(cliques[j])) >= cliqueSize-1:
                edgelist.append((i, j))

        cliqueLinks = ig.Graph(edgelist, directed=False)
        cliqueComps = cliqueLinks.components()

        numComps = len(cliqueComps)
        if numComps == 0:
            return None

        cpmCover = pd.DataFrame(index = sorted(self.G.vs["name"]), columns=range(numComps), dtype=int)
        for col in cpmCover.columns:
            cpmCover[col].values[:] = 0

        commIndex = 0
        for id, comp in enumerate(cliqueComps):
            # each component is a list of vertices
            for v in comp:
                cpmCover.loc[self.nodeNames[v], id] = 1
            commIndex+=1

        # cpmCover = cpmCover.astype(np.int8)
        # this makes sure only those communities with at least 3 modes are included.
        # the astype is necessary as results_wide init as NaNs, which are stored as floats.
        cpmCover = cpmCover.loc[:, (cpmCover.sum(axis=0) >= 3)].astype(np.int8)
        return cpmCover