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

    def compareCDMethods(self, foundcover, methods = []):
        self.overlapCliquePercolation()
        NMI = self.calculateNMI(foundcover)
        print("NMI: {}".format(NMI))

    def overlapCliquePercolation(self):
        # https://stackoverflow.com/questions/20063927/overlapping-community-detection-with-igraph-or-other-libaries

        minClique = 3
        maxClique = 3

        cliques = list(map(set, self.G.maximal_cliques(min=minClique)))
        # so, each clique is assigned an index in the list.

        for k in range(minClique,maxClique+1):
            filteredCliques = [c for c in cliques if len(c) >= k]

            edgelist = []

            # adding all the edges in each clique
            for clique in filteredCliques:
                for i, j in iter.combinations(clique, 2):
                    edgelist.append((i, j))

            # adding any edges between two cliques with sufficient overlap
            for i, j in iter.combinations(range(len(filteredCliques)), 2):
                # iter.combinations(iterable, 2) returns every possible pair of values
                # in the iterable (irrespective of order), but returns them ordered
                if len(filteredCliques[i].intersection(filteredCliques[j])) >= k-1:
                    edgelist.append((i, j))

            cliqueLinks = ig.Graph(edgelist, directed=False)
            cliqueComps = cliqueLinks.components()

            numComps = len(cliqueComps)
            if numComps == 0:
                break

            self.compareCover = pd.DataFrame(index = sorted(self.G.vs["name"]), columns=range(numComps), dtype=int)
            for col in self.compareCover.columns:
                self.compareCover[col].values[:] = 0

            commIndex = 0
            for id, comp in enumerate(cliqueComps):
                # each component is a list of vertices
                for v in comp:
                    self.compareCover.loc[self.nodeNames[v], id] = 1
                commIndex+=1

            # self.compareCover = self.compareCover.astype(np.int8)
            # this makes sure only those communities with at least 3 modes are included.
            # the astype is necessary as results_wide init as NaNs, which are stored as floats.
            self.compareCover = self.compareCover.loc[:, (self.compareCover.sum(axis=0) >= 3)].astype(np.int8)
