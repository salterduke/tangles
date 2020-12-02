import math
import numpy as np
import scipy as sp
import csv
import pandas as pd
import collections as coll
import igraph as ig
import itertools as it
import tools
import ete3
import sys

import BaseTangleSet as btang

class VertexTangleSet(btang.TangleSet):
    def __init__(self, G, job, log):
        self.G = G
        self.groundsetSize = self.G.ecount()
        self.groundset = {"{}_{}".format(self.G.vs[edge.source]['name'], self.G.vs[edge.target]['name']) for edge in self.G.es}
        super().__init__(job, log)

    def checkTangleAxioms(self, newSep):
        ### since all single edges must be on the small side
        ### if any side has more than m-2 edges, it only takes
        ### two single edges to cover the whole graph.
        if len(newSep) >= self.groundsetSize - 2:
            return False

        return super().checkTangleAxioms(newSep)



    def findNextOrderSeparations(self, k = None):
        if k is None:  ### ie, first time running
            self.TangleTree.add_feature("midsets", [])

            # self.log.tick("minimum_size_separators")
            minMidsets = self.G.minimum_size_separators()
            # self.log.tock()

            print("All midsets:")
            print([[self.G.vs[node]['name'] for node in midset] for midset in minMidsets])

            self.kmin = len(minMidsets[0])

            self.log.tick("processMidset{}".format(self.kmin))
            for midset in minMidsets:
                self.processMidset(midset)
            self.log.tock()
        else:
            print("********************************")
            print("k: {}".format(k))
            self.log.tick("processMidset{}".format(k))
            for ksub in it.combinations(self.G.vs.indices, k):
                self.processMidset(ksub)
            #### *** convert to list? so is ordered?
            self.separations[k] = list(self.separations[k])
            self.definitelySmall[k] = list(self.definitelySmall[k])
            self.log.tock()

    def processMidset(self, midset):
        def sepExists(sep, k):
            for i in range(self.kmin, k):
                if sep in self.separations[i]:
                    return True
            return False

        def defSmallExists(side, k):
            for i in range(self.kmin, k):
                if side in self.definitelySmall[i]:
                    return True
            return False

        def addToSepList(side, midset):
            namedMidset = [self.G.vs[node]['name'] for node in midset]

            sideNodes = [node for node in it.chain.from_iterable(side)]
            sideNodes.extend(namedMidset)
            sideSubGraph = self.G.subgraph(sideNodes)

            size = len(midset)

            sideEdges = [(sideSubGraph.vs[edge.source]['name'], sideSubGraph.vs[edge.target]['name']) for edge in sideSubGraph.es]
            complementSubGraph = self.G - sideEdges
            sideEdgesFormatted = frozenset({"{}_{}".format(sideSubGraph.vs[edge.source]['name'], sideSubGraph.vs[edge.target]['name']) for edge in sideSubGraph.es})
            complementEdgesFormatted = frozenset(self.groundset - sideEdgesFormatted)

            ######## *******
            if all(deg <= 2 for deg in sideSubGraph.degree()):
                if not defSmallExists(sideEdgesFormatted, size):
                    self.definitelySmall[size].append(sideEdgesFormatted)
            elif all(deg <= 2 for deg in complementSubGraph.degree()):
                if not defSmallExists(complementEdgesFormatted, size):
                    self.definitelySmall[size].append(complementEdgesFormatted)
            else:
                separation = frozenset([sideEdgesFormatted, complementEdgesFormatted])
                if not sepExists(separation, len(midset)):
                    self.separations[size].append(separation)

        Gcopy = self.G.copy()
        Gcopy.delete_vertices(midset)
        names = Gcopy.vs['name']
        newComps = []
        for cluster in Gcopy.components():
            newComps.append([names[member] for member in cluster])

        if len(newComps) >= 2:
            for partition in tools.partition(newComps, 2):
                ##### Order these *****???
                for side in partition:
                    addToSepList(side, midset)
            return True
        else:
            return False
