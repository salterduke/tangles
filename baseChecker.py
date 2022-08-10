import pandas as pd
import numpy as np
import itertools as iter

class commChecker():
    def __init__(self, nodeList):
        self.nodeNames = nodeList
        self.compareCover = pd.DataFrame(index=sorted(nodeList))
        self.nNodes = len(nodeList)

    # def makeSimMatrix(self):
    #     print("makeSimMatrix must be overridden for each child")
    #     raise NotImplementedError

    def getAveSim(self):
        #### Crack shits if haven't run simMatrix first
        return self.aveSim

    def getOverlapMetadata(self):
        #### Crack shits if haven't run simMatrix first
        return self.overlapMetadata

    def getCompareCover(self):
        #### Crack shits if haven't run simMatrix first
        return self.compareCover


    def calculateNMI(self, coverY, coverX = None):
        if coverX is None:
            coverX = self.compareCover

        try:
            coverIntersection = coverX.transpose() @ coverY
        except:
            dummy = 1

        HXgivenY = self.conditionalEntropy(coverX, coverY, coverIntersection)
        print("HXgivenY: {}".format(HXgivenY))
        HYgivenX = self.conditionalEntropy(coverY, coverX, coverIntersection.transpose())
        print("HYgivenX: {}".format(HYgivenX))
        NMI = 1 - 0.5*(HXgivenY + HYgivenX)
        return(NMI)

    def conditionalEntropy(self, coverX, coverY, coverIntersection):
        def h(p):
            # print("p: {}".format(p))
            if p == 0:
                return(0)
            else:
                return(-p*np.log2(p))   ###  ***** Check if correct base????

        N = coverX.shape[0] # num of vertices
        coverXsizes = coverX.sum(axis = 0).tolist()
        coverYsizes = coverY.sum(axis = 0).tolist()

        HXkbigYNormSum = 0

        # num of comms
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
                    # H(X,Y) - H(Y)
                    HXkYlList.append(HXkYl)
            if len(HXkYlList) == 0:
                HXkbigY = h(PX1) + h(PX0) # H(Xk)
            else:
                HXkbigY = min(HXkYlList)

            HXkbigYNorm = HXkbigY / (h(PX1) + h(PX0)) # H(Xk)
            HXkbigYNormSum += HXkbigYNorm

        try:
            HbigXbigYNorm = HXkbigYNormSum / coverX.shape[1]
        except:
            dummy = 1
        return(HbigXbigYNorm)


    def getSimilarityRatio(self, foundcover):
        # todo: integrate all of this!
        # Note: totalSim is only i vs j, not j vs i, does not include i vs i

        nComms = foundcover.shape[1]

        # similarity = \
        #     self.protChecker.makeSimMatrix(self.giantComp.vs['name'])
        # aveSim = self.protChecker.getAveSim()

        totalCommSim = 0
        numPairs = 0

        # for commIndex, commNodes in self.commLists.items():
        for commIndex in foundcover.columns.to_list():
            commNodes = foundcover.index[foundcover[commIndex]==1].tolist()
            # results_wide.index[results_wide['BoolCol']].tolist()
            if len(commNodes) >= 3:
                for pair in iter.combinations(commNodes, 2):
                    totalCommSim += self.simDF[pair[0]][pair[1]]
                    numPairs+=1
        aveCommSim = totalCommSim / numPairs
        return(aveCommSim / self.aveSim)
        # todo also add per comm values

