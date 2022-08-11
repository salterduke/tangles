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

    def calculateNMI_alt(self, mshipX, mshipY):
        keepIDs = [i for i in range(len(mshipX)) if mshipX[i] is not None and mshipY[i] is not None]
        mshipX = [mshipX[i] for i in keepIDs]
        mshipY = [mshipY[i] for i in keepIDs]
        df = pd.DataFrame(index=sorted(pd.unique(mshipX)), columns=sorted(pd.unique(mshipY)))

        denom = len(mshipX)
        for r in df.index:
            for c in df.columns:
                df.loc[r, c] = sum(1 for i in range(denom) if mshipX[i] == r and mshipY[i] == c)

        miSum = 0
        HXsum = 0
        HYsum = 0
        for r in df.index:
            for c in df.columns:
                pxy = df.loc[r, c] / denom
                px = sum(df.loc[r, :]) / denom
                py = sum(df.loc[:, c]) / denom
                if pxy != 0:
                    term = pxy * np.log2(pxy / (px * py))
                else:
                    term = 0
                miSum += term

        for r in df.index:
            px = sum(df.loc[r, :]) / denom
            if px != 0:
                HXsum+=px*np.log2(px)
            else:
                HXsum+=0

        for c in df.columns:
            py = sum(df.loc[:, c]) / denom
            if py != 0:
                HYsum+=py*np.log2(py)
            else:
                HYsum+=0

        NMI = 2*miSum / (-HXsum - HYsum)

        return NMI

    def calculateNMI(self, coverX, coverY = None, remove_none=False):
        if coverY is None:
            coverY = self.compareCover

        if remove_none:
            # remove any vertices not assigned to any comm, to match ig.compare_communities()
            coverY = coverY.loc[(coverY.sum(axis=1) != 0) & (coverX.sum(axis=1)!=0), :]
            coverX = coverX.loc[(coverY.sum(axis=1) != 0) & (coverX.sum(axis=1)!=0), :]

        HXgivenY = self.conditionalEntropy(coverX, coverY)
        print("HXgivenY: {}".format(HXgivenY))
        HYgivenX = self.conditionalEntropy(coverY, coverX)
        print("HYgivenX: {}".format(HYgivenX))
        NMI = 1 - 0.5*(HXgivenY + HYgivenX)
        return(NMI)

    def conditionalEntropy(self, coverX, coverY):
        def h(p):
            # print("p: {}".format(p))
            if p == 0:
                return(0)
            else:
                return(-p*np.log2(p))   ###  ***** Check if correct base????

        coverIntersection = coverX.transpose() @ coverY


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

