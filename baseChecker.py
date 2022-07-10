import pandas as pd
import numpy as np
import itertools as iter

class commChecker():
    def __init__(self, nodeList):
        self.realcover = pd.DataFrame(index=sorted(nodeList))
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

    def getRealCover(self):
        #### Crack shits if haven't run simMatrix first
        return self.realcover


    def calculateNMI(self, coverY):
        coverX = self.realcover

        coverIntersection = coverX.transpose() @ coverY

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
                    HXkYlList.append(HXkYl)
            if len(HXkYlList) == 0:
                HXkbigY = h(PX1) + h(PX0) # H(Xk)
            else:
                HXkbigY = min(HXkYlList)

            HXkbigYNorm = HXkbigY / (h(PX1) + h(PX0)) # H(Xk)
            HXkbigYNormSum += HXkbigYNorm

        HbigXbigYNorm = HXkbigYNormSum / coverX.shape[1]
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


    def getOverlapMI(self, foundcover):
        # todo: discuss whether this is a meaningful measure, as tangles are *not* overlapping
        # todo also make this fucker work too.
        overlapMeta = self.getOverlapMetadata()
        overlapMeta = self.getOverlapMetadata()
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
            .format(self.job['outName'], commQual, mi))
        plt.savefig('{}/{}-ScatterPlotCommVsGO.pdf'.\
            format(self.job['outputFolder'], self.job['outName']),
            bbox_inches='tight')
        plt.close()
        ##### working out NMI (Lancichinetti 2009)
