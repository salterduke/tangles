import pandas as pd


class commChecker():
    def __init__(self, nodeList):
        self.realcover = pd.DataFrame(index=sorted(nodeList))
        self.size = len(nodeList)

    def makeSimMatrix(self):
        print("makeSimMatrix must be overridden for each child")
        raise NotImplementedError

    def getAveSim(self):
        #### Crack shits if haven't run simMatrix first
        return self.aveSim

    def getOverlapMetadata(self):
        #### Crack shits if haven't run simMatrix first
        return self.overlapMetadata

    def getRealCover(self):
        #### Crack shits if haven't run simMatrix first
        return self.realcover
