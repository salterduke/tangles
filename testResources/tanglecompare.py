import pandas as pd
import numpy as np
import collections
import sys
import string



# this will have to be a little kludgy:
# a) don't think can guarantee all v's will be given explicitly in orientation file
# so pass in number of v's and make sure to use orientation file with ids, not names
# b) massive pain in the arse to take into account v merging in pre-processing
# so must ensure that both files have the same pre-processing, which means that:
# c) make sure the number of v's passed is the number *after* pre-processing

class comparerClass():

    def compareOrder(self, ord):
        ordOrig = self.origDF.loc[:,("ord={}".format(ord) in self.origDF.colums)]

    def compareTangles(self, origFile, newFile, N):

        self.origDF = pd.read_csv(origFile, sep=",")
        self.newDF = pd.read_csv(newFile, sep=",")

        self.N = N

        self.groundset = set(range(N))

        # first check column headers
        if (self.origDF.shape[1] != self.newDF.shape[1]) or not all(self.origDF.columns == self.newDF.columns):
            print("Column headers don't match")
            # print(pd.DataFrame({"orig": self.origDF.columns, "new": self.newDF.columns}) )
            print("Testing columns that do match")
            mismatch = min(self.origDF.shape[1], self.newDF.shape[1])
            # todo add better heading checking
        else:
            mismatch = len(self.origDF.columns)

        # start at col 1 as 0 is index
        tangResults = []
        for id in range(1,mismatch):
            print(id)
            origCol = self.origDF.iloc[:,id]
            origCol = origCol[origCol.notna()]

            newCol = self.newDF.iloc[:,id]
            newCol = newCol[newCol.notna()]

            origList = list(map(self.getSmallSide, origCol))
            newList = list(map(self.getSmallSide, newCol))
            colResults = np.zeros(len(origList), dtype=bool)
            for sid, sideFromLong in enumerate(origList):
                if sideFromLong in newList:
                    colResults[sid] = True
                else:
                    for sideFromShort in newList:
                        if sideFromLong.issubset(sideFromShort):
                            colResults[sid] = True
                            break
                        else:
                            compBleaves = {32, 36, 5, 8, 12, 15, 20, 22, 24, 25, 29, 30, 31}
                            # todo fix this if need other graph than compB
                            longNoLeaves = sideFromLong - compBleaves
                            longExactLeaves = sideFromLong - longNoLeaves
                            numLeaves = len(longExactLeaves)
                            if longNoLeaves.issubset(sideFromShort):
                                # ie, the main bits match, now check leaves
                                subsetDiff = sideFromShort - longNoLeaves
                                leavesInShort = {v for v in subsetDiff if (v in longExactLeaves or v < 0)}
                                # ie, either one of the specified leaves, or unnamed leaf placeholder
                                if len(leavesInShort) >= numLeaves:
                                    colResults[sid] = True
                                    break


            if not all(colResults):
                dummy=1

            tangResults.append(all(colResults))
            # colDF = pd.DataFrame({"seps": origList, "result": colResults})
            # print(colDF)
            # print(newList)
            dummy = 1

        resultsDF = pd.DataFrame({"tangle": self.origDF.columns[1:mismatch], "result": tangResults})
        print(resultsDF)
        
        # for ord in range(2, 5):
        #     self.compareOrder(ord)

    def getSmallSide(self, sideStr):

        left, right = map(str.strip, sideStr.split("/"))
        if "G\\" in left:
            complement = set(map(int, right.strip("[").strip("]").split(",")))
            smallSet = self.groundset - complement
        else:
            smallSet = set(map(int, left.strip("[").strip("]").split(",")))
        return(smallSet)



print("---------------------------------------------")
comparer = comparerClass()
# comparer.compareTangles("./orientFiles/YeastGSCompB_core-Orientations-Original.csv", "./orientFiles/YeastGSCompB_core-Orientations-smallCheck.csv", 37)
comparer.compareTangles("./orientFiles/YeastGSCompB_core-Orientations-Original.csv", "../outputDevVY/YeastGSCompB_core-Orientations.csv", 37)
# comparer.findDistSeps()
print("---------------------------------------------")

