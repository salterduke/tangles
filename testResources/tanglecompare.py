import pandas as pd
import numpy as np
import collections
import sys
import string

# this is the main function
def compareSeps(filename1, filename2):
    seps1 = pd.read_csv(filename1, sep="\t")
    seps2 = pd.read_csv(filename2, sep="\t")

    cuts1 = seps1["cut"].map(extractCuts)
    cuts2 = seps2["cut"].map(extractCuts)

    kmin1 = min(seps1["order"])
    kmax1 = max(seps1["order"])

    kmin2 = min(seps2["order"])
    kmax2 = max(seps2["order"])

    if kmin1 != kmin2:
        sys.exit("min order {} in {} not equal to min order {} in {}".format(kmin1, filename1, kmin2, filename2))

    # for k in range(kmin1, max(kmax1, kmax2)+1):
    for k in range(kmin1, min(kmax1, kmax2)+1):
        kcuts1 = set(cuts1[cuts1.map(len) == k])
        kcuts2 = set(cuts2[cuts2.map(len) == k])
        in1not2 = kcuts1 - kcuts2
        in2not1 = kcuts2 - kcuts1
        if len(in1not2) != 0:
            print("Order {} cuts:".format(k))
            print("{} Cuts in {} but not in {}".format(len(in1not2), filename1, filename2))
            printCuts(in1not2)
        if len(in2not1) != 0:
            print("Order {} cuts:".format(k))
            print("{} Cuts in {} but not in {}".format(len(in2not1), filename2, filename1))
            printCuts(in2not1)
        if len(in1not2) == 0 and len(in2not1) == 0:
            print("Order {} cuts all appear the same:".format(k))

# todo printCuts - print out nicely. Also maybe add logic for flagging matching orders or summarising
def printCuts(cuts):
    # cuts is a set
    for cut in cuts:
        print(str(cut)
              .replace("frozenset", "")
              .replace("{'", "'")
              .replace("'}","'")
              .replace("({","{")
              .replace("})","}")
        )

def extractCuts(cutStr):
    plain = cutStr.replace("[", "").replace("]", "").replace("(", "").replace(" ", "")
    plainSet = set(plain.split("),"))
    plainSet = {item.replace("'", "").replace(")","") for item in plainSet}
    setOfSets = frozenset(frozenset(item.split(",")) for item in plainSet)
    return(setOfSets)


def compareCovers(filename1, filename2):
    cover1 = pd.read_csv(filename1)
    cover2 = pd.read_csv(filename2)



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

        compBsingletons = [32, 36, 5, 8, 12, 15, 20, 22, 24, 25, 29, 30, 31]
        # this is a kludge! fix properly later

        # start at col 1 as 0 is index
        tangResults = []
        for id in range(1,mismatch):
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
                        for single in compBsingletons:
                            newSide = sideFromShort | {single}
                            if sideFromLong.issubset(newSide):
                                colResults[sid] = True
                                break
                        if colResults[sid]:
                            break
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


# todo extractEdges
def extractEdges(textStr):
    plain = textStr.replace("[", "").replace("]", "")
    plainSet = set(plain.split(", "))
    plainSet = {item.replace("'", "") for item in plainSet}
    return plainSet



def printMissed(fullData, seps):
    for sep in seps:
        sepStr = "['" + "', '".join(sorted(sep)) + "']"
        # print(sepStr)
        ord = int(fullData.loc[fullData["side1"] == sepStr]["order"])
        print("{}: {}".format(ord, sepStr))

def compare(filename1, filename2):
    tset1 = pd.read_csv(filename1)
    columns1 = list(tset1)

    tset2 = pd.read_csv(filename2)
    columns2 = list(tset2)

    if columns1 == columns2:
        print("So far so good, headers are equal")

    dict1 = collections.defaultdict(list)

    print(columns1)
    #### skip first col - just indices
    for col in columns1[1:]:
        newset = set()
        order = col.split("=")[1]
        tangle = tset1[col]
        for item in tangle:
            if str(item).find("1.00") > -1:
                newset.add(item)
        dict1[order].append(newset)
        if len(newset) == 0:
            print("Huh, no 1s")
            print(tangle)



    for col in columns2[1:]:
        order = col.split("=")[1]

        testset = set()
        tangle = tset2[col]
        for item in tangle:
            if str(item).find("1.00") > -1:
                testset.add(item)

        if len(testset) == 0:
            print("Huh, no 1s")
            print(tangle)


        if testset in dict1[order]:
            print("Yay, found {} from file2 in file1".format(col))
            # print(testset)
        else:
            print("Bummer, did not find {} from file2 in file1".format(col))
            print(testset)



print("---------------------------------------------")
comparer = comparerClass()
# comparer.compareTangles("./orientFiles/YeastGSCompB_core-Orientations-Original.csv", "./orientFiles/YeastGSCompB_core-Orientations-smallCheck.csv", 37)
comparer.compareTangles("./orientFiles/YeastGSCompB_core-Orientations-Original.csv", "../outputDevYWS/YeastGSCompB_core-Orientations.csv", 37)
# comparer.findDistSeps()
print("---------------------------------------------")

