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

    def compareTangles(self, filename1, filename2, N):

        orients = [0,0]
        orients[0] = pd.read_csv(filename1, sep=",")
        orients[1] = pd.read_csv(filename2, sep=",")


        self.groundset = set(range(N))

        # first check column headers
        if not all(orients[0].columns == orients[1].columns):
            print("Column headers don't match")
            print(pd.DataFrame({"file1": orients[0].columns, "file2": orients[1].columns}) )
            print("Testing columns that do match")
            mismatch = np.where(orients[0].columns != orients[1].columns)[0][0]
        else:
            mismatch = len(orients[0].columns) + 1

        # start at col 1 as 0 is index
        tangResults = []
        for id in range(1,mismatch-1):
            cols = [0,0]
            for fid in range(2):
                cols[fid] = orients[fid].iloc[:,id]
                cols[fid] = cols[fid][cols[fid].notna()]

            cols = sorted(cols, key=len)
            shortCol = list(map(self.getSmallSide, cols[0]))
            longCol = list(map(self.getSmallSide, cols[1]))
            colResults = np.zeros(len(longCol), dtype=bool)
            for sid, sideFromLong in enumerate(longCol):
                if sideFromLong in shortCol:
                    colResults[sid] = True
                else:
                    for sideFromShort in shortCol:
                        if sideFromLong.issubset(sideFromShort):
                            colResults[sid] = True
                            break
            tangResults.append(all(colResults))
            # colDF = pd.DataFrame({"seps": longCol, "result": colResults})
            # print(colDF)
            # print(shortCol)
            dummy = 1

        resultsDF = pd.DataFrame({"tangle": orients[0].columns[1:mismatch-1], "result": tangResults})
        print(resultsDF)

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
# comparer.compareTangles("./orientFiles/YeastGSCompB_core-Orientations-Original.csv", "./orientFiles/YeastGSCompB_core-Orientations-Supersets.csv", 37)
comparer.compareTangles("./orientFiles/YeastGSCompB_core-Orientations-Original.csv", "./orientFiles/YeastGSCompB_core-Orientations-Preclude.csv", 37)
print("---------------------------------------------")

