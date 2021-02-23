import pandas as pd
import collections
import sys

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


# todo compareTangles is not tested!
def compareTangles(filename1, filename2):
    seps1 = pd.read_csv(filename1, sep="\t")
    seps2 = pd.read_csv(filename2, sep="\t")


    groundset = extractVertices(seps1["side1"][0]) | extractVertices(seps1["side2"][0])

    small1 = seps1.loc[seps1["orientation"] == 1]
    small2 = seps1.loc[seps1["orientation"] == 2]
    small = list(small1["side1"]) + list(small2["side2"])

    # missedSeps = seps2.loc[~seps2["side1"].isin(set(seps1["side1"]))]
    missedSeps = set(seps2["side1"]) - set(seps1["side1"])
    missedSeps = [extractVertices(side) for side in missedSeps]
    print("Missed seps before: {}".format(len(missedSeps)))

    print("len old, {} len new {}".format(len(seps2), len(seps1)))

    additionalSmall = []

    ### list of sets
    smallAsSets = [extractVertices(side) for side in small]
    for i in range(len(smallAsSets)):
        for j in range(i + 1, len(smallAsSets)):
            newSep = smallAsSets[i] | smallAsSets[j]
            if newSep in missedSeps:
                missedSeps.remove(newSep)
            elif groundset - newSep in missedSeps:
                print("Removing complement 1")
                missedSeps.remove(newSep)

            additionalSmall.append(newSep)
    print("Missed seps after: {}".format(len(missedSeps)))

    for i in range(len(additionalSmall)):
        for j in range(i + 1, len(additionalSmall)):
            newSep = additionalSmall[i] | additionalSmall[j]
            if newSep in missedSeps:
                missedSeps.remove(newSep)
            elif groundset - newSep in missedSeps:
                print("Removing complement 2")
                missedSeps.remove(newSep)
    print("Missed seps after stage 2: {}".format(len(missedSeps)))
    # print("\n".join([str(sep) for sep in missedSeps]))

    copyMissedSeps = missedSeps.copy()
    # print("\n".join([str(sep) for sep in copyMissedSeps]))
    print("---------------------")
    testSet = smallAsSets + additionalSmall

    for sep in copyMissedSeps:
        # print("testing: {}".format(sep))
        for smallSep in testSet:
            if sep.issubset(smallSep):
                # print(sep)
                print("Found: {}".format(sep))
                missedSeps.remove(sep)
                break

    print("Missed seps after stage 3: {}".format(len(missedSeps)))
    # print("\n".join([str(sep) for sep in missedSeps]))
    printMissed(seps2, missedSeps)


def extractVertices(textStr):
    plain = textStr.replace("[", "").replace("]", "")
    plainSet = set(plain.split(", "))
    plainSet = {item.replace("'", "") for item in plainSet}
    return plainSet

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




# compare("YeastGSComp4-AvoidProps-Old.csv", "YeastGSComp4-AvoidProps-New.csv")
# compare("YeastGSComp4-AvoidProps-New.csv", "YeastGSComp4-AvoidProps-Old.csv")


print("---------------------------------------------")
compareSeps("../outputCutfinder/TinyEdges-SepList-CF.tsv", "../outputGHU_all/TinyEdges-SepList-GHU.tsv")
# compareSeps("../outputCutfinder/YeastGSCompA-SepList-CF.tsv", "../outputGHU_all/YeastGSCompA-SepList-GHU.tsv")
# compareSeps("../outputCutfinder/YeastGSCompB-SepList-CF.tsv", "../outputGHU_all/YeastGSCompB-SepList-GHU.tsv")
print("---------------------------------------------")

#compareSeps("../outputBacktoStupid/TinyEdges-SepList.tsv", "../outputGHU_all/TinyEdges-SepList-GHU.tsv")
#compareSeps("../outputBacktoStupid/YeastGSCompB-SepList.tsv", "../outputGHU_all/YeastGSCompB-SepList-GHU.tsv")
# compareSeps("../outputBacktoStupid/YeastGSCompA-SepList.tsv", "../outputGHU_all/YeastGSCompA-SepList-GHU.tsv")
