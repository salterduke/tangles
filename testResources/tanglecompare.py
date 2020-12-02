import csv
import pandas as pd
import collections
import sys

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

def compareSeps(filename1, filename2):
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

    print("len old, {} len new {}".format(len(seps2), len(seps1 )))

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
compareSeps("../outputDestupidification/YeastGSCompA-SepList.tsv", "../outputBacktoStupid/YeastGSCompA-SepList.tsv")
print("---------------------------------------------")
