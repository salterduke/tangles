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


# not implemented
def compareCovers(filename1, filename2):
    cover1 = pd.read_csv(filename1)
    cover2 = pd.read_csv(filename2)




print("---------------------------------------------")
compareSeps("../outputCutfinder/YeastGSCompB-SepList-CF.tsv", "../outputGHU_all/YeastGSCompB-SepList-GHU.tsv")
print("---------------------------------------------")



