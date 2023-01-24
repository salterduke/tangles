import numpy as np
import pandas as pd
# import collections as coll
# import colorsys
# from matplotlib import cm
import matplotlib.pyplot as plt
# import igraph as ig
# import itertools as iter
# import sklearn.metrics as skl
from sklearn.linear_model import LinearRegression
# import os
# import socket
import random
import seaborn as sns

sns.set()

VYfiles = [
"../outputTestVY/results2022-11-23 20.49.11.331334.csv",
"../outputTestVY/results2022-11-25 00.56.39.342710.csv",
"../outputTestVY/results2022-11-25 13.58.17.443580.csv",
"../outputTestVY/results2022-11-26 22.46.30.592230.csv",
"../outputTestVY/results2022-12-10 09.10.23.870668.csv",
"../outputTestVY/results2022-12-10 12.38.50.019390.csv",
"../outputTestVY/results2022-12-11 10.43.25.282799.csv",
"../outputTestVY/results2022-12-11 14.12.00.441343.csv",
"../outputTestVY/results2022-12-12 02.14.09.467554.csv",
"../outputTestVY/results2022-12-12 06.07.45.566873.csv",
"../outputTestVY/results2022-12-12 12.19.58.915419.csv",
"../outputTestVY/results2022-12-12 15.48.30.065434.csv",
"../outputTestVY/results2022-12-13 15.02.59.722227.csv",
"../outputTestVY/results2022-12-13 19.05.02.100191.csv",
"../outputTestVY/results2023-01-04 18.50.10.710700.csv",
"../outputTestVY/results2023-01-03 17.54.11.479864.csv",
"../outputTestVY/results2022-12-27 18.20.50.612801.csv",
"../outputTestVY/results2022-12-26 23.03.21.097428.csv",
"../outputTestVY/results2022-12-19 16.52.11.097076.csv",
"../outputTestVY/results2022-12-18 22.15.15.806828.csv"

]

YWSfiles = [
"../outputTestYWS/results2022-11-24 04.13.04.126166.csv",
"../outputTestYWS/results2022-11-25 09.06.28.692849.csv",
"../outputTestYWS/results2022-11-28 17.05.09.511741.csv",
"../outputTestYWS/results2022-11-29 19.19.01.456783.csv",
"../outputTestYWS/results2022-12-10 09.17.26.935472.csv",
"../outputTestYWS/results2022-12-10 13.08.02.281325.csv",
"../outputTestYWS/results2022-12-11 10.50.25.246702.csv",
"../outputTestYWS/results2022-12-11 14.41.15.343567.csv",
"../outputTestYWS/results2022-12-12 02.21.42.212079.csv",
"../outputTestYWS/results2022-12-12 06.38.50.450858.csv",
"../outputTestYWS/results2022-12-12 12.27.02.689464.csv",
"../outputTestYWS/results2022-12-12 16.17.42.518543.csv",
"../outputTestYWS/results2022-12-13 15.10.50.769288.csv",
"../outputTestYWS/results2022-12-13 19.37.09.424597.csv",
"../outputTestYWS/results2023-01-08 20.27.45.931906.csv",
"../outputTestYWS/results2023-01-04 01.38.36.360110.csv",
"../outputTestYWS/results2022-12-30 10.29.12.984751.csv",
"../outputTestYWS/results2022-12-27 06.33.52.986415.csv",
"../outputTestYWS/results2022-12-22 12.04.04.355848.csv",
"../outputTestYWS/results2022-12-19 05.22.53.195299.csv"
]

def readAlgData(timingFiles, algName):

    dfList = []

    for id, file in enumerate(timingFiles):
        # print("Reading ", file)
        df = pd.read_csv(file, delimiter=',', header=0, comment="#")
        df["fileID"] = id
        df["algorithm"] = algName
        dfList.append(df)


    results_wide = pd.concat(dfList, ignore_index=True)

    df1 = results_wide['tangCounts'].str.split('-', expand=True).add_prefix('Order').fillna('')
    results_wide = pd.concat([results_wide, df1], axis = 1)

    df1 = results_wide['timings'].str.split('-', expand=True).add_prefix('Time').fillna('')
    results_wide = pd.concat([results_wide, df1], axis = 1)

    df1 = results_wide['network'].str.split('_', expand=True).add_prefix('Nominal').fillna('')
    results_wide = pd.concat([results_wide, df1], axis = 1)

    # reswide_before = results_wide.copy(deep = True)
    results_wide.drop('Order4', axis=1, inplace=True)

    try:
        for col in [col for col in results_wide.columns if "Order" in str(col)]:
            results_wide[col] = results_wide.apply(lambda x: x[col][1], axis=1)
    except:
        print("moocow")

    for i in range(1,10):
        results_wide["k{}".format(i)] = np.nan

    orders = [col for col in results_wide.columns if "Order" in col]
    times = [col for col in results_wide.columns if "Time" in col]
    nOrders = len(orders)

    for rid, row in results_wide.iterrows():
        for i in range(nOrders):
            # row["k{}".format(row["Order{}".format(i)])] = row["Time{}".format(i)]
            results_wide.loc[rid, "k{}".format(row["Order{}".format(i)])] = row["Time{}".format(i)]

    results_wide = results_wide.drop(columns=["tangCounts", "timings", "Nominal0"] + orders + times)
    results_wide = results_wide.dropna(how='all', axis='columns')
    results_wide = results_wide.rename(columns = {"Nominal1": "NominalVs", "Nominal2": "NominalEs"})

    results = pd.wide_to_long(results_wide, ["k"], i=["fileID", "network"], j="order")
    results = results.rename(columns = {"k": "time"})
    results = results.reset_index()
    results.time = pd.to_numeric(results.time)
    results.Vs = pd.to_numeric(results.Vs)
    results.Es = pd.to_numeric(results.Es)
    results.NominalVs = pd.to_numeric(results.NominalVs)
    results.NominalEs = pd.to_numeric(results.NominalEs)

    return results

# ------------------------------------------------------------
def plotResults(results):
    results["nm"] = results["Vs"] * results["Es"]

    doLog = "log"
    # doLog = "linear"

    # # for vs in (20,50,100):
    for vs in (150, 200):
        singleVs = results.loc[results.NominalVs == vs].groupby(["NominalEs", "NominalVs", "order", "algorithm"])[
            "time"].mean().reset_index()
        # singleVs = singleVs[singleVs["order"] < 5]
        sns.relplot(x="NominalEs", y="time", col="order", col_wrap=2, hue="algorithm", data=singleVs)
        plt.yscale(doLog)
        plt.savefig("./Timings/Vertices_{}_{}.png".format(doLog, vs))
        # plt.savefig("./Timings/Vertices_{}.pdf".format(vs))

    for ord in (range(2, 6)):
        singleOrd = results.loc[results.order == ord].groupby(["NominalEs", "NominalVs", "order", "algorithm"])[
            "time"].mean().reset_index()
        sns.relplot(x="NominalEs", y="time", col="NominalVs", col_wrap=2, hue="algorithm", data=singleOrd)
        plt.yscale(doLog)
        plt.savefig("./Timings/Order_{}_{}.png".format(doLog, ord))

    for ord in (range(2, 6)):
        singleOrd = results.loc[results.order == ord].groupby(["NominalEs", "NominalVs", "nm", "order", "algorithm"])[
            "time"].mean().reset_index()
        sns.relplot(x="nm", y="time", col="NominalVs", col_wrap=2, hue="algorithm", data=singleOrd)
        plt.yscale(doLog)
        plt.savefig("./Timings/Order_{}_nm_{}.png".format(doLog, ord))

# ------------------------------------------------------------
def regressResults(results):
    pass

# ------------------------------------------------------------


fileLists = [VYfiles, YWSfiles]
algNames = ["VY", "YWS"]
resDFs = []

for id in (0,1):
    resDFs.append(readAlgData(fileLists[id], algNames[id]))

results = pd.concat(resDFs)

# plotResults(results)
regressResults(results)

# nomDiffs = results.groupby(["network"]).agg({
#     "NominalVs":"mean",
#     "Vs":"mean",
#     "NominalEs":"mean",
#     "Es":"mean"
# })

# shortres = results.loc[results.NominalVs.isin([20, 50, 90])]
# # shortres_noOutlier = shortres.loc[shortres.time < 1500]
#
# lowOrder = shortres.loc[shortres.order<=4].groupby(["NominalEs", "NominalVs", "order", "algorithm"])["time"].mean().reset_index()
# highOrder = shortres.loc[shortres.order>=5].groupby(["NominalEs", "NominalVs", "order", "algorithm"])["time"].mean().reset_index()
#
# sns.relplot(x="NominalEs", y="time", col="order", row="NominalVs", hue="algorithm", data=lowOrder)
# plt.savefig("./Timings/lowOrder.png")
#
# sns.relplot(x="NominalEs", y="time", col="order", row="NominalVs", hue="algorithm", data=highOrder)
# plt.savefig("./Timings/highOrder.png")
#

