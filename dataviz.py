import numpy as np
import pandas as pd
# import collections as coll
# import colorsys
# from matplotlib import cm
import matplotlib.pyplot as plt
# import igraph as ig
# import itertools as iter
# import sklearn.metrics as skl
# import os
# import socket
import random
import seaborn as sns
sns.relplot(x="NominalEs", y="time", col="order", row="NominalVs", hue="algorithm", data=resSummary)
sns.set()

VYfiles = [
"./Timings/VY_results2022-07-06 22.59.30.385284.csv",
"./Timings/VY_results2022-07-10 21.33.23.724766.csv"
]
YWSfiles = [
"./Timings/YWS_results2022-07-12 05.28.58.377122.csv",
"./Timings/YWS_results2022-07-12 21.24.05.949376.csv"
]


def readAlgData(timingFiles, algName):

    dfList = []

    for id, file in enumerate(timingFiles):
        df = pd.read_csv(file, delimiter=',', header=0, comment="#")
        df["fileID"] = id
        df["algorithm"] = algName
        dfList.append(df)


    results_wide = pd.concat(dfList)

    df1 = results_wide['tangCounts'].str.split('-', expand=True).add_prefix('Order').fillna('')
    results_wide = pd.concat([results_wide, df1], axis = 1)

    df1 = results_wide['timings'].str.split('-', expand=True).add_prefix('Time').fillna('')
    results_wide = pd.concat([results_wide, df1], axis = 1)

    df1 = results_wide['network'].str.split('_', expand=True).add_prefix('Nominal').fillna('')
    results_wide = pd.concat([results_wide, df1], axis = 1)


    for col in [col for col in results_wide.columns if "Order" in str(col)]:
        results_wide[col] = results_wide.apply(lambda x: x[col][1], axis=1)

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

fileLists = [VYfiles, YWSfiles]
algNames = ["VY", "YWS"]

resDFs = []

for id in (0,1):
    resDFs.append(readAlgData(fileLists[id], algNames[id]))

results = pd.concat(resDFs)

shortres = results.loc[results.NominalVs.isin([20, 50, 90])]
# shortres_noOutlier = shortres.loc[shortres.time < 1500]

resSummary = shortres.groupby(["NominalEs", "NominalVs", "order", "algorithm"])["time"].mean().reset_index()

sns.relplot(x="NominalEs", y="time", col="order", row="NominalVs", hue="algorithm", data=resSummary)
plt.show()


