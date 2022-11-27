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
sns.set()

VYfiles = [
'../outputTestVY/results2022-07-14 20.38.05.058843.csv',
'../outputTestVY/results2022-07-16 07.30.43.685643.csv',
'../outputTestVY/results2022-07-17 04.04.03.003248.csv',
'../outputTestVY/results2022-07-24 21.31.01.705614.csv',
'../outputTestVY/results2022-07-25 21.26.21.847385.csv'
]

YWSfiles = [
'../outputTestYWS/results2022-07-14 03.14.23.538113.csv',
'../outputTestYWS/results2022-07-15 07.42.25.115347.csv',
'../outputTestYWS/results2022-07-17 21.17.23.678404.csv',
'../outputTestYWS/results2022-07-20 02.30.23.855924.csv',
'../outputTestYWS/results2022-07-26 21.18.09.023190.csv'
]


def readAlgData(timingFiles, algName):

    dfList = []

    for id, file in enumerate(timingFiles):
        print("Reading ", file)
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

for vs in (20,50,100):
    singleVs = results.loc[results.NominalVs == vs].groupby(["NominalEs", "NominalVs", "order", "algorithm"])["time"].mean().reset_index()
    sns.relplot(x="NominalEs", y="time", row="order", hue="algorithm", data=singleVs)
    plt.savefig("./Timings/Vertices_{}.png".format(vs))

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

