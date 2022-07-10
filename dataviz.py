import numpy as np
import pandas as pd
# import collections as coll
# import colorsys
# from matplotlib import cm
# import matplotlib.pyplot as plt
# import igraph as ig
# import itertools as iter
# import sklearn.metrics as skl
# import os
# import socket
import random
import seaborn as sns

sns.set()

# VYfiles = [
# "./outputTestVY/results2022-06-19 23.01.27.550205.csv",
# "./outputTestSupersetness/results2022-07-04 05.16.59.577390.csv"
# ]

VYfiles = [
'./outputTestSupersetness/results2022-07-04 00.32.07.559015.csv',
'./outputTestSupersetness/results2022-07-04 00.35.03.977908.csv',
'./outputTestSupersetness/results2022-07-04 00.44.41.179319.csv',
'./outputTestSupersetness/results2022-07-04 00.45.42.309361.csv',
'./outputTestSupersetness/results2022-07-04 00.46.44.004117.csv',
'./outputTestSupersetness/results2022-07-04 00.51.14.681204.csv',
'./outputTestSupersetness/results2022-07-04 00.59.46.965891.csv',
'./outputTestSupersetness/results2022-07-04 01.02.07.834240.csv',
'./outputTestSupersetness/results2022-07-04 01.09.53.986792.csv',
'./outputTestSupersetness/results2022-07-04 04.29.08.020363.csv',
'./outputTestSupersetness/results2022-07-04 04.31.19.749463.csv',
'./outputTestSupersetness/results2022-07-04 04.37.17.187231.csv',
'./outputTestSupersetness/results2022-07-04 04.47.47.120816.csv',
'./outputTestSupersetness/results2022-07-04 05.16.59.577390.csv'
]

dfList = []

for id, file in enumerate(VYfiles):
    df = pd.read_csv(VYfiles[1], delimiter=',', header=0, comment="#")
    df["fileID"] = id
    dfList.append(df)


results_wide = pd.concat(dfList)

df1 = results_wide['tangCounts'].str.split('-', expand=True).add_prefix('Order').fillna('')
results_wide = pd.concat([results_wide, df1], axis = 1)

df1 = results_wide['timings'].str.split('-', expand=True).add_prefix('Time').fillna('')
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

results_wide = results_wide.drop(columns=["tangCounts", "timings"] + orders + times)
results_wide = results_wide.dropna(how='all', axis='columns')

print(results_wide.loc[0])

results = pd.wide_to_long(results_wide, ["k"], i=["fileID", "network"], j="order")

results = results.rename(columns = {"k": "time"})


# todo delete these when using real data!!!!!!
shiftVals = []
for i in range(7):
    shiftVals = shiftVals + [random.randint(0,20)]*6
results.time = pd.to_numeric(results.time) + np.random.uniform(-1,1,size=results.shape[0])
results.Vs = pd.to_numeric(results.Vs) + shiftVals
results.Es = pd.to_numeric(results.Es) + shiftVals + np.random.randint(0,10,size=results.shape[0])


sns.relplot(x="Es", y="time", col="order", row="Vs", data=results)