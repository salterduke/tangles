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
"./Timings/VY_results2022-07-06 22.59.30.385284.csv"
]

# VYfiles = [
# './outputTestSupersetness/results2022-07-04 00.32.07.559015.csv',
# './outputTestSupersetness/results2022-07-04 00.35.03.977908.csv',
# './outputTestSupersetness/results2022-07-04 00.44.41.179319.csv',
# './outputTestSupersetness/results2022-07-04 00.45.42.309361.csv',
# './outputTestSupersetness/results2022-07-04 00.46.44.004117.csv',
# './outputTestSupersetness/results2022-07-04 00.51.14.681204.csv',
# './outputTestSupersetness/results2022-07-04 00.59.46.965891.csv',
# './outputTestSupersetness/results2022-07-04 01.02.07.834240.csv',
# './outputTestSupersetness/results2022-07-04 01.09.53.986792.csv',
# './outputTestSupersetness/results2022-07-04 04.29.08.020363.csv',
# './outputTestSupersetness/results2022-07-04 04.31.19.749463.csv',
# './outputTestSupersetness/results2022-07-04 04.37.17.187231.csv',
# './outputTestSupersetness/results2022-07-04 04.47.47.120816.csv',
# './outputTestSupersetness/results2022-07-04 05.16.59.577390.csv'
# ]

dfList = []

for id, file in enumerate(VYfiles):
    df = pd.read_csv(file, delimiter=',', header=0, comment="#")
    df["fileID"] = id
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

# shiftVals = []
# for i in range(7):
#     shiftVals = shiftVals + [random.randint(0,20)]*6
# results.time = pd.to_numeric(results.time) + np.random.uniform(-1,1,size=results.shape[0])
# results.Vs = pd.to_numeric(results.Vs) + shiftVals
# results.Es = pd.to_numeric(results.Es) + shiftVals + np.random.randint(0,10,size=results.shape[0])

pd.unique(results.NominalVs)


shortres = results.loc[results.NominalVs.isin([20, 50, 90])]
shortres_noOutlier = shortres.loc[shortres.time < 1500]

shortres.shape

sns.relplot(x="NominalEs", y="time", col="order", row="NominalVs", data=shortres_noOutlier)
plt.show()