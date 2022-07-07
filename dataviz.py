# import numpy as np
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
# import random
import seaborn as sns

sns.set()

VYfiles = [
"./outputTestVY/results2022-06-19 23.01.27.550205.csv",
"./outputTestSupersetness/results2022-07-04 05.16.59.577390.csv"
]

df = pd.read_csv(VYfiles[1], delimiter=',', header=0, comment="#")

df1 = df['tangCounts'].str.split('-', expand=True).add_prefix('Order').fillna('')
df = df.join(df1)

df1 = df['timings'].str.split('-', expand=True).add_prefix('Time').fillna('')
df = df.join(df1)



for col in [col for col in df.columns if "Order" in str(col)]:
    df[col] = df.apply(lambda x: x[col][1], axis=1)

for i in range(1,11):
    df["k{}".format(i)] = ""


orders = [col for col in df.columns if "Order" in col]
times = [col for col in df.columns if "Time" in col]
nOrders = len(orders)



for rid, row in df.iterrows():
    for i in range(nOrders):
        # row["k{}".format(row["Order{}".format(i)])] = row["Time{}".format(i)]
        df.loc[rid, "k{}".format(row["Order{}".format(i)])] = row["Time{}".format(i)]

df = df.drop(columns=["tangCounts", "timings"] + orders + times)

print(df.loc[0])