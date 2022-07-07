import numpy as np
import pandas as pd
import collections as coll
import colorsys
from matplotlib import cm
import matplotlib.pyplot as plt
import igraph as ig
import itertools as iter
import sklearn.metrics as skl
import os
import socket
import random
import seaborn as sns

sns.set()

VYfiles = [
"./outputTestVY/results2022-06-19 23.01.27.550205.csv"
]

pd.read_csv(VYfiles[0], delimiter=',', header=0, comment="#")

df1 = df['tangCounts'].str.split('-', expand=True).add_prefix('Order').fillna('')
df = df.join(df1)

for col in [col for col in df.columns if "Order" in str(col)]:
    df[col] = df.apply(lambda x: x[col][1], axis=1)

for i in range(1,11):
    df[i] = ""