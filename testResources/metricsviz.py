import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()



datasets = [
    {"short": "Karate", "title": "Karate Club Members", "fname":"./outputDevYWS/Karate-CDcomparisons.csv"},
    {"short": "YeastA", "title": "Yeast PPI, Component A", "fname":"./outputDevYWS/YeastGSCompA-CDcomparisons.csv"},
    {"short": "Copperfield", "title": "Copperfield", "fname":"./outputDevYWS/Copperfield-CDcomparisons.csv"},
    {"short": "YeastB", "title": "Yeast PPI, Component B", "fname": "./outputDevYWS/YeastGSCompB_core-CDcomparisons.csv"}
]

def methodComparisons(dataset):
    compData = pd.read_csv(dataset["fname"], delimiter=',', header=0, comment="#")
    compData = compData[~compData.method.str.contains("CPM")]

    for metricVal in pd.unique(compData.metric):
        metricData = compData[compData.metric == metricVal]

        chart = sns.catplot(x="method", y="value", col="order", kind="bar", data=metricData)

        for axes in chart.axes.flat:
            _ = axes.set_xticklabels(axes.get_xticklabels(), rotation=-45)

        plt.tight_layout()
        plt.savefig("./testResources/metricGraphs/{}_{}.pdf".format(dataset["short"], metricVal))
        # plt.show()

        # looks shithouse
        # grid = sns.FacetGrid(metricData, col="order",palette="Set3")
        # grid.map(sns.barplot,'method','value')
        # plt.savefig("./testResources/test1_{}.png".format(metricVal))



# --------------------------------
for entry in datasets:
    methodComparisons(entry)


