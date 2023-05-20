import numpy as np
import pandas as pd
# import collections as coll
# import colorsys
# from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, MaxNLocator
import igraph as ig
# import itertools as iter
# import sklearn.metrics as skl
# from sklearn.linear_model import LinearRegression
# import os
# import socket
import random
import seaborn as sns
import warnings
import platform
import math

# warnings.filterwarnings( "ignore", module = "matplotlib\..*" )
warnings.filterwarnings( "ignore", message = ".*required_interactive_framework" )
warnings.filterwarnings( "ignore", message = ".*deprecated" )

sns.set()

class Grapher():
    def __init__(self, dataFolder = None, outputFolder=None):
        if dataFolder is not None:
            self.dataFolder = dataFolder

        if outputFolder is not None:
            self.outputFolder = outputFolder # overwrite if have info
        elif dataFolder is not None:
            self.outputFolder = dataFolder # set default
        else:
            self.outputFolder = "./"


    def plotNetworks(self, inherit = False):

        dataSets = {
            # "YeastA": "../../NetworkData/BioDBs/YeastPPI/YuEtAlGSCompA.csv",
            # "YeastB": "../../NetworkData/BioDBs/YeastPPI/YuEtAlGSCompB.csv",
            # "Bsubtilis": "../../NetworkData/BioDBs/HINTformatted/BacillusSubtilisSubspSubtilisStr168-htb-hq.txt",
            # "Jazz": "../../NetworkData/MediumSize/Jazz.csv",
            # "Copperfield": "../../NetworkData/MediumSize/Copperfield.csv",
            # "Iceland": "../../NetworkData/MediumSize/Iceland.csv",
            "Zebra": "../../NetworkData/MediumSize/Zebra.csv"
            # "Dolphins": "../../NetworkData/MediumSize/Dolphins.csv"
        }

        plainOnly = {
            "Karate": "../../NetworkData/SmallNWs/KarateEdges.csv",
            "Celegans": "../../NetworkData/Celegans/NeuronConnect.csv",
            "Football": "../../NetworkData/MediumSize/Football.csv"
        }

        # for datakey, graphFile in plainOnly.items():
        #     G = ig.Graph.Read_Ncol(graphFile, names=True, directed=False)
        #     G = G.connected_components().giant().simplify()
        #     self.plotSingleClustering(datakey, "plain", G, None)

        for datakey, graphFile in dataSets.items():
            print(datakey)
            G = ig.Graph.Read_Ncol(graphFile, names=True, directed=False)
            G = G.connected_components().giant().simplify()
            if inherit:
                filelabel = "Disj_inherit"
            else:
                filelabel = "Disj_negatives"
            mshipFile = "{}{}_expandedTangleComms{}.csv".format(self.dataFolder, datakey, filelabel)
            mships = pd.read_csv(mshipFile)
            mships = mships.set_index("Unnamed: 0")

            if any(mships.iloc[:, 0] == -1):
                print("{} has unassigned vs at lowest level".format(datakey))
                print(mships.loc[mships.iloc[:, 0] == -1])

            mships, nTangleComms = self.consolidateTangleComms(mships)
            for clusteringName in list(mships.columns):
                if inherit:
                    self.plotSingleClustering(datakey, clusteringName, G, mships[clusteringName], inherit = nTangleComms)
                else:
                    self.plotSingleClustering(datakey, clusteringName, G, mships[clusteringName], inherit = 0)

            self.plotSingleClustering(datakey, "plain", G, None)

    # if diff order tangles equivalent, make sure they've got the same comm ID
    def consolidateTangleComms(self, mships):
        tangOrds = [col for col in mships.columns if col.isnumeric()]
        tangDF = mships[tangOrds]

        replacementTable = dict()

        for columnID, order in enumerate(tangOrds):
            for commID in tangDF[order].unique():
                commVs = tangDF.loc[tangDF[order] == commID]
                for nextColumnID in range(columnID+1, len(tangOrds)):
                    # commVs[tangOrds[nextColumnID]]
                    uniqueComms = commVs[tangOrds[nextColumnID]].unique()
                    if len(uniqueComms) == 1 and uniqueComms[0] != -1:
                        # if there's only one value, the communities are equivalent
                        # (higher order tangle must be \subseteq lower order, so no need to check other way round)
                        higherOrdCommID = uniqueComms[0]
                        if higherOrdCommID != commID:
                            replacementTable[higherOrdCommID] = commID

        chainedReplacements = dict()
        firstGo = True # emulating do while loop
        foundValToChain = False
        # keep checking until we're sure we've got them all.
        while firstGo or foundValToChain:
            foundValToChain = False
            firstGo = False
            for old, new in replacementTable.items():
                if new in replacementTable.keys():
                    foundValToChain = True
                    chainedReplacements[old] = replacementTable[new]
            for old, new in chainedReplacements.items():
                replacementTable[old] = new
            # doing this way so not modifying while iterating

        replacementInstructions = {order: replacementTable for order in tangOrds}
        newMships = mships.replace(replacementInstructions)

        uniqueCommsEachOrd = []
        for order in tangOrds:
            uniqueCommsEachOrd = uniqueCommsEachOrd + list(newMships[order].unique())

        allUnique = np.unique(uniqueCommsEachOrd)
        allUnique.sort()

        allUnique = [comm for comm in allUnique if comm != -1]

        if max(allUnique) != len(allUnique) - 1:
            # ie, there are gaps
            # replace so numbers all contiguous
            for id, val in enumerate(allUnique):
                if id != val:
                    replacementInstructions = {order: {val:id} for order in tangOrds}
                    # print(replacementInstructions)
                    newMships = newMships.replace(replacementInstructions)

        return newMships, len(allUnique)


    def plotSingleClustering(self, dataName, clusteringName, G, clustering, inherit=0, doGrid=False, gridWidth=None):

        G.es["curved"] = 0

        visual_style = {}
        visual_style["vertex_size"] = 14
        visual_style["edge_width"] = 1

        nTangleComms = inherit

        if clustering is None:
            visual_style["vertex_color"] = "Light Sky Blue"
        else:
            if nTangleComms and clusteringName.isnumeric():
                pal = ig.drawing.colors.ClusterColoringPalette(nTangleComms)
            else:
                pal = ig.drawing.colors.ClusterColoringPalette(len(np.unique(clustering)))

            clustering.index = clustering.index.astype("str")

            if clusteringName.isnumeric():
                clusteringName = "ord" + clusteringName
            try:
                # make the colour ids non-negative, and order by frequency so biggest comm first
                commIDS = {val: id for id, val in enumerate(clustering.value_counts().index)}

                if doGrid:
                    visual_style["vertex_frame_width"] = 2
                    visual_style["vertex_frame_color"] = pal.get_many((commIDS[clustering.loc[v["name"]]] for v in G.vs))
                    visual_style["vertex_color"] = G.vs["vertex_color"]
                    visual_style["vertex_label"] = G.vs["name"]
                    visual_style["edge_label"] = list(map(int, G.es["weight"]))

                    if gridWidth is not None:
                        gridHeight = math.ceil(G.vcount() / gridWidth)
                    else:
                        gridHeight = np.sqrt(G.vcount())
                        gridWidth = gridHeight

                    for idx, v in enumerate(G.vs):
                        if clustering.loc[v["name"]] == -1:
                            visual_style["vertex_frame_color"][idx] = "white"
                else:
                    G.vs['color'] = pal.get_many((commIDS[clustering.loc[v["name"]]] for v in G.vs))
                    for v in G.vs:
                        if clustering.loc[v["name"]] == -1:
                            v['color'] = "white"
            except:
                exit("Arrg, something wrong other than keyerror")

        if inherit:
            filelabel = "-inherit"
        else:
            filelabel = "-unassigned"

        if doGrid:
            layout = G.layout_grid(gridWidth)
        else:
            layout = G.layout_kamada_kawai()


        if platform.system() == "Windows":
            filetype = "png"
        else:
            filetype = "pdf"

        # with labels:
        # visual_style["bbox"] = (0, 0, gridWidth * 50, gridHeight * 50)
        # imageFilename = "{}{}-{}{}-labelled.{}".format(self.outputFolder, dataName, clusteringName, filelabel, filetype)
        # ig.plot(G, target=imageFilename, layout=layout, **visual_style)

        # unlabelled
        visual_style["bbox"] = (0, 0, gridWidth * 30, gridHeight * 30)
        visual_style["vertex_label"] = ""
        visual_style["edge_label"] = ""
        imageFilename = "{}{}-{}{}-plain.{}".format(self.outputFolder, dataName, clusteringName, filelabel, filetype)
        ig.plot(G, target=imageFilename, layout=layout, **visual_style)



    def processCDComparisons(self, inherit, methodClass):

        # Only using LFK NMI, so remove qualifier.
        self.metricLongNames = {
            'LFK': "Overlapping NMI",
            'MGH': 'Overlapping NMI - MGH',
            'adjusted_rand': "Adjusted Rand Index",
            'nmi': "Normalised Mutual Information",
            'omega': "Omega Index of Resemblance",
            'rand': "Rand Index",
            'split-join': "Split-Join Distance",
            'vi': "Variation of Information"
        }

        self.tidyShort = {
            'LFK': "LFK",
            'MGH': 'MGH',
            'adjusted_rand': "Adj. Rand",
            'nmi': "NMI",
            'omega': "Omega",
            'rand': "Rand",
            'split-join': "Split-Join",
            'vi': "VI"
        }


        self.methodLongNames = {
            'CPM3': "CPM clique size 3",
            'CPM4': "CPM clique size 4",
            'CPM5': "CPM clique size 5",
            'CPM6': "CPM clique size 6",
            'between': "Edge-betweenness",
            'eigen': "Leading Eigenvector",
            'fastgreedy': "Fast-greedy",
            'infomap': "Infomap",
            'labelprop': "Label Propagation",
            'leiden': "Leiden Algorithm",
            'modularity': "Optimal Modularity",
            'multilevel': "Louvain Modularity",
            'spinglass': "Spinglass",
            'walktrap': "Walktrap",
            'reference': "Ref. Infomap vs Louvain"
        }
        # Note that modularity and multilevel (Louvain) have to be one after another
        # to ensure graphs put them in the same place on the x axis (since only one is included)

        dfList = []
        coverFolder = "../outputdevResults_VY/inheritComparisons/SelectResults/"

        infileLabel = methodClass
        # if inherit:
        #     infileLabel = "{}_inherit".format(methodClass)
        # else:
        #     infileLabel = "{}_negatives".format(methodClass)

        self.refValsDF = pd.read_csv("{}referenceMetrics_{}.csv".format(coverFolder, infileLabel))

        self.objDF = pd.read_csv("{}objectiveFns{}.csv".format(coverFolder, infileLabel))
        self.compDF = pd.read_csv("{}ComparisonValues{}_All.csv".format(coverFolder, infileLabel), delimiter=",", header=0)

        self.compDF["metricLongName"] = self.compDF["metric"].map(self.metricLongNames)
        self.compDF["metricShortName"] = self.compDF["metric"].map(self.tidyShort)
        self.compDF["value"] = self.compDF["value"].round(3)
        self.compDF.drop('Unnamed: 0', axis=1, inplace=True)
        self.compDF["methodLongName"] = self.compDF["method"].map(self.methodLongNames)

        self.refValsDF.drop('Unnamed: 0', axis=1, inplace=True)
        self.refValsDF["method"] = "reference"
        self.refValsDF["metricLongName"] = self.refValsDF["metric"].map(self.metricLongNames)
        self.refValsDF["metricShortName"] = self.refValsDF["metric"].map(self.tidyShort)
        self.refValsDF["methodLongName"] = self.refValsDF["method"].map(self.methodLongNames)
        self.refValsDF.drop('method1', axis=1, inplace=True)
        self.refValsDF.drop('method2', axis=1, inplace=True)

        # todo this for average CD as well
        # self.plotSimilarityAgainstobjFunc(inherit)

        for dataName in np.unique(self.compDF["dataName"]):
            self.processSingleNetwork(dataName, inherit)

    def makeNetworkOrderTable(self, df, disjoint = True, fileLabel = "table"):
        outputTableFolder = "../outputdevResults_VY/Tables/"
        if df.size == 0:
            return

        dataName = df["dataName"].iloc[0]
        order = df["order"].iloc[0]

        try:
            df = df.pivot(columns="metricShortName", index="methodLongName", values="value")
        except:
            dummy = 1
            exit("File not pivoting right")

        sortOrder = [methodLong for methodLong in self.methodLongNames.values()
                 if methodLong in np.unique(df.index)]
        df = df.loc[sortOrder]

        tableFilename = "{}{}-ord{}-{}.csv".format(outputTableFolder, dataName, order, fileLabel)
        df.index.name = "Community Detection Method"
        df.to_csv(tableFilename)

    def plotSimilarityAgainstobjFunc(self, inherit):

        self.objDF.drop('Unnamed: 0', axis=1, inplace=True)
        self.objDF = self.objDF.rename(columns={"objFunction":"method"})

        specificComps = self.compDF.loc[self.compDF.method.isin(("multilevel", "infomap"))].reset_index()
        specificComps = specificComps.set_index(["dataName", "method"])
        self.objDF = self.objDF.set_index(["dataName", "method"])

        joinedDF = specificComps.join(self.objDF, lsuffix="_metric", rsuffix="_obj")
        joinedDF.xs("infomap", level="method")

        joinedDF = joinedDF.loc[joinedDF.metric.isin(("nmi", "adjusted_rand"))]

        if inherit:
            fileLabel = "inherit"
        else:
            fileLabel = "negatives"

        for method in joinedDF.index.unique(level='method'):
            xlabel = {
                "multilevel": "Modularity",
                "infomap": "Map equation"
            }

            truncData = joinedDF.xs(method, level="method")
            # truncData = truncData.loc[truncData["order"] < 10]
            # testing:
            g = sns.relplot(x="order", y="value_metric", col="metricLongName", hue="dataName", palette="bright",
                            data=truncData)
            g.set_axis_labels("Tangle order", "Similarity metric value")
            g.set_titles(col_template='{col_name}')
            g._legend.set_title("Network")
            plt.savefig("../outputdevResults_VY/Visualisations_v2/Similarity-vs-order-{}-{}.pdf".format(method, fileLabel))


            # g = sns.relplot(x="value_obj", y="value_metric", col="metricLongName", hue="dataName", palette="bright",
            #                 data=joinedDF.xs(method, level="method"))
            g = sns.relplot(x="value_obj", y="value_metric", col="metricLongName", hue="dataName", palette="bright",
                            size="order", sizes=(20, 200), data=joinedDF.xs(method, level="method"))
            g.set_axis_labels("{} value".format(xlabel[method]), "Similarity metric value")
            g.set_titles(col_template='{col_name}')
            g._legend.set_title("Network")
            for textItem in g.legend.texts:
                if textItem._text == "dataName":
                    textItem._text = "Network"
                elif textItem._text == "order":
                    textItem._text = "Order"
                    break

            # plt.show(block=True)
            plt.savefig("../outputdevResults_VY/Visualisations_v2/Similarity-vs-objective-{}-{}.pdf".format(method, fileLabel))


    def plotSelectedMetrics(self, df, disjoint = True, fileLabel = None, onlyOrders=None):
        outputGraphsFolder = "../outputdevResults_VY/Visualisations_v2/"
        if df.size == 0:
            return

        # removing because igraph leiden fn doesn't seem to work with default params
        df = df.loc[df["method"] != "leiden"]

        # removing because don't need all modularity methods
        if any(df.method == "modularity"):
            deleteModMethods = {"eigen", "multilevel", "fastgreedy"}
        else:
            deleteModMethods = {"eigen", "fastgreedy"}
        df = df.loc[~df["method"].isin(deleteModMethods)]

        if disjoint:
            metrics = ("adjusted_rand", "nmi")
        else:
            metrics = np.unique(df["metric"])

        df = df.loc[df.metric.isin(metrics)]
        if onlyOrders:
            df = df.loc[df.order.isin(onlyOrders)]

        fg = sns.catplot(data=df, kind="bar", x="methodLongName", y="value", col="order", row="metricLongName",
                    order=[methodLong for methodLong in self.methodLongNames.values()
                        if methodLong in np.unique(df["methodLongName"])]
                    )

        # fg.set_xticklabels(rotation=-90, horizontalalignment='left', rotation_mode='anchor')
        fg.tick_params("x", direction="out", bottom=True, labelrotation=-90)
        fg.set(xlabel = "Community Detection Method")
        fg.set_titles(col_template='{col_name}', row_template='{row_name}')
        for ax in fg.axes.flat:
            title = ax.get_title()
            metric = title.split("|")[0].strip()
            order = title.split("|")[1].strip()
            ax.set_title("Order {} tangles".format(order))
            ax.set(ylabel=metric)

        fg.fig.subplots_adjust(wspace=.1)

        # fig.tight_layout(pad=1)
        # fig.align_ylabels()

        if fileLabel is not None:
            dataName = df["dataName"].iloc[0]
            outputFileName = "{}{}-{}.pdf".format(outputGraphsFolder, dataName, fileLabel)
            fg.savefig(outputFileName)
            plt.close(fg.fig)
        else:
            # actually this won't work, but eh.
            fg.show()

    def plotNetworkOrder(self, df, numMetrics = 5, disjoint = True, fileLabel = None, removeExtraneous = False):
        outputGraphsFolder = "../outputdevResults_VY/Visualisations/"
        if df.size == 0:
            return

        if removeExtraneous:
            # removing because igraph leiden fn doesn't seem to work with default params
            df = df.loc[df["method"] != "leiden"]

            # removing because don't need all modularity methods
            if "modularity" in df["method"]:
                deleteModMethods = {"eigen", "multilevel", "fastgreedy"}
            else:
                deleteModMethods = {"eigen", "fastgreedy"}
            df = df.loc[~df["method"].isin(deleteModMethods)]

        if numMetrics == 3:
            if disjoint:
                metrics = ("adjusted_rand", "split-join", "nmi")
            else:
                metrics = np.unique(df["metric"])

            fig, ax = plt.subplots(3,1, figsize=(6.4,9))
            for id, metric in enumerate(metrics):
                sns.barplot(data = df.loc[df["metric"] == metric],
                            x="methodLongName", y="value", ax=ax[id],
                            order=[methodLong for methodLong in self.methodLongNames.values()
                                   if methodLong in np.unique(df["methodLongName"])]
                            )
                ax[id].set_box_aspect(6/len(ax[id].get_xticklabels()))
                if id == 2:
                    ax[id].tick_params("x", direction="out", bottom=True)
                    # ax[id].set_xticklabels(ax[id].get_xticklabels(), rotation=-90, horizontalalignment='right', rotation_mode='anchor')
                    ax[id].set_xticklabels(ax[id].get_xticklabels(), rotation=-90)
                    ax[id].set(ylabel=self.metricLongNames[metric], xlabel = "Community Detection Method")
                else:
                    ax[id].tick_params("x", direction="out", bottom=False)
                    ax[id].set_xticklabels([])
                    ax[id].set(ylabel=self.metricLongNames[metric], xlabel="")

        elif numMetrics == 5:
            metrics = np.unique(df["metric"])
            fig, allAxes = plt.subplots(3,2, figsize=(7.5,9))
            allAxes.flat[1].set_visible(False)
            for metID, metric in enumerate(metrics):
                if metID == 1 or metID == 3:
                    id = metID + 2
                    # so that the gap is on the top of the second column, not the bottom
                else:
                    id = metID
                ax = allAxes.flat[id]
                sns.barplot(data = df.loc[df["metric"] == metric],
                            x="methodLongName", y="value", ax=ax,
                            order=[methodLong for methodLong in self.methodLongNames.values()
                                   if methodLong in np.unique(df["methodLongName"])]
                            )
                if id == 4 or id == 5: #ie, the bottom ones
                    ax.tick_params("x", direction="out", bottom=True)
                    # ax.set_xticklabels(ax.get_xticklabels(), rotation=-90, horizontalalignment='right', rotation_mode='anchor')
                    ax.set_xticklabels(ax.get_xticklabels(), rotation=-90)
                    ax.set(ylabel=self.metricLongNames[metric], xlabel = "Community Detection Method")
                else:
                    ax.tick_params("x", direction="out", bottom=False)
                    ax.set_xticklabels([])
                    ax.set(ylabel=self.metricLongNames[metric], xlabel="")

        fig.tight_layout(pad=1)
        fig.align_ylabels()

        if fileLabel is not None:
            dataName = df["dataName"].iloc[0]
            order = df["order"].iloc[0]
            outputFileName = "{}{}-ord{}-{}.pdf".format(outputGraphsFolder, dataName, order, fileLabel)
            fig.savefig(outputFileName)
            plt.close()
        else:
            fig.show()
            plt.close()

    def plotMetricsLinegraph(self, df, dataName, disjoint=True, inherit=True):

        if disjoint:
            df = df.loc[df.metric.isin(("nmi", "adjusted_rand"))]
            fileLabel = dataName + "-disjoint"
            if any(df.method == "modularity"):
                deleteModMethods = {"eigen", "multilevel", "fastgreedy", "leiden"}
            else:
                deleteModMethods = {"eigen", "fastgreedy", "leiden"}
            df = df.loc[~df["method"].isin(deleteModMethods)]
        else:
            df = df.loc[df.metric.isin(("LFK", "omega"))]
            fileLabel = dataName + "-overlap"

        if inherit:
            fileLabel = fileLabel + "-inherit"
        else:
            fileLabel = fileLabel + "-negatives"

        if disjoint:
            # p = sns.relplot(data=singleVs, x="Es", y="delay", col="sepOrder", col_wrap=2, hue="algorithm",
            #                 palette="bright", errorbar=('ci', 95), kind='line')
            g = sns.relplot(x="order", y="value", kind="line", col="metricLongName", errorbar=('ci', 95),
                            hue="methodLongName", marker="o", palette="bright", data=df)
        else:
            g = sns.relplot(x="order", y="value", kind="line", col="metricLongName",
                            hue="methodLongName", marker="o", palette="bright", data=df)
        g.set_axis_labels("Tangle order", "Similarity metric value")
        g.set_titles(col_template='{col_name}')
        g._legend.set_title("CD Method")
        ax = plt.gca()
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.savefig("../outputdevResults_VY/Visualisations_CI/{}-SimilarityLinegraph.pdf".format(fileLabel))

    def processSingleNetwork_Old(self):
        # this is the order currently used in the results discussion in latex
        selectedOrders = {
            "YeastA": (2,3,4,5),
            "YeastB": (3,4),
            "Bsubtilis": (3,),
            "Jazz": (4,),
            "Copperfield": (5,),
            "Dolphins": (3,7,8,9),
            "Zebra": (4,5,6,14),
            "Iceland": (3,4,5,6)
        }

        singleDF = self.compDF.loc[self.compDF["dataName"] == dataName]

        singleRefs = self.refValsDF.loc[self.refValsDF["dataName"] == dataName]
        meanRefs = singleRefs.groupby("metric")["value"].mean()

        disjointMethodsAll = singleDF.loc[~singleDF["method"].str.contains("CPM")]
        overlapMethodsAll = singleDF.loc[singleDF["method"].str.contains("CPM")]
        for order in np.unique(singleDF["order"]):
            for metric, val in meanRefs.items():
                rowdict =  {'order':order,
                            'method':"Reference",
                            'metric':metric,
                            'value':val,
                            'dataName':dataName,
                            'methodLongName':self.methodLongNames["Reference"],
                            'metricLongName':self.metricLongNames[metric],
                            'metricShortName':self.tidyShort[metric]
                            }
                newRow = pd.DataFrame([rowdict])
                disjointMethodsAll = pd.concat([disjointMethodsAll, newRow], ignore_index=True)

        onlyOrders = selectedOrders[dataName]
        print(dataName)
        print(onlyOrders)
        if len(onlyOrders) <= 2:
            self.plotSelectedMetrics(disjointMethodsAll, disjoint=True, fileLabel="selected", onlyOrders = onlyOrders)
            self.plotSelectedMetrics(overlapMethodsAll, disjoint=False, fileLabel="selected", onlyOrders = onlyOrders)
        elif len(onlyOrders) <= 4:
            ords1 = onlyOrders[0:2]
            ords2 = onlyOrders[2:]
            self.plotSelectedMetrics(disjointMethodsAll, disjoint=True, fileLabel="selected", onlyOrders = ords1)
            self.plotSelectedMetrics(disjointMethodsAll, disjoint=True, fileLabel="selectedPt2", onlyOrders = ords2)
            self.plotSelectedMetrics(overlapMethodsAll, disjoint=False, fileLabel="selected", onlyOrders = ords1)
            self.plotSelectedMetrics(overlapMethodsAll, disjoint=False, fileLabel="selectedPt2", onlyOrders = ords2)

        for order in np.unique(singleDF["order"]):
            disjointMethods = singleDF.loc[(~singleDF["method"].str.contains("CPM")) & (singleDF["order"] == order)]
            overlapMethods = singleDF.loc[(singleDF["method"].str.contains("CPM")) & (singleDF["order"] == order)]
            # yeah, I know this is redundant, but the dfs just above are based on singleDF
            #
            for metric, val in meanRefs.items():
                rowdict =  {'order':order,
                            'method':"Reference",
                            'metric':metric,
                            'value':val,
                            'dataName':dataName,
                            'methodLongName':self.methodLongNames["Reference"],
                            'metricLongName':self.metricLongNames[metric],
                            'metricShortName':self.tidyShort[metric]
                            }
                newRow = pd.DataFrame([rowdict])
                disjointMethods = pd.concat([disjointMethods, newRow], ignore_index=True)

            self.plotNetworkOrder(disjointMethods, numMetrics=3, disjoint=True, fileLabel="disj-3", removeExtraneous = True)
            self.plotNetworkOrder(disjointMethods, numMetrics=5, disjoint=True, fileLabel="disj-5")
            self.plotNetworkOrder(overlapMethods, numMetrics=3, disjoint=False, fileLabel="overlap")

            self.makeNetworkOrderTable(disjointMethods, fileLabel="disj")
            self.makeNetworkOrderTable(overlapMethods, fileLabel="overlap")

    def processSingleNetwork(self, dataName, inherit):

        singleDF = self.compDF.loc[self.compDF["dataName"] == dataName]

        singleRefs = self.refValsDF.loc[self.refValsDF["dataName"] == dataName]

        disjointMethodsAll = singleDF.loc[~singleDF["method"].str.contains("CPM")]
        overlapMethodsAll = singleDF.loc[singleDF["method"].str.contains("CPM")]
        for order in singleDF["order"].unique():
            for metric in singleRefs["metric"].unique():
                # 10 reps hard coded I know, but I don't care at this stage.
                # pick 10 ref metricss to add to DF
                refMetrics = singleRefs.loc[singleRefs["metric"] == metric]
                for rep in range(10):

                    rowdict =  {'order':order,
                                'method':"reference",
                                'metric':metric,
                                'value':refMetrics.iloc[rep]["value"],
                                'dataName':dataName,
                                'methodLongName':self.methodLongNames["reference"],
                                'metricLongName':self.metricLongNames[metric],
                                'metricShortName':self.tidyShort[metric]
                                }
                    newRow = pd.DataFrame([rowdict])
                    disjointMethodsAll = pd.concat([disjointMethodsAll, newRow], ignore_index=True)

        print(dataName)

        if disjointMethodsAll.size > 0:
            self.plotMetricsLinegraph(disjointMethodsAll, dataName=dataName, disjoint=True, inherit=inherit)
        if overlapMethodsAll.size > 0:
            self.plotMetricsLinegraph(overlapMethodsAll, dataName=dataName, disjoint=False, inherit=inherit)

        # todo fix these
        # self.makeNetworkOrderTable(disjointMethods, fileLabel="disj")
        # self.makeNetworkOrderTable(overlapMethods, fileLabel="overlap")

# ---------------------------------------------------------------------------------------
    def readAlgTimingData(self, timingFiles, algName):

        dfList = []

        for id, file in enumerate(timingFiles):
            # print("Reading ", file)
            df = pd.read_csv(file, delimiter=',', header=0, comment="#")
            # print(np.unique(df.Vs))
            df["fileID"] = id
            df["algorithm"] = algName

            dfList.append(df)

        results_wide = pd.concat(dfList, ignore_index=True)
        # np.unique(results_wide.groupby("network").count().Vs)
        # This was just to get an idea of what timing data I actually *have*
        # actual = results_wide.groupby("network")[("Vs", "Es")].max()
        # sns.scatterplot(data=actual, x="Vs", y="Es").set(title=algName)
        # plt.savefig("Es-vs-Vs-{}.pdf".format(algName))

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
            exit()

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
    def plotTimingResults(self, results):

        #print("Results-------------")
        #print(results.head())


        timeSummary = results.groupby(["NominalEs", "NominalVs", "order", "algorithm"])["delay"].mean().reset_index()

        timeSummary["sepOrder"] = timeSummary["order"] - 1

        # eplots = sns.relplot(x="NominalEs", y="delay", col="order", row="NominalVs", hue="algorithm", data=timeSummary)
        # eplots.set_axis_labels("Number of Edges (m)", "Delay per cut (seconds)")
        # eplots.set_titles(row_template='Vertices n={row_name}', col_template='Order {col_name} Separations')
        # eplots._legend.set_title("Algorithm")
        # plt.show(block=True)
        # # plt.savefig("./Timings/All-against-edges.pdf")

        # nplots = sns.relplot(x="NominalVs", y="delay", row="sepOrder", hue="NominalEs", col="algorithm", palette="copper_r", data=timeSummary)
        # nplots.set_axis_labels("Number of vertices (n)", "Delay per cut (seconds)")
        # nplots.set_titles(row_template='Order {row_name} Separations', col_template='Algorithm: {col_name}')
        # nplots._legend.set_title("Edges (m)")
        # # plt.show(block=True)
        # plt.savefig("./Timings/All-against-verts.pdf")
        #
        # for ord in np.unique(timeSummary["sepOrder"]):
        #     singleOrd = timeSummary.loc[timeSummary["sepOrder"] == ord]
        #     nplots = sns.relplot(x="NominalVs", y="delay", hue="NominalEs", col="algorithm", palette="copper_r", data=singleOrd)
        #     nplots.set_axis_labels("Number of vertices (n)", "Delay per cut (seconds)")
        #     nplots.set_titles(col_template='Algorithm: {col_name}')
        #     # plt.show(block=True)
        #     plt.savefig("./Timings/Ord{}-against-verts.pdf".format(ord))

        for vs in (20,50,100,150,200):
            singleVs = results.loc[results.NominalVs == vs].groupby(["NominalEs", "NominalVs", "order", "algorithm"])[
                "delay"].mean().reset_index()
            # singleVs = singleVs[singleVs["order"] < 5]
            singleVs["sepOrder"] = singleVs["order"] - 1 
            eplots = sns.relplot(x="NominalEs", y="delay", col="sepOrder", col_wrap=2, hue="algorithm", data=singleVs)
            eplots.set_axis_labels("Number of Edges (m)", "Delay per cut (seconds)")
            eplots.set_titles(col_template='Order {col_name} Separations')
            eplots._legend.set_title("Algorithm")
            # plt.show(block=True)
            plt.savefig("./Timings/Vs-{}-against-edges.pdf".format(vs))

    # ------------------------------------------------------------

    def readSepCounts(self, countFiles, algName):
        dfList = []

        for id, file in enumerate(countFiles):
            # print("Reading ", file)
            df = pd.read_csv(file, delimiter=',', header=0, comment="#")
            df["fileID"] = id
            df["algorithm"] = algName
            dfList.append(df)

        counts_wide = pd.concat(dfList, ignore_index=True)

        df1 = counts_wide['sepCounts'].str.split('-', expand=True).add_prefix('sepNumber').fillna('')
        counts_wide = pd.concat([counts_wide, df1], axis = 1)

        df1 = counts_wide['tangCounts'].str.split('-', expand=True).add_prefix('Order').fillna('')
        counts_wide = pd.concat([counts_wide, df1], axis = 1)

        counts_wide.drop('Order4', axis=1, inplace=True, errors="ignore")
        counts_wide.drop('Order5', axis=1, inplace=True, errors="ignore")
        counts_wide.drop('Order6', axis=1, inplace=True, errors="ignore")

        for col in [col for col in counts_wide.columns if "Order" in str(col)]:
            counts_wide[col] = counts_wide.apply(lambda x: x[col][1], axis=1)

        orders = [col for col in counts_wide.columns if "Order" in col]
        counts = [col for col in counts_wide.columns if "sepNumber" in col]
        nOrders = len(orders)

        for i in range(1,10):
            counts_wide["seps_at_ord{}".format(i)] = np.nan

        for rid, row in counts_wide.iterrows():
            for i in range(nOrders):
                counts_wide.loc[rid, "seps_at_ord{}".format(row["Order{}".format(i)])] = row["sepNumber{}".format(i)]

        counts_wide = counts_wide.drop(columns=["tangCounts", "timings", "Es", "Vs", "Unnamed: 0", "secs", "sepCounts"] + orders + counts)
        counts_wide = counts_wide.dropna(how='all', axis='columns')

        results = pd.wide_to_long(counts_wide, ["seps_at_ord"], i=["fileID", "network"], j="order")
        results = results.rename(columns = {"seps_at_ord": "numSeps"})
        results = results.reset_index()
        results.numSeps = pd.to_numeric(results.numSeps)

        return results


    def combineSepcounts(self, timeDF, sepDF):
        # first check that sepnums in VY and YWS match!
        sepNumsVY = sepDF.loc[sepDF.algorithm == "VY"]
        sepNumsYWS = sepDF.loc[sepDF.algorithm == "YWS"]

        YWStoadd = []
        for id, row in sepNumsVY.iterrows():
            rowYWS = sepNumsYWS.loc[
                (sepNumsYWS.network == row.network) &
                (sepNumsYWS.order == row.order)
            ]
            try:
                rowYWS = pd.Series(rowYWS.iloc[0])
            except:
                # kind of dodgy, but it should work, since the sepNums *should* be the same
                # can remove when YWS 200 finished running
                rowYWS = row
                rowYWS["algorithm"] = "YWS"
                YWStoadd.append(rowYWS)

            if row.numSeps != rowYWS.numSeps:
                print("ARRGH, numSeps doesn't match")
                print(row)
                print(rowYWS)
                exit()

        try:
            df = pd.concat(YWStoadd, axis=1).T
            sepDF = pd.concat([sepDF, df])
        except:
            dummy = 1

        sepDF = sepDF.set_index(["network", "order", "algorithm"])
        timeDF = timeDF.set_index(["network", "order", "algorithm"])

        sepDF = sepDF.drop(columns = ["fileID"])

        timeDF = timeDF.join(sepDF)
        timeDF["delay"] = timeDF.time / timeDF.numSeps

        return timeDF

    def processTimingData(self):
        VYfiles = [
            "../outputTestVY/results2022-11-23 20.49.11.331334.csv",
            "../outputTestVY/results2022-11-25 00.56.39.342710.csv",
            "../outputTestVY/results2022-11-25 13.58.17.443580.csv",
            "../outputTestVY/results2022-11-26 22.46.30.592230.csv",
            "../outputTestVY/results2022-12-18 22.15.15.806828.csv",
            "../outputTestVY/results2022-12-19 16.52.11.097076.csv",
            "../outputTestVY/results2022-12-26 23.03.21.097428.csv",
            "../outputTestVY/results2022-12-27 18.20.50.612801.csv",
            "../outputTestVY/results2023-01-03 17.54.11.479864.csv",
            "../outputTestVY/results2023-01-04 18.50.10.710700.csv",
            "../outputTestVY/results2023-01-24 00.38.06.570735.csv",
            "../outputTestVY/results2023-01-24 06.17.45.901736.csv",
            "../outputTestVY/results2023-01-24 01.55.55.688959.csv",
            "../outputTestVY/results2023-01-24 08.53.11.574409.csv",
            "../outputTestVY/results2023-01-24 04.07.53.285378.csv"
        ]
        # last five are 20 - 100

        VYfiles_3ordsonly = [
            "../outputTestVY/results2022-12-10 09.10.23.870668.csv",
            "../outputTestVY/results2022-12-10 12.38.50.019390.csv",
            "../outputTestVY/results2022-12-11 10.43.25.282799.csv",
            "../outputTestVY/results2022-12-11 14.12.00.441343.csv",
            "../outputTestVY/results2022-12-12 02.14.09.467554.csv",
            "../outputTestVY/results2022-12-12 06.07.45.566873.csv",
            "../outputTestVY/results2022-12-12 12.19.58.915419.csv",
            "../outputTestVY/results2022-12-12 15.48.30.065434.csv",
            "../outputTestVY/results2022-12-13 15.02.59.722227.csv",
            "../outputTestVY/results2022-12-13 19.05.02.100191.csv"
        ]

        YWSfiles = [
            "../outputTestYWS/results2022-11-24 04.13.04.126166.csv",
            "../outputTestYWS/results2022-11-25 09.06.28.692849.csv",
            "../outputTestYWS/results2022-11-28 17.05.09.511741.csv",
            "../outputTestYWS/results2022-11-29 19.19.01.456783.csv",
            "../outputTestYWS/results2022-12-19 05.22.53.195299.csv",
            "../outputTestYWS/results2022-12-22 12.04.04.355848.csv",
            "../outputTestYWS/results2022-12-27 06.33.52.986415.csv",
            "../outputTestYWS/results2022-12-30 10.29.12.984751.csv",
            "../outputTestYWS/results2023-01-04 01.38.36.360110.csv",
            "../outputTestYWS/results2023-01-08 20.27.45.931906.csv",
            "../outputTestYWS/results2023-01-24 01.05.40.484337.csv",
            "../outputTestYWS/results2023-01-24 07.14.30.492746.csv",
            "../outputTestYWS/results2023-01-24 02.40.56.709631.csv",
            "../outputTestYWS/results2023-01-24 09.48.47.964709.csv",
            "../outputTestYWS/results2023-01-24 04.42.06.618790.csv"
        ]
        # last five are 20 - 100

        YWSfiles_3ordsonly = [
            "../outputTestYWS/results2022-12-10 09.17.26.935472.csv",
            "../outputTestYWS/results2022-12-10 13.08.02.281325.csv",
            "../outputTestYWS/results2022-12-11 10.50.25.246702.csv",
            "../outputTestYWS/results2022-12-11 14.41.15.343567.csv",
            "../outputTestYWS/results2022-12-12 02.21.42.212079.csv",
            "../outputTestYWS/results2022-12-12 06.38.50.450858.csv",
            "../outputTestYWS/results2022-12-12 12.27.02.689464.csv",
            "../outputTestYWS/results2022-12-12 16.17.42.518543.csv",
            "../outputTestYWS/results2022-12-13 15.10.50.769288.csv",
            "../outputTestYWS/results2022-12-13 19.37.09.424597.csv"
        ]

        # VYfiles = [
        #     "../outputTestVY/timesexampleVY.csv"
        # ]
        #
        # YWSfiles = [
        #     "../outputTestYWS/timesexampleYWS.csv"
        # ]

        sepNumbersVY = [
            "../outputTestVY/results2023-01-24 21.31.56.721371.csv",
            "../outputTestVY/results2023-01-26 05.33.37.125949.csv",
            "../outputTestVY/results2023-04-13 01.54.13.023727.csv"
        ]
        sepNumbersYWS = [
            "../outputTestYWS/results2023-01-25 10.27.19.290499.csv",
            "../outputTestYWS/results2023-04-13 02.40.57.263470.csv"
        ]

        # sepNumbersVY = [
        #     "../outputTestVY/sepcountsVY.csv"
        # ]
        #
        # sepNumbersYWS = [
        #     "../outputTestYWS/sepcountsYWS.csv"
        # ]

        fileLists = [VYfiles, YWSfiles]
        sepcountLists = [sepNumbersVY, sepNumbersYWS]
        algNames = ["VY", "YWS"]
        timeDFs = []
        sepcountDFs = []

        for id in (0,1):
            timeDFs.append(self.readAlgTimingData(fileLists[id], algNames[id]))
            sepcountDFs.append(self.readSepCounts(sepcountLists[id], algNames[id]))

        timings = pd.concat(timeDFs)
        sepCounts = pd.concat(sepcountDFs)


        results = self.combineSepcounts(timings, sepCounts)

        self.plotRestrictedTimings(results)
        # self.plotTimingResults(results)

    def plotRestrictedTimings(self, results):

        bins = list(range(10, 210, 10))
        labels = list(range(20, 210, 10))
        results["Vcategory"] = pd.cut(results["Vs"], bins=bins, labels=labels)
        results = results.reset_index()
        results["sepOrder"] = results["order"] - 1

        for vs in (20,50,100,150,200):
            singleVs = results.loc[results["Vcategory"] == vs]
            p = sns.relplot(data=singleVs, x="Es", y="delay", col="sepOrder", col_wrap=2, hue="algorithm",
                            palette="bright", errorbar=('ci', 95), kind='line')
            p.set_axis_labels("Number of Edges (m)", "Delay per cut (seconds)")
            p.set_titles(col_template='Order {col_name} Separations')
            p._legend.set_title("Algorithm")
            # plt.show(block=True)
            plt.savefig("../TimingPlots/Vs-{}-against-edges.pdf".format(vs))

            # if want bars instead:
            # p = sns.relplot(data=singleVs, x="Es", y="delay", col="sepOrder", col_wrap=2, hue="algorithm",
            #     err_style="bars", errorbar=('ci', 95), kind='line', err_kws={"capsize": 3})
            # plt.show()




if __name__ == '__main__':
    # processTimingData()
    grapher = Grapher(dataFolder="../outputdevResults_VY/inheritComparisons/")

    # for inherit in (True, False):
    for inherit in (True,):
        # for methodClass in ("CPM", "Disj"):
        for methodClass in ("Disj",):
            grapher.processCDComparisons(inherit, methodClass)

    # grapher.plotNetworks(inherit=False)
    # grapher.processTimingData()

