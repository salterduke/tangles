import igraph as ig
import numpy as np
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rc('axes', labelsize=14)
mpl.rc('xtick', labelsize=12)
mpl.rc('ytick', labelsize=12)


class MNIST():
    def __init__(self):
        self.dim = 28 # change if type of image changes, ditto below
        self.numColours = 256

        try:
            from sklearn.datasets import fetch_openml
            self.mnist = fetch_openml('mnist_784', version=1, cache=True)
            self.mnist.target = self.mnist.target.astype(np.int8) # fetch_openml() returns targets as strings
            dummy = 1
        except ImportError:
            from sklearn.datasets import fetch_mldata
            self.mnist = fetch_mldata('MNIST original')

    def fetchSingleImage(self, id = None):
        if id is None:
            id = random.randrange(0, self.mnist.data.shape[0])
        # if id < 0 or id >= self.mnist.data.shape[0]:
            # some error
        imArray = self.mnist.data[id]

        # down-sample to three colours (white grey black)
        # just doing it a stupid way to start with
        self.maxColourCode = 2 # *should* be able to convert to more shades easily, this way
        # todo make smarter

        # image = imArray.reshape(28, 28)
        # plt.imshow(image, cmap=mpl.cm.binary, interpolation="nearest")
        # plt.axis("off")


        imArray = np.array(list(map(lambda x: round(x / self.numColours * self.maxColourCode), imArray)))

        graph = ig.Graph()
        graph.add_vertices(self.dim*self.dim)
        # iterate over *nodes*
        for i in range(self.dim): # rows
            for j in range(self.dim): # columns
                sourceID = i * self.dim + j
                graph.vs[sourceID]["name"] = "{}_{}".format(i,j)
                # add only edges going "right" or "down" to avoid double count
                if i < self.dim - 1:
                    # new edge
                    targetID = (i+1) * self.dim + j
                    # so if adj vs have same/similar colour edge weight is *high*
                    # ie, small cuts between *different* colours
                    edgeWeight = self.maxColourCode - np.abs(imArray[sourceID] - imArray[targetID])
                    graph.add_edge(sourceID, targetID, weight = edgeWeight)
                if j < self.dim - 1:
                    # new edge
                    targetID = (i) * self.dim + (j+1)
                    # so if adj vs have same/similar colour edge weight is *high*
                    # ie, small cuts between *different* colours
                    edgeWeight = self.maxColourCode - np.abs(imArray[sourceID] - imArray[targetID])
                    graph.add_edge(sourceID, targetID, weight = edgeWeight)

        dummy = 1
        graph.delete_edges(weight_eq=0)
        return graph

        # visual_style = {}
        # visual_style["vertex_label"] = imArray
        # color_list = ["white","grey","black"]
        # visual_style["edge_label"] = graph.es["weight"]
        # visual_style["bbox"] = (0, 0, 1000, 1000)
        #
        # visual_style["vertex_color"] = [color_list[colCode] for colCode in imArray]
        #
        # layout = graph.layout_grid()
        # ig.plot(graph, layout=layout, **visual_style)


if __name__ == '__main__':
    print("running test on one image")
    M = MNIST()
    M.fetchSingleImage(0)