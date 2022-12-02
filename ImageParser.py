import igraph as ig
import numpy as np
import random
import matplotlib as mpl
import PIL

import matplotlib.pyplot as plt
mpl.rc('axes', labelsize=14)
mpl.rc('xtick', labelsize=12)
mpl.rc('ytick', labelsize=12)


class ImageParser():
    def __init__(self):
        self.numColours = 256   # change if type of image changes
        self.iconList = ["headphones",
                         "printer",
                         "binderpage",
                         "lightbulb",
                         "flame",
                         "folder"]

        self.mnist = None

    def initMNIST(self):
        try:
            from sklearn.datasets import fetch_openml
            self.mnist = fetch_openml('mnist_784', version=1, cache=True)
            self.mnist.target = self.mnist.target.astype(np.int8) # fetch_openml() returns targets as strings
            dummy = 1
        except ImportError:
            from sklearn.datasets import fetch_mldata
            self.mnist = fetch_mldata('MNIST original')


    def fetchMNISTasARRAY(self, id=None):
        if self.mnist is None:
            self.initMNIST()    # ie, only do this if actually doing MNIST, as it takes a little while

        self.dim = 28 # change if type of image changes, ditto below
        if id is None:
            id = random.randrange(0, self.mnist.data.shape[0])
        # if id < 0 or id >= self.mnist.data.shape[0]:
            # some error
        imArray = self.mnist.data[id]

        # MNIST stored as 0 = white, 255 = black. Convert to the other way round to match standard hex codes
        imArray = np.array(list(map(lambda x: self.numColours - 1 - x, imArray)))

        A = imArray.reshape((self.dim, self.dim))
        newdim = 20
        margin = int((self.dim - newdim) / 2)
        # todo deal with odd numbers
        imArray = A[margin:(self.dim - margin), margin:(self.dim - margin)].reshape(-1)
        self.dim = newdim

        return imArray

    def fetchICONasARRAY(self, id = None):
        self.dim = 16 # change if type of image changes,

        if id >= len(self.iconList) or id < 0:
            exit("Invalid icon id code {}".format(id))

        if id is None:
            id = random.randrange(0, len(self.iconList))


        iconLocs = "/home/saltermich1/PhDThesisGdrive/Code/NetworkData/MNIST"
        # todo sort this properly later, ie machine indep, change image etc
        filename = "{}/{}.ico".format(iconLocs, self.iconList[id])
        img = PIL.Image.open(filename)

        new_image = PIL.Image.new("RGBA", img.size, "WHITE")  # Create a white rgba background
        new_image.paste(img, (0, 0), img)  # Paste the image on the background. Go to the links given below for details.

        # colArray = np.asarray(img)
        imArray = np.asarray(new_image.convert(mode="L")).reshape(-1)

        return imArray

    def fetchSingleImage(self, type="MNIST", id = None):
        if type.upper() == "MNIST":
            imArray = self.fetchMNISTasARRAY(id)
        elif type.upper() == "ICON":
            imArray = self.fetchICONasARRAY(id)
        else:
            exit("invalid image type {}. Must be MNIST or ICON.".format(type))

        # down-sample to three colours (white grey black)
        # just doing it a stupid way to start with
        self.maxColourCode = 2 # *should* be able to convert to more shades easily, this way
        # todo make smarter

        # image = imArray.reshape(28, 28)
        # plt.imshow(image, cmap=mpl.cm.binary, interpolation="nearest")
        # plt.axis("off")


        imArray = np.array(list(map(lambda x: int(np.round(x / self.numColours * self.maxColourCode)), imArray)))



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
                    try:
                        edgeWeight = self.maxColourCode - np.abs(imArray[sourceID] - imArray[targetID])
                    except:
                        print("moocow")
                    graph.add_edge(sourceID, targetID, weight = edgeWeight)
                if j < self.dim - 1:
                    # new edge
                    targetID = (i) * self.dim + (j+1)
                    # so if adj vs have same/similar colour edge weight is *high*
                    # ie, small cuts between *different* colours
                    edgeWeight = self.maxColourCode - np.abs(imArray[sourceID] - imArray[targetID])
                    graph.add_edge(sourceID, targetID, weight = edgeWeight)

        dummy = 1
#        graph.delete_edges(weight_eq=0)
        # todo check if this is valid

        # visual_style = {}
        # visual_style["vertex_label"] = range(len(imArray))
        # color_list = ["black","grey","white"]
        # visual_style["edge_label"] = graph.es["weight"]
        # visual_style["bbox"] = (0, 0, 1000, 1000)
        #
        # try:
        #     visual_style["vertex_color"] = [color_list[colCode] for colCode in imArray]
        # except:
        #     print("moocow")
        #
        # graph.es["curved"] = 0
        # layout = graph.layout_grid()
        # output = ig.plot(graph, layout=layout, **visual_style)
        # outfile = "{}.png".format(self.iconList[id])
        # output.save(outfile)


        return graph



if __name__ == '__main__':
    print("running test on one image")
    M = ImageParser()
    M.fetchSingleImage(type="ICON", id = 0)