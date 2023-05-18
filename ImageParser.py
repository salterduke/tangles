import igraph as ig
import numpy as np
import pandas as pd
import random
import matplotlib as mpl
import PIL

import matplotlib.pyplot as plt
mpl.rc('axes', labelsize=14)
mpl.rc('xtick', labelsize=12)
mpl.rc('ytick', labelsize=12)


# inFile;outName;depth;numColours;cropsize;rowoffset;coloffset
# ../NetworkData/Images/photos/monalisa0.png;monalisa0;14;4;14;2;14

class ImageParser():
    def __init__(self):
        self.numColoursOrig = 256   # change if type of image changes
        self.mnist = None
        self.dim = 40

    def initMNIST(self):
        try:
            from sklearn.datasets import fetch_openml
            self.mnist = fetch_openml('mnist_784', version=1, cache=True)
            self.mnist.target = self.mnist.target.astype(np.int8) # fetch_openml() returns targets as strings
            dummy = 1
        except ImportError:
            from sklearn.datasets import fetch_mldata
            self.mnist = fetch_mldata('MNIST original')

    def crop(self, longA, newdim, newdim2, rowstart = None, colstart = None):
        # all right, currently input image has to be square, but can be cropped not square
        # change if needed
        A = longA.reshape((self.dim, self.dim))

        if pd.isna(rowstart) or rowstart is None:
            rowstart = int((self.dim - newdim) / 2) # int truncates if odd, which is fine here
        rowend = int(rowstart + newdim)

        if rowstart < 0 or rowend > self.dim:
            exit("crack the sads, row crop outside existing border: rowstart {}, rowend {}, current dim {} ".format(rowstart, rowend, self.dim))

        if pd.isna(colstart) or colstart is None:
            colstart = int((self.dim - newdim2) / 2)  # int truncates if odd, which is fine here
        colend = int(colstart + newdim2)

        if colstart < 0 or colend > self.dim:
            exit("crack the sads, col crop outside existing border: colstart {}, colend {}, current dim {} ".format(colstart, colend, self.dim))

        imArray = A[rowstart:rowend, colstart:colend].reshape(-1)
        self.dim = newdim
        self.dim2 = newdim2
        return imArray


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
        imArray = np.array(list(map(lambda x: self.numColoursOrig - 1 - x, imArray)))

        return imArray

    def fetchIMAGEasARRAY(self, filename, transparent=False):

        img = PIL.Image.open(filename).convert("L")

        if transparent:
            new_image = PIL.Image.new("L", img.size, "WHITE")  # Create a white rgba background
            new_image.paste(img, (0, 0), img)  # Paste the image on the background. Go to the links given below for details.
            img = new_image

        # colArray = np.asarray(img)
        imArray = np.asarray(img).reshape(-1)
        if img.size[0] != img.size[1]:
            exit("crack the sads, image not square")
        self.dim = img.size[0]

        return imArray

    def fetchSingleImage(self, job):
        # provide some defaults. Default default is None
        imtype = job.get("imType", "IMAGE")
        id = job.get("MNISTid")
        numColours = job.get("numColours", 256)
        cropsize = job.get("cropsize", 16)
        rowoffset = job.get("rowoffset", None)
        coloffset = job.get("coloffset", None)
        doDiagonal = job.get("doDiagonal", False)


        crop2 = job.get("crop2", None)
        if pd.isna(crop2) or crop2 is None:
            crop2 = cropsize
            job["crop2"] = cropsize
        else:
            # don't know why it's not int already, but at this point I don't care...
            crop2 = int(crop2)
            job["crop2"] = crop2


        job["outName"] = "{}-{}cols-size{}x{}-at-{}-{}".format(job["outName"], numColours, cropsize, crop2, rowoffset, coloffset)

        if imtype == "MNIST":
            imArray = self.fetchMNISTasARRAY(id)
        elif imtype == "IMAGE":
            imArray = self.fetchIMAGEasARRAY(job.get("inFile"))
        else:
            exit("invalid image type {}. Must be MNIST or IMAGE.".format(imtype))

        imArray = self.crop(imArray, cropsize, crop2, rowoffset, coloffset)

        self.numColoursNew = numColours

        # image = imArray.reshape(28, 28)
        # plt.imshow(image, cmap=mpl.cm.binary, interpolation="nearest")
        # plt.axis("off")

        imArray = np.array(list(map(lambda x: \
                    int(np.round(x / (self.numColoursOrig-1) * (self.numColoursNew -1))), imArray)))

        graph = ig.Graph()
        graph.add_vertices(self.dim*self.dim2)

        # so if adj vs have same/similar colour edge weight is *high*
        # ie, small cuts between *different* colours
        weightExp = 1

        numSubtractFrom = self.numColoursNew - 1
        # numSubtractFrom = self.numColoursNew


        # iterate over *nodes*
        for i in range(self.dim): # rows
            for j in range(self.dim2): # columns
                sourceID = i * self.dim2 + j
                graph.vs[sourceID]["name"] = "{}_{}".format(i,j)
                # add only edges going "right" or "down" to avoid double count
                if i < self.dim - 1:
                    # add right edge
                    targetID = (i+1) * self.dim2 + j
                    # -1 is offset for 0 index
                    edgeWeight = numSubtractFrom - np.abs(imArray[sourceID] - imArray[targetID])

                    edgeWeight = edgeWeight**weightExp
                    graph.add_edge(sourceID, targetID, weight = edgeWeight)
                if j < self.dim2 - 1:
                    # add down edge
                    targetID = (i) * self.dim2 + (j+1)
                    # -1 is offset for 0 index
                    edgeWeight = numSubtractFrom - np.abs(imArray[sourceID] - imArray[targetID])

                    edgeWeight = edgeWeight**weightExp
                    graph.add_edge(sourceID, targetID, weight = edgeWeight)
                # add diagonal edges
                if doDiagonal:
                    if i < self.dim - 1 and j < self.dim2 - 1:
                        #add down right diagonal

                        targetID = (i+1) * self.dim2 + (j+1)
                        # -1 is offset for 0 index
                        edgeWeight = numSubtractFrom  - np.abs(imArray[sourceID] - imArray[targetID])

                        edgeWeight = edgeWeight**weightExp
                        graph.add_edge(sourceID, targetID, weight = edgeWeight)

                    if i > 0 and j < self.dim2 - 1:
                        #add down left diagonal
                        targetID = (i-1) * self.dim2 + (j+1)
                        # -1 is offset for 0 index
                        edgeWeight = numSubtractFrom  - np.abs(imArray[sourceID] - imArray[targetID])

                        edgeWeight = edgeWeight**weightExp
                        graph.add_edge(sourceID, targetID, weight = edgeWeight)

        labels = ""
        visual_style = {}
        if labels == "name":
            visual_style["vertex_label"] = graph.vs["name"]
        elif labels == "index":
            visual_style["vertex_label"] = graph.vs.indices
        else:
            visual_style["vertex_label"] = ""

        if labels:
            visual_style["edge_label"] = graph.es["weight"]
            visual_style["bbox"] = (0, 0, self.dim2 * 50, self.dim * 50)
        else:
            visual_style["bbox"] = (0, 0, self.dim2 * 30, self.dim * 30)

        pal = ig.GradientPalette("black", "white", self.numColoursNew)

        graph.vs["vertex_color"] = [pal[int(colCode)] for colCode in imArray]
        visual_style["vertex_color"] = graph.vs["vertex_color"]

        graph.es["curved"] = 0
        layout = graph.layout_grid(self.dim2)
        output = ig.plot(graph, layout=layout, **visual_style)
        outfile = "{}/{}_grid.png".format(job["outputFolder"], job["outName"])
        output.save(outfile)

        return graph

    def createRandomBackground(self, propBlack):

        randBits = np.random.choice((0,1), p=(propBlack, 1-propBlack), size=self.dim**2)
        dummy = 1

if __name__ == '__main__':
    print("running test on one image")
    M = ImageParser()

    M.createRandomBackground(0.5)

    # job = {"imType":"IMAGE",
    #  "MNISTid":0,
    #  "numColours":3,
    #  "cropsize":16,
    #  "rowoffset":2,
    #  "coloffset":12}
    # M.fetchSingleImage(job)

