import numpy as np
from time import time
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rc('axes', labelsize=14)
mpl.rc('xtick', labelsize=12)
mpl.rc('ytick', labelsize=12)


class MNIST():
    def __init__(self):
        try:
            from sklearn.datasets import fetch_openml
            self.mnist = fetch_openml('mnist_784', version=1, cache=True)
            self.mnist.target = mnist.target.astype(np.int8) # fetch_openml() returns targets as strings
            dummy = 1
        except ImportError:
            from sklearn.datasets import fetch_mldata
            self.mnist = fetch_mldata('MNIST original')

    def fetchSingleAsArray(self, id = None):
        if id is None:
            id = random.randrange(0, self.mnist.data.shape[0])
        # if id < 0 or id >= self.mnist.data.shape[0]:
            # some error
        imArray = self.mnist.data[id]

        # down-sample to three colours (white grey black)
        # just doing it a stupid way to start with
        # todo make smarter
        imArray = np.array(list(map(lambda x: round(x / 128), imArray))).reshape(28, 28)
        return imArray

    def fetchSingleAsGraph(self, id = None):
        singleArr = self.fetchSingleAsArray(id)
