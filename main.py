import networkCD as netCD
import numpy as np
import datetime
import pandas as pd
import igraph as ig
import logger
import os
import shutil
import multiprocessing
import platform

if __name__ == '__main__':
    np.set_printoptions(precision=3)

    configFile = "config2.txt"
    # configFile = "configConst.txt"
    testName = "Constructed"

    log = logger.logger(testName)
    copyPics = False
# ------------------------------------------------------------------------------

def runAnalysis(job):

    log.log("Job: {}".format(job['outName']))
    log.tick("{} RunAnalysis".format(job['outName']))
    jobGraph = netCD.graphCD(job, log)

    jobGraph.findOverLapCommunities()

    log.tock()

def copyPicsToLatex():
    picFolder = "./output{}".format(testName)
    latexFolder = "../../Reports_Thesis/Latex/Images"

    pics = os.listdir(picFolder)
    for pic in pics:
        if ".pdf" in pic:
            fullName = os.path.join(picFolder, pic)
            if os.path.isfile(fullName):
                shutil.copy(fullName, latexFolder)

# ------------------------------------------------------------------------------
# Set name of graph file to input, and name to use for output files
# In config file
if __name__ == '__main__':
    if platform.system() == 'Linux':
        print("Running linux")
        multiprocessing.set_start_method('forkserver')

    jobsToRun = pd.read_csv(configFile, delimiter=';', header=0, comment="#")

    for index, job in jobsToRun.iterrows():
        job['outputFolder'] = "./output{}".format(testName)
        runAnalysis(job)

    if copyPics:
        copyPicsToLatex()

    log.end()
