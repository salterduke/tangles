import networkCD as netCD
import numpy as np
import datetime
import pandas as pd
import igraph as ig
import logger
import os
import shutil

np.set_printoptions(precision=3)

configFile = "config.txt"
# testName = "BacktoStupid"
# testName = "GHU_all"
testName = "Cutfinder"

log = logger.logger(testName)
copyPics = False
# ------------------------------------------------------------------------------

def runAnalysis(job, method):

    log.log("Job: {}".format(job['outName']))
    log.tick("{} {}".format(job['outName'], method))
    jobGraph = netCD.graphCD(job, log)
    jobGraph.method = method

    jobGraph.findOverLapCommunities(method)

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
jobsToRun = pd.read_csv(configFile, delimiter=';', header=0, comment="#")

for index, job in jobsToRun.iterrows():
    job['outputFolder'] = "./output{}".format(testName)

    for method in (["tangles"]):
        runAnalysis(job, method)

if copyPics:
    copyPicsToLatex()

log.end()
