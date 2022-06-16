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

    # configFile = "configConst.txt"
    # testName = "Constructed"
    configFile = "config2.txt"
    testName = "TestCores"

    log = logger.logger(testName)
    copyPics = False
# ------------------------------------------------------------------------------

def runAnalysis(job):

    log.log("Job: {}".format(job['outName']))
    ticktoken = log.tick("{} RunAnalysis".format(job['outName']))
    jobGraph = netCD.graphCD(job, log)

    n, m, tangCounts, timings = jobGraph.findTangleComms(dep = 2, sepsOnly=True)
    secs = log.tock(ticktoken)

    # jobGraph.overLapCliquePercolation()

    return([job['outName'], n, m, secs, "-".join(map(str, tangCounts)), "-".join(map(str, timings))])

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

    jobResults = []

    for index, job in jobsToRun.iterrows():
        job['outputFolder'] = "./output{}".format(testName)
        jobres = runAnalysis(job)
        jobResults.append(jobres)

    # timing tests:
    # jobResults = []
    # for index, job in jobsToRun.iterrows():
    #     job['outputFolder'] = "./output{}".format(testName)
    #     jobres = runAnalysis(job)
    #     jobResults.append(jobres)


    if copyPics:
        copyPicsToLatex()

    resultsDF = pd.DataFrame(jobResults, columns=['network', 'Vs', 'Es', 'secs', 'tangCounts', 'timings'])

    resultsDF.to_csv("./output{}/results{}.csv".format(testName, log.timeString))

    log.end()
