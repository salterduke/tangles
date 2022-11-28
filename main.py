import networkCD as netCD
import MNIST
import numpy as np
import datetime
import pandas as pd
import igraph as ig
import logger
import os
import sys
import shutil
import multiprocessing
import platform

if __name__ == '__main__':
    np.set_printoptions(precision=3)

    configFile = "const.txt"

    if len(sys.argv) >= 2 and "dev" in sys.argv[1].lower():
        print("Dev testing")
        testName = "DevYWS" # note leaving YWS in name so alg is correctly selected later
        configFile = "config2.txt"
    elif len(sys.argv) >= 2 and "YWS" in sys.argv[1]:
        print("Running YWS")
        testName = "TestYWS"
    elif len(sys.argv) >= 2 and "VY" in sys.argv[1]:
        print("Running VY")
        testName = "TestVY"
    else:
        exit("No test type specified. Either YWS or VY")

    log = logger.logger(testName)
    copyPics = False
# ------------------------------------------------------------------------------

def runAnalysis(job):

    log.log("Job: {}".format(job['outName']))
    ticktoken = log.tick("{} RunAnalysis".format(job['outName']))
    jobGraph = netCD.graphCD(job, log)

    n, m, tangCounts, timings = jobGraph.findTangleComms(dep = 3, sepsOnly=True)
    secs = log.tock(ticktoken)

    # jobGraph.overLapCliquePercolation()

    return([job['outName'], n, m, secs, "-".join(map(str, tangCounts)), "-".join(map(str, timings))])

def runMNIST(job, MNISTdata, imageID):
    job["outName"] = "MNIST_id{}".format(imageID)
    log.log("Job: {}, image {}".format(job["outName"], imageID))
    ticktoken = log.tick("{} RunAnalysis".format(job['outName']))
    job["doMNIST"] = True
    job["MNISTid"] = imageID
    job["MNISTdata"] = MNISTdata

    jobGraph = netCD.graphCD(job, log)

    n, m, tangCounts, timings = jobGraph.findTangleComms(dep = 0, sepsOnly=False)
    secs = log.tock(ticktoken)

    # jobGraph.overLapCliquePercolation()

    return([job['outName'], n, m, secs, "-".join(map(str, tangCounts)), "-".join(map(str, timings))])

def runMadeupGraph(job, N, M):
    job["outName"] = "constructed_{}_{}".format(N,M)
    job["construct"] = True
    job["N"] = N
    job["M"] = M
    log.log("Job: {}".format(job['outName']))
    ticktoken = log.tick("{} RunAnalysis".format(job['outName']))
    jobGraph = netCD.graphCD(job, log)

    n, m, tangCounts, timings = jobGraph.findTangleComms(dep = 4, sepsOnly=False)
    secs = log.tock(ticktoken)

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

    doConstructed = False
    doMNIST = False
    if doConstructed:
        # timing tests:
        job = {'outputFolder': "./output{}".format(testName), 'testName': testName}
        jobResults = []
        # for n in (20, 50, 100):
        for n in (150, 200):
            for m in range(n + 10, 3 * n+10, 20):
                jobres = runMadeupGraph(job, n, m)
                jobResults.append(jobres)
    elif doMNIST:
        job = {'outputFolder': "./output{}".format(testName), 'testName': testName}
        jobResults = []
        MNISTdata = MNIST.MNIST()

        ids = [0]
        # todo add different ids here
        for id in ids:
            jobres = runMNIST(job, MNISTdata, id)
            jobResults.append(jobres)
    else:
        jobResults = []
        for index, job in jobsToRun.iterrows():
            job['outputFolder'] = "./output{}".format(testName)
            job['testName'] = testName
            jobres = runAnalysis(job)
            jobResults.append(jobres)


    if copyPics:
        copyPicsToLatex()

    resultsDF = pd.DataFrame(jobResults, columns=['network', 'Vs', 'Es', 'secs', 'tangCounts', 'timings'])

    resultsDF.to_csv("./output{}/results{}.csv".format(testName, log.timeString))

    log.end()
