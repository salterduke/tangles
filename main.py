import networkCD as netCD
import ImageParser
import numpy as np
import datetime
import pandas as pd
import igraph as ig
import Modules.logger as logger
import os
import sys
import shutil
import multiprocessing
import platform
import ImageParser

# print("{}: Outside EVERYTHING ___________________________".format(__name__))

if __name__ == '__main__':
    np.set_printoptions(precision=3)

    # default
    configFile = "const.txt"
    # runDepth = 4
    printTimings = False

    if len(sys.argv) >= 2 and "dev" in sys.argv[1].lower():
        # testName = "DevYWS" # note leaving YWS in name so alg is correctly selected later
        testName = sys.argv[1]
        print("Dev testing, name {}".format(testName))
        configFile = "config2.txt"
        # runDepth = 6
    elif len(sys.argv) >= 2 and "img" in sys.argv[1].lower():
        print("Dev Image testing")
        testName = "DevVY" + sys.argv[1] # note leaving VY in name so alg is correctly selected later
        configFile = "configImage.txt"
    elif len(sys.argv) >= 2 and "YWS" in sys.argv[1]:
        print("Running YWS")
        testName = "TestYWS"
        printTimings = True
    elif len(sys.argv) >= 2 and "VY" in sys.argv[1]:
        print("Running VY")
        testName = "TestVY"
        printTimings = True
    else:
        exit("No test type specified. Either YWS or VY")

    log = logger.logger(testName)
    copyPics = False
# ------------------------------------------------------------------------------

def runAnalysis(job):

    log.log("Job: {}".format(job['outName']))
    ticktoken = log.tick("{} RunAnalysis".format(job['outName']))
    jobGraph = netCD.graphCD(job, log)

    depDefault = 4
    runDepth = job.get("depth", depDefault)
    if np.isnan(runDepth):
        runDepth = depDefault

    # modified so dep is total number of orders, not total after the first one
    # n, m, tangCounts, timings, sepCounts = jobGraph.findTangleComms(dep = runDepth, sepsOnly=False)
    n, m, tangCounts, timings, sepCounts = jobGraph.findTangleComms(dep = runDepth, maxEmptyOrders=4, sepsOnly=False)


    secs = log.tock(ticktoken)


    return([job['outName'], n, m, secs, "-".join(map(str, tangCounts)), "-".join(map(str, timings)), "-".join(map(str, sepCounts))])

# def runMadeupGraph(job, N, M):
#     job["outName"] = "constructed_{}_{}".format(N,M)
#     job["construct"] = True
#     job["N"] = N
#     job["M"] = M
#     log.log("Job: {}".format(job['outName']))
#     ticktoken = log.tick("{} RunMadeUpGraph".format(job['outName']))
#     jobGraph = netCD.graphCD(job, log)
#
#     n, m, tangCounts, timings, sepCounts = jobGraph.findTangleComms(dep = 4, sepsOnly=False)
#     secs = log.tock(ticktoken)
#
#     return([job['outName'], n, m, secs, "-".join(map(str, tangCounts)), "-".join(map(str, timings)), "-".join(map(str, sepCounts))])
#

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
    print("-----------------------------------------")
    print(jobsToRun)
    print("-----------------------------------------")

    doConstructed = False
    if doConstructed:
        # timing tests:
        job = {'outputFolder': "./output{}".format(testName), 'testName': testName}
        jobResults = []
        # for n in (20, 50, 100):
        for n in (150, 200):
            for m in range(n + 10, 3 * n+10, 20):
                jobres = runMadeupGraph(job, n, m)
                jobResults.append(jobres)
    else:
        jobResults = []
        # creating the parser here, so it's a single object, and thus MNIST data doesn't have to
        # be recreated for each new job
        parser = ImageParser.ImageParser()
        for index, job in jobsToRun.iterrows():
            job['outputFolder'] = "./output{}".format(testName)
            job['testName'] = testName
            job["imParser"] = parser
            job = {key: val for key, val in job.items() if not pd.isna(val)}
            jobres = runAnalysis(job)
            jobResults.append(jobres)

    # elif doImage:
    #
    #     jobResults = []
    #     parser = ImageParser.ImageParser()
    #
    #     ids = range(1)
    #     # todo add different ids here
    #     for id in ids:
    #         job["numColours"] = 3
    #         job["cropsize"] = 16
    #         job["rowoffset"] = 2  # get right hand side, not centre
    #         job["coloffset"] = 12  # get right hand side, not centre
    #         jobres = runImage(job, parser, imageType, id)
    #         jobResults.append(jobres)


    if copyPics:
        copyPicsToLatex()

    if printTimings:
        resultsDF = pd.DataFrame(jobResults, columns=['network', 'Vs', 'Es', 'secs', 'tangCounts', 'timings', 'sepCounts'])
        resultsDF.to_csv("./output{}/results{}.csv".format(testName, log.timeString))

    log.end()
