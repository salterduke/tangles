import datetime
import os
import pandas as pd
# from collections import defaultdict

class logger():
    def __init__(self, testName):
        self.startTime = datetime.datetime.utcnow()
        logDir = "./output{}/logs".format(testName)
        if not os.path.exists(logDir):
            os.makedirs(logDir)

        self.filename = "{}/{}{}.log".format(logDir, testName, str(self.startTime).replace(":","."))
        # self.datafile = "{}/{}.csv".format(logDir, testName)
        # #### Check if exists!!!!!!
        # self.infoData = pd.read_csv(self.datafile, index_col=0)

        text = "Test: {}, time: {}\n".format(testName, self.startTime)
        with open(self.filename, 'w+') as the_file:
            the_file.write(text)
        print(text)

    def log(self, text):
        with open(self.filename, 'a') as the_file:
            the_file.write(text + "\n")
        print(text)

    def tick(self, taskName):
        self.jobstartTime = datetime.datetime.utcnow()
        self.currTaskName = taskName

    def tock(self):
        endTime = datetime.datetime.utcnow()
        seconds = (endTime - self.jobstartTime).total_seconds()
        self.log("Job {} took {} seconds".format(self.currTaskName,seconds))
        self.log("-----------------------")

    # def newRecord(self, recName):
    #     self.infoData[recName]
    #
    # def info(self, key, value):
    #     self.infoData[key] = value

    def end(self):
        endTime = datetime.datetime.utcnow()
        seconds = (endTime - self.startTime).total_seconds()
        self.log("Total Running Time {} seconds".format(seconds))
        self.log("*****************************************")
