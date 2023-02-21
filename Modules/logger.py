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

        self.timeString = str(self.startTime).replace(":",".")

        text = "Test: {}, time: {}\n".format(testName, self.startTime)
        with open(self.filename, 'w+') as the_file:
            the_file.write(text)
        print(text)
        self.tokenID = -1
        self.allStartTimes = []
        self.allTaskNames = []

    def log(self, text):
        with open(self.filename, 'a') as the_file:
            the_file.write(text + "\n")
        print(text)

    def tick(self, taskName):
        self.tokenID+=1
        self.allStartTimes.append(datetime.datetime.utcnow())
        self.allTaskNames.append(taskName)
        # todo should really check the indexes work
        return(self.tokenID)

    def tock(self, token=None):
        if token is None:
            token = self.tokenID #ie, latest tick
        endTime = datetime.datetime.utcnow()
        seconds = round((endTime - self.allStartTimes[token]).total_seconds(), 2)
        self.log("Job {} took {} seconds".format(self.allTaskNames[token],seconds))
        self.log("-----------------------")
        return(seconds)

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
