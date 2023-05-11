import Bio.UniProt.GOA as GOA
import os
from ftplib import FTP
import gzip
import pandas as pd
import numpy as np
from Modules import baseChecker as bch


class protChecker(bch.commChecker):
    def __init__(self, protList, organism = "Yeast"):

        bch.commChecker.__init__(self, protList)

        protSet = set(protList)
        ## do once as 'in' is more efficient for sets

        #### put shit here for other organisms if necessary
        org_uri = "/pub/databases/GO/goa/YEAST/goa_yeast.gaf.gz"
        org_fn = org_uri.split('/')[-1]
        data_folder = "./GOfiles/"


        ### setting shit up
        TAMFile = "./GOfiles/GO_FS.txt"
        # TAMFile = "./GOfiles/Short.txt"
        headers=["A", "B", "LCAcount"]
        self.protAliases = pd.DataFrame(index = protList,
            columns=["TAMalias", "GOalias", "GOcount"])
        self.TAMScores = pd.read_csv(TAMFile, sep="\t", header=None,
            names=headers, index_col=(0,1))

        # Check if the file exists already
        org_gaf = os.path.join(data_folder, org_fn)
        if(not os.path.isfile(org_gaf)):
            # Login to FTP server
            ebi_ftp = FTP('ftp.ebi.ac.uk')
            ebi_ftp.login() # Logs in anonymously

            # Download
            with open(org_gaf,'wb') as org_fp:
                ebi_ftp.retrbinary('RETR {}'.format(org_uri), org_fp.write)

            # Logout from FTP server
            ebi_ftp.quit()


        # File is a gunzip file, so we need to open it in this way
        with gzip.open(org_gaf, 'rt') as org_gaf_fp:
            self.org_funcs = {}  # Initialise the dictionary of functions

            # Iterate on each function using Bio.UniProt.GOA library.
            for entry in GOA.gafiterator(org_gaf_fp):
                prot = entry.pop('DB_Object_Symbol')
                if prot in protSet:
                    self.protAliases.loc[prot]["GOalias"] = prot
                    key = prot
                else:
                    syns = [syn for syn in entry['Synonym'] if syn in protSet]
                    if len(syns) > 1:
                        print("WTF, found too many syns for {}".format(prot))
                        print(syns)
                        for syn in syns:
                            self.protAliases.loc[syn]["GOalias"] = prot
                    elif len(syns) == 1:
                        self.protAliases.loc[syns[0]]["GOalias"] = prot
                        key = syns[0]
                    else:
                        ### not in protlist - therefore don't care about it!
                        continue

                ### Do not add GO term if qualifier is or contains 'NOT'
                if not any('NOT' in word for word in entry['Qualifier']):
                    if prot not in self.org_funcs:
                        self.org_funcs[prot] = {
                            'Name': entry['DB_Object_Name'],
                            'Synonym': entry['Synonym'],
                            'GoTerms': {entry['GO_ID']: entry['Aspect']}
                        }
                    else:
                        self.org_funcs[prot]['GoTerms'][entry['GO_ID']] =\
                            entry['Aspect']

                    if key == "YER110C":
                        dummy = 1

                    if entry["GO_ID"] not in self.compareCover.columns:
                        # if col not exist, init
                        self.compareCover[entry["GO_ID"]] = np.zeros(self.nNodes, dtype=int)
                    self.compareCover.loc[key, entry["GO_ID"]] = 1

                    if prot in self.TAMScores.index:
                        # print("Found {}".format(prot))
                        self.protAliases.loc[key]["TAMalias"] = prot
                    else:
                        for synonym in entry["Synonym"]:
                            if synonym in self.TAMScores.index:
                                # print("{} -> {}".format(prot, synonym))
                                self.protAliases.loc[key]["TAMalias"] = synonym
                                break

        # self.getGOcounts(protList)  # todo check what this is doing and if I need it

        # Based on Yu et al 2007 and 2008
        #########
        size = self.nNodes
        simMatrix = np.zeros((size, size))

        nProts = max(len(self.TAMScores.index.levels[0]), len(self.TAMScores.index.levels[1]))
        N = nProts * (nProts - 1) / 2

        totalSim = 0
        for i in range(size):
            # simMatrix[i][i] = 1   ### Actually, probably not... Unclear.
            for j in range(i+1, size):
                # print((protList[i], protList[j]))
                # n = TAMScores.loc[(protList[i], protList[j])]["LCAcount"]
                if pd.isnull(self.protAliases.iloc[i]["TAMalias"]) or pd.isnull(self.protAliases.iloc[j]["TAMalias"]):
                    n = N   ### quick kludge to make sim = 0
                elif self.protAliases.iloc[i]["TAMalias"] == self.protAliases.iloc[j]["TAMalias"]:
                    n = 1  ### if we have two different prots with the same synonyms
                    ## we can reasonably conclude they're similar (perhaps subunits of same prot)
                else:
                    if (self.protAliases.iloc[i]["TAMalias"], self.protAliases.iloc[j]["TAMalias"]) in self.TAMScores.index:
                        n = self.TAMScores.loc[(self.protAliases.iloc[i]["TAMalias"],
                            self.protAliases.iloc[j]["TAMalias"])]["LCAcount"]
                    elif (self.protAliases.iloc[j]["TAMalias"], self.protAliases.iloc[i]["TAMalias"]) in self.TAMScores.index:
                        n = self.TAMScores.loc[(self.protAliases.iloc[j]["TAMalias"],
                            self.protAliases.iloc[i]["TAMalias"])]["LCAcount"]
                    else:
                        n = N # makes sim = 0


                ### Ahn et al uses this cutoff
                simMatrix[i][j] = simMatrix[j][i] = 1 if n/N < 10**(-3) else 0
                totalSim+= simMatrix[i][j]
                # print("{}: {}".format(n/N, simMatrix[i][j]))

        self.simDF = pd.DataFrame(data=simMatrix,
            index=self.protAliases.index, columns=self.protAliases.index)

        self.aveSim = totalSim / (size * (size-1) / 2)



    def query(self, protID):
        print("-----------")
        print(protID)

        if protID in self.org_funcs:
            # print(self.org_funcs[protID])
            print("Found: {}".format(protID))
        else:
            for key in self.org_funcs:
                if protID in self.org_funcs[key]["Synonym"]:
                    print("Found as synonym")
                    print("{}: {}".format(protID, self.org_funcs[key]["Synonym"]))
                else:
                    print("Well fuck, not found as synonym")

            ## add to dictionary to save time on subsequent lookups
            self.org_funcs[protID] = self.org_funcs[key]

    def getGOcounts(self, protList):

        allGOTerms = set()

        for prot in self.protAliases.index:
            # if prot in self.org_funcs:
            #     self.protAliases.loc[prot]["GOcount"] = len(self.org_funcs[prot]["GoTerms"].keys())
            # else:
            #     # print("Can't find!!!!! {}".format(prot))
            #     self.protAliases.loc[prot]["GOcount"] = 0

            if not pd.isnull(self.protAliases.loc[prot]["GOalias"]):
                self.protAliases.loc[prot]["GOcount"] = len(self.org_funcs[self.protAliases.loc[prot]["GOalias"]]["GoTerms"].keys())

                ##### Finds out what/ how many GO terms in whole NW
                allGOTerms.update(self.org_funcs[self.protAliases.loc[prot]["GOalias"]]["GoTerms"].keys())
            else:
                self.protAliases.loc[prot]["GOcount"] = 0


        # self.protAliases.fillna(value = {"GOcount": 0}, inplace = True)


         #self.protAliases.apply(lambda row: len(self.org_funcs[row.GOalias]["GoTerms"].keys()), axis=1)
        # self.protAliases["GOcount"] = \
        #     [len(self.org_funcs[row.GOalias]["GoTerms"].keys()) for index, row in self.protAliases.iterrows()]

        # print("NUmber in TAM but not GOA: {}".format(TAMnotGOA))
        self.overlapMetadata = self.protAliases["GOcount"]
        self.protAliases.sort_values(by="TAMalias", inplace=True)

        ##### Finds out what/ how many GO terms in whole NW
        print("Total number of GO Terms: {} ".format(len(allGOTerms)))



