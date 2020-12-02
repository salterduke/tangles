import math
import numpy as np
import scipy as sp
import csv
import pandas as pd
import collections as coll
import igraph as ig
import itertools as it
import tools
# import treeswift
import ete3
import sys

class TangleSet():
    def __init__(self, job, log):
        self.separations = coll.defaultdict(list)
        self.definitelySmall = coll.defaultdict(list)

        self.TangleTree = ete3.Tree()
        self.TangleLists = coll.defaultdict(list)  # or set?
        self.currentTangle = 0
        self.job = job
        self.log = log
        self.nodeIndex = 1

        #### smart - creating hasse diag to keep track of maximal small sides
        # self.hasseDiagram = ig.Graph()
        self.maximalSeps = set()




    def findNextOrderSeparations(self):
        print("This must be overridden for each child")
        raise NotImplementedError

    #### delete depth later, when we're sure this works. ******
    def findAllTangles(self, depth=4):
        print("Finding All Tangles")
        # self.AvoidThreshold = threshold

        ### so that the first tangle tree starts at the root.
        self.TangleTree.add_feature("smallSides", [])

        self.findNextOrderSeparations()
        if self.kmin is None:
            print("Crack the shits, kmin not set")
            exit()
        # self.createHasseDiagram(self.kmin)
        self.kmax = self.kmin + depth

        self.TangleLists[self.kmin - 1] = [self.TangleTree]

        # ### see Evrendilek for heuristic?
        self.log.tick("kTangle{}".format(self.kmin))
        if not self.kTangle(self.kmin):
            print("No tangles exist")
            self.log.tock()
            return

        # # # ### find all tangles at each k
        for k in range(self.kmin+1, self.kmax+1):
            self.log.tick("kTangle{} Find seps".format(k))
            self.findNextOrderSeparations(k)
            self.log.tock()
            self.log.tick("kTangle{} Build Tangle".format(k))
            # self.createHasseDiagram(k)

            if not self.kTangle(k):
                print("NOT k: {}".format(k))
                self.tangleNum = k
                self.log.tock()
                ####### ********
                # return(self.commLists)
                # return(self.allBigSideProps)
                return(False)
            self.log.tock()

        self.log.end()
        ####### ********
        # return(self.commLists)
        return(self.allBigSideProps)

    def checkTangleAxioms(self, newSep):
        if newSep == set(list('ijk')):
            # print(side1, side2, newSep)
            print(self.smallSidesStack)

        ###### yes, I know this is probably not the best programming practice!!!
        if "Edge" in str(self.__class__):
            offset = 0
        elif "Vertex" in str(self.__class__):
            # because vertex tangles separate edges, and every single edge is
            # an order 2 separation, so we do this so we don't need to include
            # each edge explicitly
            # We *don't* want to do this for edge tangles
            print("Should probably worry about this when k=1 though...")
            print("Fix later..... ******")
            input("Press any key to continue")
            offset = 1
        else:
            print("New Tangle Type: Decide what you're doing here!")
            print("Define offset in checkTangleAxioms in baseTangleSet")
            exit()


        ### Axiom 2
        for side1, side2 in it.combinations(self.smallSidesStack, 2):
            # if newSep == set(list('ijk')):
            #     print(side1, side2, newSep)
            ###### ***** unnecessary, but maybe useful???
            # if newSep.issuperset(self.groundset - side1) or newSep.issuperset(self.groundset - side2):
            #     # print("**************")
            #     # print("{}\n is superset of  existing bigside".format(newSep))
            #     # print("**************")
            #     return False

            #### looks like shitty code, but might mean
            #### some of the unions can be avoided - O(m)
            double1 = side1 | newSep
            if len(double1) >= self.groundsetSize - offset:
                return False
            double2 = side2 | newSep
            if len(double2) >= self.groundsetSize - offset:
                return False
            triple = side1 | double2
            if len(triple) == self.groundsetSize:
                return False

        return True


    def kTangle(self, k):
        def sortB(sep):
            return(len(sep[1]))

        def formatSideName(side):   ##### fix this to switch on "verbose"
            # return ""
            side = str(side)
            side = side.replace("frozenset(", "").replace(")", "")
            side = side.replace("'", "")
            return side

        def addSideAsSep(side, parent, sepNum):
            # checkparent = set(['YBL014C', 'YCR042C', 'YDR145W', 'YDR167W', 'YDR392W', 'YDR448W', 'YER107C', 'YER110C', 'YER148W', 'YGL016W', 'YGL092W', 'YGL238W', 'YGL241W', 'YGR047C', 'YGR218W', 'YGR246C', 'YGR274C', 'YIL063C', 'YIL115C', 'YJL025W', 'YKL058W', 'YLR293C', 'YLR335W', 'YLR347C', 'YML007W', 'YML043C', 'YML114C', 'YMR005W', 'YMR047C', 'YMR125W', 'YMR308C', 'YNL039W', 'YNL189W', 'YOR098C', 'YOR160W', 'YOR194C', 'YPL082C', 'YPL178W'])
            # pset = set(parent.name.strip("{}").split(", "))

            if self.checkTangleAxioms(side):

                child = parent.add_child(name=formatSideName(side))
                self.nodeIndex+=1
                self.smallSidesStack.append(side)
                child.add_feature("smallSides", self.smallSidesStack.copy())
                prevBranches.append(child)


                if sepNum == numkSeps-1:
                    self.foundTangle +=1
                    self.TangleLists[k].append(child)
                    child.name = "T{}{}".format(self.currentTangle, formatSideName(side))
                    self.currentTangle += 1

                self.smallSidesStack.pop()


        self.foundTangle = 0

        numkdefSmall = len(self.definitelySmall[k])
        numkSeps = len(self.separations[k]) + numkdefSmall

        prevBranches = self.TangleLists[k-1]

        for sepNum in range(numkSeps):
            currentBranches = prevBranches
            prevBranches = []
            for truncTangle in currentBranches:
                self.smallSidesStack = truncTangle.smallSides   ###### *****
                if sepNum < numkdefSmall:
                    #### No branching
                    addSideAsSep(self.definitelySmall[k][sepNum], truncTangle, sepNum)
                elif sepNum < numkSeps:
                    ### check both sides of the separation
                    for side in self.separations[k][sepNum - numkdefSmall]:
                        addSideAsSep(side, truncTangle, sepNum)

        if self.foundTangle:
            self.finaliseAndPrint(k)
            return(True)
        else:
            return(False)


    def finaliseAndPrint(self, k):

        tangArray = pd.DataFrame()
        orientArray = pd.DataFrame()

        # self.commLists = coll.defaultdict(set)
        self.allBigSideProps = coll.defaultdict(dict)

        def printSingleTangle(tangle, tangID, tangOrder, verbose=False):

            tangOrients = []

            bigSideCounts = dict.fromkeys(self.groundset, 0)
            numSeps = len(tangle.smallSides)

            ##########################
            for i in range(numSeps):
                complement = self.groundset - tangle.smallSides[i]
                for element in complement:
                    bigSideCounts[element] += 1

                if verbose:
                    print("---------------------")
                    print("{}: order {}".format(tangID, tangOrder))
                    if len(tangle.smallSides[i]) <= len(complement):
                        print("{} / G\\*".format(set(tangle.smallSides[i])))
                    else:
                        print("G\\* / {}".format(complement))

                if len(tangle.smallSides[i]) <= len(complement):
                    tangOrients.append("{} / G\\*".format(sorted(set(tangle.smallSides[i]))))
                else:
                    tangOrients.append("G\\* / {}".format(sorted(complement)))

            bigSideProps = {element:(count/numSeps) for element, count in bigSideCounts.items()}

            # self.commLists[tangID]={element for element in self.groundset
            #     if bigSideProps[element] >= self.AvoidThreshold}

            tangHeader = "{}: ord={}".format(tangID, tangOrder)
            tangList = ["{}: {:.2f}".format(element, prop) for
                element, prop in sorted(bigSideProps.items(),
                key=lambda x:x[1], reverse=True)]
            tangArray[tangHeader] = tangList

            temparray = pd.DataFrame()
            temparray[tangHeader] = tangOrients
            nonlocal orientArray
            orientArray = pd.concat([orientArray, temparray], axis=1)

            # orientArray[tangHeader] = tangOrients

            self.allBigSideProps[tangID] = bigSideProps

        # print("**********************************")
        # print("Found {} tangles of order {}".format(self.foundTangle, k+1))
        # print("**********************************")
        tangleCounter = 0
        for i in range(self.kmin, k+1):
            for tangle in self.TangleLists[i]:
                printSingleTangle(tangle, tangleCounter, i+1, verbose=False)
                tangleCounter+=1

        AvoidOutfile = "{}/{}-AvoidProps.csv".\
        format(self.job['outputFolder'], self.job['outName'])
        tangArray.to_csv(AvoidOutfile)

        OrientsOutfile = "{}/{}-Orientations.csv".\
        format(self.job['outputFolder'], self.job['outName'])
        orientArray.to_csv(OrientsOutfile)


        with open('sepList.txt', 'w') as f:
            for i in [1,2,3]:
                f.write("sep order: {}\n".format(i))
                for item in [i]:
                    f.write("{}\n".format(item))

        # print(self.TangleTree)
        #####################
        doTreePrint = False
        if doTreePrint:
            outfile = "{}/{}-TangleTree-{}.png".\
                format(self.job['outputFolder'], self.job['outName'], k)

            style = ete3.NodeStyle()
            style["size"] = 0
            for n in self.TangleTree.traverse():
                n.set_style(style)

            ts = ete3.TreeStyle()
            # ts.show_branch_length = False
            ts.show_scale = False
            ts.show_leaf_name = False
            def my_layout(node):
                F = ete3.TextFace(node.name, tight_text=True)
                # F.fsize = 16
                F.rotation = -90
                ete3.add_face_to_node(F, node, column=0, position="branch-right")

            ts.layout_fn = my_layout
            ts.rotation = (90)

            # ts.min_leaf_separation = 10
            ts.branch_vertical_margin = 10
            # ts.allow_face_overlap = True

            try:
                self.TangleTree.render(outfile, tree_style=ts)
            except Exception as rendError:
                print(rendError)

            tidyTree = self.TangleTree.copy()
            # print(tidyTree.get_ascii(show_internal=True))
            for treenode in tidyTree.traverse():
                if "T" not in treenode.name:
                    treenode.delete(prevent_nondicotomic = False)
            ts.show_branch_length = False

            for node in tidyTree.iter_descendants():
                node.dist*=4

            ts.branch_vertical_margin = 8
            ts.scale = 360


            tidyOutfile = "{}/{}-TidyTangleTree-{}.pdf".\
                format(self.job['outputFolder'], self.job['outName'], k)
            try:
                tidyTree.render(tidyOutfile, tree_style=ts)
            except Exception as rendError:
                print(rendError)
