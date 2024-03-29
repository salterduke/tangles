import pandas as pd
import collections as coll
import itertools as it
import ete3
import numpy as np

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
        self.smallSidesList = []

    def getTangleCounts(self):
        countList = list()
        for k in range(self.kmin, self.kmax):
            # note k is SEP order, so +1 for tang order
            tangNumbers = (k+1, len(self.TangleLists[k]))
            countList.append(tangNumbers)
        return(countList)

    def findNextOrderSeparations(self):
        print("This must be overridden for each child")
        raise NotImplementedError

    #### delete depth later, when we're sure this works. ******
    def findAllTangles(self, depth=4, maxEmptyOrders=4, sepsOnly = False):
        print("Finding All Tangles")

        orderCount = 0
        emptyCount = 0

        timings = []
        sepCounts = []

        ### so that the first tangle tree starts at the root.
        self.TangleTree.add_feature("smallSides", [])

        print("-------------------------------------------")
        self.log.tick("kTangle min Find seps")
        maxdepth = depth + maxEmptyOrders
        numSeps = self.findNextOrderSeparations(None, maxdepth)
        orderCount+=1
        if self.kmin is None:
            print("Crack the shits, kmin not set")
            exit()
        timings.append(self.log.tock())
        sepCounts.append(numSeps)
        # self.kmax = self.kmin + depth

        self.TangleLists[self.kmin - 1] = [self.TangleTree]

        if not sepsOnly:
            # ### see Evrendilek for heuristic?
            self.log.tick("kTangle{} (min) Build Tangle".format(self.kmin))
            if not self.kTangle(self.kmin):
                self.log.log("No tangles exist at min order {}".format(self.kmin))
                self.log.tock()
                return(timings, sepCounts)
            else:
                self.log.tock()

        # # # ### find all tangles at each k
        k = self.kmin
        while orderCount < depth and emptyCount < maxEmptyOrders:
        # for k in range(self.kmin+1, self.kmax+1):
            print("-------------------------------------------")
            k+=1
            self.log.tick("kTangle{} Find seps".format(k))
            numSeps = self.findNextOrderSeparations(k, maxdepth)
            timings.append(self.log.tock())
            sepCounts.append(numSeps)
            # note that this is the total seps found, which could be different to the number of seps stored,
            # as not all seps are necessary

            # check if found AND STORED any separations
            kSeps = len(self.separations[k]) + len(self.definitelySmall[k])
            self.log.log("Found {} seps at order k = {}".format(kSeps, k))
            if kSeps > 0:
                orderCount+=1
            else:
                emptyCount+=1

            if not sepsOnly and kSeps > 0:
                self.log.tick("kTangle{} Build Tangle".format(k))
                if not self.kTangle(k):
                    self.log.log("No tangles at tangle order k: {}".format(k+1))
                    self.kmax = k-1  # todo: make sure not off-by-one
                    self.log.tock()
                    self.log.end()
                    return(timings, sepCounts)
                else:
                    self.log.log("Found tangles at tangle order k: {}".format(k+1))
                self.log.tock()
        self.log.end()
        return(timings, sepCounts)


    def checkTangleAxioms(self, newSep):
        # return vals are passesCheck, and addToList, for seps that are subsets so not needed
        # (both bool)

        ### Axiom 2
        # todo - reverting to storing all def small in stack to handle cases where *all* def small
        # todo - consider doing properly later
        # sepsSoFar = list(it.chain(*self.definitelySmall.values(), self.smallSidesStack))
        sepsSoFar = self.smallSidesList

        if sepsSoFar is not None:
            self.keepSeps = np.ones(len(sepsSoFar), dtype = int)

            if len(sepsSoFar) == 1:
                side1 = sepsSoFar[0]
                if newSep.issubset(side1):
                    # ie, it's okay, but don't need it
                    return True, False
                elif newSep.issuperset(side1):
                    self.keepSeps[0] = 0
                double1 = side1 | newSep
                # if len(double1) > self.groundsetSize or (-1 not in double1 and len(double1) == self.groundsetSize):
                #     return False, False
                leafCount=0
                for l in range(-1, self.leafExtent-1,-1):
                    leafCount+=sum(1 for sep in (side1, newSep) if l in sep)
                leftOut = self.groundset - double1
                if len(leftOut) == 0 or (len(leftOut) <= leafCount and leftOut.issubset(self.leaves)):
                    return False, False
                # leafCount = max(0, leafCount-1) # working out how many to disregard, but need to count one -1 if present
                # todo put this back in if not working
                # if len(double1) + leafCount >= self.groundsetSize:
                #     return False, False
            else:
                checkedSubsets = False
                # using a list so easier to double iterate
                for id1, side1 in enumerate(sepsSoFar[:len(sepsSoFar)-1]):
                    if not checkedSubsets:
                        if newSep.issubset(side1):
                            # ie, it's okay, but don't need it
                            return True, False
                        elif newSep.issuperset(side1):
                            self.keepSeps[id1] = 0

                    double1 = side1 | newSep
                    # if len(double1) > self.groundsetSize or (-1 not in double1 and len(double1) == self.groundsetSize):
                    #     return False, False
                    leafCount = 0
                    for l in range(-1, self.leafExtent - 1, -1):
                        leafCount += sum(1 for sep in (side1, newSep) if l in sep)
                    leftOut = self.groundset - double1
                    if len(leftOut) == 0 or (len(leftOut) <= leafCount and leftOut.issubset(self.leaves)):
                        return False, False
                    # leafCount = max(0,leafCount - 1)  # working out how many to disregard, but need to count one -1 if present
                    # todo put back in if not working
                    # if len(double1) + leafCount >= self.groundsetSize:
                    #     return False, False

                    for id2 in range(id1 + 1, len(sepsSoFar)):
                        side2 = sepsSoFar[id2]
                        if not checkedSubsets:
                            if newSep.issubset(side2):
                                return True, False
                            elif newSep.issuperset(side2):
                                self.keepSeps[id2] = 0
                        triple = side2 | double1
                        # if len(triple) > self.groundsetSize or (-1 not in triple and len(triple) == self.groundsetSize):
                        #     return False, False
                        leafCount = 0
                        for l in range(-1, self.leafExtent - 1, -1):
                            leafCount += sum(1 for sep in (side1, side2, newSep) if l in sep)
                        leftOut = self.groundset - triple
                        if len(leftOut) == 0 or (len(leftOut) <= leafCount and leftOut.issubset(self.leaves)):
                            return False, False
                        # leafCount = max(0, leafCount - 1)  # working out how many to disregard, but need to count one -1 if present
                        # todo put back in if not working
                        # if len(triple) + leafCount >= self.groundsetSize:
                        #     return False, False
                    checkedSubsets = True

        # return first True because we got to this point without hitting a false
        # second true, because we did not hit early True ie subset
        return True, True


    def kTangle(self, k):

        def formatSideName(side):   ##### fix this to switch on "verbose"
            side = str(side)
            side = side.replace("frozenset(", "").replace(")", "")
            side = side.replace("'", "")
            return side

        def addSideAsSep(side, parent, sepNum):

            passes, toAdd = self.checkTangleAxioms(side)
            if passes:

                if toAdd:
                    child = parent.add_child(name=formatSideName(side))
                    self.nodeIndex+=1
                    newList = [smallSide for id, smallSide in enumerate(self.smallSidesList) if self.keepSeps[id]]
                    newList.append(side)
                    child.add_feature("smallSides", newList)
                    # self.prevBranches.append(child)
                    currentBranch = child
                else:
                    currentBranch = parent

                self.prevBranches.append(currentBranch)

                if sepNum == numkSeps-1:
                    self.foundTangle +=1
                    # todo editing to make tidy tree tidier - kludge. Fix later to code with sep ids and A/B maybe?
                    # child.name = "T{}{}".format(self.currentTangle, formatSideName(side))
                    currentBranch.name = "T{}-{}".format(self.currentTangle, currentBranch.name)
                    self.TangleLists[k].append(currentBranch)
                    self.currentTangle += 1

                if toAdd:
                    # if it passes and we need to add it, ie, it's not a subset, then we DO need to check the comp
                    # so precludesComp = False
                    return False
                else:
                    # if we don't add it, it's because it's a subset, therefore the complement is precluded
                    return True
            else:
                # ret val is precludesComp. If does not pass checks, DO want to check comp
                return False

        self.foundTangle = 0

        numkdefSmall = len(self.definitelySmall[k])
        numkSeps = len(self.separations[k]) + numkdefSmall

        # loop through and find the last order which had separations
        lastfound = False
        for lastorder in range(k-1, -1, -1):
            if len(self.TangleLists[lastorder]) > 0:
                self.prevBranches = self.TangleLists[lastorder]
                # print("Lastorder: {}".format(lastorder))
                lastfound = True
                break
        if not lastfound:
            exit("crack the sads, no prev tangles found")


        # self.prevBranches = self.TangleLists[k-1]
        print("Before building {}: Len of prevBranches: {}".format(k, len(self.prevBranches)))
        # self.separations[k] = sorted(self.separations[k], key=len, reverse=True)
        # Do the most uneven separations first, as they're likely to break first
        # note that since this list only contains the smallest side of each separation
        # the smallest small side means the most uneven separation


        for sepNum in range(numkSeps):
            currentBranches = self.prevBranches
            self.prevBranches = []
            for truncTangle in currentBranches:
                # todo !!!! changes here.
                # self.smallSidesStack = truncTangle.smallSides   ###### *****
                self.smallSidesList = truncTangle.smallSides   ###### *****
                if self.smallSidesList is None:
                    self.smallSidesList = []
                if sepNum < numkdefSmall:
                    #### No branching
                    addSideAsSep(self.definitelySmall[k][sepNum], truncTangle, sepNum)
                elif sepNum < numkSeps:
                    ### check both sides of the separation
                    # NOTE: Edited so that only the side with fewest elements is stored
                    # other must be calculated
                    # for side in self.separations[k][sepNum - numkdefSmall]:
                    #     addSideAsSep(side, truncTangle, sepNum)
                    side = self.separations[k][sepNum - numkdefSmall]
                    precludesComp = addSideAsSep(side, truncTangle, sepNum)
                    if not precludesComp:
                        complement = self.groundset - side
                        addSideAsSep(complement, truncTangle, sepNum)

        if self.foundTangle:
            self.finaliseAndPrint(k)
            return(True)
        else:
            return(False)


    def finaliseAndPrint(self, k):

        orientArray = pd.DataFrame()
        orientArrayNamed = pd.DataFrame()


        def printSingleTangle(tangle, tangID, tangOrder, verbose=False):

            tangOrients = []
            tangOrientsNamed = []

            if tangle.smallSides is None:
                tangle.smallSides = []
            ##########################
            for smallSide in tangle.smallSides:
                # if smallSide == {6,7,8}:
                #     global debugBreak
                #     debugBreak = True
                #     print("Found trouble sep in printSingleTangle")

                complement = self.groundset - smallSide

                if verbose:
                    print("---------------------")
                    print("{}: order {}".format(tangID, tangOrder))
                    if len(smallSide) <= len(complement):
                        print("{} / G\\*".format(set(smallSide)))
                    else:
                        print("G\\* / {}".format(complement))

                if len(smallSide) <= len(complement):
                    tangOrients.append("{} / G\\*".format(sorted(set(smallSide))))
                    nameList = [self.names[id] if id >= 0 else str(id) for id in smallSide]
                    tangOrientsNamed.append("{} / G\\*".format(sorted(nameList)))
                    # tangOrientsNamed.append("{} / G\\*".format(sorted(list(map(self.names.__getitem__, smallSide)))))
                else:
                    tangOrients.append("G\\* / {}".format(sorted(complement)))
                    nameList = [self.names[id] if id >= 0 else str(id) for id in complement]
                    tangOrientsNamed.append("{} / G\\*".format(sorted(nameList)))
                    # tangOrientsNamed.append("G\\* / {}".format(sorted(list(map(self.names.__getitem__, complement)))))

            tangHeader = "{}: ord={}".format(tangID, tangOrder)
            temparray = pd.DataFrame()
            temparrayNamed = pd.DataFrame()
            temparray[tangHeader] = tangOrients
            temparrayNamed[tangHeader] = tangOrientsNamed
            nonlocal orientArray
            nonlocal orientArrayNamed
            orientArray = pd.concat([orientArray, temparray], axis=1)
            orientArrayNamed = pd.concat([orientArrayNamed, temparrayNamed], axis=1)

        tangleCounter = 0
        for i in range(self.kmin, k+1):
            for tangle in self.TangleLists[i]:
                printSingleTangle(tangle, tangleCounter, i+1, verbose=False)
                tangleCounter+=1

        # I think this file gets overwritten at each order. Not ideal...
        OrientsOutfile = "{}/{}-Orientations.csv".\
        format(self.job['outputFolder'], self.job['outName'])
        orientArray.to_csv(OrientsOutfile)

        OrientsOutfileNamed = "{}/{}-OrientationsNamed.csv".\
        format(self.job['outputFolder'], self.job['outName'])
        orientArrayNamed.to_csv(OrientsOutfileNamed)
        # print("OUTPUT ORIENTATIONS")

        # with open('sepList.txt', 'w') as f:
        #     for i in [1,2,3]:
        #         f.write("sep order: {}\n".format(i))
        #         for item in [i]:
        #             f.write("{}\n".format(item))

        #####################
        doTreePrint = False
        if doTreePrint:
            try:
                self.printTangleTree(k)
            except Exception as error:
                print("There was an error trying to print the tangle tree:")
                print(error)



    def printTangleTree(self, k):
        outfile = "{}/{}-TangleTree-{}.png". \
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
            print("Render doesn't work correctly in Python 3.10. Use 3.9 or lower.")
            # print(rendError)

        tidyTree = self.TangleTree.copy()
        # print(tidyTree.get_ascii(show_internal=True))
        for treenode in tidyTree.traverse():
            if "T" not in treenode.name:
                treenode.delete(prevent_nondicotomic=False)
            else:
                dummy = 1
        ts.show_branch_length = False

        for node in tidyTree.iter_descendants():
            node.dist *= 4

        ts.branch_vertical_margin = 8
        ts.scale = 360

        tidyOutfile = "{}/{}-TidyTangleTree-{}.pdf". \
            format(self.job['outputFolder'], self.job['outName'], k)
        try:
            tidyTree.render(tidyOutfile, tree_style=ts)
        except Exception as rendError:
            print("Render doesn't work correctly in Python 3.10. Use 3.9 or lower.")
            # print(rendError)


