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
        for k in range(self.kmin, self.kmax+1):
            # note k is SEP order, so +1 for tang order
            tangNumbers = (k+1, len(self.TangleLists[k]))
            countList.append(tangNumbers)
        return(countList)

    def findNextOrderSeparations(self):
        print("This must be overridden for each child")
        raise NotImplementedError

    #### delete depth later, when we're sure this works. ******
    def findAllTangles(self, depth=4, sepsOnly = False):
        print("Finding All Tangles")

        timings = []

        ### so that the first tangle tree starts at the root.
        self.TangleTree.add_feature("smallSides", [])

        print("-------------------------------------------")
        self.log.tick("kTangle min Find seps")
        self.findNextOrderSeparations(None, depth)
        if self.kmin is None:
            print("Crack the shits, kmin not set")
            exit()
        timings.append(self.log.tock())
        self.kmax = self.kmin + depth

        self.TangleLists[self.kmin - 1] = [self.TangleTree]
        self.prevBranches = [self.TangleTree]

        if not sepsOnly:
            # ### see Evrendilek for heuristic?
            self.log.tick("kTangle{} (min) Build Tangle".format(self.kmin))
            if not self.kTangle(self.kmin):
                print("No tangles exist")
                self.log.tock()
                return(timings)
            else:
                self.log.tock()

        # # # ### find all tangles at each k
        for k in range(self.kmin+1, self.kmax+1):
            print("-------------------------------------------")
            self.log.tick("kTangle{} Find seps".format(k))
            self.findNextOrderSeparations(k)
            timings.append(self.log.tock())
            if not sepsOnly:
                self.log.tick("kTangle{} Build Tangle".format(k))
                if not self.kTangle(k):
                    print("No tangles at k: {}".format(k))
                    self.kmax = k-1  # todo: make sure not off-by-one
                    self.log.tock()
                    self.log.end()
                    return(timings)
                self.log.tock()
        self.log.end()
        return(timings)


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
                    return True, False  # todo Really not sure about this!!!
                elif newSep.issuperset(side1):
                    self.keepSeps[0] = 0
                double1 = side1 | newSep
                if len(double1) >= self.groundsetSize:
                    return False, False
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
                    if len(double1) >= self.groundsetSize:
                        return False, False

                    # todo problem is enumerate starts id2 at 0!!!!!!!!
                    for id2 in range(id1 + 1, len(sepsSoFar)):
                        side2 = sepsSoFar[id2]
                    # for id2, side2 in enumerate(sepsSoFar[id1 + 1:]):
                        if not checkedSubsets:
                            if newSep.issubset(side2):
                                return True, False
                            elif newSep.issuperset(side2):
                                self.keepSeps[id2] = 0
                        triple = side2 | double1
                        if len(triple) >= self.groundsetSize:
                            return False, False
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
                    currentBranch.name = "T{}".format(self.currentTangle)
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
        # numkSeps = len(self.separations[k])


        # self.prevBranches = self.TangleLists[k-1]

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

        self.prevBranches = self.TangleLists[k]

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
            # todo !!!!! change to enum over set
            # for i in range(numSeps):
            for smallSide in tangle.smallSides:
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
                    tangOrientsNamed.append("{} / G\\*".format(sorted(list(map(self.names.__getitem__, smallSide)))))
                else:
                    tangOrients.append("G\\* / {}".format(sorted(complement)))
                    tangOrientsNamed.append("G\\* / {}".format(sorted(list(map(self.names.__getitem__, complement)))))

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

        OrientsOutfile = "{}/{}-Orientations.csv".\
        format(self.job['outputFolder'], self.job['outName'])
        orientArray.to_csv(OrientsOutfile)

        OrientsOutfileNamed = "{}/{}-OrientationsNamed.csv".\
        format(self.job['outputFolder'], self.job['outName'])
        orientArrayNamed.to_csv(OrientsOutfileNamed)
        # print("OUTPUT ORIENTATIONS")

        with open('sepList.txt', 'w') as f:
            for i in [1,2,3]:
                f.write("sep order: {}\n".format(i))
                for item in [i]:
                    f.write("{}\n".format(item))

        #####################
        doTreePrint = False
        if doTreePrint:
            self.printTangleTree(k)



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
            print(rendError)

        tidyTree = self.TangleTree.copy()
        # print(tidyTree.get_ascii(show_internal=True))
        for treenode in tidyTree.traverse():
            if "T" not in treenode.name:
                treenode.delete(prevent_nondicotomic=False)
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
            print(rendError)


