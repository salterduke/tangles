import math
import numpy as np
import scipy as sp
import csv
import pandas as pd
import collections as coll
import igraph as ig
import itertools as it
import tools
import ete3
import sys
import heapq
import copy


import BaseTangleSet as btang



class partialCut(object):
    def __init__(self, Gdir, Gdircopy, pcut=None, mincut=None):
        def unmergeVnames(mergedVids):
            newVids = set()
            # print("Merged Vids")
            # print(mergedVids)
            # todo unmerge, but check for lone vertices being assigned to the wrong side
            for vid in mergedVids:
                # print(vid)
                if len(Gdircopy.vs[vid]["name"]) > 0:
                    newVids.update(v.index for v in map(Gdir.vs.find, Gdircopy.vs[vid]["name"].split(",")))
            return newVids

        self.weight = mincut.value
        self.pcut = copy.deepcopy(pcut)

        self.mcut = [unmergeVnames(sideVs) for sideVs in mincut.partition] # todo fix for merges!

        self.cutEdges = frozenset(sorted(
            [tuple(sorted((edge["st"][0], edge["st"][1])))
             for edge in mincut.es]))
        # todo this isn't right yet!

    def __lt__(self, other):
        return self.weight < other.weight

    def __str__(self):
        # todo add more shit later?
        return "wt {} S:{}".format(self.weight, self.pcut["S"])

    def __repr__(self):
        # todo add more shit later?
        return "\nwt {} S:{}".format(self.weight, self.pcut["S"])

    # todo crack the shits if still all none?


class EdgeTangleSet(btang.TangleSet):
    def __init__(self, G, job, log):
        self.G = G
        self.groundset = set(self.G.vs["name"])
        self.groundsetSize = self.G.vcount()
        btang.TangleSet.__init__(self, job, log)
        self.cuts = set()
        self.Gdirected = self.G.as_directed()

        self.names = self.G.vs['name']


        self.lineGraph = self.G.linegraph()

        ####
        self.cutfinder = False
        self.doGH = False
        self.stupid = True


        if self.cutfinder:
            self.sepFilename = "{}/{}-SepList-CF.tsv". \
                format(job['outputFolder'], job['outName'])
        elif self.doGH:
            self.sepFilename = "{}/{}-SepList-GHU.tsv". \
                format(job['outputFolder'], job['outName'])
            self.GHTree = self.G.gomory_hu_tree()
            for v in self.GHTree.vs():
                v["name"] = self.names[v.index]
            # self.findAllMinCuts()
        else:
            self.sepFilename = "{}/{}-SepList.tsv".\
                    format(job['outputFolder'], job['outName'])


        # todo Check if need to add orientation as per git.
        text = "order\tcut\tside1\tside2\torientation\n"
        with open(self.sepFilename, 'w+') as the_file:
            the_file.write(text)



    def findNextOrderSeparations(self, k = None):
        if self.cutfinder:
            self.findNextOrderSeparationsCutFinder(k)
        elif self.doGH:
            self.findNextOrderSeparationsGH(k)
        else:
            self.findNextOrderSeparationsBrute(k)


    # todo whack a numba thingy on this guy? - probably a wrapper
    # todo https://stackoverflow.com/questions/41769100/how-do-i-use-numba-on-a-member-function-of-a-class
    def findNextOrderSeparationsCutFinder(self, k=None, maxdepth=4):
        def mergeVnames(names):
            return ",".join(names)

        def mergeVertices(pcut):
            minS = min(pcut["S"])
            minT = min(pcut["T"])

            self.Gdircopy = self.Gdirected.copy()
            # todo - can I think of a faster way of doing this?
            mapvector = [minS if v in pcut["S"] else (minT if v in pcut["T"] else v) for v in self.Gdircopy.vs.indices]

            # print("------------------------------- in merge")
            # print(pcut)
            # print(mapvector)


            self.Gdircopy.contract_vertices(mapvector, combine_attrs=dict(name=mergeVnames))
            # self.Gdircopy.contract_vertices(mapvector, combine_attrs="first")
            # self.Gdircopy.contract_vertices(mapvector)

            # print([(e.source, e.target) for e in self.Gdircopy.es])
            return minS, minT


        # todo lets def as an inner fn here, that way if I numba it, I can fairly easily chuck it out into a separate jit-able fn.
        # from Yeh, Li-Pu, Biing-Feng Wang, and Hsin-Hao Su. 2010. “Efficient Algorithms for the Problems of Enumerating Cuts by Non-Decreasing Weights.” Algorithmica. An International Journal in Computer Science 56 (3): 297–312.
        def basicPartition():
            for edge in self.Gdirected.es:
                edge["st"] = (self.Gdirected.vs[edge.source]["name"], self.Gdirected.vs[edge.target]["name"])
                # store original source target as attr, so okay when merge.

            B0 = [] # initialise the heap
            n = self.Gdirected.vcount()
            # self.Gdircopy = self.Gdirected.copy()
            self.vids = self.Gdirected.vs.indices
            # let s = node 0, (see Yeh and Wang) so start at 1
            self.U = range(1,n)  # like this so easier to change later
            S = set([0])
            # this is so we can do a Source-t cut,
            for vi in self.U:
                partcut = {
                    "S": S.copy(),
                    "T": set([vi])
                }
                # print("merging in basic")
                s, t = mergeVertices(partcut)
                mincut = self.Gdircopy.mincut(source=s, target=t)
                if len(mincut.partition) == 2:
                    heapq.heappush(B0, partialCut(self.Gdirected,self.Gdircopy,partcut,mincut))
                else:
                    # todo Is this okay? I think we don't want it on the heap if non-minimal.
                    print("More than 2 components: {}".format(mincut))
                    input("press any key to continue")
                S.add(vi) # todo I think I still add it to S though. I think.

            return(B0)

        def getNextPartCut(k):
            if len(self.partcutHeap) > 0 and self.partcutHeap[0].weight <= k:
                return heapq.heappop(self.partcutHeap)
            else:
                return None

        def resetGroupSourceSink(currentpcut):
            self.U = range(1,self.G.vcount())

            # print("merging in reset")
            mergeVertices(currentpcut)
            # todo Check this

        def extractMinPart(partial):
            resetGroupSourceSink(partial.pcut)
            SuT = partial.pcut["S"].union(partial.pcut["T"])
            self.U = [u for u in self.U if u not in SuT]
            n = self.G.vcount()
            # prevpcut = copy.deepcopy(partial.pcut)
            prevpcut = partial.pcut

            sidetoaddto = None

            for Ucounter in range(len(self.U)):
                for side in ["S", "T"]:
                    if (side == "S"):
                        otherside = "T"
                    else:  # == "T"
                        otherside = "S"

                    if minInPart({side:prevpcut[side] | set([self.U[Ucounter]]), otherside:prevpcut[otherside]}, partial.mcut):
                        sidetoaddto = side
                        # then continue

                    else:
                        # print("merging in extract min")

                        s, t, = mergeVertices(
                            {side: prevpcut[side] | set([self.U[Ucounter]]), otherside: prevpcut[otherside]})
                        # print(s,t)
                        # todo check edges correctly added by here
                        mincut = self.Gdircopy.mincut(source=s, target=t)
                        # print("----------------------------------------------")
                        # print("mincut partition")
                        # print(mincut.partition)
                        if len(mincut.partition) == 2:
                            # convert to list of sets for easier subset checking
                            if mincut.value <= self.kmax:
                                heapq.heappush(self.partcutHeap,
                                               partialCut(self.Gdirected,self.Gdircopy,{side:prevpcut[side] | set([self.U[Ucounter]]), otherside:prevpcut[otherside]}, mincut))
                                # print([(e.source, e.target) for e in self.Gdircopy.es])
                                #
                                # print("moocow")
                        else:
                            # todo Is this okay? I think we don't want it on the heap if non-minimal.
                            print("More than 2 components in extract min: {}".format(mincut))
                            input("Press any key to continue")
                if(sidetoaddto):
                    # print("side to add: ", sidetoaddto)
                    prevpcut[sidetoaddto].add(self.U[Ucounter])

        def minInPart(pcut, mcut):
            return ((pcut["S"].issubset(mcut[0]) and pcut["T"].issubset(mcut[1])) or
                    (pcut["S"].issubset(mcut[1]) and pcut["T"].issubset(mcut[0])))

        if k is None:  ### ie, first time running
            # todo check initialise heap.
            self.partcutHeap = basicPartition()
            self.TangleTree.add_feature("cutsets", [])
            # todo kmax for avoiding heap pushes.
            self.kmin = int(self.partcutHeap[0].weight)
            k = self.kmin
            self.kmax = k + maxdepth

        while (partcut := getNextPartCut(k)) != None:

            extractMinPart(partcut)
            self.addToSepList(partcut)
            # todo note that the vertex IDs *should* be the same for G and Gdirected,
            # todo - and this will break if they are not
        # print("End of fn")
        # sys.exit()


    def findNextOrderSeparationsGH(self, k=None):
        def cutIsSuperset(newCut):
            for cut in self.cuts:
                if newCut.issuperset(cut):
                    return True
            return False

        if k is None:  ### ie, first time running
            self.TangleTree.add_feature("cutsets", [])

            self.kmin = int(min(self.GHTree.es()["flow"]))
            self.kmax = int(max(self.GHTree.es()["flow"]))
            k = self.kmin

            # print("GHTree -------------------------------------------")
            # print(self.GHTree.es()["flow"])
            # print("kmin = {}".format(self.kmin))

        # print("********************************")
        # print("k: {}".format(k))

        testAllcuts = True

        for edge in self.GHTree.es.select(flow_eq=k):
            s = self.GHTree.vs[edge.source]["name"]
            t = self.GHTree.vs[edge.target]["name"]
            if testAllcuts:
                # print("About to run all_st_cuts. flow={}".format(k))
                cutnum=0
                for cut in self.Gdirected.all_st_cuts(s, t):
                    cutnum+=1
                    # print("cut {}".format(cutnum))
                    cutEdges = frozenset(sorted(
                        [tuple(sorted((self.G.vs[edge.source]["name"], self.G.vs[edge.target]["name"]))) for edge in
                         cut.es]))
                    self.addToSepListOld(cut.partition, cutEdges)

            else:
                # print("About to run all_st_MINcuts. flow={}".format(k))
                cutnum=0
                for cut in self.Gdirected.all_st_mincuts(s, t):
                    cutnum+=1
                    # print("cut {}".format(cutnum))
                    cutEdges = frozenset(sorted(
                        [tuple(sorted((self.G.vs[edge.source]["name"], self.G.vs[edge.target]["name"]))) for edge in
                         cut.es]))
                    self.addToSepListOld(cut.partition, cutEdges)

        # Adding singletons
        # todo not sure I need to do this, with doing ALL_st_cuts
        for vert in self.G.vs.select(_degree_eq=k):
            side1 = {vert.index}
            side2 = set(self.G.vs.indices) - side1
            cutEdges = frozenset(sorted([tuple(
                sorted((self.G.vs[self.G.es[edge].source]["name"], self.G.vs[self.G.es[edge].target]["name"]))) for
                                         edge in self.G.incident(vert)]))
            if not cutIsSuperset(cutEdges):
                self.addToSepListOld([side1, side2], cutEdges)


    def addToSepList(self, partial):
        def cutIsSuperset(newCut):
            for cut in self.cuts:
                if newCut.issuperset(cut):
                    return True
            return False


        def defSmallExists(side, k):
            for i in range(self.kmin, k+1):
                if side in self.definitelySmall[i]:
                    # print("Does this even happen with edge connectivity? Yes, apparently")
                    return True
            return False

        def ambigSepExists(sep, k):
            for i in range(self.kmin, k+1):
                if sep in self.separations[i]:
                    # print("Does this even happen with edge connectivity? Yes, apparently")
                    return True
            return False

        def printSepToFile(separation, cut, orientation):
            separation = sorted(separation, key=len)

            sideNodes = sorted([self.names[node] for node in separation[0]])
            complementNodes = sorted([self.names[node] for node in separation[1]])

            text = "{}\t{}\t{}\t{}\t{}\n".format(len(cut), sorted(cut), sideNodes, complementNodes, orientation)
            with open(self.sepFilename, 'a') as the_file:
                the_file.write(text)

        if cutIsSuperset(partial.cutEdges):
            return

        components = partial.mcut  # saves me refactoring.

        # smart: making sure the first side is the shortest
        components.sort(key=len)

        sideNodes = frozenset({self.names[node] for node in components[0]})
        complementNodes = frozenset({self.names[node] for node in components[1]})

        sideSubGraph = self.G.subgraph(components[0])
        complementSubGraph = self.G.subgraph(components[1])

        size = len(partial.cutEdges)

        ######## ******* do other checks for "easy" seps (do shit here)
        ### smart:  initial check for def small - add more checks here?
        # todo - check don't need existence checks with new Yeh and Wang algorithm?
        # todo - what about if both sides def small?
        # todo - add check for union of all def small sides? Only if *all* seps def orientable?
        makeTrunk = True
        if(makeTrunk):
            # if sideSubGraph.vcount() == sideSubGraph.ecount() + 1 or \
            # (all(deg <= 2 for deg in sideSubGraph.degree()) and size >= 2):
            if (all(deg <= 2 for deg in sideSubGraph.degree()) and size >= 2) or (all(deg <= 1 for deg in sideSubGraph.degree())):
                # if not defSmallExists(sideNodes, size):
                if True:
                    # print("adding small - Side")
                    self.definitelySmall[size].append(sideNodes)
                    orientation = 1
                else:
                    return
            # elif complementSubGraph.vcount() == complementSubGraph.ecount() + 1 or \
            # (all(deg <= 2 for deg in complementSubGraph.degree()) and size >= 2):
            elif (all(deg <= 2 for deg in complementSubGraph.degree()) and size >= 2) or (all(deg <= 1 for deg in complementSubGraph.degree())):
                # if not defSmallExists(complementNodes, size):
                if True:
                    # print("adding small - Complement")
                    self.definitelySmall[size].append(complementNodes)
                    orientation = 2
                else:
                    return

            else:
                ### smart - this was here before, outside if else
                separation = (sideNodes, complementNodes)
                # if not ambigSepExists(separation, size):
                if True:
                    self.separations[size].append(separation)
                    orientation = 3
                else:
                    return

        else:
            separation = (sideNodes, complementNodes)
            if not ambigSepExists(separation, size):
                self.separations[size].append(separation)
                orientation = 3
            else:
                return

        self.cuts.add(partial.cutEdges)
        printSepToFile(components, partial.cutEdges, orientation)

# #########################################################
    # this is for GH
    # todo check if brute is using this one or the other one.
    def addToSepListOld(self, components, cutEdges):
        def defSmallExists(side, k):
            for i in range(self.kmin, k+1):
                if side in self.definitelySmall[i]:
                    # print("Does this even happen with edge connectivity? Yes, apparently")
                    return True
            return False

        def ambigSepExists(sep, k):
            for i in range(self.kmin, k+1):
                if sep in self.separations[i]:
                    # print("Does this even happen with edge connectivity? Yes, apparently")
                    return True
            return False

        def printSepToFile(separation, cut, orientation):
            separation = sorted(separation, key=len)

            sideNodes = sorted([self.names[node] for node in separation[0]])
            complementNodes = sorted([self.names[node] for node in separation[1]])

            text = "{}\t{}\t{}\t{}\t{}\n".format(len(cut), sorted(cut), sideNodes, complementNodes, orientation)
            with open(self.sepFilename, 'a') as the_file:
                the_file.write(text)


        if len(components) > 2:
            print("Arrgh, too many components, crack the shits!")
            print(components)
            print("============")
            print(self.definitelySmall)
            exit()

        # add the condition in here to account for the added node in Gdirected for S,t cuts
        for side in [0,1]:
            try:
                components[side].remove(self.G.vcount())
            except:
                pass

            try:
                components[side].remove(self.G.vcount()+1)
            except:
                pass


        # smart: making sure the first side is the shortest
        components.sort(key=len)

        sideNodes = frozenset({self.names[node] for node in components[0]})
        complementNodes = frozenset({self.names[node] for node in components[1]})

        sideSubGraph = self.G.subgraph(components[0])
        complementSubGraph = self.G.subgraph(components[1])


        size = len(cutEdges)

        ######## ******* do other checks for "easy" seps (do shit here)
        ### smart:  initial check for def small - add more checks here?
        makeTrunk = True
        if(makeTrunk):
            # if sideSubGraph.vcount() == sideSubGraph.ecount() + 1 or \
            # (all(deg <= 2 for deg in sideSubGraph.degree()) and size >= 2):
            if (all(deg <= 2 for deg in sideSubGraph.degree()) and size >= 2):
                if not defSmallExists(sideNodes, size):
                    # print("adding small - Side")
                    self.definitelySmall[size].append(sideNodes)
                    orientation = 1
                else:
                    return
            # elif complementSubGraph.vcount() == complementSubGraph.ecount() + 1 or \
            # (all(deg <= 2 for deg in complementSubGraph.degree()) and size >= 2):
            elif (all(deg <= 2 for deg in complementSubGraph.degree()) and size >= 2):
                if not defSmallExists(complementNodes, size):
                    # print("adding small - Complement")
                    self.definitelySmall[size].append(complementNodes)
                    orientation = 2
                else:
                    return

            else:
                ### smart - this was here before, outside if else
                separation = (sideNodes, complementNodes)
                if not ambigSepExists(separation, size):
                    self.separations[size].append(separation)
                    orientation = 3
                else:
                    return

        else:
            separation = (sideNodes, complementNodes)
            if not ambigSepExists(separation, size):
                self.separations[size].append(separation)
                orientation = 3
            else:
                return

        # todo Note - should be able to remove existence checks with new algorithm.
        self.cuts.add(cutEdges)
        printSepToFile(components, cutEdges, orientation)




    def findNextOrderSeparationsBrute(self, k = None):
        if k is None:  ### ie, first time running
            self.TangleTree.add_feature("cutsets", [])

            minCut = self.G.mincut()
            # print("minCut")
            # print(minCut.value)

            self.kmin = int(minCut.value)
            k = self.kmin

        print("********************************")
        print("k: {}".format(k))
        for ksub in it.combinations(self.G.es, k):
            self.processCut(ksub)
        print("Brute")

    def processCut(self, edgeObjects):
        def cutIsSuperset(newCut):
            for cut in self.cuts:
                if newCut.issuperset(cut):
                    return True
            return False

        def printSepToFile(separation, cut):
            separation = sorted(separation, key=len)

            text = "{}\t{}\t{}\t{}\n".format(len(cut), sorted(cut), sorted(separation[0]), sorted(separation[1]))
            with open(self.sepFilename, 'a') as the_file:
                the_file.write(text)

        ### Kludge - if doing new way, components is not anything sensible
        # todo Fix later
        def addToSepList(components, cut):

            if self.stupid:
                sideNodes = frozenset({names[node] for node in components[0]})
                complementNodes = frozenset({names[node] for node in components[1]})
            else:
                sideNodes = frozenset(self.G.vs.select(clabel_eq=1)['name'])
                complementNodes = frozenset(self.G.vs.select(clabel_eq=2)['name'])


            ### remove this check since can't implement in new method.
            ### I don't *think* it's necessary?
            # if len(components) > 2:
            #     print("Arrgh, too many components, crack the shits!")
            #     exit()

            size = len(cut)

            ######## ******* do other checks for "easy" seps (do shit here)
            separation = frozenset([sideNodes, complementNodes])
            self.cuts.add(cut)
            ### self.cuts is only needed to check if new proposed cuts are supersets
            self.separations[size].append(separation)
            ######################################

            printSepToFile(separation, cut)

        def isCutSet(cut):
            # todo CHANGE GOES HERE
            # print(dir(self.G))

            dfs = self.lineGraph.dfsiter(0, mode=3, advanced=False)   # 3 = ALL.

            # dfs = self.G.dfsiter(0, mode=3, advanced=True)   # 3 = ALL.
            # reset labels to blank
            self.G.vs['clabel'] = 0
            # print("--------------------------")
            # print(cut)
            # print("--------------------------")

            ### keep running track of last vertex so edge oriented correct way
            u = None

            for evert in dfs:
                e = self.G.es[evert.index]

                ### Assign u and v appropriately
                if u == None:
                    ## first run
                    # print("Case None")
                    u = self.G.vs[e.source]
                    v = self.G.vs[e.target]
                    u['clabel'] = 1
                elif self.G.vs[e.source]['clabel']:
                    # print("Case Source")
                    u = self.G.vs[e.source]
                    v = self.G.vs[e.target]
                elif self.G.vs[e.target]['clabel']:
                    # print("Case Target")
                    u = self.G.vs[e.target]
                    v = self.G.vs[e.source]
                else:
                    print("Crack the poopahs, no label")
                    exit()

                tup = (v['name'],u['name'])
                if (v['name'],u['name']) in cut or (u['name'],v['name']) in cut:
                    # print("found")
                    # print(tup)
                    # print(cut)
                    newlabel = 3 - u['clabel'] ### switch between 1 and 2
                else:
                    newlabel = u['clabel']

                # print(tup, v['clabel'], u['clabel'], newlabel)

                if v['clabel'] and v['clabel'] != newlabel:
                    # print("Ruh roh - mismatch label, not cut set")
                    return(False)
                else:
                    # print("Not clabel")
                    v['clabel'] = newlabel
            # end for evert in dfs

            # for edge in self.G.es:
            #     if (self.G.vs[edge.source]["name"], self.G.vs[edge.target]["name"]) in cut:
            #         print("cut: ", self.G.vs[edge.source]["clabel"], self.G.vs[edge.target]["clabel"])
            #     else:
            #         print("Not cut: ", self.G.vs[edge.source]["clabel"], self.G.vs[edge.target]["clabel"])

            # for v in self.G.vs:
            #     print(v['name'], v["clabel"])


            ### if we've traversed the whole graph without mismatch, cutset:
            return(True)
        # end isCutSet()

        cut = frozenset(sorted([tuple(sorted((self.G.vs[edge.source]["name"], self.G.vs[edge.target]["name"]))) for edge in edgeObjects]))
        if not cutIsSuperset(cut):
            if self.stupid:
                Gcopy = self.G.copy()
                names = Gcopy.vs['name']
                Gcopy.delete_edges(edgeObjects)
                if len(Gcopy.components()) >= 2:
                    addToSepList(Gcopy.components(), cut)
            else:
                if isCutSet(cut):
                    addToSepList("Moocow", cut)
