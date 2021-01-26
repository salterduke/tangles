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

import BaseTangleSet as btang




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
        self.doGH = True
        if self.doGH:
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

        self.stupid = True


    # todo - I think this was just checking shit, and not needed.
    def findAllMinCuts(self):
        # Gdirected = self.G.as_directed()
        allMinCuts = pd.DataFrame(columns=list(self.G.vs["name"]),
                                  index=list(self.G.vs["name"]))
        treeNames = self.GHTree.vs["name"]

        moreIntersCount = 0

        for i in range(self.G.vcount()-1):
            for j in range(i+1,self.G.vcount()):
                # i = 3
                # j = 6

                s = self.G.vs[i]["name"]
                t = self.G.vs[j]["name"]

                # GHPath = self.GHTree.get_all_shortest_paths(s,t)
                GHPath = self.GHTree.get_all_shortest_paths(i,j)
                # print("---------------------------")
                # print(i, j, s, t)
                # print(GHPath)
                # GHPathVerts = [treeNames[v] for v in GHPath[0]]
                GHPathEdges = [(GHPath[0][vid], GHPath[0][vid+1])
                    for vid in range(len(GHPath[0])-1)]
                # print(GHPathEdges)
                GHAdj = self.GHTree.get_adjacency(attribute="flow")

                flowList = []
                for e in GHPathEdges:
                    # print("edge", e[0], e[1])
                    # print(treeNames[e[0]], treeNames[e[1]])
                    flow = GHAdj[e[0],e[1]]
                    # print(flow)
                    flowList.append(flow)

                minFlow = min(flowList)
                # todo Get min flow val - I would think there's an easier way

                ##### calc num cuts between i, j
                numCuts = len(self.Gdirected.all_st_mincuts(i, j))
                allMinCuts.loc[s,t] = numCuts

                numInterCuts = 0
                for cutId in range(len(flowList)):
                    if flowList[cutId] == minFlow:
                        e = GHPathEdges[cutId]
                        numInterCuts += len(self.Gdirected.all_st_mincuts(e[0],e[1]))
                        # print(self.Gdirected.all_st_mincuts(e[0],e[1]))

                # print("Total num cuts: {}".format(numCuts))
                # # print(self.Gdirected.all_st_mincuts(i, j))
                # print("Total Inter cuts: {}".format(numInterCuts))
                if numCuts > numInterCuts:
                    print("Shit!")
                    # cuts1 = self.Gdirected.all_st_mincuts(3,4)
                    # cuts2 = self.Gdirected.all_st_mincuts(10,4)
                    # print([c.cut for c in cuts1])
                    # print([c.cut for c in cuts2])
                    exit()
                elif numCuts < numInterCuts:
                    # print("More inters")
                    moreIntersCount+=1


        #### ** remove later!!!!!
        print("{} cuts had more inters".format(moreIntersCount))
        # exit()
        allMinCuts.to_csv("allMinCuts.csv")



    def findNextOrderSeparations(self, k = None):
        if self.doGH:
            self.findNextOrderSeparationsGH(k)
        else:
            self.findNextOrderSeparationsBrute(k)

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

        print("********************************")
        print("k: {}".format(k))

        testAllcuts = True

        for edge in self.GHTree.es.select(flow_eq=k):
            s = self.GHTree.vs[edge.source]["name"]
            t = self.GHTree.vs[edge.target]["name"]
            if testAllcuts:
                for cut in self.Gdirected.all_st_cuts(s, t):
                    cutEdges = frozenset(sorted(
                        [tuple(sorted((self.G.vs[edge.source]["name"], self.G.vs[edge.target]["name"]))) for edge in
                         cut.es]))
                    self.addToSepList(cut.partition, cutEdges)

            else:
                for cut in self.Gdirected.all_st_mincuts(s, t):
                    cutEdges = frozenset(sorted(
                        [tuple(sorted((self.G.vs[edge.source]["name"], self.G.vs[edge.target]["name"]))) for edge in
                         cut.es]))
                    self.addToSepList(cut.partition, cutEdges)

        # Adding singletons
        for vert in self.G.vs.select(_degree_eq=k):
            side1 = {vert.index}
            side2 = set(self.G.vs.indices) - side1
            cutEdges = frozenset(sorted([tuple(
                sorted((self.G.vs[self.G.es[edge].source]["name"], self.G.vs[self.G.es[edge].target]["name"]))) for
                                         edge in self.G.incident(vert)]))
            if not cutIsSuperset(cutEdges):
                self.addToSepList([side1, side2], cutEdges)

    # this is for GH
    # todo check if brute is using this one or the other one.
    def addToSepList(self, components, cutEdges):
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

        # smart: making sure the first side is the shortest
        components.sort(key=len)

        sideNodes = frozenset({self.names[node] for node in components[0]})
        complementNodes = frozenset({self.names[node] for node in components[1]})

        sideSubGraph = self.G.subgraph(components[0])
        complementSubGraph = self.G.subgraph(components[1])

        if len(components) > 2:
            print("Arrgh, too many components, crack the shits!")
            print(components)
            print("============")
            print(self.definitelySmall)
            exit()

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
