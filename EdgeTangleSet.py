import itertools as it
import sys
import heapq
import copy
import functools
import multiprocessing

import BaseTangleSet as btang


def mergeVnames(names):
    return ",".join(names)


def extCutIsSuperset(currCuts, newCut):
    if (set([('YGL016W', 'YLR335W'), ('YNL189W', 'YOR098C'), ('YLR347C', 'YOR098C'), ('YLR293C', 'YOR160W')]).issubset(
            newCut) or
            set([('YGL016W', 'YLR335W'), ('YOR098C', 'YNL189W'), ('YLR347C', 'YOR098C'),
                 ('YLR293C', 'YOR160W')]).issubset(newCut)):
        print(newCut)
    for cut in currCuts:
        if newCut.issuperset(cut):
            return True
    return False


def externalExtractMinPart(partcut, Gdir, kmax, currCuts):
    SuT = partcut.pcut["S"].union(partcut.pcut["T"])
    U = [u for u in range(1, Gdir.vcount()) if u not in SuT]
    newHeapCuts = []
    for uid in range(len(U)):

        newnode = U[uid]
        sidetoaddto = "T" if newnode in partcut.mcut[0] else "S"
        # note for this bit, we're adding the new node to the side that *doesn't* match the mcut
        newpartcut = {
            "S": partcut.pcut["S"].union(partcut.mcut[0].intersection(U[0:uid])),
            "T": partcut.pcut["T"].union(partcut.mcut[1].intersection(U[0:uid]))
        }
        newpartcut[sidetoaddto].add(newnode)

        Gdircopy, s, t = externalMergeVertices(Gdir, newpartcut)
        mincut = Gdircopy.mincut(source=s, target=t)

        if mincut.value <= kmax:
            newPartial = partialCut(Gdir, Gdircopy, newpartcut, mincut)
            # if True:
            if not extCutIsSuperset(currCuts, newPartial.cutEdges):
                newHeapCuts.append(newPartial)
                # note that partialCut takes care of the vids re the adjusted graph
    return newHeapCuts


def externalMinPartBranch(uid, Gdir, partcut):
    if 0 not in partcut.mcut[0]:
        print(partcut.mcut)
        sys.exit("0 not in side 0 in mcut!")

    newnode = tangset.U[uid]
    sidetoaddto = "T" if newnode in partcut.mcut[0] else "S"
    # note for this bit, we're adding the new node to the side that *doesn't* match the mcut
    newpartcut = {
        "S": partcut.pcut["S"].union(partcut.mcut[0].intersection(tangset.U[0:uid])),
        "T": partcut.pcut["T"].union(partcut.mcut[1].intersection(tangset.U[0:uid]))
    }
    newpartcut[sidetoaddto].add(newnode)

    Gdircopy, s, t = externalMergeVertices(tangset, newpartcut)
    mincut = Gdircopy.mincut(source=s, target=t)
    if len(mincut.partition) == 2:
        if mincut.value <= tangset.kmax:
            return partialCut(tangset.Gdirected, Gdircopy, newpartcut, mincut)
            # note that partialCut takes care of the vids re the adjusted graph
        else:
            return None     # todo check if this is best way to handle
    else:
        # todo Is this okay? I think we don't want it on the heap if non-minimal.
        print("More than 2 components: {}".format(mincut))
        input("press any key to continue")


def externalMergeVertices(Gdir, pcut):
    minS = min(pcut["S"])
    minT = min(pcut["T"])

    Gdircopy = Gdir.copy()
    # todo - can I think of a faster way of doing this?
    mapvector = [minS if v in pcut["S"] else (minT if v in pcut["T"] else v) for v in Gdircopy.vs.indices]

    Gdircopy.contract_vertices(mapvector, combine_attrs=dict(name=mergeVnames))
    return Gdircopy, minS, minT


def externalBasicPartitionBranch(uid, tangset):
    partcut = {
        "S": set(tangset.U[0:uid]).union([0]),
        "T": set([tangset.U[uid]])
    }
    # print(multiprocessing.current_process().name) # todo - if this works, switch from tangset to Gdir here as well
    Gdircopy, s, t = externalMergeVertices(tangset.Gdirected, partcut)
    mincut = Gdircopy.mincut(source=s, target=t)
    if len(mincut.partition) == 2:
        return partialCut(tangset.Gdirected, Gdircopy, partcut, mincut)
        # note that partialCut takes care of the vids re the adjusted graph
    else:
        # todo Is this okay? I think we don't want it on the heap if non-minimal.
        print("More than 2 components: {}".format(mincut))
        input("press any key to continue")


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

        self.mcut = [unmergeVnames(sideVs) for sideVs in mincut.partition]

        if 0 not in self.mcut[0]:
            self.mcut[0], self.mcut[1] = self.mcut[1], self.mcut[0]
            print("O not in side 0. New mcut:")
            print(self.mcut)

        self.cutEdges = frozenset(sorted(
            [tuple(sorted((edge["st"][0], edge["st"][1])))
             for edge in mincut.es]))

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

        self.lineGraph = None  # only used if stupid, easier to check for None and then initialise than check if exists

        ####
        self.cutfinder = True
        self.doGH = True
        self.stupid = False


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



    def findNextOrderSeparations(self, k = None, maxdepth=4):
        if self.cutfinder:
            self.findNextOrderSeparationsCutFinder(k, maxdepth)
        elif self.doGH:
            self.findNextOrderSeparationsGH(k)
        else:
            self.findNextOrderSeparationsBrute(k)


    # todo whack a numba thingy on this guy? - probably a wrapper
    # todo https://stackoverflow.com/questions/41769100/how-do-i-use-numba-on-a-member-function-of-a-class
    def findNextOrderSeparationsCutFinder(self, k=None, maxdepth=4):

        # from Yeh, Li-Pu, Biing-Feng Wang, and Hsin-Hao Su. 2010. “Efficient Algorithms for the Problems of Enumerating Cuts by Non-Decreasing Weights.” Algorithmica. An International Journal in Computer Science 56 (3): 297–312.
        def basicPartition(pool):
            print("Basic Partition")
            for edge in self.Gdirected.es:
                edge["st"] = (self.Gdirected.vs[edge.source]["name"], self.Gdirected.vs[edge.target]["name"])
                # store original source target as attr, so okay when merge.

            n = self.Gdirected.vcount()
            self.vids = self.Gdirected.vs.indices
            # let s = node 0, (see Yeh and Wang) so start at 1
            self.U = range(1,n)  # like this so easier to change later

            self.partcutHeap = [partcut for partcut in pool.map(functools.partial(externalBasicPartitionBranch, tangset=self), range(len(self.U))) if partcut.weight <= self.kmax]
            print("Size of partcutHeap after basic: {}".format(len(self.partcutHeap)))
            heapq.heapify(self.partcutHeap)

        with multiprocessing.Pool() as pool:
            if k is None:  ### ie, first time running
                # ---------------------------------
                # slight waste of time, but saves memory by being able to avoid putting things on heap
                self.kmin = int(self.Gdirected.mincut().value)
                k = self.kmin
                self.kmax = k + maxdepth
                # ---------------------------------
                basicPartition(pool)
                self.TangleTree.add_feature("cutsets", [])

                # todo note that the vertex IDs *should* be the same for G and Gdirected,
                # todo - and this will break if they are not


            while len(self.partcutHeap) > 0 and self.partcutHeap[0].weight <= k:
                # I know this looks stupid, but partcutHeap gets modified by this loop
                # and we want to repeat until no more relevant partcuts in heap
                partcutList = []
                while len(self.partcutHeap) > 0 and self.partcutHeap[0].weight <= k:
                    newpartcut = heapq.heappop(self.partcutHeap)
                    partcutList.append(newpartcut)
                    self.addToSepList(newpartcut)

                results = pool.map(functools.partial(externalExtractMinPart, Gdir=self.Gdirected, kmax=self.kmax, currCuts = self.cuts), partcutList)
                origSize = len(self.partcutHeap)
                pcutCount = 0
                for pcut in [item for subresults in results for item in subresults]:  # todo make sure returns work okay
                    pcutCount+=1
                    if pcut.weight <= self.kmax:
                        heapq.heappush(self.partcutHeap, pcut)
                sizediff = len(self.partcutHeap) - origSize

                print("{} partcuts calculated {}, added {} more, newsize {}".format(len(partcutList), pcutCount, sizediff, len(self.partcutHeap)))


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
        # components.sort(key=len)  # fucks some stuff up # todo not important, but see if fixable

        sideNodes = frozenset({self.names[node] for node in components[0]})
        complementNodes = frozenset({self.names[node] for node in components[1]})


        size = len(partial.cutEdges)

        ######## ******* do other checks for "easy" seps (do shit here)
        ### smart:  initial check for def small - add more checks here?
        # todo - what about if both sides def small?
        # todo - add check for union of all def small sides? Only if *all* seps def orientable?
        makeTrunk = True
        if(makeTrunk):
            if (len(components[0])==1) or (max(self.G.degree(components[0])) <= 2 and size >= 2) or (max(self.G.degree(components[0])) <= 1):
                self.definitelySmall[size].append(sideNodes)
                orientation = 1
            elif (len(components[1])==1) or (max(self.G.degree(components[1])) <= 2 and size >= 2) or (max(self.G.degree(components[1])) <= 1):
                self.definitelySmall[size].append(complementNodes)
                orientation = 2
            else:
                separation = (sideNodes, complementNodes)
                self.separations[size].append(separation)
                orientation = 3

        else:
            separation = (sideNodes, complementNodes)
            self.separations[size].append(separation)
            orientation = 3

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

            if not self.lineGraph:
                self.lineGraph = self.G.linegraph()
            dfs = self.lineGraph.dfsiter(0, mode=3, advanced=False)   # 3 = ALL.

            self.G.vs['clabel'] = 0

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

                if (v['name'],u['name']) in cut or (u['name'],v['name']) in cut:
                    newlabel = 3 - u['clabel'] ### switch between 1 and 2
                else:
                    newlabel = u['clabel']


                if v['clabel'] and v['clabel'] != newlabel:
                    return(False)
                else:
                    v['clabel'] = newlabel
            # end for evert in dfs

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
