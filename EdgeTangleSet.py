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
    for cut in currCuts:
        if newCut.issuperset(cut):
            return True
    return False

class HO_cutfinder(object):
    def __init__(self, Gdir):
        # modified init
        pass

    def findmin(self):
        pass

    def relabel(self):
        pass

    def selectSink(self):
        pass


def HO_mincut(Gdir):
    cutfinder = HO_cutfinder(Gdir)




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

            if True:
            # todo - this is the check that speeds it up lots, but doesn't seem valid.
            # todo - are there alternate checks that are valid?
            # if not extCutIsSuperset(currCuts, newPartial.cutEdges):
                newHeapCuts.append(newPartial)
                # note that partialCut takes care of the vids re the adjusted graph
    return newHeapCuts


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

        self.cutEdges = frozenset(
            [tuple(sorted((edge["st"][0], edge["st"][1])))
             for edge in mincut.es])
        # note to self. frozenset doesn't maintain order. No need to sort here - possibly need to sort later, before output.
        dummy=1


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

            self.partcutHeap = [partcut for partcut in
                                pool.map(functools.partial(externalBasicPartitionBranch, tangset=self),
                                         range(len(self.U))) if partcut.weight <= self.kmax]
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

                # todo - check if need currCuts,
                results = pool.map(functools.partial(externalExtractMinPart, Gdir=self.Gdirected, kmax=self.kmax, currCuts = self.cuts), partcutList)
                origSize = len(self.partcutHeap)
                pcutCount = 0
                for pcut in [item for subresults in results for item in subresults]:  # todo make sure returns work okay
                    pcutCount+=1

                    if pcut.weight <= self.kmax:
                        heapq.heappush(self.partcutHeap, pcut)
                    else:
                        print("this got through")
                sizediff = len(self.partcutHeap) - origSize

                print("{} partcuts calculated {}, added {} more, newsize {}".format(len(partcutList), pcutCount, sizediff, len(self.partcutHeap)))

# this is the current working version - I think!
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
