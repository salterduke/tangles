import heapq
import functools
import multiprocessing
import sys
import numpy as np
import igraph as ig
import collections as coll
import itertools as iter
import BaseTangleSet as btang
import Modules.tools as tools
from Modules.tools import mergeVnames

def extractComponents(mcut, vcount):
    # probably there's a quicker way of doing this, but just get it working for now.
    mcutList = list('{0:b}'.format(mcut).zfill(vcount)[::-1])
    comps = [set(), set()]
    for i in range(len(mcutList)):
        comps[int(mcutList[i])].add(i)
    return (comps)



def extCutIsSuperset(currCuts, newCut):
    for cut in currCuts:
        if newCut.issuperset(cut):
            return True
    return False



def externalExtractMinPart(partcut, Gdir, kmax):
    longpcut = extractComponents(partcut.pcut["binrep"], partcut.pcut["binlen"])

    SuT = longpcut[0].union(longpcut[1])
    U = [u for u in range(1, Gdir.vcount()) if u not in SuT]
    newHeapCuts = []

    components = extractComponents(partcut.mcut, Gdir.vcount())

    for uid in range(len(U)):

        newnode = U[uid]
        sidetoaddto = 1 if newnode in components[0] else 0
        # note for this bit, we're adding the new node to the side that *doesn't* match the mcut

        newpcut = [
            longpcut[0].union(components[0].intersection(U[0:uid])),
            longpcut[1].union(components[1].intersection(U[0:uid]))
        ]
        newpcut[sidetoaddto].add(newnode)

        Gdircopy, s, t = externalMergeVertices(Gdir, newpcut)
        mincut = Gdircopy.mincut(source=s, target=t)

        if mincut.value <= kmax:
            pcutlen = uid + len(SuT) + 1
            newPartial = partialCut(Gdir, Gdircopy, mincut, newpcut, pcutlen)  # todo again, do I need a + 1?

            newHeapCuts.append(newPartial)
            # note that partialCut takes care of the vids re the adjusted graph
    return newHeapCuts


def externalMergeVertices(Gdir, pcut):
    minS = min(pcut[0])
    minT = min(pcut[1])

    Gdircopy = Gdir.copy()
    # todo - can I think of a faster way of doing this?
    mapvector = [minS if v in pcut[0] else (minT if v in pcut[1] else v) for v in Gdircopy.vs.indices]

    Gdircopy.contract_vertices(mapvector, combine_attrs=dict(name=mergeVnames))
    return Gdircopy, minS, minT


# NOTE can probably just pass the few bits of tangset we need, but this isn't the slow bit, so eh...
def externalBasicPartitionBranch(uid, tangset):
    newpcut = [
        set(tangset.U[0:uid]).union([0]),
        set([tangset.U[uid]])
    ]
    Gdircopy, s, t = externalMergeVertices(tangset.Gdirected, newpcut)
    mincut = Gdircopy.mincut(source=s, target=t)
    dummy = 1
    if len(mincut.partition) == 2:
        return partialCut(tangset.Gdirected, Gdircopy, mincut, newpcut, uid+2) # todo note, 2 for: node 0, and 0 indexing
        # note that partialCut takes care of the vids re the adjusted graph
    else:
        # todo Is this okay? I think we don't want it on the heap if non-minimal.
        print("More than 2 components: {}".format(mincut))
        input("press any key to continue")


class partialCut(object):
    def __init__(self, Gdir, Gdircopy, mincut=None, longpcut=None, pcutlen=0):
        def getMcutShort(mincut):
            component1 = unmergeVnames(mincut.partition[1])  # only need [1] because only 1s add to binary val
            # sum(c << i for i, c in enumerate(mincut.membership))
            return(sum(1 << i for i in component1))

        def unmergeVnames(mergedVids):

            newVids = set()
            # todo unmerge, but check for lone vertices being assigned to the wrong side
            for vid in mergedVids:
                # print(vid)
                if len(Gdircopy.vs[vid]["name"]) > 0:
                    try:
                        newVids.update(v.index for v in map(Gdir.vs.find, Gdircopy.vs[vid]["name"].split(",")))
                    except BaseException as err:
                        print("error here in unmergeVnames. ")

            return newVids

        self.weight = mincut.value
        self.pcut = {
            "binrep": sum(1 << i for i in longpcut[1]),
            "binlen": pcutlen
        }

        self.mcut = getMcutShort(mincut)
        self.cutEdges = frozenset([edge.index for edge in mincut.es])

    def __lt__(self, other):
        return self.weight < other.weight

    def __str__(self):
        # todo add more shit later?
        return "wt {} binrep:{}, binlen: {}".format(self.weight, self.pcut["binrep"], self.pcut["binlen"])

    def __repr__(self):
        # todo add more shit later?
        return "\nwt {} binrep:{}, binlen: {}".format(self.weight, self.pcut["binrep"], self.pcut["binlen"])

class EdgeTangleSet(btang.TangleSet):
    def __init__(self, G, job, log):

        print("Initialising VY Edge Tangle Set")

        # Takes every "twig" hanging off the edges of the main graph, and condenses it to a single node as a "stub"
        # since a twig is always "small", don't need to deal with all parts of the twig.
        self.G = tools.pruneToStubs(G)

        self.leafExtent = 0
        self.leaves = set(self.G.vs.select(_degree_eq=1).indices)

        self.groundset = set(self.G.vs.indices)
        self.groundsetSize = self.G.vcount()
        btang.TangleSet.__init__(self, job, log)
        self.cuts = set()
        self.Gdirected = self.G.as_directed()

        self.names = self.G.vs['name']

        # ig.plot(self.G)
        #
        self.sepFilename = "{}/{}-SepList-CF.tsv". \
            format(job['outputFolder'], job['outName'])

        text = "order\tcut\tside1\tside2\torientation\n"
        with open(self.sepFilename, 'w+') as the_file:
            the_file.write(text)

    # todo whack a numba thingy on this guy? - probably a wrapper
    # todo https://stackoverflow.com/questions/41769100/how-do-i-use-numba-on-a-member-function-of-a-class
    def findNextOrderSeparations(self, k=None, maxdepth=4):

        # from Yeh, Li-Pu, Biing-Feng Wang, and Hsin-Hao Su. 2010. “Efficient Algorithms for the Problems of Enumerating Cuts by Non-Decreasing Weights.” Algorithmica. An International Journal in Computer Science 56 (3): 297–312.
        def basicPartition(pool):

            n = self.Gdirected.vcount()
            self.vids = self.Gdirected.vs.indices
            # let s = node 0, (see Yeh and Wang) so start at 1
            self.U = range(1,n)  # like this so easier to change later

            self.partcutHeap = [partcut for partcut in
                                pool.map(functools.partial(externalBasicPartitionBranch, tangset=self),
                                         range(len(self.U))) if partcut.weight <= self.kmax]
            heapq.heapify(self.partcutHeap)

        self.numkSeps = 0

        with multiprocessing.Pool() as pool:
            if k is None:  ### ie, first time running
                # slight waste of time, but saves memory by being able to avoid putting things on heap
                self.kmin = int(self.Gdirected.mincut().value)
                k = self.kmin
                self.kmax = k + maxdepth - 1
                basicPartition(pool)
                # self.TangleTree.add_feature("cutsets", [])

            # todo remove when done debugging
            # externalExtractMinPart(
            #     partialCut(S=np.array([1, 1, 1, 1, 1, 0, 0, 0, 0], dtype=bool),
            #                T=np.array([0, 0, 0, 0, 0, 0, 1, 0, 0], dtype=bool),
            #                Tstar=np.array([0, 0, 0, 0, 0, 1, 1, 1, 1], dtype=bool),
            #                weight=1.0),
            #     self.Gdirected, self.kmax
            # )
            # dummy = 1
            # exit()

            while len(self.partcutHeap) > 0 and self.partcutHeap[0].weight <= k:
                # I know this looks stupid, but partcutHeap gets modified by this loop
                # and we want to repeat until no more relevant partcuts in heap
                partcutList = []
                while len(self.partcutHeap) > 0 and self.partcutHeap[0].weight <= k:
                    newpartcut = heapq.heappop(self.partcutHeap)
                    partcutList.append(newpartcut)
                    self.addToSepList(newpartcut)

                results = pool.map(functools.partial(externalExtractMinPart, Gdir=self.Gdirected, kmax=self.kmax), partcutList)
                for partcut in [item for subresults in results for item in subresults]:  # todo make sure returns work okay
                    heapq.heappush(self.partcutHeap, partcut)

        return self.numkSeps

    def addToSepList(self, partial):
        def convertLeaves(newcomp):
            compLeaves = self.leaves & newcomp

            leafNumber = -1
            if len(newcomp) > 1 and len(compLeaves) > 0:
                for leaf in compLeaves:
                    if self.G.neighbors(leaf)[0] not in newcomp:
                        # ie, it's a separate lone leaf, not attached
                        # newcomp.remove(leaf)
                        # newcomp.add(-1)
                        newcomp = newcomp | {leafNumber}
                        newcomp = newcomp - {leaf}
                        leafNumber-=1
                self.leafExtent = min(self.leafExtent, leafNumber)
            return newcomp

        def sideIsDefSmall(side, sep_k):
            if (len(side) == 1) or (self.G.maxdegree(side) <= 2 and sep_k >= 2):
                return True
            # else:
            #     return False

            if (len(side) == 2) and (self.G.maxdegree(side) <= sep_k):
                # ie, each of two v's makes the small side of some sep,
                # and the union of two small sides is small if the k of each <= k of the union
                return True

            # todo this might be extendable to maxdeg <= sep_k?
            if (self.G.maxdegree(side) <= sep_k) and (sep_k >= 3):
                subG = self.G.induced_subgraph(side)
                comps = subG.clusters()
                if len(comps) == 1 and comps.subgraph(0).is_tree(mode="all"):
                    return True
                elif len(comps) > 1:
                    isTree = [1 if comp.is_tree(mode = "all") else 0 for comp in comps.subgraphs() ]
                    if all(isTree):
                        return True
                    elif sum(isTree) == len(comps) - 1:
                        notTree = comps.subgraph(isTree.index(0))
                        notTreeNames = notTree.vs["name"]
                        notTree.delete_vertices(notTree.vs.select(_degree_eq=1))
                        if notTree.maxdegree() <= 2:
                            # if this higher order sep is small this way, the subcomponent sep must be small
                            # therefore terminate any tangles that contradict this
                            notTreeIDs = set(self.G.vs.select(name_in=notTreeNames).indices)
                            comp = self.groundset - notTreeIDs
                            # todo removed this as seems to be error?
                            # try:
                            #     self.prevBranches = [branch for branch in self.prevBranches if comp not in branch.smallSides]
                            #     for branch in self.prevBranches:
                            #         branch.smallSides.remove(notTreeIDs)
                            #         # is not necessary, as is subset of this larger side.
                            # except:
                            #     # if sepsOnly, prevBranches is not defined, so error.
                            #     # but the above is only necessary when building tangles,
                            #     # so can ignore
                            #     pass
                            return True


        def addDefSmall(newcomp, newsize):

            toAdd = True
            for size in range(self.kmin, newsize+1):
                for comp in self.definitelySmall[size]:
                    if len(newcomp) < len(comp) and newcomp.issubset(comp):
                        toAdd = False
                        break
                    elif len(comp) < len(newcomp) and comp.issubset(newcomp):
                        self.definitelySmall[size].remove(comp)
                        # todo I know this isn't the most efficient. Not sure how to improve
                if not toAdd:
                    break
                    # first break on breaks inner loop.
                    # Don't need to continue checking at all.

            if toAdd:
                self.definitelySmall[newsize].append(newcomp)
            # todo, note that all seps are written to file, even if not tested.

            # todo I think this needs to not happen. Might be okay to do just for the same size?
            # for size in range(self.kmin, newsize+1):
            #     self.definitelySmall[size] = [comp for comp in self.definitelySmall[size]\
            #                                   if not newcomp.issuperset(comp)]
            # todo: possible version?
            # self.definitelySmall[newsize] = [comp for comp in self.definitelySmall[newsize]\
            #                               if not newcomp.issuperset(comp)]
            # self.definitelySmall[newsize].append(newcomp)

        def printSepToFile(components, cut, orientation):
            sideNodes = sorted([self.names[node] for node in components[0]])
            complementNodes = sorted([self.names[node] for node in components[1]])

            cutLong = []
            for eid in cut:
                edge = self.Gdirected.es[eid]
                eString = "('{}', '{}')".format(self.Gdirected.vs[edge.source]["name"], self.Gdirected.vs[edge.target]["name"])
                cutLong.append(eString)

            text = "{}\t{}\t{}\t{}\t{}\n".format(len(cut), sorted(cutLong), sideNodes, complementNodes, orientation)
            text = text.replace('\"', '')
            with open(self.sepFilename, 'a') as the_file:
                the_file.write(text)

        self.numkSeps+=1 # note, adding this here, as not all seps are added to defSmall so can't just get list len after

        components = extractComponents(partial.mcut, self.Gdirected.vcount())
        components = sorted(components, key=len)

        # size = len(partial.cutEdges)
        sep_k = int(partial.weight)

        ######## ******* do other checks for "easy" seps (do shit here)
        ### smart:  initial check for def small - add more checks here?
        # todo - what about if both sides def small?
        if sideIsDefSmall(components[0], sep_k):
            addDefSmall(convertLeaves(components[0]), sep_k)
            orientation = 1
        elif sideIsDefSmall(components[1], sep_k):
            addDefSmall(convertLeaves(components[1]), sep_k)
            orientation = 2
        else:
            # Note: edited so only adding the shortest side.
            self.separations[sep_k].append(convertLeaves(components[0]))
            orientation = 3

        # need to do this because we need to check for superset-ness
        # self.cuts.add(partial.cutEdges)
        # todo consider if can remove cutEdges?
        if "constructed" not in self.sepFilename:
            printSepToFile(components, partial.cutEdges, orientation)
