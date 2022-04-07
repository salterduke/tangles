import heapq
import functools
import multiprocessing
import sys
import numpy as np
import igraph as ig
import collections as coll
import itertools as iter


import BaseTangleSet as btang

def extractComponents(mcut, vcount):
    # probably there's a quicker way of doing this, but just get it working for now.
    mcutList = list('{0:b}'.format(mcut).zfill(vcount)[::-1])
    comps = [set(), set()]
    for i in range(len(mcutList)):
        comps[int(mcutList[i])].add(i)
    return (comps)


def mergeVnames(names):
    return ",".join(names)


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

            if True:
            # todo - this is the check that speeds it up lots, but doesn't seem valid.
            # todo - are there alternate checks that are valid?
            # if not extCutIsSuperset(currCuts, newPartial.cutEdges):
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
    def __init__(self, pcut, mask, mcut, weight):
        # expects bool arrays, converts to ints for storage
        self.pcut = self.encode(pcut)
        self.mask = self.encode(mask)
        self.mcut = self.encode(mcut)
        self.weight = weight


    def encode(self, A):
        return(sum(int(val) << idx for idx, val in enumerate(A)))

    def decode(self):
        pass

class partialCut_Old(object):
    def __init__(self, Gdir, Gdircopy, mincut=None, longpcut=None):
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
                    newVids.update(v.index for v in map(Gdir.vs.find, Gdircopy.vs[vid]["name"].split(",")))
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

class HaoOrlin():
    #### Partitions the partial cut P(s, \empty), and returns all elements <= kmax
    #### Assumes G is directed
    # todo add check for directed
    def __init__(self, H):
        self.H = H
        self.H.es["flow"] = 0
        if not self.H.is_weighted():
            self.H.es["weight"] = 1
        self.adj = np.array(self.H.get_adjacency(attribute = "weight").data)
        self.N = self.H.vcount() # for convenience
        self.pcutList = []

    def initFor_s(self, s):
        # s is the *index* of the vertex (at this stage, assuming it works)
        self.Dset = []
        J = self.H.successors(s)
        # todo this could be more efficient with get_eids and a comprehension, but let's just get it working first, hey?
        for j in J:
            eid = self.H.get_eid(s, j)
            self.H.es[eid]["flow"] = self.H.es[eid]["weight"]

        # self.Dset.append({s})
        self.Dset.append(1 << s)
        self.Dmax = 0

        # each vertex is represented by a single bit, with all the bits combined stored as an integer
        self.W = np.ones(self.N, dtype = bool) # using 1 and 0 for true and false in subsequent lines
        self.W[s] = 0

        self.S = np.zeros(self.N, dtype = bool)
        self.S[s] = 1

        self.t = 1  # todo write checks s != t      # todo also change to T for extract min?
        self.d = np.ones(self.H.vcount(), dtype=np.int16) # todo consider dtype
        self.d[self.t] = 0
        self.dCount = coll.Counter({1:(self.N-1), 0:1}) # ie, there are N-1 nodes at dist 1, 1 node at dist 0
        # todo make sure to update dCount every time d is updated!

        # todo *** maybe take out if not debugging?
        self.H.es["label"] = self.H.es["flow"]
        self.H.es["curved"] = 0


    # if findall == False, find only the min, else find all <= kmax
    def HOfindCuts(self, s, kmax, findall = True):
        self.initFor_s(s)
        self.kmax = kmax

        while sum(self.S) < self.N:
            while self.existsActiveNode():
                i = self.activeNode
                if self.existsAdmissableEdge(i):
                    self.pushFlow()
                else:
                    self.relabel(i)
            self.updateCutList()
            dummy = 1
            self.selectSink()

    def updateCutList(self):
        # pcut is S, T = {t} (for basic partition)
        # completion (mcut) is D, W (S, T)
        weight = self.adj[~self.W,:][:,self.W].sum()
        if weight <= self.kmax:
            # todo edit later when T is more than one vertex
            self.T = np.zeros(self.N, dtype=bool)
            self.T[self.t] = 1

            self.pcutList.append(
                partialCut(pcut = self.T,
                           mask = self.S | self.T,
                           mcut = self.W,
                           weight = weight)
            )

    def existsActiveNode(self):
        # for v in iter.chain(self.activeNode, self.W - {self.t, self.activeNode}):
        # for v in iter.chain(self.activeNode, [idx for idx, val in enumerate(bin(self.W)[:1:-1]) if val == '1' and idx not in [self.t, self.activeNode]]):
        for v in [idx for idx, val in enumerate(self.W) if
                  val == 1 and idx != self.t]:
            # todo update when T can be more than one node
            self.excess_i = sum(self.H.es[self.H.incident(v, mode="in")]["flow"]) - \
            sum(self.H.es[self.H.incident(v, mode="out")]["flow"])
            if self.excess_i > 0:
                self.activeNode = v
                return True
        self.activeNode = -1
        return False


    def existsAdmissableEdge(self, i):
        J = self.H.successors(i)
        # for j in self.W.intersection(J):
        for j in [idx for idx, val in enumerate(self.W) if val == 1 and idx in J]:
            if self.d[i] == self.d[j] + 1:
                self.ij = self.H.get_eid(i, j)
                self.ji = self.H.get_eid(j, i)
                self.r_ij = self.H.es[self.ij]["weight"] - self.H.es[self.ij]["flow"] + self.H.es[self.ji]
                if  self.r_ij > 0:
                    return True
        return False

    def pushFlow(self):
        delta = min(self.r_ij, self.excess_i)

        # reverse any backflow before pushing forward flow
        if self.H.es[self.ji]["flow"] > 0:
            delta_1 = min(self.H.es[self.ji]["flow"], delta)
            self.H.es[self.ji]["flow"] -= delta_1
            delta -= delta_1
        self.H.es[self.ij]["flow"] += delta


    def relabel(self, i):
        di = self.d[i]
        if self.dCount[di] == 1:
            # R = {j for j in self.W if self.d[j] >= di}
            # R = sum(1 << idx for idx, val in enumerate(self.W) if val == '1' and self.d[idx] >= di)
            R = sum(1 << idx for idx, val in enumerate(self.W) if val == 1 and self.d[idx] >= di)
            self.Dmax+=1
            self.Dset.append(R)
            self.W = self.W - R
        else:
            J = self.H.successors(i)
            j_in_W  = [idx for idx, val in enumerate(self.W) if val == 1 and idx in J]
            # j_in_W = self.W.intersection(arcs)
            if len(j_in_W) == 0:
                R = 1 << i
                self.Dmax += 1
                self.Dset.append(R)
                self.W = self.W - R
            else:
                jDists = {self.d[j] for j in j_in_W if (self.H.es[self.H.get_eid(i, j)]["weight"] - self.H.es[self.H.get_eid(i, j)]["flow"] + self.H.es[self.H.get_eid(j, i)]["flow"] ) > 0}
                self.d[i] = min(jDists) + 1
                # todo update dcount!!!!
        # this func done. Does it work? Who can say?

    def selectSink(self):
        pass
        # todo: up to here #####


class EdgeTangleSet(btang.TangleSet):
    def __init__(self, G, job, log):

        # todo check this change
        shell = np.array(G.coreness())
        # remove 1-shell completely:
        # G.delete_vertices(np.where(shell == 1)[0])

        # leave "stubs"
        del_Vs = []
        for vid in np.where(shell == 1)[0]:
            to_del = True
            for nb in G.neighbors(vid):
                if shell[nb] > 1:
                    to_del = False
            if to_del:
                del_Vs.append(vid)

        G.delete_vertices(del_Vs)

        # todo change ends here
        self.G = G


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
        # def basicPartitionPreYWS(pool):
        #
        #     n = self.Gdirected.vcount()
        #     self.vids = self.Gdirected.vs.indices
        #     # let s = node 0, (see Yeh and Wang) so start at 1
        #     self.U = range(1,n)  # like this so easier to change later?
        #
        #     self.partcutHeap = [partcut for partcut in
        #                         pool.map(functools.partial(externalBasicPartitionBranch, tangset=self),
        #                                  range(len(self.U))) if partcut.weight <= self.kmax]
        #     heapq.heapify(self.partcutHeap)



        def basicPartition(pool):

            n = self.Gdirected.vcount()
            self.vids = self.Gdirected.vs.indices
            # let s = node 0, (see Yeh and Wang) so start at 1
            HO = HaoOrlin(self.Gdirected)
            self.partcutHeap = HO.HOfindCuts(s=0, kmax=self.kmax)

            heapq.heapify(self.partcutHeap)


        with multiprocessing.Pool() as pool:
            if k is None:  ### ie, first time running
                # slight waste of time, but saves memory by being able to avoid putting things on heap
                self.kmin = int(self.Gdirected.mincut().value)
                k = self.kmin
                self.kmax = k + maxdepth
                basicPartition(pool)
                # self.TangleTree.add_feature("cutsets", [])

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

            # do the singletons that are in the middle of the graph, so that the cut removing them is actually a composition of cuts.

    def addToSepList(self, partial):
        def addDefSmall(newcomp, newsize):

            for size in range(self.kmin, newsize):
                self.definitelySmall[size] = [comp for comp in self.definitelySmall[size]\
                                              if not newcomp.issuperset(comp)]
            self.definitelySmall[newsize].append(newcomp)

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


        components = extractComponents(partial.mcut, self.Gdirected.vcount())
        components = sorted(components, key=len)

        # size = len(partial.cutEdges)
        size = int(partial.weight)

        ######## ******* do other checks for "easy" seps (do shit here)
        ### smart:  initial check for def small - add more checks here?
        # todo - what about if both sides def small?
        if (len(components[0])==1) or (max(self.G.degree(components[0])) <= 2 and size >= 2) or (max(self.G.degree(components[0])) <= 1):
            addDefSmall(components[0], size)
            # self.definitelySmall[size].append(components[0])
            orientation = 1
        elif (len(components[1])==1) or (max(self.G.degree(components[1])) <= 2 and size >= 2) or (max(self.G.degree(components[1])) <= 1):
            addDefSmall(components[1], size)
            # self.definitelySmall[size].append(components[1])
            orientation = 2
        else:
            # Note: edited so only adding the shortest side.
            self.separations[size].append(components[0])
            orientation = 3

        # need to do this because we need to check for superset-ness
        # self.cuts.add(partial.cutEdges)
        # todo consider if can remove cutEdges?
        printSepToFile(components, partial.cutEdges, orientation)
