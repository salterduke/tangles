import heapq
import functools
import multiprocessing
import sys
import numpy as np
import igraph as ig
import collections as coll
import itertools as iter
import BaseTangleSet as btang


def mergeVnames(names):
    return ",".join(names)

def externalExtractMinPart(partcut, Gdir, kmax):
    # new partcut is pcut, mcut, weight, mask

    newPartcutList = []

    HO = HaoOrlin(Gdir)
    print("------------------- After creation")
    print(Gdir.vs.indices)
    print(HO.G.vs.indices)
    HO.initForPartial(partcut)
    print("------------------- before side 0")
    print(Gdir.vs.indices)
    print(HO.G.vs.indices)
    if HO.getSide(sideID = 0, partial = partcut):  # sets .H to be just the source side of partcut
        HO.HOfindCuts(kmax=kmax)
        newPartcutList = HO.getInducedCuts(sideID = 0, oldpart = partcut)

    print("------------------- before side 1")
    print(HO.G.vs.indices)
    if HO.getSide(sideID = 1, partial = partcut):  # sets .H to be just the target side of partcut
        HO.HOfindCuts(kmax=kmax)
        newPartcutList = newPartcutList + HO.getInducedCuts(sideID = 1, oldpart = partcut)

    print("------------------- after side 1")
    print(HO.G.vs.indices)

    return newPartcutList



def externalExtractMinPartOld(partcut, Gdir, kmax):
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


def externalMergeVerticesOld(Gdir, pcut):
    minS = min(pcut[0])
    minT = min(pcut[1])

    Gdircopy = Gdir.copy()
    # todo - can I think of a faster way of doing this?
    mapvector = [minS if v in pcut[0] else (minT if v in pcut[1] else v) for v in Gdircopy.vs.indices]

    Gdircopy.contract_vertices(mapvector, combine_attrs=dict(name=mergeVnames))
    return Gdircopy, minS, minT


# NOTE can probably just pass the few bits of tangset we need, but this isn't the slow bit, so eh...
# NOTE not being used for YWS
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

    def __init__(self, S, T, Tstar, weight):
        # expects bool arrays, converts to ints for storage
        # note that int bits are right-to-left, bools are left-to-right
        self.pcut = self.encode(T)
        self.mask = self.encode(S | T)
        self.mcut = self.encode(Tstar)
        self.weight = weight
        self.bitlen = len(S) # todo error checking etc.
        # todo
        # since bitlen changes when v's contracted, have to store for each separate object rather than globally
        # this takes up extra space. At some point, if too much mem used, look into converting to array of ints?
        # Will this reduce total overhead?

    def encode(self, Arr):
        # note that int bits are right-to-left, bools are left-to-right
        return(sum(int(val) << idx for idx, val in enumerate(Arr)))

    def decode(self, bits):
        # note that int bits are right-to-left, bools are left-to-right
        # -1 reverses the order, :1 leaves out the "0b". bitlen saves time by allowing pre-allocation.
        return (np.fromiter((int(i) for i in bin(bits)[:1:-1].ljust(self.bitlen, "0")), dtype=bool, count=self.bitlen))

    def getS(self, asString = False):
        if asString:
            return(bin(self.mask ^ self.pcut)[:1:-1].ljust(self.bitlen, "0"))
            # todo check this works
        else:
            return self.decode(self.mask) ^ self.decode(self.pcut)
            # todo check if ^ works inside decode brackets
            # todo check for mask false but pcut true? Should be impossible?

    def getT(self, asString = False):
        if asString:
            return(bin(self.pcut)[:1:-1].ljust(self.bitlen, "0"))
        else:
            return self.decode(self.pcut)

    def getSstar(self, asString = False):
        # not sure if this will ever get used
        if asString:
            # xor with all 1s
            return bin(self.mcut ^ (2 ** (self.bitlen)) - 1)[:1:-1].ljust(self.bitlen, "0")
            # todo check this later.
        else:
            return ~self.decode(self.mcut)

    def getTstar(self, asString = False):
        if asString:
            return bin(self.mcut)[:1:-1].ljust(self.bitlen, "0")
        else:
            return self.decode(self.mcut)

    def mapvector(self):
        # these take up more memory, so this only done *after* popped out of heap
        self.S = self.getS()
        self.T = self.getT()

        # print("in initForPartial")
        # print(self)
        #
        minS = 0
        if not self.S[0]:
            print("What the shit? node 0 not in getS")
            exit()

        minT = self.T.argmax()
        if not self.T[minT]:
            print("What the shit? no nodes in getT")
            exit()

        # todo - can I think of a faster way of doing this?
        # note - this assumes that G.vs.indices == range(bitlen)
        return [minS if self.S[v] else (minT if self.T[v] else v) for v in range(self.bitlen)], minS, minT

    def matchesPartition(self, partn):
        # partn is list of lists.

        self.Tstar = self.getTstar()

        if 0 in partn[1]:
            # so ordered as getSstar, getTstar
            partn.reverse()

        for sideID in (0, 1):
            for vid in partn[sideID]:
                if self.Tstar[vid] != sideID:
                    print("Oops, partition doesn't match mincut")
                    print("Partition", partn)
                    print(self)
                    return False

        return True

    def __lt__(self, other):
        return self.weight < other.weight

    def __str__(self):
        # todo add more shit later?
        return "wt {} S:{}, T: {}, Tstar: {}".format(self.weight, self.getS(asString=True), self.getT(asString=True), self.getTstar(asString=True)) + " wt {} pcut:{}, mask: {}, mcut: {}".format(self.weight, self.pcut, self.mask, self.mcut)


    def __repr__(self):
        # todo add more shit later?
        return "\nwt {} S:{}, T: {}, Tstar: {}".format(self.weight, self.getS(asString=True), self.getT(asString=True), self.getTstar(asString=True)) + " wt {} pcut:{}, mask: {}, mcut: {}".format(self.weight, self.pcut, self.mask, self.mcut)


    def components(self):
        # returns [getSstar, getTstar] as sets of vertex indices
        if not hasattr(self, "compList"):
            # todo add hasattr check to other "get" functions
            self.compList = [set(), set()]
            Tstar = self.getTstar()
            for idx, elem in enumerate(Tstar):
                self.compList[elem].add(idx)

        return self.compList


class HaoOrlin():
    #### Partitions the partial cut P(s, \empty), and returns all elements <= kmax
    #### Assumes G is directed
    # todo add check for directed
    def __init__(self, G):
        self.G = G.copy()   # I think this is necessary so we don't fuck up the original?
        if not self.G.is_weighted():
            self.G.es["weight"] = 1
        self.Gadj = np.array(self.G.get_adjacency(attribute = "weight").data)
        self.H = self.G     # G is the full graph, H is the working version. for basic partition, they're the same. later will be only getS* or getT*

    def initFor_s(self, s):
        # s is the *index* of the vertex (at this stage, assuming it works)
        self.s = s

        self.H.es["flow"] = 0

        if self.H.ecount() > 0:
            self.adj = np.array(self.H.get_adjacency(attribute = "weight").data)
        else:
            self.adj = np.array(self.H.get_adjacency(attribute=None).data)

        self.N = self.H.vcount() # for convenience
        self.partcutList = []
        self.flowlist = []


        self.Dset = []
        J = self.H.successors(s)
        # todo this could be more efficient with get_eids and a comprehension, but let's just get it working first, hey?
        for j in J:
            eid = self.H.get_eid(s, j)
            self.H.es[eid]["flow"] = self.H.es[eid]["weight"]

        # each vertex is represented by a single bit, with all the bits combined stored as an integer
        self.W = np.ones(self.N, dtype = bool) # using 1 and 0 for true and false in subsequent lines
        self.W[s] = 0

        self.S = np.zeros(self.N, dtype = bool)
        self.S[s] = 1

        # self.Dset.append(1 << s)
        self.Dset.append(self.S.copy())
        self.Dmax = 0


        # todo sort out what happens when N = 1, also check that N is updated properly.
        self.t = 1  # todo write checks s != t      # todo also change to getT for extract min?
        self.d = np.ones(self.H.vcount(), dtype=np.int16) # todo consider dtype
        print("--------------")
        print(self.H.vcount(), self.t, self.N)
        print("--------------**")
        self.d[self.t] = 0
        self.dCount = coll.Counter({1:(self.N-1), 0:1}) # ie, there are N-1 nodes at dist 1, 1 node at dist 0
        # todo make sure to update dCount every time d is updated!

        # todo *** maybe take out if not debugging?
        self.H.es["label"] = self.H.es["flow"]
        self.H.es["curved"] = 0

    def initForPartial(self, partcut):

        mapvector, minS, minT = partcut.mapvector()


        self.G.vs["name"] = [str(id) for id in self.G.vs.indices]
        self.G.contract_vertices(mapvector, combine_attrs=dict(name=mergeVnames)) # todo consider replacing with a lambda
        self.G.simplify(combine_edges="sum")
        self.G.delete_vertices([v.index for v in self.G.vs if v["name"] == ''] )


        self.mergedNodes = [minS, minT]

        # get flow and check it matches up
        flow = self.G.maxflow(minS, minT, capacity="weight")
        if flow.value != partcut.weight:
            print("What the shit? flow and cut weights don't match up")
            exit()

        if not partcut.matchesPartition(flow.partition):
            exit()


        for ij_idx in self.G.es.indices:
            i = self.G.es[ij_idx].source
            j = self.G.es[ij_idx].target
            ji_idx = self.G.get_eid(j, i)
            self.G.es[ij_idx]["weight"] = self.G.es[ij_idx]["weight"] - flow.flow[ij_idx] + flow.flow[ji_idx]
            # update weights to get residual graph
        self.partcut = partcut

    def getSide(self, sideID, partial):
        sourceName = self.G.vs[self.mergedNodes[sideID]]["name"]
        print("---------------------")
        print(sourceName, sideID)
        print(partial.components())
        print(self.G.vs.indices)
        print(self.G.vs["name"])
        if self.G.vcount() < 5:
            print("moocow")
        print("---------------------**")


        self.H = self.G.induced_subgraph(self.G.vs.select(
            lambda v: {int(vid) for vid in v["name"].split(",")}.issubset(partial.components()[sideID])).indices)
        if self.H.vcount() <= 1:
            return False

        self.s = self.H.vs.find(sourceName).index   # in case indices get all stuffed about, can't just use minS, minT?

        if sideID == 1:
            # ie, side T - need to get graph transpose

            A = np.array(self.H.get_adjacency(attribute = "weight").data)
            newH = ig.Graph.Weighted_Adjacency(A.T.tolist(), attr="weight") # transpose
            newH.vs["name"] = self.H.vs["name"] # copy names
            self.H = newH
            # I don't think there's an easier way to do it?

        self.initFor_s(self.s)
        return True

    def getInducedCuts(self, sideID, oldpart):
        newList = []
        for sidepart in self.partcutList:
            Sside = self.expandLabels(sidepart.getS())
            Tside = self.expandLabels(sidepart.getT())
            Sstarside = self.expandLabels(sidepart.getSstar())
            Tstarside = self.expandLabels(sidepart.getTstar())
            # Note, can't do this within partial cut, as that doesn't store the graphs

            if sideID == 0:
                # ie, side S
                Snew = Sside
                Tnew = oldpart.getT() | Tside
                Sstarnew = Sstarside
                Tstarnew = oldpart.getTstar() | Tstarside
            elif sideID == 1:
                # ie, side T - note that these are based on edge-reversed version of graph
                Snew = oldpart.getSstar() | Tside
                Tnew = Sside
                Sstarnew = oldpart.getSstar() | Tstarside
                Tstarnew = Sstarside
            else:
                print("Crack the shits, invalid sideID, ", sideID)
                exit()

            # todo check that Sstar + Tstar is everything, and disjoint
            weight = self.Gadj[Sstarnew, :][:, Tstarnew].sum()
            # note using full adjacency, not just part

            newList.append(partialCut(
                S = Snew,
                T = Tnew,
                Tstar = Tstarnew,
                weight = weight
            ))
        return newList

    def expandLabels(self, valArray):
        labels = [longlabel.split(",") for longlabel in
                   self.H.vs[(idx for idx, val in enumerate(valArray) if val == 1)]["name"] ]
        vids = [int(label) for label in iter.chain.from_iterable(labels)]
        newArray = np.zeros(self.G.vcount(), dtype=bool)
        print("----------------------------------")
        print(self.G.vcount())
        print(labels)
        print(vids)
        newArray[vids] = 1
        return newArray
        # could probably combine into one command, but that seems hard to read!


    def HOfindCuts(self, kmax):
        self.kmax = kmax

        doneFlag = False

        # Not checking for getS == N here, as done in selectSink, and no point doing sum() twice
        while not doneFlag:
            while self.existsActiveNode():
                i = self.activeNode
                if self.existsAdmissableEdge(i):
                    self.pushMaxFlow(i, self.j)
                else:
                    self.relabel(i)
            self.updateCutList()
            doneFlag = self.selectSink()


    def updateCutList(self):
        # pcut is getS, getT = {t} (for basic partition)
        # completion (mcut) is D, W (getS, getT)
        weight = self.adj[~self.W,:][:,self.W].sum()
        if weight <= self.kmax:
            # todo edit later when getT is more than one vertex
            self.T = np.zeros(self.N, dtype=bool)
            self.T[self.t] = 1

            self.partcutList.append(
                partialCut(S = self.S,
                           T = self.T,
                           Tstar = self.W,
                           weight = weight)
            )
            self.flowlist.append(self.H.es["flow"]) # todo remove when done debugging

    def existsActiveNode(self):
        # for v in iter.chain(self.activeNode, self.W - {self.t, self.activeNode}):
        # for v in iter.chain(self.activeNode, [idx for idx, val in enumerate(bin(self.W)[:1:-1]) if val == '1' and idx not in [self.t, self.activeNode]]):
        for v in [idx for idx, val in enumerate(self.W) if
                  val == 1 and idx != self.t]:
            # todo update when getT can be more than one node
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
                ij = self.H.get_eid(i, j)
                ji = self.H.get_eid(j, i)
                r_ij = self.H.es[ij]["weight"] - self.H.es[ij]["flow"] + self.H.es[ji]["flow"]
                # this ends up getting recalculated when doing pushflow, but it might be safer that way
                # # since pushflow doesn't always follow this method, so making r_ij global and relying
                # on it seems dubious
                # todo consider alternatives
                if  r_ij > 0:
                    self.j = j
                    return True
                # todo consider returning -1 if not exists, j otherwise? Better practice but fiddly
        return False

    def pushMaxFlow(self, a, b):

        ab = self.H.get_eid(a, b)
        ba = self.H.get_eid(b, a)

        r_ab = self.H.es[ab]["weight"] - self.H.es[ab]["flow"] + self.H.es[ba]["flow"]

        excess_a = sum(self.H.es[self.H.incident(a, mode="in")]["flow"]) - \
                        sum(self.H.es[self.H.incident(a, mode="out")]["flow"])

        delta = min(r_ab, excess_a)
        # if r_ab > excess_a and a == self.t:
        #     # ie, if doing selectSink
        #     print("Not enough excess")
        #     exit()
        # todo see if this is okay!!

        # reverse any backflow before pushing forward flow
        if self.H.es[ba]["flow"] > 0:
            delta_1 = min(self.H.es[ba]["flow"], delta)
            self.H.es[ba]["flow"] -= delta_1
            delta -= delta_1
        self.H.es[ab]["flow"] += delta
        # todo do I need an error code?

    def relabel(self, i):
        di = self.d[i]
        if self.dCount[di] == 1:
            # R = {j for j in self.W if self.d[j] >= di}
            # R = sum(1 << idx for idx, val in enumerate(self.W) if val == '1' and self.d[idx] >= di)
            R = np.fromiter((1 if val == 1 and self.d[idx] >= di else 0 for idx, val in enumerate(self.W) ),
                            dtype=bool, count=self.N)
            self.Dmax+=1
            self.Dset.append(R)
            self.W = self.W ^ R
        else:
            J = self.H.successors(i)
            j_in_W  = [idx for idx, val in enumerate(self.W) if val == 1 and idx in J]
            # j_in_W = self.W.intersection(arcs)
            if len(j_in_W) == 0:
                # R = 1 << i
                R = np.zeros(self.N, dtype=bool)
                R[i] = 1
                self.Dmax += 1
                self.Dset.append(R)
                self.W = self.W ^ R
            else:
                oldDist = self.d[i]
                minDist = min({self.d[j] for j in j_in_W if (self.H.es[self.H.get_eid(i, j)]["weight"] - self.H.es[self.H.get_eid(i, j)]["flow"] + self.H.es[self.H.get_eid(j, i)]["flow"] ) > 0})
                self.d[i] = minDist + 1
                self.dCount[oldDist]-=1
                self.dCount[minDist + 1]+=1

        # this func done. Does it work? Who can say?

    def selectSink(self):
        self.W[self.t] = 0
        self.S[self.t] = 1
        self.Dset[0][self.t] = 1
        if sum(self.S) == self.N:
            return True # ie, done = True

        neighbours = self.H.successors(self.t)
        for k in neighbours:
            if self.S[k] == 0:  # ie, k != getS
                self.pushMaxFlow(self.t, k)

        if not np.any(self.W):
            self.W = self.Dset[self.Dmax]
            self.Dmax-=1

        # select j in W such that d[j] is minimum
        # todo is there a nice way of doing this?
        jlist = list((idx for idx, val in enumerate(self.W) if val == 1))
        self.t = jlist[np.argmin(self.d[jlist])]
        dummy = 1
        return False

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

        def basicPartition():

            # self.vids = self.Gdirected.vs.indices
            # let s = node 0, (see Yeh and Wang) so start at 1
            HO = HaoOrlin(self.Gdirected)
            HO.initFor_s(0)
            HO.HOfindCuts(kmax=self.kmax)
            self.partcutHeap = HO.partcutList
            heapq.heapify(self.partcutHeap)


        with multiprocessing.Pool() as pool:
            if k is None:  ### ie, first time running
                # slight waste of time, but saves memory by being able to avoid putting things on heap
                self.kmin = int(self.Gdirected.mincut(capacity = "weight").value)
                k = self.kmin
                self.kmax = k + maxdepth
                basicPartition()  # not done with pool (at least at this stage)
                # self.TangleTree.add_feature("cutsets", [])

            while len(self.partcutHeap) > 0 and self.partcutHeap[0].weight <= k:
                # I know this looks stupid, but partcutHeap gets modified by this loop
                # and we want to repeat until no more relevant partcuts in heap
                partcutList = []
                while len(self.partcutHeap) > 0 and self.partcutHeap[0].weight <= k:
                    newpartcut = heapq.heappop(self.partcutHeap)
                    partcutList.append(newpartcut)
                    self.addToSepList(newpartcut)

                # todo remove after debugging this part!!!!
                # res = externalExtractMinPart(partcutList[0], Gdir=self.Gdirected, kmax=self.kmax)
                # print(res)
                # print("moocow")
                # exit()

                results = pool.map(functools.partial(externalExtractMinPart, Gdir=self.Gdirected, kmax=self.kmax), partcutList)
                for partcut in [item for subresults in results for item in subresults]:  # todo make sure returns work okay
                    heapq.heappush(self.partcutHeap, partcut)

            # do the singletons that are in the middle of the graph, so that the cut removing them is actually a composition of cuts.

    def addToSepList(self, partcut):
        def addDefSmall(newcomp, newsize):

            # todo I think this needs to not happen. Might be okay to do just for the same size?
            # for size in range(self.kmin, newsize):
            #     self.definitelySmall[size] = [comp for comp in self.definitelySmall[size]\
            #                                   if not newcomp.issuperset(comp)]
            # todo: possible version?
            # self.definitelySmall[newsize] = [comp for comp in self.definitelySmall[newsize]\
            #                               if not newcomp.issuperset(comp)]
            self.definitelySmall[newsize].append(newcomp)

        def printSepToFile(cutweight, components, orientation):
            sideNodes = sorted([self.names[node] for node in components[0]])
            complementNodes = sorted([self.names[node] for node in components[1]])

            text = "{}\t{}\t{}\t{}\n".format(cutweight, sideNodes, complementNodes, orientation)
            text = text.replace('\"', '')
            with open(self.sepFilename, 'a') as the_file:
                the_file.write(text)

        components = partcut.components()
        components = sorted(components, key=len)

        size = int(partcut.weight)

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
        # todo seems okay, but unsure what the issue mentioned on slack was.

        printSepToFile(size, components, orientation)
