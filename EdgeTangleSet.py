import heapq
import functools
import multiprocessing
import sqlite3

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

# next two defs just stubs...
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
            newPartial = partialCut(Gdir, Gdircopy, fromDB=False, mincut=mincut, longpcut=newpcut, pcutlen=pcutlen)  # todo again, do I need a + 1?

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
def externalBasicPartitionBranch(uid, U, Gdirected):
    newpcut = [
        set(U[0:uid]).union([0]),
        set([U[uid]])
    ]
    Gdircopy, s, t = externalMergeVertices(Gdirected, newpcut)
    mincut = Gdircopy.mincut(source=s, target=t)

    if len(mincut.partition) == 2:
        return partialCut(Gdirected, Gdircopy, fromDB=False, mincut=mincut, longpcut=newpcut, pcutlen=uid+2) # todo note, 2 for: node 0, and 0 indexing
        # note that partialCut takes care of the vids re the adjusted graph

class partialCut(object):
    def __init__(self, Gdir=None, Gdircopy=None, fromDB=False, mincut=None, longpcut=None, pcutlen=0, dbRow=None):
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

        if fromDB:
            self.weight = dbRow[1]
            self.pcut = {
                "binrep": dbRow[2],
                "binlen": dbRow[3]
            }
            self.mcut = dbRow[4]
            self.cutEdges = frozenset(map(int, dbRow[5].split(",")))

        else:
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

        self.sepFilename = "{}/{}-SepList-CF.tsv". \
            format(job['outputFolder'], job['outName'])

        # todo Check if need to add orientation as per git.
        text = "order\tcut\tside1\tside2\torientation\n"
        with open(self.sepFilename, 'w+') as the_file:
            the_file.write(text)

        # todo;: do I want to put other stuff in a db?
        self.partcutConn = sqlite3.connect("partcutDB.sqlite")  # does this give rel path?
        partcutCursor = self.partcutConn.cursor()
        createTableStr = """ CREATE TABLE IF NOT EXISTS partcuts (
                            id integer PRIMARY KEY,
                            weight integer,
                            pcutrep integer,
                            pcutlen integer,
                            mcut integer,
                            edgeCuts text
                        ); """

        partcutCursor.execute(createTableStr)


    def __del__(self):
        pass
        # todo: this needs to delete all the partcuts from self.partcutDB


    # todo whack a numba thingy on this guy? - probably a wrapper
    # todo https://stackoverflow.com/questions/41769100/how-do-i-use-numba-on-a-member-function-of-a-class
    def findNextOrderSeparations(self, k=None, maxdepth=4):

        # from Yeh, Li-Pu, Biing-Feng Wang, and Hsin-Hao Su. 2010. “Efficient Algorithms for the Problems of Enumerating Cuts by Non-Decreasing Weights.” Algorithmica. An International Journal in Computer Science 56 (3): 297–312.
        def basicPartition(pool):
            print("Basic Partition")

            n = self.Gdirected.vcount()
            self.vids = self.Gdirected.vs.indices
            # let s = node 0, (see Yeh and Wang) so start at 1
            self.U = range(1,n)  # like this so easier to change later

            partcutList = [partcut for partcut in
                                pool.map(functools.partial(externalBasicPartitionBranch, U=self.U, Gdirected=self.Gdirected),
                                         range(len(self.U))) if partcut.weight <= self.kmax]
            # todo: consider changing to imap? to avoid double iteration? Mess with later
            # https://stackoverflow.com/questions/26520781/multiprocessing-pool-whats-the-difference-between-map-async-and-imap

            for partcut in partcutList:
                if partcut.weight <= self.kmax:
                    self.addtoPartcutDB(partcut)

        with multiprocessing.Pool() as pool:
            if k is None:  ### ie, first time running
                # ---------------------------------
                # slight waste of time, but saves memory by being able to avoid putting things on heap
                self.kmin = int(self.Gdirected.mincut().value)
                k = self.kmin
                self.kmax = k + maxdepth
                # ---------------------------------
                basicPartition(pool)
                # self.TangleTree.add_feature("cutsets", [])

                # todo note that the vertex IDs *should* be the same for G and Gdirected,
                # todo - and this will break if they are not

            # ----------------------------------------------------------------------------
            selectString = """ SELECT * FROM partcuts WHERE weight=? """
            deleteString = """ DELETE FROM partcuts WHERE id=? """
            checksizeString = """ SELECT COUNT(*) FROM partcuts WHERE weight=? """

            cur1 = self.partcutConn.cursor()
            while(cur1.execute(checksizeString, (k,)).fetchone()[0]):
                # I know this looks stupid, but partcutHeap gets modified by this loop
                # and we want to repeat until no more relevant partcuts in heap

                partcutList = []
                cur2 = self.partcutConn.cursor()
                cur2.execute(selectString, (k,))

                for row in cur2:
                    cur3 = self.partcutConn.cursor()
                    cur3.execute(deleteString, (row[0],))
                    # def __init__(self, Gdir, Gdircopy, fromDB=FALSE, mincut=None, longpcut=None, pcutlen=0, dbRow=None):
                    newpartcut = partialCut(fromDB=True, dbRow=row)
                    partcutList.append(newpartcut)
                    self.addToSepList(newpartcut)
                self.partcutConn.commit()

            # todo: convert to db - get *all* of k cuts, store, iter through and add in mem?
            # while len(self.partcutHeap) > 0 and self.partcutHeap[0].weight <= k:
            #     # I know this looks stupid, but partcutHeap gets modified by this loop
            #     # and we want to repeat until no more relevant partcuts in heap
            #     partcutList = []
            #     while len(self.partcutHeap) > 0 and self.partcutHeap[0].weight <= k:
            #         newpartcut = heapq.heappop(self.partcutHeap)
            #         partcutList.append(newpartcut)
            #         self.addToSepList(newpartcut)

                results = pool.map(functools.partial(externalExtractMinPart, Gdir=self.Gdirected, kmax=self.kmax), partcutList)
                for partcut in [item for subresults in results for item in subresults]:  # todo make sure returns work okay
                    self.addtoPartcutDB(partcut)

    def addtoPartcutDB(self, partcut):
        sql = """ INSERT INTO partcuts(weight,pcutrep,pcutlen,mcut,edgeCuts) VALUES(?,?,?,?,?) """

        dbTuple = (partcut.weight,
                   partcut.pcut["binrep"],
                   partcut.pcut["binlen"],
                   partcut.mcut,
                   ",".join(map(str, partcut.cutEdges)))
        partcutCursor = self.partcutConn.cursor()
        partcutCursor.execute(sql, dbTuple)
        self.partcutConn.commit()

    def addToSepList(self, partial):
        def cutIsSuperset(newCut):
            for cut in self.cuts:
                if newCut.issuperset(cut):
                    return True
            return False

        def printSepToFile(separation, cut, orientation):
            sideNodes = sorted([self.names[node] for node in separation[0]])
            complementNodes = sorted([self.names[node] for node in separation[1]])

            cutLong = []
            for eid in cut:

                edge = self.Gdirected.es[eid]
                eString = "('{}', '{}')".format(self.Gdirected.vs[edge.source]["name"], self.Gdirected.vs[edge.target]["name"])
                cutLong.append(eString)

            text = "{}\t{}\t{}\t{}\t{}\n".format(len(cut), sorted(cutLong), sideNodes, complementNodes, orientation)
            text = text.replace('\"', '')
            with open(self.sepFilename, 'a') as the_file:
                the_file.write(text)

        # def extractComponents(mcut):
        #     # probably there's a quicker way of doing this, but just get it working for now.
        #     mcutList = list('{0:b}'.format(mcut).zfill(self.Gdirected.vcount())[::-1])
        #     comps = [list(),list()]
        #     for i in len(mcutList):
        #         comps[int(mcutList[i])] = i
        #     return(comps)

        if cutIsSuperset(partial.cutEdges):
            return

        components = extractComponents(partial.mcut, self.Gdirected.vcount())
        components = sorted(components, key=len)

        sideNodes = frozenset({self.names[node] for node in components[0]})
        complementNodes = frozenset({self.names[node] for node in components[1]})

        size = len(partial.cutEdges)

        ######## ******* do other checks for "easy" seps (do shit here)
        ### smart:  initial check for def small - add more checks here?
        # todo - what about if both sides def small?
        # todo - add check for union of all def small sides? Only if *all* seps def orientable?
        if (len(components[0])==1) or (max(self.G.degree(components[0])) <= 2 and size >= 2) or (max(self.G.degree(components[0])) <= 1):
            self.definitelySmall[size].append(sideNodes)
            orientation = 1
        elif (len(components[1])==1) or (max(self.G.degree(components[1])) <= 2 and size >= 2) or (max(self.G.degree(components[1])) <= 1):
            self.definitelySmall[size].append(complementNodes)
            orientation = 2
        else:
            # Note: edited so only adding the shortest side.
            # separation = (sideNodes, complementNodes)
            # self.separations[size].append(separation)
            self.separations[size].append(sideNodes)
            orientation = 3

        # need to do this because we need to check for superset-ness
        self.cuts.add(partial.cutEdges)
        printSepToFile(components, partial.cutEdges, orientation)
