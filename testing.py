import igraph as ig
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# import random
from sklearn.linear_model import LinearRegression
import sys
from scipy.stats import entropy

def l2(x):

    if x == 0:
        return x
    else:
        return np.log2(x)


def le(x):
    if x == 0:
        return x
    else:
        return np.log(x)

def MI_e(P):
    # todo note these have dims hard coded!!!
    total = 0
    for i in range(3):
        for j in range(2):
            # print(i,j,P[i,j],P.sum(axis=1)[i],P.sum(axis=0)[j],P[i,j]*l2(P[i,j] / (P.sum(axis=1)[i] * P.sum(axis=0)[j])))
            total+= P[i,j]*le(P[i,j] / (P.sum(axis=1)[i] * P.sum(axis=0)[j]))
    return total


def MI(P):
    # todo note these have dims hard coded!!!
    total = 0
    for i in range(3):
        for j in range(2):
            # print(i,j,P[i,j],P.sum(axis=1)[i],P.sum(axis=0)[j],P[i,j]*l2(P[i,j] / (P.sum(axis=1)[i] * P.sum(axis=0)[j])))
            total+= P[i,j]*l2(P[i,j] / (P.sum(axis=1)[i] * P.sum(axis=0)[j]))
    return total


def NMI(N):
    top = 0
    bleft = 0
    bright = 0
    for i in range(3):
        bleft += N.sum(axis=1)[i]*l2(N.sum(axis=1)[i]/N.sum())
        for j in range(2):
            top += N[i,j]*l2(N[i,j]*N.sum()/(N.sum(axis=1)[i]*N.sum(axis=0)[j]))

    for j in range(2):
        bright += N.sum(axis=0)[j] * l2(N.sum(axis=0)[j] / N.sum())

    return(-2*top/(bleft+bright))

NMI(N)


def NMI_e(N):
    top = 0
    bleft = 0
    bright = 0
    for i in range(3):
        bleft += N.sum(axis=1)[i]*le(N.sum(axis=1)[i]/N.sum())
        for j in range(2):
            top += N[i,j]*le(N[i,j]*N.sum()/(N.sum(axis=1)[i]*N.sum(axis=0)[j]))

    for j in range(2):
        bright += N.sum(axis=0)[j] * le(N.sum(axis=0)[j] / N.sum())

    return(-2*top/(bleft+bright))


NMI(N1)
NMI_e(N1)

N = np.array([3,0,2,0,0,4]).reshape((3,2))
P = N/9

N1 = np.array([2,1,1,0,3,2]).reshape((3,2))
P1 = N1/9


I = MI(P)
HC2 = entropy(P.sum(axis=0), base=2)
HC1 = entropy(P.sum(axis=1), base=2)

I_e = MI_e(P)
HC2_e = entropy(P.sum(axis=0)) # base e
HC1_e = entropy(P.sum(axis=1))

myVI_e = HC1_e + HC2_e - 2*I_e

myVI_e - VI

m1 = [0, 0, 0, 1, 1, 2, 2, 2, 2]
m2 = [0, 0, 0, 0, 0, 1, 1, 1, 1]

VI = ig.compare_communities(m1, m2, method = "vi")
NMI_ig = ig.compare_communities(m1, m2, method = "nmi")

(HC1 + HC2 - VI)/2

myVI = HC1 + HC2 - 2*I
myVI/VI
VI/myVI

# G = ig.Graph.Read_GML("../NetworkData/MediumSize/adjnoun.gml")
#
# visual_style = {}
# visual_style["vertex_color"] = "white"
# visual_style["vertex_label"] = G.vs.indices
# ig.plot(G, **visual_style)
#
# G = G.simplify()
#
# G.write_ncol("../NetworkData/MediumSize/Copperfield.csv", names="label", weights=None)
#
# def longp(l):
#     print("\n".join(map(str, l)))


# print(sys.argv[0])
# print(sys.argv[1])
# print(len(sys.argv))




# minDist = min({self.d[j] for j in j_in_W
#                if (self.H.es[self.H.get_eid(i, j)]["weight"] - self.H.es[self.H.get_eid(i, j)]["flow"] + self.H.es[self.H.get_eid(j, i)]["flow"] ) > 0})
#
#
#
# print(ig.__version__)
#
# results = []
#
# for v in range(10,100,10):
#     for e in range(v+10, 3*v, 10):
#         e_actual = 0
#         fudge = 0
#         print("--------------------------")
#         while e_actual < e:
#             p = (e - 3.52663056*v +201.35998895)/(206.97060932) + fudge
#             if p >= 1:
#                 p = 0.99999999
#             if p <= 0:
#                 p = 0.00000001
#             G = ig.Graph.Forest_Fire(v, p)
#             fudge += ((e - e_actual)/(e*100) * (1 - p))
#             # print(p, G.vcount(), G.ecount())
#         #     results.append((p, G.vcount(), G.ecount()))
#             results.append((G.vcount(), e, G.ecount()))
#             e_actual = G.ecount()
#             print(p, fudge, (G.vcount(), e, G.ecount()))
#
# df = pd.DataFrame(results, columns = ["V","e","E"])
# # summ = results_wide.groupby(["p", "V"]).to_numpy().mean()

# p = summ.index.array.reshape(-1,1)
# V = summ.V.array.reshape(-1,1)
# E = summ.E.array.reshape(-1,1)
# pV = np.hstack((p,V))
# regr = LinearRegression()
# fit = regr.fit(pV,E)
#
# print(fit.intercept_, fit.coef_)



# summ.plot(y="E")
# plt.show()




#
#
# G = ig.Graph.Famous("Zachary")
#
# g= ig.Graph.Full(n=4)
# t = g.gomory_hu_tree(flow="label")
# t.es["curved"] = 0
#
# e = t.es[0]
#
#
# # t = g.gomory_hu_tree()
visual_style = {}
# visual_style["vertex_color"] = "white"
visual_style["vertex_label"] =
ig.plot(t, **visual_style)
#
# ig.plot(t)

#
# # todo this is just scratch stuff - delete
#
# cmake -D CMAKE_C_COMPILER=gcc-4.2 -D CMAKE_CXX_COMPILER=g++-4.2 path/to/your/source
#
#
# gname = "C:/Users/mrousset/Documents/PhDThesisLaptop/Code/NetworkData/BioDBs/YeastPPI/YuEtAlGSCompB.csv"
# # gname = "C:/Users/mrousset/Documents/PhDThesisLaptop/Code/NetworkData/SmallNWs/TinyEdges.csv"
# G = ig.Graph.Read_Ncol(gname, names=True, directed=False)
#
# # origVs = pd.DataFrame({"name": G.vs["name"],
# #                        "deg": G.degree()})
#
# shell = np.array(G.coreness())
#
# del_Vs = []
# for vid in np.where(shell==1)[0]:
#     to_del = True
#     for nb in G.neighbors(vid):
#         if shell[nb] > 1:
#             to_del = False
#     if to_del:
#         del_Vs.append(vid)
#
# G.delete_vertices(del_Vs)
#
#
# # ig.plot(G)
#
# GH = G.gomory_hu_tree(flow = "label")
# # GH = G.gomory_hu_tree()
#
# # GH.contract_vertices([0,0,0,1,2,3,3,3,3], combine_attrs="concat")
#
#
# G = ig.Graph.Forest_Fire(50,0.75);G.layout_fruchterman_reingold();ig.plot(G)

# visual_style = {}
# visual_style["vertex_size"] = 40
# visual_style["vertex_label"] = G.vs.indices
# visual_style["vertex_color"] = "white"
# visual_style["vertex_label_size"] = 8
# # visual_style["edge_width"] = [1 + 2 * int(is_formal) for is_formal in g.es["is_formal"]]
# # visual_style["layout"] = layout
# visual_style["bbox"] = (300, 300)
#
# # ig.plot(t, **visual_style)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# # fbase = "C:/Users/mrousset/Documents/PhDThesisLaptop/Code/NetworkData/Constructed/"
# # fullSize = 2016
# # edgesToDel = 252
# #
# #
# # G = ig.Graph.Full(64)
# # fname = fbase + "Graph_v_64_es_2016.ncol"
# # G.write_ncol(fname, None, None)
# #
# # for i in range(1, 8):
# #     delCount = 0
# #     print("Graph {}".format(i))
# #     while(delCount < edgesToDel):
# #         Es = list(G.es.indices)
# #         newEs = random.sample(Es, len(Es)-1)
# #         newG = G.subgraph_edges(newEs, delete_vertices = False)
# #         if (newG.is_connected()) :
# #             delCount += 1
# #             G = newG
# #     fname = "{}Graph_v_64_es_{}.ncol".format(fbase, G.ecount())
# #     print("---------------------")
# #     print("Number of nodes: {}".format(G.vcount()))
# #     print("Number of edges: {}".format(G.ecount()))
# #     G.write_ncol(fname, None, None)
# #
# # edgesToDel = 21
# #
# # for i in range(1, 10    ):
# #     delCount = 0
# #     print("Graph {}".format(i))
# #     while(delCount < edgesToDel):
# #         Es = list(G.es.indices)
# #         newEs = random.sample(Es, len(Es)-1)
# #         newG = G.subgraph_edges(newEs, delete_vertices = False)
# #         if (newG.is_connected()) :
# #             delCount += 1
# #             G = newG
# #     fname = "{}Graph_v_64_es_{}.ncol".format(fbase, G.ecount())
# #     print("---------------------")
# #     print("Number of nodes: {}".format(G.vcount()))
# #     print("Number of edges: {}".format(G.ecount()))
# #     G.write_ncol(fname, None, None)
#
#
# # # import sqlite3
# # import pandas as pd
# #
# # longlist = ("YuEtAlCCSB.csv","YuEtAlGSComp4.csv","YuEtAlGSCompB.csv","YuEtAlRRS.csv","YuEtAlGS.csv","YuEtAlGSCompA.csv","YuEtAlY2H.csv","YuEtAlGSComp3.csv","YuEtAlPRS.csv")
# # shortlist = ("GT.E-Map.PMID-16269340.YeastNet.v3.138gene.209link.txt",)
# #
# # shortlist = ("YeastNet.csv",)
# #
# # # for file in shortlist:
# # # # for file in longlist:
# # # #     results_wide = pd.read_csv(file, sep="\t", header=None)
# # #     results_wide = pd.read_csv(file, sep=" ")
# # #     print(results_wide.head())
# # #     # print("---------------")
# # #     # print(file)
# # #     m = len(results_wide.index)
# # #     Vs = pd.concat([results_wide.iloc[:,0], results_wide.iloc[:,1]])
# # #     n = len(Vs.unique())
# # #     print(file, m, n)
# # #     # results_wide.iloc[:,0:2].to_csv("YeastNet.csv", index=False, sep = " ", header = ["A", "B"])
# #
# #
# # print((364+18+128+16+24)*25000/(1024*1024))
# #
# #
# #
# #
# #
# #
# #
# # #
# # # partcutConn = sqlite3.connect("partcutDB.sqlite")  # does this give rel path?
# # # partcutCursor = partcutConn.cursor()
# # # createTableStr = """ CREATE TABLE IF NOT EXISTS partcuts (
# # #                     id integer PRIMARY KEY,
# # #                     weight integer,
# # #                     pcutrep integer,
# # #                     pcutlen integer,
# # #                     mcut integer,
# # #                     edgeCuts text
# # #                 ); """
# # #
# # # partcutCursor.execute(createTableStr)
# # #
# # # sql = """ INSERT INTO partcuts(weight,pcutrep,pcutlen,mcut,edgeCuts) VALUES(?,?,?,?,?) """
# # #
# # #
# # # partcutCursor = partcutConn.cursor()
# # # # partcutCursor.execute(sql, (1,1,3,4,"123123"))
# # # # partcutCursor.execute(sql, (1,2,3,4,"123123"))
# # # # partcutCursor.execute(sql, (1,3,3,4,"123123"))
# # # # partcutCursor.execute(sql, (1,4,3,4,"123123"))
# # # # partcutConn.commit()
# # #
# # # selectAll = """ SELECT * FROM partcuts """
# # # selectString = """ SELECT * FROM partcuts WHERE weight=? """
# # # deleteString = """ DELETE FROM partcuts WHERE id=? """
# # # checksizeString = """ SELECT COUNT(*) FROM partcuts WHERE weight=? """
# # #
# # # cur1 = partcutConn.cursor()
# # # # cur1.execute(selectString, (1,))
# # #
# # # partcutList = []
# # #
# # # # print(cur1.execute(checksizeString, (1,)).fetchone()[0])
# # # # print(cur1.fetchone()[0])
# # #
# # # cur1.execute(selectAll)
# # # for row in cur1:
# # #     print(row)
# # #     # cur2 = self.partcutConn.cursor()
# # #     # cur2.execute(deleteString, (c[0],))
# # #
# # # # self.partcutConn.commit()


