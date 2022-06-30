import igraph as ig
# import numpy as np
# import pandas as pd
# import random


minDist = min({self.d[j] for j in j_in_W
               if (self.H.es[self.H.get_eid(i, j)]["weight"] - self.H.es[self.H.get_eid(i, j)]["flow"] + self.H.es[self.H.get_eid(j, i)]["flow"] ) > 0})



print(ig.__version__)

#
#
# G = ig.Graph.Famous("Zachary")
#
g= ig.Graph.Full(n=4)
t = g.gomory_hu_tree(flow="label")
t.es["curved"] = 0

e = t.es[0]


# t = g.gomory_hu_tree()
visual_style = {}
# visual_style["vertex_color"] = "white"
visual_style["vertex_label"] = range(4)
# ig.plot(t, **visual_style)

ig.plot(t)

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
G = ig.Graph.Forest_Fire(50,0.75);G.layout_fruchterman_reingold();ig.plot(G)
for p in range(70,100,1):
    G = ig.Graph.Forest_Fire(50, p/100);print(p, G.vcount(), G.ecount())

visual_style = {}
visual_style["vertex_size"] = 40
visual_style["vertex_label"] = t.vs["name"]
visual_style["vertex_color"] = "white"
visual_style["vertex_label_size"] = 8
visual_style["edge_width"] = [1 + 2 * int(is_formal) for is_formal in g.es["is_formal"]]
visual_style["layout"] = layout
visual_style["bbox"] = (300, 300)
#
# # ig.plot(t, **visual_style)
#
#
#
[e for e in y if e not in x]
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
# # # #     df = pd.read_csv(file, sep="\t", header=None)
# # #     df = pd.read_csv(file, sep=" ")
# # #     print(df.head())
# # #     # print("---------------")
# # #     # print(file)
# # #     m = len(df.index)
# # #     Vs = pd.concat([df.iloc[:,0], df.iloc[:,1]])
# # #     n = len(Vs.unique())
# # #     print(file, m, n)
# # #     # df.iloc[:,0:2].to_csv("YeastNet.csv", index=False, sep = " ", header = ["A", "B"])
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
