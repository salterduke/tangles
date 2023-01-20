import igraph as ig
import matplotlib.pyplot as plt

outfolder = "C:/Users/mrousset/Documents/PhDThesisLaptop/Code/Tangles/testResources/diagrams/"

# --------------------------------------------------------------------
# tinyEdges
def drawGraphExample():
    gFile = "C:/Users/mrousset/Documents/PhDThesisLaptop/Code/NetworkData/SmallNWs/TinyEdges.csv"
    layout = [(0, 0), (0, 2), (1, 1), (2, 1), (2, 0), (3, 1), (4, 0), (5, 1), (4, 2)]
    G = ig.Graph.Read_Ncol(gFile, names=True, directed=False)

    visual_style = {}
    visual_style["vertex_color"] = "white"
    visual_style["vertex_label"] = G.vs["name"]
    visual_style["bbox"] = (0, 0, 500, 250)

    outfile = outfolder + "tinyedges.png"

    g = ig.plot(G, layout=layout, **visual_style)
    # g = ig.plot(G, layout=layout, target=outfile, **visual_style)
    # g.show()
    g.save(outfile)

# --------------------------------------------------------------------
# branch decomp of tinyEdges

def drawBranchDecomp():
    bDecompFile = "C:/Users/mrousset/Documents/PhDThesisLaptop/Code/NetworkData/SmallNWs/TinyBranchDecomp.csv"

    G = ig.Graph.Read_Ncol(bDecompFile, names=True, directed=False)
    G.es["curved"] = 0

    layout = G.layout_reingold_tilford()
    # layout = [(0,0), (1,1), (0,2), (2,1), (2,2), (3,1), (3,2), (4,1), (4,2), (5,1), (5,2), (6,1), (6,2), (7,1), (7,2), (8,1), (8,2), (9,1), (9,2), (10,1), (11,2), (11,0)]

    visual_style = {}
    visual_style["vertex_color"] = "white"
    visual_style["vertex_label"] = list(n if not n.isnumeric() else "" for n in G.vs["name"] )
    visual_style["vertex_size"] = list(25 if not n.isnumeric() else 1 for n in G.vs["name"] )
    visual_style["bbox"] = (0, 0, 500, 250)
    visual_style["edge_label"] = list(" {}".format(int(w)) if w>0 else ""  for w in G.es["weight"])
    # visual_style["edge_label_dist"] = 50
    # visual_style["edge_background"] = "white"


    g = ig.plot(G, layout=layout, **visual_style)
    # g.show()
    outfile = outfolder + "tinybranchdecomp.png"
    g.save(outfile)

# --------------------------------------------------------------------
# carving decomp of tinyEdges

def drawCarving():
    cDecompFile = "C:/Users/mrousset/Documents/PhDThesisLaptop/Code/NetworkData/SmallNWs/TinyCarvingDecomp.csv"

    G = ig.Graph.Read_Ncol(cDecompFile, names=True, directed=False)
    G.es["curved"] = 0

    layout = G.layout_reingold_tilford()
    # layout = [(0,0), (1,1), (0,2), (2,1), (2,2), (3,1), (3,2), (4,1), (4,2), (5,1), (5,2), (6,1), (6,2), (7,1), (7,2), (8,1), (8,2), (9,1), (9,2), (10,1), (11,2), (11,0)]

    visual_style = {}
    visual_style["vertex_color"] = "white"
    visual_style["vertex_label"] = list(n if not n.isnumeric() else "" for n in G.vs["name"] )
    visual_style["vertex_size"] = list(25 if not n.isnumeric() else 1 for n in G.vs["name"] )
    # visual_style["bbox"] = (0, 0, 500, 130)
    visual_style["edge_label"] = list("{}".format(int(w)) if w>0 else "" for w in G.es["weight"])
    # visual_style["edge_label_dist"] = 50
    # visual_style["edge_background"] = "white"

    # sugi_layers = [4, 3, 4, 2, 4, 1, 3, 4, 4, 0, 1, 4, 2, 4, 3, 4, 4]
    # layout = G.layout_sugiyama()

    # layout = G.layout_reingold_tilford(root=0)
    layout = G.layout_reingold_tilford()
    visual_style["bbox"] = (0, 0, 500, 250)
    g = ig.plot(G, layout=layout, **visual_style)
    g.show()

    outfile = outfolder + "tinycarvingdecomp.png"
    g.save(outfile)

# --------------------------------------------------------------------
# tree diagram for VY, YWS
def drawCutFinderTree():
    treeFile = "C:/Users/mrousset/Documents/PhDThesisLaptop/Code/NetworkData/SmallNWs/ForDiagrams/treeDiagVY.csv"

    G = ig.Graph.Read_Ncol(treeFile, names=True, directed=False)
    G.es["curved"] = 0

    layout = G.layout_reingold_tilford()
    # layout = [(0,0), (1,1), (0,2), (2,1), (2,2), (3,1), (3,2), (4,1), (4,2), (5,1), (5,2), (6,1), (6,2), (7,1), (7,2), (8,1), (8,2), (9,1), (9,2), (10,1), (11,2), (11,0)]

    visual_style = {}
    visual_style["vertex_color"] = "white"
    visual_style["vertex_label"] = list(n if not n.isnumeric() else "" for n in G.vs["name"] )
    visual_style["vertex_size"] = list(25 if not n.isnumeric() else 1 for n in G.vs["name"] )
    # visual_style["bbox"] = (0, 0, 500, 130)
    visual_style["edge_label"] = list("{}".format(int(w)) if w>0 else "" for w in G.es["weight"])
    # visual_style["edge_label_dist"] = 50
    # visual_style["edge_background"] = "white"

    # sugi_layers = [4, 3, 4, 2, 4, 1, 3, 4, 4, 0, 1, 4, 2, 4, 3, 4, 4]
    # layout = G.layout_sugiyama()

    # layout = G.layout_reingold_tilford(root=0)
    layout = G.layout_reingold_tilford()
    visual_style["bbox"] = (0, 0, 500, 250)
    g = ig.plot(G, layout=layout, **visual_style)
    g.show()

    outfile = outfolder + "cutfinder.png"
    g.save(outfile)


# _____________________________________________

drawCutFinderTree()

# drawGraphExample()
# drawBranchDecomp()
# drawCarving()