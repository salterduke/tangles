import igraph as ig
import matplotlib.pyplot as plt

basefolder = "/home/saltermich1/PhDThesisGdrive/"
# basefolder = "C:/Users/mrousset/Documents/PhDThesisLaptop/"
outfolder = basefolder + "Code/Tangles/testResources/diagrams/"

# --------------------------------------------------------------------
# tinyEdges
def drawGraphExample(format = "png"):
    # gFile = "C:/Users/mrousset/Documents/PhDThesisLaptop/Code/NetworkData/SmallNWs/TinyEdges.csv"
    gFile = basefolder + "Code/NetworkData/SmallNWs/TinyEdges.csv"
    layout = [(0, 0), (0, 2), (1, 1), (2, 1), (2, 0), (3, 1), (4, 0), (5, 1), (4, 2)]
    G = ig.Graph.Read_Ncol(gFile, names=True, directed=False)

    visual_style = {}
    visual_style["vertex_color"] = "white"
    visual_style["vertex_label"] = G.vs["name"]
    # visual_style["vertex_label"] = G.vs.indices
    visual_style["bbox"] = (0, 0, 500, 250)

    outfile = outfolder + "tinyedges." + format

    ig.plot(G, target=outfile, layout=layout, **visual_style)
    # g = ig.plot(G, layout=layout, target=outfile, **visual_style)
    # g.show()

# --------------------------------------------------------------------
# branch decomp of tinyEdges

def drawBranchDecomp(format = "png"):
    # bDecompFile = "C:/Users/mrousset/Documents/PhDThesisLaptop/Code/NetworkData/SmallNWs/TinyBranchDecomp.csv"
    bDecompFile = basefolder + "Code/NetworkData/SmallNWs/TinyBranchDecomp.csv"

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
    outfile = outfolder + "tinybranchdecomp." + format
    g.save(outfile)

# --------------------------------------------------------------------
# carving decomp of tinyEdges

def drawCarving(format = "png"):
    # cDecompFile = "C:/Users/mrousset/Documents/PhDThesisLaptop/Code/NetworkData/SmallNWs/TinyCarvingDecomp.csv"
    cDecompFile = basefolder + "Code/NetworkData/SmallNWs/TinyCarvingDecomp.csv"

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
    # g.show()

    outfile = outfolder + "tinycarvingdecomp." + format
    g.save(outfile)

# --------------------------------------------------------------------
# tree diagram for VY, YWS
def drawCutFinderTree(format = "png"):
    # treeFile = "C:/Users/mrousset/Documents/PhDThesisLaptop/Code/NetworkData/SmallNWs/ForDiagrams/treeDiagVY.csv"
    treeFile = basefolder + "Code/NetworkData/SmallNWs/ForDiagrams/treeDiagVY.csv"

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

    outfile = outfolder + "cutfinder." + format
    g.save(outfile)


# _____________________________________________

# drawCutFinderTree()

# note that it appears that saving as pdf *doesn't* mean it's saved as a vector format,
# it's just an image within a pdf. Similarly for svg. (Unless it's done as target = outfile?)
outformat = "pdf"

drawGraphExample(format = outformat)
drawBranchDecomp(format = outformat)
drawCarving(format = outformat)
