import numpy as np
import functools

def plotG(G, outfile):

    visual_style = {}
    visual_style["vertex_label"] = G.vs.indices

    visual_style["edge_label"] = G.es["weight"]
    visual_style["bbox"] = (0, 0, 14 * 50, 14 * 50)

    pal = ig.GradientPalette("black", "white", 3)
    visual_style["vertex_color"] = G.vs["vertex_color"]

    G.es["curved"] = 0
    layout = G.layout_grid()
    output = ig.plot(G, layout=layout, **visual_style)
    output.save(outfile)

def ifelse(cond, valifTrue, valifFalse):
    if cond:
        return valifTrue
    else:
        return valifFalse

# returns a list with dim the num of cols in df1, recording whether that col exists in df2
# sorts the rows by index to ensure comparisons correct
def matchCols(df1, df2):
    df1 = df1.sort_index(axis="index")
    df2 = df2.sort_index(axis="index")
    return [any([df1[j].equals(df2[i]) for i in df2.columns]) for j in df1.columns]


def pruneToStubs(G):
    G.simplify(combine_edges="sum")

    # Takes every "twig" hanging off the edges of the main graph, and condenses it to a single node as a "stub"
    # since a twig is always "small", don't need to deal with all parts of the twig.
    outershell = G.induced_subgraph(np.where(np.array(G.coreness()) == 1)[0])
    twigs = outershell.clusters()

    mapvector = G.vs.indices
    for twig in twigs:
        vnames = outershell.vs[twig]["name"]
        vids_inG = G.vs.select(name_in=vnames).indices
        for i in vids_inG:
            mapvector[i] = min(vids_inG)
            # probably a oneliner way of doing this, but eh.

    G.contract_vertices(mapvector, combine_attrs=dict(name=functools.partial(mergeVnames, sep=";")))
    G.simplify(combine_edges="sum")
    G.delete_vertices([v.index for v in G.vs if v["name"] == ''])
    return G

def mergeVnames(names, sep=","):
    return sep.join(names)
