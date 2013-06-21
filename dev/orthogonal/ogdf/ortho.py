import networkx, snappy, spherogram, os, sys, re, collections

def random_link():
    return spherogram.DTcodec(snappy.HTLinkExteriors.random().DT_code()).link()

def link_to_gml_file(L, filename='graphs/link.gml'):
    # Have to add dummy vertices in the twist regions
    G = networkx.MultiGraph()
    verts_to_int = { v:i for i, v in enumerate(L.vertices) }
    verts = sorted(verts_to_int.values())
    edges = collections.Counter(tuple(sorted( [verts_to_int[e.tail], verts_to_int[e.head]])) for e in L.edges)
    assert set(edges.values()).issubset(set([1,2]))
    for (a,b), m in edges.items():
        if m == 2:
            c = len(verts)
            verts.append(c)
            edges[ (a,c) ] = 1
            edges[ (c, b) ] = 1
    for a,b in edges:
        G.add_edge( a, b)
    networkx.write_gml(G, filename)

def test_orthogonal_drawing(L):
    link_to_gml_file(L)
    print os.popen("./ortho").read()[:-1]
    os.system('open -a "Google Chrome" graphs/link.svg')
