def orient_tree(T, root):
    "Returns a digraph by orienting T outwards from root."
    if not T.is_tree():
        raise Exception("Input must be a tree.")
    if not T.edges():
        verts = T.vertices()
        v = verts[0]
        return graph.DiGraph({v:list()})
    H = T.copy()
    H.delete_vertex(root)
    root_edges = list()
    for e in T.edges_incident(root):
        if e[0] == root:
            root_edges.append(e)
        else:
            root_edges.append((e[1],e[0],e[2]))
    pieces = list() #Orient the subtrees gotten by deleting root.
    for e in root_edges:
        comp = H.subgraph(H.connected_component_containing_vertex(e[1]))
        pieces.append(orient_tree(comp,e[1]))
    all_edges = root_edges
    for P in pieces:
        all_edges = all_edges+P.edges()
    return graph.DiGraph(all_edges)

def _crossing_to_edge(K,G,c):
    "Returns the edge of G corresponding to crossing c, where G is either the black graph or white graph."
    edge = list()
    verts = G.vertices()
    black_faces = list() #Strip off the corners (0, 1, 2 or 3) from the crossings
    for v in verts:
        black_faces.append([x[0] for x in v])
    for i in range(len(black_faces)):
        if c in black_faces[i]:
            edge.append(verts[i])
    edge.append(c)
    if len(edge) != 3:
        raise Exception("Did not find two faces incident to c="+repr(c))
    return tuple(edge)

def _dual_spanning_tree(K,T):
    "Returns the spanning tree for the white graph dual to the spanning tree T for the black graph. Here, dual means the edges do not intersect, i.e., no common crossings."
    edges = list()
    G = K.white_graph()
    crossings_to_ignore = [e[2] for e in T.edges()]
    crossings_to_use = [c for c in K.crossings if not c in crossings_to_ignore]
    return Graph([K._crossing_to_edge(G,c) for c in crossings_to_use])

def Kauffman_states(K):
    "Returns the set of Kauffman states for the Alexander polynomial, corresponding to the spanning trees in the black graph. Returns a list of dictionaries, with keys crossings and values faces in the knot projection."
    G = K.black_graph()
    trees = G.spanning_trees()
    marked_edge = ((K.crossings[0],0), K.crossings[0].adjacent[0]) #We (arbitrarily) mark an edge in K.
    states = list()
    for T in trees:
        #Find the root of T. 
        for v in T.vertices():
            if marked_edge[0] in v:
                root = v
        #Orient T out from the root.
        oT = orient_tree(T,root)
        #Find the planar dual spanning tree for the white graph
        dT = _dual_spanning_tree(K,T)
        #Find the root of dT. 
        for v in dT.vertices():
            if marked_edge[0] in v:
                droot = v
        odT = orient_tree(dT,droot)
        state = dict()
        for e in oT.edges():
            state[e[2]] = e[1]
        for e in odT.edges():
            state[e[2]] = e[1]
        states.append(state)
    return states
