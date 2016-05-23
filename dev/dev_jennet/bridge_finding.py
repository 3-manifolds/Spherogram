"""
Used to find all bridges of a graph G in an efficient way.
Follows the solution to 4.1.36 in Sedgewick's _Algorithms_ (4th ed.)
"""

def find_bridges(G):
    global cnt
    global low # low[v] is the lowest preorder number of any vertex connected to v
    global preorder
    global bridges

    cnt = 0    
    low = dict()
    pre = dict()
    bridges = []

    verts = G.vertices()
    low = dict([(v,-1) for v in verts])
    preorder = dict([(v,-1) for v in verts])

    for v in verts:
        if preorder[v] == -1:
            _recursive_bridge_finding(G, v, v)

    return bridges

def _recursive_bridge_finding(G, u, v):
    global cnt
    global low
    global preorder
    global bridges
  
    cnt += 1
    preorder[v] = cnt
    low[v] = preorder[v];
    for w in G.neighbor_iterator(v):
        if (preorder[w] == -1):
            _recursive_bridge_finding(G, v, w)
            low[v] = min(low[v], low[w])
            if (low[w] == preorder[w]):
                #print(repr(v) + "-" + repr(w) + " is a bridge")
                bridges.append( (min(v,w), max(v,w), None) )
        # update low number, ignoring reverse of edge leading to v
        # if there are multiple edges from u to v, ignore only one
        elif (u != w or connections(G,u,v) > 1):
            low[v] = min(low[v], preorder[w])
    
def connections(G,a,b):
    return len([x for x in G.edges() if x[0] == a and x[1] == b] + [x for x in G.edges() if x[0] == b and x[1] == a])
    
