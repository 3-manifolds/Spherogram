"""
Wrapper for Boyer's (C) planarity algorithm.
"""

cdef extern from "c/graph.h":

    ctypedef struct vertexRec:
        int  link[2]
        int  index
        unsigned int flags
    ctypedef vertexRec * vertexRecP

    ctypedef struct edgeRec:
        int  link[2]
        int  neighbor
        unsigned int flags
    ctypedef edgeRec * edgeRecP

    ctypedef struct baseGraphStructure:
        vertexRecP V
        edgeRecP E
        int N
        int M
    ctypedef baseGraphStructure * graphP
    
    cdef int OK, EMBEDFLAGS_PLANAR, NOTOK, NONEMBEDDABLE
    
    cdef graphP gp_New()
    cdef void gp_Free(graphP *pGraph)
    cdef int gp_InitGraph(graphP theGraph, int N)
    cdef int gp_AddEdge(graphP theGraph, int u, int ulink, int v, int vlink)
    cdef int gp_Embed(graphP theGraph, int embedFlags)
    cdef int gp_SortVertices(graphP theGraph)
    
def planar(fatgraph):
    if len(fatgraph.edges) == 0:
        return True, None
    cdef graphP theGraph
    cdef int status

    theGraph = gp_New()
    status = gp_InitGraph(theGraph, len(fatgraph.vertices))
    if status != OK:
        raise RuntimeError("gp_InitGraph status is not ok.")
    vertices = list(fatgraph.vertices)
    for edge in fatgraph.edges:
        start, end = edge
        m, n = vertices.index(start), vertices.index(end)
        if m != n:    # remove loops
            status = gp_AddEdge(theGraph, m, 0, n, 0)
            if status == NOTOK:
                raise RuntimeError("gp_AddEdge status is not ok.")
            elif status == NONEMBEDDABLE:
                return False, None
    status = gp_Embed(theGraph, EMBEDFLAGS_PLANAR)
    gp_SortVertices(theGraph)
    if status == NOTOK:
        raise RuntimeError("not ok.")
    if status == NONEMBEDDABLE:
        return False, None
    embedding = {}
    for i from 0 <= i < theGraph.N:
        adjacency_list = []
        j = theGraph.V[i].link[0]    # the first edge
        last = theGraph.V[i].link[1] # the last edge
        while j >= 0:
            adjacency_list.append(vertices[theGraph.E[j].neighbor])
            j = theGraph.E[j].link[0] # the next edge
        embedding[vertices[i]] = adjacency_list
    gp_Free(&theGraph)
    return True, embedding
