#cython: language_level=3

from libc.stdlib cimport malloc, free
import random

cdef extern from 'PMdef.h':
    ctypedef struct pmSize:
        char m, b  #  map and basic map type
        long e, v, f  # edges, vertices, faces 
        long r, g, d  # red and black vertices, max degree
        long t        # tolerance on e
        long *dgArr   # pt on vertex list

    ctypedef struct pmMethod:
        char core, pic
        long seed
        char verbose

    ctypedef struct pmMemory:
        char dTree
        long sTree, rTree, gTree, sWrd, sEdge, sVtx, sLeaf

    ctypedef struct pm_vertex:
        pm_edge* root
        pm_vertex* next
        long mark
        short type
        long label
        
    ctypedef struct pm_edge:
        pm_vertex* c_from "from"
        pm_vertex* face
        pm_edge* prev
        pm_edge* next
        pm_edge* oppo
        long mark
        short type
        long label

    ctypedef struct pmMap:
        pm_edge *root
        long e, v, f, i

    void set_pmRandom_callback(long (*function)(long))

cdef extern from 'PMplanmap.h':
    extern int pmMemoryInit(pmSize *S, pmMethod *Meth, pmMemory *M)
    extern int pmPlanMap(pmSize *S, pmMethod *Meth,
                         pmMemory *M, pmMap *Map)
    extern int pmFreeMap(pmMap *Map)

cdef extern from 'stats.h':
    extern long pmStatGauss(pmMap *Map)

# We want Planarmap to use the same pseudo-random number generator as Python.

cdef long randrange_callback(long n):
    return random.randrange(n) + 1

set_pmRandom_callback(randrange_callback)


# The main function

def random_map(num_vertices, int edge_connectivity=4,
               int num_link_comps=0, int max_tries=100):
    """
    Use Gilles Schaeffer's "Planarmap program" to generate
    a random 4-valent planar graph with the given number
    of vertices.

    The "edge_connectivity" parameter can be any of 2, 4, or 6.
    Recall that a graph is k-edge-connected if removing any set
    of *less than* k edges disconnects the graph.  In particular,
    for 4-valent graphs, being 2-connected is just the same as
    being connected.  In particular, a 2-connected graph can
    (and frequently do) have looped edges.

    Under the hood, it uses Python's pseudo-random number
    generator.
    """
    cdef pmSize size
    cdef pmMethod method
    cdef pmMemory memory
    cdef pmMap the_map
    cdef pm_edge *edge
    cdef pm_vertex *vert

    if edge_connectivity==2:
        size.m, size.b = 4, 4
    elif edge_connectivity==4:
        size.m, size.b = 5, 5
    elif edge_connectivity==6:
        size.m, size.b = 6, 5
    else:
        raise ValueError("Invalid edge_connectivity parameter")

    size.v = num_vertices
    size.e, size.f, size.r, size.g, size.d = 0, 0, 0, 0, 0
    size.t = -1
    size.dgArr = NULL

    method.core, method.pic = 0, 0
    method.verbose = 0
    pmMemoryInit(&size, &method, &memory)
    pmPlanMap(&size, &method, &memory, &the_map)
    if num_link_comps > 0:
        tries = 0
        while pmStatGauss(&the_map) != num_link_comps and tries < max_tries:
            pmFreeMap(&the_map)
            pmPlanMap(&size, &method, &memory, &the_map)
            tries += 1

        if tries==max_tries:
            return 

    ans = []
    vert = the_map.root.c_from
    while vert != NULL:
        edges_at_vert = []
        edge = vert.root
        while edge != vert.root.prev:
            edges_at_vert.append(edge.label)
            edge = edge.next
        edges_at_vert.append(edge.label)
        ans.append( (vert.label, tuple(edges_at_vert)) )
        vert = vert.next
    pmFreeMap(&the_map)
    return ans
