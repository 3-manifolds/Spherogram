try:
    from sage.all import *
except ImportError:
    pass
    
import spherogram

G = spherogram.Graph()
G.add_vertex(0)
G.add_vertex('a')
G.add_vertex('c')
G.add_edge(0, 'a')
G.add_edge(0, 'a')
G.add_edge('c', 'a')
G.add_edge('c', 0)
G.add_edge('c', 'c')

def convert(G):
    S = Graph(multiedges=True)
    S.add_vertices(G.vertices)
    for e in G.edges:
        S.add_edge(e.ends[0], e.ends[1], repr(e))
    return S

#print G
R = G.reduced()
#print R
print(G.is_planar())
print(R.embedding())
