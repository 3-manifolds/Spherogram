from tangle_patch import *
import snappy
from sage.graphs.graph import Graph

links = snappy.HTLinkExteriors(crossings=12, knots_vs_links='knots')

"""
for M in links:
    L = M.link()
    mutants = all_mutants(L)
    ids = [str(K.exterior().identify()[-1]) for K in mutants]
    print(ids)
"""

"""
mutation_graph = Graph()
i = 0
for M in links:
    i += 1
    print('%d / %d' % (i, len(links)))
    L = M.link()
    try:
        M_id = str(M.identify()[-1])
    except:
        continue
    mutants = all_mutants(L)
    different_from_original = []
    for mutant in mutants:
        mutant_ext = mutant.exterior()
        try:
            is_isometric = mutant_ext.is_isometric_to(M)
        except:
            continue
        if not mutant_ext.is_isometric_to(M):
            different_from_original.append(mutant)
    try:
        ids = [str(K.exterior().identify()[-1]) for K in different_from_original]
    except:
        continue
    for other_id in ids:
        mutation_graph.add_edge(M_id,other_id)
"""        

def mutant_neighborhood_graph(link):
    G = Graph()
    isosig = link.exterior().isometry_signature(of_link=True)
    isosig_to_link_dict = {isosig: link.PD_code()}
    links_boundary = [link]
    while links_boundary:
        new_links_boundary = []
        for L in links_boundary:
            L_iso = L.exterior().isometry_signature(of_link=True)
            mutants = all_mutants(L)
            for M in mutants:
                isosig = M.exterior().isometry_signature(of_link=True)
                if isosig not in isosig_to_link_dict:
                    new_links_boundary.append(M)
                    isosig_to_link_dict[isosig] = M.PD_code()
                G.add_edge(L_iso, isosig)
        links_boundary = new_links_boundary

    return G, isosig_to_link_dict

mutation_graph = Graph()
isosig_dict = {}
i = 0
for link in links:
    print('%d / %d' % (i, len(links)))
    i += 1
    try:
        isosig = link.isometry_signature(of_link=True)
    except:
        print('isosig failed')
        continue
    if isosig in mutation_graph:
        continue
    G, d = mutant_neighborhood_graph(link.link())
    print(len(G.edges()))
    mutation_graph.add_edges(G.edges())
    for iso in d:
        isosig_dict[iso] = d[iso]

with open('mutation_graph_knots12.txt', 'w') as f:
    f.write(str(mutation_graph.edges()))

with open('isosig_to_PD_dict12.txt', 'w') as g:
    g.write(str(isosig_dict))

