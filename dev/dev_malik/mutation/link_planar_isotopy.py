import spherogram
from spherogram.links.planar_isotopy import link_isosig
from spherogram.links.random_links import map_to_link, random_map
from random import choice
from collections import Counter


def random_rooted_link(size, edge_conn=2):
    K = map_to_link(random_map(size, edge_conn))
    root = choice(K.crossing_strands())
    return K, root


def isosig(self, root=None, over_or_under=False):
    return link_isosig(self, root, over_or_under)


spherogram.Link.isosig = isosig

num_samples = 600
size = 4
edge_conn = 4

isosigs = []
isosig_to_link = {}
for i in range(num_samples):
    print(i)
    K, root = random_rooted_link(size,edge_conn)
    isosig = K.isosig(root,over_or_under=True)
    isosig_to_link[isosig]=(K,root)
    isosigs.append(isosig)

C = Counter(isosigs)
print(C.most_common(10))
