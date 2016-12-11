"""
Generating random knots starting with randomly generated planar maps.

Code written as part of a 2013 Illinois Geometry Lab project by

Nathan Dunfield, Alexander Ehrenberg, Shiladitya Bhattacharyya, and Dongming Lei

Details to hopefullly appear in some paper or other.  I wouldn't hold
your breath, though.
"""

import os, sys, re, gzip, random
from .. import graphs
from . import links, twist
from spherogram.planarmap import random_map as raw_random_map


class LinkGenerationError(Exception):
    def __init__(self, num_tries):
        message = "Didn't generate requested link after %d tries." % num_tries
        Exception.__init__(self, message)

def random_map(num_verts, edge_conn_param=4,
               num_link_comps=0, max_tries=100):
    """
    Returns a dictionary of endpoints of edges in the form:

    (signed edge) -> (vertex, position)
    """
    if isinstance(num_verts, list):
        data = num_verts
    else:
        data = raw_random_map(num_verts, edge_conn_param,
                              num_link_comps, max_tries)
        if data is None:
            raise LinkGenerationError(max_tries)
        

    vertex_adjacencies = []
    for vertex, adjacencies in data:
        if random.randrange(2):
            adjacencies = adjacencies[1:] + adjacencies[:1]
        vertex_adjacencies.append(adjacencies)

    edge_adjacencies = dict()
    for v, edges in enumerate(vertex_adjacencies):
        for i, edge in enumerate(edges):
            edge_adjacencies[edge] = (v, i)
    return edge_adjacencies
        
def map_to_link(map):
    num_edges = len(map) // 2
    crossings = [links.Crossing() for i in range(num_edges//2)]
    for e in range(1, num_edges+1):
        (a, i), (b,j) = map[e], map[-e]
        crossings[a][i] = crossings[b][j]
    ans = links.Link(crossings, check_planarity=False)
    return ans

def self_crossings(component):
    comp_set = set(component)
    return  [ce for ce in component if ce.other() in comp_set]
    
def longest_components(link, num_components):
    self_crosses = [(len(self_crossings(comp)), comp) for comp in link.link_components]
    self_crosses.sort(reverse=True)
    return [x[1] for x in self_crosses[:num_components]]

def simplified_prime_pieces(link, simplify_fun):
    ans = []
    cur = link.deconnect_sum(True)
    while len(cur):
        L = cur.pop()
        if simplify_fun(L):
            cur += L.deconnect_sum(True)
        else:
            ans.append(L)
    return ans

def largest_prime_piece(link, simplify_fun):
    pieces = simplified_prime_pieces(link, simplify_fun)
    if len(pieces) == 0:
        return links.Link([])
    return max(pieces, key=lambda L:len(L.crossings))


def random_link(crossings,
                num_components = 'any', 
                initial_map_gives_link = False,
                alternating = False,
                consistent_twist_regions = False,
                simplify = 'basic',
                prime_decomposition = True,
                return_all_pieces = False, 
                max_tries=100):
    """
    Generates a random link from a model that starts with a random
    4-valent planar graph sampled with the uniform distribution by
    Schaeffer's `PlanarMap program.
    <http://www.lix.polytechnique.fr/~schaeffe/PagesWeb/PlanarMap/index-en.html>`_ 

    The ``crossings`` argument specifies the number of vertices of the
    initial planar graph G; the number of crossing in the returned knot
    will typically be less. The meanings of the optional arguments are as
    follows:

    1. ``num_components``: The number of components of the returned link.
       The link naively associated to G may have too few or too many
       components. The former situation is resolved by picking another G,
       and the latter by either

       a. Taking the sublink consisting of the components with the largest
          self-crossing numbers.

       b. Resampling G until the desired number of components is achieved;
          this can take a very long time as the expected number of
          components associated to G grows linearly in the number of
          vertices.

       When the argument ``initial_map_gives_link`` is ``False`` the
       program does (a) and this is the default behavior. If you want (b)
       set this argument to ``True``.

       To get the entire link associated to G, set ``num_components`` to
       ```any```, which is also the default.

    2. The 4-valent vertices of G are turned into crossings by flipping a
       fair coin. If you want the unique alternating diagram associated to
       G, pass ``alternating = True``.  If you want there to be no
       obvious Type II Reidemeister moves, pass
       ``consistent_twist_regions = False``.

    3. ``simplify``: Whether and how to try to reduce the number of
       crossings of the link via Reidemeister moves using the method
       ``Link.simplify``.  For no simplification, set ``simplify = None``;
       otherwise set ``simplify`` to be the appropriate mode for
       ``Link.simplify``, for example ``basic`` (the default), ``level``,
       or ``global``.

    4. ``prime_decomposition``:  The initial link generated from G may not
       be prime (and typically isn't if ``initial_map_gives_link`` is
       ``False``). When set (the default), the program undoes any connect
       sums that are "diagrammatic obvious", simplifies the result, and
       repeats until pieces are "diagrammatically prime".  If
       ``return_all_pieces`` is ``False`` (the default) then only the
       largest (apparently) prime component is returned; otherwise all
       summands are returned as a list.


    Some examples:
    
    >>> K = random_link(25, num_components=1, initial_map_gives_link=True, alternating=True)
    >>> K
    <Link: 1 comp; 25 cross>

    >>> L= random_link(30, consistent_twist_regions=True, simplify = 'global')
    >>> isinstance(random_link(30, return_all_pieces=True), list)
    True
    """

    # This means no trivial loops.  PlanarMap accepts 6, which means
    # no bigons, but this is unbearably slow.  
    edge_conn_param = 4

    # Generate the initial link    
    if num_components == 'any':
        plane_map = random_map(crossings, edge_conn_param)
        link = map_to_link(plane_map)
    elif initial_map_gives_link:
        plane_map = random_map(crossings, edge_conn_param,
                               num_components, max_tries)
        link = map_to_link(plane_map)
    else:
        for i in range(max_tries):
            plane_map = random_map(crossings, edge_conn_param)
            link = map_to_link(plane_map)
            if len(link.link_components) >= num_components:
                break
        
        comps = longest_components(link, num_components)
        link = link.sublink(comps)
        if len(link.link_components) != num_components:
            raise LinkGenerationError(max_tries)


    # Adjust the currently random crossings to match the request

    if alternating:
        link = link.alternating()
        return link

    if consistent_twist_regions:
        twist.make_twist_regions_consistent(link)

    # Initial simplification, if any.

    def simplify_func(link):
        if isinstance(simplify, dict):
            return link.simplify(**simplify)
        elif isinstance(simplify, str):
            return link.simplify(mode=simplify)
        return False
            
    simplify_func(link)
            
    # Pull into "prime" pieces, if requested.  
    
    if prime_decomposition:
        if return_all_pieces:
            return simplified_prime_pieces(link, simplify_func)
        else:
            return largest_prime_piece(link, simplify_func)
    else:
        return link 

#def random_knot(crossings, **kwargs):
#    return random_link(crossings, num_components=1, **kwargs)
#
#def new_random_knot(crossings, prime_decomposition=True):
#    return random_knot(crossings, edge_conn_param=4,
#                initial_map_gives_knot=True, 
#                return_all_pieces=False, 
#                prime_decomposition=prime_decomposition,
#                consistent_twist_regions=True)

