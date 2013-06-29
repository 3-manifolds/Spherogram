"""
The *bridge number* of a planar link diagram D is

b(D) = min { # of maxima of h on D }

where h is a generic height function on R^2. Here, we compute b(D) via
linear programming, using ideas from

[DP] Didimo and Pizzonia, Upward Embeddings and Orientations of
Undirected Planar Graphs. <http://jgaa.info/getPaper?id=68>

To each corner of a face of D is classified into "large", "flat", or
"small" (which is the default) as per [DP], which correspond to integers
2, 1, and 0 respectively.  The geometric conditions are

* Every vertex has either a single "large" angle or two "flat" ones which
are opposite.

* The sum of the corner types around every face is degree - 2, with the
exception of the exterior face (which is degree + 2).
"""

import spherogram, snappy, random
from spherogram import DTcodec, RationalTangle
from spherogram.links.links import CrossingStrand



def bridge_LP(link):
    """
    An integer linear program which computes the bridge number of the given
    link diagram.   
    """
    LP = MixedIntegerLinearProgram(maximization=False)

    # Near a crossing, the level sets of the height function are either
    # horizontal or vertical, and so there are two flat corners which
    # must be opposite.  
    hor_cross = LP.new_variable(binary=True)
    vert_cross = LP.new_variable(binary=True)

    # We add a dummy vertex in the middle of each edge so it can bend.
    # The corners here must either both be flat or just one is large.  For
    # the corner to large, *both* flat_edge and large_edge must be 1.  
    flat_edge = LP.new_variable(binary=True)
    large_edge = LP.new_variable(binary=True)

    # Exactly one complementary region is the exterior one.
    exterior = LP.new_variable(binary=True)
    faces = link.faces()
    LP.add_constraint(sum(exterior[i] for i, F in enumerate(faces)) == 1)

    # Now for the main constraints 
    for c in link.crossings:
        LP.add_constraint(hor_cross[c] + vert_cross[c] == 1)
        for ce in c.entry_points():
            s = CrossingStrand(c, ce.entry_point)
            t = s.opposite()
            LP.add_constraint( flat_edge[s] == flat_edge[t] )
            LP.add_constraint( flat_edge[s] + large_edge[s] + large_edge[t] == 1 )

    for i, face in enumerate(faces):
        eqn = 0
        for cs in face:
            flat = hor_cross if cs.entry_point % 2 == 0 else vert_cross
            eqn += flat[cs.crossing] + flat_edge[cs] + 2*large_edge[cs]
        LP.add_constraint(eqn == (2*len(face) - 2 + 4*exterior[i]))

    LP.set_objective(sum(large_edge.values()))
    bridge = int(LP.solve())
    assert bridge % 2 == 0
    return bridge//2, LP.get_values([hor_cross, vert_cross, flat_edge, large_edge, exterior])



def test1():
    for M in snappy.HTLinkExteriors(num_cusps=1):
        if M.volume() < 2.0:
            continue
        L = link_from_manifold(M)
        b = bridge_LP(L)[0]
        tb = M.is_two_bridge()
        if b == 2:
            assert tb
        if tb and b > 2 or b > 3:
            print M, b, M.is_two_bridge()

unknot = RationalTangle(1).numerator_closure()
hopf = RationalTangle(2).numerator_closure()
trefoil = DTcodec([(4,6,2)]).link()
big_knot = DTcodec([(4, 12, 14, 22, 20, 2, 28, 24, 6, 10, 26, 16, 8, 18)]).link()
big_link = DTcodec([(8, 12, 16), (18, 22, 24, 20), (4, 26, 14, 2, 10, 6)]).link()


def link_from_manifold(manifold):
    return DTcodec(manifold.DT_code()).link()

def random_16():
    n = random.randrange(1, 1008907)
    return snappy.Manifold('16n%d' % n)

def test_16(N):
    for i in xrange(N):
        M = random_16()
        L = link_from_manifold(random_16())
        print M.name(), bridge_LP(L)[0]
