from __future__ import print_function
from . import links, tangles
Crossing, Link, RationalTangle, IdentityBraid = links.Crossing, links.Link, tangles.RationalTangle, tangles.IdentityBraid
import os, sys, re

#----- Some basic tests, constructing links by hand -------

def figure8():
    a, b, c, d = [Crossing(x) for x in 'abcd']
    a[0] = c[1]
    a[1] = d[0]
    a[2] = b[1]
    a[3] = b[0]
    b[2] = d[3]
    b[3] = c[2]
    c[3] = d[2]
    c[0] = d[1]
    return Link([a,b,c,d])

def punct_torus():
    a = Crossing('a')
    a[0] = a[2]
    a[1] = a[3]
    return Link([a], check_planarity=False)

def whitehead():
    a, b, c, d, e =  crossings = [Crossing(x) for x in 'abcde']
    a[0] = b[3]
    a[1] = d[0]
    a[2] = d[3]
    a[3] = c[0]
    b[0] = c[3]
    b[1] = e[2]
    b[2] = e[1]
    c[1] = d[2]
    c[2] = e[3]
    d[1] = e[0]
    return Link(crossings)
 
def basic_test():
    K, W, T = figure8(), whitehead(), punct_torus()
    print( K.is_planar(), W.is_planar(), punct_torus().is_planar() )
    print( K.PD_code(True) )
    print( K.DT_code(True) , K.peer_code())
    print( W.PD_code(True) )
    print( W.DT_code(True) , K.peer_code())


# ----- checking the SnapPy link exterior code ------------

def knot(fractions):
    if len(fractions) == 1:
        return RationalTangle(*fractions[0]).denominator_closure()
    else:
        A, B, C = [RationalTangle(*f) for f in fractions]
        T = A + B + C
        return T.numerator_closure()

def some_knots():
    from . import hyperbolic_montesinos 
    return [ (K, knot(fractions)) for K, fractions in hyperbolic_montesinos.knots] 

def exterior_test():
    try:
        import snappy
    except ImportError:
        print("SnapPy not installed, skipping link exterior test.")

    print(figure8().exterior().volume(), whitehead().exterior().volume())

    C, Id = RationalTangle(1), IdentityBraid(1)
    x = C | Id
    y = Id | (-C)
    print((x*y*x*y).denominator_closure().exterior().volume())

    for name, K in some_knots():
        M0, M1 = K.exterior(), snappy.Manifold(K.DT_code(True))
        N0 = snappy.LinkExteriors.identify(M0)
        N1 = snappy.LinkExteriors.identify(M1)
        assert N0.name() == name and N1.name() == name

    print("Checked two different ways of building 167 hyperbolic knots.")


def alexander_polynomial_test():
    """
    Export a bunch of Montesinos knots as PD diagrams and
    compute the Alexander polynomials via KnotTheory and
    also Sage.  Make sure they match.
    """
    
    try:
        sys.path.append('/Users/dunfield/work/python')
        import nsagetools
    except ImportError:
        print("Skipping this test as you're not within Sage and/or are not Nathan.")
        return 

    from sage.all import mathematica, PolynomialRing, ZZ

    mathematica.execute('AppendTo[$Path, "/pkgs"]')
    mathematica.execute('<<KnotTheory`')
    mathematica.execute('Off[KnotTheory::loading]')
    mathematica.execute('Off[KnotTheory::credits]')
    mathematica.execute('MyAlex[L_] := CoefficientList[t^100*Alexander[L][t], t]')

    def alex_by_KnotTheory(L):
        p = mathematica.MyAlex(L.PD_code(True)).sage()
        i = min( [i for i,c in enumerate(p) if c != 0])
        R = PolynomialRing(ZZ, 'a')
        return R(p[i:])

    def alex_match(L):
        p1 = alex_by_KnotTheory(L)
        p2 = nsagetools.alexander_polynomial(L.exterior())
        return 0 in [p1 - p2, p1 + p2]

    for name, K in some_knots():
        assert alex_match(K)

    print("Checked Alexander polynomials via KnotTheory for these 167 knots.")

if __name__ == '__main__':
    basic_test()
    exterior_test()
    alexander_polynomial_test()
