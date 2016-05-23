import spherogram, snappy
from sage.all import Link as SageLink
from sage.all import LaurentPolynomialRing, PolynomialRing, ZZ, var
from sage.symbolic.ring import SymbolicRing

def alexander_poly_of_sage(knot):
    p = knot.alexander_polynomial()
    ans = p.polynomial_construction()[0]
    if ans.leading_coefficient() < 0:
        ans = -ans
    return ans

def jones_polynomial_of_sage(knot):
    """
    To match our conventions (which seem to agree with KnotAtlas), we
    need to swap q and 1/q.
    """
    q = LaurentPolynomialRing(ZZ, 'q').gen()
    p = knot.jones_polynomial(skein_normalization=True)
    exps = p.exponents()
    assert all(e % 4 == 0 for e in exps)
    return sum(c*q**(-e//4) for c, e in zip(p.coefficients(), exps))

def test_knot(snappy_manifold):
    M = snappy_manifold
    assert M.num_cusps() == 1

    U = M.link()
    T = U.sage_link()

    assert [list(x) for x in U.PD_code(min_strand_index=1)] == T.pd_code()
    assert U.alexander_polynomial() == alexander_poly_of_sage(T)
    assert U.signature() == T.signature()
    assert U.jones_polynomial() == jones_polynomial_of_sage(T)

    T_alt = SageLink(U.braid_word(as_sage_braid=True))
    assert U.signature() == T_alt.signature()
    U_alt = spherogram.Link(T.braid())
    assert U.signature() == U_alt.signature()
    
    
    #q = var('q')
    #print (SR(K_us.jones_poly()) - K_them.jones_polynomial(variab='q').subs(q=1/q)).simplify()
    #print K_us.signature(), K_them.signature()


    
M = snappy.Manifold('K11n42')
knot_manifolds = snappy.HTLinkExteriors(cusps=1)

