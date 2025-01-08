import spherogram
import snappy
import numpy as np
import mpmath

from sage.all import PolynomialRing, LaurentPolynomialRing, RR, ZZ, RealField, ComplexField, matrix, arccos, exp


def is_palindromic(p):
    coeffs = p.coefficients(sparse=False)
    return coeffs == list(reversed(coeffs))

def compact_form(p):
    """
    Given a palindromic polynomial p(x) of degree 2n, return a
    polynomial g so that p(x) = x^n g(x + 1/x).
    """
    assert is_palindromic(p) and p.degree() % 2 == 0
    R = p.parent()
    x = R.gen()
    f, g = p, R(0)

    while f != 0:
        c = f.leading_coefficient()
        assert f.degree() % 2 == 0
        d = f.degree()//2
        g += c*x**d
        f = (f - c*(x**2 + 1)**d)
        if f != 0:
            e = min(f.exponents())
            assert e > 0 and f % x**e == 0
            f = f // x**e

    # Double check
    L = LaurentPolynomialRing(R.base_ring(), repr(x))
    y = L.gen()
    assert p == y**(p.degree() // 2) * g(y + 1/y)
    return g


def roots_on_unit_circle(poly, prec=53):
    """
    For a palindromic polynomial p(x) of even degree, return all the
    roots on the unit circle in the form

        (argument/(2 pi), multiplicity)

    """
    assert is_palindromic(poly) and poly.degree() % 2 == 0
    assert poly.parent().is_exact()
    # Deal with these corner cases at some point if needed.
    assert poly(1) != 0 and poly(-1) != 0

    ans = []
    RR = RealField(prec)
    pi = RR.pi()
    g = compact_form(poly)
    for f, e in g.factor():
        roots = [r for r in f.roots(RR, False) if -2 < r < 2]
        args = [arccos(r/2)/(2*pi) for r in roots]
        args += [1-a for a in args]
        ans += [ (arg, e) for arg in args]
    return sorted(ans)


def alexander_poly_from_seifert(V):
    R = PolynomialRing(ZZ, 't')
    t = R.gen()
    poly = (t*V - V.transpose()).determinant()
    if poly.leading_coefficient() < 0:
        poly = -poly
    e = min(poly.exponents())
    if e > 0:
        poly = poly // t**e
    return poly

def signature_via_numpy(A):
    CC = ComplexField(53)
    A = A.change_ring(CC).numpy()
    assert np.linalg.norm(A - A.conjugate().transpose()) < 10e-10
    eigs = np.linalg.eigh(A)[0]
    smallest = min(np.abs(eigs))
    assert smallest > 1e-5
    return np.sum(eigs > 0) - np.sum(eigs < 0)


def signature_function_of_integral_matrix(V, prec=53):
    """
    Computes the signature function sigma of V via numerical methods.
    Returns two lists, the first representing a partition of [0, 1]:

         x_0 = 0 < x_1 < x_2 < ... < x_n = 1

    and the second list consisting of the values [v_0, ... , v_(n-1)]
    of sigma on the interval (x_i, x_(i+1)).  Currently, the value of
    sigma *at* x_i is not computed.
    """
    poly = alexander_poly_from_seifert(V)
    RR = RealField(prec)
    CC = ComplexField(prec)
    pi = RR.pi()
    I = CC.gen()
    partition = [RR(0)] + [a for a, e in roots_on_unit_circle(poly, prec)] + [RR(1)]
    n = len(partition) - 1
    values = []
    for i in range(n):
        omega = exp((2*pi*I)*(partition[i] + partition[i + 1])/2)
        A = (1 - omega)*V + (1 - omega.conjugate())*V.transpose()
        values.append(signature_via_numpy(A))

    assert list(reversed(values)) == values
    return partition, values


def signature_function(L, prec=53):
    """
    Computes the signature function sigma of K via numerical methods.
    Returns two lists, the first representing a partition of [0, 1]:

         x_0 = 0 < x_1 < x_2 < ... < x_n = 1

    and the second list consisting of the values [v_0, ... , v_(n-1)]
    of sigma on the interval (x_i, x_(i+1)).  Currently, the value of
    sigma *at* x_i is not computed.
    """
    V = L.seifert_matrix()
    return signature_function_of_integral_matrix(V)


def basic_knot_test():
    # All passed!
    for M in snappy.HTLinkExteriors(cusps=1):
        print(M.name())
        R = PolynomialRing(ZZ, 't')
        K = M.link()
        V = matrix(K.seifert_matrix())
        p0 = R(M.alexander_polynomial())
        p1 = R(K.alexander_polynomial())
        p2 = alexander_poly_from_seifert(V)
        assert p0 == p1 == p2
        partition, values = signature_function_of_integral_matrix(V)
        n = len(values)
        assert n % 2 == 1
        m = (n - 1)//2
        assert K.signature() == values[m]


def invertible_representative(seifert_matrix):
    """
    For many calculations in the algebraic concordance group, one
    needs to replace the initial Seifert matrix with an invertible
    matrix that represents the same class.  See::

      Levine, Invariants of knot cobordism. Invent Math 8, 98–110 (1969).
      https://doi.org/10.1007/BF01404613

    or the summary in::

      Livingston, The algebraic concordance order of a knot.
      https://arxiv.org/abs/0806.3068

    To keep the sizes of the entries of the final matrix under
    control, we use LLL at various points.  We also use a "block
    matrix" form of Levine's lemma to reduce the amount of recursion.

    Warning: this code has not been extensively tested.
    """
    V = seifert_matrix
    n = V.nrows()
    assert V.base_ring() == ZZ
    assert (V - V.transpose()).det().abs() == 1
    assert (V + V.transpose()).rank() == n

    # Start with the sublattice L that is left-orthogonal to
    # everything, i.e. the v with v*V = 0.

    L = V.left_kernel_matrix()
    k = L.nrows()
    if k == 0:
        return V
    L = L.LLL()  # The basis is the rows.

    # Now look at the left-orthogonal complement of L, that is,
    #
    #    U = {u with u * V * v = 0 for all v in L}

    U = matrix(V * L.transpose()).left_kernel_matrix()
    assert U * V * L.transpose() == 0

    # L is a submodule of U, and we now find a basis of U that extends
    # the one for L.

    L_in_U = U.solve_left(L).change_ring(ZZ)
    E = extend_to_primitive(L_in_U).change_ring(ZZ)
    assert E.det().abs() == 1

    # Another basis for U, this time with the basis for L for the
    # first vectors

    U = E * U
    assert U[:k, :] == L
    # cleanup
    U = L.stack(U[k: , :].LLL())

    # Extend U to a Z-basis of Z^n where the new vectors give the
    # standard dual basis of Hom(L, Z) via a -> (x -> a*V*x).

    R = extend_to_primitive(U).change_ring(ZZ)
    A = R[U.nrows(): , :].LLL()
    F = A * V * L.transpose()
    assert F.det().abs() == 1
    A = F.inverse().change_ring(ZZ) * A

    C = L.stack(A).stack(U[k: , :])
    assert C.det().abs() == 1

    # We apply the block matrix generalization of the lemma in Levine's
    # Inventiones 1969 paper on page 101, to see that B is equivalent
    # in G^Q and hence G^Z. The block version follows from applying the
    # original repeatedly.

    W = C * V * C.transpose()
    assert W[:k, :] == 0 and W[k:2*k, 0:k] == 1 and W[2*k: , :k] == 0
    B = W[2*k:, 2*k:]
    assert (B - B.transpose()).det().abs() == 1
    return invertible_representative(B)


if __name__ == '__main__':
    K = spherogram.Link('K12n123')
    V = matrix(ZZ, K.seifert_matrix())
    basic_knot_test()
