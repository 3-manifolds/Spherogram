"""

Computing the Jones polynomial following the conventions of Jake
Rasmussen's lectures at PCMI, Summer 2019:

https://www.dpmms.cam.ac.uk/~jar60/PCMINotes.pdf
https://www.youtube.com/watch?v=--l-XOhDXOU
https://www.youtube.com/watch?v=3VwcGHycyAE




"""

from ..sage_helper import _within_sage
from . import exhaust
from .links_base import Link

if _within_sage:
    from sage.all import ZZ, LaurentPolynomialRing, PerfectMatchings, PerfectMatching
    R = LaurentPolynomialRing(ZZ, 'q')
    q = R.gen()
else:
    pass


def num_Pn(n):
    """
    An element of Jake's P_{0, n} of planar tangles from 0 points to n
    points will be a PerfectMatching of [0,..,n-1] where
    is_noncrossing is True.
    """
    return len([m for m in PerfectMatchings(n) if m.is_noncrossing()])


def insert_cup(matching, i):
    """
    Insert a new adjacent matching which joins i and i + 1.

    sage: m = PerfectMatching([(0, 1), (2, 5), (3, 4)])
    sage: insert_cup(m, 0)
    [(0, 1), (2, 3), (4, 7), (5, 6)]
    sage: insert_cup(m, 1)
    [(0, 3), (1, 2), (4, 7), (5, 6)]
    sage: insert_cup(m, 6)
    [(0, 1), (2, 5), (3, 4), (6, 7)]
    """
    assert len(matching.base_set()) >= i

    def shift(pair):
        return [a if a < i else a + 2 for a in pair]
    return PerfectMatching([shift(pair) for pair in matching] + [(i, i + 1)])


def cap_off(matching, i):
    """
    Merge i and i + 1 with a cap.  Returns a new matching and whether
    or not a circle was created.

    sage: m = PerfectMatching([(0, 5), (1, 4), (2, 3)])
    sage: cap_off(m, 2)
    ([(0, 3), (1, 2)], True)
    sage: cap_off(m, 3)
    ([(0, 3), (1, 2)], False)
    """
    def shift(a):
        return a if a < i else a - 2

    def follow(a):
        b = matching.partner(a)
        if b == i:
            b = matching.partner(i + 1)
        elif b == i + 1:
            b = matching.partner(i)
        return b

    n = len(matching.base_set())
    circle = matching.partner(i) == i + 1
    new_match = set()
    for a in range(n):
        if a not in {i, i + 1}:
            u, v = shift(a), shift(follow(a))
            if u > v:
                u, v = v, u
            new_match.add((u, v))
    return PerfectMatching(new_match), circle


class VElement():
    """
    An element of some V_{0, n} which is the free R-module on P_{0, n}

    sage: m = PerfectMatching([(0, 1), (3, 4), (2, 5)])
    sage: v1 = VElement(m)
    sage: v1
    (1)*[(0, 1), (2, 5), (3, 4)]
    sage: v2 = (q + q**-1)*v1
    sage: v2
    (q^-1 + q)*[(0, 1), (2, 5), (3, 4)]
    sage: v3 = q* VElement(PerfectMatching([(5, 0), (4, 3), (1, 2)]))
    sage: v1 + v2 + v3
    (q^-1 + 1 + q)*[(0, 1), (2, 5), (3, 4)] + (q)*[(0, 5), (1, 2), (3, 4)]
    sage: v2.insert_cup(6)
    (q^-1 + q)*[(0, 1), (2, 5), (3, 4), (6, 7)]
    sage: (v1 + v2 + v3).cap_off(1)
    (q^-1 + 2 + q + q^2)*[(0, 3), (1, 2)]
    """
    def __init__(self, spec=None):
        self.dict = dict()
        if spec is None:
            spec = PerfectMatching([])
        if isinstance(spec, dict):
            self.dict = spec
        if isinstance(spec, PerfectMatching):
            assert spec.is_noncrossing()
            self.dict[spec] = R.one()

    def __rmul__(self, other):
        if other in R:
            other = R(other)
            return VElement({m: other * c for m, c in self.dict.items()})

    def __add__(self, other):
        if isinstance(other, VElement):
            ans_dict = self.dict.copy()
            for matching, coeff in other.dict.items():
                cur_coeff = self.dict.get(matching, R.zero())
                ans_dict[matching] = cur_coeff + coeff
        return VElement(ans_dict)

    def __repr__(self):
        matchings = sorted(self.dict)
        terms = ['(%s)*%s' % (self.dict[m], m) for m in matchings]
        if len(terms) == 0:
            return '0'
        return ' + '.join(terms)

    def insert_cup(self, i):
        """
        Insert an new matching at (i, i + 1)
        """
        return VElement({insert_cup(m, i): c for m, c in self.dict.items()})

    def cap_off(self, i):
        ans_dict = dict()
        for matching, coeff in self.dict.items():
            new_matching, has_circle = cap_off(matching, i)
            cur_coeff = ans_dict.get(new_matching, R.zero())
            if has_circle:
                coeff = (q + q**-1) * coeff
            ans_dict[new_matching] = cur_coeff + coeff
        return VElement(ans_dict)

    def cap_then_cup(self, i):
        return self.cap_off(i).insert_cup(i)

    def add_positive_crossing(self, i):
        return self + (-q) * self.cap_then_cup(i)

    def add_negative_crossing(self, i):
        return self.cap_then_cup(i) + (-q) * self

    def is_multiple_of_empty_pairing(self):
        return len(self.dict) == 1 and (PerfectMatching([]) in self.dict)


def kauffman_bracket(link):
    """
    sage: L = Link('T(2, 3)')
    sage: kauffman_bracket(L)
    q^-2 + 1 + q^2 - q^6

    sage: U4 = Link(braid_closure=[1, -1, 2, -2, 3, -3])
    sage: kauffman_bracket(U4)
    -q^-1 - 4*q - 6*q^3 - 4*q^5 - q^7

    sage: U3 = Link([])
    sage: U3.unlinked_unknot_components = 3
    sage: kauffman_bracket(U3) == (q + q**-1)**3
    True
    """
    ans = VElement()
    if isinstance(link, exhaust.MorseEncoding):
        encoded = link
    else:
        exhaustion = exhaust.MorseExhaustion(link)
        encoded = exhaust.MorseEncoding(exhaustion)
    for event in encoded:
        if event.kind == 'cup':
            ans = ans.insert_cup(event.min)
        elif event.kind == 'cap':
            ans = ans.cap_off(event.min)
        else:
            assert event.kind == 'cross'
            if event.a < event.b:
                ans = ans.add_positive_crossing(event.min)
            else:
                ans = ans.add_negative_crossing(event.min)
    assert ans.is_multiple_of_empty_pairing()
    return ans.dict[PerfectMatching([])]


def jones_polynomial(link, normalized=True):
    bracket = kauffman_bracket(link)
    if normalized:
        factor = q + q**-1
        norm_bracket = bracket // factor
        assert norm_bracket * factor == bracket
    else:
        norm_bracket = bracket

    signs = [c.sign for c in link.crossings]
    n_minus, n_plus = signs.count(-1), signs.count(1)
    return (-1)**n_minus * q**(n_plus - 2 * n_minus) * norm_bracket


def test_links(N):
    import snappy
    for i in range(N):
        M = snappy.HTLinkExteriors.random()
        L = M.link()
        print(M.name(), test_one_link(L))


if __name__ == '__main__':
    import doctest
    print(doctest.testmod())
