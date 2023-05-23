"""
When used within Sage, the Link class gains many methods to compute
standard invariants.  Much of this code was contributed by Robert Lipshitz
and Jennet Dickinson.
"""

from . import links_base, alexander
from .links_base import CrossingStrand, Crossing
from ..sage_helper import _within_sage, sage_method, sage_pd_clockwise

deprecation_warnings_issued = set()

if _within_sage:
    from sage.matrix.constructor import matrix
    from sage.groups.free_group import FreeGroup
    import sage.graphs.graph as graph
    from sage.groups.braid import Braid, BraidGroup
    from sage.rings.rational_field import QQ
    from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
    from sage.quadratic_forms.quadratic_form import QuadraticForm
    try:
        from sage.knots.knot import Knot as SageKnot
        from sage.knots.link import Link as SageLink
    except ImportError:  # Sage older than 7.2
        SageKnot, SageLink = type(None), type(None)
    import functools
else:
    pass


def normalize_alex_poly(p, t):
    """
    Normalize the sign of the leading coefficient and make all
    exponents positive, then return as an ordinary rather than Laurent
    polynomial.
    """
    if len(t) == 1:
        p = p * (t[0]**(-min(p.exponents())))
        if p.coefficients()[-1] < 0:
            p = -p
        p, e = p.polynomial_construction()
        assert e == 0
        return p

    max_degree = max(sum(x) for x in p.exponents())
    highest_monomial_exps = [x for x in p.exponents() if sum(x) == max_degree]
    leading_exponents = max(highest_monomial_exps)
    leading_monomial = functools.reduce(lambda x, y: x * y,
                                        [t[i]**(leading_exponents[i])
                                         for i in range(len(t))])
    l = p.monomial_coefficient(leading_monomial)

    if l < 0:
        p = -p

    for i, ti in enumerate(t):
        min_exp = min(x[i] for x in p.exponents())
        p = p * (ti**(-min_exp))

    R = p.parent()
    p = R.polynomial_ring()(p)
    return p


def sage_braid_as_int_word(braid):
    """
    Convert a Sage Braid to a word.
    """
    # Could simplify using braid.Tietze().

    G = braid.parent()
    gen_idx = {g: i + 1 for i, g in enumerate(G.gens())}
    ans = []
    for g, e in braid.syllables():
        if e > 0:
            ans += e * [gen_idx[g]]
        else:
            ans += abs(e) * [-gen_idx[g]]

    n = G.ngens()
    if max(abs(a) for a in ans) < n:
        ans += [n, -n]
    return ans


extra_docstring = """
    You can also convert to and from SageMath braid and link types,
    see the documentation for the "sage_link" method for details.
"""


class Link(links_base.Link):
    __doc__ = links_base.Link.__doc__ + extra_docstring

    def __init__(self, crossings=None, braid_closure=None, check_planarity=True, build=True):
        if _within_sage:
            if isinstance(crossings, Braid):
                assert braid_closure is None
                braid_closure = crossings
                crossings = None
            if isinstance(braid_closure, Braid):
                braid_closure = sage_braid_as_int_word(braid_closure)
            if crossings is not None and isinstance(crossings, (SageKnot, SageLink)):
                # Sage's PD_code lists strands *clockwise* not our
                # *anticlockwise* prior to Sage 10.1.
                crossings = crossings.pd_code()
                if sage_pd_clockwise:
                    crossings = [[x[0], x[3], x[2], x[1]] for x in crossings]

        links_base.Link.__init__(self, crossings, braid_closure, check_planarity, build)

    @sage_method
    def linking_matrix(self):
        """
        Calculates the linking number for each pair of link components.

        Returns a linking matrix, in which the (i,j)th component is the
        linking number of the ith and jth link components.
        """
        mat = [[0 for i in range(len(self.link_components))]
               for j in range(len(self.link_components))]
        for n1, comp1 in enumerate(self.link_components):
            for n2, comp2 in enumerate(self.link_components):
                tally = [[0 for m in range(len(self.crossings))]
                         for n in range(2)]
                if comp1 != comp2:
                    for i, c in enumerate(self.crossings):
                        for x1 in comp1:
                            if x1[0] == c:
                                tally[0][i] += 1
                        for x2 in comp2:
                            if x2[0] == c:
                                tally[1][i] += 1
                for k, c in enumerate(self.crossings):
                    if (tally[0][k] == 1 and tally[1][k] == 1):
                        mat[n1][n2] += 0.5 * (c.sign)
                mat[n1][n2] = int(mat[n1][n2])
        return mat

    @sage_method
    def knot_group(self):
        """
        Computes the knot group using the Wirtinger presentation.

        Returns a finitely presented group::

           sage: K = Link('3_1')
           sage: G = K.knot_group()
           sage: type(G)
           <class 'sage.groups.finitely_presented.FinitelyPresentedGroup_with_category'>
        """
        n = len(self.crossings)
        F = FreeGroup(n)
        rels = []
        pieces = self._pieces()

        for z in self.crossings:
            for m, p in enumerate(pieces):
                for t, q in enumerate(p):
                    if q[0] == z:
                        if t == 0:
                            j = m
                        elif t == len(p) - 1:
                            i = m
                        else:
                            k = m
            i += 1
            j += 1
            k += 1
            if z.sign > 0:
                r = F([-k, i, k, -j])
            if z.sign < 0:
                r = F([k, i, -k, -j])
            rels.append(r)

        return F / rels

    @sage_method
    def alexander_matrix(self, mv=True):
        """
        Returns the Alexander matrix of the link::

            sage: L = Link('3_1')
            sage: A = L.alexander_matrix()
            sage: A                           # doctest: +SKIP
            ([   -1     t 1 - t]
            [1 - t    -1     t]
            [    t 1 - t    -1], [t, t, t])

            sage: L = Link([(4,1,3,2),(1,4,2,3)])
            sage: A = L.alexander_matrix()
            sage: A                           # doctest: +SKIP
            ([      -1 + t1^-1 t1^-1*t2 - t1^-1]
            [t1*t2^-1 - t2^-1       -1 + t2^-1], [t2, t1])
        """
        comp = len(self.link_components)
        if comp < 2:
            mv = False

        G = self.knot_group()
        num_gens = len(G.gens())

        L_g = LaurentPolynomialRing(QQ, [f'g{i+1}' for i in range(num_gens)])
        g = list(L_g.gens())

        if mv:
            L_t = LaurentPolynomialRing(QQ, [f't{i+1}' for i in range(comp)])
            t = list(L_t.gens())

            # determine the component to which each variable corresponds
            g_component = [c.strand_components[2] for c in self.crossings]
            for i, gci in enumerate(g_component):
                g[i] = t[gci]

        else:
            L_t = LaurentPolynomialRing(QQ, 't')
            t = L_t.gen()
            g = [t] * len(g)

        B = G.alexander_matrix(g)

        return (B, g)

    @sage_method
    def alexander_poly(self, *args, **kwargs):
        """
        Please use the "alexander_polynomial" method instead.
        """
        if 'alexander_poly' not in deprecation_warnings_issued:
            deprecation_warnings_issued.add('alexander_poly')
            print('Deprecation Warning: use "alexander_polynomial" instead of "alexander_poly".')
        return self.alexander_polynomial(*args, **kwargs)

    @sage_method
    def alexander_polynomial(self, multivar=True, v='no', method='default',
                             norm=True, factored=False):
        """
        Calculates the Alexander polynomial of the link.

        For links with one component,
        can evaluate the alexander polynomial at v::

            sage: K = Link('4_1')
            sage: K.alexander_polynomial()
            t^2 - 3*t + 1
            sage: K.alexander_polynomial(v=[4])
            5

            sage: K = Link('L7n1')
            sage: K.alexander_polynomial(norm=False)
            t1^-1*t2^-1 + t1^-2*t2^-4

        The default algorithm for *knots* is Bar-Natan's super-fast
        tangle-based algorithm.  For links, we apply Fox calculus to a
        Wirtinger presentation for the link::

            sage: L = Link('K13n123')
            sage: L.alexander_polynomial() == L.alexander_polynomial(method='wirtinger')
            True
        """
        # sign normalization still missing, but when "norm=True" the
        # leading coefficient with respect to the first variable is made
        # positive.
        if method == 'snappy':
            try:
                return self.exterior().alexander_polynomial()
            except ImportError:
                raise RuntimeError('the method "snappy" for '
                                   'alexander_polynomial requires SnapPy')

        comp = len(self.link_components)
        if comp < 2:
            multivar = False

        # If single variable, use the super-fast method of Bar-Natan.
        if comp == 1 and method == 'default' and norm:
            p = alexander.alexander(self)
        else:  # Use a simple method based on the Wirtinger presentation.
            if method not in ['default', 'wirtinger']:
                raise ValueError("Available methods are 'default' and 'wirtinger'")

            if multivar:
                L = LaurentPolynomialRing(QQ, [f't{i+1}' for i in range(comp)])
                t = list(L.gens())
            else:
                L = LaurentPolynomialRing(QQ, 't')
                t = [L.gen()]

            M = self.alexander_matrix(mv=multivar)
            C = M[0]
            m = C.nrows()
            n = C.ncols()
            if n > m:
                k = m - 1
            else:
                k = n - 1

            subMatrix = C[0: k, 0: k]
            p = subMatrix.determinant()
            if p == 0:
                return 0
            if multivar:
                t_i = M[1][-1]
                p = (p.factor()) / (t_i - 1)
                p = p.expand()

            if norm:
                p = normalize_alex_poly(p, t)

        if v != 'no':
            return p(*v)

        if multivar and factored:  # it's easier to view this way
            return p.factor()
        return p

    def knot_floer_homology(self, prime=2, complex=False):
        """
        Uses Zoltán Szabó's HFK Calculator to compute the knot Floer
        homology.  This also gives the Seifert genus, whether the knot
        fibers, etc:

        >>> K = Link('K3a1')
        >>> K.knot_floer_homology()    # doctest: +NORMALIZE_WHITESPACE
        {'L_space_knot': True,
         'epsilon': 1,
         'fibered': True,
         'modulus': 2,
         'nu': 1,
         'ranks': {(-1, -2): 1, (0, -1): 1, (1, 0): 1},
         'seifert_genus': 1,
         'tau': 1,
         'total_rank': 3}

        The homology itself is encoded by 'ranks', with the form::

          (Alexander grading, Maslov grading): dimension

        For example, here is the Conway knot, which has Alexander
        polynomial 1 and genus 3:

        >>> L = Link('K11n34')
        >>> ranks = L.knot_floer_homology()['ranks']
        >>> [(a, m) for a, m in ranks if a == 3]
        [(3, 3), (3, 4)]
        >>> ranks[3, 3], ranks[3, 4]
        (1, 1)

        Computation is done over F_2 by default, other primes less
        than 2^15 can be used instead via the optional "prime"
        parameter.

        If the parameter `complex` is set to True, then the simplified
        "UV = 0" knot Floer chain complex is returned. This complex is
        computed over the ring F[U,V]/(UV = 0), where F is the integers
        mod the chosen prime; this corresponds to only the horizontal and
        vertical arrows in the full knot Floer complex. The complex is
        specified by:

        * generators: a dictionary from the generator names to their
          (Alexander, Maslov) gradings.  The number of generators is
          equal to the total_rank.

        * differential: a dictionary whose value on (a, b) is an integer
          specifying the coefficient on the differential from generator a
          to generator b, where only nonzero differentials are
          recorded. (The coefficient on the differential is really an
          element of F[U,V]/(UV = 0), but the power of U or V can be
          recovered from the gradings on a and b so only the element of F
          is recorded.)

        For example, to compute the vertical differential, whose homology
        is HFhat(S^3), you can do:

        >>> data = L.knot_floer_homology(prime=31991, complex=True)
        >>> gens, diff = data['generators'], data['differentials']
        >>> vert = {(i,j):diff[i, j] for i, j in diff
        ...                          if gens[i][1] == gens[j][1] + 1}

        sage: from sage.all import matrix, GF
        sage: M = matrix(GF(31991), len(gens), len(gens), vert, sparse=True)
        sage: M*M == 0
        True
        sage: M.right_kernel().rank() - M.rank()
        1
        """
        import knot_floer_homology
        if len(self.link_components) + self.unlinked_unknot_components > 1:
            raise ValueError('Only works for knots, this has more components')
        if len(self.link_components) == 0 and self.unlinked_unknot_components == 1:
            return Link(braid_closure=[1, 1, -1]).knot_floer_homology()
        return knot_floer_homology.pd_to_hfk(self, prime=prime, complex=complex)

    @sage_method
    def black_graph(self):
        """
        Returns the black graph of K.

        If the black graph is disconnected (which can only happen for
        a split link diagram), returns one connected component. The
        edges are labeled by the crossings they correspond to.

        Example::

            sage: K=Link('5_1')
            sage: K.black_graph()
            Subgraph of (): Multi-graph on 2 vertices

        WARNING: While there is also a "white_graph" method, it need
        not be the case that these two graphs are complementary in the
        expected way.
        """
        faces = []
        for x in self.faces():
            l = []
            for y in x:
                l.append((y[0], y[1]))
                l.append((y[0], (y[1] + 1) % 4))
            faces.append(l)

        coords = list()
        for i in range(len(faces) - 1):
            for j in range(i + 1, len(faces)):
                a = set(faces[i])
                b = set(faces[j])
                s = a.union(b)
                for x in range(len(self.crossings)):
                    total = set(self.crossings[x][i] for i in range(4))
                    if total.issubset(s):
                        coords.append((tuple(faces[i]), tuple(faces[j]),
                                       self.crossings[x]))  # label by the crossing.

        G = graph.Graph(coords, multiedges=True)
        component = G.connected_components(sort=False)[1]
        return G.subgraph(component)

    @sage_method
    def white_graph(self):
        """
        Return the white graph of a non-split link projection.

        This method generates a multigraph whose vertices correspond
        to the faces of the diagram, with an edge joining two
        vertices whenever the corresponding faces contain opposite
        corners at some crossing.  To avoid hashability issues, the
        vertex corresponding to a face is the index of the face in the
        list returned by Link.faces().

        According to the conventions of "Gordon, C. McA. and
        Litherland, R. A, 'On the signature of a link', Inventiones
        math. 47, 23-69 (1978)", in a checkerboard coloring of a link
        diagram the unbounded region is always the first white region.
        Of course, the choice of which region is unbounded is
        arbitrary; it is just a matter of which region on S^2 contains
        the point at infinity.  In this method an equivalent arbitrary
        choice is made by just returning the second component of the
        multigraph, as determined by Graph.connected_components().
        (Empirically, the second component tends to be smaller than
        the first.)

        Note that this may produce a meaningless result in the case of
        a split link diagram.  Consequently if the diagram is split,
        i.e if the multigraph has more than 2 components, a ValueError
        is raised::

            sage: K=Link('5_1')
            sage: K.white_graph()
            Subgraph of (): Multi-graph on 2 vertices

        WARNING: While there is also a "black_graph" method, it need
        not be the case that these two graphs are complementary in the
        expected way.
        """
        # Map corners (i.e. CrossingStrands) to faces.
        face_of = {corner: n for n, face in enumerate(self.faces())
                   for corner in face}

        # Create the edges, labeled with crossing and sign.
        edges = []
        for c in self.crossings:
            edges.append((face_of[CrossingStrand(c, 0)],
                          face_of[CrossingStrand(c, 2)],
                          {'crossing': c, 'sign': 1}))
            edges.append((face_of[CrossingStrand(c, 1)],
                          face_of[CrossingStrand(c, 3)],
                          {'crossing': c, 'sign': -1}))

        # Build the graph.
        G = graph.Graph(edges, multiedges=True)
        components = G.connected_components()
        if len(components) > 2:
            raise ValueError('The link diagram is split.')
        return G.subgraph(components[1])

    @sage_method
    def goeritz_matrix(self, return_graph=False):
        """
        Call self.white_graph() and return the Goeritz matrix of the result.
        If the return_graph flag is set, also return the graph::

            sage: K=Link('4_1')
            sage: abs(K.goeritz_matrix().det())
            5
        """
        G = self.white_graph()
        V = G.vertices(sort=True)
        N = len(V)
        m = matrix(N, N)
        vertex = {v: n for n, v in enumerate(V)}
        for e in G.edges(sort=False):
            i, j = vertex[e[0]], vertex[e[1]]
            m[(i, j)] = m[(j, i)] = m[(i, j)] + e[2]['sign']
        for i in range(N):
            m[(i, i)] = -sum(m.column(i))
        m = m.delete_rows([0]).delete_columns([0])
        return (m, G) if return_graph else m

    @sage_method
    def signature(self, new_convention=True):
        """
        Returns the signature of the link, computed from the Goeritz matrix using
        the algorithm of Gordon and Litherland::

            sage: K = Link('4a1')
            sage: K.signature()
            0
            sage: L = Link('9^3_12')
            sage: Lbar = L.mirror()
            sage: L.signature() + Lbar.signature()
            0
            sage: M = Link(braid_closure=[1, 2, 1, 2, 1, 2, 1, 2])
            sage: M.signature()
            -6

        SnapPy 3.0 switched the sign convention for the signature so
        that "positive knots have negative signatures".  You can
        recover the previous default by::

          sage: L = Link('3a1')
          sage: L.signature()
          -2
          sage: L.signature(new_convention=False)
          2
        """
        m, G = self.goeritz_matrix(return_graph=True)
        correction = sum(e['sign'] for _, _, e in G.edges(sort=False)
                         if e['sign'] == e['crossing'].sign)
        ans = QuadraticForm(QQ, m).signature() + correction
        if new_convention:
            ans = -ans
        return ans

    @sage_method
    def _colorability_matrix(self):
        """Auxiliary function used by determinant."""
        edges = self._pieces()
        m = matrix(len(self.crossings), len(edges))
        for c in self.crossings:
            for i in range(4):
                for s in edges:
                    if (c, i) in s:
                        ind = edges.index(s)
                        if i % 2:
                            m[(self.crossings.index(c), ind)] += 1
                        else:
                            m[(self.crossings.index(c), ind)] -= 1
                        break
        return m

    @sage_method
    def determinant(self, method='goeritz'):
        """
        Returns the determinant of the link, a non-negative integer.

        Possible methods are 'wirt', using the Wirtinger presentation; 'goeritz',
        using the Goeritz matrix, and 'color', using the 'colorability matrix', or
        anything else, to compute the Alexander polynomial at -1.  Example::

            sage: K = Link( [(4,1,5,2),(6,4,7,3),(8,5,1,6),(2,8,3,7)] )  # Figure 8 knot
            sage: K.determinant()
            5
        """
        if method == 'color':
            M = self._colorability_matrix()
            size = len(self.crossings) - 1
            N = matrix(size, size)
            for i in range(size):
                for j in range(size):
                    N[(i, j)] = M[(i + 1, j + 1)]
            return abs(N.determinant())
        if method == 'goeritz':
            return abs(self.goeritz_matrix().determinant())
        return abs(self.alexander_polynomial(multivar=False, v=[-1], norm=False))

    @sage_method
    def morse_number(self, solver='GLPK'):
        """
        The *Morse number* of a planar link diagram D is

            m(D) = min { # of maxima of h on D }

        where h is a height function on R^2 which is generic on D; alternatively,
        this is the minimum number of cups/caps in a `MorseLink presentation
        <http://katlas.math.toronto.edu/wiki/MorseLink_Presentations>`_
        of the diagram D.  The Morse number is very closely related to the more
        traditional bridge number.   Examples::

            sage: K = Link('5_2')
            sage: K.morse_number()
            2
            sage: Link('6^3_2').morse_number()
            3
        """
        from . import morse
        return morse.morse_via_LP(self, solver)[0]

    @sage_method
    def morse_diagram(self):
        """
        Returns a MorseLinkDiagram of this link diagram, that is a choice
        of height function which realizes the Morse number::

            sage: L = Link('L8n2')
            sage: D = L.morse_diagram()
            sage: D.morse_number == L.morse_number()
            True
            sage: D.is_bridge()
            True
            sage: B = D.bridge()
            sage: len(B.bohua_code())
            64
        """
        from . import morse
        return morse.MorseLinkDiagram(self)

    @sage_method
    def jones_polynomial(self, variable=None, new_convention=True):
        """
        Returns the Jones polynomial of the link, following the
        conventions of https://arxiv.org/abs/math/0201043

        In particular, it obeys the oriented skein relation::

          q^2 V(L-) - q^-2 V(L+) = (q - q^-1) V(L0)

        and V(n-component unlink) = (q + q^-1)^(n-1).

        WARNING: The default conventions changed in SnapPy 3.0.  You
        can recover the old conventions as illustrated below::

          sage: L = Link('8_5')
          sage: J = L.jones_polynomial(); J
          1 - q^2 + 3*q^4 - 3*q^6 + 3*q^8 - 4*q^10 + 3*q^12 - 2*q^14 + q^16
          sage: Jold = L.jones_polynomial(new_convention=False); Jold
          1 - q + 3*q^2 - 3*q^3 + 3*q^4 - 4*q^5 + 3*q^6 - 2*q^7 + q^8

        Here are the values one unlinks with 4 and 5 components::

          sage: U4 = Link(braid_closure=[1, -1, 2, -2, 3, -3])
          sage: U5 = Link(braid_closure=[1, -1, 2, -2, 3, -3, 4, -4])
          sage: U4.jones_polynomial().factor()
          (q^-3) * (1 + q^2)^3
          sage: U5.jones_polynomial().factor()
          (q^-4) * (1 + q^2)^4
          sage: U4.jones_polynomial(new_convention=False).factor()
          (-q^-2) * (1 + q)^3
          sage: U5.jones_polynomial(new_convention=False).factor()
          (q^-2) * (1 + q)^4

        """
        from . import jones, jones_old

        if new_convention:
            J = jones.jones_polynomial(self, normalized=True)
        else:
            if len(self.link_components) != 1:
                J = jones_old.Jones_poly(self, new_convention=False)
            else:
                J = jones.jones_polynomial(self, normalized=True)
                R = J.parent()
                q = R.gen()
                terms = [J[e] * q**(e // 2) for e in J.exponents()]
                J = sum(terms, R(0))

        if variable is not None:
            J = J(variable)
        return J

    def seifert_matrix(self):
        """
        Returns the Seifert matrix of the link::

            sage: L = Link('K10n11')
            sage: A = L.seifert_matrix()
            sage: alex = L.alexander_polynomial()
            sage: t = alex.parent().gen()
            sage: B = t*A - A.transpose()
            sage: t**4 * alex == -B.det()
            True

        Uses the algorithm described in

        J. Collins, "An algorithm for computing the Seifert matrix of a link
        from a braid representation." (2007).

        after first making the link isotopic to a braid closure.
        """
        from . import seifert
        ans = seifert.seifert_matrix(self)
        if _within_sage:
            ans = matrix(ans)
        return ans

    def braid_word(self, as_sage_braid=False):
        """
        Return a list of integers which defines a braid word whose closure is the
        given link.  The natural numbers 1, 2, 3, etc are the generators and the
        negatives are the inverses.

        >>> L = Link('K6a2')
        >>> word = L.braid_word()
        >>> B = Link(braid_closure=word)
        >>> B.exterior().identify()    # doctest: +SNAPPY
        [m289(0,0), 6_2(0,0), K5_19(0,0), K6a2(0,0)]

        Within Sage, you can get the answer as an element of the
        appropriate BraidGroup and also check our earlier work::

            sage: Link('K6a2').braid_word(as_sage_braid=True)
            (s0*s1^-1)^2*s0^2
            sage: L.signature(), B.signature()
            (-2, -2)

        Implementation follows P. Vogel, "Representation of links by
        braids, a new algorithm".
        """
        from . import seifert
        word = seifert.braid_word(self)
        if as_sage_braid:
            if not _within_sage:
                raise ValueError('Requested Sage braid outside of Sage.')
            n = max(abs(a) for a in word) + 1
            word = BraidGroup(n)(word)
        return word

    @sage_method
    def sage_link(self):
        """
        Convert to a SageMath Knot or Link::

           sage: L = Link('K10n11')   # Spherogram link
           sage: K = L.sage_link(); K
           Knot represented by 10 crossings
           sage: L.alexander_polynomial()/K.alexander_polynomial()  # Agree up to units
           -t^3
           sage: L.signature(), K.signature()
           (-4, -4)

        Can also go the other way::

           sage: L = Link('K11n11')
           sage: M = Link(L.sage_link())
           sage: L.signature(), M.signature()
           (-2, -2)

        Can also take a braid group perspective.

            sage: B = BraidGroup(4)
            sage: a, b, c = B.gens()
            sage: Link(braid_closure=(a**-3) * (b**4) * (c**2) * a * b * c )
            <Link: 2 comp; 12 cross>
            sage: L = Link(a * b * c); L
            <Link: 1 comp; 3 cross>
            sage: S = L.sage_link(); S
            Knot represented by 3 crossings
            sage: Link(S)
            <Link: 1 comp; 3 cross>
        """
        if SageKnot is None:
            raise ValueError('Your SageMath does not seem to have a native link type')
        sage_type = SageKnot if len(self.link_components) == 1 else SageLink
        # Sage's PD_code lists strands *clockwise* not our *anticlockwise*.
        code = [[x[0], x[3], x[2], x[1]] for x in self.PD_code(min_strand_index=1)]
        return sage_type(code)

    @sage_method
    def _sage_(self):
        """
        A quick test:

            sage: L = Link('K13n100')
            sage: L._sage_()
            Knot represented by 13 crossings
        """
        return self.sage_link()


class ClosedBraid(Link):
    """
    This is a convenience class for constructing closed braids.

    The constructor accepts either a single argument, which should be a list of
    integers to be passed to the Link constructor as the braid_closure
    parameter, or one or more integer arguments which will be packaged as a list
    and used as the braid_closure parameter.

    >>> B = ClosedBraid(1,-2,3)
    >>> B
    ClosedBraid(1, -2, 3)
    >>> B = ClosedBraid([1,-2,3]*3)
    >>> B
    ClosedBraid(1, -2, 3, 1, -2, 3, 1, -2, 3)
    """
    def __init__(self, *args, **kwargs):
        if args and 'braid_closure' not in kwargs:
            if len(args) == 1:
                self.braid_word = kwargs['braid_closure'] = tuple(args[0])
                args = ()
            elif isinstance(args[0], int):
                self.braid_word = kwargs['braid_closure'] = args
                args = ()
        Link.__init__(self, *args, **kwargs)

    def __repr__(self):
        return 'ClosedBraid%s' % str(self.braid_word)
