"""
When used within Sage, the Link class gains many methods to compute
standard invariants.  Much of this code was contributed by Robert Lipshitz
and Jennet Dickinson.  
"""

from . import links_base
from .links_base import Crossing, Strand, CrossingStrand
from ..sage_helper import _within_sage, sage_method

if _within_sage:
    import sage.all
    import sage.graphs.graph
    from sage.matrix.constructor import matrix
    from sage.symbolic.ring import SR
    from sage.groups.free_group import FreeGroup
    import sage.graphs.graph as graph
    from sage.symbolic.ring import var
    from sage.groups.braid import Braid, BraidGroup
    from sage.all import ZZ, QQ
    from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
    from sage.graphs.graph import Graph 
    from sage.quadratic_forms.quadratic_form import QuadraticForm
else:
    pass 

def normalize_alex_poly(p,t):
    # Normalize the sign of the leading coefficient and make all exponents
    #positive


    if len(t) == 1:
        p = p*(t[0]**(-min(p.exponents())))
        if p.coefficients()[-1]<0:
            p = -p
        return p
   

    max_degree = max([sum(x) for x in p.exponents()])
    highest_monomial_exps = [x for x in p.exponents() if sum(x)==max_degree]
    leading_exponents = max(highest_monomial_exps)
    leading_monomial = reduce(lambda x,y: x*y,[t[i]**(leading_exponents[i]) for i in range(len(t))])
    l = p.monomial_coefficient(leading_monomial)
    
    if l < 0:
        p = -p

    for i in range(len(t)):
        exps = map(lambda x: x[i], p.exponents())
        min_exp = min(exps)
        p = p*(t[i]**(-min_exp))
    return p

def braidword_to_crossings(braidword):
    """
    Initialize a link from a sage braid word.
    """
    braidgens = list(braidword.parent().gens())
    strands = [Strand('s'+repr(i)) for i in range(len(braidgens)+1)]
 
    s = [(x,0) for x in strands] # start                                                               
    l = [(x,1) for x in strands] # loose ends     
                                                     
    braidsylls = []
    for b in braidword.syllables():
        if b[1]>0:
            braidsylls = braidsylls+[b[0]]*abs(b[1])
        if b[1]<0:
            braidsylls = braidsylls+[b[0]**-1]*abs(b[1])
    xings = [0]*len(braidsylls)

    for i, b in enumerate(braidsylls): 
        # for each syllable, there is a single crossing               
        label = "x"+repr(i)
        xings[i] = Crossing(label)
        for j, a in enumerate(braidgens): 
            # j tells us which two strands are crossing                

            if b == a: # if crossing is negative  
                xings[i][1] = l[j][0][l[j][1]]
                xings[i][0] = l[j+1][0][l[j+1][1]]
                l[j]   = (xings[i],2)
                l[j+1] = (xings[i],3)
                    
            if b**-1 == a: # if crossing is positive      
                xings[i][0] = l[j][0][l[j][1]]
                xings[i][3] = l[j+1][0][l[j+1][1]]
                l[j]   = (xings[i],1)
                l[j+1] = (xings[i],2)
                    
    for i in range(len(s)):
        s[i][0][s[i][1]] = l[i][0][l[i][1]]

    crossings = xings+strands
    return crossings


class Link(links_base.Link):
    __doc__ = links_base.Link.__doc__
    
    def __init__(self, crossings, check_planarity=True, build=True):
        if _within_sage:
            crossings = self._from_braid(crossings)
        
        links_base.Link.__init__(self, crossings, check_planarity, build)

    @sage_method
    def _from_braid(self, crossings):
        """
         sage: B = BraidGroup(4)
         sage: a, b, c = B.gens()
         sage: Link( (a**-3) * (b**4) * (c**2) * a * b * c )
         <Link: 2 comp; 12 cross>
         """
        if isinstance(crossings, Braid):
            crossings = braidword_to_crossings(crossings)
        return crossings

    @sage_method
    def linking_matrix(self):
        """
        Calcluates the linking number for each pair of link components.
        Returns a linking matrix, in which the (i,j)th component is the
        linking number of the ith and jth link components.
        """
        matrix = [ [0 for i in range(len(self.link_components)) ] for j in range(len(self.link_components)) ]
        for n1, comp1 in enumerate(self.link_components):
            for n2, comp2 in enumerate(self.link_components):
                tally = [ [0 for m in range(len(self.crossings)) ] for n in range(2) ]
                if not (comp1 == comp2):
                    for i, c in enumerate(self.crossings):
                        for x1 in comp1:
                            if x1[0] == c:
                                tally[0][i] += 1
                        for x2 in comp2:
                            if x2[0] == c:
                                tally[1][i] += 1
                for k, c in enumerate(self.crossings):
                    if (tally[0][k] == 1 and tally[1][k] == 1):
                        matrix[n1][n2] += 0.5 * (c.sign)
        for i1, m1 in enumerate(matrix):
            for i2, m2 in enumerate(m1):
                matrix[i1][i2] = int(matrix[i1][i2])
        return matrix

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
        g = list(F.gens())
        rels = []
        pieces = self.pieces()

        for z in self.crossings:
            for m, p in enumerate(pieces):
                for t, q in enumerate(p):
                    if q[0] == z:
                        if t == 0:
                            j = m
                        elif t == len(p)-1:
                            i = m
                        else:
                            k = m
            i+=1; j+=1; k+=1
            if z.sign > 0:
                r = F([-k,i,k,-j])
            if z.sign < 0:
                r = F([k,i,-k,-j])
            rels.append(r)

        G = F/rels
        return G

    @sage_method
    def alexander_matrix(self, mv=True):
        """
        Returns the Alexander matrix of the link::

            sage: L = Link('3_1')
            sage: L.alexander_matrix()
            ([       -1      t^-1 -t^-1 + 1]
            [-t^-1 + 1        -1      t^-1]
            [     t^-1 -t^-1 + 1        -1], [t, t, t])

            sage: L = Link([(4,1,3,2),(1,4,2,3)])
            sage: L.alexander_matrix()    #doctest: +SKIP
            ([ t1 - 1 -t2 + 1]
            [-t1 + 1  t2 - 1], [t2, t1])
        """
        comp = len(self.link_components)
        if comp < 2:
            mv = False

        G = self.knot_group()
        num_gens = len(G.gens())

        L_g = LaurentPolynomialRing(QQ,['g%d' % (i+1) for i in range(num_gens)])
        g = list(L_g.gens())

        if(mv):
            L_t = LaurentPolynomialRing(QQ,['t%d' % (i+1) for i in range(comp)])
            t = list(L_t.gens())
            #determine the component to which each variable corresponds
            g_component = [c.strand_components[2] for c in self.crossings]
            for i in range(len(g)):
                g[i] = t[g_component[i]]
            
        else:
            L_t = LaurentPolynomialRing(QQ,'t')
            t = L_t.gen()
            g = [t]*len(g)


        B = G.alexander_matrix(g)

        return (B,g)

    @sage_method
    def alexander_poly(self, multivar=True, v='no', method='wirt', norm = True, factored = False):
        """
        Calculates the Alexander polynomial of the link. For links with one component,
        can evaluate the alexander polynomial at v::

            sage: K = Link('4_1')
            sage: K.alexander_poly()
            1 - 3*t + t^2
            sage: K.alexander_poly(v=[4])
            5
            
            sage: K = Link('L7n1')
            sage: K.alexander_poly(norm=False)
            -t2 - t1*t2^-2
        """

        # sign normalization still missing, but when "norm=True" the
        # leading coefficient with respect to the first variable is made
        # positive. 
        if method == 'snappy':
            try:
                return self.exterior().alexander_polynomial()
            except ImportError:
                raise RuntimeError('this method for alexander_poly '+no_snappy_msg)
        else:
            comp = len(self.link_components)
            if comp < 2:
                multivar = False

            if(multivar):
                L = LaurentPolynomialRing(QQ,['t%d' % (i+1) for i in range(comp)])
                t = list(L.gens())
            else:
                L = LaurentPolynomialRing(QQ,'t')
                t = [L.gen()]

            M = self.alexander_matrix(mv=multivar)
            C = M[0]
            m = C.nrows()
            n = C.ncols()
            if n>m:
                k = m-1
            else:
                k = n-1
                
            subMatrix = C[0:k,0:k]
            p = subMatrix.determinant()
            if p == 0: return 0
            if multivar:
                t_i = M[1][-1]
                p = (p.factor())/(t_i-1)
                p = p.expand()

            if(norm):
                p = normalize_alex_poly(p,t)

            if v != 'no':
                dict1 = dict([(t[i],v[i]) for i in range(len(t))])
                return p.subs(dict1)
                
            if multivar and factored: # it's easier to view this way
                return p.factor()
            else:
                return p
        
    def _edge_sign(K, edge):
        "Returns the sign (+/- 1) associated to given edge in the black graph."
        crossing = edge[2]
        if set(((crossing,0),(crossing,1))).issubset(set(edge[0])) or set(((crossing,0),(crossing,1))).issubset(set(edge[1])):
            return +1
        return -1
        
    @sage_method
    def black_graph(self):
        """
        Returns the black graph of K. If the black graph is disconnected
        (which can only happen for a split link diagram), returns one
        connected component. The edges are labeled by the crossings
        they correspond to.  Example::

            sage: K=Link('5_1')                                                                                
            sage: K.black_graph()
            Subgraph of (): Multi-graph on 2 vertices
        """

        faces = []
        for x in self.faces():
            l = []
            for y in x:
                l.append((y[0],y[1]))
                l.append((y[0],(y[1]+1)%4))
            faces.append(l)

        coords=list()
        for i in range(len(faces)-1):
            for j in range (i+1, len(faces)):
                a=set(faces[i])
                b=set(faces[j])
                s=a.union(b)
                for x in range(len(self.crossings)):
                    crossings=[self.crossings[x][0],self.crossings[x][1],self.crossings[x][2],self.crossings[x][3]]
                    total=set(crossings)
                    if total.issubset(s):
                        coords.append((tuple(faces[i]),tuple(faces[j]),self.crossings[x])) #label by the crossing.

        G=graph.Graph(coords, multiedges=True)
        component=G.connected_components()[1]
        G=G.subgraph(component)
        return G

    @sage_method      
    def old_goeritz_matrix(self):
        """
        Finds the black graph of a knot, and from that, returns the Goeritz
        matrix of that knot::
        
            sage: K=Link('4_1')
            sage: K.goeritz_matrix().det()
            5
        """
        g=self.black_graph()
        l=g.vertices()
        m=matrix(len(g.vertices()),len(g.vertices()))
        for e in g.edges():
            i=l.index(e[0])
            j=l.index(e[1])
            m[(i,j)]=m[(i,j)]+self._edge_sign(e)
            m[(j,i)]=m[(j,i)]+self._edge_sign(e)
        for i in range(len(g.vertices())):
            m[(i,i)]=sum(m.column(i))*(-1)
        m=m.delete_rows([0])
        m=m.delete_columns([0])
        return m

    @sage_method
    def old_signature(self):
        """
        Returns the signature of the link.  Examples::
        
            sage: K = Link('4a1')                                                                               
            sage: K.signature()          
            0
            sage: L = Link('9^3_12')
            sage: Lbar = L.mirror()
            sage: L.signature() + Lbar.signature()
            0
        """
        answer=0
        G = self.black_graph()
        for e in G.edges():
            v = e[2]
            if self._edge_sign(e) == v.sign:
                answer = answer + v.sign
        m=self.old_goeritz_matrix()
        vals=m.eigenvalues()
        pos=0
        neg=0
        for v in vals:
            if v>0:
                pos+=1
            elif v<0:
                neg+=1
        sig=pos-neg+answer
        return sig

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
        is raised.

            sage: K=Link('5_1')           
            sage: K.white_graph()
            Subgraph of (): Multi-graph on 2 vertices
        """
        # Map corners (i.e. CrossingStrands) to faces.
        face_of = dict((corner, n) for n, face in enumerate(self.faces()) for corner in face) 

        # Create the edges, labeled with crossing and sign.
        edges = []
        for c in self.crossings:
            edges.append((face_of[CrossingStrand(c, 0)], face_of[CrossingStrand(c, 2)],
                          {'crossing':c, 'sign':1}))
            edges.append((face_of[CrossingStrand(c, 1)], face_of[CrossingStrand(c, 3)],
                          {'crossing':c, 'sign':-1}))
                
        # Build the graph.
        G = graph.Graph(edges, multiedges=True)
        components = G.connected_components()
        if len(components) > 2:
            raise ValueError, 'The link diagram is split.'
        return G.subgraph(components[1])

    @sage_method      
    def goeritz_matrix(self, return_graph=False):
        """
        Call self.white_graph() and return the Goeritz matrix of the result.
        If the return_graph flag is set, also return the graph.
        
            sage: K=Link('4_1')
            sage: abs(K.goeritz_matrix().det())
            5
        """
        G = self.white_graph()
        V = G.vertices()
        N = len(V)
        m = matrix(N,N)
        vertex = dict((v, n) for n, v in enumerate(V))
        for e in G.edges():
            i, j = vertex[e[0]], vertex[e[1]]
            m[(i,j)] = m[(j,i)] = m[(i,j)] + e[2]['sign']
        for i in xrange(N):
            m[(i,i)] = -sum(m.column(i))
        m = m.delete_rows([0]).delete_columns([0])
        if return_graph:
            return m, G
        else:
            return m

    @sage_method
    def signature(self):
        """
        Returns the signature of the link, computed from the Goeritz matrix using
        the algorithm of Gordon and Litherland.

            sage: K = Link('4a1')            
            sage: K.signature()          
            0
            sage: L = Link('9^3_12')
            sage: Lbar = L.mirror()
            sage: L.signature() + Lbar.signature()
            0
        """
        m, G = self.goeritz_matrix(return_graph=True)
        correction = sum([e[2]['sign'] for e in G.edges()
                          if e[2]['sign'] == e[2]['crossing'].sign])
        return QuadraticForm(QQ, m).signature() + correction

    @sage_method
    def _colorability_matrix(self):
        """Auxiliary function used by determinant."""
        edges = self.pieces()
        m=matrix(len(self.crossings), len(edges))
        for c in self.crossings:
            for i in range(4):
                for s in edges:
                    if (c,i) in s:
                        ind=edges.index(s)
                        if i==1 or i==3:
                            m[(self.crossings.index(c),ind)]+=1
                        if i==2 or i==0:
                            m[(self.crossings.index(c),ind)]-=1
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
        if method=='color':
            M = self._colorability_matrix()
            size = len(self.crossings)-1
            N = matrix(size, size)
            for i in range(size):
                for j in range(size):
                    N[(i,j)]=M[(i+1,j+1)]
            return abs(N.determinant())
        elif method=='goeritz':
            return abs(self.goeritz_matrix().determinant())
        else:
            return abs(self.alexander_poly(multivar=False, v=[-1], norm = False))

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
    def jones_poly(self, variable=None):
        """
        Returns the Jones polynomial of the link::

            sage: L = Link('8_5')
            sage: L.jones_poly()
            1 - q + 3*q^2 - 3*q^3 + 3*q^4 - 4*q^5 + 3*q^6 - 2*q^7 + q^8
        """
        from . import jones
        return jones.Jones_poly(self, variable)
