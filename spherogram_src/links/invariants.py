"""
When used within Sage, the Link class gains many methods to compute
standard invariants.  Much of this code was contributed by Robert Lipshitz
and Jennet Dickinson.  
"""

from . import links
from .links import Crossing, Strand

import sage.all
import sage.graphs.graph
from sage.matrix.constructor import matrix
from sage.symbolic.ring import SR
from sage.groups.free_group import FreeGroup
import sage.graphs.graph as graph
from sage.symbolic.ring import var
from sage.groups.braid import Braid, BraidGroup
from sage.all import ZZ

def normalize_alex_poly(p,t):
    # Normalize the sign of the leading coefficient 
    l = p
    for v in t:
        l = l.leading_coefficient(v)
    if l < 0:
        p = -p            
    for i in range(0,len(t)):
        exps = [x[1] for x in p.coeffs(t[i])]
        a = max(exps)
        b = min(exps)
        c = -1*(a+b)/ZZ(2)
        p = p*(t[i]**c)
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


class Link(links.Link):
    """
    Links are made from Crossings.  The general model is that of a PD
    diagram as described in

    http://katlas.org/wiki/Planar_Diagrams

    See the file "doc.pdf" for the conventions, which can be accessed
    via "spherogram.pdf_docs()", and the file "test.py" for some
    examples of creating links.

    Here are two ways of creating the figure-8 knot, first via a PD code

    >>> K1 = Link([[8,3,1,4],[2,6,3,5],[6,2,7,1],[4,7,5,8]])

    and by directly gluing up Crossings:

    >>> a, b, c, d = [Crossing(x) for x in 'abcd']
    >>> a[0], a[1], a[2], a[3] = c[1], d[0], b[1], b[0]
    >>> b[2], b[3] = d[3], c[2]
    >>> c[3], c[0] = d[2], d[1]
    >>> K2 = Link([a,b,c,d])

    Some families of named links are available, such a torus knots

    >>> Link('T(4, 2)')
    <Link: 2 comp; 6 cross>

    and if you have SnapPy installed also the Rolfsen and Hoste-Thistlethwaite
    tables:

    >>> Link('8_20')
    <Link 8_20: 1 comp; 8 cross>
    >>> Link('K12a123')
    <Link K12a123: 1 comp; 12 cross>
    >>> Link('L12n123')
    <Link L12n123: 2 comp; 12 cross>

    You can also create links from a Sage braid:

    >>> B = BraidGroup(4)
    >>> a, b, c = B.gens()
    >>> Link( (a**-3) * (b**4) * (c**2) * a * b * c )
    <Link: 2 comp; 12 cross>
    """
    def __init__(self, crossings, check_planarity=True, build=True):
        if isinstance(crossings, Braid):
            crossings = braidword_to_crossings(crossings)
        if isinstance(crossings, str) and crossings.startswith('T('):
            p, q = map(int, crossings[2:-1].split(','))
            B = BraidGroup(p)
            b = B([i+1 for i in range(0,p-1)]*q)
            crossings = braidword_to_crossings(b)
        
        links.Link.__init__(self, crossings, check_planarity, build)

    
    def linking_matrix(self):
        """
        Calcluates the linking number for each pair of link components.
        
        Returns a linking matrix, in which the (i,j)th component is the linking                                    
        number of the ith and jth link components.
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

    def knot_group(self):
        """
        Computes the knot group using the Wirtinger presentation. 
        Returns a finitely presented group.
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

    def alexander_matrix(self, mv=True):
        """
        Returns the alexander matrix of self.

        >>> L = Link('3_1')
        >>> L.alexander_matrix()
        ([      -1 -1/t + 1      1/t]
        [     1/t       -1 -1/t + 1]
        [-1/t + 1      1/t       -1], [t, t, t])
        >>> L = Link([(4,1,3,2),(1,4,2,3)])
        >>> L.alexander_matrix()  # doctest: +SKIP
        ([ t1 - 1 -t2 + 1]
        [-t1 + 1  t2 - 1], [t2, t1])
        """
        comp = len(self.link_components)
        if comp < 2:
            mv = False

        G = self.knot_group()
        B = G.alexander_matrix()
        g = list(var('g%d' % (i+1)) for i in range(len(G.gens())))

        if(mv):
            t = list(var('t%d' % (i+1)) for i in range(0,comp))
        else:
            t = [var('t')]*comp

        import sage.symbolic.relation as rel

        eq = [y(g) == 1 for y in G.relations()]
        solns = rel.solve(eq,g, solution_dict=True)
        r = list(set([solns[0][h] for h in g]))
        dict1 = {r[i]:t[i] for i in range(len(t))}

        for i in range(len(g)):
            g[i] = g[i].subs(solns[0]).subs(dict1)

        m = B.nrows()
        n = B.ncols()
        C = matrix(SR,m,n)

        for i in range(n):
            for j in range(m):
                for k in B[j,i].terms():
                    x = k.leading_item()
                    C[j,i] = C[j,i] + x[1]*x[0](g)

        return (C,g)
    
    def alexander_poly(self, multivar=True, v='no', method='wirt', norm = True):
        """
        Calculates the alexander polynomial of self. For links with one component,
        can evaluate the alexander polynomial at v.

        >>> K = Link('4_1')
        >>> K.alexander_poly()
        t + 1/t - 3
        >>> K.alexander_poly(v=[4])
        5/4

        >>> K = Link('L7n1')
        >>> K.alexander_poly(norm=False)
        (t1*t2^3 + 1)/(t1*t2^4)
        """

        # sign normalization still missing, but when "norm=True" the
        # leading coefficient with respect to the first variable is made
        # positive. 
        if method == 'snappy':
            try:
                import snappy
                return snappy.snap.alexander_polynomial(self.exterior())
            except ImportError:
                raise RunTimeError('this method for alexander_poly '+no_snappy_msg)
        else:
            comp = len(self.link_components)
            if comp < 2:
                multivar = False

            if(multivar):
                t = list(var('t%d' % (i+1)) for i in range(0,comp))
            else:
                t = [var('t')]

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

            if multivar:
                t_i = M[1][-1]
                p = (p.factor())/(t_i-1)

            if(norm):
                p = normalize_alex_poly(p.expand(),t)

            if v != 'no':
                dict1 = {t[i]:v[i] for i in range(len(t))}
                return p.subs(dict1)
                
            if multivar: # it's easier to view this way
                return p.factor()
            else:
                return p.expand()

    def _edge_sign(K, edge):
        "Returns the sign (+/- 1) associated to given edge in the black graph."
        crossing = edge[2]
        if set(((crossing,0),(crossing,1))).issubset(set(edge[0])) or set(((crossing,0),(crossing,1))).issubset(set(edge[1])):
            return +1
        return -1
        
    def _find_crossing(self, e):
        """
        Auxiliary function used by signature to find which
        crossing corresponds to an edge in the black graph.
        """
        a=set(e[0])
        b=set(e[1])
        s=a.union(b)
        for x in range(len(self.crossings)):
            crossings=[self.crossings[x][0],self.crossings[x][1],self.crossings[x][2],self.crossings[x][3]]
            total=set(crossings)
            if total.issubset(s):
                return self.crossings[x]

    def black_graph(self):
        """
        Returns the black graph of K. If the black graph is disconnected (which
        can only happen for a split link diagram), returns one connected component.

        The edges are labeled by the crossings they correspond to.

        >>> K=Link('5_1')                                                                                
        >>> K.black_graph()
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

        G=graph.Graph(coords)
        component=G.connected_components()[1]
        G=G.subgraph(component)
        return G


            
    def goeritz_matrix(self):
        """
        Finds the black graph of a knot, and from that, returns the Goeritz matrix of that knot.
       
        >>> K=Link('4_1')
        >>> K.goeritz_matrix().det()
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

    def signature(self):
        """
        Returns the signature of the link.       
        >>> K = Link('4a1')                                                                               
        >>> K.signature()                                                                              
        0
        >>> L = Link('9^3_12')
        >>> Lbar = L.mirror()
        >>> L.signature() + Lbar.signature()
        0
        """
        answer=0
        G = self.black_graph()
        for e in G.edges():
            v = e[2]
            if self._edge_sign(e) == v.sign:
                answer = answer + v.sign
        m=self.goeritz_matrix()
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

    def colorability_matrix(self):
        """Auxiliary function used by determinant. Returns 'colorability matrix'."""
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

    def determinant(self, method='goeritz'):
        """Returns the determinant of the knot K, a non-negative integer.                

        Possible methods are 'wirt', using the Wirtinger presentation; 'goeritz',
        using the Goeritz matrix, and 'color', using the 'colorability matrix', or
        anything else, to compute the Alexander polynomial at -1.
        
        >>> K = Link( [(4,1,5,2),(6,4,7,3),(8,5,1,6),(2,8,3,7)] )  # Figure 8 knot
        >>> K.determinant()                                                        
        5
        """
        if method=='color':
            M = self.colorability_matrix()
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