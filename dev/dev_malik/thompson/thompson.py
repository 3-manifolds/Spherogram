import matplotlib.pyplot as plt
import itertools
from spherogram.graphs import FatGraph
from random import randint

class DyadicRational:
    def __init__(self, numerator, denom_exp):
        self.numerator = numerator
        self.denom_exp = denom_exp
        self.reduce()

    def reduce(self):
        if self.numerator == 0:
            return
        while self.numerator % 2 == 0:
            self.numerator = self.numerator // 2
            self.denom_exp -= 1

    def children(self):
        l = DyadicRational(2*self.numerator-1, self.denom_exp+1)
        r = DyadicRational(2*self.numerator+1, self.denom_exp+1)
        return (l,r)

    def parent(self):
        if self.numerator % 4 == 1:
            return DyadicRational(self.numerator+1,self.denom_exp)
        else:
            return DyadicRational(self.numerator-1,self.denom_exp)

    def __repr__(self):
        return '%s/2^%s' %(self.numerator,self.denom_exp)

    def __eq__(self,otherdyadic):
        if self.numerator == 0:
            return otherdyadic.numerator == 0
        else:
            numeq = (self.numerator == otherdyadic.numerator)
            deneq = (self.denom_exp == otherdyadic.denom_exp)
            return (numeq and deneq)

    def __lt__(self,otherdyadic):
        m = max(self.denom_exp,otherdyadic.denom_exp)
        diffself, diffother = (m-self.denom_exp),(m-otherdyadic.denom_exp)
        return (self.numerator * (2**diffself)) < (otherdyadic.numerator * (2**diffother))

    def float(self):
        return float(self.numerator) / (2 ** self.denom_exp)

    def __mul__(self,other):
        return DyadicRational(self.numerator * other.numerator, self.denom_exp + other.denom_exp)

    def __add__(self,other):
        m = max(self.denom_exp,other.denom_exp)
        diffself, diffother = (m-self.denom_exp),(m-other.denom_exp)
        new_num = (self.numerator * (2**diffself)) + (other.numerator * (2**diffother))
        return DyadicRational(new_num,m)

    def __sub__(self,other):
        m = max(self.denom_exp,other.denom_exp)
        diffself, diffother = (m-self.denom_exp),(m-other.denom_exp)
        new_num = (self.numerator * (2**diffself)) - (other.numerator * (2**diffother))
        return DyadicRational(new_num,m)

    def __div__(self,other):
        if self.numerator % other.numerator == 0:
            return DyadicRational(self.numerator/other.numerator,self.denom_exp-other.denom_exp)
        else:
            raise Exception('Quotient not a dyadic rational')

class TreeSequence:
    def __init__(self, s):
        self.dyadic_seq = s
        self._make_tree_sequence()

    def _make_tree_sequence(self):
        #prevent duplicates
        #self.dyadic_seq = list(set(self.dyadic_seq))

        #add necessary values
        zero = DyadicRational(0,1)
        one = DyadicRational(1,0)
        half = DyadicRational(1,1)
        if zero not in self.dyadic_seq:
            self.dyadic_seq.insert(0,zero)
        if one not in self.dyadic_seq:
            self.dyadic_seq.append(one)
        if half not in self.dyadic_seq:
            self.dyadic_seq.append(half)

        #add all parents
        while not self._has_all_parents():
            for d in self.dyadic_seq:
                if d not in [zero,one,half]:
                    if d.parent() not in self.dyadic_seq:
                        self.dyadic_seq.append(d.parent())
        self.dyadic_seq.sort()


    def _has_all_parents(self):
        for d in self.dyadic_seq:
            if d not in [DyadicRational(0,1), DyadicRational(1,1), DyadicRational(1,0)]:
                if not self.find_dyadic(d.parent()):
                    return False
        return True

    def _is_ordered(self):
        for i in range(self.length()-1):
            if not self.dyadic_seq[i] < self.dyadic_seq[i+1]:
                return False
        return True

    def length(self):
        return len(self.dyadic_seq)

    def find_dyadic(self,d):
        for dyadic in self.dyadic_seq:
            if d == dyadic:
                return True

    def delete_dyadic(self,position):
        self.dyadic_seq.pop(position)

    def __repr__(self):
        reps = [repr(d) for d in self.dyadic_seq]
        return ' < '.join(reps)

    def merge(self,other):
        dyadic_seq = []
        i = 0
        j = 0
        while i < self.length() or j < other.length():
            if self.dyadic_seq[i] == other.dyadic_seq[j]:
                dyadic_seq.append(self.dyadic_seq[i])
                i += 1
                j += 1
            elif self.dyadic_seq[i] < other.dyadic_seq[j]:
                dyadic_seq.append(self.dyadic_seq[i])
                i += 1
            else:
                dyadic_seq.append(other.dyadic_seq[j])
                j += 1
        return TreeSequence(dyadic_seq)

    def __contains__(self, dyadic_rational):
        for dyadic in self.dyadic_seq:
            if dyadic_rational == dyadic:
                return True

    def tree_vertices(self, start_node=DyadicRational(1, 1)):
        if start_node not in self:
            raise Exception('Starting node not in TreeSequence')
        left_child, right_child = start_node.children()
        left_in, right_in = left_child in self, right_child in self
        if not left_in and not right_in:
            return [left_child, right_child]
        if left_in and not right_in:
            L = self.tree_vertices(start_node=left_child)
            L.append(right_child)
            return L
        if not left_in and right_in:
            L = [left_child]
            L.extend(self.tree_vertices(start_node=right_child))
            return L
        if left_in and right_in:
            L = self.tree_vertices(start_node=left_child)
            L.extend(self.tree_vertices(start_node=right_child))
            return L

    def exponents(self, start_node=DyadicRational(1, 1)):
        vertices = self.tree_vertices()
        exps = []
        for v in vertices:
            exp = 0
            while v.numerator % 4 == 1 and v.denom_exp > 1:
                v = v.parent()
                exp += 1
            if ((2**v.denom_exp) - v.numerator) == 1:
                exp -= 1
            exps.append(max(exp,0))
        return exps

    #Take left child until you leave the treesequence, return first vertex outside TreeSequence
    #Analogous for right child, returns pair
    def geodesic_endpoints(self,start_vertex):
        if start_vertex not in self:
            return start_vertex, start_vertex
        left_child, right_child = start_vertex.children()
        left_still_in_tree = left_child in self
        right_still_in_tree = right_child in self
        while left_still_in_tree:
            left_child = left_child.children()[0]
            left_still_in_tree = left_child in self
        while right_still_in_tree:
            right_child = right_child.children()[1]
            right_still_in_tree = right_child in self
        return left_child, right_child

class TreeMap:

    def __init__(self,TreeSeq1, TreeSeq2):
        self.break_points = TreeSeq1
        self.images = TreeSeq2
        if not TreeSeq1.length() == TreeSeq2.length():
            raise Exception("Sequences have different lengths")

    def size(self):
        return self.break_points.length()

    def graph(self):
        x = [d.float() for d in self.break_points.dyadic_seq]
        y = [d.float() for d in self.images.dyadic_seq]
        plt.plot(x,y)
        plt.show()

    def reduce(self):
        for i in reversed(range(self.size()-2)):
            if self.slope(i) == self.slope(i+1):
                if self.break_points.dyadic_seq[i+1] != DyadicRational(1,1):
                    if self.images.dyadic_seq[i+1] != DyadicRational(1,1):
                        self.break_points.delete_dyadic(i+1)
                        self.images.delete_dyadic(i+1)

    def slope(self, position):
        rise = self.images.dyadic_seq[position+1] - self.images.dyadic_seq[position]
        run = self.break_points.dyadic_seq[position+1] - self.break_points.dyadic_seq[position]
        return rise/run

    def __repr__(self):
        return repr(self.break_points)+'\n-->\n' + repr(self.images)

    def evaluate(self,dyadic_ratl):
        if dyadic_ratl == DyadicRational(0,1):
            return DyadicRational(0,1)
        if dyadic_ratl == DyadicRational(1,0):
            return DyadicRational(1,0)
        i = 0
        while self.break_points.dyadic_seq[i] < dyadic_ratl:
            i += 1
        #i is now right endpoint of interval in which dyadic_ratl falls
        slope = self.slope(i-1)

        left_endpoint,left_image = self.break_points.dyadic_seq[i-1],self.images.dyadic_seq[i-1]

        diff = dyadic_ratl-left_endpoint

        return left_image+(diff*slope)

    def inverse(self):
        return TreeMap(self.images,self.break_points)

    def exponents(self):
        return self.break_points.exponents(), self.images.exponents()

    def is_identity(self):
        br_exps = self.break_points.exponents()
        im_exps = self.images.exponents()
        for i in range(len(br_exps)):
            if br_exps[i] != im_exps[i]:
                return False
        return True

    #switched self and other to make ordering better
    def __mul__(self,other):
        inv = other.inverse()
        inv_breaks = TreeSequence([inv.evaluate(d) for d in self.break_points.dyadic_seq])
        break_points = other.break_points.merge(inv_breaks)
        images = TreeSequence([self.evaluate(other.evaluate(d)) for d in break_points.dyadic_seq])
        return TreeMap(break_points,images)

    def __pow__(self,exp):
        result = TreeMap(TreeSequence([]),TreeSequence([])) #identity map
        if exp > 0:
            for i in range(exp):
                result *= self
        if exp < 0:
            inv = self.inverse()
            for i in range(-exp):
                result *= inv
        return result

    def __eq__(self,other):
        return (self*(other.inverse())).is_identity()

    def planar_graph(self):
        # Make two graphs, merge at end
        G1 = FatGraph()
        G2 = FatGraph()
        vertices_copy = None #save vertices labels for later
        inverted_vertices_copy = None
        for tree_seq, inverted in [(self.break_points,False), (self.images,True)]:
            G = G1
            if inverted:
                G = G2
            dyadics = tree_seq.dyadic_seq[1:-1] #trim off 0 and 1
            geodesic_endpts = [[tree_seq.geodesic_endpoints(d.children()[0]),d] for d in dyadics]

            vertices = tree_seq.tree_vertices()
            if not inverted:
                vertices_copy = [str(v) for v in vertices]
            else:
                inverted_vertices_copy = [str(v) for v in vertices]
            edges = []
            for triple in geodesic_endpts:
                l, r = triple[0]
                d = triple[1]
                r_ind = vertices.index(r)
                r = vertices[vertices.index(r)+1]
                edges.append( [(l,r),d] )

            edges_copy = edges[:]

            #grouping by left vertex
            grouped_edges = []
            for triple in edges:
                l = triple[0][0]
                group = []
                for other_triple in edges_copy:
                    if other_triple[0][0] == l:
                        group.append(other_triple)
                for other_triple in group: #remove grouped elements from bank of elements
                    edges_copy.remove(other_triple)
                #sort by biggest denom first
                group.sort(key=lambda x: x[1].denom_exp, reverse=True)
                if group:
                    grouped_edges.append(group)

            #switch out denominators for indices
            #prep to add edge

            for group in grouped_edges:
                for n, triple in enumerate(group):
                    l, r = triple[0]
                    G.add_edge( (str(l),n), (str(r),-1), 0)

            #switch -1's for farthest point counterclockwise
            for edge in G.edges:
                for v in edge.incident_to():
                    edge.set_slot(v,edge.slot(v) % len(G.incidence_dict[v]))

            #if bottom side (images) switch all positions to negative
            if inverted:
                for edge in G.edges:
                    edge.slots[0] = -edge.slots[0]-1
                    edge.slots[1] = -edge.slots[1]-1

        #Merge graphs, ordering vertices left to right
        graph_nodes = ['v'+str(i) for i in range(len(vertices_copy))]
        G = FatGraph()
        edges = G1.edges
        inv_edges = G2.edges
        while edges:
            edge = edges.pop()
            vert1, vert2 = edge.incident_to()
            v1G, v2G = graph_nodes[vertices_copy.index(vert1)], graph_nodes[vertices_copy.index(vert2)]
            G.add_edge( (v1G,edge.slot(vert1)), (v2G,edge.slot(vert2)), 0)
        while inv_edges:
            edge = inv_edges.pop()
            vert1, vert2 = edge.incident_to()
            v1G, v2G = graph_nodes[inverted_vertices_copy.index(vert1)], graph_nodes[inverted_vertices_copy.index(vert2)]
            G.add_edge( (v1G,edge.slot(vert1)), (v2G,edge.slot(vert2)), 1)
        return G


class ABWord:

    def __init__(self, s, a, b):
        for letter in s:
            if letter not in 'abAB':
                raise Exception('Not a valid word')
        self.word_string = s
        self.a = a
        self.b = b
        self.A = a**(-1)
        self.B = b**(-1)
        self.identity = a*(a**(-1))

    def evaluate(self):
        product = self.identity
        for letter in self.word_string:
            product = product * (eval('self.' + letter))
        return product

    def __repr__(self):
        return self.word_string

    def is_reduced(self):
        for i in range(len(self.word_string) - 1):
            u, v = self.word_string[i: i + 2]
            if u == 'a' and v == 'A':
                return False
            if u == 'A' and v == 'a':
                return False
            if u == 'b' and v == 'B':
                return False
            if u == 'B' and v == 'b':
                return False
        return True

    def reduce(self):
        while not self.is_reduced():
            self.word_string = self.word_string.replace('aA','')
            self.word_string = self.word_string.replace('Aa','')
            self.word_string = self.word_string.replace('bB','')
            self.word_string = self.word_string.replace('Bb','')

    def __mul__(self,other):
        return ABWord(''.join([self.word_string, other.word_string]), self.a, self.b)


def random_word(complexity):
    lets = 'abAB'
    word = lets[randint(0,3)]
    for i in range(complexity-1):
        last_let = word[i].swapcase()
        new_lets = lets.replace(last_let,'')
        word = ''.join([word, new_lets[randint(0, 2)]])
    return word

a = TreeMap( TreeSequence([DyadicRational(3,2)]) , TreeSequence([DyadicRational(1,2)]) )
b = TreeMap( TreeSequence([DyadicRational(7,3)]) , TreeSequence([DyadicRational(5,3)]) )
A = a.inverse()
B = b.inverse()

def random_sequence(length, bound):
    return [randint(0,bound) for i in range(length)]

num_gens = 100

xs = [a,b]
x = A*b*a
for i in range(num_gens-2):
    xs.append(x)
    x = A*x*a
x_invs = [i.inverse() for i in xs]

def exponents_to_tree_map(exps_domain, exps_range):
    result = a*A
    inv = a*A
    for n, exp in enumerate(exps_range):
        result = result * (xs[n]**exp)
    for n, exp in enumerate(exps_domain):
        inv = inv * (xs[n]**exp)
    return result*(inv.inverse())


def random_xs(length, complexity):
    forbidden_let = (None, None)
    word = ''
    result = a*A
    for i in range(length):
        inverse = randint(0,1)
        next_pos = randint(0,complexity)
        while (next_pos,inverse) == forbidden_let:
            inverse == randint(0,1)
            next_pos = randint(0,complexity)
        if inverse:
            result *= x_invs[next_pos]
            forbidden_let = (next_pos, 1-inverse)
            word += 'X'+str(next_pos)+' '
        else:
            result *= xs[next_pos]
            forbidden_let = (next_pos,1-inverse)
            word += 'x'+str(next_pos)+' '
    return result, word


#letters = 'abAB'
#letters_list = [letters for i in range(10)]
#i=0
#for elt in itertools.product(*letters_list):
    #print(i)
    #i += 1

#    s= ''.join(elt)
#    word = ABWord(s,a,b)
#    word_ev = word.evaluate()
#    word_ev.reduce()
#    if word.is_reduced() and word_ev.size() <= 3:
#        print(word)
#        print(word_ev)
