from .graphs import ReducedGraph, Digraph, Poset

class Alphabet():
    """
    An Alphabet translates between integers and strings.
    Call as a function to go from string to integer; use getitem
    to go from integer to string.
    """
    def __init__(self, identity, pos, neg, separator=''):
        pos, neg = list(pos), list(neg)
        neg.reverse()
        self.strings = [identity] + pos + [identity] + neg
        self.separator = separator
        self.size = 1 + len(pos)

    def __getitem__(self, index):
        return self.strings[index]
        
    def __call__(self, string):
        n = self.strings.index(string)
        return n if n < self.size else n - 2*self.size

    def spell(self, int_list):
        """
        Convert a sequence of integers to a string.
        """
        return self.separator.join([self[x] for x in int_list])

abc = Alphabet('1',
               'abcdefghijklmnopqrstuvwxyz',
               'ABCDEFGHIJKLMNOPQRSTUVWXYZ')

class Word(list):
    """
    A word in a free group is represented as a list of non-zero
    integers.  Inverse corresponds to negation.  The optional alphabet
    controls how these lists are displayed to the user.
    """
    def __init__(self, word, alphabet=abc):
        list.__init__(self)
        if isinstance(word, str):
            for char in word:
                self.append(alphabet(char))
        else:
            self.extend(word)
        self.letters = set(map(abs,self))
        self.alphabet = alphabet
        self.cancel()

    def __mul__(self, other):
        product = Word(self + other, alphabet=self.alphabet)
        
        product.cancel()
        return product

    def __invert__(self):
        inverse = map(operator.neg, self)
        inverse.reverse()
        return Word(inverse, alphabet=self.alphabet)

    def __repr__(self):
        return self.alphabet.spell(self)

    def cancel(self):
        done = False
        while not done:
            n, done = 1, True
            while n < len(self):
                if self[n] == -self[n-1]:
                    done = False
                    self.pop(n-1)
                    self.pop(n-1)
                else:
                    n += 1

class CyclicWord(Word):
    
    def cancel(self):
        Word.cancel(self)
        while len(self) > 1 and self[0] == -self[-1]:
            self.pop()
            self.pop(0)
            
class Presentation:
    """
    A Presentation contains a list of CyclicWords as relators and
    a list of (integer) letters as generators.  The generators are
    implicitly defined by their appearance in some relator.  (To add a
    free generator x, include the relation xX.)  Instantiate with a
    list of objects that can be used to instantiate a CyclicWord,
    i.e. strings or lists of non-zero integers.
    """ 
    def __init__(self, relator_list, alphabet=abc):
        self.alphabet = alphabet
        self.relators = []
        self.generators = set()
        if not isinstance(relator_list, list):
            raise ValueError, 'Please provide a list of relators.'
        for r in relator_list:
            W = CyclicWord(r, alphabet)
            if len(W) > 0:
                self.relators.append(W)
            if isinstance(r,str):
                self.generators.update([abs(alphabet(l)) for l in r])
            else:
                self.generators.update(abs(r))
        self.build_reduced_whitehead_graph()

    def __repr__(self):
        generators = ', '.join([self.alphabet[x] for x in self.generators])
        return 'generators: [%s]\nrelators: %s'%(generators, self.relators)

    def __len__(self):
        return sum( [len(W) for W in self.relators] )
    
    def build_reduced_whitehead_graph(self):
        self.whitehead = Wh = ReducedGraph()
        for letter in self.generators:
            Wh.add_vertex(letter)
            Wh.add_vertex(-letter)
        for relator in self.relators:
            for n in range(-1, len(relator)-1):
                Wh.add_edge(relator[n], -relator[n+1])

    def find_reducers(self):
        reducers = []
        levels = []
        for x in self.generators:
            cut = self.whitehead.one_min_cut(x, -x)
            valence = self.whitehead.valence(x)
            length_change = cut['size'] - valence
            if length_change < 0:
                reducers.append( (length_change, x, cut['set']) )
            elif length_change == 0:
                levels.append( (x, cut) )
        reducers.sort(key=lambda x: x[0])
        return reducers, levels

    def whitehead_move(self, a, cut_set):
        """
        Perform a Whitehead move (T-transformation) on the presentation.
        """
        new_relators = []
        for relator in self.relators:
            new_relator = []
            for x in relator:
                # Be careful with these minus signs!
                if -x not in cut_set and x != a:
                    new_relator.append(a)
                new_relator.append(x)
                if x not in cut_set and x != -a:
                    new_relator.append(-a)
            W = CyclicWord(new_relator, self.alphabet)
            new_relators.append(W)
        self.relators = new_relators
        self.build_reduced_whitehead_graph()

    def shorten(self):
        """
        Apply maximally reducing Whitehead moves until the minimal
        length is reached.

        >>> P = Presentation(['AAAAABBAACCC', 'AAABBBCCCCC', 'AABDCCBD'])
        >>> P.whitehead.is_planar()
        False
        >>> P.shorten()
        [AAAAABBAACCC, AAABBBCCCCC, AABDCCBD]
        [AAAAABBAACCC, AAABBBCCCCC, AADCCD]
        >>> P.whitehead.is_planar()
        True
        >>> P = Presentation(['xyyxyyxy', 'xyy'])
        >>> P.shorten()
        [xyyxyyxy, xyy]
        [xyxyx, xy]
        [xyy, y]
        [xy, y]
        [x, y]
        """
        starting_length = len(self)
        print self.relators
        while True:
            reducers, levels = self.find_reducers()
            if not reducers:
                return starting_length - len(self)
            reduction, a, cut_set = reducers[0]
            self.whitehead_move(a, cut_set)
            print self.relators

    def find_level_transformations(self):
        """
        Return a list of non-trivial level transformations.

        >>> P = Presentation(['AABCaBacAcbabC'])
        >>> for x, X in P.find_level_transformations():
        ...   P = Presentation(['AABCaBacAcbabC'])
        ...   P.whitehead_move(x,X)
        ...   print P, len(P)
        ... 
        generators: [a, b, c]
        relators: [ABCaaBcAAcbabC] 14
        generators: [a, b, c]
        relators: [AABCBaaccbabCA] 14
        generators: [a, b, c]
        relators: [ABCaBaaccbbACA] 14
        generators: [a, b, c]
        relators: [AACaBacBAcabbC] 14
        generators: [a, b, c]
        relators: [AABaBcacAbaCbC] 14
        """

#        For each generator x we find a minimal (x,x^-1)-cut.  We then
#        construct a digraph D from all of the saturated edges of the
#        associated maximal flow, using their flow directions.  (Both
#        directions may occur, giving rise to cycles of length 2.)  For
#        each edge which is not saturated by the flow, we add a cycle of
#        length 2.  We then form the DAG of strong components of D, and
#        its associated poset.  The level transformations returned have
#        the form (x,X) where X is the union of all strong components
#        in the transitive closure of a single strong component, and
#        where neither X nor its complement has size 1.

        reducers, levels = self.find_reducers()
        result = []
        if reducers:
            raise ValueError, 'Presentation is not minimal.'
        for generator, cut in levels:
            edges = set()
            for weight, path in cut['paths']:
                for vertex, edge in path:
                    edges.add( (vertex, edge(vertex)) )
            for edge in cut['unsaturated']:
                x, y = edge
                edges.add( (x,y) )
                edges.add( (y,x) )
            D = Digraph(edges)
            P = Poset(D.component_DAG())
            for x in P:
                C = P.closure([x])
                if 1 < len(C) < len(P) - 1:
                    cut_set = set()
                    for x in C:
                        cut_set |= x
                    result.append( (generator, cut_set) )
        return result
