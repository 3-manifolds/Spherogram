from .graphs import ReducedGraph, Digraph, Poset
from collections import deque
import operator

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
        if int_list:
            return self.separator.join([self[x] for x in int_list])
        else:
            return self[0]
        
abc = Alphabet('1',
               'abcdefghijklmnopqrstuvwxyz',
               'ABCDEFGHIJKLMNOPQRSTUVWXYZ')

ABC = Alphabet('1',
               'ABCDEFGHIJKLMNOPQRSTUVWXYZ',
               'abcdefghijklmnopqrstuvwxyz')

class Word(deque):
    """
    A word in a free group is represented as a deque of non-zero
    integers.  Inverse corresponds to negation.  The optional alphabet
    controls how these lists are displayed to the user.
    """
    def __init__(self, word, alphabet=ABC):
        deque.__init__(self)
        if isinstance(word, str):
            for char in word:
                self.append(alphabet(char))
        else:
            self.extend(word)
        self.letters = set(map(abs,self))
        self.alphabet = alphabet
        self.cancel()

    def __mul__(self, other):
        product = Word(list(self) + list(other), alphabet=self.alphabet)
        product.cancel()
        return product

    def __invert__(self):
        inverse = [-x for x in reversed(self)]
        return Word(inverse, alphabet=self.alphabet)

    def __repr__(self):
        return self.alphabet.spell(self)

    def cancel(self):
        done = False
        while not done:
            done, size, n = True, len(self), 1
            self.rotate(-1)
            while len(self) > 0 and n < size:
                if self[0] == -self[-1]:
                    done = False
                    self.popleft()
                    self.pop()
                else:
                    self.rotate(-1)
                n += 1

class CyclicWord(Word):
    
    def cancel(self):
        Word.cancel(self)
        while len(self) > 0 and self[0] == -self[-1]:
            self.popleft()
            self.pop()

    def __mul__(self):
        raise ValueError, 'Cyclic words cannot be multiplied.'

    def __invert__(self):
        inverse = [-x for x in reversed(self)]
        return CyclicWord(inverse, alphabet=self.alphabet)
    
    def spun(self, start=0):
        """
        Generator for letters in cyclic order, starting at start.
        """
        N = len(self)
        for n in xrange(start, start+N):
            yield self[n%N]

    def invert(self):
        """
        Invert this cyclic word in place.
        """
        # This is builtin, starting from python 2.7
        for n in xrange(len(self)/2):
            self[n], self[-1-n] = self[-1-n], self[n]
        map(operator.neg, self)

    def shuffle(self, perm_dict={}):
        """
        Permute generators according to the supplied dictionary.
        Keys must be positive integers.  Values may be negative.
        The set of keys must equal the set of values up to sign.
        """
        abs_image = set( map(operator.abs, perm_dict.values()) )
        if set(perm_dict.keys()) != abs_image:
            raise ValueError, 'Not a permutation!'
        for n in xrange(len(self)):
            x = self[n]
            self[n] = perm_dict.get(x,x) if x > 0 else -perm_dict.get(-x,-x)
            
    def powers(self, start=0):
        """
        Return a list of pairs (letter, power) for the exponential
        representation of the spun word, beginning at start.
        """
        result = []
        last_letter = self[start]
        count = 0
        for letter in self.spun(start):
            if letter == last_letter:
                count += 1
            else:
                result.append( (last_letter, count))
                count = 1
                last_letter = letter
        result.append( (last_letter, count) )
        return result
    
    def complexity(self, size, ordering=[], spin=0):
        """
        Returns the complexity of the word relative to an extension of
        the ordering, and returns the extended ordering.  The
        lexicographical complexity is a list of integers, representing
        the ranking of each letter in the ordering.  The size is the
        total number of generators, which may be larger than the number
        of distinct generators that appear in the word.  If x appears
        in the ordering with rank r, then x^-1 has rank size + r.
        Unordered generators are added to the ordering as they are
        encountered in the word.
        """
        the_ordering = list(ordering)
        complexity = []
        for letter in self.spun(spin):
            if letter in the_ordering:
                complexity.append(the_ordering.index(letter))
            elif -letter in the_ordering:
                complexity.append(size + the_ordering.index(-letter))
            else:
                complexity.append(len(the_ordering))
                the_ordering.append(letter)
        return Complexity(complexity), the_ordering
                                  
    def minima(self, size, ordering=[]):
        """
        Return the minimal complexity of all rotations and inverted
        rotations, and a list of the words and orderings that realize
        the minimal complexity.
        """
        least = Complexity([])
        minima = []
        for inverted in (False, True):
            for n in xrange(len(self)):
                complexity, O = self.complexity(size, ordering, spin=n)
                if complexity < least:
                    least = complexity
                    minima = [ (CyclicWord(self.spun(n)), O) ]
                elif complexity == least:
                    minima.append( (CyclicWord(self.spun(n)), O)  )
            self.invert()
        return least, minima
    
class Complexity(list):
    def __lt__(self, other):
        if len(self) > len(other):
            return True
        else:
            return list.__lt__(self, other)
    def __gt__(self, other):
        if len(self) < len(other):
            return True
        else:
            return list.__gt__(self, other)
    def __le__(self, other):
        return not self < other

    def __gt__(self, other):
        return not self > other
        
    def __cmp__(self, other):
        result = cmp(len(other), len(self))
        if result == 0:
            result = cmp(self, other)
        return result

        
class WhiteheadMove:
    """
    Holds the data describing a Whitehead move.
    """
    def __init__(self, letter, cut_set, generators, alphabet):
        self.letter = letter
        self.cut_set = cut_set
        self.generators = generators
        self.alphabet = alphabet

    def __repr__(self):
        subs = []
        for x in self.generators:
            sub = '%s -> '%self.alphabet[x]
            # Be careful with these minus signs!
            if -x not in self.cut_set and x != self.letter:
                sub += self.alphabet[self.letter]
            sub += self.alphabet[x]
            if x not in self.cut_set and x != -self.letter:
                sub += self.alphabet[-self.letter]
            subs.append(sub)
        return ', '.join(subs)
            
class Presentation:
    """
    A Presentation contains a list of CyclicWords as relators and
    a list of (integer) letters as generators.  The generators are
    implicitly defined by their appearance in some relator.  (To add a
    free generator x, include the relation xX.)  Instantiate with a
    list of objects that can be used to instantiate a CyclicWord,
    i.e. strings or lists of non-zero integers.
    """ 
    def __init__(self, relator_list, alphabet=ABC):
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
            print WhiteheadMove(a, cut_set, self.generators, self.alphabet)
            self.whitehead_move(a, cut_set)
            print self.relators

    def level_transformations(self):
        """
        Generator for non-trivial level transformations.

        >>> P = Presentation(['AABCaBacAcbabC'])
        >>> for x, X in P.level_transformations():
        ...   P = Presentation(['AABCaBacAcbabC'])
        ...   P.whitehead_move(x, X)
        ...   print P, len(P)
        ...
        generators: [A, B, C]
        relators: [ABCaaBacAcbbAC] 14
        generators: [A, B, C]
        relators: [ABCaaBcAAcbabC] 14
        generators: [A, B, C]
        relators: [AABCaBaaccbbAC] 14
        generators: [A, B, C]
        relators: [AAABCBaaccbabC] 14
        generators: [A, B, C]
        relators: [AACaBacBAcabbC] 14
        generators: [A, B, C]
        relators: [AABaBcacAbaCbC] 14
        """

#        For each generator x we find one minimal (x,x^-1)-cut.  We
#        then construct a digraph D from all of the saturated edges of
#        the associated maximal flow, using their flow directions.
#        (Both directions may occur, giving rise to cycles of length
#        2.)  For each edge which is not saturated by the flow, we add
#        a cycle of length 2.  We then form the DAG of strong
#        components of D, and its associated poset.  The transitively
#        closed subsets determine all cuts, which must be level
#        transformations.  The generator yields all level
#        transformations of the form (x,X) where where neither X nor
#        its complement has size 1.

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
            for subset in P.closed_subsets():
                if 1 < len(subset) < len(P)-1:
                    yield generator, frozenset.union(*subset)

