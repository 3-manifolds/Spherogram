from __future__ import print_function
from .graphs import ReducedGraph, Digraph, Poset
from collections import deque
import operator

class Alphabet(object):
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

    def __pow__(self, exponent):
        if exponent < 0:
            inverse = self.__invert__()
            return inverse**(-exponent)
        product = Word(exponent*list(self), alphabet=self.alphabet)
        product.cancel()
        return product

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

    def syllables(self):
        if len(self) == 0:
            return []
        ans, curr = [], None
        for x in self:
            g = abs(x)
            e = 1 if x > 0 else -1
            if g == curr:
                count += e
            else:
                if curr != None:
                    ans.append( (curr, count) )
                curr, count = g, e

        ans.append( (curr, count) )
        return ans

    def verbose_string(self, separator=' * '):
        ans = []
        alpha = self.alphabet
        for g, e in self.syllables():
            part = alpha[g] if e == 1 else alpha[g] + '^' + repr(e)
            ans.append(part)
        return separator.join(ans)
        

class CyclicWord(Word):
    
    def cancel(self):
        Word.cancel(self)
        while len(self) > 0 and self[0] == -self[-1]:
            self.popleft()
            self.pop()

    def __mul__(self):
        raise ValueError('Cyclic words cannot be multiplied.')

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

    def rewrite(self, ordering):
        seq = []
        for letter in self:
            if letter in ordering:
                seq.append(1 + ordering.index(letter))
            else:
                seq.append(-1 - ordering.index(-letter))
        return CyclicWord(seq)
                    
    def shuffle(self, perm_dict={}):
        """
        Permute generators according to the supplied dictionary.
        Keys must be positive integers.  Values may be negative.
        The set of keys must equal the set of values up to sign.
        """
        abs_image = set( map(operator.abs, perm_dict.values()) )
        if set(perm_dict.keys()) != abs_image:
            raise ValueError('Not a permutation!')
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
        for word in (self, ~self):
            for n in xrange(len(self)):
                complexity, Xordering = word.complexity(size, ordering, spin=n)
                if complexity < least:
                    least = complexity
                    minima = [ (CyclicWord(word.spun(n)), Xordering) ]
                elif complexity == least:
                    minima.append( (CyclicWord(word.spun(n)), Xordering)  )
        return least, minima
    
class Complexity(list):
    def __lt__(self, other):
        if len(self) == len(other):
            return list.__lt__(self, other)
        else:
            return len(self) > len(other)

    def __gt__(self, other):
        if len(self) == len(other):
            return list.__gt__(self, other)
        else:
            return len(other) > len(self)

    def __le__(self, other):
        return not self > other

    def __ge__(self, other):
        return not self < other
        
    def __cmp__(self, other):
        result = cmp( len(other), len(self) )
        if result == 0:
            result = cmp(self, other)
        return result

        
class WhiteheadMove(object):
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
            if x == self.letter or x == -self.letter:
                sub += self.alphabet[x]
            else:
                # Be careful with these minus signs!
                if -x not in self.cut_set and x != self.letter:
                    sub += self.alphabet[self.letter]
                sub += self.alphabet[x]
                if x not in self.cut_set and x != -self.letter:
                    sub += self.alphabet[-self.letter]
            subs.append(sub)
        return ', '.join(subs)
            
class Presentation(object):
    """
    A Presentation contains a list of CyclicWords as relators and a
    list of (integer) letters as generators.  The generators are
    implicitly defined by their appearance in some relator.  Other
    generators can be provided in an optional list.  Instantiate with
    a list of objects that can be used to instantiate a CyclicWord,
    i.e. strings or lists of non-zero integers.
    """ 
    def __init__(self, relators, generators=[], alphabet=ABC):
        self.alphabet = alphabet
        self.relators = []
        self.generators = set(generators)
        if isinstance(relators, str):
            raise ValueError('Please provide a sequence of relators.')
        for r in relators:
            W = CyclicWord(r, alphabet)
            if len(W) > 0:
                self.relators.append(W)
            if isinstance(r, str):
                self.generators.update([abs(alphabet(l)) for l in r])
            else:
                self.generators.update([abs(l) for l in r])

    def __repr__(self):
        generators = ', '.join([self.alphabet[x] for x in self.generators])
        return 'generators: [%s]\nrelators: %s'%(generators, self.relators)

    def __len__(self):
        return sum( [len(W) for W in self.relators] )

    def __eq__(self, other):
        return (self.relators == other.relators and
                self.generators == other.generators)
    
    def whitehead_graph(self):
        Wh = ReducedGraph()
        vertex_dict = {}
        for letter in self.generators:
            vertex_dict[letter] = letter
            vertex_dict[-letter] = -letter
            Wh.add_vertex(vertex_dict[letter])
            Wh.add_vertex(vertex_dict[-letter])
        for relator in self.relators:
            for n in range(-1, len(relator)-1):
                Wh.add_edge(vertex_dict[relator[n]], vertex_dict[-relator[n+1]])
        return Wh

    def find_reducers(self):
        whitehead = self.whitehead_graph()
        reducers = []
        levels = []
        for x in self.generators:
            cut = whitehead.one_min_cut(x, -x)
            valence = whitehead.multi_valence(x)
            length_change = cut['size'] - valence
            if length_change < 0:
                reducers.append( (length_change, x, cut['set']) )
            elif length_change == 0:
                levels.append( (x, cut) )
        reducers.sort(key=lambda x: x[0])
        return reducers, levels

    def whitehead_move(self, a, cut_set):
        """
        Return a presentation obtained by performing a Whitehead move
        (T-transformation) on this presentation.
        """
        new_relators = []
        for relator in self.relators:
            new_relator = []
            for x in relator:
                if x == a or x == -a:
                    new_relator.append(x)
                else:
                    # Be careful with these minus signs!
                    if -x not in cut_set:
                        new_relator.append(a)
                    new_relator.append(x)
                    if x not in cut_set:
                        new_relator.append(-a)
            W = CyclicWord(new_relator, self.alphabet)
            new_relators.append(W)
        return Presentation(new_relators, self.generators)


    def shorten(self):
        """
        Apply maximally reducing Whitehead moves until the minimal
        length is reached.  Return the resulting minimal presentation. 

        >>> P = Presentation(['AAAAABBAACCC', 'AAABBBCCCCC', 'AABDCCBD'])
        >>> P.whitehead_graph().is_planar()
        False
        >>> S = P.shorten()
        >>> print(S)
        generators: [A, B, C, D]
        relators: [AAAAABBAACCC, AAABBBCCCCC, AADCCD]
        >>> S.whitehead_graph().is_planar()
        True
        >>> P = Presentation(['xyyxyyxy', 'xyy'])
        >>> P.shorten()
        generators: [X, Y]
        relators: [x, y]
        """
        result = Presentation(self.relators, self.generators)
        starting_length = len(result)
        #print(result.relators)
        while True:
            reducers, levels = result.find_reducers()
            if not reducers:
                return result
            reduction, a, cut_set = reducers[0]
            #w = WhiteheadMove(a, cut_set, self.generators, self.alphabet)
            #print(w)
            result = result.whitehead_move(a, cut_set)
            #print(result.relators)
        return result

    def level_transformations(self):
        """
        Generator for non-trivial level transformations.

        >>> P = Presentation(['AABCaBacAcbabC'])
        >>> relators = []
        >>> for x, X in P.level_transformations():
        ...   P = Presentation(['AABCaBacAcbabC'])
        ...   P = P.whitehead_move(x, X)
        ...   relators += P.relators
        >>> sorted(relators)
        [AAABCBaaccbabC, AABaBcacAbaCbC, AABCaBaaccbbAC, AACaBacBAcabbC, ABCaaBcAAcbabC, ABCaaBacAcbbAC]
        """

#        For each generator x we find one minimal (x,x^-1)-cut.  We
#        then construct a digraph D from all of the saturated edges of
#        the associated maximal flow, using their flow directions.
#        (Both directions may occur, giving rise to cycles of length
#        2.)  For each edge which is not saturated by the flow, we add
#        a cycle of length 2.  We then form the DAG of strong
#        components of D, and its associated poset.  The transitively
#        closed subsets determine all cuts, which must be level
#        transformations.  This generator yields all level
#        transformations of the form (x,X) where neither X nor
#        its complement has size 1.

        reducers, levels = self.find_reducers()
        result = []
        if reducers:
             raise ValueError('Presentation is not minimal.')
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

    def level_orbit(self, verbose=False):
        """
        Generator for all presentations obtained from this one by
        length preserving Whitehead moves. Does a depth first search
        of the orbit.  If the verbose flag is True, yields a tuple:
        (parent, Whitehead move, result, canonical presentation).
        """
        S = self.signature()
        queue = deque([ (None, None, self, S) ])
        seen = set([S])
        while queue:
            parent, move, pres, sig = queue.popleft()
            for a, A in pres.level_transformations():
                P = Presentation(pres.relators, pres.generators)
                P = P.whitehead_move(a,A)
                signature = P.signature()
                if signature not in seen:
                    WM = WhiteheadMove(a,A, pres.generators, self.alphabet),
                    queue.append( (pres, WM, P, signature) )
                    seen.add(signature)
            if verbose:
                yield (parent, move, pres, Presentation(*sig))
            else:
                yield Presentation(*sig)

    def signature(self):
        """
        Return the relators of a canonized presentation as a tuple
        of tuples.  The result is hashable, but can be used to
        generate a canonical presentation equivalent to this one.
        """
        queue = deque()
        P = Presentation([], self.generators)
        queue.append(CanonizeNode(P, self.relators))
        done = False
        while True:
            if len(queue[0].remaining) > 0:
                left = queue.popleft()
                for child in left.children():
                    queue.append(child)
            else:
                break
        relators = queue[0].presentation.relators
        ordering = queue[0].ordering
        generators = tuple( range(1, len(self.generators) + 1 ) )
        return tuple([tuple(R.rewrite(ordering)) for R in relators]), generators

    def magma_string(self):
        gens = sorted([self.alphabet[g] for g in self.generators])
        ans = 'Group<' + ', '.join(gens) + ' | '
        ans += ', '.join([R.verbose_string() for R in self.relators])
        return ans + '>'
                          

class CanonizeNode(object):
    def __init__(self, presentation, remaining, ordering=[]):
        self.presentation = presentation
        self.generators = presentation.generators
        self.relators = presentation.relators
        self.remaining = remaining
        self.ordering = ordering

    def __repr__(self):
        return '%s\n%s'%(self.presentation, self.ordering)
    
    def children(self):
        childlist = []
        least = Complexity()
        length = len(self.generators)
        for relator in self.remaining:
            complexity, minima = relator.minima(length, self.ordering)
            if complexity > least:
                continue
            if complexity < least:
                least = complexity
                childlist = []
            for minimum in minima:
                word, ordering = minimum
                relators = self.relators + [word]
                remaining = list(self.remaining)
                remaining.remove(relator)
                P = Presentation(relators, generators=self.generators)
                childlist.append(CanonizeNode(P, remaining, ordering))
        return childlist
