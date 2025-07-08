from collections import deque, Counter
import operator
import networkx as nx


def all_descendants(direct_acyclic_graph, nodes):
    ans = set(nodes)
    for n in nodes:
        ans.update(nx.descendants(direct_acyclic_graph, n))
    return ans


def nonempty_transitively_closed_subets(direct_acyclic_graph):
    for A in nx.antichains(direct_acyclic_graph):
        if len(A):
            yield all_descendants(direct_acyclic_graph, A)


class Alphabet:
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
        return n if n < self.size else n - 2 * self.size

    def spell(self, int_list):
        """
        Convert a sequence of integers to a string.
        """
        if int_list:
            return self.separator.join(self[x] for x in int_list)
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
        self.letters = set(map(abs, self))
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
        product = Word(exponent * list(self), alphabet=self.alphabet)
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
                if curr is not None:
                    ans.append((curr, count))
                curr, count = g, e

        ans.append((curr, count))
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

    def __mul__(self, other):
        raise ValueError('Cyclic words cannot be multiplied.')

    def __invert__(self):
        inverse = [-x for x in reversed(self)]
        return CyclicWord(inverse, alphabet=self.alphabet)

    def spun(self, start=0):
        """
        Generator for letters in cyclic order, starting at start.
        """
        N = len(self)
        for n in range(start, start + N):
            yield self[n % N]

    def invert(self):
        """
        Invert this cyclic word in place.
        """
        # This is builtin, starting from python 2.7
        for n in range(len(self) // 2):
            self[n], self[-1 - n] = self[-1 - n], self[n]
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
        abs_image = set(map(operator.abs, perm_dict.values()))
        if set(perm_dict) != abs_image:
            raise ValueError('Not a permutation!')
        for n in range(len(self)):
            x = self[n]
            self[n] = perm_dict.get(x, x) if x > 0 else -perm_dict.get(-x, -x)

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
                result.append((last_letter, count))
                count = 1
                last_letter = letter
        result.append((last_letter, count))
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
            for n in range(len(self)):
                complexity, Xordering = word.complexity(size, ordering, spin=n)
                if complexity < least:
                    least = complexity
                    minima = [(CyclicWord(word.spun(n)), Xordering)]
                elif complexity == least:
                    minima.append((CyclicWord(word.spun(n)), Xordering))
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
            sub = '%s -> ' % self.alphabet[x]
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


class Presentation:
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
        generators = ', '.join(self.alphabet[x] for x in self.generators)
        return 'generators: [%s]\nrelators: %s' % (generators, self.relators)

    def __len__(self):
        return sum(len(W) for W in self.relators)

    def __eq__(self, other):
        return (self.relators == other.relators and
                self.generators == other.generators)

    def whitehead_graph(self):
        transitions = []
        for relator in self.relators:
            for n in range(-1, len(relator) - 1):
                v, w = relator[n], -relator[n + 1]
                if v > w:
                    v, w = w, v
                transitions.append((v, w))

        Wh = nx.Graph()
        Wh.add_nodes_from(self.generators)
        Wh.add_nodes_from(-g for g in self.generators)
        for (u, v), m in Counter(transitions).items():
            Wh.add_edge(u, v, multiplicity=m)
        return Wh

    def find_reducers(self):
        Wh = self.whitehead_graph()
        reducers = []
        levels = []
        for x in self.generators:
            cut_size, partition = nx.minimum_cut(Wh, x, -x,
                                                 capacity='multiplicity')
            valence = Wh.degree(x, weight='multiplicity')
            length_change = cut_size - valence
            if length_change < 0:
                reducers.append((length_change, x, partition[0]))
        reducers.sort()
        return reducers

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
        >>> nx.is_planar(P.whitehead_graph())
        False
        >>> S = P.shorten()
        >>> print(S)
        generators: [A, B, C, D]
        relators: [AAAAABBAACCC, AAABBBCCCCC, AADCCD]
        >>> nx.is_planar(S.whitehead_graph())
        True
        >>> P = Presentation(['xyyxyyxy', 'xyy'])
        >>> P.shorten()
        generators: [X, Y]
        relators: [x, y]
        """
        result = Presentation(self.relators, self.generators)
        while True:
            reducers = result.find_reducers()
            if not reducers:
                break
            reduction, a, cut_set = reducers[0]
            result = result.whitehead_move(a, cut_set)
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

        # For each generator x we find a maximal (x,x^-1)-flow.  We
        # then construct the residual graph for the flow, which is the
        # diagraph D with a directed edge from y to z if there is
        # remaining capacity in that direction along the edge.
        # (Both directions may occur, giving rise to cycles of length 2.)
        # For each edge which is not saturated by the flow, we add a
        # cycle of length 2.
        #
        # We then form the condensation DAG C of strong components of
        # D. The transitively closed subsets of D determine all cuts,
        # which must be level transformations.  This generator yields
        # all level transformations of the form (x, X) where neither X
        # nor its complement has size 1.

        Wh = self.whitehead_graph()
        for x in self.generators:
            valence = Wh.degree(x, weight='multiplicity')
            flow_val, flow = nx.maximum_flow(Wh, x, -x, capacity='multiplicity')
            if flow_val < valence:
                raise ValueError('Presentation is not minimal.')

            D = nx.DiGraph()
            for u, v, m in Wh.edges(data='multiplicity'):
                f = flow[u][v]
                if f < m:
                    D.add_edge(u, v)

                f = flow[v][u]
                if f < m:
                    D.add_edge(v, u)

            C = nx.condensation(D)
            members = C.nodes('members')
            for B in nonempty_transitively_closed_subets(C):
                closure = set.union(*[members[b] for b in B])
                if 1 < len(closure) < Wh.number_of_nodes() - 1:
                    yield x, closure

    def level_orbit(self, verbose=False):
        """
        Generator for all presentations obtained from this one by
        length preserving Whitehead moves. Does a depth first search
        of the orbit.  If the verbose flag is True, yields a tuple:
        (parent, Whitehead move, result, canonical presentation).
        """
        S = self.signature()
        queue = deque([(None, None, self, S)])
        seen = {S}
        while queue:
            parent, move, pres, sig = queue.popleft()
            for a, A in pres.level_transformations():
                P = Presentation(pres.relators, pres.generators)
                P = P.whitehead_move(a, A)
                signature = P.signature()
                if signature not in seen:
                    WM = WhiteheadMove(a, A, pres.generators, self.alphabet),
                    queue.append((pres, WM, P, signature))
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
        while True:
            if len(queue[0].remaining) > 0:
                left = queue.popleft()
                for child in left.children():
                    queue.append(child)
            else:
                break
        relators = queue[0].presentation.relators
        ordering = queue[0].ordering
        generators = tuple(range(1, len(self.generators) + 1))
        return tuple(tuple(R.rewrite(ordering)) for R in relators), generators

    def magma_string(self):
        gens = sorted(self.alphabet[g] for g in self.generators)
        ans = 'Group<' + ', '.join(gens) + ' | '
        ans += ', '.join(R.verbose_string() for R in self.relators)
        return ans + '>'


class CanonizeNode:
    def __init__(self, presentation, remaining, ordering=[]):
        self.presentation = presentation
        self.generators = presentation.generators
        self.relators = presentation.relators
        self.remaining = remaining
        self.ordering = ordering

    def __repr__(self):
        return '%s\n%s' % (self.presentation, self.ordering)

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
