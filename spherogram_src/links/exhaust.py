from .links_base import Strand, Crossing, Link
import random
import collections


def insert_space(point_dict, i):
    """
    Insert two new points which will be labeled i and i + 1.
        """
    ans = dict()
    for a, v in point_dict.items():
        if a >= i:
            a = a + 2
        ans[a] = v
    return ans


def remove_space(point_dict, i):
    """
    Remove the points i and i + 1
    """
    ans = dict()
    for a, v in point_dict.items():
        if a < i:
            ans[a] = v
        if a > i + 1:
            ans[a - 2] = v
    return ans


class Event():
    """
    There are three kinds of events:

    * ('cup', a, b, c=None): Introduce a new oriented strand with
      endpoints a and b oriented from a to b.  The optional argument c
      is the index of the link component in the final link.

    * ('cap', a, b): Join these strands by a cap.  The order of a and
      b is not important.

    * ('cross', a, b): Cross strand a *over* strand b.
    """

    def __init__(self, kind, a, b, c=None):
        if kind not in ['cup', 'cap', 'cross']:
            raise ValueError('Event kind must be cup, cap, or cross')
        self.kind = kind
        if abs(a - b) != 1 or a < 0 or b < 0:
            raise ValueError('Invalid strand info')
        self.a, self.b, self.c = a, b, c
        self.min = min(a, b)
        self.max = max(a, b)

    def __repr__(self):
        ans = (self.kind, self.a, self.b)
        if self.c is not None:
            ans += (self.c,)
        return repr(ans)


class MorseEncoding():
    """
    A MorseEncoding is a concrete encoding of a Morse diagram of an
    oriented link in the sprit of:

    http://katlas.org/wiki/MorseLink_Presentations

    as a sequence of events.  At each stage, there are n active strand
    ends labeled in order [0, 1, ... , n - 1].  The events can be
    specified either as instances of Event or as tuples as shown
    below.

    >>> events = 2*[('cup', 0, 1)] + 2*[('cross', 1, 2)] + 2*[('cross', 3, 2)]
    >>> events += [('cap', 1, 2), ('cap', 0, 1)]
    >>> me = MorseEncoding(events)
    >>> me.events[1:3]
    [('cup', 0, 1), ('cross', 1, 2)]
    >>> me.width
    4
    >>> me_copy = MorseEncoding(me.events)
    >>> L = me.link()
    >>> L.exterior().identify()[0]  #doctest: +SNAPPY
    m004(0,0)
    >>> U = MorseEncoding([('cup', 0, 1), ('cap', 0, 1)]).link()
    >>> U.exterior().fundamental_group().num_generators()  #doctest: +SNAPPY
    1
    """

    def __init__(self, events):
        self.events = []
        for event in events:
            if not isinstance(event, Event):
                event = Event(*event)
            self.events.append(event)
        self._check_and_set_width()

    def _check_and_set_width(self):
        width = 0
        max_width = 0
        for event in self.events:
            if event.kind == 'cup':
                assert event.max < width + 2
                width += 2
                max_width = max(width, max_width)
            elif event.kind == 'cap':
                assert event.max < width
                width += -2
            elif event.kind == 'cross':
                assert event.max < width
        assert width == 0
        self.width = max_width

    def link(self):
        active = dict()
        crossings = []
        for event in self.events:
            if event.kind == 'cup':
                S = Strand()
                crossings.append(S)
                active = insert_space(active, event.min)
                active[event.a] = S[0]
                active[event.b] = S[1]
            elif event.kind == 'cap':
                S = Strand()
                crossings.append(S)
                S[0] = active[event.a]
                S[1] = active[event.b]
                active = remove_space(active, event.min)
            elif event.kind == 'cross':
                C = Crossing()
                crossings.append(C)
                if event.a < event.b:
                    C[3] = active[event.a]
                    C[0] = active[event.b]
                    active[event.a] = C[2]
                    active[event.b] = C[1]
                else:
                    C[3] = active[event.a]
                    C[2] = active[event.b]
                    active[event.a] = C[0]
                    active[event.b] = C[1]
        return Link(crossings)

    def __iter__(self):
        return self.events.__iter__()


class BiDict():
    """
    A bijective mapping from range(n) to a set of hashable non-integers.

    >>> bd = BiDict({0:'a', 1:'b', 2:'c', 3:'d'})
    >>> bd[0], bd['b']
    ('a', 1)
    >>> bd._check()
    True
    >>> bd.insert_space(4, 2)
    >>> bd._check()
    False
    >>> bd[4] = 'e'
    >>> bd['f'] = 5
    >>> bd._check()
    True
    >>> bd
    BiDict({0: 'a', 1: 'b', 2: 'c', 3: 'd', 4: 'e', 5: 'f'})
    >>> bd.insert_space(3, 2)
    >>> bd[3], bd['y'] = 'x', 4
    >>> bd
    BiDict({0: 'a', 1: 'b', 2: 'c', 3: 'x', 4: 'y', 5: 'd', 6: 'e', 7: 'f'})
    >>> bd.pop(1); bd.pop('c'); bd.pop('e'); bd
    BiDict({0: 'a', 1: 'x', 2: 'y', 3: 'd', 4: 'f'})
    >>> (0 in bd, 'x' in bd, 'z' in bd, 10 in bd)
    (True, True, False, False)
    >>> bd.insert_space(-1, 2); bd[0] = 'u'; bd[1] = 'v'; bd
    BiDict({0: 'u', 1: 'v', 2: 'a', 3: 'x', 4: 'y', 5: 'd', 6: 'f'})
    >>> sorted(bd.values())
    ['a', 'd', 'f', 'u', 'v', 'x', 'y']
    """

    def __init__(self, int_to_set_dict):
        self.n = n = len(int_to_set_dict)
        assert sorted(int_to_set_dict) == list(range(n))
        self.int_to_set = int_to_set_dict
        self.set_to_int = {v: k for k, v in int_to_set_dict.items()}

    def __getitem__(self, index):
        if isinstance(index, int):
            return self.int_to_set[index]
        elif index in self.set_to_int:
            return self.set_to_int[index]

    def __setitem__(self, index, value):
        if isinstance(index, int):
            assert not isinstance(value, int)
            the_int, the_set = index, value
        else:
            assert isinstance(value, int)
            the_int, the_set = value, index
        assert 0 <= the_int < self.n
        self.int_to_set[the_int] = the_set
        self.set_to_int[the_set] = the_int

    def __contains__(self, value):
        return (value in self.set_to_int) or (value in self.int_to_set)

    def pop(self, index):
        """
        Remove the given item and shift the indices downwards when necessary
        to close the gap.
        """
        if isinstance(index, int):
            the_int = index
            the_set = self.int_to_set[the_int]
        else:
            the_set = index
            the_int = self.set_to_int[the_set]
        self.int_to_set.pop(the_int)
        self.set_to_int.pop(the_set)

        def shift(j):
            return j if j < the_int else j - 1
        self.int_to_set = {shift(k): v for k, v in self.int_to_set.items()}
        self.set_to_int = {k: shift(v) for k, v in self.set_to_int.items()}
        self.n += -1

    def insert_space(self, i, n):
        """
        Shift indices upwards when necessary so that the n slots

        i, i + 1, ... , i + n - 1

        are unassigned.
        """
        assert isinstance(i, int)

        def shift(j):
            return j if j < i else j + n
        self.n = self.n + n
        self.int_to_set = {shift(k): v for k, v in self.int_to_set.items()}
        self.set_to_int = {k: shift(v) for k, v in self.set_to_int.items()}

    def values(self):
        return self.int_to_set.values()

    def __repr__(self):
        itos = self.int_to_set
        keys = sorted(itos)
        items = ['%d: %s' % (k, repr(itos[k])) for k in keys]
        return 'BiDict({' + ', '.join(items) + '})'

    def __len__(self):
        return len(self.int_to_set)

    def _check(self):
        return sorted(self.int_to_set) == list(range(self.n))


def is_range(L):
    """
    >>> is_range([2, 3, 4]), is_range([2, 3, 5])
    (True, False)
    """
    return L == list(range(min(L), max(L) + 1))


class Frontier(BiDict):
    def overlap_indices(self, crossing):
        neighbors = [cs.opposite() for cs in crossing.crossing_strands()]
        return [self[ns] for ns in neighbors if ns in self]

    def overlap_is_consecutive(self, crossing):
        overlap = self.overlap_indices(crossing)
        return len(overlap) > 0 and is_range(overlap)

    def biggest_all_consecutive_overlap(self):
        """
        Return a random crossing from among those with the maximal possible
        overlap.
        """
        overlap_indices = collections.defaultdict(list)
        for i, cs in self.int_to_set.items():
            overlap_indices[cs.opposite()[0]].append(i)
        possible = []
        for crossing, overlap in overlap_indices.items():
            overlap = sorted(overlap)
            if is_range(overlap):
                possible.append((len(overlap), min(overlap), crossing))
        max_overlap = max(possible)[0]
        good_choices = [pos for pos in possible if pos[0] == max_overlap]
        return random.choice(good_choices)


class MorseExhaustion():
    """
    An exhaustion of a link where crossings are added in one-by-one
    so that the resulting tangle is connected at every stage.

    Starting at the given crossing, it uses a greedy algorithm to try
    to minimize the sizes of the frontiers of the intermediate tangles.

    The answer is returned as a sequence of events describing a
    MorseEncoding.

    If no initial crossing is specified, one is chosen at random.

    >>> L = Link('L2a1')
    >>> mexhaust = MorseExhaustion(L, L.crossings[0])
    >>> mexhaust
    [('cup', 0, 1), ('cup', 0, 1), ('cross', 1, 2), ('cross', 1, 2), ('cap', 0, 1), ('cap', 0, 1)]
    >>> me = MorseEncoding(mexhaust)
    >>> me.link().exterior().fundamental_group().relators() #doctest: +SNAPPY
    ['abAB']

    >>> K = Link([[0, 0, 1, 1]])  # Unknot
    >>> MorseExhaustion(K)
    [('cup', 0, 1), ('cup', 0, 1), ('cross', 1, 2), ('cap', 0, 1), ('cap', 0, 1)]
    """
    def __init__(self, link, crossing=None):
        events = []
        if link.crossings:
            if crossing is None:
                crossing = random.choice(link.crossings)
            crossings = [crossing]
            events = [('cup', 0, 1), ('cup', 0, 1), ('cross', 1, 2)]
            css = crossing.crossing_strands()
            frontier = Frontier({0: css[3], 1: css[2], 2: css[1], 3: css[0]})
            frontier_lengths = [4]
            if len(link.crossings) == 1:
                events += [('cap', 0, 1), ('cap', 0, 1)]
        else: # Only unlinked unknotted components
            crossings = []
            frontier_lengths = []

        while len(crossings) < len(link.crossings):
            overlap, i, C = frontier.biggest_all_consecutive_overlap()
            cs = frontier[i]
            cs_opp = cs.opposite()
            assert C not in crossings
            crossings.append(C)
            if overlap == 1:
                i = frontier[cs]
                events.append(('cup', i + 1, i + 2))
                if cs_opp.strand_index in {1, 3}:
                    events.append(('cross', i, i + 1))
                else:
                    events.append(('cross', i + 1, i))
                frontier.insert_space(i + 1, 2)
                for s in range(3):
                    frontier[i + s] = cs_opp.rotate(-(s + 1))
            elif overlap == 2:
                if cs_opp.strand_index in {1, 3}:
                    events.append(('cross', i, i + 1))
                else:
                    events.append(('cross', i + 1, i))
                for s in range(2):
                    frontier[i + s] = cs_opp.rotate(-(s + 1))
            elif overlap == 3:
                if cs_opp.strand_index in {1, 3}:
                    events.append(('cross', i, i + 1))
                else:
                    events.append(('cross', i + 1, i))
                events.append(('cap', i + 1, i + 2))
                frontier.pop(i + 2)
                frontier.pop(i)
                frontier[i] = cs_opp.rotate(-1)
            else:
                assert overlap == 4
                if cs_opp.rotate().strand_index in {1, 3}:
                    events.append(('cross', i + 1, i + 2))
                else:
                    events.append(('cross', i + 2, i + 1))
                events += 2 * [('cap', i, i + 1)]
                for a in range(4):
                    frontier.pop(i)
            assert frontier_lengths[-1] + 4 - 2 * overlap == len(frontier)
            assert frontier._check()
            frontier_lengths.append(len(frontier))

        c = link.unlinked_unknot_components
        events += c*[('cup', 0, 1), ('cap', 0, 1)]
        frontier_lengths += c*[2, 0]

        self.link = link
        self.crossings = crossings
        self.frontier_lengths = frontier_lengths
        self.events = events
        self.width = max(frontier_lengths) // 2

    def __repr__(self):
        return repr(self.events)

    def __iter__(self):
        return self.events.__iter__()


def good_exhaustion(link, max_tries=20):
    """
    Uses a random search to try to find an Exhaustion with small width.

    >>> ge = good_exhaustion(Link('K4a1'))
    >>> ge.width
    2
    """
    E_best = None
    crossings = list(link.crossings)
    tries = 0
    while tries < max_tries:
        random.shuffle(crossings)
        for C in crossings:
            E = MorseExhaustion(link, C)
            if E_best is None or E.width < E_best.width:
                E_best = E
            tries += 1
            if tries >= max_tries:
                break
    return E_best


def test_morse_machine(link):
    E = link.exterior()
    exhaust = MorseExhaustion(link, link.crossings[0])
    encoding = MorseEncoding(exhaust)
    new_link = encoding.link()
    new_E = new_link.exterior()
    return E.is_isometric_to(new_E)


def test_many(N):
    for i in range(N):
        M = snappy.HTLinkExteriors().random()
        if M.solution_type() == 'all tetrahedra positively oriented':
            assert test_morse_machine(M.link())


def morse_encoding_from_zs_hfk(link):
    """
    >>> me = morse_encoding_from_zs_hfk(Link('K8n1'))
    >>> me.width
    6
    """
    from knot_floer_homology import pd_to_morse
    pd = 'PD: ' + repr(link.PD_code())
    zs_morse = pd_to_morse(pd)
    ans = MorseEncoding(zs_morse['events'])
    assert zs_morse['girth'] == ans.width
    return ans


def orient_pres_isometric(A, B):
    for iso in A.is_isometric_to(B, True):
        mat = iso.cusp_maps()[0]
        if mat.det() == 1:
            return True
    return False


def test_zs_hfk(crossings, how_many):
    from networkx.algorithms import approximation
    for _ in range(how_many):
        K = spherogram.random_link(crossings, num_components=1,
                                   initial_map_gives_link=True,
                                   consistent_twist_regions=True)

        E = K.exterior()
        if not E.solution_type().startswith('all tet'):
            continue

        exhaust = good_exhaustion(K, 100)
        encoding0 = MorseEncoding(exhaust)
        encoding1 = morse_encoding_from_zs_hfk(K)

        M0 = encoding0.link().exterior()
        M1 = encoding1.link().exterior()
        assert orient_pres_isometric(E, M0)
        print(encoding0.width, encoding1.width, orient_pres_isometric(E, M1))


if __name__ == '__main__':
    import doctest
    print(doctest.testmod())
