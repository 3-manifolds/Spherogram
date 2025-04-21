"""
Implements the upper bounds on bridge number from:

1. R. Blair, A. Kjuchukova, and E. Pfaff,
   The Plain Sphere Number of a Link, Preprint 2025.
   https://arxiv.org/abs/2504.10517

2. R. Blair, A. Kjuchukova, R. Velazquez, and P. Villanueva,
   Wirtinger systems of generators of knot groups,
   Comm. in Analysis and Geometry, 28 (2020), 243-262.
   https://dx.doi.org/10.4310/CAG.2020.v28.n2.a2
"""


from .simplify import dual_graph_as_nx
from .links_base import CrossingStrand, Link
import networkx as nx
from packaging.version import Version


# Finding simple_cycles in *undirected* graphs was only added in
# NetworkX 3.1, circa March 2023.  For earlier vesions of NetworkX, we
# use the local backport.

if Version(nx.__version__) >= Version('3.1'):
    simple_cycles = nx.simple_cycles
else:
    from .nx_simple_cycles import simple_cycles


def number_strands(link):
    """
    Numbers the strands in the diagram and returns a dictionary that
    maps each CrossingStrand or CrossingEntryPoint to the
    corresponding strand number.

    >>> K = Link('K8n1')
    >>> len(set(number_strands(K).values()))
    8
    """
    to_meridian = dict()
    for i, C in enumerate(link.crossings):
        cs = CrossingStrand(C, 0)
        to_meridian[cs] = i
        cs = cs.opposite()
        while cs.strand_index != 2:
            to_meridian[cs] = i
            cs = cs.rotate(2)
            to_meridian[cs] = i
            cs = cs.opposite()
        to_meridian[cs] = i

    return to_meridian


class MeridianRelation:
    """
    For a relator consisting of a product of meridians, a MeridianRelation
    records:

    1) The strand number of each (unoriented) meridian that appears.

    2) Which meridians appear at least twice.

    >>> r = MeridianRelation([1, 1, 2, 3, 4])
    >>> r.applies({1, 2, 4})
    (True, 3)
    >>> r.applies({2, 3, 4, 5})
    (False, None)
    >>> r.applies({1})
    (False, None)
    """
    def __init__(self, relation):
        self.meridians = set(relation)
        self.appears_twice = {r for r in relation if relation.count(r) > 1}
        self.num_meridians = len(self.meridians)

    def applies(self, current_meridians):
        """
        The relation applies when there is a meridan M in it that
        appears only once, and the other meridians in the relation are
        in current_meridians.

        Returns (relation can be applied, None or the new meridian M)
        """
        if len(current_meridians) < self.num_meridians - 1:
            return False, None

        outside = None
        for m in self.meridians:
            if m not in current_meridians:
                if (m in self.appears_twice) or (outside is not None):
                    return False, None
                outside = m
        return True, outside


def meridian_relations(link):
    """
    Identifies all simple cycles in the dual graph to the given planar
    diagram and returns the corresponding MeridianRelations, sorted by
    length.

    >>> K = Link('9_42')
    >>> rels = meridian_relations(K)
    >>> len(rels)
    74
    """
    to_meridian = number_strands(link)
    G = dual_graph_as_nx(link)
    relations = []
    for cycle in simple_cycles(G):
        relation = []
        n = len(cycle)
        for i in range(n):
            a = cycle[i]
            b = cycle[(i + 1) % n]
            meridians = [to_meridian[cs]
                         for cs in G.edges[a, b]['interface'].values()]
            assert len(set(meridians)) == 1
            m = meridians[0]
            relation.append(m)

        relations.append(relation)

    relations = [MeridianRelation(rel)
                 for rel in sorted(relations, key=lambda x:(len(x), x))]

    return relations


def crossing_relations(link):
    """
    Returns the meridian relations corresponding to the crossings of
    the diagram.

    >>> K = Link('9_42')
    >>> len(crossing_relations(K))
    9
    """
    to_meridian = number_strands(link)
    relations = []
    for C in link.crossings:
        rel = [to_meridian[cs] for cs in C.crossing_strands()]
        relations.append(rel)
    return [MeridianRelation(rel) for rel in relations]


class MeridianSet:
    """
    A set of meridians generated from a "seed" of such, together with a
    collection of meridian relations.
    """
    def __init__(self, total_meridians, relations, meridians=None, seed=None):
        self.n = total_meridians
        self.remaining_relations = relations[:]
        self.exhausted_relations = []
        if meridians is not None and seed is not None:
            self.meridians = set(meridians)
            self.seed = set(seed)
        else:
            self.seed = set()
            self.meridians = set()

    def add_to_seed(self, meridian):
        """
        Add a new meridian to the seed and apply meridian relations to
        expand as much as possible.
        """
        self.seed.add(meridian)
        self.meridians.add(meridian)
        self.saturate()

    def done(self):
        return len(self.meridians) == self.n

    def expand(self):
        """
        Helper method for saturate: expand the current set of
        meridians with a single pass through the remaining relations.
        """
        start_size = len(self.meridians)
        unused = []
        for rel in self.remaining_relations:
            progress, new_meridian = rel.applies(self.meridians)
            if new_meridian is not None:
                self.meridians.add(new_meridian)
            if not progress:
                unused.append(rel)
        self.remaining_relations = unused
        return len(self.meridians) > start_size

    def saturate(self):
        """
        Expand the current set of meridians as much as possible using
        the relations.
        """
        start_size = len(self.meridians)
        while self.expand():
            pass
        return len(self.meridians) > start_size

    def __repr__(self):
        return f'<MS of size {len(self.meridians)} out of {self.n}  with seed {sorted(self.seed)}>'

    def descendents(self):
        """
        Return all MeridianSets created by adding an unused meridian to
        the current seed.

        To avoid duplication, we only add new meridians who whose
        index is larger than any currently in the seed.
        """
        start = 0 if len(self.seed) == 0 else max(self.seed) + 1
        for m in range(start, self.n):
            if m not in self.meridians:
                child = MeridianSet(self.n, self.remaining_relations,
                                    self.meridians, self.seed)
                child.add_to_seed(m)
                yield child


def bridge_upper_bound(link, method='plain sphere', return_meridians=False):
    """
    Computes an upper bound on the bridge number of the given link.
    By default, it computes the plain sphere number rho(D) of the
    diagram D from [BKP2025]_, but with ``method='wirtinger'``, it
    computes the weaker Wirtinger number omega(D) from [BKVV2020]_.
    The latter is significantly faster to compute.

    If ``return_meridians=True`` is set, returns a list of meridians
    generating the knot group, where each meridian is specified by the
    index of the crossing for which it is the 0 input strand.

    >>> K = Link('K14n1527')
    >>> bridge_upper_bound(K)
    3
    >>> bridge_upper_bound(K, method='wirtinger')
    4
    >>> bridge_upper_bound(K, return_meridians=True)
    [2, 3, 9]

    .. [BKP2025] R. Blair, A. Kjuchukova, and E. Pfaff,
                 *The Plain Sphere Number of a Link.*
                 https://arxiv.org/abs/2504.10517

    .. [BKVV2020] R. Blair, A. Kjuchukova, R. Velazquez, and P. Villanueva,
                  *Wirtinger systems of generators of knot groups.*
                  https://dx.doi.org/10.4310/CAG.2020.v28.n2.a2
    """
    n = len(link.crossings)
    if method == 'plain sphere':
        rels = meridian_relations(link)
    elif method == 'wirtinger':
        rels = crossing_relations(link)
    else:
        raise ValueError("Available methods are 'plain sphere' and 'wirtinger'")

    empty = MeridianSet(n, rels)
    current = [empty]
    next_round = []
    while True:
        for S in current:
            for D in S.descendents():
                if D.done():
                    if return_meridians:
                        return sorted(D.seed)
                    else:
                        return len(D.seed)
                next_round.append(D)
        current = sorted(next_round, key=lambda x:len(x.meridians), reverse=True)
        next_round = []


if __name__ == '__main__':
    import doctest
    print(doctest.testmod())
