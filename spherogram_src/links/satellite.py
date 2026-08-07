from typing import Literal

from .links_base import Crossing, CrossingStrand, Link
from .tangles import Tangle, IdentityBraid, BraidTangle, RationalTangle, join_strands


def _blackboard_cable_tangle(link: Link, n: int) -> Tangle:
    """
    Return an unbuilt tangle constructed by cutting open `link` at the 0th strand of
    the first crossing in `link.crossings`, and `n`-cabling that tangle in a
    blackboard-framed way.
    """
    new_crossings = {
        (x, i, j): Crossing(label=str((x, i, j)))
        for x in link.crossings
        for i in range(n)
        for j in range(n)
    }

    for x in link.crossings:
        for i in range(n):
            for j in range(n - 1):
                new_crossings[(x, i, j)][2] = new_crossings[(x, i, j + 1)][0]
        for i in range(n - 1):
            for j in range(n):
                new_crossings[(x, i, j)][1] = new_crossings[(x, i + 1, j)][3]

    ccw_indices = {
        0: [(i, 0) for i in range(n)],
        1: [(n - 1, i) for i in range(n)],
        2: [(n - 1 - i, n - 1) for i in range(n)],
        3: [(0, n - 1 - i) for i in range(n)],
    }

    entry_crossing, entry_index = link.crossing_entries()[0]
    bottom_points = [
        (new_crossings[(entry_crossing, *ccw_indices[entry_index][i])], entry_index)
        for i in range(n)
    ]
    exit_crossing, exit_index = entry_crossing.adjacent[entry_index]
    top_points = [
        (new_crossings[(exit_crossing, *ccw_indices[exit_index][-i - 1])], exit_index)
        for i in range(n)
    ]

    for x0, i0 in link.crossing_entries()[1:]:
        x1, i1 = x0.adjacent[i0]
        for i in range(n):
            new_entry = (x0, *ccw_indices[i0][i])
            new_exit = (x1, *ccw_indices[i1][-i - 1])
            new_crossings[new_entry][i0] = new_crossings[new_exit][i1]

    return Tangle(
        boundary=n,
        crossings=list(new_crossings.values()),
        entry_points=bottom_points + top_points,
        build=False,
    )


def _pq_braid(p: int, q: int) -> Tangle:
    """
    Return a `BraidTangle` with index `p` which closes to a (p, q)-torus knot
    """
    if q == 0:
        return IdentityBraid(p)
    if q < 0:
        return -_pq_braid(p, -q)
    braid_word = list(range(p - 1, 0, -1)) * q
    return BraidTangle(braid_word, p)


def satellite(knot: Link, pattern: Tangle, twists: int = 0) -> Link:
    """
    Returns the satellite with the specified companion, pattern, and number of full
    twists.

    The pattern is specified by a `Tangle` with the same number of top and bottom
    boundary strands. Identifying the top and bottom strands gives the pattern knot
    in the solid torus. The resulting knot inherits its orientation from `pattern`,
    so that satelliting by a tangle with one strand oriented upwards gives back the
    companion. If the tangle does not have coherent orientations, so that the sign
    on the bottom strand is the opposite of the corresponding top strand, no
    orientation behavior is guaranteed.

    The denominator closure of `pattern` may be a link of more than one component.
    If so, satelliting by `pattern` gives a link of the same number of components.
    """
    if len(knot.link_components) + knot.unlinked_unknot_components > 1:
        raise ValueError("Only defined for knots, this has more components")
    if pattern.boundary[0] != pattern.boundary[1]:
        raise ValueError(
            "Invalid pattern, numbers of strands on the top and bottom must be equal"
        )

    pattern = pattern.copy()
    width = pattern.boundary[0]

    if knot.unlinked_unknot_components == 1:
        for i in range(width):
            join_strands(pattern.adjacent[i], pattern.adjacent[width + i])

        satellite_link = type(knot)(pattern.crossings + pattern.boundary_strands)
        satellite_link._rebuild(same_components_and_orientations=True)
        return satellite_link

    cable_tangle = _blackboard_cable_tangle(knot, width)
    twist_tangle = _pq_braid(width, width * (twists - knot.writhe()))

    for i in range(width):
        join_strands(cable_tangle.adjacent[i], twist_tangle.adjacent[width + i])
        join_strands(twist_tangle.adjacent[i], pattern.adjacent[width + i])
        join_strands(pattern.adjacent[i], cable_tangle.adjacent[width + i])

    # get start_orientations for the satellite from the pattern
    # first, from closed components of the pattern
    open_cpts = set(s.strand_component for s in pattern.boundary_strands)
    pattern_start_orns = pattern._start_orientations()
    closed_start_orns: list[CrossingStrand] = [
        so for i, so in enumerate(pattern_start_orns) if i not in open_cpts
    ]

    # then, from open components of the pattern
    # taking only one CrossingStrand from each component in the denominator closure
    cpts_in_closure: list[set[int]] = []
    for i in range(pattern.boundary[0]):
        bottom_index: int = pattern.boundary_strands[i].strand_component
        bottom_cpt = next(
            (j for j, cpt in enumerate(cpts_in_closure) if bottom_index in cpt), None
        )

        top_index: int = pattern.boundary_strands[width + i].strand_component
        top_cpt = next(
            (j for j, cpt in enumerate(cpts_in_closure) if top_index in cpt), None
        )

        if bottom_cpt is not None and top_cpt is not None:
            if bottom_cpt != top_cpt:
                cpts_in_closure[bottom_cpt].update(cpts_in_closure[top_cpt])
                del cpts_in_closure[top_cpt]
        elif bottom_cpt is not None:
            cpts_in_closure[bottom_cpt].add(top_index)
        elif top_cpt is not None:
            cpts_in_closure[top_cpt].add(bottom_index)
        else:
            cpts_in_closure.append({bottom_index, top_index})

    open_start_orns: list[CrossingStrand] = []
    for cpt_set in cpts_in_closure:
        strand_index = None
        for i in range(pattern.boundary[0]):
            if pattern.boundary_strands[i].strand_component in cpt_set:
                strand_index = i
                break
        assert strand_index is not None

        bd_strand = CrossingStrand(*cable_tangle.adjacent[width + strand_index])
        strand = bd_strand.rotate().opposite()
        if pattern.boundary_signs[strand_index] == 1:
            strand = strand.rotate(2)

        open_start_orns.append(strand)

    # an unbuilt link for orienting the crossings following the pattern
    intermediate_link = Link(
        cable_tangle.crossings
        + cable_tangle.boundary_strands
        + twist_tangle.crossings
        + twist_tangle.boundary_strands
        + pattern.crossings
        + pattern.boundary_strands,
        check_planarity=False,
        build=False,
    )

    for c in intermediate_link.crossings:
        c._clear()

    intermediate_link._orient_crossings(closed_start_orns + open_start_orns)

    return type(knot)(intermediate_link.crossings)


def cable(knot: Link, p: int, q: int) -> Link:
    """
    Returns the (`p`, `q`)-cable of `knot`.

    If `p` and `q` are coprime, the result is a knot. Otherwise, it may be a link
    of more than one component.
    """
    if len(knot.link_components) + knot.unlinked_unknot_components != 1:
        raise ValueError("Only defined for knots, this does not have one component")
    if knot.unlinked_unknot_components == 1:
        return type(knot)(f"T({p},{q})")

    tangle = _blackboard_cable_tangle(knot, p) * _pq_braid(p, q - p * knot.writhe())
    tangle.make_upward()
    for i in range(p):
        join_strands(tangle.adjacent[i], tangle.adjacent[p + i])

    return type(knot)(tangle.crossings + tangle.boundary_strands)


def whitehead_double(
    knot: Link,
    clasp_sign: Literal["positive", "negative"] = "positive",
    twists: int = 0,
) -> Link:
    """
    Returns the Whitehead double of `knot`.

    The sign of the clasp of the pattern can be specified as `"positive"` or
    `"negative"`. You may add some number of full twists with the `twists`
    parameter.
    """
    if clasp_sign not in ["positive", "negative"]:
        raise ValueError("Invalid clasp sign")

    pattern = RationalTangle(-2 if clasp_sign == "positive" else 2)
    return satellite(knot, pattern, twists)
