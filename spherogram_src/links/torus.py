"""
Torus knots
"""
from .links import Crossing, Link


def torus_knot(name: str) -> Link:
    """
    Return a `(p,q)`-torus knot, as an instance of the ``Link`` class.

    If exactly one of `p` and `q` is negative, it returns the mirror of
    `T(|p|, |q|)`.
    """
    p, q = map(int, name[2:-1].split(','))

    if p == 0 or q == 0:
        raise ValueError("torus_knot(p,q) requires non zero p and q")
    to_mirror = p * q < 0
    p, q = abs(p), abs(q)
    if p == 2:
        our_crossings = [Crossing(i) for i in range(q)]
        if q > 1:
            # set up conditions true for all two strand situations if q > 1
            our_crossings[0][0] = our_crossings[q - 1][1]
            our_crossings[0][3] = our_crossings[q - 1][2]
            our_crossings[0][1] = our_crossings[1][0]
            our_crossings[0][2] = our_crossings[1][3]
            # set up in between crossings
            for i in range(1, q - 1):
                our_crossings[i][1] = our_crossings[i + 1][0]
                our_crossings[i][2] = our_crossings[i + 1][3]
            link = Link(our_crossings)
            return link.mirror() if to_mirror else link

        our_crossings[0][0] = our_crossings[0][1]
        our_crossings[0][2] = our_crossings[0][3]
        link = Link(our_crossings)
        return link.mirror() if to_mirror else link

    our_crossings = {(i, j): Crossing((i, j))
                     for i in range(q) for j in range(p - 1)}

    # set up connecting ends
    our_crossings[(0, 0)][3] = our_crossings[(q - 1, 0)][2]
    our_crossings[(0, p - 2)][0] = our_crossings[(q - 1, p - 2)][1]

    # middle strands of connecting ends
    for i in range(p - 2):
        our_crossings[(0, i)][0] = our_crossings[(q - 1, i + 1)][2]

    if q > 1:
        # set up side connections
        for i in range(q - 1):
            our_crossings[(i, 0)][2] = our_crossings[(i + 1, 0)][3]
            our_crossings[(i, p - 2)][1] = our_crossings[(i + 1, p - 2)][0]

        # set up connections between crossings
        for i in range(p - 2):  # was p - 3
            for j in range(q):
                our_crossings[(j, i)][1] = our_crossings[(j, i + 1)][3]

        for i in range(1, p - 1):  # was p - 2
            for j in range(q - 1):
                our_crossings[(j, i)][2] = our_crossings[(j + 1, i - 1)][0]

    link = Link(list(our_crossings.values()))
    return link.mirror() if to_mirror else link
