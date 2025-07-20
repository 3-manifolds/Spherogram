"""

Given a collection of hyperbolic links, select one from each isotopy
class, modulo mirror image.

"""

import snappy
import collections
from snappy import Manifold


def pos_tets(manifold):
    return manifold.solution_type() == 'all tetrahedra positively oriented'


def degen_tets(manifold):
    decent = ['all tetrahedra positively oriented',
              'contains negatively oriented tetrahedra']
    return manifold.solution_type() not in decent


def are_isometric_as_links(A, B):
    if degen_tets(A) or degen_tets(B):
        return False
    try:
        isos = A.is_isometric_to(B, return_isometries=True)
        return any(iso.extends_to_link() for iso in isos)
    except RuntimeError:
        return False


def link_isotopy_classes(links_with_manifolds):
    """
    >>> mflds = [Manifold('L13n5325'), Manifold('L13n5834')]
    >>> links_w_mflds = [(M.link(), M) for M in mflds]
    >>> len(link_isotopy_classes(links_w_mflds))  # These are actually duplicates!
    1
    >>> mflds = snappy.HTLinkExteriors[:3]
    >>> links_w_mflds = [(M.link(), M) for M in mflds]
    >>> links = link_isotopy_classes(links_w_mflds)
    >>> [pos_tets(M) for L, M in links]
    [True, True, False]
    """
    if len(links_with_manifolds) == 0:
        return []
    isom_sigs = dict()
    troublesome = []
    for (L, M) in links_with_manifolds:
        try:
            sig = M.isometry_signature(of_link=True)
            if sig not in isom_sigs:
                isom_sigs[sig] = (L, M)
        except RuntimeError:
            troublesome.append((L, M))

    return list(isom_sigs.values()) + troublesome


def coarse_link_isotopy_classes(links_with_manifolds):
    data = collections.defaultdict(list)
    for (L, M) in links_with_manifolds:
        poly = L.sage_link().homfly_polynomial()
        data[poly].append((L, M))

    def crossings(link_w_mfld):
        return len(link_w_mfld[0].crossings)
    
    ans = [min(pairs, key=crossings) for pairs in data.values()]
    ans.sort(key=crossings)
    return ans
