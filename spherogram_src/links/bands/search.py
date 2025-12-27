from ..links import Link
from .merge_links import (link_isotopy_classes,
                          are_isometric_as_links,
                          are_same_triangulations_and_links)
from .core import Band, add_one_band, banded_links, normalize_crossing_labels


def linking_nums_all_zero(link):
    """
    >>> linking_nums_all_zero(Link('L2a1'))  # Hopf
    False
    >>> linking_nums_all_zero(Link('L5a1'))  # Whitehead
    True
    """
    return all(x == 0 for row in link.linking_matrix() for x in row)


def is_unlink_exterior(manifold):
    """
    If there is clearly a hyperbolic structure present, we don't even
    look at the fundamental group as it takes some time to compute and
    simplify.
    """
    if manifold.solution_type() == 'all tetrahedra positively oriented':
        return False
    return manifold.fundamental_group().num_relators() == 0


def could_be_strongly_slice(link):
    """
    Checks several obstructions for a link L to be strongly slice,
    that is that all components of L simultaneously bound embedded
    disks in the 4-ball.  Returns False when L has been shown to not
    be strongly slice::

      sage: L = Link('L2a1')
      sage: could_be_strongly_slice(L)  #doctest: +SNAPPY
      False
    """
    if not linking_nums_all_zero(link):
        return False
    if link.signature() != 0:
        return False
    if not link.exterior().fox_milnor_test():
        return False
    return True


def remove_reidemeister_I(link):
    """
    Do Reidemeister I moves until none are possible.

    Adapted from simplify.basic_simplify
    """
    from spherogram.links.simplify import reidemeister_I

    to_visit = set(link.crossings)
    eliminated = set()
    while to_visit:
        crossing = to_visit.pop()
        elim, changed = reidemeister_I(link, crossing)
        assert not elim.intersection(changed)
        eliminated.update(elim)
        to_visit.difference_update(elim)
        to_visit.update(changed)

    success = len(eliminated) > 0

    # Redo the strand labels (used for DT codes)
    if success:
        component_starts = []
        for component in link.link_components:
            assert len(component) > 0
            if len(component) > 1:
                a, b = component[:2]
            else:
                a = component[0]
                b = a.next()
            if a.strand_label() % 2 == 0:
                component_starts.append(a)
            else:
                component_starts.append(b)
        link._build_components(component_starts)
    return success

def are_same_link(L0, L1, tries=10):
    """
    Try to determine whether the given links are isotopic *or* one is
    isotopic to the mirror of the other.

    If it returns True, then the links are definitely equivalent.  If
    it returns False, then they may still be the same.

    When both links are hyperbolic, this is function is robust.
    Otherwise, whether it succeeds is pretty stochastic.
    """
    L0, L1 = L0.copy(), L1.copy()
    for i in range(tries):
        if (L0.PD_code() == L1.PD_code() or
            are_isometric_as_links(L0, L1) or
            are_same_triangulations_and_links(L0, L1)):
            return True

        L0.backtrack(30)
        L0.simplify()
        L0.simplify('global')
        L1.backtrack(30)
        L1.simplify()
        L1.simplify('global')

    return False


def verify_ribbon_to_unknot(link, certificate):
    """
    >>> L = Link('K12n553')
    >>> cert = [[(3,1,4,0),(1,12,2,13),(11,2,12,3),(4,15,5,16),(14,5,15,6),
    ...          (6,13,7,14),(18,8,19,7),(8,18,9,17),(22,10,23,9),(10,22,11,21),
    ...          (16,19,17,20),(20,23,21,0)],'2c02_0_0',[(1,8,2,9),(7,2,8,3),
    ...          (3,6,4,7),(11,5,0,4),(5,11,6,10),(9,0,10,1)],'0900_0_-1','unknot']
    >>> verify_ribbon_to_unknot(L, cert)  #doctest: +SNAPPY
    True

    TODO: Sometimes doesn't work when some intermediate link is
    non-hyperbolic.  in this case a connected sum of two trefoils.
    The below should return True but might fail occasionally.

    >>> L = Link([(0,28,1,27),(10,1,11,2),(2,11,3,12),(36,3,37,4),(19,5,20,4),(5,23,6,22),(15,6,16,7),(34,8,35,7),(8,25,9,26),(26,9,27,10),(12,36,13,35),(13,20,14,21),(21,14,22,15),(23,17,24,16),(17,31,18,30),(31,19,32,18),(24,33,25,34),(28,0,29,37),(29,33,30,32)])
    >>> cert = [[(0,28,1,27),(10,1,11,2),(2,11,3,12),(36,3,37,4),(19,5,20,4),(5,23,6,22),(15,6,16,7),(34,8,35,7),(8,25,9,26),(26,9,27,10),(12,36,13,35),(13,20,14,21),(21,14,22,15),(23,17,24,16),(17,31,18,30),(31,19,32,18),(24,33,25,34),(28,0,29,37),(29,33,30,32)],'483d_0_0',[(5,2,6,3),(3,6,4,7),(1,4,2,5),(7,11,8,10),(11,9,0,8),(9,1,10,0)],'1603_0_1','unknot']
    >>> verify_ribbon_to_unknot(L, cert)  #doctest: +SNAPPY
    True
    """
    import snappy
    if certificate[-1] != 'unknot':
        try:
            snappy.RibbonLinks[certificate[-1]]
        except KeyError:
            raise ValueError('ribbon_cert does not end with known ribbon link')

    L = link.copy()
    E = L.exterior()
    for i in range(0, len(certificate) - 1, 2):
        L_cert = Link(certificate[i])
        E_cert = L_cert.exterior()

        if not are_same_link(L, L_cert):
            return False

        band = Band(certificate[i + 1])
        L = add_one_band(L_cert, band)
        # If L already has the correct PD code, don't try to simplify
        # as that can only mess things up.
        if i + 2 < len(certificate):
            next_link_PD_code = certificate[i + 2]
            if L.PD_code() != next_link_PD_code:
                L.simplify('global')
        else:
            L.simplify('global')
        L.unlinked_unknot_components = 0
        E = L.exterior()

    if certificate[-1] != 'unknot':
        E_cert = snappy.RibbonLinks[certificate[-1]]
        return are_isometric_as_links(E, E_cert)

    return len(L.link_components) == 0 or is_unlink_exterior(E)


def ribbon_concordant_links(link_or_manifold,
                            max_bands=1,
                            max_twists=2,
                            max_band_len=None,
                            paths='shortest',
                            stop_at_unlink=True,
                            only_return_unlink=False,
                            R1_only=False,
                            R1_R2_only=False,
                            certify=False,
                            print_progress=False,
                            use_ribbon_link_cache=True,
                            filter_for_plausibly_slice=True):
    """
    For the input link L0 returns all links obtained by ribbon
    concordances using bands generated by banded_links with the given
    parameters.  Assuming we're filtering for plausibly slice knots,
    this requires Sage::

      sage: L = Link('K6a3')
      sage: ans = ribbon_concordant_links(L, max_twists=1, certify=True)   #doctest: +SNAPPY
      sage: len(ans), ans['unknot'][1:]                                    #doctest: +SNAPPY
      (1, ['1201_0_1', 'unknot'])
      sage: M = Link('K12n553')
      sage: len(ribbon_concordant_links(M, max_bands=1, max_band_len=2, use_ribbon_link_cache=False))   #doctest: +SNAPPY
      2
      sage: ribbon_concordant_links(M, max_bands=2, max_band_len=2, use_ribbon_link_cache=False)        #doctest: +SNAPPY
      ['unknot']
      sage: ans = ribbon_concordant_links(M, max_bands=1, max_band_len=2, certify=True)   #doctest: +SNAPPY
      sage: ans['unknot'][-1]  #doctest: +SNAPPY
      'ribbon_1_6_e73be35b'

    If certify is True, it returns a dictionary whose keys are the
    links and whose values are the history of the bands and
    intermediate links.

    """
    import snappy

    def log_progress(message):
        if print_progress:
            print(message)

    if isinstance(link_or_manifold, Link):
        link = link_or_manifold.copy()
        manifold = link.exterior()
    else:
        manifold = link_or_manifold
        link = manifold.link()

        manifold = manifold.copy()

    normalize_crossing_labels(link)

    old_links = [(link, manifold)]
    certificates = {link:[link.PD_code()]}
    for i in range(max_bands):
        log_progress(f'Adding band {i} to {len(old_links)} link(s)')
        new_links = []
        for L0, E0 in old_links:
            log_progress(f'  starting {L0}\n{L0.PD_code()}\n')
            for L, spec in banded_links(L0, max_twists, max_band_len, paths):
                # Skip links with the wrong linking numbers
                if filter_for_plausibly_slice and not linking_nums_all_zero(L):
                    continue

                if R1_only:
                    remove_reidemeister_I(L)
                elif R1_R2_only:
                    L.simplify('basic')
                else:
                    L.simplify('global')
                normalize_crossing_labels(L)

                # cap off any trivial components
                L.unlinked_unknot_components = 0

                if len(L.link_components) == 0:
                    if 'unknot' not in certificates:
                        past = certificates[L0]
                        certificates['unknot'] = past + [spec, 'unknot']

                else:
                    E = L.exterior()
                    if is_unlink_exterior(E):
                        if 'unknot' not in certificates:
                            past = certificates[L0]
                            certificates['unknot'] = past + [spec, 'unknot']
                    elif only_return_unlink and i == max_band_len - 1:
                        pass
                    elif ((not filter_for_plausibly_slice or could_be_strongly_slice(L)) and
                          not are_isometric_as_links(E, E0)):
                        matched = False
                        if use_ribbon_link_cache and E.solution_type(enum=True) in {1, 2}:
                            F = snappy.RibbonLinks.identify(E, extends_to_link=True)
                            if F:
                                matched = True
                                if 'unknot' not in certificates:
                                    past = certificates[L0]
                                    certificates['unknot'] = past + [spec, F.name()]

                        if not matched:
                            past = certificates[L0]
                            certificates[L] = past + [spec, L.PD_code()]
                            new_links.append((L, E))

                if stop_at_unlink and 'unknot' in certificates:
                    if certify:
                        return {'unknot':certificates['unknot']}
                    else:
                        return ['unknot']


        old_links = link_isotopy_classes(new_links)

    final_links = [L for L, E in old_links]
    if 'unknot' in certificates:
        final_links.append('unknot')

    if certify:
        final_links = {L:certificates[L] for L in final_links}

    return final_links
