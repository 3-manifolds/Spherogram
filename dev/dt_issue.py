import snappy, spherogram, plink
import spherogram.links.orthogonal
from nsnappytools import appears_hyperbolic
from sage.all import *


links = snappy.HTLinkExteriors()
for i in range(1):
    dt = links.random().DT_code()
    dtm = [tuple([-d for d in comp]) for comp in dt]
    dtm = dtm
    for d in [dt, dtm]:
        dtc = spherogram.DTcodec(d)
        L = dtc.link()
        assert d==L.DT_code()

dt = [(8, -10, 14), (12, -20, 18, 22, 4, 24, -2, 6, 26, 16)]
dtc = spherogram.DTcodec(dt)
L = dtc.link()

def asymmetric_link_DTs(N):
    found = 0
    while found < N:
        M = snappy.HTLinkExteriors.random()
        if appears_hyperbolic(M) and M.symmetry_group().order() == 1:
            found += 1
            yield M.DT_code()

def matches_peripheral_data(isom):
    Id = matrix(ZZ, [[1,0], [0,1]])
    cusp_perm = isom.cusp_images()
    cusp_maps = isom.cusp_maps()
    n = len(cusp_perm)
    return cusp_perm == range(n) and cusp_maps == (n*[Id])

def manifolds_match(M0, M1):
    isoms = M0.is_isometric_to(M1, True)
    assert len(isoms) == 1
    return matches_peripheral_data(isoms[0])
            
def test_DT(dt, M2=None):
    if M2 is None:
        M2 = snappy.Manifold()
    dtc = spherogram.DTcodec(dt)
    L = dtc.link()
    M0, M1 = dtc.exterior(), L.exterior()
    L.view(M2.LE)
    #M2.LE.sorted_components()
    M2.LE.callback()
    return manifolds_match(M0, M1) and manifolds_match(M1, M2)

def fill_vols(M):
    vols = []
    n = M.num_cusps()
    for i in range(n):
        N = M.without_hyperbolic_structure()
        fill = n*[(1,0)]
        fill[i] = (0,0)
        N.dehn_fill(fill)
        N = N.filled_triangulation().with_hyperbolic_structure()
        print(N.num_cusps())
        vols.append(N.volume())
    return vols

def plink_component_order_test(N):
    M2 = snappy.Manifold()
    return set([test_DT(dt, M2) for dt in asymmetric_link_DTs(N)])

