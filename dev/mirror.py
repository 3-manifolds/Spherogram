from spherogram import Link, Crossing
from spherogram.links.links import CrossingEntryPoint
import snappy

L = Link('9^3_12')
#L = Link('4a1')

def mirror(link):
    """
    Create the mirror image of the link, preserving link orientations and
    component order.
    """
    # Basically, we just mirror every crossing, but the particular data 
    # structures used make this a little involved.
    
    new_crossings = dict()
    for C in link.crossings:
        C_new = Crossing(label=C.label)
        C_new.sign = -C.sign
        new_crossings[C] = C_new

    def convert(C, c):
        """
        Go from a crossing to the mirror, which requires the a rotation of the
        entry points; the direction of rotation depends on the sign of
        the crossing.
        """
        return new_crossings[C], (c + C.sign) % 4
        
    for A in link.crossings:
        entry_points = [CEP.entry_point for CEP in A.entry_points()]
        for a in entry_points:
            B, b = A.adjacent[a]
            B_new, b_new = convert(B, b)
            B_new[b_new] = convert(A, a)

    new_link = Link(new_crossings.values(),
                    check_planarity=False, build=False)

    # Build the link components, starting in the same place as
    # the original link.
    
    component_starts = []
    for component in link.link_components:
        C_new, c_new = convert(*component[0])
        component_starts.append( CrossingEntryPoint(C_new, c_new) )

    new_link._build_components(component_starts)
    return new_link
        

def test_link(L):
    Lbar = mirror(L)
    N = L.exterior()
    Nbar = Lbar.exterior()
    assert len([i for i in N.is_isometric_to(Nbar, True)
                if i.cusp_images() == range(len(L.link_components))]) > 0
    assert L.signature() + Lbar.signature() == 0
    assert [len(C) for C in L.link_components] == [len(C) for C in Lbar.link_components]

def basic_test():
    for M in snappy.LinkExteriors():
        if M.solution_type().startswith('all tetra'):
            L = Link(M.name())
            test_link(L)

        

