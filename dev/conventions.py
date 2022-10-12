"""
Various Trefoils:

Start with what spherogram calls positive crossings and stack
vertically with C at the top and A at the bottom.

OzSz and Lickorish agree that this is a +1 crossing

This is the *right-handed* trefoil according to Jake Rasmussen and
Wikipedia.

>>> A, B, C = Crossing(), Crossing(), Crossing()
>>> A[2], A[1] = B[3], B[0]
>>> B[2], B[1] = C[3], C[0]
>>> C[2], C[1] = A[3], A[0]
>>> T0 = Link([A, B, C])
>>> T0.signature()
-2
>>> T0.jones_polynomial()
q^2 + q^6 - q^8
>>> E = T0.exterior()
>>> E.dehn_fill((1, 1))
>>> gap(E.fundamental_group()).Order()
120


As per OzSz, HFK in Alexander degree 0 is supported in Maslov grading
(signature/2), which is the case now:


>>> HFK = T0.knot_floer_homology()
>>> HFK['tau']
1
>>> HFK['ranks']
{(-1, -2): 1, (0, -1): 1, (1, 0): 1}

More generally "Positive knots have negative signature" according to
Jake, Rolfsen, and Przytycki's article with that title:

https://arxiv.org/pdf/0905.0922.pdf

Some other ways of getting the trefoil:


>>> T1 = RationalTangle(1, 3).denominator_closure()
>>> T1.signature()
2

>>> T2 = RationalTangle(3, 1).numerator_closure()
>>> T2.signature()
-2

>>> T3 = Link('T(3, 2)')
>>> T3.signature()
-2

>>> T4 = Link(braid_closure=[1, 1, 1])
>>> T4.signature()
-2

PD code from KnotInfo, which gives -2 as the signature.

>>> T5 = Link([[1,5,2,4],[3,1,4,6],[5,3,6,2]])
>>> T5.signature()
-2

>>> T6 = Link('K3a1')
>>> T6.signature()
-2


Versus Sage
===========

>>> b = [1, 1, 1, 2, 1, 3, 1, 2, 2, 3, 3]
>>> K = Link(braid_closure=b)
>>> pd = K.PD_code(min_strand_index=1)
>>> L = Link(pd)
>>> K.signature(), L.signature()
(-6, -6)
>>> K.exterior().is_isometric_to(L.exterior())
True

>>> import sage.all
>>> B4 = sage.all.BraidGroup(4)
>>> w = B4(b)
>>> SK = sage.all.Link(w)
>>> SK.signature()
-6

>>> SB = Link(braid_closure=w)
>>> SB.signature()
-6

"""

import snappy
from spherogram import Crossing, Link, RationalTangle
import doctest
from sage.all import gap

print(doctest.testmod())
