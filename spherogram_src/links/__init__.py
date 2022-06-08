import os
import sys

from .links import Crossing, Strand, Link, ClosedBraid
from .tangles import Tangle, CapTangle, CupTangle, RationalTangle, ZeroTangle, InfinityTangle, MinusOneTangle, OneTangle, IntegerTangle, IdentityBraid, ComponentTangle, join_strands
from . import orthogonal
from .random_links import random_link

Link.view = orthogonal.orthogonal_draw


def pdf_docs():
    "Open the PDF docs for links and tangles using the default application"
    file = __path__[0] + os.sep + "doc.pdf"
    if sys.platform.startswith('darwin'):
        command = 'open'
    elif sys.platform.startswith('win'):
        command = 'start'
    else:
        command = 'xdg-open'
    os.system(command + ' ' + file)


__all__ = ['Crossing', 'Strand', 'Link', 'ClosedBraid',
           'Tangle', 'CapTangle', 'CupTangle', 'RationalTangle',
           'ZeroTangle', 'InfinityTangle', 'MinusOneTangle', 'OneTangle', 'IntegerTangle',
           'IdentityBraid', 'join_strands',
           'pdf_docs', 'random_link']
