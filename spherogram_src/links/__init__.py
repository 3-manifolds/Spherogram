from .links import Crossing, Strand, Link
from .tangles import Tangle, RationalTangle, ZeroTangle, InfinityTangle, IdentityBraid, join_strands
from . import orthogonal 
import os, sys
Link.view = orthogonal.orthogonal_draw
from .random_links import random_knot

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

__all__ = ['Crossing', 'Strand', 'Link', 'Tangle', 'RationalTangle', 'ZeroTangle', 'InfinityTangle', 'IdentityBraid', 'join_strands', 'pdf_docs', 'random_knot']
