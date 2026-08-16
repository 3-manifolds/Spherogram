"""

Submodule for doing band moves to construct ribbon concordances and
search for ribbon disks.  Code originally developed by Nathan Dunfield
and Sherry Gong for their paper::

  Ribbon concordances and slice obstructions: experiments and examples
  https://arXiv.org/abs/2512.21825

"""

from . import merge_links
from . import core
from . import search
from .core import Band, normalize_crossing_labels
from .search import verify_ribbon_to_unknot



