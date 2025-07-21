"""

Submodule for doing band moves to construct ribbon concordances and
seach for ribbon disks.  Code orignally developed by Nathan Dunfield
and Sherry Gong for their paper::

  Ribbon concordances and slice obstructions: experiments and examples
  https://arXiv.org/abs/FILLIN

"""

from . import merge_links
from . import core
from . import search
from .core import Band
from .search import verify_ribbon_to_unknot



