from .graphs import *
from .presentations import *
from .links import *
from .codecs import *
# try:
#     import snappy
# except ImportError:
#     pass

# Make the module version number easily accessible.
from . import version as _version


def version():
    return _version.version


__version__ = version()


__all__ = ['ABC', 'ClosedBraid', 'CyclicList', 'CyclicWord',
           'DTcodec', 'Digraph',
           'DirectedEdge', 'DirectedMultiEdge', 'Edge', 'FatEdge',
           'FatGraph', 'Graph',
           'Link', 'MultiEdge', 'Poset', 'Presentation',
           'Crossing', 'Strand', 'WhiteheadMove',
           'Word', 'random_link',
           # from spherogram.links.tangles:
           'Tangle', 'CapTangle', 'CupTangle', 'RationalTangle',
           'ZeroTangle', 'InfinityTangle', 'MinusOneTangle', 'OneTangle', 'IntegerTangle',
           'IdentityBraid', 'ComponentTangle', 'join_strands']
