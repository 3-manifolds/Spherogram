from .presentations import *
from .links import *
from .codecs import *
from .links.reshetikhin_turaev import DictLaurentPolynomial, RMatrix, colored_links_gould_R_matrices, colored_jones_R_matrices, prefactor_colored_jones

# Make the module version number easily accessible.
from . import version as _version


def version():
    return _version.version


__version__ = version()


__all__ = ['ABC', 'ClosedBraid', 'CyclicList', 'CyclicWord',
           'DTcodec', 'Digraph',
           'DirectedEdge', 'DirectedMultiEdge', 'Edge', 'FatEdge',
           'FatGraph', 'Graph',
           'Link', 'MultiEdge', 'Presentation',
           'Crossing', 'Strand', 'WhiteheadMove',
           'Word', 'random_link',
           # from spherogram.links.tangles:
           'Tangle', 'CapTangle', 'CupTangle', 'RationalTangle',
           'ZeroTangle', 'InfinityTangle', 'MinusOneTangle', 'OneTangle', 'IntegerTangle',
           'IdentityBraid', 'BraidTangle', 'ComponentTangle', 'join_strands',
           'DictLaurentPolynomial', 'RMatrix', 'colored_links_gould_R_matrices', 'colored_jones_R_matrices', 'prefactor_colored_jones']
