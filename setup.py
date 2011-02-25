import setuptools
import pkg_resources
import sys, os, glob
sys.path.remove(os.path.realpath(os.curdir))

# Hack to patch setuptools so that it treats Cython
# as a replacement for pyrex.

from distutils.core import Extension as _Extension
from setuptools.dist import _get_unpatched
_Extension = _get_unpatched(_Extension)

try:
    from Cython.Distutils import build_ext
except ImportError:
    have_cython = False
else:
    have_cython = True


class Extension(_Extension):
    """
    This modified version of setuptools Extension allows us
    to use Cython instead of pyrex.  If Cython is not installed
    on this system, it will assume that a Cython-generated .c
    file is present in the distribution.
    """
    if not have_cython:
        # convert .pyx extensions to .c
        def __init__(self,*args,**kw):
            _Extension.__init__(self,*args,**kw)
            sources = []
            for s in self.sources:
                if s.endswith('.pyx'):
                    sources.append(s[:-3]+'c')
                else:
                    sources.append(s)
            self.sources = sources

class Library(Extension):
    """Just like a regular Extension, but built as a library instead"""

import sys, distutils.core, distutils.extension
distutils.core.Extension = Extension
distutils.extension.Extension = Extension
if 'distutils.command.build_ext' in sys.modules:
    sys.modules['distutils.command.build_ext'].Extension = Extension
# End of hack

from setuptools import setup, Command
from pkg_resources import load_entry_point

# The planarity extension

planarity_dir = ['planarity-read-only/c']
planarity_extra_objects = glob.glob('planarity-read-only/c/*.o')

Planarity = Extension(
    name = "spherogram.planarity",
    sources = ["planarity.pyx"], 
    include_dirs = planarity_dir, 
    extra_objects = planarity_extra_objects,
)

setup( name = "spherogram",
#       version = version,
       zip_safe = False,
       install_requires = [],
       dependency_links = [],
       packages = ["spherogram"],
       ext_modules = [Planarity],
       cmdclass =  {'build_ext': build_ext},
       entry_points = {},
       author = "Marc Culler and Nathan Dunfield and John Berge",
       author_email = "culler@math.uic.edu, nmd@illinois.edu",
       description = "Spherical diagrams for 3-manifold topology",
       license = "GPL",
       keywords = "graphs, presentations",
       url = "",
       )

