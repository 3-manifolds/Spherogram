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

try:
    import sage.all
    ext_modules = []
except ImportError:    
    planarity_dir = ['planarity_src/planarity-read-only/c']
    planarity_extra_objects = glob.glob('planarity_src/planarity-read-only/c/*.o')
    if not os.path.exists(planarity_dir[0]) and 'clean' not in sys.argv:
        os.chdir('planarity_src')
        os.system('sh build_planarity.sh')
        os.chdir('..')
        planarity_extra_objects = glob.glob('planarity_src/planarity-read-only/c/*.o')
        if len(planarity_extra_objects) == 0:
            print("NOTE: Need to run 'build_planarity.sh' script in 'planarity_src' before this module can be built.")
            sys.exit()

    Planarity = Extension(
        name = 'spherogram.planarity',
        sources = ['planarity_src/planarity.pyx'], 
        include_dirs = planarity_dir, 
        extra_objects = planarity_extra_objects,
        )

    ext_modules = [Planarity]


class clean(Command):
    user_options = []
    def initialize_options(self):
        pass 
    def finalize_options(self):
        pass
    def run(self):
        os.system('rm -rf build dist')
        os.system('rm -rf spherogram*.egg-info')
        os.system('rm -rf planarity_src/planarity-read-only')
        os.system('rm -f planarity_src/planarity.c')


# Main module 

setup( name = 'spherogram',
       zip_safe = False,
       install_requires = [],
       dependency_links = [],
       packages = ['spherogram', 'spherogram.links'],
       package_dir = {'spherogram' : 'spherogram_src'},
       package_data = {'spherogram.links'  :  ['doc.pdf']}, 
       ext_modules = ext_modules,
       cmdclass =  {'build_ext': build_ext, 'clean':clean},
       entry_points = {},
       author = 'Marc Culler and Nathan Dunfield and John Berge',
       author_email = 'culler@math.uic.edu, nmd@illinois.edu',
       description = 'Spherical diagrams for 3-manifold topology',
       license = 'GPL',
       keywords = 'graphs, presentations',
       url = '',
       )

