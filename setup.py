import sys, os, glob
from setuptools import setup, Command
from distutils.extension import Extension
from Cython.Distutils import build_ext

# The planarity extension

try:
    import sage.libs
    ext_modules = []
except ImportError:    
    planarity_dir = 'planarity_src/c/'
    planarity_ui_sources = glob.glob(planarity_dir + 'planarity*.c')
    planarity_sources = [file for file in glob.glob('planarity_src/c/*.c')
                         if not file in planarity_ui_sources]
    
    Planarity = Extension(
        name = 'spherogram.planarity',
        sources = ['planarity_src/planarity.pyx'] + planarity_sources, 
        include_dirs = [planarity_dir], 
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
        os.system('rm -f planarity_src/planarity.c')

# Main module

version='1.3'
with open('version.txt','w') as output:
    output.write(version)

setup( name = 'spherogram',
       version = version,
       zip_safe = False,
       install_requires = ['networkx>=1.3'],
       dependency_links = [],
       packages = ['spherogram', 'spherogram.links', 'spherogram.codecs', 'spherogram.dev', 'spherogram.dev.dev_jennet'],
       package_dir = {'spherogram' : 'spherogram_src', 'spherogram.dev':'dev'},
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
