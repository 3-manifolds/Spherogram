import sys, os, glob
from setuptools import setup, Command, Extension

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
        sources = ['planarity_src/planarity.c'] + planarity_sources, 
        include_dirs = [planarity_dir], 
        )

    ext_modules = [Planarity]


# The planarmap extension

pmap_dir = 'planarmap_src/'
pmap_src_dir = pmap_dir + '/src/'
pmap_src_files = [pmap_src_dir + file for file in
                  ['PMdef.c', 'PMplanmap.c', 'PMenlight.c',
                   'PMconjugation.c', 'PMextract.c', 'stats.c']]

Planarmap = Extension(
    name = 'spherogram.planarmap',

    sources =  [pmap_dir + 'planarmap.c'] + pmap_src_files, 
    include_dirs = [pmap_src_dir]
    )

ext_modules.append(Planarmap)


# If have Cython, check that .c files are up to date:

try:
    from Cython.Build import cythonize
    if 'clean' not in sys.argv:
        targets = ['planarity_src/planarity.pyx', pmap_dir + 'planarmap.pyx']
        targets = [file for file in targets if os.path.exists(file)]
        cythonize(targets)
except ImportError:
    pass 


# A real clean

class clean(Command):
    user_options = []
    def initialize_options(self):
        pass 
    def finalize_options(self):
        pass
    def run(self):
        os.system('rm -rf build dist')
        os.system('rm -rf spherogram*.egg-info')
        
# Main module

exec(open('spherogram_src/version.py').read())

# Get long description from README
long_description = open('README').read()
long_description = long_description.split('==\n')[1]
long_description = long_description.split('\nDeveloped')[0]

setup( name = 'spherogram',
       version = version,
       install_requires = ['networkx>=1.3', 'decorator'],
       dependency_links = [],
       packages = ['spherogram', 'spherogram.links',
                   'spherogram.links.test', 'spherogram.codecs',
                   'spherogram.dev', 'spherogram.dev.dev_jennet'],
       package_dir = {'spherogram' : 'spherogram_src', 'spherogram.dev':'dev'},
       package_data = {'spherogram.links'  :  ['doc.pdf']}, 
       ext_modules = ext_modules,
       cmdclass =  {'clean':clean},
       zip_safe = False,

       description= 'Spherical diagrams for 3-manifold topology', 
       long_description = long_description,
       author = 'Marc Culler and Nathan M. Dunfield',
       author_email = 'culler@uic.edu, nathan@dunfield.info',
       license='GPLv2+',
       url = 'https://bitbucket.org/t3m/spherogram',
       classifiers = [
           'Development Status :: 5 - Production/Stable',
           'Intended Audience :: Science/Research',
           'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
           'Operating System :: OS Independent',
           'Programming Language :: C',
           'Programming Language :: Python',
           'Topic :: Scientific/Engineering :: Mathematics',
        ],
        keywords = 'knot, link, SnapPy',
)

