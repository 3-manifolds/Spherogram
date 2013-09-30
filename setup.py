import sys, os, glob
from setuptools import setup, Command
from distutils.extension import Extension
from Cython.Distutils import build_ext

# The planarity extension

try:
    import sage.libs
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

version='1.2'
with open('version.txt','w') as output:
    output.write(version)

setup( name = 'spherogram',
       version = version,
       zip_safe = False,
       install_requires = ['networkx>=1.3'],
       dependency_links = [],
       packages = ['spherogram', 'spherogram.links', 'spherogram.codecs'],
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
