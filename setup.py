import os, shutil, sys, sysconfig, subprocess
from glob import glob
from distutils.util import get_platform
from setuptools import setup, Command, Extension


# Defensive linker flags for Linux:
if sys.platform.startswith('linux'):
    extra_link_args=['-Wl,-Bsymbolic-functions', '-Wl,-Bsymbolic']
else:
    extra_link_args=[]


# The planarity extension

try:
    import sage.libs
    ext_modules = []
except ImportError:    
    planarity_dir = 'planarity_src/c/'
    planarity_ui_sources = glob(planarity_dir + 'planarity*.c')
    planarity_sources = [file for file in glob('planarity_src/c/*.c')
                         if not file in planarity_ui_sources]

    if sys.platform.startswith('win'):
        extra_compile_args = ['-D_CRT_SECURE_NO_WARNINGS']
    else:
        extra_compile_args = []
    
    Planarity = Extension(
        name = 'spherogram.planarity',
        sources = ['planarity_src/planarity.c'] + planarity_sources, 
        include_dirs = [planarity_dir],
        extra_compile_args = extra_compile_args,
        extra_link_args = extra_link_args
        )

    ext_modules = [Planarity]

# A real clean

class SpherogramClean(Command):
    user_options = []
    def initialize_options(self):
        pass 
    def finalize_options(self):
        pass
    def run(self):
        junkdirs = (glob('build/lib*') +
                    glob('build/bdist*') +
                    glob('snappy*.egg-info') +
                    ['__pycache__', os.path.join('python', 'doc')]
        )
        for dir in junkdirs:
            try:
                shutil.rmtree(dir)
            except OSError:
                pass
        junkfiles = glob('*/*.so*') + glob('*/*.pyc') + glob('*/*.c') 
        for file in junkfiles:
            try:
                os.remove(file)
            except OSError:
                pass
        
class SpherogramTest(Command):
    user_options = []
    def initialize_options(self):
        pass 
    def finalize_options(self):
        pass
    def run(self):
        build_lib_dir = os.path.join(
            'build',
            'lib.{platform}-{version_info[0]}.{version_info[1]}'.format(
                platform=sysconfig.get_platform(),
                version_info=sys.version_info)
        )
        sys.path.insert(0, build_lib_dir)
        from spherogram.test import run_all_tests
        sys.exit(run_all_tests())

def check_call(args):
    try:
        subprocess.check_call(args)
    except subprocess.CalledProcessError:
        executable = args[0]
        command = [a for a in args if not a.startswith('-')][-1]
        raise RuntimeError(command + ' failed for ' + executable)
        
class SpherogramRelease(Command):
    user_options = [('install', 'i', 'install the release into each Python')]
    def initialize_options(self):
        self.install = False
    def finalize_options(self):
        pass
    def run(self):
        if os.path.exists('build'):
            shutil.rmtree('build')
        if os.path.exists('dist'):
            shutil.rmtree('dist')

        pythons = os.environ.get('RELEASE_PYTHONS', sys.executable).split(',')
        for python in pythons:
            check_call([python, 'setup.py', 'build'])
            check_call([python, 'setup.py', 'test'])
            if sys.platform.startswith('linux'):
                plat = get_platform().replace('linux', 'manylinux1')
                plat = plat.replace('-', '_')
                check_call([python, 'setup.py', 'bdist_wheel', '-p', plat])
                check_call([python, 'setup.py', 'bdist_egg'])
            else:
                check_call([python, 'setup.py', 'bdist_wheel'])

            if self.install:
                check_call([python, 'setup.py', 'install'])

        # Build sdist using the *first* specified Python
        check_call([pythons[0], 'setup.py', 'sdist'])

        # Double-check the Linux wheels
        if sys.platform.startswith('linux'):
            for name in os.listdir('dist'):
                if name.endswith('.whl'):
                    subprocess.check_call(['auditwheel', 'repair', os.path.join('dist', name)])

# The planarmap extension

pmap_dir = 'planarmap_src/'
pmap_src_dir = pmap_dir + '/src/'
pmap_src_files = [pmap_src_dir + file for file in
                  ['PMdef.c', 'PMplanmap.c', 'PMenlight.c',
                   'PMconjugation.c', 'PMextract.c', 'stats.c']]

Planarmap = Extension(
    name = 'spherogram.planarmap',
    sources =  [pmap_dir + 'planarmap.c'] + pmap_src_files, 
    include_dirs = [pmap_src_dir],
    extra_link_args = extra_link_args
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

# Main module

exec(open('spherogram_src/version.py').read())

# Get long description from README
long_description = open('README').read()
long_description = long_description.split('==\n')[1]
long_description = long_description.split('\nDeveloped')[0]

if sys.version_info >= (2, 7):
    networkx = 'networkx>=1.3'
else:
    networkx = 'networkx<1.10'

setup( name = 'spherogram',
       version = version,
       install_requires = [networkx, 'decorator', 'future'],
       dependency_links = [],
       packages = ['spherogram', 'spherogram.links',
                   'spherogram.links.test', 'spherogram.codecs',
                   'spherogram.dev', 'spherogram.dev.dev_jennet'],
       package_dir = {'spherogram' : 'spherogram_src', 'spherogram.dev':'dev'},
       package_data = {'spherogram.links'  :  ['doc.pdf']}, 
       ext_modules = ext_modules,
       cmdclass =  {'clean': SpherogramClean,
                    'test': SpherogramTest,
                    'release': SpherogramRelease,
       },
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

