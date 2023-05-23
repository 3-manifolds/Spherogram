import os
import shutil
import sys
import sysconfig
import subprocess
from glob import glob
from setuptools import setup, Command, Extension


# Defensive linker flags for Linux:
if sys.platform.startswith('linux'):
    extra_link_args=['-Wl,-Bsymbolic-functions', '-Wl,-Bsymbolic']
else:
    extra_link_args=[]

# Install a custpom build environ for MSVC, if necessary.
if sys.platform == 'win32':
    try:
        # To customize your msvc environment create msvc_env.py
        # and define the dict dev_env to be a copy of the os.environ
        # that you have when you run python in the Command Prompt for
        # your installation of Visual Studio.
        from msvc_env import dev_env
        os.environ.update(dev_env)
        MSVC_extra_objects = []    
    except:
        pass

# The planarity extension

try:
    import sage.libs
    ext_modules = []
except ImportError:
    planarity_dir = 'planarity_src/c/'
    planarity_ui_sources = glob(planarity_dir + 'planarity*.c')
    planarity_sources = [file for file in glob('planarity_src/c/*.c')
                         if file not in planarity_ui_sources]

    if sys.platform.startswith('win'):
        extra_compile_args = [
            '-D_CRT_SECURE_NO_WARNINGS',
            r'-IC:\Program Files (x86)\Windows Kits\10\Include\10.0.17134.0\ucrt',
            r'-IC:\Program Files (x86)\Windows Kits\10\Include\10.0.17134.0\shared'
        ]
        extra_link_args = [
            r'/LIBPATH:C:\Program Files (x86)\Windows Kits\10\Lib\10.0.17134.0\um\x64',
            r'/LIBPATH:C:\Program Files (x86)\Windows Kits\10\Lib\10.0.17134.0\ucrt\x64',
        ]
    else:
        extra_compile_args = []
        extra_link_args = []

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
        sys.path.insert(0, os.path.abspath(build_lib_dir))
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
            if self.install:
                check_call([python, 'setup.py', 'pip_install'])
            else:
                check_call([python, 'setup.py', 'bdist_wheel'])

        # Build sdist using the *first* specified Python
        check_call([pythons[0], 'setup.py', 'sdist'])

        # Double-check the Linux wheels
        if sys.platform.startswith('linux'):
            for name in os.listdir('dist'):
                if name.endswith('.whl'):
                    subprocess.check_call(['auditwheel', 'repair', os.path.join('dist', name)])

class SpherogramPipInstall(Command):
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        python = sys.executable
        check_call([python, 'setup.py', 'bdist_wheel'])
        egginfo = 'spherogram.egg-info'
        if os.path.exists(egginfo):
            shutil.rmtree(egginfo)
        wheels = glob('dist' + os.sep + '*.whl')
        new_wheel = max(wheels, key=os.path.getmtime)
        check_call([python, '-m', 'pip', 'uninstall', '-y', 'spherogram'])
        check_call([python, '-m', 'pip', 'install', '--upgrade',
                    '--upgrade-strategy', 'only-if-needed',
                    new_wheel])

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
    needed_c_files = [pmap_dir + 'planarmap.c', 'planarity_src/planarity.c']
    if not all(os.path.exists(file) for file in needed_c_files):
        print("Error: The Python package Cython is missing so can't build Spherogram.")
        sys.exit(1)

# Main module

exec(open('spherogram_src/version.py').read())

# Get long description from README
long_description = open('README.rst').read()
long_description = long_description.split('==\n')[1]
long_description = long_description.split('\nDeveloped')[0]

install_requires = ['decorator', 'networkx',
                    'snappy_manifolds>=1.1.2', 'knot_floer_homology>=1.2']

setup( name = 'spherogram',
       version = version,
       install_requires = install_requires,
       python_requires = '>=3',
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
                    'pip_install':SpherogramPipInstall,
       },
       zip_safe = False,

       description= 'Spherical diagrams for 3-manifold topology',
       long_description = long_description,
       author = 'Marc Culler and Nathan M. Dunfield',
       author_email = 'culler@marc-culler.info, nathan@dunfield.info',
       license='GPLv2+',
       url = 'https://github.com/3-manifolds/Spherogram',
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
