from distutils.core import setup
from distutils.extension import Extension
import os
import sys
import platform

openmm_dir = '@OPENMM_DIR@'
rigidbodyplugin_header_dir = '@RIGIDBODYPLUGIN_HEADER_DIR@'
rigidbodyplugin_library_dir = '@RIGIDBODYPLUGIN_LIBRARY_DIR@'

# setup extra compile and link arguments on Mac
extra_compile_args = []
extra_link_args = []

if platform.system() == 'Darwin':
    extra_compile_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7']
    extra_link_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7', '-Wl', '-rpath', openmm_dir+'/lib']

extension = Extension(name='_rigidbodyplugin',
                      sources=['RigidBodyPluginWrapper.cpp'],
                      libraries=['OpenMM', 'RigidBodyPlugin'],
                      include_dirs=[os.path.join(openmm_dir, 'include'), rigidbodyplugin_header_dir],
                      library_dirs=[os.path.join(openmm_dir, 'lib'), rigidbodyplugin_library_dir],
                      extra_compile_args=extra_compile_args,
                      extra_link_args=extra_link_args
                     )

setup(name='rigidbodyplugin',
      version='1.0',
      py_modules=['rigidbodyplugin'],
      ext_modules=[extension],
     )
