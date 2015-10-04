#!/usr/bin/env python3
import setuptools #needed to enable develop
from numpy.distutils.core import setup,Extension
from os.path import join
from glob import glob

fortranfiles=glob('*.f')

root='./'

fortranpaths = [join(root,f) for f in fortranfiles]
#%%
with open('README.rst','r') as f:
	long_description = f.read()

ext=[Extension(name='lcpfct',
               sources=fortranpaths,
               f2py_options=['--quiet'],
#               extra_f77_compile_args=['-finit-local-zero']
)]
               #include_dirs=[root],
               #library_dirs=[root])]

#%% install
setup(name='lcpfct',
      version='0.1',
	 description='Python wrapper for LCPFCT model',
	 long_description=long_description,
	 author='Michael Hirsch',
	 url='https://github.com/scienceopen/lcpfct',
      packages=['lcpfct'],
      ext_modules=ext,
	  install_requires=['pandas','matplotlib'],
      )
