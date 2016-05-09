#!/usr/bin/env python3
import subprocess
import setuptools #needed to enable develop
from os.path import join
from glob import glob

try:
    subprocess.run(['conda','install','--yes','--file','requirements.txt'])
except Exception as e:
    print('you will need to install packages in requirements.txt')
#%%
fortranfiles=glob('*.f')
#%%
with open('README.rst','r') as f:
	long_description = f.read()
#%%
from numpy.distutils.core import setup,Extension

ext=[Extension(name='lcpfct',
               sources=fortranfiles,
               f2py_options=['--quiet'],
               extra_f77_compile_args=['-Wno-unused-label']
)]
               #include_dirs=[root],
               #library_dirs=[root])]

#%% install
setup(name='pylcpfct',
      version='0.1',
	 description='Python wrapper for LCPFCT model',
	 long_description=long_description,
	 author='Michael Hirsch',
	 url='https://github.com/scienceopen/pylcpfct',
      packages=['pylcpfct'],
      ext_modules=ext,
	  install_requires=[],
      )
