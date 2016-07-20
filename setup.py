#!/usr/bin/env python
import subprocess
import setuptools #needed to enable develop

try:
    subprocess.call(['conda','install','--file','requirements.txt'])
except Exception as e:
    pass
#%%
from numpy.distutils.core import setup,Extension

ext=[Extension(name='lcpfct',
               sources=['lcpfct.f','gasdyn.f'],
               f2py_options=['--quiet'],
               extra_f77_compile_args=['-Wno-unused-label']
               ),
Extension(name='shock',
               sources=['shock.f','gasdyn.f','lcpfct.f'],
               f2py_options=['--quiet'],
               extra_f77_compile_args=['-Wno-unused-label']
               ),
Extension(name='fast2d',
               sources=['fast2d.f','gasdyn.f','lcpfct.f'],
               f2py_options=['--quiet'],
               extra_f77_compile_args=['-Wno-unused-label']
               )
    ]

               #include_dirs=[root],
               #library_dirs=[root])]

#%% install
setup(name='pylcpfct',
	 description='Python wrapper for LCPFCT model',
	 author='Michael Hirsch',
	 url='https://github.com/scienceopen/pylcpfct',
      ext_modules=ext,
	  install_requires=[],
    packages=['pylcpfct'],
      )
