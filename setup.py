#!/usr/bin/env python
req = ['nose','pandas','numpy','matplotlib']
import pip
try:
    import conda.cli
    conda.cli.main('install',*req)
except Exception as e:
    pip.main(['install'] + req)
# %%
import setuptools #needed to enable develop

# %%
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
      packages=['pylcpfct'],
	  description='Python wrapper for LCPFCT model',
	  author='Michael Hirsch, Ph.D.',
	  url='https://github.com/scivision/pylcpfct',
      ext_modules=ext,
      )
