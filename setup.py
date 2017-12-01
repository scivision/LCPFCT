#!/usr/bin/env python
install_requires = ['pandas','numpy','matplotlib']
tests_require=['nose','coveralls']
#
from setuptools import find_packages
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
      packages=find_packages(),
	  description='Python wrapper for LCPFCT model',
	  author='Michael Hirsch, Ph.D.',
	  url='https://github.com/scivision/pylcpfct',
      ext_modules=ext,
      install_requires=install_requires,
      python_requires='>=2.7',
      tests_require=tests_require,
      extras_require={'tests':tests_require},
      )
