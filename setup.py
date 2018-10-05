#!/usr/bin/env python
import setuptools  # noqa: F401
from numpy.distutils.core import setup, Extension
from pathlib import Path
import os


if os.name == 'nt':
    sfn = Path(__file__).parent / 'setup.cfg'
    stxt = sfn.read_text()
    if '[build_ext]' not in stxt:
        with sfn.open('a') as f:
            f.write("[build_ext]\ncompiler = mingw32")

R = Path(__file__).parent / 'src'

ext = [Extension(name='lcpfctfort',
                 sources=[str(R/'lcpfct.f'), str(R/'gasdyn.f')],
                 extra_f77_compile_args=['-Wno-unused-label']
                 ),
       Extension(name='shock',
                 sources=[str(R/'shock.f'), str(R/'gasdyn.f'), str(R/'lcpfct.f')],
                 extra_f77_compile_args=['-Wno-unused-label']
                 ),
       Extension(name='fast2d',
                 sources=[str(R/'fast2d.f'), str(R/'gasdyn.f'), str(R/'lcpfct.f')],
                 extra_f77_compile_args=['-Wno-unused-label']
                 )
       ]

setup(ext_modules=ext)
