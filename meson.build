project('lcpfct', 'fortran', default_options : ['default_library=static'])

fc = meson.get_compiler('fortran')
if fc.get_id() == 'gcc'
  add_global_arguments('-Wno-unused-label', '-Wline-truncation', language : 'fortran')
endif

lcpfct = library('lcpfct', 'src/gasdyn.f', 'src/lcpfct.f')

fast2d = executable('fast2d', 'src/RunFast2d.f90', 'src/fast2d.f', link_with : lcpfct)

shock = executable('shock', 'src/RunShock.f90', 'src/shock.f', link_with : lcpfct)

convect = executable('convect', 'src/convect.f', link_with : lcpfct)

diaphragm = executable('diaphragm', 'src/diaphragm.f', link_with : lcpfct)
