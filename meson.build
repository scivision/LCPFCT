project('lcpfct', 'fortran',
  default_options : ['default_library=static', 'buildtype=release'])

fc = meson.get_compiler('fortran')

if fc.get_id() == 'gcc'
  add_project_arguments('-mtune=native', '-Werror=array-bounds', '-Wno-unused-label', '-Werror=line-truncation', language: 'fortran')
endif

lcpfct = library('lcpfct', 'src/gasdyn.f', 'src/lcpfct.f')

fast2d_exe = executable('fast2d', 'src/RunFast2d.f90', 'src/fast2d.f',
  link_with: lcpfct)
test('Shock2d', fast2d_exe)

shock_exe = executable('shock', 'src/RunShock.f90', 'src/shock.f',
  link_with: lcpfct)
test('Shock1d', shock_exe)

convect_exe = executable('convect', 'src/convect.f',
  link_with: lcpfct)
test('Convection', convect_exe)

diaphragm_exe = executable('diaphragm', 'src/diaphragm.f',
  link_with: lcpfct)
test('Diaphragm',  diaphragm_exe)
