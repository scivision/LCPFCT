cmake_minimum_required(VERSION 3.14...3.21)

project(lcpfct LANGUAGES Fortran)

enable_testing()


if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  add_compile_options(-mtune=native -Werror=array-bounds
                    -Wall -Wno-unused-label -Werror=line-truncation)
endif()

add_library(lcpfct src/gasdyn.f src/lcpfct.f)

add_executable(fast2d src/RunFast2d.f90 src/fast2d.f)
target_link_libraries(fast2d lcpfct)
add_test(NAME Shock2d COMMAND fast2d)

add_executable(shock src/RunShock.f90 src/shock.f)
target_link_libraries(shock lcpfct)
add_test(NAME Shock1d COMMAND shock)

add_executable(convect src/convect.f)
target_link_libraries(convect lcpfct)
add_test(NAME Convection COMMAND convect)

add_executable(diaphragm src/diaphragm.f)
target_link_libraries(diaphragm lcpfct)
add_test(NAME Diaphragm COMMAND diaphragm)
