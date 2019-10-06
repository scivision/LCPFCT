cmake_minimum_required(VERSION 3.1)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Debug or Release")
endif()
project(lcpfct Fortran)
enable_testing()


if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  add_compile_options(-mtune=native -Werror=array-bounds
                    -Wall -Wextra -Wpedantic -Wno-unused-label -Wline-truncation)

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Flang)

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
