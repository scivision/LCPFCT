cmake_minimum_required(VERSION 3.7)
project(lcpfct Fortran)

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  add_compile_options(-g -O0)
else()
  add_compile_options(-O3)
endif()

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
  add_compile_options(-mtune=native -Werror=array-bounds -fexceptions
                    -Wall -Wextra -Wpedantic -Wno-unused-label -Wline-truncation)
  if(UNIX AND NOT APPLE)
    add_compile_options(-fstack-protector-all)
  endif()
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)

elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Flang)
  link_libraries(-static-flang-libs)
endif()

add_library(lcpfct src/gasdyn.f src/lcpfct.f)

add_executable(fast2d src/RunFast2d.f90 src/fast2d.f)
target_link_libraries(fast2d lcpfct)

add_executable(shock src/RunShock.f90 src/shock.f)
target_link_libraries(shock lcpfct)

add_executable(convect src/convect.f)
target_link_libraries(convect lcpfct)

add_executable(diaphragm src/diaphragm.f)
target_link_libraries(diaphragm lcpfct)
