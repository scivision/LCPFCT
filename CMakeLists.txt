cmake_minimum_required(VERSION 3.10)
project(lcpfct Fortran)

add_compile_options(-mtune=native -g -O3)

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)
  add_compile_options(-Werror=array-bounds -fexceptions
                    -Wall -Wextra -Wpedantic -Wno-unused-label -Wline-truncation)
  if(NOT WIN32)
    add_compile_options(-fstack-protector-all)
  endif()
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)

elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Flang)
  link_libraries(-static-flang-libs)
endif()

add_library(lcpfct gasdyn.f lcpfct.f)

add_executable(fast2d RunFast2d.f90 fast2d.f)
target_link_libraries(fast2d lcpfct)

add_executable(shock RunShock.f90 shock.f)
target_link_libraries(shock lcpfct)

add_executable(convect convect.f)
target_link_libraries(convect lcpfct)

add_executable(diaphragm diaphragm.f)
target_link_libraries(diaphragm lcpfct)
