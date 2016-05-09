cmake_minimum_required(VERSION 2.8.12)
project( lcpfct Fortran )

add_compile_options(-mtune=native -Werror=array-bounds -fexceptions -fstack-protector-all
                    -Wall -Wextra -Wpedantic -Wno-unused-label)

add_library(lcpfct gasdyn.f lcpfct.f)

add_executable(fast2d RunFast2d.f90 fast2d.f)
target_link_libraries(fast2d lcpfct)

add_executable(shock RunShock.f90 shock.f)
target_link_libraries(shock lcpfct)
