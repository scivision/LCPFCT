implicit none

integer, parameter :: NPT = 202, &
                      NX = 50   ! Number of cells in the computational domain

real PYOUT(NPT*NX,6)

call shock(NX,PYOUT)

!write (*,*) PYOUT

end program
