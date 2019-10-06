implicit none

integer, parameter :: MAXTSTP = 801, NGRID = 64

real pyRho(NGRID,NGRID,MAXTSTP), pyVR(NGRID,NGRID,MAXTSTP), &
pyVZ(NGRID,NGRID,MAXTSTP), pyErg(NGRID,NGRID,MAXTSTP)

call fast2d(pyRho,pyVR,pyVZ,pyErg)

!write (*,*) PYOUT(50,1,:)

end program
