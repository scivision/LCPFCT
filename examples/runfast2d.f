        Program RUNFAST2D
        
        Implicit None
        
        Integer,Parameter :: NPT = 202, MAXTSTP = 801, NGRID = 64
         
         Real PYOUT(NGRID,9,MAXTSTP)
         
         call fast2d(PYOUT)
         
c         write (*,*) PYOUT(50,1,:)
         
         stop ' end of fast2d '
         End Program 
