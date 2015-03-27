c Michael Hirsch 2015.

        Program runshock
        
        Implicit None
         
         Integer NPT, NX
         Parameter ( NPT = 202,  
     &                NX = 50  ) ! Number of cells in the computational domain
         Real PYOUT(NPT*NX,6)
     
        Call shock(NX,PYOUT)
        
c        Write (*,*) PYOUT
        
        Stop ' End of program SHOCK '
        End Program
