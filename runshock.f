c Michael Hirsch 2015.
c Examples of using NRL LCPFCT with Python, wrapping demo program SHOCK
c
c Step 1: compile other subroutines, not directly accessed from Python
c gfortran -c lcpfct.f gasdyn.f
c
c Step 2: 
c f2py3 -m shock -c shock.f


        Program runshock
        
        Implicit none
         
         Integer NPT, NX
         Parameter ( NPT = 202,  
     &                NX = 50  ) ! Number of cells in the computational domain
         Real RHON_PY(NPT*NX,6)

         Real      RHON(NPT),  RVXN(NPT),  RVTN(NPT),  ERGN(NPT), RELAX
         Real      RHO_IN,     PRE_IN,     VEL_IN,     GAMMA0
         Real      RHOAMB,     PREAMB,     VELAMB,     GAMMAM

         Common    / ARRAYS / RHON,   RVXN,   RVTN,   ERGN,   RELAX, 
     &                        RHO_IN, PRE_IN, VEL_IN, GAMMA0, 
     &                        RHOAMB, PREAMB, VELAMB, GAMMAM
     
        Call shock(NX,RHON_PY)
        
c        Write (*,*) RHON_PY
        
        Stop ' End of program RUNSHOCK '
        End
