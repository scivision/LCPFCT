C=======================================================================

          Subroutine FAST2D(PYRHO,PYVR,PYVZ,PYERG)

C-----------------------------------------------------------------------
c
c  BURSTING DIAPHRAGM "MUZZLE FLASH" - LCPFCT TEST # 4       August 1992
c
c  The problem begins with 1000:1 pressure and 100:1 density ratios 
c  across a diaphragm inside a solid cylindrical barrel.  The ideal wall
c  of the barrel is 10 cells thick (1.0 cm) with its inner radius given
c  as 1.5 cm and its outer radius of 2.5 cm.  The run starts at time 
c  t = 0.0 when the diaphragm at interface J = 11 (inside the barrel) is
c  ruptured.  The flow then expands upward in a 1D manner, spilling out 
c  of the barrel in a 2D flow which eventually reaches the boundaries at
c  R = 4.0 cm and Z = 4.0 cm where a very simple extrapolative outflow
c  condition is expressed through the LCPFCT boundary conditions values.
c  The outflow condition used here includes a slow relaxation to ambient 
c  conditions far from the origin.
c
C-----------------------------------------------------------------------

         Implicit  NONE
         
         Logical doplot
         Parameter ( doplot=.false.)
         

         Integer   NPT, I, J, IJ, NR, NZ, MAXSTP
         Parameter ( NPT = 202, NR = 64, NZ=64 , MAXSTP=801)

c         Real,   Intent(OUT)   ::   PYOUT(NR,9,MAXSTP) 
      Real,Intent(OUT) :: PYRHO(NR,NZ,MAXSTP),PYVR(NR,NZ,MAXSTP),
     &                    PYVZ(NR,NZ,MAXSTP),PYERG(NR,NZ,MAXSTP) 

         Integer   NRP, IALFR,     BC_AXIS, BC_WALL, BC_OUTF
         Integer   NZP, IALFZ,        IPRINT
         Integer   ICIN, ICOUT, JCTOP,      ISTEP
         Real      DR,         DZ,         DT,         TIME
         Real      COURANT,    DTNEW,      VTYPICAL,   RELAX
         Real      DTMAX,      VZMAX,      R(NPT),     Z(NPT)
         Real      RHO(NR,NZ), RVR(NR,NZ), RVZ(NR,NZ), ERG(NR,NZ)
         Real      DIN(9)

         Real      RHO_IN,     PRE_IN,     VEL_IN,     GAMMA0
         Real      RHOAMB,     PREAMB,     VELAMB,     GAMMAM
         Real      RHON(NPT),  RVRN(NPT),  RVTN(NPT),  ERGN(NPT)
         Common    / ARRAYS / RHON,   RVRN,   RVTN,   ERGN,   RELAX, 
     &                         RHO_IN, PRE_IN, VEL_IN, GAMMA0,
     &                         RHOAMB, PREAMB, VELAMB, GAMMAM

 1000 Format ('1', /, '   LCPFCT Test # 4 - FAST2D Barrel Explosion:',
     1        '  Step =', I4, /, 5X, I3, ' x', I3, ' Uniform Grid.',
     2        '  Time=', 1PE12.4, ' and   DT =', E12.4 )
 1005 Format (' After step ', I5, '  TIME = ', 1PE12.4,
     &           ' and timestep DT = ', E12.4 )
 1010 Format (1X, /, ' Fluid variables on selected lines ',
     1           /, ' I, J', '    RHO axis VZ     PRE ',
     2                      '     RHO top  VR     PRE ',
     3                      '     RHO wall VZ     PRE ', / )
 1011 Format ( I3, 1X, 3(1X, F8.2, F8.2, F8.2) ) 

c  The 2D barrel explosion program control parameters are specified.
c  (Change here to run other cases) . . .
C-----------------------------------------------------------------------
c      NR      =  64  ! Number of cells in the first (radial R) direction
      IALFR   =   2  ! Sets cylindrical coordinates in the R direction
      DR      =  0.1 ! Cell size (e.g., cm) in the radial direction
c      NZ      =  64  ! Number of cells in the second (axial Z) direction
      IALFZ   =   1  ! Sets Cartesian coordinates in the Z direction
      DZ      =  0.1 ! Cell size (e.g., cm) in the axial direction
      BC_AXIS =   1  ! Cylindrical axis set as an impermeable wall
      BC_OUTF =   2  ! Outer boundaries set as extrapolative outflow
      BC_WALL =   1  ! Walls of the barrel set as an ideal solid wall
c      MAXSTP  = 801  ! Maximum number of timesteps of length DT
      IPRINT  =  25  ! Initial frequency of validation printout results
      COURANT =  0.4 ! Approximate maximum Courant number allowed
      DTMAX = 2.0E-7 ! Maximum timestep allowed in the computation
      DT    = 1.0E-9 ! Initial (small guess) for starting timestep

c  Initialize the test problem geometry, a cylindrical shell JCTOP cells 
c  high in Z (indexed by J) which extends from the left of cell ICIN to 
c  the right of cell ICOUT in X (indexed by I).  
C-----------------------------------------------------------------------
      ICIN    =  16  ! Number of the innermost radial cell in the barrel
      ICOUT   =  25  ! Number of the outermost radial cell in the barrel
      JCTOP   =  20  ! Number of the uppermost axial cell in the barrel
      GAMMA0 = 1.4   ! Gas constant
      RHOAMB = 0.00129       ! Initialization and relaxation BCN = 2
      PREAMB = 1.013E+6      ! Initialization and relaxation BCN = 2
      VELAMB = 0.0           ! Initialization and relaxation BCN = 2
      RHO_IN = 100.0*RHOAMB  ! Initialization and relaxation BC1 = 2
      PRE_IN = 1000.0*PREAMB ! Initialization and relaxation BC1 = 2
      VEL_IN = 0.0           ! Initialization and relaxation BC1 = 2
      RELAX =  0.002 ! Relaxation rate, used when BC1 or BCN = 2
      GAMMAM = GAMMA0 - 1.0

c  Determine the cell interface locations, here a uniform grid . . .
C-----------------------------------------------------------------------
      NRP = NR + 1
      NZP = NZ + 1
      Do I = 1, NRP
       R(I) = DR*FLOAT(I-1)
      End do
      Do J = 1, NZP
       Z(J) = DZ*FLOAT(J-1)
      End do
c  Fill the arrays with air at STP and behind the diaphragm increase the
c  density by 100 to 1 and the pressure by 1000 to 1 . . .
C-----------------------------------------------------------------------
      Do J = 1, NZ
        Do  I = 1, NR
         RHO(I,J) = RHOAMB
         RVR(I,J) = 0.0
         RVZ(I,J) = 0.0
         ERG(I,J) = PREAMB/GAMMAM
        End do
      End do
      Do J = 1, JCTOP/2
        Do I = 1, ICIN - 1
         ERG(I,J) = PRE_IN/GAMMAM
         RHO(I,J) = RHO_IN
        End do
      End do
c  Mark the unused cells inside the cylindrical 'barrel' so they will
c  show up distinctly compared to ambient values in the plots.  This
c  has no effect as the simulation does not access these values . . .
C-----------------------------------------------------------------------
      Do J = 1, JCTOP
        Do I = ICIN, ICOUT
         ERG(I,J) = 20.0*ERG(I,J)
         RHO(I,J) = 20.0*RHO(I,J)
        End do
      End do
c  Begin loop over the timesteps . . .
C-----------------------------------------------------------------------
      TIME = 0.0
      Do 9999 ISTEP = 1, MAXSTP

c  Compute the next timestep based on a 'Courant' number COURANT . . .
C-----------------------------------------------------------------------
         VZMAX = 0.0
         Do J = 1, NZ
           Do I = 1, NR
            VTYPICAL = ERG(I,J)/RHO(I,J) 
            VZMAX = MAX ( VTYPICAL, VZMAX ) 
           End do
         End do
         VZMAX = SQRT ( VZMAX )
         DTNEW = COURANT*MIN(DR,DZ)/VZMAX
         DT = MIN ( DTMAX, 1.25*DT, DTNEW )

c  The results are printed when required . . .
C-----------------------------------------------------------------------
          If (doplot) Then
            Do IJ = 1, 64
               DIN(1) = RHO(1,IJ)/RHOAMB
               DIN(2) = 0.01*RVZ(1,IJ)/RHO(1,IJ)
               DIN(3) = (GAMMAM/PREAMB)*(ERG(1,IJ) - 0.5*
     &                  (RVR(1,IJ)**2 + RVZ(1,IJ)**2)/RHO(1,IJ))
               DIN(4) = RHO(IJ,64)/RHOAMB
               DIN(5) = 0.01*RVR(IJ,64)/RHO(IJ,64)
               DIN(6) = (GAMMAM/PREAMB)*(ERG(IJ,64) - 0.5*
     &                  (RVR(IJ,64)**2 + RVZ(IJ,64)**2)/RHO(IJ,64))
               DIN(7) = RHO(64,IJ)/RHOAMB
               DIN(8) = 0.01*RVZ(64,IJ)/RHO(64,IJ)
               DIN(9) = (GAMMAM/PREAMB)*(ERG(64,IJ) - 0.5*
     &                  (RVR(64,IJ)**2 + RVZ(64,IJ)**2)/RHO(64,IJ))
c                write (*,*) ISTEP
c                PYOUT(IJ,:,ISTEP) = DIN
            End do  
          End If

c  Integrate the fluid equations in the radial direction (indexed by I). 
c  The outer boundary condition at interface I = NR+1 is an extra- 
c  polation from the interior cell values with a slow relaxation to
c  the known distant ambient conditions . . .
C-----------------------------------------------------------------------
         Call RESIDIFF ( 0.999 )
         Call MAKEGRID ( R, R, 1, NRP, IALFR )
         Do J = 1, NZ

c  Pick up the data from the 2D arrays in the radial direction, setting
c  the temporary, compact 1D arrays for GASDYN . . .
C-----------------------------------------------------------------------
            Do I = 1, NR
               RHON(I) = RHO(I,J)
               RVRN(I) = RVR(I,J)
               RVTN(I) = RVZ(I,J)
               ERGN(I) = ERG(I,J)
            End do
c  Integrate along the radials inside and outside the cylinder  . . .
C-----------------------------------------------------------------------
            If ( J .le. JCTOP ) Then
               Call GASDYN ( 1, ICIN-1, BC_AXIS, BC_WALL, DT )
               Call GASDYN ( ICOUT+1, NR, BC_WALL, BC_OUTF, DT)

c  Integrate along the radials (indexing in I) above the cylinder
c  which reach from the axis to the outer boundary . ..
C-----------------------------------------------------------------------
            Else
               Call GASDYN ( 1, NR, BC_AXIS, BC_OUTF, DT )
            End If

c  Put the data back into the 2D arrays in the radial direction . . .    
C-----------------------------------------------------------------------
            Do I = 1, NR
               RHO(I,J) = RHON(I)
               RVR(I,J) = RVRN(I)
               RVZ(I,J) = RVTN(I)
               ERG(I,J) = ERGN(I)
            End do   
         End do     ! End loop integrating the NZ rows.

c  Integrate along the axials (indexing in J) which reach from the 
c  lower active J cell (1 or JCTOP+1) to the upper boundary.  The
c  upper boundary condition at interface J = NZ+1 (BCN = 2 ) is an 
c  extrapolation from the interior cell values with a slow relaxation 
c  to the known distant ambient conditions . . .
C-----------------------------------------------------------------------
         Call MAKEGRID ( Z, Z, 1, NZP, IALFZ )
         Do I = 1, NR

c  Pick up the data from the 2D arrays in the axial direction, setting
c  the temporary, compact 1D arrays for GASDYN . . .
C-----------------------------------------------------------------------
            Do J = 1, NZ
               RHON(J) = RHO(I,J)
               RVTN(J) = RVR(I,J)
               RVRN(J) = RVZ(I,J)
               ERGN(J) = ERG(I,J)
            End do
c  Integrate along the axials either from the lower solid boundary at 
c  interface J = 1 or from the top of the barrel at J = 21 for cells
c  with I = ICIN to ICOUT . . .
C-----------------------------------------------------------------------
            If ( I.ge.ICIN .and. I.le.ICOUT ) Then
               Call GASDYN ( JCTOP+1, NZ, BC_WALL, BC_OUTF, DT)
            Else
               Call GASDYN ( 1, NZ, BC_WALL, BC_OUTF, DT )
            End If 

c  Put the data back into the 2D arrays in the axial direction . . .    
C-----------------------------------------------------------------------
            Do J = 1, NZ
               RHO(I,J) = RHON(J)
               RVR(I,J) = RVTN(J)
               RVZ(I,J) = RVRN(J)
               ERG(I,J) = ERGN(J)
            End do   
         End do     ! End loop integrating the NR columns.

         TIME = TIME + DT
         
         PYRHO(:,:,ISTEP) = RHO
         PYVR(:,:,ISTEP)  = RVR
         PYVZ(:,:,ISTEP)  = RVZ
         PYERG(:,:,ISTEP) = ERG

 9999 End do        ! End of the timestep loop.

      Return
      End Subroutine FAST2D

C=======================================================================
