      Program DIAPHRAGM

C-----------------------------------------------------------------------
c
c     1D BURSTING DIAPHRAGM PROBLEM - LCPFCT Test # 3        August 1992
c
c  This program runs a very simple 1D bursting diaphragm test problem 
c  using LCPFCT and its utility routines.  The fluid is ideal and
c  inviscid with constant GAMMA = 1.667.  The end walls are reflecting 
c  so the fluid should be totally contained in the domain.
c
C-----------------------------------------------------------------------
         use, intrinsic :: iso_fortran_env,only: stdout=>output_unit
         implicit none
       
         Integer   NPT, ALPHA, BC1, BCN,   MAXSTP,     IPRINT
         Parameter ( NPT = 202 )
         Integer   NX, NXP,    MX, ISTEP,  JSTEP, I,   LOUT
         Real      TIME,       DELTAX,     DELTAT
         REAL      XGRID(NPT), XNEXT(NPT), DXOFF,      DX_OSC 
         Real      RHOSUM,     RHVSUM,     PRESUM,     ERGSUM
         Real      VNEW(NPT),  PNEW(NPT),  TNEW(NPT)
         Real      VELX,       VXPAND,     SCALEG
         Real      ERG_IN,     ERGAMB,     RELAX

         Real      RHO_IN,     PRE_IN,     VEL_IN,     GAMMA0
         Real      RHOAMB,     PREAMB,     VELAMB,     GAMMAM
         Real      RHON(NPT),  RVXN(NPT),  RVTN(NPT),  ERGN(NPT)
         Common    / ARRAYS / RHON,   RVXN,   RVTN,   ERGN,   RELAX, 
     &                        RHO_IN, PRE_IN, VEL_IN, GAMMA0, 
     &                        RHOAMB, PREAMB, VELAMB, GAMMAM

 1000 Format ( '1', /, '    LCPFCT Test # 3 - Bursting Diaphragm:',
     &         '   Step =', I4, '   NX =', I3, '   DT =', F6.3, / )
 1001 Format ( 2X, I3, 6F12.5 )
 1002 Format ( '    I    Density  Temperature   Pressure  ',
     &    '  Velocity     Energy    Interfaces', / )
 1003 Format ( '0 Conservation Sums ', /, 5X, F12.5, 12X, 3F12.5, / )

c  The Bursting Diaphragm run control parameters are specified here.
c  (change for other cases) . . .
C-----------------------------------------------------------------------
      DELTAX = 1.0  ! Cell size 
      DELTAT = 0.05 ! Timestep for the calculation 
      ALPHA  =   1  ! (1 = Cartesian, 2 = Cylindrical, 3 = Spherical)
      LOUT   =  stdout  ! Logical unit number of printed output device
      NX     = 100  ! Number of cells in the computational domain
      MX     =  60  ! Number of cells initialized behind the shock
      MAXSTP = 1601 ! Maximum number of timesteps of length DELTAT
      IPRINT =  50  ! Frequency of intermediate result printouts

c  Initialize the variables in Common for use in GASDYN . . .
C-----------------------------------------------------------------------
      RHO_IN =  1.0    ! Initial mass density behind the diaphragm
      PRE_IN = 10.0    ! Initial (higher) pressure behind the diaphragm
      VEL_IN =  0.0    ! Initial velocity
      VXPAND =  3.1    ! Characteristic system expansion velocity 
      DX_OSC =  0.125  ! Amplitude of the grid jiggling after step 200
      GAMMA0 = 1.66667 ! Gas constant
      PREAMB =  1.0    ! Ambient (unshocked) pressure on the right
      RHOAMB =  1.0    ! Density of the unshocked fluid on the right
      VELAMB =  0.0    ! Initial velocity
      RELAX  =  0.002  ! Relaxation rate, used when BC1, BCN = 2
      GAMMAM = GAMMA0 - 1.0
      ERGAMB = PREAMB/(GAMMA0 - 1.0) + 0.5*RHOAMB*VELAMB**2
      ERG_IN = PRE_IN/(GAMMA0 - 1.0) + 0.5*RHO_IN*VEL_IN**2

c  Set up the fluid variables with the diaphragm at interface MX+1 . . .
C-----------------------------------------------------------------------
      Do I = MX+1, NX
         RHON(I) = RHOAMB
         RVXN(I) = RHOAMB*VELAMB
         RVTN(I) = 0.0
         ERGN(I) = ERGAMB
      Enddo
      Do I = 1, MX
         RHON(I) = RHO_IN
         RVXN(I) = RHO_IN*VEL_IN
         RVTN(I) = 0.0
         ERGN(I) = ERG_IN
      Enddo
         
c  Begin loop over timesteps . . .
C-----------------------------------------------------------------------
      BC1 = 1
      BCN = 1
      Call RESIDIFF ( 0.998 )

      TIME = 0.0
      Do 9999 ISTEP = 1, MAXSTP

c  Define the cell interface locations and physical variables.  The grid
c  is expanded at the rate VELX at I = NXP after step 201 (as the shock 
c  approaches the boundary) by making XNEXT at the end of the timestep
c  proportionately larger than XGRID at the beginning of the step.  The 
c  system length is then renormalized to its original size, capturing 
c  the similarity solution by equating the grid expansion to the shock 
c  velocity.  A small jiggle is added to this systematic expansion to
c  show the added flexibility of the continuity solver.
C-----------------------------------------------------------------------
         NXP = NX + 1
         If ( ISTEP .gt. 200 ) Then
            VELX = VXPAND
            DXOFF = DX_OSC
            If ( MOD(ISTEP,2) .eq. 0 ) DXOFF = - DX_OSC
         Else
            VELX = 0.0
            DXOFF = 0.0
         End If

         Do I = 1, NXP
            XGRID(I)  = FLOAT(I-MX-1)*DELTAX
            XNEXT(I) = XGRID(I)  
         Enddo
         SCALEG = (VELX*DELTAT + XGRID(NXP))/XGRID(NXP)
         Do I = 3, NX-2
            XNEXT(I)   = (XNEXT(I) + DXOFF)*SCALEG
            XGRID(I)   = XGRID(I) - DXOFF
         Enddo
         Call MAKEGRID ( XGRID, XNEXT, 1, NXP, ALPHA )

c  The results are printed when required . . .
C-----------------------------------------------------------------------
         If ( MOD(ISTEP-1, IPRINT) .eq. 0) Then
            JSTEP = ISTEP - 1
            IPRINT = 2*(ISTEP - 1)
            If ( ISTEP .lt. 200 ) IPRINT = 50 
            Write ( LOUT, 1000 ) JSTEP, NX, DELTAT
            Do I = 1, NX
               VNEW(I) = RVXN(I)/RHON(I)
               PNEW(I) = GAMMAM*(ERGN(I) - 0.5*RVXN(I)*VNEW(I))
               TNEW(I) = PNEW(I)/RHON(I)
            Enddo
            Write (  LOUT, 1002 )
            Write (  LOUT, 1001 ) ( I, RHON(I), TNEW(I), PNEW(I), 
     &           VNEW(I), ERGN(I), XGRID(I), I = 1, NX )
            Call CONSERVE (RHON, 1, NX, RHOSUM)
            Call CONSERVE (PNEW, 1, NX, PRESUM)
            Call CONSERVE (RVXN, 1, NX, RHVSUM)
            Call CONSERVE (ERGN, 1, NX, ERGSUM)
            Write (  LOUT, 1003 ) RHOSUM, PRESUM, RHVSUM, ERGSUM
         End If

c  The FCT integration of the continuity equations is performed . . .
C-----------------------------------------------------------------------
         Call GASDYN ( 1, NX, BC1, BCN, DELTAT )
         TIME = TIME + DELTAT

 9999 Continue    ! End of the loop over timesteps.

      End Program
