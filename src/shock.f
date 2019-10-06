C=======================================================================

      Subroutine SHOCK(NX,PYOUT)
C-----------------------------------------------------------------------
c
c     PROGRESSING 1D GASDYNAMIC SHOCK - LCPFCT Test # 2      August 1992
c
c     This program runs a simple 1D gasdynamic shock through a uniform
c     grid using LCPFCT and its utility routines. The fluid is ideal and
c     inviscid with constant GAMMA0 = 1.4.  The boundary conditions are
c     specified external values on both ends of the system.
c
C-----------------------------------------------------------------------

         Implicit NONE

         Logical,parameter ::doplot=.false.
         Integer   ALPHA, BC1, BCN,   MAXSTP, IPRINT
         Integer, Intent(IN) ::   NX
         Integer,Parameter :: NPT = 202
         Real, Intent(OUT)   ::   PYOUT(NPT*NX,6)

         Integer   NXP,    MX, ISTEP,  JSTEP, I,   LOUT
         Real      MACH,       V0,         DELTAX,     DELTAT
         Real      RHOSUM,     RHVSUM,     PRESUM,     ERGSUM
         Real      VNEW(NPT),  PNEW(NPT),  TNEW(NPT),  XINT(NPT)
         Real      CSAMB,      ERG_IN,     ERGAMB,     RELAX
         Real      TIME

         Real      RHO_IN,     PRE_IN,     VEL_IN,     GAMMA0
         Real      RHOAMB,     PREAMB,     VELAMB,     GAMMAM
         Real      RHON(NPT),  RVXN(NPT),  RVTN(NPT),  ERGN(NPT)
         Common    / ARRAYS / RHON,   RVXN,   RVTN,   ERGN,   RELAX,
     &                        RHO_IN, PRE_IN, VEL_IN, GAMMA0,
     &                        RHOAMB, PREAMB, VELAMB, GAMMAM

 1000 Format ( '1    LCPFCT Test # 2 - Progressing Shock:',
     &         '   Step =', I4, '   NX =', I3, '   DT =', F6.3 )
 1001 Format (2X, I3, 1P6E13.5)
 1002 Format ('    I    Density  Temperature   Pressure  ',
     &    '  Velocity     Energy    Interfaces')
 1003 Format ('0 Conservation Sums ', /, 5X, 1PE12.4, 12X, 3E12.4, /)

c  The Progressing Shock test run control parameters are specified here.
c  (change for other cases) . . .
C-----------------------------------------------------------------------
      MACH   = 5.0  ! Mach number of the incoming ambient flow
      V0     = 3.0  ! Shock speed in the lab frame
      DELTAX = 1.0  ! Cell size
      DELTAT = 0.05 ! Timestep for the calculation
      ALPHA  = 1    ! (1 = Cartesian, 2 = Cylindrical, 3 = Spherical)
      LOUT   =  11  ! Logical unit number of printed output device
      MX     =  10  ! Number of cells initialized behind the shock
      MAXSTP = 201  ! Maximum number of timesteps of length DELTAT
      IPRINT =  25  ! Frequency of intermediate result printouts

c  Initialize the variables in Common for use in GASDYN . . .
C-----------------------------------------------------------------------
      GAMMA0 = 1.4  ! Gas constant
      PREAMB = 1.0  ! Ambient (unshocked) pressure on the right
      RHOAMB = 1.0  ! Density of the unshocked fluid on the right
      RELAX  = 0.0  ! Relaxation rate, not used when BC1, BCN = 2
      GAMMAM = GAMMA0 - 1.0

c  The Rankine-Hugoniot conditions are set for boundaries . . .
C-----------------------------------------------------------------------
      CSAMB  = SQRT (GAMMA0*PREAMB/RHOAMB)
      VELAMB = -MACH*CSAMB
      VEL_IN = VELAMB*(GAMMAM + 2.0/MACH**2)/(GAMMA0 + 1.0)
      RHO_IN = RHOAMB*VELAMB/VEL_IN
      PRE_IN = PREAMB - RHO_IN*VEL_IN**2 + RHOAMB*VELAMB**2
      VELAMB = VELAMB + V0
      VEL_IN = VEL_IN + V0
      ERGAMB = PREAMB/(GAMMA0 - 1.0) + 0.5*RHOAMB*VELAMB**2
      ERG_IN = PRE_IN/(GAMMA0 - 1.0) + 0.5*RHO_IN*VEL_IN**2

c  Define the cell interface locations and physical variables . . .
C-----------------------------------------------------------------------
      NXP = NX + 1
      Do I = 1, NXP
          XINT(I) = FLOAT(I-1)*DELTAX
      End Do
      Do I = MX+1, NX
         RHON(I) = RHOAMB
         RVXN(I) = RHOAMB*VELAMB
         RVTN(I) = 0.0
         ERGN(I) = ERGAMB
      End Do
      Do I = 1, MX
         RHON(I) = RHO_IN
         RVXN(I) = RHO_IN*VEL_IN
         RVTN(I) = 0.0
         ERGN(I) = ERG_IN
      End Do
c one time write of column headers
      If (doplot) Then
          Write (  LOUT, 1002 )
      End if
c  Begin loop over timesteps . . .
C-----------------------------------------------------------------------
      BC1 = 4
      BCN = 4
      Call RESIDIFF ( 1.000 )
      Call MAKEGRID ( XINT, XINT, 1, NXP, ALPHA )

      Do ISTEP = 1, MAXSTP

c  The results are printed when required . . .
C-----------------------------------------------------------------------
c         If ( MOD(ISTEP-1, IPRINT) .eq. 0) Then

            JSTEP = ISTEP - 1
            Do I = 1, NX
               VNEW(I) = RVXN(I)/RHON(I)
               PNEW(I) = GAMMAM*(ERGN(I) - 0.5*RVXN(I)*VNEW(I))
               TNEW(I) = PNEW(I)/RHON(I)
            End do

            Call CONSERVE (RHON, 1, NX, RHOSUM)
            Call CONSERVE (PNEW, 1, NX, PRESUM)
            Call CONSERVE (RVXN, 1, NX, RHVSUM)
            Call CONSERVE (ERGN, 1, NX, ERGSUM)
            PRESUM = PRESUM/GAMMAM
           If (doplot) Then
            Write ( LOUT, 1000 ) JSTEP, NX, DELTAT
            Write (  LOUT, 1001 ) ( I, RHON(I), TNEW(I), PNEW(I),
     &           VNEW(I), ERGN(I), XINT(I), I = 1, NX )
            Write (  LOUT, 1003 ) RHOSUM, PRESUM, RHVSUM, ERGSUM
           End If
c         End If

c  The FCT integration of the continuity equations is performed . . .
C-----------------------------------------------------------------------
         Call GASDYN ( 1, NX, BC1, BCN, DELTAT )
         TIME = TIME + DELTAT
c let's pass an array out to Python!
         PYOUT((NX*(ISTEP-1)+1):(NX*ISTEP),1) = RHON(1:NX)
         PYOUT((NX*(ISTEP-1)+1):(NX*ISTEP),2) = TNEW(1:NX)
         PYOUT((NX*(ISTEP-1)+1):(NX*ISTEP),3) = PNEW(1:NX)
         PYOUT((NX*(ISTEP-1)+1):(NX*ISTEP),4) = VNEW(1:NX)
         PYOUT((NX*(ISTEP-1)+1):(NX*ISTEP),5) = ERGN(1:NX)
         PYOUT((NX*(ISTEP-1)+1):(NX*ISTEP),6) = XINT(1:NX)
      End do

      End Subroutine SHOCK

