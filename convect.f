      Program CONVECT

C-----------------------------------------------------------------------
c
c  CONSTANT VELOCITY CONVECTION - LCPFCT Test # 1            August 1992
c
c  This program runs three periodic convection problems using LCPFCT and
c  the FCT utility routines .  The three profiles are the square wave, a
c  semicirle, and a Gaussian peak. The velocity is constant in space and
c  time and the grid is kept stationary.
c
C-----------------------------------------------------------------------
         use, intrinsic :: iso_fortran_env,only: stdout=>output_unit
         Implicit none

         Integer    NPT,           NX,            NXP
         Parameter  ( NPT = 202 )
         Logical    USE_LCP
         Integer    ISTEP,         JSTEP,         I
         Integer    MAXSTP,        IPRINT,        LOUT
         Real       DX, DT,        VELX,          TIME, XCELL
         Real       CSQUARE,       CCIRCLE,       CGAUSSP
         Real       ISQUARE,       ICIRCLE,       IGAUSSP
         Real       ESQUARE,       ECIRCLE,       EGAUSSP
         Real       ASQUARE(NPT),  ACIRCLE(NPT),  AGAUSSP(NPT)
         Real       SQUARE(NPT),   CIRCLE(NPT),   GAUSSP(NPT)
         Real       XINT(NPT),     VINT(NPT)

 1000 Format ( '1',/,'    LCPFCT Test #1 - Constant V Convection:',
     &         '   step =', I4, ' and TIME =', F7.3, /, 10X, 
     &         'with DX =', F6.3, ' DT =', F6.3, ' and VX =', F6.3, /)
 1001 Format ( '    I     X(I)     Square    exact     Circle ',
     &         '   exact   Gaussian    exact' ) 
 1002 Format ( I5, 7F10.5 )
 1003 Format ( 1X, /, ' Conserved sums', 6F10.5 )
 1004 Format (        ' Absolute error', F10.5, 10X, F10.5, 10X, F10.5 )

c  The Constant Velocity Convection control parameters are initialized.
c  (change here for other cases) . . .
C-----------------------------------------------------------------------
      USE_LCP = .true. ! Use the LCPFCT routine rather than CNVFCT
      NX      =   50   ! Number of cells in the periodic system
      DX      =  1.0   ! Cell size
      DT      =  0.2   ! Timestep for the calculation
      VELX    =  1.0   ! Constant X velocity, VELX*DT/DX = 0.2
      MAXSTP  =  501   ! Number of timesteps, two cycles of the system
      LOUT    = stdout   ! Logical unit number of printed output device
      IPRINT  =  125   ! Printout frequency, fluid moves 25 cells

c  The grid, velocity, and three density profiles are initialized . . .
C-----------------------------------------------------------------------
      NXP = NX + 1
      Do I = 1, NXP
         XINT(I) = FLOAT(I-1)*DX
         VINT(I) = VELX
      Enddo

      Call PROFILE ( 1, TIME, SQUARE, XINT, NX, NXP, VELX )
      Call PROFILE ( 2, TIME, CIRCLE, XINT, NX, NXP, VELX )
      Call PROFILE ( 3, TIME, GAUSSP, XINT, NX, NXP, VELX )

c  Set residual diffusion, grid, and velocity factors in LCPFCT . . .
C-----------------------------------------------------------------------
      Call RESIDIFF ( 1.0000 )
      Call MAKEGRID ( XINT, XINT, 1, NXP, 1 )
      Call VELOCITY ( VINT, 1, NXP, DT )

c  Begin loop over timesteps . . .
C-----------------------------------------------------------------------
      TIME = 0.0
      Do 9999 ISTEP = 1, MAXSTP

c  Results are printed as required . . .
C-----------------------------------------------------------------------
         If ( MOD(ISTEP-1,IPRINT) .eq. 0 ) Then
            JSTEP = ISTEP - 1
            Write ( LOUT, 1000 ) JSTEP, TIME, DX, DT, VELX
            Write ( LOUT, 1001 )
            Call PROFILE ( 1, TIME, ASQUARE, XINT, NX,NXP, VELX )
            Call PROFILE ( 2, TIME, ACIRCLE, XINT, NX,NXP, VELX )
            Call PROFILE ( 3, TIME, AGAUSSP, XINT, NX,NXP, VELX )
            ESQUARE = 0.0
            ECIRCLE = 0.0
            EGAUSSP = 0.0
            Do I = 1, NX
               XCELL = XINT(I) + 0.5*DX
               ESQUARE = ESQUARE + ABS(SQUARE(I) - ASQUARE(I))
               ECIRCLE = ECIRCLE + ABS(CIRCLE(I) - ACIRCLE(I))
               EGAUSSP = EGAUSSP + ABS(GAUSSP(I) - AGAUSSP(I))
               Write ( LOUT, 1002 )  I, XCELL, SQUARE(I), ASQUARE(I), 
     &            CIRCLE(I), ACIRCLE(I), GAUSSP(I), AGAUSSP(I)
            Enddo

            Call CONSERVE ( SQUARE, 1, NX, CSQUARE )
            Call CONSERVE ( CIRCLE, 1, NX, CCIRCLE )
            Call CONSERVE ( GAUSSP, 1, NX, CGAUSSP )
            If ( ISTEP .eq. 1 ) Then
               ISQUARE = CSQUARE
               ICIRCLE = CCIRCLE
               IGAUSSP = CGAUSSP
            End If
            Write ( LOUT, 1003) CSQUARE, ISQUARE, CCIRCLE, ICIRCLE,
     &                          CGAUSSP, IGAUSSP
            ESQUARE = ESQUARE/CSQUARE
            ECIRCLE = ECIRCLE/CCIRCLE
            EGAUSSP = EGAUSSP/CGAUSSP
            Write ( LOUT, 1004) ESQUARE, ECIRCLE, EGAUSSP
         End If

c    Advance the densities one timestep using LCPFCT and CNVFCT . . .
C-----------------------------------------------------------------------
         If ( USE_LCP ) Then
            Call LCPFCT ( SQUARE, SQUARE, 1, NX, 0.,0.,0.,0., .true.)
            Call LCPFCT ( CIRCLE, CIRCLE, 1, NX, 0.,0.,0.,0., .true.)
            Call LCPFCT ( GAUSSP, GAUSSP, 1, NX, 0.,0.,0.,0., .true.)
         Else
            Call CNVFCT ( SQUARE, SQUARE, 1, NX, 0.,0.,0.,0., .true.)
            Call CNVFCT ( CIRCLE, CIRCLE, 1, NX, 0.,0.,0.,0., .true.)
            Call CNVFCT ( GAUSSP, GAUSSP, 1, NX, 0.,0.,0.,0., .true.)
         End If
         TIME = TIME + DT

 9999 Continue           !  End of the timestep loop.

      End program

C=======================================================================

      Subroutine PROFILE ( TYPE, TIME, ARRAY, XINT, NX, NXP, VX )

C-----------------------------------------------------------------------
c
c  This subroutine computes three different analytic density profiles
c  depending on the value of TYPE . . .
c
c     TYPE = 1   Square wave profile between 0.0 and HEIGHT
c     TYPE = 2   Semicircular (elliptical) profile from 0.0 to HEIGHT
c     TYPE = 3   Gaussian peak profile between 0.0 and HEIGHT
c
c  The profiles are presented on a periodic domain NX cells long and a
c  crude integration is done within each cell to better approximate the
c  curved functions and to give an analytic approximation accounting for
c  convection across a partial cell.  
c
C-----------------------------------------------------------------------

      Implicit   NONE
      Integer    TYPE, NX, NXP, I, K
      Real       ARRAY(NX), XINT(NXP), TIME,  VX, SYSLEN, ARG
      Real       HEIGHT, X0, WIDTH,    XLEFT, XK, XRIGHT, XCENT
      Data       HEIGHT, X0, WIDTH / 1.0, 20.0, 10.0 /

      SYSLEN = XINT(NXP)
      Go To ( 100, 200, 300 ), TYPE

c  Compute the profile of the square wave . . .
C-----------------------------------------------------------------------
  100    XLEFT = (X0 - WIDTH) + VX*TIME
  101    If ( XLEFT .gt. SYSLEN ) Then
            XLEFT = XLEFT - SYSLEN 
            Go To 101
         End If
  102    If ( XLEFT .lt. 0.0 ) Then
            XLEFT = XLEFT + SYSLEN 
            Go To 102
         End If
         XRIGHT = XLEFT + 2.0*WIDTH

c  Loop over the cells in the numerical profile to be determined . . .
         Do 120 I = 1, NX
            ARRAY(I) = 0.0
            Do 110 K = 1, 10
               XK = XINT(I) + 0.1*(FLOAT(K)-0.5)*(XINT(I+1) - XINT(I))
               If ( XK .gt. XLEFT .and. XK .lt. XRIGHT ) Then
                  ARRAY(I) = ARRAY(I) + 0.1*HEIGHT
               Else 
                  XK = XK + SYSLEN 
                  If ( XK .gt. XLEFT .and. XK .lt. XRIGHT ) Then
                     ARRAY(I) = ARRAY(I) + 0.1*HEIGHT
                  End If
               End If
  110       Continue
  120    Continue
      Return

c  Compute the profile of the semicircle density hump . . .
C-----------------------------------------------------------------------
  200    XLEFT = (X0 - WIDTH) + VX*TIME
  201    If ( XLEFT .gt. SYSLEN ) Then
            XLEFT = XLEFT - SYSLEN 
            Go To 201
         End If
         XRIGHT = XLEFT + 2.0*WIDTH
  202    If ( XLEFT .lt. 0.0 ) Then
            XLEFT = XLEFT + SYSLEN 
            Go To 202
         End If
         XRIGHT = XLEFT + 2.0*WIDTH

c  Loop over the cells in the numerical profile to be determined . . .
         Do 220 I = 1, NX
            ARRAY(I) = 0.0
            Do 210 K = 1, 10
               XK = XINT(I) + 0.1*(FLOAT(K)-0.5)*(XINT(I+1) - XINT(I))
               If ( XK .gt. XLEFT .and. XK .lt. XRIGHT ) Then
                  XCENT = XLEFT + WIDTH
                  ARRAY(I) = ARRAY(I) + 0.1*HEIGHT*
     &               SQRT ( 1.0 - ((XK - XCENT)/WIDTH)**2 )
               Else 
                  XK = XK + SYSLEN 
                  If ( XK .gt. XLEFT .and. XK .lt. XRIGHT ) Then
                     XCENT = XLEFT + WIDTH
                     ARRAY(I) = ARRAY(I) + 0.1*HEIGHT*
     &                 SQRT ( 1.0 - ((XK - XCENT)/WIDTH)**2 )
                  End If
               End If
  210       Continue
  220    Continue
      Return

c  Compute the profile of the Gaussian density hump  . . .
C-----------------------------------------------------------------------
  300    XCENT = X0 + VX*TIME
  301    If ( XCENT .gt. SYSLEN ) Then
            XCENT = XCENT - SYSLEN 
            Go To 301
         End If
  302    If ( XCENT .lt. 0.0 ) Then
            XCENT = XCENT + SYSLEN 
            Go To 302
         End If

c  Loop over the cells in the numerical profile to be determined . . .
         Do 320 I = 1, NX
            ARRAY(I) = 0.0
            Do 310 K = 1, 10
               XK = XINT(I) + 0.1*(FLOAT(K)-0.5)*(XINT(I+1) - XINT(I))
               If ( XK .gt. (XCENT + 0.5*SYSLEN) ) XK = XK - SYSLEN
               If ( XK .lt. (XCENT - 0.5*SYSLEN) ) XK = XK + SYSLEN
               ARG = 4.0*((XK - XCENT)/WIDTH)**2
               ARRAY(I) = ARRAY(I) + 0.1*HEIGHT/EXP(AMIN1(30.0,ARG))
  310       Continue
  320    Continue
      End subroutine profile

