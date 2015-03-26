C=======================================================================

      Subroutine GASDYN ( K1, KN, BC1, BCN, DT )

C-----------------------------------------------------------------------
c
c  This routine integrates the gasdynamic equations using the momentum
c  component RVRN as the direction of integration and the momentum RVTN 
c  as the transverse direction.  In 2D models the two directions of 
c  integration are chosen by exchanging RVRN and RVTN in Common.

c  K1  . . .   Index of the integration's first cell       
c  KN  . . .   Index of the integration's last cell        
c  BC1 . . .   Indicates boundary condition on integration at K1
c  BCN . . .   Indicates boundary condition on integration at K1
c  DT  . . .   Timestep for the integrations of this step
c
C-----------------------------------------------------------------------

         Implicit  NONE
         Integer   NPT, K1, K1P, BC1, BCN, K, KN, KNP, IT
         Parameter ( NPT = 202 )
         Logical   PBC
         Real      SBC1,  SRV1,  SBCN,     SRVN,  VRHO1, VRHON
         Real      VRVR1, VRVRN, VRVT1,    VRVTN, VERG1, VERGN
         Real      MPINT(NPT), VEL(NPT),   UNIT(NPT),  ZERO(NPT)
         Real      RHOO(NPT),  RVRO(NPT),  RVTO(NPT),  ERGO(NPT)
         Real      VINT(NPT),  PRE(NPT),   MPVINT(NPT)
         Real      DTSUB,      DT,         RELAX
         Data      UNIT / NPT*1.0 /,       ZERO / NPT*0.0 /

         Real      RHO_IN,     PRE_IN,     VEL_IN,     GAMMA0
         Real      RHOAMB,     PREAMB,     VELAMB,     GAMMAM
         Real      RHON(NPT),  RVRN(NPT),  RVTN(NPT),  ERGN(NPT)
         Common    / ARRAYS / RHON,   RVRN,   RVTN,   ERGN,   RELAX,
     &                        RHO_IN, PRE_IN, VEL_IN, GAMMA0,
     &                        RHOAMB, PREAMB, VELAMB, GAMMAM

c  Prepare for the time integration. Index K is either I or J depending
c  on the definitions of RVRN and RVTN. Copies of the physical variable 
c  are needed to recover values for the whole step integration . . .
C-----------------------------------------------------------------------
      KNP = KN + 1
      K1P = K1 + 1
      PBC = .false.
      If ( BC1.eq.3 .OR. BCN.eq.3 ) PBC = .true.
      Do K = K1, KN
         RHOO(K) = RHON(K)
         RVRO(K) = RVRN(K)
         RVTO(K) = RVTN(K)
         ERGO(K) = ERGN(K)
      End do

c  Integrate first the half step then the whole step . . .
C-----------------------------------------------------------------------
      Do 500 IT = 1, 2
         DTSUB = 0.5*DT*FLOAT(IT)

         Do K = K1, KN
            VEL(K) = RVRN(K)/RHON(K) 
            PRE(K) = GAMMAM*(ERGN(K) 
     &             - 0.5*(RVRN(K)**2 + RVTN(K)**2)/RHON(K)) 
         End do
c  Calculate the interface velocities and pressures as weighted values
c  of the cell-centered values computed just above . . .
C-----------------------------------------------------------------------
         Do K = K1+1, KN
            VINT(K)  =  0.5*(VEL(K) + VEL(K-1))
            MPINT(K) = -0.5*(PRE(K) + PRE(K-1))
            MPVINT(K) = -0.5*(PRE(K)*vel(k) + PRE(K-1)*vel(k-1))
         End do

c  The unweighted interface averages can be computed as follows . . .


c  Call the FCT utility routines and set the boundary conditions. Other
c  boundary conditions could be added for inflow, outflow, etc . . .
c  BC1, BCN = 1  => ideal solid wall or axis boundary condition 
c  BC1, BCN = 2  => an extrapolative outflow boundary condition 
c  BC1, BCN = 3  => periodic boundary conditions . . .
c  BC1, BCN = 4  => specified boundary values (e.g. shock tube problem)
C-----------------------------------------------------------------------
         Go To ( 310, 320, 330, 340 ), BC1
  310       VINT(K1)   = 0.0
            MPINT(K1)  = - PRE(K1)
            MPVINT(K1) = 0.0
            Go To 350
  320       VINT(K1)   = VEL(K1)*(1.0 - RELAX) 
            MPINT(K1)  = - PRE(K1)*(1.0 - RELAX) - RELAX*PRE_IN
            MPVINT(K1) = MPINT(K1)*VINT(K1)
            Go To 350
  330       MPVINT(K1) = 1.0/( RHON(K1) + RHON(KN) )
            VINT(K1)   = (VEL(K1)*RHON(KN)+VEL(KN)*RHON(K1)) *MPVINT(K1)
            MPINT(K1)  = -(PRE(K1)*RHON(KN)+PRE(KN)*RHON(K1))*MPVINT(K1)
            MPVINT(K1) = -(PRE(K1)*VEL(K1)*RHON(KN) 
     &                + PRE(KN)*VEL(KN)*RHON(K1))*MPVINT(K1)
            Go To 350
  340       VINT(K1)   = VEL_IN
            MPINT(K1)  = - PRE_IN
            MPVINT(K1) = - PRE_IN*VEL_IN

  350    Go To ( 410, 420, 430, 440 ), BCN
  410       VINT(KNP)   = 0.0
            MPINT(KNP)  = - PRE(KN)
            MPVINT(KNP) = 0.0
            Go To 450
  420       VINT(KNP)   = VEL(KN)*(1.0 - RELAX)
            MPINT(KNP)  = - PRE(KN)*(1.0 - RELAX) - RELAX*PREAMB
            MPVINT(KNP) = MPINT(KNP)*VINT(KNP)
            Go To 450
  430       VINT(KNP)   = VINT(K1)
            MPINT(KNP)  = MPINT(K1)
            MPVINT(KNP) = MPVINT(K1)
            Go To 450
  440       VINT(KNP)   = VELAMB
            MPINT(KNP)  = - PREAMB
            MPVINT(KNP) = - PREAMB*VELAMB
  450    Continue

c  The velocity dependent FCT coefficients are set and the boundary 
c  condition calculations are completed.  Here the periodic boundary 
c  conditions require no action as (S)lope and (V)alue boundary value
c  specifiers are ignored in LCPFCT when PBC = .true.
C-----------------------------------------------------------------------
         Call VELOCITY ( VINT, K1, KNP, DTSUB )

         Go To ( 510, 520, 550, 540 ), BC1
  510       Call ZEROFLUX ( K1 )
            SBC1  = 1.0
            SRV1  = -1.0
            VRHO1 = 0.0
            VRVR1 = 0.0
            VRVT1 = 0.0
            VERG1 = 0.0
            Go To 550
  520       Call ZERODIFF ( K1 )
            SBC1  = 1.0 - RELAX
            SRV1  = 1.0 - RELAX
            VRHO1 = RELAX*RHO_IN
            VRVR1 = 0.0
            VRVT1 = 0.0
            VERG1 = RELAX*PRE_IN/GAMMAM
            Go To 550
  540       SBC1  = 0.0
            SRV1  = 0.0
            VRHO1 = RHO_IN
            VRVR1 = RHO_IN*VEL_IN
            VRVT1 = 0.0
            VERG1 = PRE_IN/GAMMAM + 0.5*RHO_IN*VEL_IN**2

  550    Go To ( 610, 620, 650, 640 ), BCN
  610       Call ZEROFLUX ( KNP )
            SBCN = 1.0
            SRVN = -1.0
            VRHON = 0.0
            VRVRN = 0.0
            VRVTN = 0.0
            VERGN = 0.0
            Go To 650
  620       Call ZERODIFF ( KNP )
            SBCN  = 1.0 - RELAX
            SRVN  = 1.0 - RELAX
            VRHON = RELAX*RHOAMB
            VRVRN = 0.0
            VRVTN = 0.0
            VERGN = RELAX*PREAMB/GAMMAM
            Go To 650
  640       SBCN  = 0.0
            SRVN  = 0.0
            VRHON = RHOAMB
            VRVRN = RHOAMB*VELAMB
            VRVTN = 0.0
            VERGN = PREAMB/GAMMAM + 0.5*RHOAMB*VELAMB**2
  650    Continue

c  Integrate the continuity equations using LCPFCT . . .  
C-----------------------------------------------------------------------
         Call LCPFCT ( RHOO, RHON, K1,KN, SBC1,VRHO1, SBCN,VRHON, PBC )

         Call SOURCES( K1,KN, DTSUB, 5, UNIT, MPINT, 
     &                                        MPINT(K1),  MPINT(KNP)  )    

         Call LCPFCT ( RVRO, RVRN, K1,KN, SRV1,VRVR1, SRVN,VRVRN, PBC )

         Call LCPFCT ( RVTO, RVTN, K1,KN, SBC1,VRVT1, SBCN,VRVTN, PBC )

         Call SOURCES( K1,KN, DTSUB, 4, UNIT, MPVINT, 
     &                                        MPVINT(K1), MPVINT(KNP) )    

         Call LCPFCT ( ERGO, ERGN, K1,KN, SBC1,VERG1, SBCN,VERGN, PBC )

  500 Continue       ! End of halfstep-wholestep loop.
      Return
      End

C=======================================================================
