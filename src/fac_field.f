      SUBROUTINE fac_field (A,theta0,dTheta,mode,isEven,isSmoothed,
     .     dPhi,xKappa,
     .     tilt,X,Y,Z,
     .     BX,BY,BZ)            !   SEE NB# 6, P.60
C
C  CALCULATES COMPONENTS  OF REGION 1/2 FIELD IN SPHERICAL COORDS.  DERIVED FROM THE S/R DIPDEF2C (WHICH
C    DOES THE SAME JOB, BUT INPUT/OUTPUT THERE WAS IN SPHERICAL COORDS, WHILE HERE WE USE CARTESIAN ONES)
C
C   INPUT:  NUMB=1 (2) FOR REGION 1 (2) CURRENTS
C           MODE=1 YIELDS SIMPLE SINUSOIDAL MLT VARIATION, WITH MAXIMUM CURRENT AT DAWN/DUSK MERIDIAN
C     WHILE MODE=2 YIELDS THE SECOND HARMONIC.
C
C
      IMPLICIT NONE
      
C     Inputs
      REAL*8  A(30)
      REAL*8  theta0
      REAL*8  dTheta
      INTEGER mode
      LOGICAL isEven
      LOGICAL isSmoothed
      REAL*8  dPhi
      REAL*8  xKappa
      REAL*8  tilt
      REAL*8  X,Y,Z
      
C     Outputs
      REAL*8  BX,BY,BZ      
      
c     Internal variables      
      INTEGER M
      REAL*8  XS,YS,ZS
      REAL*8  Xsc,Ysc,Zsc
      REAL*8  BXS,BYAS,BZS
      
C  (1) DPHI:   HALF-DIFFERENCE (IN RADIANS) BETWEEN DAY AND NIGHT LATITUDE OF FAC OVAL AT IONOSPHERIC ALTITUDE;
C              TYPICAL VALUE: 0.06
C  (2) B:      AN ASYMMETRY FACTOR AT HIGH-ALTITUDES;  FOR B=0, THE ONLY ASYMMETRY IS THAT FROM DPHI
C              TYPICAL VALUES: 0.35-0.70
C  (3) RHO_0:  A FIXED PARAMETER, DEFINING THE DISTANCE RHO, AT WHICH THE LATITUDE SHIFT GRADUALLY SATURATES AND
C              STOPS INCREASING
C              ITS VALUE WAS ASSUMED FIXED, EQUAL TO 7.0.
C  (4) XKAPPA: AN OVERALL SCALING FACTOR, WHICH CAN BE USED FOR CHANGING THE SIZE OF THE F.A.C. OVAL
C
      BX=0.d0
      BY=0.d0
      BZ=0.d0
      
      Xsc=X*xKappa
      Ysc=Y*xKappa
      Zsc=Z*xKappa

      M=MODE

C     Stretch the position
      call stretchpos(tilt, dPhi, Xsc,Ysc,Zsc, XS,YS,ZS)
C     Evaluate the field
      CALL twocones (A,theta0,dTheta,M,isEven,isSmoothed,XS,Ysc,ZS,
     .     BXS,BYAS,BZS)
C     Stretch the field
      call stretchfield(BXS,BYAS,BZS, BX,BY,BZ)
      
      BX=BX*xKappa              ! SCALING
      BY=BY*xKappa
      BZ=BZ*xKappa
      
      RETURN
      END
c
C=========================================================================
c
      SUBROUTINE twocones (A,theta0,dTheta,mode,isEven,isSmoothed,
     .     X,Y,Z, BX,BY,BZ)
C
C   ADDS FIELDS FROM TWO CONES (NORTHERN AND SOUTHERN), WITH A PROPER SYMMETRY OF THE CURRENT AND FIELD,
C     CORRESPONDING TO THE REGION 1 BIRKELAND CURRENTS. (SEE NB #6, P.58).
C
      IMPLICIT NONE
      
C     Inputs
      REAL*8  A(30)
      REAL*8  theta0
      REAL*8  dTheta
      INTEGER mode
      LOGICAL isEven
      LOGICAL isSmoothed
      REAL*8 x,y,z

c     Outputs
      REAL*8 bx,by,bz

c     Internal variables      
      REAL*8 bxn,byn,bzn
      REAL*8 bxs,bys,bzs

      if (isSmoothed .eqv. .false.) then
         call one_cone (A,theta0,dTheta,mode,isEven,X,Y,Z,BXN,BYN,BZN)
         call one_cone (A,theta0,dTheta,mode,isEven,X,Y,-Z,BXS,BYS,BZS)
      else
         call one_cone_smooth (A,theta0,dTheta,mode,isEven,X,Y,Z,
     .        BXN,BYN,BZN)
         call one_cone_smooth (A,theta0,dTheta,mode,isEven,X,Y,-Z,
     .        BXS,BYS,BZS)         
      endif
      
      BX=BXN-BXS
      BY=BYN-BYS
      BZ=BZN+BZS

      RETURN
      END

