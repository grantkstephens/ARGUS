      SUBROUTINE deformbirkpos(A, r,theta,phi, rS,thetaS,phiS)

      IMPLICIT NONE

C     Inputs
      REAL*8 A(31)
      REAL*8 r,theta,phi

c     Outputs
      REAL*8 rS,thetaS,phiS

c     Constants
      REAL*8 DR
      PARAMETER (DR=1.D-6)        !   JUST FOR NUMERICAL DIFFERENTIATION
      REAL*8 DT
      PARAMETER (DT=1.D-6)

c     Functions
      REAL*8 R_S,THETA_S
      
c     Internal variables
      LOGICAL initialized /.false./
      REAL*8 RC
      REAL*8 DRSDR,DRSDT,DTSDR,DTSDT,STSST,RSR
      
c     saved
c     save xs,ys,zs
      save initialized
      save RC
      save DRSDR,DRSDT,DTSDR,DTSDT,STSST,RSR
      
C     Inputs
      REAL*8 br,btheta,bphi
      
c     Outputs
      REAL*8 brS,bthetaS,bphiS

c     Internal variables
      
      initialized=.true.
      RC = R ! copy R so we can save it for use later
C
C   MAKE THE DEFORMATION OF COORDINATES:
C
       rS=R_S(A,R,THETA)
       thetaS=THETA_S(A,R,THETA)
       phiS=PHI

C     
C   NOW TRANSFORM B{R,T,F}_AST BY THE DEFORMATION TENSOR:
C
C      FIRST OF ALL, FIND THE DERIVATIVES:
C
       DRSDR=(R_S(A,R+DR,THETA)-R_S(A,R-DR,THETA))/(2.D0*DR)
       DRSDT=(R_S(A,R,THETA+DT)-R_S(A,R,THETA-DT))/(2.D0*DT)
       DTSDR=(THETA_S(A,R+DR,THETA)-THETA_S(A,R-DR,THETA))/(2.D0*DR)
       DTSDT=(THETA_S(A,R,THETA+DT)-THETA_S(A,R,THETA-DT))/(2.D0*DT)
       STSST=DSIN(THETAS)/DSIN(THETA)
       RSR=RS/R
       
      return
      
      entry deformbirkfield(br,btheta,bphi, brS,bthetaS,bphiS)
      
c     check that warppos was called before warpfield
      if (initialized .eqv. .false.) then
         stop 'deformbirkfield must be called after '//
     .        'deformbirkpos to initialize the fields'
      endif
      
      brS     =-RSR/RC*STSST*btheta*DRSDT !   NB#6, P.43    BRAST DOES NOT ENTER HERE
      bthetaS = RSR*STSST*btheta*DRSDR !               (SINCE IT IS ZERO IN OUR CASE)
      bphiS   = RSR*bphi*(DRSDR*DTSDT-DRSDT*DTSDR)

      return
      
      end
C
C=====================================================================================
      DOUBLE PRECISION FUNCTION R_S(A,R,THETA)
      IMPLICIT NONE
C     Inputs
      REAL*8 A(31)
      REAL*8 r,theta
C
      R_S=R+A(2)/R+A(3)*R/DSQRT(R**2+A(11)**2)+A(4)*R/(R**2+A(12)**2)
     *+(A(5)+A(6)/R+A(7)*R/DSQRT(R**2+A(13)**2)+A(8)*R/(R**2+A(14)**2))*
     * DCOS(THETA)
     *+(A(9)*R/DSQRT(R**2+A(15)**2)+A(10)*R/(R**2+A(16)**2)**2)
     * *DCOS(2.D0*THETA)
C
      RETURN
      END
C
C-----------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION THETA_S(A,R,THETA)
      IMPLICIT NONE
C     Inputs
      REAL*8 A(31)
      REAL*8 r,theta
c
      THETA_S=THETA+(A(17)+A(18)/R+A(19)/R**2
     *                +A(20)*R/DSQRT(R**2+A(27)**2))*DSIN(THETA)
     * +(A(21)+A(22)*R/DSQRT(R**2+A(28)**2)
     *                +A(23)*R/(R**2+A(29)**2))*DSIN(2.D0*THETA)
     * +(A(24)+A(25)/R+A(26)*R/(R**2+A(30)**2))*DSIN(3.D0*THETA)
C
      RETURN
      END
