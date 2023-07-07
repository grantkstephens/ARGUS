c
C-------------------------------------------------------------------------
C
      SUBROUTINE one_cone(A,theta0,dTheta,M,isEven, x,y,z,
     .     bx,by,bz)
c
c  RETURNS FIELD COMPONENTS FOR A DEFORMED CONICAL CURRENT SYSTEM, FITTED TO A BIOSAVART FIELD
c    BY SIM_14.FOR.  HERE ONLY THE NORTHERN CONE IS TAKEN INTO ACCOUNT.
c

      IMPLICIT NONE

C     Inputs
      REAL*8  A(30)
      REAL*8  theta0
      REAL*8  dTheta
      INTEGER M
      LOGICAL isEven
      REAL*8  x,y,z
      
C     Outputs
      REAL*8  bx,by,bz

C     Internal variables      
      REAL*8 rho2,rho
      REAL*8 r,theta,phi
      REAL*8 rS,thetaS,phiS
      REAL*8 bThetaC,bPhiC
      REAL*8 br,btheta,bphi
      REAL*8 S,C,SF,CF,BE
      
      rho2=x**2+y**2
      rho=DSQRT(rho2)
      
      r=DSQRT(rho2+z**2)
      theta=DATAN2(rho,z)
      phi=DATAN2(y,x)
C     
C     MAKE THE DEFORMATION OF COORDINATES:
C     
      call deformbirkpos(A, r,theta,phi, rS,thetaS,phiS)
      
C
C   CALCULATE FIELD COMPONENTS AT THE NEW POSITION (ASTERISKED):
C     
      CALL conical (M,theta0,dTheta,isEven,
     .     rS,thetaS,phiS, bThetaC,bPhiC) !   MODE #M

C     Scale the field values by 800.0
      bThetaC=bThetaC*800.0d0
      bPhiC=bPhiC*800.0d0

C     Now deform the field
      call deformbirkfield(0.0d0,bThetaC,bPhiC, br,btheta,bphi)

C     Convert the field from spherical to Cartesian
      S=rho/r
      C=z/r
      SF=y/rho
      CF=x/rho
      
      BE=br*S+btheta*C

C     Scale by the amplitude coefficient
      bx=A(1)*(BE*CF-bphi*SF)
      by=A(1)*(BE*SF+bphi*CF)
      bz=A(1)*(BR*C-btheta*S)
      
      RETURN
      END
