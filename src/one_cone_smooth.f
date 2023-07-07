c
C-------------------------------------------------------------------------
C
      SUBROUTINE ONE_CONE_SMOOTH(A,theta0,dTheta,M,isEven, x,y,z,
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
      
c     Internal variables
      INTEGER n
      REAL*8 ddTheta,theta00
      REAL*8 bxD(5),byD(5),bzD(5)
      REAL*8 rho2,rho
      REAL*8 r,theta,phi
      REAL*8 rS,thetaS,phiS
      REAL*8 bThetaC,bPhiC
      REAL*8 br,btheta,bphi
      REAL*8 S,C,SF,CF,BE
      
c      
c  ////////SMOOTHING FAC modules/////////
c
      DO n=1,5
         ddTheta=0.2d0*dTheta
         theta00=theta0-2.d0*ddTheta+(n-1)*ddTheta !  dynamic THETA0
      
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
         CALL conical (M,theta00,dTheta,isEven,
     .        rS,thetaS,phiS, bThetaC,bPhiC) !   MODE #M
         
         bThetaC=bThetaC*800.0d0
         bPhiC=bPhiC*800.0d0
      
         call deformbirkfield(0.0d0,bThetaC,bPhiC, br,btheta,bphi)
         
         S=rho/r
         C=z/r
         SF=y/rho
         CF=x/rho
      
         BE=br*S+btheta*C
         
         bxD(n)=A(1)*(BE*CF-bphi*SF)
         byD(n)=A(1)*(BE*SF+bphi*CF)
         bzD(n)=A(1)*(BR*C-btheta*S)
         
      ENDDO
  
      bx=bxD(1)+bxD(2)+bxD(3)+bxD(4)+bxD(5)
      by=byD(1)+byD(2)+byD(3)+byD(4)+byD(5)
      bz=bzD(1)+bzD(2)+bzD(3)+bzD(4)+bzD(5)
       
      bx=bx/5.d0
      by=by/5.d0
      bz=bz/5.d0
c
c   /////////END OF SMOOTHING FAC modules////////////
c
       RETURN
       END
     
