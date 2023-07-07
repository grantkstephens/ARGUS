C==========================================================================
C     May 2019, G.K.Stephens, this subroutine was adabpted to the Saturn
C     magnetic field model.
C      
C     July 2017, G.K.Stephens, this routine updated to incorporate Jay Albert's
C     improvements to the Bessel function evaluation. The Bessel function
C     values are now precomputed and  passed into this rather than computed 
C     inside this routine.
C      
      SUBROUTINE tailsheet_sym (n,rho0,D0,AJM,
     .     x,y,z,
     .     bx,by,bz)

      IMPLICIT NONE

C     Inputs      
      INTEGER n
      REAL*8 rho0,D0
      REAL*8 x,y,z
      REAL*8 AJM(0:*)

c     Outputs
      REAL*8 bx,by,bz

c     Internal variables
      REAL*8 rho,cosPhi,sinPhi,kn,knZ,zDist,j0,j1,ex
      
c      COMMON /TAIL/ D  ! THE COMMON BLOCKS FORWARDS TAIL SHEET THICKNESS
c      DIMENSION AJM(0:5),AJMD(0:5)
C-----------------------------------------------------------------------------------
C
c      RNOT=20.0        !    This can be replaced by introducing them

      bx=0.d0
      by=0.d0
      bz=0.d0
      
      rho=DSQRT(x*x+y*y)
      cosPhi=x/rho
      sinPhi=y/rho
C
      kn=n/rho0
      knZ=kn*z

      zDist=DSQRT(z*z+D0*D0)
C
      j0=AJM(0)
      j1=AJM(1)
c     July 2017, G.K.Stephens, Bessel functions are now passed in.
c      RJ0=bessj0(RKMR)
c      RJ1=bessj1(RKMR)
      ex=DEXP(kn*zDist)
C
      bx=knZ*j1*cosPhi/zDist/ex
      by=knZ*j1*sinPhi/zDist/ex
      bz=kn*j0/ex
C
C    CALCULATION OF THE MAGNETOTAIL CURRENT CONTRIBUTION IS FINISHED
C     
      RETURN
      END
