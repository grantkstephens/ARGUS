C==========================================================================
C     May 2019, G.K.Stephens, this subroutine was adabpted to the Saturn
C     magnetic field model.
C      
C     July 2017, G.K.Stephens, this routine updated to incorporate Jay Albert's
C     improvements to the Bessel function evaluation. The Bessel function
C     values are now precomputed and  passed into this rather than computed 
C     inside this routine.
      SUBROUTINE tailsheet_asym (isEven,m,n,rho0,D0,AJM,AJMD,
     .     x,y,z,
     .     bx,by,bz)

      IMPLICIT NONE

C     Inputs
      LOGICAL isEven
      INTEGER M,N
      REAL*8 rho0,D0
      REAL*8 x,y,z
      REAL*8 AJM(0:*),AJMD(0:*)
      
c     Outputs
      REAL*8 bx,by,bz

c     Internal variables
      REAL*8 rho,phi,cosPhi,sinPhi,cosMPhi,sinMPhi,
     .     kn,knZ,knRho,zDist,j0,j1,ex,
     .     bRho,bPhi
      
C-----------------------------------------------------------------------------------
C
c      RNOT=20.0     !    Rho_0 - scale parameter along the tail axis
c      DLTK=1.0      !    step in Km
c
c -----------------------------------------------------------------------------------

      bx=0.d0
      by=0.d0
      bz=0.d0
      
      rho=DSQRT(x*x+y*y)
      cosPhi=x/rho
      sinPhi=y/rho

      phi=DATAN2(y,x)
      cosMPhi=DCOS(m*phi)
      sinMPhi=DSIN(m*phi)
C
      kn=n/rho0
C
      knRho=kn*rho
C
      zDist=DSQRT(z*z+D0*D0)
C
      ex=DEXP(kn*zDist)
c
c     July 2017, G.K.Stephens, Jm is now passed in, not computed internally
c ---- calculating Jm and its derivatives ------
c
c                   if(m.gt.2) then
c                   AJM=bessj(m,RKMR)
c                   AJM1=bessj(m-1,RKMR)
c                   AJMD=AJM1-m*AJM/RKMR
c                   else
c                  --------------------
c                   if(m.eq.2) then
c                   AJM=bessj(2,RKMR)
c                   AJM1=bessj1(RKMR)
c                   AJMD=AJM1-m*AJM/RKMR
c                   else
c                  --------------------
c                   AJM=bessj1(RKMR)
c                   AJM1=bessj0(RKMR)
c                   AJMD=AJM1-AJM/RKMR
c                  --------------------
c                   end if
c                  --------------------
c                   endif
c -----------------------------------------
c     
      if(isEven.eqv. .TRUE.) then
c -----------------------------------------
c calculating symmetric modes
c -----------------------------------------
         bRho = -m*sinMPhi*z*AJMD(m)/zDist/ex
         bPhi = -m*m*cosMPhi*z*AJM(m)/knRho/zDist/ex
         bz   = +m*sinMPhi*AJM(m)/ex
      else
c -----------------------------------------
c calculating asymmetric modes
c -----------------------------------------
c     
         bRho = +m*cosMPhi*z*AJMD(m)/zDist/ex
         bPhi = -m*m*sinMPhi*z*AJM(m)/knRho/zDist/ex
         bz   = -m*cosMPhi*AJM(m)/ex
      end if
c
c --- transformation from cylindrical ccordinates to GSM ---
c
      bx=bRho*cosPhi-bPhi*sinPhi
      by=bRho*sinPHi+bPhi*cosPhi
C
C    CALCULATION OF THE MAGNETOTAIL CURRENT CONTRIBUTION IS FINISHED
C
      RETURN
      END
