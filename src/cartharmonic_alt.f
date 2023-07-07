      SUBROUTINE cartharmonic_alt (isIeven,isKeven,iEnd,kEnd,
     .     pis,pks,aiks,
     .     x,y,z,
     .     bx,by,bz)
c      
c     Computes the field for the following Cartesian scalar potential
c     solution to Laplace's equation:
c     Uik = aik*exp(x*sqrt(pi^2+pk^2))*sin(pi*y)*sin(pk*z)
c     
c     Inputs:
c     logical isIeven: Determines the trig parity for the Y terms
c       (odd=sine, even=cosine).
c     logical isKeven: Determines the trig parity for the Z terms
c       (odd=sine, even=cosine).
c     integer iEnd: The ending index (or size) of the Y expansion.
c     integer kEnd: The ending index (or size) of the Z expansion.
c     real*8 pis: An array containing the nonlinear parameters
c       associated with the Y expansion.
c     real*8 pks: An array containing the nonlinear parameters
c       associated with the Z expansion.
c     real*8 aiks: An array containing the linear scaling coefficients.
c     real*8 x,y,z: The location to evaluate the value of the magnetic
c       field.
c      
c     Outputs:
c     real*8 bx,by,bz: The magnetc field evaluated at the supplied
c       location.
c      
      IMPLICIT NONE

C     Inputs
      LOGICAL isIeven,isKeven
      INTEGER iEnd,kEnd
      REAL*8 pis(iEnd)
      REAL*8 pks(kEnd)
      REAL*8 aiks(iEnd,kEnd)
      REAL*8 x,y,z
      
c     Outputs
      REAL*8 bx,by,bz

c     Internal variables
      integer i,k
      REAL*8 pi,sinYpi,cosYpi,pk,sqrtP,exp,sinZpk,cosZpk,aik 

      DO i=1,Iend
         
         pi = pis(i)

         if(isIeven .eqv. .TRUE.) then         
            sinYpi = DCOS(pi*y)
            cosYpi =-DSIN(pi*y)
         else
            sinYpi = DSIN(pi*y)
            cosYpi = DCOS(pi*y)
         end if
         
         DO k=1,Kend

            pk = pks(k)
            
            sqrtP = DSQRT(pi**2+pk**2)

            exp = DEXP(x*sqrtP)

            if(isKeven .eqv. .TRUE.) then         
               sinZpk = DCOS(pk*z)
               cosZpk =-DSIN(pk*z)
            else
               sinZpk = DSIN(pk*z)
               cosZpk = DCOS(pk*z)
            end if

            aik = aiks(i,k)

            if (k .eq. Kend) then

               bx = bx - aik*exp*sinYpi*
     .              (sqrtP*z*cosZpk+sinZpk*pk*(x+1.0d0/sqrtP))
               by = by - aik*exp*pi*cosYpi*
     .              (z*cosZpk+x*pk*sinZpk/sqrtP)
               bz = bz - aik*exp*sinYpi*
     .              (cosZpk*(1.0d0+x*pk*pk/sqrtP)-z*pk*sinZpk)
               
            else
               bx = bx - aik*exp*sqrtP*sinYpi*sinZpk
               by = by - aik*exp*pi*cosYpi*sinZpk
               bz = bz - aik*exp*pk*sinYpi*cosZpk
            end if
               
         ENDDO
      ENDDO
      
      RETURN
      END
