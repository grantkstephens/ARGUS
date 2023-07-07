C==========================================================================
C     July 2017, G.K.Stephens, this routine updated to incorporate Jay Albert's
C     improvements to the Bessel function evaluation.
      SUBROUTINE cylharmonic(isEven, TS,
     .     x,y,z,
     .     fx,fy,fz)            ! modified SHTBNORM_E
      IMPLICIT NONE

C     Inputs
      LOGICAL isEven
      REAL*8 TS(80)
      REAL*8 x,y,z

c     Outputs
      REAL*8 fx,fy,fz

c     Internal variables
      INTEGER Mend,Nend
      integer m,n,L1
      REAL*8 rhoI
      REAL*8 phi
      REAL*8 CHZ,SHZ
      REAL*8 HX,HX1,HX2
      REAL*8 HY,HY1,HY2
      REAL*8 HZ
      REAL*8 DPDX,DPDY
      REAL*8 cosMPhi,sinMPhi
      REAL*8 AKN,AKNR,AKNRI
      REAL*8 rho,knRho, AJM(0:14),AJMD(0:14)
      REAL*8 AK(5)
            
      AK(1)=TS(76)
      AK(2)=TS(77)
      AK(3)=TS(78)
      AK(4)=TS(79)
      AK(5)=TS(80)

      phi=DATAN2(Y,X)
      rho=dsqrt(x*x+y*y)
      
      if(RHO.lt.1.D-8) then
         rhoI=1.D8
      else
         rhoI=1.D0/rho
      endif
      DPDX=-Y*rhoI*rhoI
      DPDY= X*rhoI*rhoI
      
      fx=0.d0
      fy=0.d0
      fz=0.d0

      Mend = 14
      Nend = 5
      
      DO n=1,Nend
         AKN=dabs(AK(n))
         AKNR=AKN*rho
         if(AKNR.lt.1.D-8) then
            AKNRI=1.D8
         else
            AKNRI=1.D0/AKNR
         endif
         
         CHZ=dcosh(Z*AKN)
         SHZ=dsinh(Z*AKN)
         
         call bessJJ(14,AKNR, AJM) !!! get all n in one call
         DO m=1,14
            AJMD(m)=AJM(m-1)-m*AJM(m)*AKNRI
         ENDDO
         AJMD(0)=-AJM(1)

         DO m=0,Mend

            cosMPhi=DCOS(m*phi)
            sinMPhi=DSIN(m*phi)
            
            if(isEven .eqv. .true.) then
               HX1=-m*DPDX*cosMPhi*SHZ*AJM(m)
               HX2=-AKN*X*rhoI*sinMPhi*SHZ*AJMD(m)
               HX=HX1+HX2
            
               HY1=-m*DPDY*cosMPhi*SHZ*AJM(m)
               HY2=-AKN*Y*rhoI*sinMPhi*SHZ*AJMD(m)
               HY=HY1+HY2
            
               HZ=-AKN*sinMPhi*CHZ*AJM(m)
            else
               HX1=m*DPDX*sinMPhi*SHZ*AJM(m)
               HX2=-AKN*X*RHOI*cosMPhi*SHZ*AJMD(m)
               HX=HX1+HX2

               HY1=m*DPDY*sinMPhi*SHZ*AJM(m)
               HY2=-AKN*Y*RHOI*cosMPhi*SHZ*AJMD(m)
               HY=HY1+HY2
               
               HZ=-AKN*cosMPhi*CHZ*AJM(m)
            endif
            
            L1=n+5*m
            
            FX=FX+HX*TS(L1)
            FY=FY+HY*TS(L1)
            FZ=FZ+HZ*TS(L1)
         ENDDO
      ENDDO
      RETURN
      END
