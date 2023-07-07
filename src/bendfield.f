      SUBROUTINE bendpos(tilt, RH0, x,y,z, xs,ys,zs)

      IMPLICIT NONE

C     Inputs
      REAL*8 tilt
      REAL*8 RH0
      REAL*8 x,y,z
      
c     Outputs
      REAL*8 xs,ys,zs

c     Constants
      REAL*8 RH2
      PARAMETER (RH2=-5.2d0)
      INTEGER IEPS
      PARAMETER (IEPS=3)
      
c     Internal variables
      LOGICAL initialized /.false./
      REAL*8 SPS,CPS
      REAL*8 R2,R,ZR,RH
      REAL*8 DRHDR,DRHDZ
      REAL*8 RRH,F,DFDR,DFDRH
      REAL*8 SPSAS,CPSAS
      REAL*8 FACPS,PSASX,PSASY,PSASZ
      REAL*8 DXASDX,DXASDY,DXASDZ,DZASDX,DZASDY,DZASDZ
      REAL*8 FAC1,FAC2,FAC3
      
c     saved
c     save xs,ys,zs
      save initialized
      save SPS,CPS
      save R2,R,ZR,RH
      save DRHDR,DRHDZ
      save RRH,F,DFDR,DFDRH
      save SPSAS,CPSAS
      save FACPS,PSASX,PSASY,PSASZ
      save DXASDX,DXASDY,DXASDZ,DZASDX,DZASDY,DZASDZ
      save FAC1,FAC2,FAC3
      
C     Inputs
      REAL*8 bx,by,bz
      
c     Outputs
      REAL*8 bxs,bys,bzs

c     Internal variables
      
      initialized=.true.
      
C
C  RH0,RH1,RH2, AND IEPS CONTROL THE TILT-RELATED DEFORMATION OF THE TAIL FIELD
C
      SPS=DSIN(tilt)
      CPS=DSQRT(1.D0-SPS**2)
      
      R2=X**2+Y**2+Z**2
      R=SQRT(R2)
      ZR=Z/R
      RH=RH0+RH2*ZR**2
      
      DRHDR=-ZR/R*2.D0*RH2*ZR
      DRHDZ= 2.D0*RH2*ZR/R
C
      RRH=R/RH
      F=1.D0/(1.D0+RRH**IEPS)**(1.D0/IEPS)
      DFDR=-RRH**(IEPS-1)*F**(IEPS+1)/RH
      DFDRH=-RRH*DFDR
c
      SPSAS=SPS*F
      CPSAS=DSQRT(1.D0-SPSAS**2)
C
      xs=X*CPSAS-Z*SPSAS
      ys=y
      zs=X*SPSAS+Z*CPSAS
C
      FACPS=SPS/CPSAS*(DFDR+DFDRH*DRHDR)/R
      PSASX=FACPS*X
      PSASY=FACPS*Y
      PSASZ=FACPS*Z+SPS/CPSAS*DFDRH*DRHDZ
C
      DXASDX=CPSAS-zs*PSASX
      DXASDY=-zs*PSASY
      DXASDZ=-SPSAS-zs*PSASZ
      DZASDX=SPSAS+xs*PSASX
      DZASDY=xs*PSASY
      DZASDZ=CPSAS+xs*PSASZ
      
      FAC1=DXASDZ*DZASDY-DXASDY*DZASDZ
      FAC2=DXASDX*DZASDZ-DXASDZ*DZASDX
      FAC3=DZASDX*DXASDY-DXASDX*DZASDY

      return
      
      entry bendfield(bx,by,bz, bxs,bys,bzs)
      
c     check that warppos was called before warpfield
      if (initialized .eqv. .false.) then
         stop 'bendfield must be called after '//
     .        'bendpos to initialize the fields'
      endif
      
      bxs=bx*DZASDZ-bz*DXASDZ+by*FAC1
      bys=by*FAC2
      bzs=bz*DXASDX-bx*DZASDX+by*FAC3
         
      return
      
      end
