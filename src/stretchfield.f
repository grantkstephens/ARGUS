      SUBROUTINE stretchpos(tilt, dPhi, x,y,z, xs,ys,zs)

      IMPLICIT NONE

C     Inputs
      REAL*8 tilt,PS
      REAL*8 dPhi
      REAL*8 x,y,z

c     Outputs
      REAL*8 xs,ys,zs

c     Constants
! THESE PARAMETERS CONTROL DAY-NIGHT ASYMMETRY OF F.A.C., AS FOLLOWS:      
C  (2) B:      AN ASYMMETRY FACTOR AT HIGH-ALTITUDES;  FOR B=0, THE ONLY ASYMMETRY IS THAT FROM DPHI
C              TYPICAL VALUES: 0.35-0.70
C  (3) RHO_0:  A FIXED PARAMETER, DEFINING THE DISTANCE RHO, AT WHICH THE LATITUDE SHIFT GRADUALLY SATURATES AND
C              STOPS INCREASING
C              ITS VALUE WAS ASSUMED FIXED, EQUAL TO 7.0.      
! parameters of the tilt-dependent deformation of the untilted F.A.C. field      
      REAL*8 BETA
      PARAMETER (BETA=0.9d0)
      REAL*8 RH
      PARAMETER (RH=10.d0)
      REAL*8 EPS
      PARAMETER (EPS=3.d0)
      REAL*8 B
      PARAMETER (B=0.5) ! not D.P.
      REAL*8 RHO_0
      PARAMETER (RHO_0=7.0) ! not D.P.
      
c     Internal variables
      LOGICAL initialized /.false./
      REAL*8 RHO,Rsc,RHO2,PHI
      REAL*8 SPHIC,CPHIC
      REAL*8 BRACK,R1RH,PSIAS
      REAL*8 PHIS,DPHISPHI,DPHISRHO,DPHISDY
      REAL*8 SPHICS,CPHICS
      
c     saved
c     save xs,ys,zs
      save initialized
      save RHO,Rsc,RHO2,PHI
      save SPHIC,CPHIC
      save BRACK,R1RH,PSIAS
      save PHIS,DPHISPHI,DPHISRHO,DPHISDY
      save SPHICS,CPHICS
      
C     Inputs
      REAL*8 bx,by,bz
      
c     Outputs
      REAL*8 bxs,bys,bzs

c     Internal variables
      REAL*8 BPHI_S,BPHIAS,BRHO_S,BRHOAS
      
      initialized=.true.
      PS=tilt

      RHO=DSQRT(x**2+z**2)

      Rsc=DSQRT(x**2+y**2+z**2)                                 !  SCALED
      RHO2=RHO_0**2

      IF (x.EQ.0.D0.AND.z.EQ.0.D0) THEN
         PHI=0.D0
      ELSE
         PHI=DATAN2(-z,x)  !  FROM CARTESIAN TO CYLINDRICAL (RHO,PHI,Y)
      ENDIF

      SPHIC=DSIN(PHI)
      CPHIC=DCOS(PHI)  !  "C" means "CYLINDRICAL", TO DISTINGUISH FROM SPHERICAL PHI

      BRACK=DPHI+B*RHO2/(RHO2+1.D0)*(RHO**2-1.D0)/(RHO2+RHO**2)
      R1RH=(Rsc-1.D0)/RH
      PSIAS=BETA*PS/(1.D0+R1RH**EPS)**(1.D0/EPS)

      PHIS=PHI-BRACK*DSIN(PHI) -PSIAS
      DPHISPHI=1.D0-BRACK*DCOS(PHI)
      DPHISRHO=-2.D0*B*RHO2*RHO/(RHO2+RHO**2)**2 *DSIN(PHI)
     *   +BETA*PS*R1RH**(EPS-1.D0)*RHO/(RH*Rsc*
     *   (1.D0+R1RH**EPS)**(1.D0/EPS+1.D0))
      DPHISDY= BETA*PS*R1RH**(EPS-1.D0)*y/(RH*Rsc*
     *   (1.D0+R1RH**EPS)**(1.D0/EPS+1.D0))

      SPHICS=DSIN(PHIS)
      CPHICS=DCOS(PHIS)

      XS= RHO*CPHICS
      YS=y
      ZS=-RHO*SPHICS
      
      return
      
      entry stretchfield(bx,by,bz, bxs,bys,bzs)
      
c     check that warppos was called before warpfield
      if (initialized .eqv. .false.) then
         stop 'deformbirkfield must be called after '//
     .        'deformbirkpos to initialize the fields'
      endif

      BRHOAS=bx*CPHICS-bz*SPHICS
      BPHIAS=-bx*SPHICS-bz*CPHICS

      BRHO_S=BRHOAS*DPHISPHI                          
      BPHI_S=(BPHIAS-RHO*(by*DPHISDY+BRHOAS*DPHISRHO))
      bys=by*DPHISPHI

      bxs=BRHO_S*CPHIC-BPHI_S*SPHIC
      bzs=-BRHO_S*SPHIC-BPHI_S*CPHIC

      return
      
      end
