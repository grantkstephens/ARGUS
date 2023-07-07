     
      SUBROUTINE warppos(tilt, G,TW, x,y,z, xs,ys,zs)
C
C   CALCULATES GSM COMPONENTS OF THE WARPED FIELD FOR TWO TAIL UNIT MODES.
C   THE WARPING DEFORMATION IS IMPOSED ON THE UNWARPED FIELD, COMPUTED
C   BY THE S/R "UNWARPED".  THE WARPING PARAMETERS WERE TAKEN FROM THE
C   RESULTS OF GEOTAIL OBSERVATIONS (TSYGANENKO ET AL. [1998]).
C   NB # 6, P.106, OCT 12, 2000.
C
      IMPLICIT NONE
      
C     Inputs
      REAL*8 tilt
      REAL*8 G,TW
      REAL*8 x,y,z
      
c     Outputs
      REAL*8 xs,ys,zs

c     Internal variables
      LOGICAL initialized /.false./
      REAL*8 dGdX,XL,dXLdX
      REAL*8 rho,rho2
      REAL*8 phi,cosPhi,sinPhi
      REAL*8 F,cosF,sinF
      REAL*8 SPS
      REAL*8 RR4L4
      REAL*8 dFdPhi,dFdRho, dFDx

c     saved
c     save xs,ys,zs
      save initialized
      save dGdX,XL,dXLdX
      save rho,rho2
      save phi,cosPhi,sinPhi
      save F,cosF,sinF
      save SPS
      save RR4L4
      save dFdPhi,dFdRho, dFDx
      
C     Inputs
      REAL*8 bx,by,bz
      
c     Outputs
      REAL*8 bxs,bys,bzs

c     Internal variables
      REAL*8 BRHO_AS,BPHI_AS, BRHO_S,BPHI_S

      initialized=.true.
      
      dGdX=0.D0
      XL=20.D0
      dXLdX=0.D0

      SPS=DSIN(tilt)
      rho2=Y**2+Z**2
      rho=DSQRT(rho2)

      IF (Y.EQ.0.D0.AND.Z.EQ.0.D0) THEN
       phi=0.D0
       cosPhi=1.D0
       sinPhi=0.D0
      ELSE
       phi=DATAN2(Z,Y)
       cosPhi=Y/rho
       sinPhi=Z/rho
      ENDIF

      RR4L4=rho/(RHO2**2+XL**4)

      F=phi+G*rho2*RR4L4*cosPhi*SPS +TW*(X/10.D0)
      
      dFdPhi=1.D0-G*rho2*RR4L4*sinPhi*SPS
      dFdRho=G*RR4L4**2*(3.D0*XL**4-rho2**2)*cosPhi*SPS
      
      dFDx=RR4L4*cosPhi*SPS*(dGdX*RHO2-G*rho*RR4L4*4.D0*XL**3*DXLDX)
     *  +TW/10.D0        !  THE LAST TERM DESCRIBES THE IMF-INDUCED TWISTING (ADDED 04/21/06)

      cosF=DCOS(F)
      sinF=DSIN(F)
            
      xs=x
      ys=rho*cosF
      zs=rho*sinF
c      YAS=rho*CF
c      ZAS=rho*SF
      
      return

      entry warpfield(bx,by,bz, bxs,bys,bzs)

c     check that warppos was called before warpfield
      if (initialized .eqv. .false.) then
         stop 'warpfield must be called after '//
     .         'warppos to initialize the fields'
      endif

      BRHO_AS =  by*cosF+bz*sinF
      BPHI_AS = -by*sinF+bz*cosF
            
      BRHO_S = BRHO_AS*dFdPhi
      BPHI_S = BPHI_AS-rho*(BX*DFDX+BRHO_AS*dFdRho)
      
      bxs=bx*DFDPHI
      bys=BRHO_S*cosPhi-BPHI_S*sinPhi
      bzs=BRHO_S*sinPhi+BPHI_S*cosPhi
      
      return
      
      end
