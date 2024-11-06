!
!© 2003 The Johns Hopkins University Applied Physics Laboratory LLC
!
      subroutine store_sst19_static(staticDir)

      IMPLICIT NONE
      
C     Inputs      
      CHARACTER*256 staticDir
      CHARACTER*256 variableDir
      INTEGER iYear,iDayOfYear,iHour,min
      CHARACTER*256 filename

C     Constants
      INTEGER    Mend,Nend
      PARAMETER (Mend = 6)
      PARAMETER (Nend = 8)
      REAL*8     RE,mu0
      PARAMETER (RE = 6371.2d0)
      PARAMETER (mu0 = 1.25663706127E-6)
      
C     Internal variables
C     time-static variables
      LOGICAL init_static /.false./
      REAL*8 TSS(80,Nend)
      REAL*8 TSO(80,Mend,Nend)
      REAL*8 TSE(80,Mend,Nend)   
c     time-variable variables
      LOGICAL init_variable /.false./
      REAL*8 A_dipsh
      REAL*8 A_eq(Nend+2*Mend*Nend)
      REAL*8 A_eq_P(Nend+2*Mend*Nend)
      REAL*8 A_eq_TCS(Nend+2*Mend*Nend)
      REAL*8 A_eq_P_TCS(Nend+2*Mend*Nend)
      REAL*8 A_R11am,A_R12am,A_R11sm,A_R12sm
      REAL*8 A_R11a, A_R12a, A_R11s, A_R12s
      REAL*8 A_R21am,A_R22am,A_R21sm,A_R22sm
      REAL*8 A_R21a, A_R22a, A_R21s, A_R22s
      REAL*8 D,DTCS,RH0,G
      REAL*8 xKappa1,xKappa2
      REAL*8 TW
      REAL*8 Pdyn
      REAL*8 tilt
      
c     saved
      save TSS,TSO,TSE
      save A_dipsh
      save A_eq, A_eq_P
      save A_eq_TCS, A_eq_P_TCS
      save A_R11am,A_R12am,A_R11sm,A_R12sm
      save A_R11a, A_R12a, A_R11s, A_R12s
      save A_R21am,A_R22am,A_R21sm,A_R22sm
      save A_R21a, A_R22a, A_R21s, A_R22s
      save D,DTCS,RH0,G
      save xKappa1,xKappa2
      save TW
      save Pdyn
      save tilt
      
C     Inputs
      REAL*8 tiltIn,PdynIn
      REAL*8 x,y,z

c     Outputs
      REAL*8 bx,by,bz
      REAL*8 jx,jy,jz

c     Internal variables
      REAL*8  BXCF,BYCF,BZCF
      REAL*8  BXTS1(Nend),BXTO1(Mend,Nend),BXTE1(Mend,Nend)
      REAL*8  BYTS1(Nend),BYTO1(Mend,Nend),BYTE1(Mend,Nend)
      REAL*8  BZTS1(Nend),BZTO1(Mend,Nend),BZTE1(Mend,Nend)
      REAL*8  BXTS2(Nend),BXTO2(Mend,Nend),BXTE2(Mend,Nend)
      REAL*8  BYTS2(Nend),BYTO2(Mend,Nend),BYTE2(Mend,Nend)
      REAL*8  BZTS2(Nend),BZTO2(Mend,Nend),BZTE2(Mend,Nend)
      REAL*8  BX11am,BY11am,BZ11am
      REAL*8  BX12am,BY12am,BZ12am
      REAL*8  BX11sm,BY11sm,BZ11sm
      REAL*8  BX12sm,BY12sm,BZ12sm
      REAL*8  BX11a,BY11a,BZ11a
      REAL*8  BX12a,BY12a,BZ12a
      REAL*8  BX11s,BY11s,BZ11s
      REAL*8  BX12s,BY12s,BZ12s
      REAL*8  BX21am,BY21am,BZ21am
      REAL*8  BX22am,BY22am,BZ22am
      REAL*8  BX21sm,BY21sm,BZ21sm
      REAL*8  BX22sm,BY22sm,BZ22sm
      REAL*8  BX21a,BY21a,BZ21a
      REAL*8  BX22a,BY22a,BZ22a
      REAL*8  BX21s,BY21s,BZ21s
      REAL*8  BX22s,BY22s,BZ22s
c     Internal variables for current
      REAL*8  dr,Acurr
      REAL*8  xp,yp,zp
      REAL*8  xm,ym,zm
      REAL*8  bx_xp,by_xp,bz_xp
      REAL*8  bx_xm,by_xm,bz_xm
      REAL*8  bx_yp,by_yp,bz_yp
      REAL*8  bx_ym,by_ym,bz_ym
      REAL*8  bx_zp,by_zp,bz_zp
      REAL*8  bx_zm,by_zm,bz_zm
      
      call read_sst19_static(staticDir,
     .     TSS,TSO,TSE)
      init_static=.true.      
      return

      entry store_sst19_variable(filename)
      call read_sst19_variable(filename,
     .     A_dipsh, A_eq,A_eq_P,A_eq_TCS,A_eq_P_TCS,
     .     A_R11am,A_R12am,A_R11sm,A_R12sm,
     .     A_R11a, A_R12a, A_R11s, A_R12s,
     .     A_R21am,A_R22am,A_R21sm,A_R22sm,
     .     A_R21a, A_R22a, A_R21s, A_R22s,
     .     D,DTCS, RH0,G, xKappa1,xKappa2, TW,
     .     Pdyn,tilt)
      init_variable=.true.
      return
      

      entry store_sst19_variable_datetime(variableDir,
     .     iYear,iDayOfYear,iHour,min)
      call read_sst19_variable_datetime(variableDir,
     .     iYear,iDayOfYear,iHour,min,
     .     A_dipsh, A_eq,A_eq_P,A_eq_TCS,A_eq_P_TCS,
     .     A_R11am,A_R12am,A_R11sm,A_R12sm,
     .     A_R11a, A_R12a, A_R11s, A_R12s,
     .     A_R21am,A_R22am,A_R21sm,A_R22sm,
     .     A_R21a, A_R22a, A_R21s, A_R22s,
     .     D,DTCS, RH0,G, xKappa1,xKappa2, TW,
     .     Pdyn,tilt)
      init_variable=.true.
      return

      
      entry sst19_field_tilt(tiltIn,x,y,z, bx,by,bz)
      
c     check that the static coefficients were loaded
      if (init_static .eqv. .false.) then
         stop 'the static shielding coefficients must be '//
     .        'stored first by calling store_sst19_static'
      endif      

c     check that the time-variable variables were loaded
      if (init_variable .eqv. .false.) then
         stop 'the time-variable parameters and coefficients '//
     .        'must be stored first by calling store_sst19_variable'
      endif
      
      call sst19_field_vect(TSS,TSO,TSE,
     .     A_dipsh, A_eq,A_eq_P, A_eq_TCS,A_eq_P_TCS,
     .     A_R11am,A_R12am,A_R11sm,A_R12sm,
     .     A_R11a, A_R12a, A_R11s, A_R12s,
     .     A_R21am,A_R22am,A_R21sm,A_R22sm,
     .     A_R21a, A_R22a, A_R21s, A_R22s,
     .     D,DTCS, RH0,G, xKappa1,xKappa2, TW,
     .     tiltIn,Pdyn,X,Y,Z,
     .     BX,BY,BZ)
      
      return
      
      
      entry sst19_field(x,y,z, bx,by,bz)
      
c     check that the static coefficients were loaded
      if (init_static .eqv. .false.) then
         stop 'the static shielding coefficients must be '//
     .        'stored first by calling store_sst19_static'
      endif      

c     check that the time-variable variables were loaded
      if (init_variable .eqv. .false.) then
         stop 'the time-variable parameters and coefficients '//
     .        'must be stored first by calling store_sst19_variable'
      endif
      
      call sst19_field_vect(TSS,TSO,TSE,
     .     A_dipsh, A_eq,A_eq_P, A_eq_TCS,A_eq_P_TCS,
     .     A_R11am,A_R12am,A_R11sm,A_R12sm,
     .     A_R11a, A_R12a, A_R11s, A_R12s,
     .     A_R21am,A_R22am,A_R21sm,A_R22sm,
     .     A_R21a, A_R22a, A_R21s, A_R22s,
     .     D,DTCS, RH0,G, xKappa1,xKappa2, TW,
     .     tilt,Pdyn,X,Y,Z,
     .     BX,BY,BZ)
      
      return

      
      entry sst19_current(x,y,z, jx,jy,jz)
      
c     check that the static coefficients were loaded
      if (init_static .eqv. .false.) then
         stop 'the static shielding coefficients must be '//
     .        'stored first by calling store_sst19_static'
      endif      

c     check that the time-variable variables were loaded
      if (init_variable .eqv. .false.) then
         stop 'the time-variable parameters and coefficients '//
     .        'must be stored first by calling store_sst19_variable'
      endif
      
      dr = 0.01d0

c     converts to current density in units of nA/m^2
      Acurr = 1.0d0/(2.0d0*dr)/(mu0*RE*1000.0d0)

c     shift the input coord. by the finite difference
      xp = x+dr
      xm = x-dr
      yp = y+dr
      ym = y-dr
      zp = z+dr
      zm = z-dr

c     evaluate the field at the displaced coords.
      call sst19_field_vect(TSS,TSO,TSE,
     .     A_dipsh, A_eq,A_eq_P, A_eq_TCS,A_eq_P_TCS,
     .     A_R11am,A_R12am,A_R11sm,A_R12sm,
     .     A_R11a, A_R12a, A_R11s, A_R12s,
     .     A_R21am,A_R22am,A_R21sm,A_R22sm,
     .     A_R21a, A_R22a, A_R21s, A_R22s,
     .     D,DTCS, RH0,G, xKappa1,xKappa2, TW,
     .     tilt,Pdyn,xp,y,z, bx_xp,by_xp,bz_xp)
      call sst19_field_vect(TSS,TSO,TSE,
     .     A_dipsh, A_eq,A_eq_P, A_eq_TCS,A_eq_P_TCS,
     .     A_R11am,A_R12am,A_R11sm,A_R12sm,
     .     A_R11a, A_R12a, A_R11s, A_R12s,
     .     A_R21am,A_R22am,A_R21sm,A_R22sm,
     .     A_R21a, A_R22a, A_R21s, A_R22s,
     .     D,DTCS, RH0,G, xKappa1,xKappa2, TW,
     .     tilt,Pdyn,xm,y,z, bx_xm,by_xm,bz_xm)
      call sst19_field_vect(TSS,TSO,TSE,
     .     A_dipsh, A_eq,A_eq_P, A_eq_TCS,A_eq_P_TCS,
     .     A_R11am,A_R12am,A_R11sm,A_R12sm,
     .     A_R11a, A_R12a, A_R11s, A_R12s,
     .     A_R21am,A_R22am,A_R21sm,A_R22sm,
     .     A_R21a, A_R22a, A_R21s, A_R22s,
     .     D,DTCS, RH0,G, xKappa1,xKappa2, TW,
     .     tilt,Pdyn,x,yp,z, bx_yp,by_yp,bz_yp)
      call sst19_field_vect(TSS,TSO,TSE,
     .     A_dipsh, A_eq,A_eq_P, A_eq_TCS,A_eq_P_TCS,
     .     A_R11am,A_R12am,A_R11sm,A_R12sm,
     .     A_R11a, A_R12a, A_R11s, A_R12s,
     .     A_R21am,A_R22am,A_R21sm,A_R22sm,
     .     A_R21a, A_R22a, A_R21s, A_R22s,
     .     D,DTCS, RH0,G, xKappa1,xKappa2, TW,
     .     tilt,Pdyn,x,ym,z, bx_ym,by_ym,bz_ym)
      call sst19_field_vect(TSS,TSO,TSE,
     .     A_dipsh, A_eq,A_eq_P, A_eq_TCS,A_eq_P_TCS,
     .     A_R11am,A_R12am,A_R11sm,A_R12sm,
     .     A_R11a, A_R12a, A_R11s, A_R12s,
     .     A_R21am,A_R22am,A_R21sm,A_R22sm,
     .     A_R21a, A_R22a, A_R21s, A_R22s,
     .     D,DTCS, RH0,G, xKappa1,xKappa2, TW,
     .     tilt,Pdyn,x,y,zp, bx_zp,by_zp,bz_zp)
      call sst19_field_vect(TSS,TSO,TSE,
     .     A_dipsh, A_eq,A_eq_P, A_eq_TCS,A_eq_P_TCS,
     .     A_R11am,A_R12am,A_R11sm,A_R12sm,
     .     A_R11a, A_R12a, A_R11s, A_R12s,
     .     A_R21am,A_R22am,A_R21sm,A_R22sm,
     .     A_R21a, A_R22a, A_R21s, A_R22s,
     .     D,DTCS, RH0,G, xKappa1,xKappa2, TW,
     .     tilt,Pdyn,x,y,zm, bx_zm,by_zm,bz_zm)

c     now compute the finite difference curl and convert to nA/m^2
      jx = Acurr*( (bz_yp-bz_ym) - (by_zp-by_zm) )
      jy = Acurr*( (bx_zp-bx_zm) - (bz_xp-bz_xm) )
      jz = Acurr*( (by_xp-by_xm) - (bx_yp-bx_ym) )
      
      return
      
      end

      
      subroutine sst19_field_vect(TSS,TSO,TSE,
     .     A_dipsh, A_eq,A_eq_P, A_eq_TCS,A_eq_P_TCS,
     .     A_R11am,A_R12am,A_R11sm,A_R12sm,
     .     A_R11a, A_R12a, A_R11s, A_R12s,
     .     A_R21am,A_R22am,A_R21sm,A_R22sm,
     .     A_R21a, A_R22a, A_R21s, A_R22s,
     .     D,DTCS, RH0,G, xKappa1,xKappa2, TW,
     .     tilt,Pdyn,X,Y,Z,
     .     BX,BY,BZ)

      IMPLICIT  NONE

C     Constants
      INTEGER    Mend,Nend
      PARAMETER (Mend = 6)
      PARAMETER (Nend = 8)

C     Inputs      
      REAL*8 TSS(80,Nend)
      REAL*8 TSO(80,Mend,Nend)
      REAL*8 TSE(80,Mend,Nend)
      REAL*8 A_dipsh
      REAL*8 A_eq(Nend+2*Mend*Nend)
      REAL*8 A_eq_P(Nend+2*Mend*Nend)
      REAL*8 A_eq_TCS(Nend+2*Mend*Nend)
      REAL*8 A_eq_P_TCS(Nend+2*Mend*Nend)
      REAL*8 A_R11am,A_R12am,A_R11sm,A_R12sm
      REAL*8 A_R11a, A_R12a, A_R11s, A_R12s
      REAL*8 A_R21am,A_R22am,A_R21sm,A_R22sm
      REAL*8 A_R21a, A_R22a, A_R21s, A_R22s
      REAL*8 D,DTCS             ! TAIL SHEET THICKNESS
      REAL*8 RH0,G
      REAL*8 xKappa1,xKappa2    !  SCALING FACTORS FOR BIRKELAND CURRENTS
      REAL*8 TW
      REAL*8 tilt,Pdyn
      REAL*8 X,Y,Z

C     Outputs
      REAL*8  BX,BY,BZ
      
c     Internal variables
      REAL*8  BXCF,BYCF,BZCF
      REAL*8  BXTS1(Nend),BXTO1(Mend,Nend),BXTE1(Mend,Nend)
      REAL*8  BYTS1(Nend),BYTO1(Mend,Nend),BYTE1(Mend,Nend)
      REAL*8  BZTS1(Nend),BZTO1(Mend,Nend),BZTE1(Mend,Nend)
      REAL*8  BXTS2(Nend),BXTO2(Mend,Nend),BXTE2(Mend,Nend)
      REAL*8  BYTS2(Nend),BYTO2(Mend,Nend),BYTE2(Mend,Nend)
      REAL*8  BZTS2(Nend),BZTO2(Mend,Nend),BZTE2(Mend,Nend)
      REAL*8  BX11am,BY11am,BZ11am
      REAL*8  BX12am,BY12am,BZ12am
      REAL*8  BX11sm,BY11sm,BZ11sm
      REAL*8  BX12sm,BY12sm,BZ12sm
      REAL*8  BX11a,BY11a,BZ11a
      REAL*8  BX12a,BY12a,BZ12a
      REAL*8  BX11s,BY11s,BZ11s
      REAL*8  BX12s,BY12s,BZ12s
      REAL*8  BX21am,BY21am,BZ21am
      REAL*8  BX22am,BY22am,BZ22am
      REAL*8  BX21sm,BY21sm,BZ21sm
      REAL*8  BX22sm,BY22sm,BZ22sm
      REAL*8  BX21a,BY21a,BZ21a
      REAL*8  BX22a,BY22a,BZ22a
      REAL*8  BX21s,BY21s,BZ21s
      REAL*8  BX22s,BY22s,BZ22s

      call sst19_field_expansion(TSS,TSO,TSE,
     .     A_dipsh, A_eq,A_eq_P, A_eq_TCS,A_eq_P_TCS,
     .     A_R11am,A_R12am,A_R11sm,A_R12sm,
     .     A_R11a, A_R12a, A_R11s, A_R12s,
     .     A_R21am,A_R22am,A_R21sm,A_R22sm,
     .     A_R21a, A_R22a, A_R21s, A_R22s,
     .     D,DTCS, RH0,G, xKappa1,xKappa2, TW,
     .     tilt,Pdyn,X,Y,Z,
     .     BXCF,BYCF,BZCF,
     .     BXTS1,BYTS1,BZTS1,BXTO1,BYTO1,BZTO1,BXTE1,BYTE1,BZTE1,
     .     BXTS2,BYTS2,BZTS2,BXTO2,BYTO2,BZTO2,BXTE2,BYTE2,BZTE2,
     .     BX11am,BY11am,BZ11am,
     .     BX12am,BY12am,BZ12am,
     .     BX11sm,BY11sm,BZ11sm,
     .     BX12sm,BY12sm,BZ12sm,
     .     BX11a,BY11a,BZ11a,
     .     BX12a,BY12a,BZ12a,
     .     BX11s,BY11s,BZ11s,
     .     BX12s,BY12s,BZ12s,
     .     BX21am,BY21am,BZ21am,
     .     BX22am,BY22am,BZ22am,      
     .     BX21sm,BY21sm,BZ21sm,
     .     BX22sm,BY22sm,BZ22sm,      
     .     BX21a,BY21a,BZ21a,
     .     BX22a,BY22a,BZ22a,
     .     BX21s,BY21s,BZ21s,
     .     BX22s,BY22s,BZ22s,
     .     BX,BY,BZ)

      return
      end

      
      subroutine sst19_field_expansion(TSS,TSO,TSE,
     .     A_dipsh, A_eq,A_eq_P, A_eq_TCS,A_eq_P_TCS,
     .     A_R11am,A_R12am,A_R11sm,A_R12sm,
     .     A_R11a, A_R12a, A_R11s, A_R12s,
     .     A_R21am,A_R22am,A_R21sm,A_R22sm,
     .     A_R21a, A_R22a, A_R21s, A_R22s,
     .     D,DTCS, RH0,G, xKappa1,xKappa2, TW,
     .     tilt,Pdyn,X,Y,Z,
     .     BXCF,BYCF,BZCF,
     .     BXTS1,BYTS1,BZTS1,BXTO1,BYTO1,BZTO1,BXTE1,BYTE1,BZTE1,
     .     BXTS2,BYTS2,BZTS2,BXTO2,BYTO2,BZTO2,BXTE2,BYTE2,BZTE2,
     .     BX11am,BY11am,BZ11am,
     .     BX12am,BY12am,BZ12am,
     .     BX11sm,BY11sm,BZ11sm,
     .     BX12sm,BY12sm,BZ12sm,
     .     BX11a,BY11a,BZ11a,
     .     BX12a,BY12a,BZ12a,
     .     BX11s,BY11s,BZ11s,
     .     BX12s,BY12s,BZ12s,
     .     BX21am,BY21am,BZ21am,
     .     BX22am,BY22am,BZ22am,      
     .     BX21sm,BY21sm,BZ21sm,
     .     BX22sm,BY22sm,BZ22sm,      
     .     BX21a,BY21a,BZ21a,
     .     BX22a,BY22a,BZ22a,
     .     BX21s,BY21s,BZ21s,
     .     BX22s,BY22s,BZ22s,
     .     BX,BY,BZ)

      IMPLICIT  NONE

C     Constants
      INTEGER    Mend,Nend
      PARAMETER (Mend = 6)
      PARAMETER (Nend = 8)

C     Inputs      
      REAL*8 TSS(80,Nend)
      REAL*8 TSO(80,Mend,Nend)
      REAL*8 TSE(80,Mend,Nend)
      REAL*8 A_dipsh
      REAL*8 A_eq(Nend+2*Mend*Nend)
      REAL*8 A_eq_P(Nend+2*Mend*Nend)
      REAL*8 A_eq_TCS(Nend+2*Mend*Nend)
      REAL*8 A_eq_P_TCS(Nend+2*Mend*Nend)
      REAL*8 A_R11am,A_R12am,A_R11sm,A_R12sm
      REAL*8 A_R11a, A_R12a, A_R11s, A_R12s
      REAL*8 A_R21am,A_R22am,A_R21sm,A_R22sm
      REAL*8 A_R21a, A_R22a, A_R21s, A_R22s
      REAL*8 D,DTCS             ! TAIL SHEET THICKNESS
      REAL*8 RH0,G
      REAL*8 xKappa1,xKappa2    !  SCALING FACTORS FOR BIRKELAND CURRENTS
      REAL*8 TW
      REAL*8 tilt,Pdyn
      REAL*8 X,Y,Z

C     Outputs
      REAL*8  BXCF,BYCF,BZCF
      REAL*8  BXTS1(Nend),BXTO1(Mend,Nend),BXTE1(Mend,Nend)
      REAL*8  BYTS1(Nend),BYTO1(Mend,Nend),BYTE1(Mend,Nend)
      REAL*8  BZTS1(Nend),BZTO1(Mend,Nend),BZTE1(Mend,Nend)
      REAL*8  BXTS2(Nend),BXTO2(Mend,Nend),BXTE2(Mend,Nend)
      REAL*8  BYTS2(Nend),BYTO2(Mend,Nend),BYTE2(Mend,Nend)
      REAL*8  BZTS2(Nend),BZTO2(Mend,Nend),BZTE2(Mend,Nend)
      REAL*8  BX11am,BY11am,BZ11am
      REAL*8  BX12am,BY12am,BZ12am
      REAL*8  BX11sm,BY11sm,BZ11sm
      REAL*8  BX12sm,BY12sm,BZ12sm
      REAL*8  BX11a,BY11a,BZ11a
      REAL*8  BX12a,BY12a,BZ12a
      REAL*8  BX11s,BY11s,BZ11s
      REAL*8  BX12s,BY12s,BZ12s
      REAL*8  BX21am,BY21am,BZ21am
      REAL*8  BX22am,BY22am,BZ22am
      REAL*8  BX21sm,BY21sm,BZ21sm
      REAL*8  BX22sm,BY22sm,BZ22sm
      REAL*8  BX21a,BY21a,BZ21a
      REAL*8  BX22a,BY22a,BZ22a
      REAL*8  BX21s,BY21s,BZ21s
      REAL*8  BX22s,BY22s,BZ22s
      REAL*8  BX,BY,BZ
      
c     Internal variables
C      
      REAL*8 XX,YY,ZZ
      REAL*8 CFX,CFY,CFZ
      REAL*8 TX,TY,TZ
      REAL*8 bFacX,bFacY,bFacZ
C
      REAL*8 Pdyn_0,P_factor,Xappa,Xappa3,coeff
C
      INTEGER m,n,ind,indS
C      

      Pdyn_0=2.d0               !   AVERAGE PRESSURE USED FOR NORMALIZATION
      
      Xappa=(Pdyn/Pdyn_0)**0.155   !   0.155 is the value obtained in TS05
      Xappa3=Xappa**3
      
C     Rescale the position vector using the pressure scaling
      XX=X*Xappa  ! pressure scaling has been reinstated here
      YY=Y*Xappa
      ZZ=Z*Xappa
C     Evaluate the shielded dipole field
      call dipoleshield(tilt, XX,YY,ZZ, CFX,CFY,CFZ) !  DIPOLE SHIELDING FIELD
      BXCF=CFX  *Xappa3
      BYCF=CFY  *Xappa3
      BZCF=CFZ  *Xappa3
C     Evaluate the deformed equatorial expansion field
      CALL deformedsheet_with_tcs (Mend,Nend, D,DTCS, RH0,G,TW,
     .     TSS,TSO,TSE,         !  TAIL FIELD (THREE MODES)
     .     tilt,XX,YY,ZZ,      
     .     BXTS1,BYTS1,BZTS1, BXTO1,BYTO1,BZTO1, BXTE1,BYTE1,BZTE1,
     .     BXTS2,BYTS2,BZTS2, BXTO2,BYTO2,BZTO2, BXTE2,BYTE2,BZTE2)
C     Evaluate the field-aligned current field, note only 4 of the 8 FAC expansions are used in TS07
      CALL fac_total_sst19 (xKappa1,xKappa2,
     .     tilt,XX,YY,ZZ,
     .     BX11am,BY11am,BZ11am,
     .     BX12am,BY12am,BZ12am,
     .     BX11sm,BY11sm,BZ11sm,
     .     BX12sm,BY12sm,BZ12sm,
     .     BX11a,BY11a,BZ11a,
     .     BX12a,BY12a,BZ12a,
     .     BX11s,BY11s,BZ11s,
     .     BX12s,BY12s,BZ12s,
     .     BX21am,BY21am,BZ21am,
     .     BX22am,BY22am,BZ22am,      
     .     BX21sm,BY21sm,BZ21sm,
     .     BX22sm,BY22sm,BZ22sm,      
     .     BX21a,BY21a,BZ21a,
     .     BX22a,BY22a,BZ22a,
     .     BX21s,BY21s,BZ21s,
     .     BX22s,BY22s,BZ22s) ! (TWO MODES FOR R1s AND TWO MODES FOR R2s)
c                                                           (but we actually use from here only R2s modes)
C
C-----------------------------------------------------------
C
C    NOW, ADD UP ALL THE COMPONENTS:
C
      TX=0.d0
      TY=0.d0
      TZ=0.d0

C --- New tail structure -------------

      P_factor=DSQRT(Pdyn/Pdyn_0)-1.d0

      ind=1
      
      DO n=1,Nend
         
         coeff=A_eq(ind)+A_eq_P(ind)*P_factor

         TX=TX+coeff*BXTS1(n)    !   2 - 6  &  47 - 51
         TY=TY+coeff*BYTS1(n)
         TZ=TZ+coeff*BZTS1(n)

         ind=ind+1
      ENDDO
      
      DO n=1,Nend
         DO  m=1,Mend

            coeff=A_eq(ind)+A_eq_P(ind)*P_factor
            
            TX=TX+coeff*BXTO1(m,n) !   7 -26  &  52 - 71
            TY=TY+coeff*BYTO1(m,n)
            TZ=TZ+coeff*BZTO1(m,n)
            
            indS=ind+Mend*Nend
            coeff=A_eq(indS)+A_eq_P(indS)*P_factor
            
            TX=TX+coeff*BXTE1(m,n) !   27 -46  &  72 - 91
            TY=TY+coeff*BYTE1(m,n)
            TZ=TZ+coeff*BZTE1(m,n)
            
            ind=ind+1
           
         ENDDO
      ENDDO

      ind=1
      
      DO n=1,Nend
         
         coeff=A_eq_TCS(ind)+A_eq_P_TCS(ind)*P_factor

         TX=TX+coeff*BXTS2(n)    !   2 - 6  &  47 - 51
         TY=TY+coeff*BYTS2(n)
         TZ=TZ+coeff*BZTS2(n)

         ind=ind+1
      ENDDO
      
      DO n=1,Nend
         DO  m=1,Mend

            coeff=A_eq_TCS(ind)+A_eq_P_TCS(ind)*P_factor
            
            TX=TX+coeff*BXTO2(m,n) !   7 -26  &  52 - 71
            TY=TY+coeff*BYTO2(m,n)
            TZ=TZ+coeff*BZTO2(m,n)
            
            indS=ind+Mend*Nend
            coeff=A_eq_TCS(indS)+A_eq_P_TCS(indS)*P_factor
            
            TX=TX+coeff*BXTE2(m,n) !   27 -46  &  72 - 91
            TY=TY+coeff*BYTE2(m,n)
            TZ=TZ+coeff*BZTE2(m,n)
            
            ind=ind+1
           
         ENDDO
      ENDDO

      bFacX=A_R11am*BX11am+A_R12am*BX12am+A_R11sm*BX11sm+A_R12sm*BX12sm+
     .      A_R11a *BX11a +A_R12a *BX12a +A_R11s *BX11s +A_R12s *BX12s +
     .      A_R21am*BX21am+A_R22am*BX22am+A_R21sm*BX21sm+A_R22sm*BX22sm+
     .      A_R21a *BX21a +A_R22a *BX22a +A_R21s *BX21s +A_R22s *BX22s
      bFacY=A_R11am*BY11am+A_R12am*BY12am+A_R11sm*BY11sm+A_R12sm*BY12sm+
     .      A_R11a *BY11a +A_R12a *BY12a +A_R11s *BY11s +A_R12s *BY12s +
     .      A_R21am*BY21am+A_R22am*BY22am+A_R21sm*BY21sm+A_R22sm*BY22sm+
     .      A_R21a *BY21a +A_R22a *BY22a +A_R21s *BY21s +A_R22s *BY22s
      bFacZ=A_R11am*BZ11am+A_R12am*BZ12am+A_R11sm*BZ11sm+A_R12sm*BZ12sm+
     .      A_R11a *BZ11a +A_R12a *BZ12a +A_R11s *BZ11s +A_R12s *BZ12s +
     .      A_R21am*BZ21am+A_R22am*BZ22am+A_R21sm*BZ21sm+A_R22sm*BZ22sm+
     .      A_R21a *BZ21a +A_R22a *BZ22a +A_R21s *BZ21s +A_R22s *BZ22s
      
C
c   -----------------------------------------------------------
C
C   AND WE HAVE THE TOTAL EXTERNAL FIELD.
C      
      BX=A_dipsh*BXCF+TX+bFacX      
      BY=A_dipsh*BYCF+TY+bFacY
      BZ=A_dipsh*BZCF+TZ+bFacZ

      return
      end
