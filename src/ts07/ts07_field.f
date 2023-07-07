      subroutine store_ts07_static(staticDir)

      IMPLICIT NONE
      
C     Inputs      
      CHARACTER*256 staticDir
      CHARACTER*256 variableDir
      INTEGER iYear,iDayOfYear,iHour,min
      CHARACTER*256 filename

C     Constants
      INTEGER    Mend,Nend
      PARAMETER (Mend = 4)
      PARAMETER (Nend = 5)
      
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
      REAL*8 A_R11a,A_R12a,A_R21a,A_R21s
      REAL*8 D,RH0,G
      REAL*8 xKappa1,xKappa2
      REAL*8 TW
      REAL*8 Pdyn
      REAL*8 tilt
      
c     saved
      save TSS,TSO,TSE
      save A_dipsh
      save A_eq, A_eq_P
      save A_R11a,A_R12a,A_R21a,A_R21s
      save D,RH0,G
      save xKappa1,xKappa2
      save TW
      save Pdyn
      save tilt
      
C     Inputs
      REAL*8 tiltIn,PdynIn
      REAL*8 x,y,z

c     Outputs
      REAL*8 bx,by,bz      

c     Internal variables
      REAL*8  BXCF,BYCF,BZCF
      REAL*8  BXTS(Nend),BXTO(Mend,Nend),BXTE(Mend,Nend)
      REAL*8  BYTS(Nend),BYTO(Mend,Nend),BYTE(Mend,Nend)
      REAL*8  BZTS(Nend),BZTO(Mend,Nend),BZTE(Mend,Nend)
      REAL*8  BXR11a,BYR11a,BZR11a, BXR12a,BYR12a,BZR12a
      REAL*8  BXR21a,BYR21a,BZR21a, BXR21s,BYR21s,BZR21s
      
      call read_ts07_static(staticDir,
     .     TSS,TSO,TSE)
      init_static=.true.      
      return

      entry store_ts07_variable(filename)
      call read_ts07_variable(filename,
     .     A_dipsh, A_eq,A_eq_P, A_R11a,A_R12a,A_R21a,A_R21s,
     .     D,RH0,G, xKappa1,xKappa2, TW,
     .     Pdyn,tilt)
      init_variable=.true.
      return
      

      entry store_ts07_variable_datetime(variableDir,
     .     iYear,iDayOfYear,iHour,min)
      call read_ts07_variable_datetime(variableDir,
     .     iYear,iDayOfYear,iHour,min,
     .     A_dipsh, A_eq,A_eq_P, A_R11a,A_R12a,A_R21a,A_R21s,
     .     D,RH0,G, xKappa1,xKappa2, TW,
     .     Pdyn,tilt)    
      init_variable=.true.
      return

      
      entry ts07_field_tilt(tiltIn,x,y,z, bx,by,bz)
      
c     check that the static coefficients were loaded
      if (init_static .eqv. .false.) then
         stop 'the static shielding coefficients must be '//
     .        'stored first by calling store_ts07_static'
      endif      

c     check that the time-variable variables were loaded
      if (init_variable .eqv. .false.) then
         stop 'the time-variable parameters and coefficients '//
     .        'must be stored first by calling store_ts07_variable'
      endif
      
      call ts07_field_expansion(TSS,TSO,TSE,
     .     A_dipsh, A_eq,A_eq_P,
     .     A_R11a,A_R12a,A_R21a,A_R21s,
     .     D, RH0,G, xKappa1,xKappa2, TW,
     .     tiltIn,Pdyn,X,Y,Z,
     .     BXCF,BYCF,BZCF,
     .     BXTS,BYTS,BZTS,BXTO,BYTO,BZTO,BXTE,BYTE,BZTE,
     .     BXR11a,BYR11a,BZR11a,BXR12a,BYR12a,BZR12a,
     .     BXR21a,BYR21a,BZR21a,BXR21s,BYR21s,BZR21s,
     .     BX,BY,BZ)
      
      return
      
      
      entry ts07_field(x,y,z, bx,by,bz)
      
c     check that the static coefficients were loaded
      if (init_static .eqv. .false.) then
         stop 'the static shielding coefficients must be '//
     .        'stored first by calling store_ts07_static'
      endif      

c     check that the time-variable variables were loaded
      if (init_variable .eqv. .false.) then
         stop 'the time-variable parameters and coefficients '//
     .        'must be stored first by calling store_ts07_variable'
      endif
      
      call ts07_field_expansion(TSS,TSO,TSE,
     .     A_dipsh, A_eq,A_eq_P,
     .     A_R11a,A_R12a,A_R21a,A_R21s,
     .     D, RH0,G, xKappa1,xKappa2, TW,
     .     tilt,Pdyn,X,Y,Z,
     .     BXCF,BYCF,BZCF,
     .     BXTS,BYTS,BZTS,BXTO,BYTO,BZTO,BXTE,BYTE,BZTE,
     .     BXR11a,BYR11a,BZR11a,BXR12a,BYR12a,BZR12a,
     .     BXR21a,BYR21a,BZR21a,BXR21s,BYR21s,BZR21s,
     .     BX,BY,BZ)
      
      return
      
      end
      

      subroutine ts07_field_expansion(TSS,TSO,TSE,
     .     A_dipsh, A_eq,A_eq_P,
     .     A_R11a,A_R12a,A_R21a,A_R21s,
     .     D, RH0,G, xKappa1,xKappa2, TW,
     .     tilt,Pdyn,X,Y,Z,
     .     BXCF,BYCF,BZCF,
     .     BXTS,BYTS,BZTS,BXTO,BYTO,BZTO,BXTE,BYTE,BZTE,
     .     BXR11a,BYR11a,BZR11a,BXR12a,BYR12a,BZR12a,
     .     BXR21a,BYR21a,BZR21a,BXR21s,BYR21s,BZR21s,
     .     BX,BY,BZ)
C
C   IOPGEN - GENERAL OPTION FLAG:  IOPGEN=0 - CALCULATE TOTAL FIELD
C                                  IOPGEN=1 - DIPOLE SHIELDING ONLY
C                                  IOPGEN=2 - TAIL FIELD ONLY
C                                  IOPGEN=3 - BIRKELAND FIELD ONLY
C                                  IOPGEN=4 - RING CURRENT FIELD ONLY
C      
      IMPLICIT  NONE

C     Constants
      INTEGER    Mend,Nend
      PARAMETER (Mend = 4)
      PARAMETER (Nend = 5)

C     Inputs      
      REAL*8 TSS(80,Nend)
      REAL*8 TSO(80,Mend,Nend)
      REAL*8 TSE(80,Mend,Nend)
      REAL*8 A_dipsh
      REAL*8 A_eq(Nend+2*Mend*Nend)
      REAL*8 A_eq_P(Nend+2*Mend*Nend)
      REAL*8 A_R11a,A_R12a,A_R21a,A_R21s
      REAL*8 D                  ! TAIL SHEET THICKNESS
      REAL*8 RH0,G
      REAL*8 xKappa1,xKappa2    !  SCALING FACTORS FOR BIRKELAND CURRENTS
      REAL*8 TW
      REAL*8 tilt,Pdyn
      REAL*8 X,Y,Z

C     Outputs
      REAL*8  BXCF,BYCF,BZCF
      REAL*8  BXTS(Nend),BXTO(Mend,Nend),BXTE(Mend,Nend)
      REAL*8  BYTS(Nend),BYTO(Mend,Nend),BYTE(Mend,Nend)
      REAL*8  BZTS(Nend),BZTO(Mend,Nend),BZTE(Mend,Nend)
      REAL*8  BXR11a,BYR11a,BZR11a, BXR12a,BYR12a,BZR12a
      REAL*8  BXR21a,BYR21a,BZR21a, BXR21s,BYR21s,BZR21s
      REAL*8  BX,BY,BZ
      
c     Internal variables
C
      REAL*8  BXR22a,BYR22a,BZR22a
      REAL*8  BXR11s,BYR11s,BZR11s, BXR12s,BYR12s,BZR12s
      REAL*8  BXR22s,BYR22s,BZR22s
C      
      REAL*8 XX,YY,ZZ
      REAL*8 CFX,CFY,CFZ
      REAL*8 TX,TY,TZ
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
      CALL deformedsheet (Mend,Nend, D, RH0,G,TW, TSS,TSO,TSE, !  TAIL FIELD (THREE MODES)
     .     tilt,XX,YY,ZZ,      
     .     BXTS,BYTS,BZTS, BXTO,BYTO,BZTO, BXTE,BYTE,BZTE)
C     Evaluate the field-aligned current field, note only 4 of the 8 FAC expansions are used in TS07
      CALL fac_total_ts07 (xKappa1,xKappa2,
     .     tilt,XX,YY,ZZ,
     .     BXR11a,BYR11a,BZR11a,BXR12a,BYR12a,BZR12a, ! BIRKELAND FIELD (TWO MODES FOR R1 AND TWO MODES FOR R2)
     .     BXR21a,BYR21a,BZR21a,BXR22a,BYR22a,BZR22a, 
     .     BXR11s,BYR11s,BZR11s,BXR12s,BYR12s,BZR12s, ! "SYMMETRIC" BIRKELAND FIELD
     .     BXR21s,BYR21s,BZR21s,BXR22s,BYR22s,BZR22s) ! (TWO MODES FOR R1s AND TWO MODES FOR R2s)
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

         TX=TX+coeff*BXTS(n)    !   2 - 6  &  47 - 51
         TY=TY+coeff*BYTS(n)
         TZ=TZ+coeff*BZTS(n)

         ind=ind+1
      ENDDO
      
      DO n=1,Nend
         DO  m=1,Mend

            coeff=A_eq(ind)+A_eq_P(ind)*P_factor
            
            TX=TX+coeff*BXTO(m,n) !   7 -26  &  52 - 71
            TY=TY+coeff*BYTO(m,n)
            TZ=TZ+coeff*BZTO(m,n)
            
            indS=ind+Mend*Nend
            coeff=A_eq(indS)+A_eq_P(indS)*P_factor
            
            TX=TX+coeff*BXTE(m,n) !   27 -46  &  72 - 91
            TY=TY+coeff*BYTE(m,n)
            TZ=TZ+coeff*BZTE(m,n)
            
            ind=ind+1
           
         ENDDO
      ENDDO

C
c   -----------------------------------------------------------
C
C   AND WE HAVE THE TOTAL EXTERNAL FIELD.
C      
      BX=A_dipsh*BXCF+TX+
     *     A_R11a*BXR11a+A_R12a*BXR12a+A_R21a*BXR21a+A_R21s*BXR21s
      
      BY=A_dipsh*BYCF+TY+
     *     A_R11a*BYR11a+A_R12a*BYR12a+A_R21a*BYR21a+A_R21s*BYR21s
      
      BZ=A_dipsh*BZCF+TZ+
     *     A_R11a*BZR11a+A_R12a*BZR12a+A_R21a*BZR21a+A_R21s*BZR21s
      
      return
      end
