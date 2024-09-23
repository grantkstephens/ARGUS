!
!© 2023 The Johns Hopkins University Applied Physics Laboratory LLC
!
      subroutine read_sst19_static(staticDir,
     .     TSS,TSO,TSE)

      IMPLICIT NONE
      
C     Inputs      
      CHARACTER*256 staticDir

C     Constants
      INTEGER    Mend,Nend
      PARAMETER (Mend = 6)
      PARAMETER (Nend = 8)
      
C     Outputs
      REAL*8 TSS(80,Nend)
      REAL*8 TSO(80,Mend,Nend)
      REAL*8 TSE(80,Mend,Nend)
      
C     Internal variables
      INTEGER IREAD,KREAD,KK
      CHARACTER*256 FILENAME
      
      
c--------------------------------------------------------------------------
c//////////////////////////////////////////////////////////////////////////
c     The following reads the static coefficients that are used for the Equatorial Field
c     shielding fields. These only need to be loaded once as they are common for all times.
      DO IREAD=1,Nend
        WRITE(filename,'(A,A,I0.2,A)')
     *        TRIM(staticDir),'/tailamebhr_',IREAD,'.par'
      OPEN (UNIT=1,FILE=FILENAME)
      READ (1,200) (TSS(KK,IREAD),KK=1,80)
 200  FORMAT(G17.10)
      ENDDO
      CLOSE(1)

      DO IREAD=1,Nend
         DO KREAD=1,Mend
      WRITE(filename,'(A,A,I0.2,A,I0.2,A)')
     *        TRIM(staticDir),'/tailamhr_o_',IREAD,'_',KREAD,'.par'
      OPEN (UNIT=1,FILE=FILENAME)
      READ (1,200) (TSO(KK,KREAD,IREAD),KK=1,80)
      ENDDO
      ENDDO
      CLOSE(1)

      DO IREAD=1,Nend
        DO KREAD=1,Mend
        WRITE(filename,'(A,A,I0.2,A,I0.2,A)')
     *        TRIM(staticDir),'/tailamhr_e_',IREAD,'_',KREAD,'.par'
      OPEN (UNIT=1,FILE=FILENAME)
      READ (1,200) (TSE(KK,KREAD,IREAD),KK=1,80)
      ENDDO
      ENDDO
      CLOSE(1)                                        

      return
      end

      
      subroutine read_sst19_variable_datetime(variableDir,
     .     iYear,iDayOfYear,iHour,min,
     .     A_dipsh, A_eq,A_eq_P, A_eq_TCS,A_eq_P_TCS,
     .     A_R11am,A_R12am,A_R11sm,A_R12sm,
     .     A_R11a, A_R12a, A_R11s, A_R12s,
     .     A_R21am,A_R22am,A_R21sm,A_R22sm,
     .     A_R21a, A_R22a, A_R21s, A_R22s,
     .     D,DTCS, RH0,G, xKappa1,xKappa2, TW,
     .     Pdyn,tilt)

      IMPLICIT NONE

C     Inputs      
      CHARACTER*256 variableDir
      INTEGER iYear,iDayOfYear,iHour,min

C     Constants
      INTEGER   Mend,Nend
      PARAMETER    (Mend = 6)
      PARAMETER    (Nend = 8)      
      
C     Outputs
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
      REAL*8  Pdyn
      REAL*8  tilt

C     Internal variables      
      CHARACTER*256 filename
      
c      WRITE(filename,'(A,I0.4,A,I0.3,A,I0.4,A,I0.3,A,I0.2,A,I0.2,A)')
c     *     TRIM(variableDir),
c     *     iYear,'_',iDayOfYear,
c     .     '/',iYear,'_',iDayOfYear,'_',iHour,'_',min,'.par'
      WRITE(filename,'(A,A,I0.4,A,I0.3,A,I0.2,A,I0.2,A)')
     .     TRIM(variableDir),'/',
     .     iYear,'_',iDayOfYear,'_',iHour,'_',min,'.par'
      
      
      call read_sst19_variable(filename,
     .     A_dipsh, A_eq,A_eq_P,A_eq_TCS,A_eq_P_TCS,
     .     A_R11am,A_R12am,A_R11sm,A_R12sm,
     .     A_R11a, A_R12a, A_R11s, A_R12s,
     .     A_R21am,A_R22am,A_R21sm,A_R22sm,
     .     A_R21a, A_R22a, A_R21s, A_R22s,
     .     D,DTCS, RH0,G, xKappa1,xKappa2, TW,
     .     Pdyn,tilt)
      
      return
      end
      
      
      subroutine read_sst19_variable(filename,
     .     A_dipsh, A_eq,A_eq_P,A_eq_TCS,A_eq_P_TCS,
     .     A_R11am,A_R12am,A_R11sm,A_R12sm,
     .     A_R11a, A_R12a, A_R11s, A_R12s,
     .     A_R21am,A_R22am,A_R21sm,A_R22sm,
     .     A_R21a, A_R22a, A_R21s, A_R22s,
     .     D,DTCS, RH0,G, xKappa1,xKappa2, TW,
     .     Pdyn,tilt)

      IMPLICIT NONE

C     Inputs      
      CHARACTER*256 filename

C     Constants
      INTEGER   Mend,Nend
      PARAMETER    (Mend = 6)
      PARAMETER    (Nend = 8)     
      
C     Outputs
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
      REAL*8  Pdyn
      REAL*8  tilt
      
C     Internal variables      
      INTEGER    Ntot
      PARAMETER (Ntot=440)
      REAL*8 params(NTOT)
      INTEGER IREAD
      
      REAL*8  QValue
      REAL*8  Brms
      INTEGER M_in
      INTEGER N_in
      
c     Now read variable parameter and coefficients file, the variable parameters and coefficients 
c     are stored in the array PARAMS. Additionaly, the dynamic pressure and the dipole tilt angle
c     are stored in PDYN and TILT respectively. All of these are variable or time-dependent values.
c     The other values stored in the file are not currently being used.
      OPEN (UNIT=1,FILE=filename,action='read') ! open the filed
      READ (1,100) (PARAMS(IREAD),IREAD=1,NTOT) ! read the variable coefficients and parameters
      READ (1,101) QValue ! the Q factor, related to chi squared, a measure of the goodness of fit
      READ (1,101) Brms 
      READ (1,102) M_in ! the number of azimuthal expansions in the equatorial field module (M=4)
      READ (1,102) N_in ! the number of radial expansions in the equatorial field module (N=5)
      READ (1,101) Pdyn ! the dynamic pressure for this time
      READ (1,101) tilt ! the dipole tilt for this time
 100  FORMAT(G15.6)                                         
 101  FORMAT(7x,G15.6)                                            
 102  FORMAT(7x,I15)                                            
      CLOSE(1)
      
      A_dipsh=params(1)

      A_eq=  params(2:105)
      A_eq_P=params(106:209)
      A_eq_TCS=  params(210:313)
      A_eq_P_TCS=params(314:417)
      
      A_R11am=params(418)
      A_R12am=params(419)
      A_R11sm=params(420)
      A_R12sm=params(421)

      A_R11a=params(422)
      A_R12a=params(423)
      A_R11s=params(424)
      A_R12s=params(425)
      
      A_R21am=params(426)
      A_R22am=params(427)
      A_R21sm=params(428)
      A_R22sm=params(429)
      
      A_R21a=params(430)
      A_R22a=params(431)
      A_R21s=params(432)
      A_R22s=params(433)
      
      D=      params(434)
      DTCS=   params(435)
      RH0=    params(436)
      G=      params(437)
      xKappa1=params(438)
      xKappa2=params(439)
      TW=     params(440) !   THIS PARAMETER CONTROLS THE IMF-INDUCED TWISTING (ADDED 04/21/06)

      return
      end
