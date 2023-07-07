      subroutine read_ts07_static(staticDir,
     .     TSS,TSO,TSE)

      IMPLICIT NONE
      
C     Inputs      
      CHARACTER*256 staticDir

C     Constants
      INTEGER    Mend,Nend
      PARAMETER (Mend = 4)
      PARAMETER (Nend = 5)
      
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
        WRITE(filename,'(A,A,I0,A)')
     *        TRIM(staticDir),'/tailamebhr',IREAD,'.par'
      OPEN (UNIT=1,FILE=FILENAME)
      READ (1,200) (TSS(KK,IREAD),KK=1,80)
 200  FORMAT(G17.10)
      ENDDO
      CLOSE(1)

      DO IREAD=1,Nend
         DO KREAD=1,Mend
      WRITE(filename,'(A,A,I0,I0,A)')
     *        TRIM(staticDir),'/tailamhr_o_',IREAD,KREAD,'.par'
      OPEN (UNIT=1,FILE=FILENAME)
      READ (1,200) (TSO(KK,KREAD,IREAD),KK=1,80)
      ENDDO
      ENDDO
      CLOSE(1)

      DO IREAD=1,Nend
        DO KREAD=1,Mend
        WRITE(filename,'(A,A,I0,I0,A)')
     *        TRIM(staticDir),'/tailamhr_e_',IREAD,KREAD,'.par'
      OPEN (UNIT=1,FILE=FILENAME)
      READ (1,200) (TSE(KK,KREAD,IREAD),KK=1,80)
      ENDDO
      ENDDO
      CLOSE(1)                                        

      return
      end

      
      subroutine read_ts07_variable_datetime(variableDir,
     .     iYear,iDayOfYear,iHour,min,
     .     A_dipsh, A_eq,A_eq_P, A_R11a,A_R12a,A_R21a,A_R21s,
     .     D,RH0,G, xKappa1,xKappa2, TW,
     .     Pdyn,tilt)

      IMPLICIT NONE

C     Inputs      
      CHARACTER*256 variableDir
      INTEGER iYear,iDayOfYear,iHour,min

C     Constants
      INTEGER   Mend,Nend
      PARAMETER    (Mend = 4)
      PARAMETER    (Nend = 5)      
      
C     Outputs
      REAL*8 A_dipsh
      REAL*8 A_eq(Nend+2*Mend*Nend)
      REAL*8 A_eq_P(Nend+2*Mend*Nend)
      REAL*8 A_R11a,A_R12a,A_R21a,A_R21s
      REAL*8 D,RH0,G
      REAL*8 xKappa1,xKappa2
      REAL*8 TW
      REAL*8  Pdyn
      REAL*8  tilt

C     Internal variables      
      CHARACTER*256 filename
      
      WRITE(filename,'(A,I0.4,A,I0.3,A,I0.4,A,I0.3,A,I0.2,A,I0.2,A)')
     *     TRIM(variableDir),
     *     iYear,'_',iDayOfYear,
     .     '/',iYear,'_',iDayOfYear,'_',iHour,'_',min,'.par'
      
      call read_ts07_variable(filename,
     .     A_dipsh, A_eq,A_eq_P, A_R11a,A_R12a,A_R21a,A_R21s,
     .     D,RH0,G, xKappa1,xKappa2, TW,
     .     Pdyn,tilt)
      
      return
      end
      
      
      subroutine read_ts07_variable(filename,
     .     A_dipsh, A_eq,A_eq_P, A_R11a,A_R12a,A_R21a,A_R21s,
     .     D,RH0,G, xKappa1,xKappa2, TW,
     .     Pdyn,tilt)

      IMPLICIT NONE

C     Inputs      
      CHARACTER*256 filename

C     Constants
      INTEGER   Mend,Nend
      PARAMETER    (Mend = 4)
      PARAMETER    (Nend = 5)
      
C     Outputs
      REAL*8 A_dipsh
      REAL*8 A_eq(Nend+2*Mend*Nend)
      REAL*8 A_eq_P(Nend+2*Mend*Nend)
      REAL*8 A_R11a,A_R12a,A_R21a,A_R21s
      REAL*8 D,RH0,G
      REAL*8 xKappa1,xKappa2
      REAL*8 TW
      REAL*8  Pdyn
      REAL*8  tilt
      
C     Internal variables      
      INTEGER    Ntot
      PARAMETER (Ntot=101)
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

      A_eq=  params(2:46)
      A_eq_P=params(47:91)
      
      A_R11a=params(92)
      A_R12a=params(93)
      A_R21a=params(94)
      A_R21s=params(95)

      D=      params(96)
      RH0=    params(97)
      G=      params(98)
      xKappa1=params(99)
      xKappa2=params(100)
      TW=     params(101) !   THIS PARAMETER CONTROLS THE IMF-INDUCED TWISTING (ADDED 04/21/06)
      
      return
      end
