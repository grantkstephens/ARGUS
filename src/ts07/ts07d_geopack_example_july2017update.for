!======================================================================================================
! July 2017: G.K.Stephens (Grant.Stephens@jhuapl.edu)
! This program demonstrates how to setup and evaluate the Tsyganenko and Sitnov [2007] (TS07D)
! empirical geomagnetic field model and how to setup and evaluate the IGRF-12 model included
! in Tsyganenko's GEOPACK package.
!
! This July 2017 update incorporates an upgraded version of the Bessel function evaluator,
! provided by Jay Albert, which significantly speeds up the model. We thank Jay Albert for these
! contributions. 
!
! Additionally, the model subroutine was updated to be a double precision subroutine.
! To indicate the update, the subroutine was renamed from EXTMODEL to TS07D_JULY_2017.
! 
! To compile this code with GFortran, use:
! gfortran ts07d_geopack_example_july2017update.for -o ts07d_geopack_example_july2017update
! and to run: ./ts07d_geopack_example_july2017update
!======================================================================================================
      P r o g r a m   TS07D
      implicit none 
c----------------------------------------------------------------------
c     Inputs to the model
      INTEGER   IOPT
      REAL*8  PARMOD(10), XGSM,YGSM,ZGSM      

c     Input parameters that will be passed through the common block
      INTEGER    NTOT
      PARAMETER (NTOT=101)
      REAL*8 PDYN,PARAMS(NTOT),TSS,TSO,TSE
      COMMON /PARAM/ PARAMS
      COMMON /INPUT/  PDYN
      COMMON /TSS/ TSS(80,5)
      COMMON /TSO/ TSO(80,5,4)
      COMMON /TSE/ TSE(80,5,4)

c     Set up the Geopack common block
      REAL*8 AAA,BBB,PSI ! the dipole tilt
      COMMON /GEOPACK1/ AAA(15),PSI,BBB(18)

c     Output parameters for TS07D and IGRF models
      REAL*8 BXGSM,BYGSM,BZGSM, HXGSM,HYGSM,HZGSM

c     TODO: change the below to paths to your local paths
c     STAICDIR is the directory that contains the static shielding coefficients
      CHARACTER*256 STATICDIR
      PARAMETER    (STATICDIR = '/Users/stephgk1/Downloads/TAIL_PAR/')
c     VARIABLEDIR is the directory where the variable coefficient directories will be
      CHARACTER*256 VARIABLEDIR
      PARAMETER(VARIABLEDIR = '/Users/stephgk1/sst19_fortran/')

c     The meta data included in the .par file
      REAL*8 COEFF_Q, COEFF_B_RMS,TILT
      INTEGER M_INX,N_INX

c
      CHARACTER*256 FILENAME
      REAL*8 VXGSE,VYGSE,VZGSE,PI
      INTEGER KK,IREAD,KREAD
      INTEGER IYEAR,IDAY,IHOUR,MIN,ISEC

c--------------------------------------------------------------------------
c//////////////////////////////////////////////////////////////////////////
c     The following reads the static coefficients that are used for the Equatorial Field
c     shielding fields. These only need to be loaded once as they are common for all times.
      DO 1001 IREAD=1,5
        WRITE(filename,'(A,A,I0,A)')
     *        TRIM(STATICDIR),'/tailamebhr',IREAD,'.par'
      OPEN (UNIT=1,FILE=FILENAME)
      READ (1,200) (TSS(KK,IREAD),KK=1,80)
 200  FORMAT(G17.10)
 1001 CLOSE(1)

      DO 1002 IREAD=1,5
        DO 1003 KREAD=1,4
        WRITE(filename,'(A,A,I0,I0,A)')
     *        TRIM(STATICDIR),'/tailamhr_o_',IREAD,KREAD,'.par'
      OPEN (UNIT=1,FILE=FILENAME)
      READ (1,200) (TSO(KK,IREAD,KREAD),KK=1,80)
 1003   CONTINUE
 1002 CLOSE(1)

      DO 1004 IREAD=1,5
        DO 1005 KREAD=1,4
        WRITE(filename,'(A,A,I0,I0,A)')
     *        TRIM(STATICDIR),'/tailamhr_e_',IREAD,KREAD,'.par'
      OPEN (UNIT=1,FILE=FILENAME)
      READ (1,200) (TSE(KK,IREAD,KREAD),KK=1,80)
 1005   CONTINUE
 1004 CLOSE(1)                                         

      PRINT *, '   SHIELDING COEFFICIENTS HAS BEEN READ INTO RAM'
c//////////////////////////////////////////////////////////////////////////
c--------------------------------------------------------------------------


c--------------------------------------------------------------------------
c//////////////////////////////////////////////////////////////////////////
c     Now, initialize the TS07D and IGRF models for the supplied date and time. If multiple times
c     are wanted, the following must be recalled for the new date and time.

c     The date and time is specified by the year, day of year (1-366), hour (0-23), minute (0-59),
c     and the second (0-59), note that the second field is not used for accessing the variable
c     parameter and coefficients file.
      IYEAR = 2005
      IDAY = 253
      IHOUR = 0
      MIN = 0
      ISEC = 0

c     Use the date and time fields to construct the variable parameter and coefficients file.
c     Or, if preferred, you could hardcode the path as:
c     filename = '/Users/username/Downloads/2015_075/2015_075_12_35.par'
      WRITE(filename,'(A,I0.4,A,I0.3,A,I0.4,A,I0.3,A,I0.2,A,I0.2,A)')
     *     TRIM(VARIABLEDIR),
     *     IYEAR,'_',IDAY,'/',IYEAR,'_',IDAY,'_',IHOUR,'_',MIN,'.par'

c     Now read variable parameter and coefficients file, the variable parameters and coefficients 
c     are stored in the array PARAMS. Additionaly, the dynamic pressure and the dipole tilt angle
c     are stored in PDYN and TILT respectively. All of these are variable or time-dependent values.
c     The other values stored in the file are not currently being used.
      OPEN (UNIT=1,FILE=filename,action='read') ! open the filed
      READ (1,100) (PARAMS(IREAD),IREAD=1,NTOT) ! read the variable coefficients and parameters
      READ (1,101) COEFF_Q ! the Q factor, related to chi squared, a measure of the goodness of fit
      READ (1,101) COEFF_B_RMS 
      READ (1,102) M_INX ! the number of azimuthal expansions in the equatorial field module (M=4)
      READ (1,102) N_INX ! the number of radial expansions in the equatorial field module (N=5)
      READ (1,101) PDYN ! the dynamic pressure for this time
      READ (1,101) TILT ! the dipole tilt for this time
 100  FORMAT(G15.6)                                         
 101  FORMAT(7x,G15.6)                                            
 102  FORMAT(7x,I15)                                            
      CLOSE(1)

C     To get the IGRF model in GSM coordinates, set the solar wind velocity vector to the following
      VXGSE=-400.d0
      VYGSE=   0.d0
      VZGSE=   0.d0
c     initializes the Geopack package for the supplied date and time
      call RECALC_08 (IYEAR,IDAY,IHOUR,MIN,ISEC,VXGSE,VYGSE,VZGSE)
c//////////////////////////////////////////////////////////////////////////
c--------------------------------------------------------------------------


c--------------------------------------------------------------------------
c//////////////////////////////////////////////////////////////////////////
c     Now the TS07D and the IGRF model can be evaluated

c     or just code in the values, e.g. XGSM = 5.0d0, units are in Earth radii (1 RE=6371.2 KM).
c     PRINT *, '  Enter S/C position: XGSM,YGSM,ZGSM '
c     READ *, XGSM,YGSM,ZGSM
      XGSM = -5.0d0
      YGSM = +4.2d0
      ZGSM = +1.2d0
      
c     Evaluate the IGRF model in GSM coordinates
      CALL IGRF_GSW_08 (XGSM,YGSM,ZGSM,HXGSM,HYGSM,HZGSM)

      IOPT = 0 ! the whole external magnetic field model (almost always what you want)
c     Now evaluate the TS07D model, note that we have two variables holding the dipole tilt,
c     PSI which comes from Geopack and TILT which comes from the variables file, which
c     ultimately also was computed using Geopack. So use either one.
      CALL TS07D_JULY_2017 (IOPT,PARMOD,PSI,XGSM,YGSM,ZGSM,
     *     BXGSM,BYGSM,BZGSM)

      PRINT *,'     Main field:'
      PRINT *,' HXGSM,HYGSM,HZGSM=',HXGSM,HYGSM,HZGSM,' nT'
      PRINT *,' '
      PRINT *,'   External field:'
      PRINT *,' BXGSM,BYGSM,BZGSM=',BXGSM,BYGSM,BZGSM,' nT'
      PRINT *,' '
      PI=4.D0*DATAN(1.D0)
      PRINT *,'  Geodipole tilt (degrees):', PSI*180d0/PI, ' deg.'
c//////////////////////////////////////////////////////////////////////////
c--------------------------------------------------------------------------

c     The program is complete
      END

C
c***********************************************************************
c


C
c***********************************************************************
c
      SUBROUTINE TS07D_JULY_2017 (IOPT,PARMOD,PS,X,Y,Z,BX,BY,BZ)     !  A DOUBLE-PRECISION SUBROUTINE
C
C  July 2017, G.K.Stephens, This routine was updated to be a double precision subroutine.
C  To indicate the update, the subroutine was renamed from EXTMODEL to TS07D_JULY_2017.
C  Additionally, this July 2017 update incorporates an upgraded version of the Bessel function evaluator,
C  provided by Jay Albert, which significantly speeds up the model. We thank Jay Albert for these
C  contributions.
C
C  This subroutine computes and returns the compoents of the external magnetic field vector in the
C  GSM coordinate system due to extraterrestrial currents, as described in the following references.
C  To compute the total magnetic field, the magnetic field vector from the Earth's internal sources
C  must be added. To compute this see Dr. Tsyganenko's most recent implementation of the IGRF
C  model in the GEOPACK package (Geopack-2008_dp.for).
C
C  References: (1) Tsyganenko, N. A., and M. I. Sitnov (2007), Magnetospheric configurations from a 
C                  high-resolution data-based magneticfield model, J. Geophys. Res., 112,A06225,
C                  doi:10.1029/2007JA012260.
C
C              (2) Sitnov, M. I., N. A. Tsyganenko, A. Y. Ukhorskiy, B. J. Anderson, H. Korth,
C                  A. T. Y. Lui, and P. C. Brandt (2010),Empirical modeling of a CIR‐driven magnetic
C                  storm, J. Geophys. Res., 115, A07231, doi:10.1029/2009JA015169.
C
C  Inputs:
C     IOPT - An option flag that allows to break out the magnetic field vector contributions from the
C        individual modules.
C        IOPT=0 - The total external magnetic field (most common)
C        IOPT=1 - The field due to only the shielding of the dipole field on the magnetopause
C        IOPT=2 - The field due to only the equatorial currents and its shielding field
C        IOPT=3 - The field due to only the Birkeland currents and its shielding field
C     PS - The Geodipole tilt angle in radians.
C     PARMOD - A 10-element array, in this model this input is not used and will have no impact on the
C        model evaluation. It is kept here because the TRACE_08 routine in the Geopack package requires
C        a consistent signature with other empirical magnetic field models that do use this parameter.
C     X,Y,Z - The Cartesian Geocentric position where the model will be evaluated in the GSM coordinate
C        system in units of Earth Radii (1 RE=6371.2 KM).
C
C  Common Block Inputs:
C     /INPUT/ PDYN - The solar wind dynamic pressure in nanoPascals (nPa)
C     /PARAM/ A - An 101-element array containing the variable (time-dependent) coefficients and
C        parameters used in evaluating the model
C     /TSS/ TSS - An 80x5-element array containing the static (time-independent) coefficients that are
C        used to shield the symmetric equatorial expansions
C     /TSO/ TSO - An 80x5x4-element array containing the static (time-independent) coefficients that are
C        used to shield the ODD axis-symmetric equatorial expansions
C     /TSE/ TSE - An 80x5x4-element array containing the static (time-independent) coefficients that are
C        used to shield the EVEN axis-symmetric equatorial expansions 
C
C  Outputs:
C     BX,BY,BZ - the evaluated magnetic field vector in the GSM coordinate system in units of nanoTesla (nT)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8    PARMOD(10)

      PARAMETER (NTOT=101)

      COMMON /INPUT/  PDYN
      DIMENSION BXTS(5),BXTO(5,4),BXTE(5,4)
      DIMENSION BYTS(5),BYTO(5,4),BYTE(5,4)
      DIMENSION BZTS(5),BZTO(5,4),BZTE(5,4)
      COMMON /PARAM/ A(NTOT)

      CALL EXTERN (IOPT,A,NTOT,PS,PDYN,X,Y,Z,BXCF,BYCF,BZCF,
     *BXTS,BYTS,BZTS,BXTO,BYTO,BZTO,BXTE,BYTE,BZTE,
     *BXR11,BYR11,BZR11,BXR12,BYR12,BZR12,BXR21a,BYR21a,BZR21a,BXR21s,
     *BYR21s,BZR21s,BX,BY,BZ)

      RETURN
C
      END
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      SUBROUTINE EXTERN (IOPGEN,A,NTOT,
     *  PS,PDYN,X,Y,Z,
     *  BXCF,BYCF,BZCF,
     *  BXTS,BYTS,BZTS,BXTO,BYTO,BZTO,BXTE,BYTE,BZTE,
     *  BXR11,BYR11,BZR11,BXR12,BYR12,BZR12,
     *  BXR21a,BYR21a,BZR21a,BXR21s,BYR21s,BZR21s,
     *  BX,BY,BZ)
C
C   IOPGEN - GENERAL OPTION FLAG:  IOPGEN=0 - CALCULATE TOTAL FIELD
C                                  IOPGEN=1 - DIPOLE SHIELDING ONLY
C                                  IOPGEN=2 - TAIL FIELD ONLY
C                                  IOPGEN=3 - BIRKELAND FIELD ONLY
C                                  IOPGEN=4 - RING CURRENT FIELD ONLY
C
      IMPLICIT  REAL * 8  (A - H, O - Z)
C
      DIMENSION A(NTOT)
C
      DIMENSION BXTS(5),BXTO(5,4),BXTE(5,4)
      DIMENSION BYTS(5),BYTO(5,4),BYTE(5,4)
      DIMENSION BZTS(5),BZTO(5,4),BZTE(5,4)
C
      COMMON /TAIL/ D  ! THE COMMON BLOCK FORWARDS TAIL SHEET THICKNESS
      COMMON /BIRKPAR/ XKAPPA1,XKAPPA2   !  SCALING FACTORS FOR BIRKELAND CURRENTS
      COMMON /G/ G,TW
      COMMON /RH0/ RH0
C
      DATA A0_A,A0_S0,A0_X0 /34.586D0,1.1960D0,3.4397D0/   !   SHUE ET AL. PARAMETERS
      DATA DSIG /0.005D0/, RH2 /-5.2D0/

      XAPPA=(PDYN/2.D0)**0.155   !   0.155 is the value obtained in TS05
      XAPPA3=XAPPA**3

      D=      A(96)
      RH0=    A(97)
      G=      A(98)
      XKAPPA1=A(99)
      XKAPPA2=A(100)
      TW=     A(101)       !   THIS PARAMETER CONTROLS THE IMF-INDUCED TWISTING (ADDED 04/21/06)
c
      XX=X*XAPPA  ! pressure scaling has been reinstated here
      YY=Y*XAPPA
      ZZ=Z*XAPPA

c     print *,XAPPA,PDYN
C
      SPS=DSIN(PS)
c
      X0=A0_X0/XAPPA   ! pressure scaling has been reinstated, even though these parameters are not used in this code
      AM=A0_A/XAPPA    ! pressure scaling has been reinstated, even though these parameters are not used in this code
      S0=A0_S0
C
C   CALCULATE THE IMF CLOCK ANGLE:
C
C        IF (BYIMF.EQ.0.D0.AND.BZIMF.EQ.0.D0) THEN
C            THETA=0.D0
C         ELSE
C            THETA=DATAN2(BYIMF,BZIMF)
C            IF (THETA.LE.0.D0) THETA=THETA+6.283185307D0
C        ENDIF
C
C       CT=COS(THETA)
C       ST=SIN(THETA)
C       YS=Y*CT-Z*ST
C       ZS=Z*CT+Y*ST
C
C       STHETAH=SIN(THETA/2.)**2
C
C  CALCULATE "IMF" COMPONENTS OUTSIDE THE MAGNETOPAUSE LAYER (HENCE BEGIN WITH "O")
C  THEY ARE NEEDED ONLY IF THE POINT (X,Y,Z) IS WITHIN THE TRANSITION MAGNETOPAUSE LAYER
C  OR OUTSIDE THE MAGNETOSPHERE:
C
C      FACTIMF=A(24)+A(25)*STHETAH
C
C      OIMFX=0.D0
C      OIMFY=BYIMF*FACTIMF
C      OIMFZ=BZIMF*FACTIMF
c
C =====================================================================
C  THIS FRAGMENT (BETWEEN THE ===== LINES) DISABLES THE CALCULATION OF THE MAGNETOPAUSE POSITION
C  IT SHOULD BE USED ONLY FOR THE FITTING (WE ASSUME THAT NO POINTS FROM THE SHEATH ARE PRESENT
C  IN THE DATASET, WHICH ITSELF IS STILL A QUESTION).
C
C  REMOVE IT IN THE FINAL VERSION.

C      SIGMA=0.D0
C      GOTO 1111
C======================================================================
c
C      R=SQRT(X**2+Y**2+Z**2)
C      XSS=X
C      ZSS=Z

C  1   XSOLD=XSS      !   BEGIN ITERATIVE SEARCH OF UNWARPED COORDS (TO FIND SIGMA)
C      ZSOLD=ZSS

C      RH=RH0+RH2*(ZSS/R)**2
C      SINPSAS=SPS/(1.D0+(R/RH)**3)**0.33333333D0
C      COSPSAS=DSQRT(1.D0-SINPSAS**2)
C      ZSS=X*SINPSAS+Z*COSPSAS
C      XSS=X*COSPSAS-Z*SINPSAS
C      DD=DABS(XSS-XSOLD)+DABS(ZSS-ZSOLD)
C      IF (DD.GT.1.D-6) GOTO 1
C                                END OF ITERATIVE SEARCH
C      RHO2=Y**2+ZSS**2
C      ASQ=AM**2
C      XMXM=AM+XSS-X0
C      IF (XMXM.LT.0.) XMXM=0. ! THE BOUNDARY IS A CYLINDER TAILWARD OF X=X0-AM
C      AXX0=XMXM**2
C      ARO=ASQ+RHO2
C      SIGMA=DSQRT((ARO+AXX0+SQRT((ARO+AXX0)**2-4.*ASQ*AXX0))/(2.*ASQ))

C==================================================================
C 1111 CONTINUE  !!!!!!!!!!!!  REMOVE IN THE FINAL VERSION
C==================================================================

C
C   NOW, THERE ARE THREE POSSIBLE CASES:
C    (1) INSIDE THE MAGNETOSPHERE   (SIGMA
C    (2) IN THE BOUNDARY LAYER
C    (3) OUTSIDE THE MAGNETOSPHERE AND B.LAYER
C       FIRST OF ALL, CONSIDER THE CASES (1) AND (2):
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C      IF (SIGMA.LT.S0+DSIG) THEN  !  CASES (1) OR (2); CALCULATE THE MODEL FIELD
C                                   (WITH THE POTENTIAL "PENETRATED" INTERCONNECTION FIELD):
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IF (IOPGEN.LE.1) THEN

C     print *,XX,YY,ZZ
         CALL SHLCAR3X3(XX,YY,ZZ,PS,CFX,CFY,CFZ)         !  DIPOLE SHIELDING FIELD
         BXCF=CFX  *XAPPA3
         BYCF=CFY  *XAPPA3
         BZCF=CFZ  *XAPPA3
      ELSE
         BXCF=0.D0
         BYCF=0.D0
         BZCF=0.D0
      ENDIF                                              !  DONE

      IF (IOPGEN.EQ.0.OR.IOPGEN.EQ.2) THEN
      CALL DEFORMED (PS,XX,YY,ZZ,                !  TAIL FIELD (THREE MODES)
     *   BXTS,BYTS,BZTS,BXTO,BYTO,BZTO,BXTE,BYTE,BZTE)
       ELSE
        DO 11 K=1,5

        BXTS(K)=0.D0
        BYTS(K)=0.D0
        BZTS(K)=0.D0

 11     CONTINUE

        DO 12 K=1,5

          DO 13 L=1,4

          BXTO(K,L)=0.D0
          BYTO(K,L)=0.D0
          BZTO(K,L)=0.D0

          BXTE(K,L)=0.D0
          BYTE(K,L)=0.D0
          BZTE(K,L)=0.D0

 13       CONTINUE

 12   CONTINUE

      ENDIF

      IF (IOPGEN.EQ.0.OR.IOPGEN.EQ.3) THEN

         CALL BIRK_TOT (PS,XX,YY,ZZ,BXR11,BYR11,BZR11,BXR12,BYR12,
     *   BZR12,BXR21a,BYR21a,BZR21a,BXR22a,BYR22a,BZR22a)    !   BIRKELAND FIELD (TWO MODES FOR R1 AND TWO MODES FOR R2)

         CALL BIRTOTSY (PS,XX,YY,ZZ,BXR11s,BYR11s,BZR11s,BXR12s,
     *   BYR12s,BZR12s,BXR21s,BYR21s,BZR21s,BXR22s,BYR22s,BZR22s)    !   "SYMMETRIC" BIRKELAND FIELD
C                                                                        (TWO MODES FOR R1s AND TWO MODES FOR R2s)
c                                                                        (but we actually use from here only R2s modes)
      ELSE
         BXR11=0.D0
         BYR11=0.D0
         BZR11=0.D0
         BXR12=0.D0
         BYR12=0.D0
         BZR12=0.D0
         BXR21a=0.D0
         BYR21a=0.D0
         BZR21a=0.D0
         BXR21s=0.D0
         BYR21s=0.D0
         BZR21s=0.D0
      ENDIF
C
C-----------------------------------------------------------
C
C    NOW, ADD UP ALL THE COMPONENTS:

      A_R11=A(92)
      A_R12=A(93)
      A_R21a=A(94)
      A_R21s=A(95)

      TX=0.D0
      TY=0.D0
      TZ=0.D0

C --- New tail structure -------------


        PDYN_0=2.D0   !   AVERAGE PRESSURE USED FOR NORMALIZATION

        P_FACTOR=DSQRT(PDYN/PDYN_0)-1.D0


      IND=1

        DO 911 K=1,5
         IND=IND+1
           TX=TX+(A(IND)+A(IND+45)*P_FACTOR)*BXTS(K)    !   2 - 6  &  47 - 51
           TY=TY+(A(IND)+A(IND+45)*P_FACTOR)*BYTS(K)
           TZ=TZ+(A(IND)+A(IND+45)*P_FACTOR)*BZTS(K)
 911     CONTINUE


        DO 912 K=1,5

          DO 913 L=1,4

          IND=IND+1

           TX=TX+(A(IND)+A(IND+45)*P_FACTOR)*BXTO(K,L)  !   7 -26  &  52 - 71
           TY=TY+(A(IND)+A(IND+45)*P_FACTOR)*BYTO(K,L)
           TZ=TZ+(A(IND)+A(IND+45)*P_FACTOR)*BZTO(K,L)

           TX=TX+(A(IND+20)+A(IND+65)*P_FACTOR)*BXTE(K,L) !   27 -46  &  72 - 91
           TY=TY+(A(IND+20)+A(IND+65)*P_FACTOR)*BYTE(K,L)
           TZ=TZ+(A(IND+20)+A(IND+65)*P_FACTOR)*BZTE(K,L)

 913       CONTINUE

 912   CONTINUE
      
      BBX=A(1)*BXCF+TX+
     * A_R11*BXR11+A_R12*BXR12+A_R21a*BXR21a+A_R21s*BXR21s

      BBY=A(1)*BYCF+TY+
     * A_R11*BYR11+A_R12*BYR12+A_R21a*BYR21a+A_R21s*BYR21s

      BBZ=A(1)*BZCF+TZ+
     * A_R11*BZR11+A_R12*BZR12+A_R21a*BZR21a+A_R21s*BZR21s
C
c   -----------------------------------------------------------
C
C   AND WE HAVE THE TOTAL EXTERNAL FIELD.
C
C  NOW, LET US CHECK WHETHER WE HAVE THE CASE (1). IF YES - WE ARE DONE:
C
C      IF (SIGMA.LT.S0-DSIG) THEN    !  (X,Y,Z) IS INSIDE THE MAGNETOSPHERE

C       BX=BBX
C       BY=BBY
C       BZ=BBZ
C                     ELSE           !  THIS IS THE MOST COMPLEX CASE: WE ARE INSIDE
C                                             THE INTERPOLATION REGION
C       FINT=0.5*(1.-(SIGMA-S0)/DSIG)
C       FEXT=0.5*(1.+(SIGMA-S0)/DSIG)
C
C       CALL DIPOLE (PS,X,Y,Z,QX,QY,QZ)
C       BX=(BBX+QX)*FINT+OIMFX*FEXT -QX
C       BY=(BBY+QY)*FINT+OIMFY*FEXT -QY
C       BZ=(BBZ+QZ)*FINT+OIMFZ*FEXT -QZ
c
C        ENDIF  !   THE CASES (1) AND (2) ARE EXHAUSTED; THE ONLY REMAINING
C                      POSSIBILITY IS NOW THE CASE (3):
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C        ELSE
C                CALL DIPOLE (PS,X,Y,Z,QX,QY,QZ)
C                BX=OIMFX-QX
C                BY=OIMFY-QY
C                BZ=OIMFZ-QZ
C        ENDIF
C

       BX=BBX
       BY=BBY
       BZ=BBZ

      RETURN
      END
C
C XXXXXXXXXXXXXXXXXXXXXXXXXXX11/15/05 16:06 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
         SUBROUTINE  SHLCAR3X3(X,Y,Z,PS,BX,BY,BZ)
C
C   THIS S/R RETURNS THE SHIELDING FIELD FOR THE EARTH'S DIPOLE,
C   REPRESENTED BY  2x3x3=18 "CARTESIAN" HARMONICS, tilted with respect
C   to the z=0 plane (see NB#4, p.74-74)
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  The 36 coefficients enter in pairs in the amplitudes of the "cartesian"
c    harmonics (A(1)-A(36).
c  The 14 nonlinear parameters (A(37)-A(50) are the scales Pi,Ri,Qi,and Si
C   entering the arguments of exponents, sines, and cosines in each of the
C   18 "Cartesian" harmonics  PLUS TWO TILT ANGLES FOR THE CARTESIAN HARMONICS
C       (ONE FOR THE PSI=0 MODE AND ANOTHER FOR THE PSI=90 MODE)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
      IMPLICIT  REAL * 8  (A - H, O - Z)
C
      DIMENSION A(50)
      DATA A/-901.2327248,895.8011176,817.6208321,-845.5880889,
     *-83.73539535,86.58542841,336.8781402,-329.3619944,-311.2947120,
     *308.6011161,31.94469304,-31.30824526,125.8739681,-372.3384278,
     *-235.4720434,286.7594095,21.86305585,-27.42344605,-150.4874688,
     *2.669338538,1.395023949,-.5540427503,-56.85224007,3.681827033,
     *-43.48705106,5.103131905,1.073551279,-.6673083508,12.21404266,
     *4.177465543,5.799964188,-.3977802319,-1.044652977,.5703560010,
     *3.536082962,-3.222069852,9.620648151,6.082014949,27.75216226,
     *12.44199571,5.122226936,6.982039615,20.12149582,6.150973118,
     *4.663639687,15.73319647,2.303504968,5.840511214,.8385953499E-01,
     *.3477844929/
C
         P1=A(37)
         P2=A(38)
         P3=A(39)
         R1=A(40)
         R2=A(41)
         R3=A(42)
         Q1=A(43)
         Q2=A(44)
         Q3=A(45)
         S1=A(46)
         S2=A(47)
         S3=A(48)

         T1  =A(49)
         T2  =A(50)
C
         CPS=DCOS(PS)
         SPS=DSIN(PS)
         S2PS=2.D0*CPS      !   MODIFIED HERE (INSTEAD OF SIN(3*PS) I TRY SIN(2*PS)

C
           ST1=DSIN(PS*T1)
           CT1=DCOS(PS*T1)
           ST2=DSIN(PS*T2)
           CT2=DCOS(PS*T2)

C     print *,X,Z

            X1=X*CT1-Z*ST1

C         print *,'X1=',X1

            Z1=X*ST1+Z*CT1
            X2=X*CT2-Z*ST2
            Z2=X*ST2+Z*CT2
C
C
c  MAKE THE TERMS IN THE 1ST SUM ("PERPENDICULAR" SYMMETRY):
C
C       I=1
C
        SQPR= DSQRT(1.D0/P1**2+1.D0/R1**2)
        CYP = DCOS(Y/P1)
        SYP = DSIN(Y/P1)
        CZR = DCOS(Z1/R1)
        SZR = DSIN(Z1/R1)
c       print *,X1
        EXPR= DEXP(SQPR*X1)
        FX1 =-SQPR*EXPR*CYP*SZR
        HY1 = EXPR/P1*SYP*SZR
        FZ1 =-EXPR*CYP/R1*CZR
        HX1 = FX1*CT1+FZ1*ST1
        HZ1 =-FX1*ST1+FZ1*CT1

        SQPR= DSQRT(1.D0/P1**2+1.D0/R2**2)
        CYP = DCOS(Y/P1)
        SYP = DSIN(Y/P1)
        CZR = DCOS(Z1/R2)
        SZR = DSIN(Z1/R2)
        EXPR= DEXP(SQPR*X1)
        FX2 =-SQPR*EXPR*CYP*SZR
        HY2 = EXPR/P1*SYP*SZR
        FZ2 =-EXPR*CYP/R2*CZR
        HX2 = FX2*CT1+FZ2*ST1
        HZ2 =-FX2*ST1+FZ2*CT1

        SQPR= DSQRT(1.D0/P1**2+1.D0/R3**2)
        CYP = DCOS(Y/P1)
        SYP = DSIN(Y/P1)
        CZR = DCOS(Z1/R3)
        SZR = DSIN(Z1/R3)
        EXPR= DEXP(SQPR*X1)
        FX3 =-EXPR*CYP*(SQPR*Z1*CZR+SZR/R3*(X1+1.D0/SQPR))
        HY3 = EXPR/P1*SYP*(Z1*CZR+X1/R3*SZR/SQPR)
        FZ3 =-EXPR*CYP*(CZR*(1.D0+X1/R3**2/SQPR)-Z1/R3*SZR)
        HX3 = FX3*CT1+FZ3*ST1
        HZ3 =-FX3*ST1+FZ3*CT1
C
C       I=2:
C
        SQPR= DSQRT(1.D0/P2**2+1.D0/R1**2)
        CYP = DCOS(Y/P2)
        SYP = DSIN(Y/P2)
        CZR = DCOS(Z1/R1)
        SZR = DSIN(Z1/R1)
        EXPR= DEXP(SQPR*X1)
        FX4 =-SQPR*EXPR*CYP*SZR
        HY4 = EXPR/P2*SYP*SZR
        FZ4 =-EXPR*CYP/R1*CZR
        HX4 = FX4*CT1+FZ4*ST1
        HZ4 =-FX4*ST1+FZ4*CT1

        SQPR= DSQRT(1.D0/P2**2+1.D0/R2**2)
        CYP = DCOS(Y/P2)
        SYP = DSIN(Y/P2)
        CZR = DCOS(Z1/R2)
        SZR = DSIN(Z1/R2)
        EXPR= DEXP(SQPR*X1)
        FX5 =-SQPR*EXPR*CYP*SZR
        HY5 = EXPR/P2*SYP*SZR
        FZ5 =-EXPR*CYP/R2*CZR
        HX5 = FX5*CT1+FZ5*ST1
        HZ5 =-FX5*ST1+FZ5*CT1

        SQPR= DSQRT(1.D0/P2**2+1.D0/R3**2)
        CYP = DCOS(Y/P2)
        SYP = DSIN(Y/P2)
        CZR = DCOS(Z1/R3)
        SZR = DSIN(Z1/R3)
        EXPR= DEXP(SQPR*X1)
        FX6 =-EXPR*CYP*(SQPR*Z1*CZR+SZR/R3*(X1+1.D0/SQPR))
        HY6 = EXPR/P2*SYP*(Z1*CZR+X1/R3*SZR/SQPR)
        FZ6 =-EXPR*CYP*(CZR*(1.D0+X1/R3**2/SQPR)-Z1/R3*SZR)
        HX6 = FX6*CT1+FZ6*ST1
        HZ6 =-FX6*ST1+FZ6*CT1
C
C       I=3:
C
        SQPR= DSQRT(1.D0/P3**2+1.D0/R1**2)
        CYP = DCOS(Y/P3)
        SYP = DSIN(Y/P3)
        CZR = DCOS(Z1/R1)
        SZR = DSIN(Z1/R1)
        EXPR= DEXP(SQPR*X1)
        FX7 =-SQPR*EXPR*CYP*SZR
        HY7 = EXPR/P3*SYP*SZR
        FZ7 =-EXPR*CYP/R1*CZR
        HX7 = FX7*CT1+FZ7*ST1
        HZ7 =-FX7*ST1+FZ7*CT1

        SQPR= DSQRT(1.D0/P3**2+1.D0/R2**2)
        CYP = DCOS(Y/P3)
        SYP = DSIN(Y/P3)
        CZR = DCOS(Z1/R2)
        SZR = DSIN(Z1/R2)
        EXPR= DEXP(SQPR*X1)
        FX8 =-SQPR*EXPR*CYP*SZR
        HY8 = EXPR/P3*SYP*SZR
        FZ8 =-EXPR*CYP/R2*CZR
        HX8 = FX8*CT1+FZ8*ST1
        HZ8 =-FX8*ST1+FZ8*CT1

        SQPR= DSQRT(1.D0/P3**2+1.D0/R3**2)
        CYP = DCOS(Y/P3)
        SYP = DSIN(Y/P3)
        CZR = DCOS(Z1/R3)
        SZR = DSIN(Z1/R3)
        EXPR= DEXP(SQPR*X1)
        FX9 =-EXPR*CYP*(SQPR*Z1*CZR+SZR/R3*(X1+1.D0/SQPR))
        HY9 = EXPR/P3*SYP*(Z1*CZR+X1/R3*SZR/SQPR)
        FZ9 =-EXPR*CYP*(CZR*(1.D0+X1/R3**2/SQPR)-Z1/R3*SZR)
        HX9 = FX9*CT1+FZ9*ST1
        HZ9 =-FX9*ST1+FZ9*CT1


       A1=A(1)+A(2)*CPS
       A2=A(3)+A(4)*CPS
       A3=A(5)+A(6)*CPS
       A4=A(7)+A(8)*CPS
       A5=A(9)+A(10)*CPS
       A6=A(11)+A(12)*CPS
       A7=A(13)+A(14)*CPS
       A8=A(15)+A(16)*CPS
       A9=A(17)+A(18)*CPS
       BX=A1*HX1+A2*HX2+A3*HX3+A4*HX4+A5*HX5+A6*HX6+A7*HX7+A8*HX8+A9*HX9
       BY=A1*HY1+A2*HY2+A3*HY3+A4*HY4+A5*HY5+A6*HY6+A7*HY7+A8*HY8+A9*HY9
       BZ=A1*HZ1+A2*HZ2+A3*HZ3+A4*HZ4+A5*HZ5+A6*HZ6+A7*HZ7+A8*HZ8+A9*HZ9


c  MAKE THE TERMS IN THE 2ND SUM ("PARALLEL" SYMMETRY):
C
C       I=1
C
       SQQS= DSQRT(1.D0/Q1**2+1.D0/S1**2)
       CYQ = DCOS(Y/Q1)
       SYQ = DSIN(Y/Q1)
       CZS = DCOS(Z2/S1)
       SZS = DSIN(Z2/S1)
       EXQS= DEXP(SQQS*X2)
       FX1 =-SQQS*EXQS*CYQ*CZS *SPS
       HY1 = EXQS/Q1*SYQ*CZS   *SPS
       FZ1 = EXQS*CYQ/S1*SZS   *SPS
       HX1 = FX1*CT2+FZ1*ST2
       HZ1 =-FX1*ST2+FZ1*CT2

       SQQS= DSQRT(1.D0/Q1**2+1.D0/S2**2)
       CYQ = DCOS(Y/Q1)
       SYQ = DSIN(Y/Q1)
       CZS = DCOS(Z2/S2)
       SZS = DSIN(Z2/S2)
       EXQS= DEXP(SQQS*X2)
       FX2 =-SQQS*EXQS*CYQ*CZS *SPS
       HY2 = EXQS/Q1*SYQ*CZS   *SPS
       FZ2 = EXQS*CYQ/S2*SZS   *SPS
       HX2 = FX2*CT2+FZ2*ST2
       HZ2 =-FX2*ST2+FZ2*CT2

       SQQS= DSQRT(1.D0/Q1**2+1.D0/S3**2)
       CYQ = DCOS(Y/Q1)
       SYQ = DSIN(Y/Q1)
       CZS = DCOS(Z2/S3)
       SZS = DSIN(Z2/S3)
       EXQS= DEXP(SQQS*X2)
       FX3 =-SQQS*EXQS*CYQ*CZS *SPS
       HY3 = EXQS/Q1*SYQ*CZS   *SPS
       FZ3 = EXQS*CYQ/S3*SZS   *SPS
       HX3 = FX3*CT2+FZ3*ST2
       HZ3 =-FX3*ST2+FZ3*CT2
C
C       I=2
C
       SQQS= DSQRT(1.D0/Q2**2+1.D0/S1**2)
       CYQ = DCOS(Y/Q2)
       SYQ = DSIN(Y/Q2)
       CZS = DCOS(Z2/S1)
       SZS = DSIN(Z2/S1)
       EXQS= DEXP(SQQS*X2)
       FX4 =-SQQS*EXQS*CYQ*CZS *SPS
       HY4 = EXQS/Q2*SYQ*CZS   *SPS
       FZ4 = EXQS*CYQ/S1*SZS   *SPS
       HX4 = FX4*CT2+FZ4*ST2
       HZ4 =-FX4*ST2+FZ4*CT2

       SQQS= DSQRT(1.D0/Q2**2+1.D0/S2**2)
       CYQ = DCOS(Y/Q2)
       SYQ = DSIN(Y/Q2)
       CZS = DCOS(Z2/S2)
       SZS = DSIN(Z2/S2)
       EXQS= DEXP(SQQS*X2)
       FX5 =-SQQS*EXQS*CYQ*CZS *SPS
       HY5 = EXQS/Q2*SYQ*CZS   *SPS
       FZ5 = EXQS*CYQ/S2*SZS   *SPS
       HX5 = FX5*CT2+FZ5*ST2
       HZ5 =-FX5*ST2+FZ5*CT2

       SQQS= DSQRT(1.D0/Q2**2+1.D0/S3**2)
       CYQ = DCOS(Y/Q2)
       SYQ = DSIN(Y/Q2)
       CZS = DCOS(Z2/S3)
       SZS = DSIN(Z2/S3)
       EXQS= DEXP(SQQS*X2)
       FX6 =-SQQS*EXQS*CYQ*CZS *SPS
       HY6 = EXQS/Q2*SYQ*CZS   *SPS
       FZ6 = EXQS*CYQ/S3*SZS   *SPS
       HX6 = FX6*CT2+FZ6*ST2
       HZ6 =-FX6*ST2+FZ6*CT2
C
C       I=3
C
       SQQS= DSQRT(1.D0/Q3**2+1.D0/S1**2)
       CYQ = DCOS(Y/Q3)
       SYQ = DSIN(Y/Q3)
       CZS = DCOS(Z2/S1)
       SZS = DSIN(Z2/S1)
       EXQS= DEXP(SQQS*X2)
       FX7 =-SQQS*EXQS*CYQ*CZS *SPS
       HY7 = EXQS/Q3*SYQ*CZS   *SPS
       FZ7 = EXQS*CYQ/S1*SZS   *SPS
       HX7 = FX7*CT2+FZ7*ST2
       HZ7 =-FX7*ST2+FZ7*CT2

       SQQS= DSQRT(1.D0/Q3**2+1.D0/S2**2)
       CYQ = DCOS(Y/Q3)
       SYQ = DSIN(Y/Q3)
       CZS = DCOS(Z2/S2)
       SZS = DSIN(Z2/S2)
       EXQS= DEXP(SQQS*X2)
       FX8 =-SQQS*EXQS*CYQ*CZS *SPS
       HY8 = EXQS/Q3*SYQ*CZS   *SPS
       FZ8 = EXQS*CYQ/S2*SZS   *SPS
       HX8 = FX8*CT2+FZ8*ST2
       HZ8 =-FX8*ST2+FZ8*CT2

       SQQS= DSQRT(1.D0/Q3**2+1.D0/S3**2)
       CYQ = DCOS(Y/Q3)
       SYQ = DSIN(Y/Q3)
       CZS = DCOS(Z2/S3)
       SZS = DSIN(Z2/S3)
       EXQS= DEXP(SQQS*X2)
       FX9 =-SQQS*EXQS*CYQ*CZS *SPS
       HY9 = EXQS/Q3*SYQ*CZS   *SPS
       FZ9 = EXQS*CYQ/S3*SZS   *SPS
       HX9 = FX9*CT2+FZ9*ST2
       HZ9 =-FX9*ST2+FZ9*CT2

       A1=A(19)+A(20)*S2PS
       A2=A(21)+A(22)*S2PS
       A3=A(23)+A(24)*S2PS
       A4=A(25)+A(26)*S2PS
       A5=A(27)+A(28)*S2PS
       A6=A(29)+A(30)*S2PS
       A7=A(31)+A(32)*S2PS
       A8=A(33)+A(34)*S2PS
       A9=A(35)+A(36)*S2PS

       BX=BX+A1*HX1+A2*HX2+A3*HX3+A4*HX4+A5*HX5+A6*HX6+A7*HX7+A8*HX8
     *   +A9*HX9
       BY=BY+A1*HY1+A2*HY2+A3*HY3+A4*HY4+A5*HY5+A6*HY6+A7*HY7+A8*HY8
     *   +A9*HY9
       BZ=BZ+A1*HZ1+A2*HZ2+A3*HZ3+A4*HZ4+A5*HZ5+A6*HZ6+A7*HZ7+A8*HZ8
     *   +A9*HZ9
C
       RETURN
       END
c
c############################################################################
c
C
      SUBROUTINE DEFORMED (PS,X,Y,Z,
     *   BXS,BYS,BZS,BXO,BYO,BZO,BXE,BYE,BZE)
C
C    CALCULATES GSM COMPONENTS OF 104 UNIT-AMPLITUDE TAIL FIELD MODES,
C    TAKING INTO ACCOUNT BOTH EFFECTS OF DIPOLE TILT:
C    WARPING IN Y-Z (DONE BY THE S/R WARPED) AND BENDING IN X-Z (DONE BY THIS SUBROUTINE)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION BXS(5),BXO(5,4),BXE(5,4)
      DIMENSION BYS(5),BYO(5,4),BYE(5,4)
      DIMENSION BZS(5),BZO(5,4),BZE(5,4)
C
      DIMENSION BXASS(5),BXASO(5,4),BXASE(5,4)
      DIMENSION BYASS(5),BYASO(5,4),BYASE(5,4)
      DIMENSION BZASS(5),BZASO(5,4),BZASE(5,4)
C
      COMMON /RH0/ RH0
      DATA RH2,IEPS /-5.2D0,3/
C
C  RH0,RH1,RH2, AND IEPS CONTROL THE TILT-RELATED DEFORMATION OF THE TAIL FIELD
C
      SPS=DSIN(PS)
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
      XAS=X*CPSAS-Z*SPSAS
      ZAS=X*SPSAS+Z*CPSAS
C
      FACPS=SPS/CPSAS*(DFDR+DFDRH*DRHDR)/R
      PSASX=FACPS*X
      PSASY=FACPS*Y
      PSASZ=FACPS*Z+SPS/CPSAS*DFDRH*DRHDZ
C
      DXASDX=CPSAS-ZAS*PSASX
      DXASDY=-ZAS*PSASY
      DXASDZ=-SPSAS-ZAS*PSASZ
      DZASDX=SPSAS+XAS*PSASX
      DZASDY=XAS*PSASY
      DZASDZ=CPSAS+XAS*PSASZ
      FAC1=DXASDZ*DZASDY-DXASDY*DZASDZ
      FAC2=DXASDX*DZASDZ-DXASDZ*DZASDX
      FAC3=DZASDX*DXASDY-DXASDX*DZASDY
C
C     DEFORM:
C
      CALL WARPED(PS,XAS,Y,ZAS,
     *   BXASS,BYASS,BZASS,BXASO,BYASO,BZASO,BXASE,BYASE,BZASE)
C
C --- New tail structure -------------

        DO 11 K=1,5

        BXS(K)=BXASS(K)*DZASDZ-BZASS(K)*DXASDZ+BYASS(K)*FAC1
        BYS(K)=BYASS(K)*FAC2
        BZS(K)=BZASS(K)*DXASDX-BXASS(K)*DZASDX+BYASS(K)*FAC3

 11     CONTINUE

        DO 12 K=1,5

          DO 13 L=1,4

          BXO(K,L)=BXASO(K,L)*DZASDZ-BZASO(K,L)*DXASDZ
     *            +BYASO(K,L)*FAC1
          BYO(K,L)=BYASO(K,L)*FAC2
          BZO(K,L)=BZASO(K,L)*DXASDX-BXASO(K,L)*DZASDX
     *            +BYASO(K,L)*FAC3

          BXE(K,L)=BXASE(K,L)*DZASDZ-BZASE(K,L)*DXASDZ
     *            +BYASE(K,L)*FAC1
          BYE(K,L)=BYASE(K,L)*FAC2
          BZE(K,L)=BZASE(K,L)*DXASDX-BXASE(K,L)*DZASDX
     *            +BYASE(K,L)*FAC3

 13       CONTINUE

 12   CONTINUE

C ------------------------------------
C
      RETURN
      END
C
C------------------------------------------------------------------
c
C
      SUBROUTINE WARPED (PS,X,Y,Z,
     *   BXS,BYS,BZS,BXO,BYO,BZO,BXE,BYE,BZE)
C
C   CALCULATES GSM COMPONENTS OF THE WARPED FIELD FOR TWO TAIL UNIT MODES.
C   THE WARPING DEFORMATION IS IMPOSED ON THE UNWARPED FIELD, COMPUTED
C   BY THE S/R "UNWARPED".  THE WARPING PARAMETERS WERE TAKEN FROM THE
C   RESULTS OF GEOTAIL OBSERVATIONS (TSYGANENKO ET AL. [1998]).
C   NB # 6, P.106, OCT 12, 2000.
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION BXS(5),BXO(5,4),BXE(5,4)
      DIMENSION BYS(5),BYO(5,4),BYE(5,4)
      DIMENSION BZS(5),BZO(5,4),BZE(5,4)
C
      DIMENSION BX_ASS(5),BX_ASO(5,4),BX_ASE(5,4)
      DIMENSION BY_ASS(5),BY_ASO(5,4),BY_ASE(5,4)
      DIMENSION BZ_ASS(5),BZ_ASO(5,4),BZ_ASE(5,4)
C
      COMMON /G/ G,TW
      DGDX=0.D0
      XL=20.D0
      DXLDX=0.D0

      SPS=DSIN(PS)
      RHO2=Y**2+Z**2
      RHO=DSQRT(RHO2)

      IF (Y.EQ.0.D0.AND.Z.EQ.0.D0) THEN
       PHI=0.D0
       CPHI=1.D0
       SPHI=0.D0
      ELSE
       PHI=DATAN2(Z,Y)
       CPHI=Y/RHO
       SPHI=Z/RHO
      ENDIF

      RR4L4=RHO/(RHO2**2+XL**4)

      F=PHI+G*RHO2*RR4L4*CPHI*SPS +TW*(X/10.D0)
      DFDPHI=1.D0-G*RHO2*RR4L4*SPHI*SPS
      DFDRHO=G*RR4L4**2*(3.D0*XL**4-RHO2**2)*CPHI*SPS
      DFDX=RR4L4*CPHI*SPS*(DGDX*RHO2-G*RHO*RR4L4*4.D0*XL**3*DXLDX)
     *  +TW/10.D0        !  THE LAST TERM DESCRIBES THE IMF-INDUCED TWISTING (ADDED 04/21/06)

      CF=DCOS(F)
      SF=DSIN(F)
      YAS=RHO*CF
      ZAS=RHO*SF

      CALL UNWARPED (X,YAS,ZAS,
     *   BX_ASS,BY_ASS,BZ_ASS,
     *   BX_ASO,BY_ASO,BZ_ASO,
     *   BX_ASE,BY_ASE,BZ_ASE)
C
        DO 11 K=1,5
C ------------------------------------------- Deforming symmetric modules
      BRHO_AS =  BY_ASS(K)*CF+BZ_ASS(K)*SF
      BPHI_AS = -BY_ASS(K)*SF+BZ_ASS(K)*CF

      BRHO_S = BRHO_AS*DFDPHI
      BPHI_S = BPHI_AS-RHO*(BX_ASS(K)*DFDX+BRHO_AS*DFDRHO)

        BXS(K)=BX_ASS(K)*DFDPHI
        BYS(K)=BRHO_S*CPHI-BPHI_S*SPHI
        BZS(K)=BRHO_S*SPHI+BPHI_S*CPHI

 11     CONTINUE

        DO 12 K=1,5

          DO 13 L=1,4
C -------------------------------------------- Deforming odd modules
      BRHO_AS =  BY_ASO(K,L)*CF+BZ_ASO(K,L)*SF
      BPHI_AS = -BY_ASO(K,L)*SF+BZ_ASO(K,L)*CF

      BRHO_S = BRHO_AS*DFDPHI
      BPHI_S = BPHI_AS-RHO*(BX_ASO(K,L)*DFDX+BRHO_AS*DFDRHO)

          BXO(K,L)=BX_ASO(K,L)*DFDPHI
          BYO(K,L)=BRHO_S*CPHI-BPHI_S*SPHI
          BZO(K,L)=BRHO_S*SPHI+BPHI_S*CPHI
C ------------------------------------------- Deforming even modules
      BRHO_AS =  BY_ASE(K,L)*CF+BZ_ASE(K,L)*SF
      BPHI_AS = -BY_ASE(K,L)*SF+BZ_ASE(K,L)*CF

      BRHO_S = BRHO_AS*DFDPHI
      BPHI_S = BPHI_AS-RHO*(BX_ASE(K,L)*DFDX+BRHO_AS*DFDRHO)

          BXE(K,L)=BX_ASE(K,L)*DFDPHI
          BYE(K,L)=BRHO_S*CPHI-BPHI_S*SPHI
          BZE(K,L)=BRHO_S*SPHI+BPHI_S*CPHI

 13       CONTINUE

 12   CONTINUE

      RETURN
      END

C==========================================================================
C     July 2017, G.K.Stephens, this routine updated to incorporate Jay Albert's
C     improvements to the Bessel function evaluation.
       SUBROUTINE UNWARPED (X,Y,Z,BXS,BYS,BZS,BXO,BYO,BZO,BXE,BYE,BZE)

C    CALCULATES GSM COMPONENTS OF THE SHIELDED FIELD OF 45 TAIL MODES WITH UNIT
C    AMPLITUDES,  WITHOUT ANY WARPING OR BENDING.  NONLINEAR PARAMETERS OF THE MODES
C    ARE FORWARDED HERE VIA A COMMON BLOCK /TAIL/.
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /TSS/ TSS(80,5)
      COMMON /TSO/ TSO(80,5,4)
      COMMON /TSE/ TSE(80,5,4)

      COMMON /TAIL/ D0

      DIMENSION BXS(5),BXO(5,4),BXE(5,4)
      DIMENSION BYS(5),BYO(5,4),BYE(5,4)
      DIMENSION BZS(5),BZO(5,4),BZE(5,4)

      DIMENSION AJM(0:5),AJMD(0:5)
      RNOT=20.0     !    Rho_0 - scale parameter along the tail axis

C
C --- New tail structure -------------
C     
C
      DO 11 K=1,5

         RHO=DSQRT(X*X+Y*Y)

         RKMR = DBLE(K)*RHO/RNOT

C     July 2017, G.K.Stephens, all the Bessel functions are now evaluated first,
C     and passed into the subroutines
         call bessJJ(5,RKMR, AJM) !!! get all n in one call
         DO 3 m=1,5
            AJMD(m)=AJM(m-1)-m*AJM(m)/RKMR
 3       CONTINUE
         AJMD(0)=-AJM(1)

         CALL TAILSHT_S (K,X,Y,Z,AJM,BXSK,BYSK,BZSK)
         CALL SHTBNORM_S (K,X,Y,Z,HXSK,HYSK,HZSK)

         BXS(K)=BXSK+HXSK
         BYS(K)=BYSK+HYSK
         BZS(K)=BZSK+HZSK

         DO 13 L=1,4
            
            CALL TAILSHT_OE (1,K,L,X,Y,Z,AJM,AJMD,BXOKL,BYOKL,BZOKL)
            CALL SHTBNORM_O (  K,L,X,Y,Z,HXOKL,HYOKL,HZOKL)

            BXO(K,L)=BXOKL+HXOKL
            BYO(K,L)=BYOKL+HYOKL
            BZO(K,L)=BZOKL+HZOKL

            CALL TAILSHT_OE (0,K,L,X,Y,Z,AJM,AJMD,BXEKL,BYEKL,BZEKL)
            CALL SHTBNORM_E (  K,L,X,Y,Z,HXEKL,HYEKL,HZEKL)

            BXE(K,L)=BXEKL+HXEKL
            BYE(K,L)=BYEKL+HYEKL
            BZE(K,L)=BZEKL+HZEKL
            
 13      CONTINUE
         
 11   CONTINUE
      
      RETURN
      END

C==========================================================================
C     July 2017, G.K.Stephens, this routine updated to incorporate Jay Albert's
C     improvements to the Bessel function evaluation. The Bessel function
C     values are now precomputed and  passed into this rather than computed 
C     inside this routine.
      SUBROUTINE TAILSHT_S (M,X,Y,Z,AJM,BX,BY,BZ)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /TAIL/ D  ! THE COMMON BLOCKS FORWARDS TAIL SHEET THICKNESS
      DIMENSION AJM(0:5),AJMD(0:5)
C-----------------------------------------------------------------------------------
C
      RNOT=20.0        !    This can be replaced by introducing them
      DLTK=1.0         !    through the above common block


      RHO=DSQRT(X*X+Y*Y)
      CSPHI=X/RHO
      SNPHI=Y/RHO
C
      DKM=1.D0+(M-1)*DLTK
      RKM=DKM/RNOT
C
      RKMZ=RKM*Z
      RKMR=RKM*RHO

      ZD=DSQRT(Z*Z+D*D)
C
      RJ0=AJM(0)
      RJ1=AJM(1)
c     July 2017, G.K.Stephens, Bessel functions are now passed in.
c      RJ0=bessj0(RKMR)
c      RJ1=bessj1(RKMR)
      REX=DEXP(RKM*ZD)
C
      BX=RKMZ*RJ1*CSPHI/ZD/REX
      BY=RKMZ*RJ1*SNPHI/ZD/REX
      BZ=RKM*RJ0/REX
C
C    CALCULATION OF THE MAGNETOTAIL CURRENT CONTRIBUTION IS FINISHED
C
      RETURN
      END

C==========================================================================
C     July 2017, G.K.Stephens, this routine updated to incorporate Jay Albert's
C     improvements to the Bessel function evaluation.
      SUBROUTINE  SHTBNORM_S(K,X,Y,Z,FX,FY,FZ) ! modified SHTBNORM_S
      IMPLICIT  REAL * 8  (A - H, O - Z)
      DIMENSION AK(5)
      DIMENSION AJM(0:14),AJMD(0:14)
      COMMON /TSS/ TSS(80,5)
      
      AK(1)=TSS(76,K)
      AK(2)=TSS(77,K)
      AK(3)=TSS(78,K)
      AK(4)=TSS(79,K)
      AK(5)=TSS(80,K)
      
      phi=DATAN2(Y,X)
      RHO=dsqrt(X*X+Y*Y)
      if (RHO.lt.1.D-8) then
         RHOI=1.D8
      else
         RHOI=1.D0/RHO
      endif
      DPDX=-Y*RHOI*RHOI
      DPDY= X*RHOI*RHOI
      
      FX=0.D0
      FY=0.D0
      FZ=0.D0
      
      DO 2 n=1,5
         AKN=dabs(AK(n))
         AKNR=AKN*RHO
         if (AKNR.lt.1.D-8) then
            AKNRI=1.D8
         else
            AKNRI=1.D0/AKNR
         endif
         
         CHZ=dcosh(Z*AKN)
         SHZ=dsinh(Z*AKN)
           
         call bessJJ(14,AKNR, AJM) !!! get all n in one call
         DO 3 m=1,14
            AJMD(m)=AJM(m-1)-m*AJM(m)*AKNRI
 3       CONTINUE
         AJMD(0)=-AJM(1)
         
         DO 4 m=0,14
            CMP=DCOS(m*phi)
            SMP=DSIN(m*phi)
            
            HX1=m*DPDX*SMP*SHZ*AJM(m)
            HX2=-AKN*X*RHOI*CMP*SHZ*AJMD(m)
            HX=HX1+HX2
            HY1=m*DPDY*SMP*SHZ*AJM(m)
            HY2=-AKN*Y*RHOI*CMP*SHZ*AJMD(m)
            HY=HY1+HY2
            HZ=-AKN*CMP*CHZ*AJM(m)
            
            L=n+5*m
            FX=FX+HX*TSS(L,K)
            FY=FY+HY*TSS(L,K)
            FZ=FZ+HZ*TSS(L,K)
 4       CONTINUE
 2    CONTINUE
      RETURN
      END

C==========================================================================
C     July 2017, G.K.Stephens, this routine updated to incorporate Jay Albert's
C     improvements to the Bessel function evaluation. The Bessel function
C     values are now precomputed and  passed into this rather than computed 
C     inside this routine.
      SUBROUTINE TAILSHT_OE (IEVO,MK,M,X,Y,Z,AJM,AJMD,BX,BY,BZ)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /TAIL/ D0  ! THE COMMON BLOCKS FORWARDS TAIL SHEET THICKNESS
      DIMENSION AJM(0:5),AJMD(0:5)
C-----------------------------------------------------------------------------------
C
      RNOT=20.0     !    Rho_0 - scale parameter along the tail axis
      DLTK=1.0      !    step in Km
c
c -----------------------------------------------------------------------------------

      RHO=DSQRT(X*X+Y*Y)
C
      CSPHI=X/RHO
      SNPHI=Y/RHO
C
      phi=DATAN2(Y,X)
      CSMPHI=DCOS(m*phi)
      SNMPHI=DSIN(m*phi)
C
      DKM=1.D0+(MK-1)*DLTK
      RKM=DKM/RNOT
C
      RKMZ=RKM*Z
      RKMR=RKM*RHO
C
      ZD=DSQRT(Z*Z+D0*D0)
C
      REX=DEXP(RKM*ZD)
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
      if(ievo.eq.0) then
c -----------------------------------------
c calculating symmetric modes
c -----------------------------------------
c
      BRO=-M*SNMPHI*Z*AJMD(m)/ZD/REX
      BPHI=-M*M*CSMPHI*Z*AJM(m)/RKMR/ZD/REX
      BZ=M*SNMPHI*AJM(m)/REX
c
c -----------------------------------------
      else
c -----------------------------------------
c calculating asymmetric modes
c -----------------------------------------
c
      BRO=M*CSMPHI*Z*AJMD(m)/ZD/REX
      BPHI=-M*M*SNMPHI*Z*AJM(m)/RKMR/ZD/REX
      BZ=-M*CSMPHI*AJM(m)/REX
c
c -----------------------------------------
      end if
c
c --- transformation from cylindrical ccordinates to GSM ---
c
      BX=BRO*CSPHI-BPHI*SNPHI
      BY=BRO*SNPHI+BPHI*CSPHI
C
C    CALCULATION OF THE MAGNETOTAIL CURRENT CONTRIBUTION IS FINISHED
C
      RETURN
      END

C==========================================================================
C     July 2017, G.K.Stephens, this routine updated to incorporate Jay Albert's
C     improvements to the Bessel function evaluation.
      SUBROUTINE  SHTBNORM_O(K,L,X,Y,Z,FX,FY,FZ) ! modified SHTBNORM_O
      IMPLICIT  REAL * 8  (A - H, O - Z)
      DIMENSION AK(5)
      DIMENSION AJM(0:14),AJMD(0:14)
      COMMON /TSO/ TSO(80,5,4)
      
      AK(1)=TSO(76,K,L)
      AK(2)=TSO(77,K,L)
      AK(3)=TSO(78,K,L)
      AK(4)=TSO(79,K,L)
      AK(5)=TSO(80,K,L)
      
      phi=DATAN2(Y,X)
      RHO=dsqrt(X*X+Y*Y)
      if(RHO.lt.1.D-8) then
         RHOI=1.D8
      else
         RHOI=1.D0/RHO
      endif
      DPDX=-Y*RHOI*RHOI
      DPDY= X*RHOI*RHOI

      FX=0.D0
      FY=0.D0
      FZ=0.D0
      
      DO 2 n=1,5
         AKN=dabs(AK(n))
         AKNR=AKN*RHO
         if(AKNR.lt.1.D-8) then
            AKNRI=1.D8
         else
            AKNRI=1.D0/AKNR
         endif
         
         CHZ=dcosh(Z*AKN)
         SHZ=dsinh(Z*AKN)

         call bessJJ(14,AKNR, AJM) !!! get all n in one call
         DO 3 m=1,14
            AJMD(m)=AJM(m-1)-m*AJM(m)*AKNRI
 3       CONTINUE
         AJMD(0)=-AJM(1)
         
         DO 4 m=0,14
            CMP=DCOS(m*phi)
            SMP=DSIN(m*phi)

            HX1=m*DPDX*SMP*SHZ*AJM(m)
            HX2=-AKN*X*RHOI*CMP*SHZ*AJMD(m)
            HX=HX1+HX2
            HY1=m*DPDY*SMP*SHZ*AJM(m)
            HY2=-AKN*Y*RHOI*CMP*SHZ*AJMD(m)
            HY=HY1+HY2
            HZ=-AKN*CMP*CHZ*AJM(m)

            L1=n+5*m
            FX=FX+HX*TSO(L1,K,L)
            FY=FY+HY*TSO(L1,K,L)
            FZ=FZ+HZ*TSO(L1,K,L)
 4       CONTINUE
 2    CONTINUE
      RETURN
      END

C==========================================================================
C     July 2017, G.K.Stephens, this routine updated to incorporate Jay Albert's
C     improvements to the Bessel function evaluation.
      SUBROUTINE  SHTBNORM_E(K,L,X,Y,Z,FX,FY,FZ) ! modified SHTBNORM_E
      IMPLICIT  REAL * 8  (A - H, O - Z)
      DIMENSION AK(5)
      DIMENSION AJM(0:14),AJMD(0:14)
      COMMON /TSE/ TSE(80,5,4)
      
      AK(1)=TSE(76,K,L)
      AK(2)=TSE(77,K,L)
      AK(3)=TSE(78,K,L)
      AK(4)=TSE(79,K,L)
      AK(5)=TSE(80,K,L)

      phi=DATAN2(Y,X)
      RHO=dsqrt(X*X+Y*Y)
      if(RHO.lt.1.D-8) then
         RHOI=1.D8
      else
         RHOI=1.D0/RHO
      endif
      DPDX=-Y*RHOI*RHOI
      DPDY= X*RHOI*RHOI
      
      FX=0.D0
      FY=0.D0
      FZ=0.D0

      DO 2 n=1,5
         AKN=dabs(AK(n))
         AKNR=AKN*RHO
         if(AKNR.lt.1.D-8) then
            AKNRI=1.D8
         else
            AKNRI=1.D0/AKNR
         endif
         
         CHZ=dcosh(Z*AKN)
         SHZ=dsinh(Z*AKN)
         
         call bessJJ(14,AKNR, AJM) !!! get all n in one call
         DO 3 m=1,14
            AJMD(m)=AJM(m-1)-m*AJM(m)*AKNRI
 3       CONTINUE
         AJMD(0)=-AJM(1)

         DO 4 m=0,14
            CMP=DCOS(m*phi)
            SMP=DSIN(m*phi)

            HX1=-m*DPDX*CMP*SHZ*AJM(m)
            HX2=-AKN*X*RHOI*SMP*SHZ*AJMD(m)
            HX=HX1+HX2
            HY1=-m*DPDY*CMP*SHZ*AJM(m)
            HY2=-AKN*Y*RHOI*SMP*SHZ*AJMD(m)
            HY=HY1+HY2
            HZ=-AKN*SMP*CHZ*AJM(m)
            
            L1=n+5*m
            FX=FX+HX*TSE(L1,K,L)
            FY=FY+HY*TSE(L1,K,L)
            FZ=FZ+HZ*TSE(L1,K,L)
 4       CONTINUE
 2    CONTINUE
      RETURN
      END

C==========================================================================
C     July 2017, G.K.Stephens, this routine updated to incorporate Jay Albert's
C     improvements to the Bessel function evaluation.
      SUBROUTINE bessJJ(n,x, bessJ)
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension bessJ(0:n) ! bessJ holds J0 to Jn
      PARAMETER (IACC=40,BIGNO=1.D10,BIGNI=1.D-10)
      logical iseven
      ax=Dabs(x)

      tox=2.D0/ax
c     start at some large m, larger than the desired n, multiply by 2 to ensure
c     m starts at an even number
      m=2*((n+int(Dsqrt(Dfloat(IACC*n))))/2)

      evnsum=0.D0! keeps track of the sum of the even Js (J0+J2+J4+...)
      iseven = .false.

c     we set the value of Jm to some arbitrary value, here Jm=1, after the loop
c     is done, the values will be normalized using the sum
      bjp=0.D0
      bj =1.D0

c     initialize to zero
      do i=0,n
         bessJ(i)=0.
      enddo

      do 12 j=m,1,-1

c     the next value int the recursion relation J_n-1 = (2*n/x)*Jn - J_n+1
         bjm=j*tox*bj-bjp
         bjp=bj! decrement so shift J_n+1 ot Jn
         bj=bjm! decrement so shift J_n ot J_n-1

c     if the value gets too large, shift the decimal of everything by 10 places
         if (Dabs(bj).gt.BIGNO) then
            bj =bj *BIGNI
            bjp=bjp*BIGNI
            evnsum=evnsum*BIGNI
	    do i=j+1,n
	       bessJ(i)=bessJ(i)*BIGNI
            enddo
         endif

         if(iseven)evnsum=evnsum+bj ! only sum over the even Jns
         iseven=.not.iseven

         if(j.le.n) bessJ(j)=bjp ! Jj(x)

12    continue

c     sum is currently the sum of all the evens
c     use Miller's algorithm for Bessel functions which uses the identity: 
c     1.0 = 2.0*sum(J_evens) - J0, thus the quantity (2.0*sum(J_evens) - J0)
c     is used as a normalization factor
      bnorm=2.D0*evnsum-bj

c     normalize the Bessel functions
      do i=1,n
         bessJ(i)=bessJ(i)/bnorm
      enddo
      bessJ(0)=bj/bnorm ! J0(x)

c     Apply Jn(-x) = (-1)^n * Jn(x)
      if (x .lt. 0.D0) then
         do i=1,n,2
            bessJ(i)=-bessJ(i)
         enddo
      endif

      return
      END


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      SUBROUTINE BIRK_TOT (PS,X,Y,Z,BX11,BY11,BZ11,BX12,BY12,BZ12,
     *                          BX21,BY21,BZ21,BX22,BY22,BZ22)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SH11(86),SH12(86),SH21(86),SH22(86)
      COMMON /BIRKPAR/ XKAPPA1,XKAPPA2   !  INPUT PARAMETERS, SPECIFIED FROM A MAIN PROGRAM
      COMMON /DPHI_B_RHO0/ DPHI,B,RHO_0,XKAPPA ! PARAMETERS, CONTROL DAY-NIGHT ASYMMETRY OF F.A.C.

      DATA SH11/46488.84663,-15541.95244,-23210.09824,-32625.03856,
     *-109894.4551,-71415.32808,58168.94612,55564.87578,-22890.60626,
     *-6056.763968,5091.368100,239.7001538,-13899.49253,4648.016991,
     *6971.310672,9699.351891,32633.34599,21028.48811,-17395.96190,
     *-16461.11037,7447.621471,2528.844345,-1934.094784,-588.3108359,
     *-32588.88216,10894.11453,16238.25044,22925.60557,77251.11274,
     *50375.97787,-40763.78048,-39088.60660,15546.53559,3559.617561,
     *-3187.730438,309.1487975,88.22153914,-243.0721938,-63.63543051,
     *191.1109142,69.94451996,-187.9539415,-49.89923833,104.0902848,
     *-120.2459738,253.5572433,89.25456949,-205.6516252,-44.93654156,
     *124.7026309,32.53005523,-98.85321751,-36.51904756,98.88241690,
     *24.88493459,-55.04058524,61.14493565,-128.4224895,-45.35023460,
     *105.0548704,-43.66748755,119.3284161,31.38442798,-92.87946767,
     *-33.52716686,89.98992001,25.87341323,-48.86305045,59.69362881,
     *-126.5353789,-44.39474251,101.5196856,59.41537992,41.18892281,
     *80.86101200,3.066809418,7.893523804,30.56212082,10.36861082,
     *8.222335945,19.97575641,2.050148531,4.992657093,2.300564232,
     *.2256245602,-.05841594319/

      DATA SH12/210260.4816,-1443587.401,-1468919.281,281939.2993,
     *-1131124.839,729331.7943,2573541.307,304616.7457,468887.5847,
     *181554.7517,-1300722.650,-257012.8601,645888.8041,-2048126.412,
     *-2529093.041,571093.7972,-2115508.353,1122035.951,4489168.802,
     *75234.22743,823905.6909,147926.6121,-2276322.876,-155528.5992,
     *-858076.2979,3474422.388,3986279.931,-834613.9747,3250625.781,
     *-1818680.377,-7040468.986,-414359.6073,-1295117.666,-346320.6487,
     *3565527.409,430091.9496,-.1565573462,7.377619826,.4115646037,
     *-6.146078880,3.808028815,-.5232034932,1.454841807,-12.32274869,
     *-4.466974237,-2.941184626,-.6172620658,12.64613490,1.494922012,
     *-21.35489898,-1.652256960,16.81799898,-1.404079922,-24.09369677,
     *-10.99900839,45.94237820,2.248579894,31.91234041,7.575026816,
     *-45.80833339,-1.507664976,14.60016998,1.348516288,-11.05980247,
     *-5.402866968,31.69094514,12.28261196,-37.55354174,4.155626879,
     *-33.70159657,-8.437907434,36.22672602,145.0262164,70.73187036,
     *85.51110098,21.47490989,24.34554406,31.34405345,4.655207476,
     *5.747889264,7.802304187,1.844169801,4.867254550,2.941393119,
     *.1379899178,.06607020029/

      DATA SH21/162294.6224,503885.1125,-27057.67122,-531450.1339,
     *84747.05678,-237142.1712,84133.61490,259530.0402,69196.05160,
     *-189093.5264,-19278.55134,195724.5034,-263082.6367,-818899.6923,
     *43061.10073,863506.6932,-139707.9428,389984.8850,-135167.5555,
     *-426286.9206,-109504.0387,295258.3531,30415.07087,-305502.9405,
     *100785.3400,315010.9567,-15999.50673,-332052.2548,54964.34639,
     *-152808.3750,51024.67566,166720.0603,40389.67945,-106257.7272,
     *-11126.14442,109876.2047,2.978695024,558.6019011,2.685592939,
     *-338.0004730,-81.99724090,-444.1102659,89.44617716,212.0849592,
     *-32.58562625,-982.7336105,-35.10860935,567.8931751,-1.917212423,
     *-260.2023543,-1.023821735,157.5533477,23.00200055,232.0603673,
     *-36.79100036,-111.9110936,18.05429984,447.0481000,15.10187415,
     *-258.7297813,-1.032340149,-298.6402478,-1.676201415,180.5856487,
     *64.52313024,209.0160857,-53.85574010,-98.52164290,14.35891214,
     *536.7666279,20.09318806,-309.7349530,58.54144539,67.45226850,
     *97.92374406,4.752449760,10.46824379,32.91856110,12.05124381,
     *9.962933904,15.91258637,1.804233877,6.578149088,2.515223491,
     *.1930034238,-.02261109942/

      DATA SH22/-131287.8986,-631927.6885,-318797.4173,616785.8782,
     *-50027.36189,863099.9833,47680.20240,-1053367.944,-501120.3811,
     *-174400.9476,222328.6873,333551.7374,-389338.7841,-1995527.467,
     *-982971.3024,1960434.268,297239.7137,2676525.168,-147113.4775,
     *-3358059.979,-2106979.191,-462827.1322,1017607.960,1039018.475,
     *520266.9296,2627427.473,1301981.763,-2577171.706,-238071.9956,
     *-3539781.111,94628.16420,4411304.724,2598205.733,637504.9351,
     *-1234794.298,-1372562.403,-2.646186796,-31.10055575,2.295799273,
     *19.20203279,30.01931202,-302.1028550,-14.78310655,162.1561899,
     *.4943938056,176.8089129,-.2444921680,-100.6148929,9.172262228,
     *137.4303440,-8.451613443,-84.20684224,-167.3354083,1321.830393,
     *76.89928813,-705.7586223,18.28186732,-770.1665162,-9.084224422,
     *436.3368157,-6.374255638,-107.2730177,6.080451222,65.53843753,
     *143.2872994,-1028.009017,-64.22739330,547.8536586,-20.58928632,
     *597.3893669,10.17964133,-337.7800252,159.3532209,76.34445954,
     *84.74398828,12.76722651,27.63870691,32.69873634,5.145153451,
     *6.310949163,6.996159733,1.971629939,4.436299219,2.904964304,
     *.1486276863,.06859991529/

C ====   LEAST SQUARES FITTING ONLY:
C       BX11=0.D0
C       BY11=0.D0
C       BZ11=0.D0
C       BX12=0.D0
C       BY12=0.D0
C       BZ12=0.D0
C       BX21=0.D0
C       BY21=0.D0
C       BZ21=0.D0
C       BX22=0.D0
C       BY22=0.D0
C       BZ22=0.D0
C===================================

      XKAPPA=XKAPPA1        !  FORWARDED IN BIRK_1N2
      X_SC=XKAPPA1-1.1D0    !  FORWARDED IN BIRK_SHL

      CALL BIRK_1N2 (1,1,PS,X,Y,Z,FX11,FY11,FZ11)           !  REGION 1, MODE 1
      CALL BIRK_SHL (SH11,PS,X_SC,X,Y,Z,HX11,HY11,HZ11)
      BX11=FX11+HX11
      BY11=FY11+HY11
      BZ11=FZ11+HZ11

      CALL BIRK_1N2 (1,2,PS,X,Y,Z,FX12,FY12,FZ12)           !  REGION 1, MODE 2
      CALL BIRK_SHL (SH12,PS,X_SC,X,Y,Z,HX12,HY12,HZ12)
      BX12=FX12+HX12
      BY12=FY12+HY12
      BZ12=FZ12+HZ12

      XKAPPA=XKAPPA2        !  FORWARDED IN BIRK_1N2
      X_SC=XKAPPA2-1.0D0    !  FORWARDED IN BIRK_SHL

      CALL BIRK_1N2 (2,1,PS,X,Y,Z,FX21,FY21,FZ21)           !  REGION 2, MODE 1
      CALL BIRK_SHL (SH21,PS,X_SC,X,Y,Z,HX21,HY21,HZ21)
      BX21=FX21+HX21
      BY21=FY21+HY21
      BZ21=FZ21+HZ21

      CALL BIRK_1N2 (2,2,PS,X,Y,Z,FX22,FY22,FZ22)           !  REGION 2, MODE 2
      CALL BIRK_SHL (SH22,PS,X_SC,X,Y,Z,HX22,HY22,HZ22)
      BX22=FX22+HX22
      BY22=FY22+HY22
      BZ22=FZ22+HZ22

      RETURN
      END
C
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
      SUBROUTINE BIRK_1N2 (NUMB,MODE,PS,X,Y,Z,BX,BY,BZ)        !   SEE NB# 6, P.60
C
C  CALCULATES COMPONENTS  OF REGION 1/2 FIELD IN SPHERICAL COORDS.  DERIVED FROM THE S/R DIPDEF2C (WHICH
C    DOES THE SAME JOB, BUT INPUT/OUTPUT THERE WAS IN SPHERICAL COORDS, WHILE HERE WE USE CARTESIAN ONES)
C
C   INPUT:  NUMB=1 (2) FOR REGION 1 (2) CURRENTS
C           MODE=1 YIELDS SIMPLE SINUSOIDAL MLT VARIATION, WITH MAXIMUM CURRENT AT DAWN/DUSK MERIDIAN
C     WHILE MODE=2 YIELDS THE SECOND HARMONIC.
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A11(31),A12(31),A21(31),A22(31)
      COMMON /MODENUM/ M
      COMMON /DTHETA/ DTHETA

      COMMON /DPHI_B_RHO0/ DPHI,B,RHO_0,XKAPPA ! THESE PARAMETERS CONTROL DAY-NIGHT ASYMMETRY OF F.A.C., AS FOLLOWS:

C  (1) DPHI:   HALF-DIFFERENCE (IN RADIANS) BETWEEN DAY AND NIGHT LATITUDE OF FAC OVAL AT IONOSPHERIC ALTITUDE;
C              TYPICAL VALUE: 0.06
C  (2) B:      AN ASYMMETRY FACTOR AT HIGH-ALTITUDES;  FOR B=0, THE ONLY ASYMMETRY IS THAT FROM DPHI
C              TYPICAL VALUES: 0.35-0.70
C  (3) RHO_0:  A FIXED PARAMETER, DEFINING THE DISTANCE RHO, AT WHICH THE LATITUDE SHIFT GRADUALLY SATURATES AND
C              STOPS INCREASING
C              ITS VALUE WAS ASSUMED FIXED, EQUAL TO 7.0.
C  (4) XKAPPA: AN OVERALL SCALING FACTOR, WHICH CAN BE USED FOR CHANGING THE SIZE OF THE F.A.C. OVAL
C
      DATA BETA,RH,EPS/0.9D0,10.D0,3.D0/ ! parameters of the tilt-dependent deformation of the untilted F.A.C. field

      DATA A11/.1618068350,-.1797957553,2.999642482,-.9322708978,
     *-.6811059760,.2099057262,-8.358815746,-14.86033550,.3838362986,
     *-16.30945494,4.537022847,2.685836007,27.97833029,6.330871059,
     *1.876532361,18.95619213,.9651528100,.4217195118,-.08957770020,
     *-1.823555887,.7457045438,-.5785916524,-1.010200918,.01112389357,
     *.09572927448,-.3599292276,8.713700514,.9763932955,3.834602998,
     *2.492118385,.7113544659/
      DATA A12/.7058026940,-.2845938535,5.715471266,-2.472820880,
     *-.7738802408,.3478293930,-11.37653694,-38.64768867,.6932927651,
     *-212.4017288,4.944204937,3.071270411,33.05882281,7.387533799,
     *2.366769108,79.22572682,.6154290178,.5592050551,-.1796585105,
     *-1.654932210,.7309108776,-.4926292779,-1.130266095,-.009613974555,
     *.1484586169,-.2215347198,7.883592948,.02768251655,2.950280953,
     *1.212634762,.5567714182/
      DATA A21/.1278764024,-.2320034273,1.805623266,-32.37241440,
     *-.9931490648,.3175085630,-2.492465814,-16.21600096,.2695393416,
     *-6.752691265,3.971794901,14.54477563,41.10158386,7.912889730,
     *1.258297372,9.583547721,1.014141963,.5104134759,-.1790430468,
     *-1.756358428,.7561986717,-.6775248254,-.04014016420,.01446794851,
     *.1200521731,-.2203584559,4.508963850,.8221623576,1.779933730,
     *1.102649543,.8867880020/
      DATA A22/.4036015198,-.3302974212,2.827730930,-45.44405830,
     *-1.611103927,.4927112073,-.003258457559,-49.59014949,.3796217108,
     *-233.7884098,4.312666980,18.05051709,28.95320323,11.09948019,
     *.7471649558,67.10246193,.5667096597,.6468519751,-.1560665317,
     *-1.460805289,.7719653528,-.6658988668,.2515179349E-05,
     *.02426021891,.1195003324,-.2625739255,4.377172556,.2421190547,
     *2.503482679,1.071587299,.7247997430/

      B=0.5
      RHO_0=7.0

      M=MODE
      IF (NUMB.EQ.1) THEN
          DPHI=0.055D0
          DTHETA=0.06D0
      ENDIF

      IF (NUMB.EQ.2) THEN
          DPHI=0.030D0
          DTHETA=0.09D0
      ENDIF

      Xsc=X*XKAPPA
      Ysc=Y*XKAPPA
      Zsc=Z*XKAPPA
      RHO=DSQRT(Xsc**2+Zsc**2)

      Rsc=DSQRT(Xsc**2+Ysc**2+Zsc**2)                                 !  SCALED
      RHO2=RHO_0**2

      IF (Xsc.EQ.0.D0.AND.Zsc.EQ.0.D0) THEN
         PHI=0.D0
      ELSE
         PHI=DATAN2(-Zsc,Xsc)  !  FROM CARTESIAN TO CYLINDRICAL (RHO,PHI,Y)
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
      DPHISDY= BETA*PS*R1RH**(EPS-1.D0)*Ysc/(RH*Rsc*
     *   (1.D0+R1RH**EPS)**(1.D0/EPS+1.D0))

      SPHICS=DSIN(PHIS)
      CPHICS=DCOS(PHIS)

      XS= RHO*CPHICS
      ZS=-RHO*SPHICS

      IF (NUMB.EQ.1) THEN
        IF (MODE.EQ.1) CALL TWOCONES (A11,XS,Ysc,ZS,BXS,BYAS,BZS)
        IF (MODE.EQ.2) CALL TWOCONES (A12,XS,Ysc,ZS,BXS,BYAS,BZS)
      ELSE
        IF (MODE.EQ.1) CALL TWOCONES (A21,XS,Ysc,ZS,BXS,BYAS,BZS)
        IF (MODE.EQ.2) CALL TWOCONES (A22,XS,Ysc,ZS,BXS,BYAS,BZS)
      ENDIF

      BRHOAS=BXS*CPHICS-BZS*SPHICS
      BPHIAS=-BXS*SPHICS-BZS*CPHICS

      BRHO_S=BRHOAS*DPHISPHI                             *XKAPPA        ! SCALING
      BPHI_S=(BPHIAS-RHO*(BYAS*DPHISDY+BRHOAS*DPHISRHO)) *XKAPPA
      BY_S=BYAS*DPHISPHI                                 *XKAPPA

      BX=BRHO_S*CPHIC-BPHI_S*SPHIC
      BY=BY_S
      BZ=-BRHO_S*SPHIC-BPHI_S*CPHIC

      RETURN
      END
c
C=========================================================================
c
      SUBROUTINE TWOCONES (A,X,Y,Z,BX,BY,BZ)
C
C   ADDS FIELDS FROM TWO CONES (NORTHERN AND SOUTHERN), WITH A PROPER SYMMETRY OF THE CURRENT AND FIELD,
C     CORRESPONDING TO THE REGION 1 BIRKELAND CURRENTS. (SEE NB #6, P.58).
C

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(31)

      CALL ONE_CONE (A,X,Y,Z,BXN,BYN,BZN)
      CALL ONE_CONE (A,X,-Y,-Z,BXS,BYS,BZS)
      BX=BXN-BXS
      BY=BYN+BYS
      BZ=BZN+BZS

      RETURN
      END
c
C-------------------------------------------------------------------------
C
      SUBROUTINE ONE_CONE(A,X,Y,Z,BX,BY,BZ)
c
c  RETURNS FIELD COMPONENTS FOR A DEFORMED CONICAL CURRENT SYSTEM, FITTED TO A BIOSAVART FIELD
c    BY SIM_14.FOR.  HERE ONLY THE NORTHERN CONE IS TAKEN INTO ACCOUNT.
c

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(31)

      COMMON /DTHETA/ DTHETA
      COMMON /MODENUM/ M

      DATA DR,DT/1.D-6,1.D-6/  !   JUST FOR NUMERICAL DIFFERENTIATION

      THETA0=A(31)

      RHO2=X**2+Y**2
      RHO=DSQRT(RHO2)
      R=DSQRT(RHO2+Z**2)
      THETA=DATAN2(RHO,Z)
      PHI=DATAN2(Y,X)
C
C   MAKE THE DEFORMATION OF COORDINATES:
C
       RS=R_S(A,R,THETA)
       THETAS=THETA_S(A,R,THETA)
       PHIS=PHI
C
C   CALCULATE FIELD COMPONENTS AT THE NEW POSITION (ASTERISKED):
C
       CALL FIALCOS (RS,THETAS,PHIS,BTAST,BFAST,M,THETA0,DTHETA)    !   MODE #M
C
C   NOW TRANSFORM B{R,T,F}_AST BY THE DEFORMATION TENSOR:
C
C      FIRST OF ALL, FIND THE DERIVATIVES:
C
       DRSDR=(R_S(A,R+DR,THETA)-R_S(A,R-DR,THETA))/(2.D0*DR)
       DRSDT=(R_S(A,R,THETA+DT)-R_S(A,R,THETA-DT))/(2.D0*DT)
       DTSDR=(THETA_S(A,R+DR,THETA)-THETA_S(A,R-DR,THETA))/(2.D0*DR)
       DTSDT=(THETA_S(A,R,THETA+DT)-THETA_S(A,R,THETA-DT))/(2.D0*DT)
       STSST=DSIN(THETAS)/DSIN(THETA)
       RSR=RS/R

       BR     =-RSR/R*STSST*BTAST*DRSDT                 !   NB#6, P.43    BRAST DOES NOT ENTER HERE
       BTHETA = RSR*STSST*BTAST*DRSDR                  !               (SINCE IT IS ZERO IN OUR CASE)
       BPHI   = RSR*BFAST*(DRSDR*DTSDT-DRSDT*DTSDR)

       S=RHO/R
       C=Z/R
       SF=Y/RHO
       CF=X/RHO

       BE=BR*S+BTHETA*C

       BX=A(1)*(BE*CF-BPHI*SF)
       BY=A(1)*(BE*SF+BPHI*CF)
       BZ=A(1)*(BR*C-BTHETA*S)

       RETURN
       END
C
C=====================================================================================
      DOUBLE PRECISION FUNCTION R_S(A,R,THETA)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(31)
C
      R_S=R+A(2)/R+A(3)*R/DSQRT(R**2+A(11)**2)+A(4)*R/(R**2+A(12)**2)
     *+(A(5)+A(6)/R+A(7)*R/DSQRT(R**2+A(13)**2)+A(8)*R/(R**2+A(14)**2))*
     * DCOS(THETA)
     *+(A(9)*R/DSQRT(R**2+A(15)**2)+A(10)*R/(R**2+A(16)**2)**2)
     * *DCOS(2.D0*THETA)
C
      RETURN
      END
C
C-----------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION THETA_S(A,R,THETA)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(31)
c
      THETA_S=THETA+(A(17)+A(18)/R+A(19)/R**2
     *                +A(20)*R/DSQRT(R**2+A(27)**2))*DSIN(THETA)
     * +(A(21)+A(22)*R/DSQRT(R**2+A(28)**2)
     *                +A(23)*R/(R**2+A(29)**2))*DSIN(2.D0*THETA)
     * +(A(24)+A(25)/R+A(26)*R/(R**2+A(30)**2))*DSIN(3.D0*THETA)
C
      RETURN
      END
C
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
      SUBROUTINE FIALCOS(R,THETA,PHI,BTHETA,BPHI,N,THETA0,DT)
C
C  CONICAL MODEL OF BIRKELAND CURRENT FIELD; BASED ON THE OLD S/R FIALCO (OF 1990-91)
C  SEE THE OLD NOTEBOOK 1985-86-88, NOTE OF MARCH 5, BUT HERE BOTH INPUT AND OUTPUT ARE IN SPHERICAL CDS.

C  BTN, AND BPN ARE THE ARRAYS OF BTHETA AND BPHI (BTN(i), BPN(i) CORRESPOND TO i-th MODE).
C   ONLY FIRST  N  MODE AMPLITUDES ARE COMPUTED (N<=10).
C    THETA0 IS THE ANGULAR HALF-WIDTH OF THE CONE, DT IS THE ANGULAR H.-W. OF THE CURRENT LAYER

C   NOTE:  BR=0  (BECAUSE ONLY RADIAL CURRENTS ARE PRESENT IN THIS MODEL)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION  BTN(10),BPN(10),CCOS(10),SSIN(10)

      SINTE=DSIN(THETA)
      RO=R*SINTE
      COSTE=DCOS(THETA)
      SINFI=DSIN(PHI)
      COSFI=DCOS(PHI)
      TG=SINTE/(1.D0+COSTE)   !        TAN(THETA/2)
      CTG=SINTE/(1.D0-COSTE)  !        CTG(THETA/2)
C
C
      TETANP=THETA0+DT
      TETANM=THETA0-DT
      IF(THETA.LT.TETANM) GOTO 1
      TGP=DTAN(TETANP*0.5D0)
      TGM=DTAN(TETANM*0.5D0)
      TGM2=TGM*TGM
      TGP2=TGP*TGP
  1   CONTINUE

      COSM1=1.D0
      SINM1=0.D0
      TM=1.D0
      TGM2M=1.D0
      TGP2M=1.D0

      DO 2 M=1,N
      TM=TM*TG
      CCOS(M)=COSM1*COSFI-SINM1*SINFI
      SSIN(M)=SINM1*COSFI+COSM1*SINFI
      COSM1=CCOS(M)
      SINM1=SSIN(M)
      IF(THETA.LT.TETANM) THEN
      T=TM
      DTT=0.5D0*M*TM*(TG+CTG)
      DTT0=0.D0
      ELSE IF(THETA.LT.TETANP) THEN
      TGM2M=TGM2M*TGM2
      FC=1.D0/(TGP-TGM)
      FC1=1.D0/(2*M+1)
      TGM2M1=TGM2M*TGM
      TG21=1.D0+TG*TG
      T=FC*(TM*(TGP-TG)+FC1*(TM*TG-TGM2M1/TM))
      DTT=0.5D0*M*FC*TG21*(TM/TG*(TGP-TG)-FC1*(TM-TGM2M1/(TM*TG)))
      DTT0=0.5D0*FC*((TGP+TGM)*(TM*TG-FC1*(TM*TG-TGM2M1/TM))+
     * TM*(1.D0-TGP*TGM)-(1.D0+TGM2)*TGM2M/TM)
      ELSE
      TGP2M=TGP2M*TGP2
      TGM2M=TGM2M*TGM2
      FC=1.D0/(TGP-TGM)
      FC1=1.D0/(2*M+1)
      T=FC*FC1*(TGP2M*TGP-TGM2M*TGM)/TM
      DTT=-T*M*0.5D0*(TG+CTG)
      ENDIF

      BTN(M)=M*T*CCOS(M)/RO
  2   BPN(M)=-DTT*SSIN(M)/R

      BTHETA=BTN(N) *800.
      BPHI  =BPN(N) *800.

      RETURN
      END
C
C-------------------------------------------------------------------------
C
C
         SUBROUTINE BIRK_SHL (A,PS,X_SC,X,Y,Z,BX,BY,BZ)
C
         IMPLICIT  REAL * 8  (A - H, O - Z)
         DIMENSION A(86)
C
         CPS=DCOS(PS)
         SPS=DSIN(PS)

         S3PS=2.D0*CPS
C
         PST1=PS*A(85)
         PST2=PS*A(86)

         ST1=DSIN(PST1)
         CT1=DCOS(PST1)
         ST2=DSIN(PST2)
         CT2=DCOS(PST2)

         X1=X*CT1-Z*ST1
         Z1=X*ST1+Z*CT1
         X2=X*CT2-Z*ST2
         Z2=X*ST2+Z*CT2
C
         L=0
         GX=0.D0
         GY=0.D0
         GZ=0.D0
C
         DO 1 M=1,2     !    M=1 IS FOR THE 1ST SUM ("PERP." SYMMETRY)
C                          AND M=2 IS FOR THE SECOND SUM ("PARALL." SYMMETRY)
             DO 2 I=1,3
                  P=A(72+I)
                  Q=A(78+I)
                  CYPI=DCOS(Y/P)
                  CYQI=DCOS(Y/Q)
                  SYPI=DSIN(Y/P)
                  SYQI=DSIN(Y/Q)
C
                DO 3 K=1,3
                   R=A(75+K)
                   S=A(81+K)
                   SZRK=DSIN(Z1/R)
                   CZSK=DCOS(Z2/S)
                   CZRK=DCOS(Z1/R)
                   SZSK=DSIN(Z2/S)
                     SQPR=DSQRT(1.D0/P**2+1.D0/R**2)
                     SQQS=DSQRT(1.D0/Q**2+1.D0/S**2)
                        EPR=DEXP(X1*SQPR)
                        EQS=DEXP(X2*SQQS)
C
                  DO 4 N=1,2  ! N=1 IS FOR THE FIRST PART OF EACH COEFFICIENT
C                                AND N=2 IS FOR THE SECOND ONE

                    DO 5 NN=1,2 !   NN = 1,2 FURTHER SPLITS THE COEFFICIENTS INTO 2 PARTS,
C                                         TO TAKE INTO ACCOUNT THE SCALE FACTOR DEPENDENCE

                    IF (M.EQ.1) THEN
                         FX=-SQPR*EPR*CYPI*SZRK
                         FY=EPR*SYPI*SZRK/P
                         FZ=-EPR*CYPI*CZRK/R
                       IF (N.EQ.1) THEN
                         IF (NN.EQ.1) THEN
                          HX=FX
                          HY=FY
                          HZ=FZ
                         ELSE
                          HX=FX*X_SC
                          HY=FY*X_SC
                          HZ=FZ*X_SC
                         ENDIF
                       ELSE
                         IF (NN.EQ.1) THEN
                          HX=FX*CPS
                          HY=FY*CPS
                          HZ=FZ*CPS
                         ELSE
                          HX=FX*CPS*X_SC
                          HY=FY*CPS*X_SC
                          HZ=FZ*CPS*X_SC
                         ENDIF
                       ENDIF

                     ELSE                            !   M.EQ.2
                         FX=-SPS*SQQS*EQS*CYQI*CZSK
                         FY=SPS/Q*EQS*SYQI*CZSK
                         FZ=SPS/S*EQS*CYQI*SZSK
                       IF (N.EQ.1) THEN
                        IF (NN.EQ.1) THEN
                          HX=FX
                          HY=FY
                          HZ=FZ
                        ELSE
                          HX=FX*X_SC
                          HY=FY*X_SC
                          HZ=FZ*X_SC
                        ENDIF
                       ELSE
                        IF (NN.EQ.1) THEN
                         HX=FX*S3PS
                         HY=FY*S3PS
                         HZ=FZ*S3PS
                        ELSE
                         HX=FX*S3PS*X_SC
                         HY=FY*S3PS*X_SC
                         HZ=FZ*S3PS*X_SC
                        ENDIF
                       ENDIF
                  ENDIF
       L=L+1

       IF (M.EQ.1) THEN
       HXR=HX*CT1+HZ*ST1
       HZR=-HX*ST1+HZ*CT1
       ELSE
       HXR=HX*CT2+HZ*ST2
       HZR=-HX*ST2+HZ*CT2
       ENDIF

       GX=GX+HXR*A(L)
       GY=GY+HY *A(L)
  5    GZ=GZ+HZR*A(L)

  4   CONTINUE
  3   CONTINUE
  2   CONTINUE
  1   CONTINUE

      BX=GX
      BY=GY
      BZ=GZ

      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
      SUBROUTINE BIRTOTSY (PS,X,Y,Z,BX11,BY11,BZ11,BX12,BY12,BZ12,
     *                          BX21,BY21,BZ21,BX22,BY22,BZ22)
C
C   THIS S/R IS ALMOST IDENTICAL TO BIRK_TOT, BUT IT IS FOR THE SYMMETRIC MODE, IN WHICH
C     J_parallel IS AN EVEN FUNCTION OF Ygsm.
C
C
C      IOPBS -  BIRKELAND FIELD MODE FLAG:
C         IOPBS=0 - ALL COMPONENTS
C         IOPBS=1 - REGION 1, MODES 1 & 2 (SYMMETRIC !)
C         IOPBS=2 - REGION 2, MODES 1 & 2 (SYMMETRIC !)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SH11(86),SH12(86),SH21(86),SH22(86)
      COMMON /BIRKPAR/ XKAPPA1,XKAPPA2   !  INPUT PARAMETERS, SPECIFIED FROM A MAIN PROGRAM
C                                            (JOINT WITH  BIRK_TOT  FOR THE ANTISYMMETRICAL MODE)
C
      COMMON /DPHI_B_RHO0/ DPHI,B,RHO_0,XKAPPA ! PARAMETERS, CONTROL DAY-NIGHT ASYMMETRY OF F.A.C.
c
      DATA SH11/ 4956703.683,-26922641.21,-11383659.85,29604361.65,
     *-38919785.97,70230899.72,34993479.24,-90409215.02,30448713.69,
     *-48360257.19,-35556751.23,57136283.60,-8013815.613,30784907.86,
     * 13501620.50,-35121638.52,50297295.45,-84200377.18,-46946852.58,
     * 107526898.8,-39003263.47,59465850.17,47264335.10,-68892388.73,
     * 3375901.533,-9181255.754,-4494667.217,10812618.51,-17351920.97,
     * 27016083.00,18150032.11,-33186882.96,13340198.63,-19779685.30,
     * -17891788.15,21625767.23,16135.32442,133094.0241,-13845.61859,
     *-79159.98442,432.1215298,-85438.10368,1735.386707,41891.71284,
     * 18158.14923,-105465.8135,-11685.73823,62297.34252,-10811.08476,
     * -87631.38186,9217.499261,52079.94529,-68.29127454,56023.02269,
     * -1246.029857,-27436.42793,-11972.61726,69607.08725,7702.743803,
     * -41114.36810,12.08269108,-21.30967022,-9.100782462,18.26855933,
     * -7.000685929,26.22390883,6.392164144,-21.99351743,2.294204157,
     * -16.10023369,-1.344314750,9.342121230,148.5493329,99.79912328,
     * 70.78093196,35.23177574,47.45346891,58.44877918,139.8135237,
     * 91.96485261,6.983488815,9.055554871,19.80484284,2.860045019,
     * .8213262337E-01,-.7962186676E-05/
c
      DATA SH12/-1210748.720,-52324903.95,-14158413.33,19426123.60,
     * 6808641.947,-5138390.983,-1118600.499,-4675055.459,2059671.506,
     * -1373488.052,-114704.4353,-1435920.472,1438451.655,61199067.17,
     *  16549301.39,-22802423.47,-7814550.995,5986478.728,1299443.190,
     *  5352371.724,-2994351.520,1898553.337,203158.3658,2270182.134,
     * -618083.3112,-25950806.16,-7013783.326,9698792.575,3253693.134,
     * -2528478.464,-546323.4095,-2217735.237,1495336.589,-914647.4222,
     * -114374.1054,-1200441.634,-507068.4700,1163189.975,998411.8381,
     * -861919.3631,5252210.872,-11668550.16,-4113899.385,6972900.950,
     * -2546104.076,7704014.310,2273077.192,-5134603.198,256205.7901,
     * -589970.8086,-503821.0170,437612.8956,-2648640.128,5887640.735,
     *  2074286.234,-3519291.144,1283847.104,-3885817.147,-1145936.942,
     *  2589753.651,-408.7788403,1234.054185,739.8541716,-965.8068853,
     *  3691.383679,-8628.635819,-2855.844091,5268.500178,-1774.372703,
     *  5515.010707,1556.089289,-3665.434660,204.8672197,110.7748799,
     *  87.36036207,5.522491330,31.06364270,73.57632579,281.5331360,
     *  140.3461448,17.07537768,6.729732641,4.100970449,2.780422877,
     *  .8742978101E-01,-.1028562327E-04/
c
      DATA SH21/-67763516.61,-49565522.84,10123356.08,51805446.10,
     * -51607711.68,164360662.1,-4662006.024,-191297217.6,-7204547.103,
     * 30372354.93,-750371.9365,-36564457.17,61114395.65,45702536.50,
     * -9228894.939,-47893708.68,47290934.33,-149155112.0,4226520.638,
     * 173588334.5,7998505.443,-33150962.72,832493.2094,39892545.84,
     * -11303915.16,-8901327.398,1751557.110,9382865.820,-9054707.868,
     * 27918664.50,-788741.7146,-32481294.42,-2264443.753,9022346.503,
     * -233526.0185,-10856269.53,-244450.8850,1908295.272,185445.1967,
     * -1074202.863,41827.75224,-241553.7626,-20199.12580,123235.6084,
     * 199501.4614,-1936498.464,-178857.4074,1044724.507,121044.9917,
     * -946479.9247,-91808.28803,532742.7569,-20742.28628,120633.2193,
     * 10018.49534,-61599.11035,-98709.58977,959095.1770,88500.43489,
     * -517471.5287,-81.56122911,816.2472344,55.30711710,-454.5368824,
     * 25.74693810,-202.5007350,-7.369350794,104.9429812,58.14049362,
     * -685.5919355,-51.71345683,374.0125033,247.9296982,159.2471769,
     * 102.3151816,15.81062488,34.99767599,133.0832773,219.6475201,
     * 107.9582783,10.00264684,7.718306072,25.22866153,5.013583103,
     * .8407754233E-01,-.9613356793E-05/
c
      DATA SH22/-43404887.31,8896854.538,-8077731.036,-10247813.65,
     * 6346729.086,-9416801.212,-1921670.268,7805483.928,2299301.127,
     * 4856980.170,-1253936.462,-4695042.690,54305735.91,-11158768.10,
     * 10051771.85,12837129.47,-6380785.836,12387093.50,1687850.192,
     * -10492039.47,-5777044.862,-6916507.424,2855974.911,7027302.490,
     * -26176628.93,5387959.610,-4827069.106,-6193036.589,2511954.143,
     * -6205105.083,-553187.2984,5341386.847,3823736.361,3669209.068,
     * -1841641.700,-3842906.796,281561.7220,-5013124.630,379824.5943,
     * 2436137.901,-76337.55394,548518.2676,42134.28632,-281711.3841,
     * -365514.8666,-2583093.138,-232355.8377,1104026.712,-131536.3445,
     *  2320169.882,-174967.6603,-1127251.881,35539.82827,-256132.9284,
     * -19620.06116,131598.7965,169033.6708,1194443.500,107320.3699,
     * -510672.0036,1211.177843,-17278.19863,1140.037733,8347.612951,
     * -303.8408243,2405.771304,174.0634046,-1248.722950,-1231.229565,
     * -8666.932647,-754.0488385,3736.878824,227.2102611,115.9154291,
     * 94.34364830,3.625357304,64.03192907,109.0743468,241.4844439,
     * 107.7583478,22.36222385,6.282634037,27.79399216,2.270602235,
     * .8708605901E-01,-.1256706895E-04/
c
      XKAPPA=XKAPPA1        !  FORWARDED IN BIR1N2SY
      X_SC=XKAPPA1-1.1D0    !  FORWARDED IN BIRSH_SY

      CALL BIR1N2SY (1,1,PS,X,Y,Z,FX11,FY11,FZ11)           !  REGION 1, MODE 1
      CALL BIRSH_SY (SH11,PS,X_SC,X,Y,Z,HX11,HY11,HZ11)

      BX11=FX11+HX11
      BY11=FY11+HY11
      BZ11=FZ11+HZ11

      CALL BIR1N2SY (1,2,PS,X,Y,Z,FX12,FY12,FZ12)           !  REGION 1, MODE 2
      CALL BIRSH_SY (SH12,PS,X_SC,X,Y,Z,HX12,HY12,HZ12)

      BX12=FX12+HX12
      BY12=FY12+HY12
      BZ12=FZ12+HZ12

      XKAPPA=XKAPPA2        !  FORWARDED IN BIR1N2SY
      X_SC=XKAPPA2-1.0D0    !  FORWARDED IN BIRSH_SY

      CALL BIR1N2SY (2,1,PS,X,Y,Z,FX21,FY21,FZ21)           !  REGION 2, MODE 1
      CALL BIRSH_SY (SH21,PS,X_SC,X,Y,Z,HX21,HY21,HZ21)

      BX21=FX21+HX21
      BY21=FY21+HY21
      BZ21=FZ21+HZ21

      CALL BIR1N2SY (2,2,PS,X,Y,Z,FX22,FY22,FZ22)           !  REGION 2, MODE 2
      CALL BIRSH_SY (SH22,PS,X_SC,X,Y,Z,HX22,HY22,HZ22)

      BX22=FX22+HX22
      BY22=FY22+HY22
      BZ22=FZ22+HZ22

      RETURN
      END
C
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      SUBROUTINE BIR1N2SY (NUMB,MODE,PS,X,Y,Z,BX,BY,BZ)        !   SEE NB# 6, P.60 and NB#7, P.35-...
C
C   THIS CODE IS VERY SIMILAR TO BIRK_1N2, BUT IT IS FOR THE "SYMMETRICAL" MODE, IN WHICH J_parallel
C     IS A SYMMETRIC (EVEN) FUNCTION OF Ygsm
C
C  CALCULATES COMPONENTS  OF REGION 1/2 FIELD IN SPHERICAL COORDS.  DERIVED FROM THE S/R DIPDEF2C (WHICH
C    DOES THE SAME JOB, BUT INPUT/OUTPUT THERE WAS IN SPHERICAL COORDS, WHILE HERE WE USE CARTESIAN ONES)
C
C   INPUT:  NUMB=1 (2) FOR REGION 1 (2) CURRENTS
C           MODE=1 YIELDS SIMPLE SINUSOIDAL MLT VARIATION, WITH MAXIMUM CURRENT AT DAWN/DUSK MERIDIAN
C     WHILE MODE=2 YIELDS THE SECOND HARMONIC.
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A11(31),A12(31),A21(31),A22(31)
      COMMON /MODENUM/ M
      COMMON /DTHETA/ DTHETA

      COMMON /DPHI_B_RHO0/ DPHI,B,RHO_0,XKAPPA ! THESE PARAMETERS CONTROL DAY-NIGHT ASYMMETRY OF F.A.C., AS FOLLOWS:

C  (1) DPHI:   HALF-DIFFERENCE (IN RADIANS) BETWEEN DAY AND NIGHT LATITUDE OF FAC OVAL AT IONOSPHERIC ALTITUDE;
C              TYPICAL VALUE: 0.06
C  (2) B:      AN ASYMMETRY FACTOR AT HIGH-ALTITUDES;  FOR B=0, THE ONLY ASYMMETRY IS THAT FROM DPHI
C              TYPICAL VALUES: 0.35-0.70
C  (3) RHO_0:  A FIXED PARAMETER, DEFINING THE DISTANCE RHO, AT WHICH THE LATITUDE SHIFT GRADUALLY SATURATES AND
C              STOPS INCREASING
C              ITS VALUE WAS ASSUMED FIXED, EQUAL TO 7.0.
C  (4) XKAPPA: AN OVERALL SCALING FACTOR, WHICH CAN BE USED FOR CHANGING THE SIZE OF THE F.A.C. OVAL
C
      DATA BETA,RH,EPS/0.9D0,10.D0,3.D0/ ! parameters of the tilt-dependent deformation of the untilted F.A.C. field

      DATA A11/.1618068350,-.1797957553,2.999642482,-.9322708978,
     *-.6811059760,.2099057262,-8.358815746,-14.86033550,.3838362986,
     *-16.30945494,4.537022847,2.685836007,27.97833029,6.330871059,
     *1.876532361,18.95619213,.9651528100,.4217195118,-.08957770020,
     *-1.823555887,.7457045438,-.5785916524,-1.010200918,.01112389357,
     *.09572927448,-.3599292276,8.713700514,.9763932955,3.834602998,
     *2.492118385,.7113544659/
      DATA A12/.7058026940,-.2845938535,5.715471266,-2.472820880,
     *-.7738802408,.3478293930,-11.37653694,-38.64768867,.6932927651,
     *-212.4017288,4.944204937,3.071270411,33.05882281,7.387533799,
     *2.366769108,79.22572682,.6154290178,.5592050551,-.1796585105,
     *-1.654932210,.7309108776,-.4926292779,-1.130266095,-.009613974555,
     *.1484586169,-.2215347198,7.883592948,.02768251655,2.950280953,
     *1.212634762,.5567714182/
      DATA A21/.1278764024,-.2320034273,1.805623266,-32.37241440,
     *-.9931490648,.3175085630,-2.492465814,-16.21600096,.2695393416,
     *-6.752691265,3.971794901,14.54477563,41.10158386,7.912889730,
     *1.258297372,9.583547721,1.014141963,.5104134759,-.1790430468,
     *-1.756358428,.7561986717,-.6775248254,-.04014016420,.01446794851,
     *.1200521731,-.2203584559,4.508963850,.8221623576,1.779933730,
     *1.102649543,.8867880020/
      DATA A22/.4036015198,-.3302974212,2.827730930,-45.44405830,
     *-1.611103927,.4927112073,-.003258457559,-49.59014949,.3796217108,
     *-233.7884098,4.312666980,18.05051709,28.95320323,11.09948019,
     *.7471649558,67.10246193,.5667096597,.6468519751,-.1560665317,
     *-1.460805289,.7719653528,-.6658988668,.2515179349E-05,
     *.02426021891,.1195003324,-.2625739255,4.377172556,.2421190547,
     *2.503482679,1.071587299,.7247997430/

      B=0.5
      RHO_0=7.0

      M=MODE
      IF (NUMB.EQ.1) THEN
          DPHI=0.055D0
          DTHETA=0.06D0
      ENDIF

      IF (NUMB.EQ.2) THEN
          DPHI=0.030D0
          DTHETA=0.09D0
      ENDIF

      Xsc=X*XKAPPA
      Ysc=Y*XKAPPA
      Zsc=Z*XKAPPA
      RHO=DSQRT(Xsc**2+Zsc**2)

      Rsc=DSQRT(Xsc**2+Ysc**2+Zsc**2)                                 !  SCALED
      RHO2=RHO_0**2

      IF (Xsc.EQ.0.D0.AND.Zsc.EQ.0.D0) THEN
         PHI=0.D0
      ELSE
         PHI=DATAN2(-Zsc,Xsc)  !  FROM CARTESIAN TO CYLINDRICAL (RHO,PHI,Y)
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
      DPHISDY= BETA*PS*R1RH**(EPS-1.D0)*Ysc/(RH*Rsc*
     *   (1.D0+R1RH**EPS)**(1.D0/EPS+1.D0))

      SPHICS=DSIN(PHIS)
      CPHICS=DCOS(PHIS)

      XS= RHO*CPHICS
      ZS=-RHO*SPHICS

      IF (NUMB.EQ.1) THEN
        IF (MODE.EQ.1) CALL TWOCONSS (A11,XS,Ysc,ZS,BXS,BYAS,BZS)
        IF (MODE.EQ.2) CALL TWOCONSS (A12,XS,Ysc,ZS,BXS,BYAS,BZS)
      ELSE
        IF (MODE.EQ.1) CALL TWOCONSS (A21,XS,Ysc,ZS,BXS,BYAS,BZS)
        IF (MODE.EQ.2) CALL TWOCONSS (A22,XS,Ysc,ZS,BXS,BYAS,BZS)
      ENDIF

      BRHOAS=BXS*CPHICS-BZS*SPHICS
      BPHIAS=-BXS*SPHICS-BZS*CPHICS

      BRHO_S=BRHOAS*DPHISPHI                             *XKAPPA        ! SCALING
      BPHI_S=(BPHIAS-RHO*(BYAS*DPHISDY+BRHOAS*DPHISRHO)) *XKAPPA
      BY_S=BYAS*DPHISPHI                                 *XKAPPA

      BX=BRHO_S*CPHIC-BPHI_S*SPHIC
      BY=BY_S
      BZ=-BRHO_S*SPHIC-BPHI_S*CPHIC

      RETURN
      END
c

C=========================================================================
c
      SUBROUTINE TWOCONSS (A,X,Y,Z,BX,BY,BZ)
C
C   DIFFERS FROM TWOCONES:  THIS S/R IS FOR THE "SYMMETRIC" MODE OF BIRKELAND CURRENTS IN THAT
C                           HERE THE FIELD IS ROTATED BY 90 DEGS FOR M=1 AND BY 45 DEGS FOR M=2
C
C   ADDS FIELDS FROM TWO CONES (NORTHERN AND SOUTHERN), WITH A PROPER SYMMETRY OF THE CURRENT AND FIELD,
C     CORRESPONDING TO THE REGION 1 BIRKELAND CURRENTS. (SEE NB #6, P.58).
C

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(31)
      COMMON /MODENUM/ M
      DATA HSQR2/0.707106781D0/

      IF (M.EQ.1) THEN   !   ROTATION BY 90 DEGS
       XAS = Y
       YAS =-X
      ELSE               !   ROTATION BY 45 DEGS
       XAS = (X+Y)*HSQR2
       YAS = (Y-X)*HSQR2
      ENDIF

      CALL ONE_CONE (A,XAS,YAS,Z,BXN,BYN,BZN)
      CALL ONE_CONE (A,XAS,-YAS,-Z,BXS,BYS,BZS)

      BXAS=BXN-BXS
      BYAS=BYN+BYS
      BZ=BZN+BZS

      IF (M.EQ.1) THEN   !   ROTATION BY 90 DEGS
        BX =-BYAS
        BY = BXAS
      ELSE
        BX=(BXAS-BYAS)*HSQR2
        BY=(BXAS+BYAS)*HSQR2
      ENDIF

      RETURN
      END
c
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
         SUBROUTINE BIRSH_SY (A,PS,X_SC,X,Y,Z,BX,BY,BZ)
C
C   THIS S/R IS QUITE SIMILAR TO BIRK_SHL, BUT IT IS FOR THE SYMMETRIC MODE OF BIRKELAND CURRENT FIELD
C     AND FOR THAT REASON THE FIELD COMPONENTS HAVE A DIFFERENT KIND OF SYMMETRY WITH RESPECT TO Y_gsm
C
         IMPLICIT  REAL * 8  (A - H, O - Z)
         DIMENSION A(86)
C
         CPS=DCOS(PS)
         SPS=DSIN(PS)

         S3PS=2.D0*CPS
C
         PST1=PS*A(85)
         PST2=PS*A(86)

         ST1=DSIN(PST1)
         CT1=DCOS(PST1)
         ST2=DSIN(PST2)
         CT2=DCOS(PST2)

         X1=X*CT1-Z*ST1
         Z1=X*ST1+Z*CT1
         X2=X*CT2-Z*ST2
         Z2=X*ST2+Z*CT2
C
         L=0
         GX=0.D0
         GY=0.D0
         GZ=0.D0
C
         DO 1 M=1,2     !    M=1 IS FOR THE 1ST SUM ("PERP." SYMMETRY)
C                          AND M=2 IS FOR THE SECOND SUM ("PARALL." SYMMETRY)
             DO 2 I=1,3
                  P=A(72+I)
                  Q=A(78+I)
                  CYPI=DCOS(Y/P)
                  CYQI=DCOS(Y/Q)
                  SYPI=DSIN(Y/P)
                  SYQI=DSIN(Y/Q)
C
                DO 3 K=1,3
                   R=A(75+K)
                   S=A(81+K)
                   SZRK=DSIN(Z1/R)
                   CZSK=DCOS(Z2/S)
                   CZRK=DCOS(Z1/R)
                   SZSK=DSIN(Z2/S)
                     SQPR=DSQRT(1.D0/P**2+1.D0/R**2)
                     SQQS=DSQRT(1.D0/Q**2+1.D0/S**2)
                        EPR=DEXP(X1*SQPR)
                        EQS=DEXP(X2*SQQS)
C
                  DO 4 N=1,2  ! N=1 IS FOR THE FIRST PART OF EACH COEFFICIENT
C                                AND N=2 IS FOR THE SECOND ONE

                    DO 5 NN=1,2 !   NN = 1,2 FURTHER SPLITS THE COEFFICIENTS INTO 2 PARTS,
C                                         TO TAKE INTO ACCOUNT THE SCALE FACTOR DEPENDENCE

                    IF (M.EQ.1) THEN
                         FX= SQPR*EPR*SYPI*SZRK
                         FY=EPR*CYPI*SZRK/P
                         FZ= EPR*SYPI*CZRK/R
                       IF (N.EQ.1) THEN
                         IF (NN.EQ.1) THEN
                          HX=FX
                          HY=FY
                          HZ=FZ
                         ELSE
                          HX=FX*X_SC
                          HY=FY*X_SC
                          HZ=FZ*X_SC
                         ENDIF
                       ELSE
                         IF (NN.EQ.1) THEN
                          HX=FX*CPS
                          HY=FY*CPS
                          HZ=FZ*CPS
                         ELSE
                          HX=FX*CPS*X_SC
                          HY=FY*CPS*X_SC
                          HZ=FZ*CPS*X_SC
                         ENDIF
                       ENDIF

                     ELSE                            !   M.EQ.2
                         FX= SPS*SQQS*EQS*SYQI*CZSK
                         FY= SPS/Q*EQS*CYQI*CZSK
                         FZ=-SPS/S*EQS*SYQI*SZSK
                       IF (N.EQ.1) THEN
                        IF (NN.EQ.1) THEN
                          HX=FX
                          HY=FY
                          HZ=FZ
                        ELSE
                          HX=FX*X_SC
                          HY=FY*X_SC
                          HZ=FZ*X_SC
                        ENDIF
                       ELSE
                        IF (NN.EQ.1) THEN
                         HX=FX*S3PS
                         HY=FY*S3PS
                         HZ=FZ*S3PS
                        ELSE
                         HX=FX*S3PS*X_SC
                         HY=FY*S3PS*X_SC
                         HZ=FZ*S3PS*X_SC
                        ENDIF
                       ENDIF
                  ENDIF
       L=L+1

       IF (M.EQ.1) THEN
       HXR=HX*CT1+HZ*ST1
       HZR=-HX*ST1+HZ*CT1
       ELSE
       HXR=HX*CT2+HZ*ST2
       HZR=-HX*ST2+HZ*CT2
       ENDIF

       GX=GX+HXR*A(L)
       GY=GY+HY *A(L)
  5    GZ=GZ+HZR*A(L)

  4   CONTINUE
  3   CONTINUE
  2   CONTINUE
  1   CONTINUE

      BX=GX
      BY=GY
      BZ=GZ

      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
       SUBROUTINE DIPOLE (PS,X,Y,Z,BX,BY,BZ)
C
C      A DOUBLE PRECISION ROUTINE
C
C  CALCULATES GSM COMPONENTS OF A GEODIPOLE FIELD WITH THE DIPOLE MOMENT
C  CORRESPONDING TO THE EPOCH OF 2000.
C
C----INPUT PARAMETERS:
C     PS - GEODIPOLE TILT ANGLE IN RADIANS,
C     X,Y,Z - GSM COORDINATES IN RE (1 RE = 6371.2 km)
C
C----OUTPUT PARAMETERS:
C     BX,BY,BZ - FIELD COMPONENTS IN GSM SYSTEM, IN NANOTESLA.
C
C  LAST MODIFICATION: JAN. 5, 2001. THE VALUE OF THE DIPOLE MOMENT WAS UPDATED TO 2000.
C    AND A "SAVE" STATEMENT HAS BEEN ADDED, TO AVOID POTENTIAL PROBLEMS WITH SOME
C    FORTRAN COMPILERS
C
C  WRITTEN BY: N. A. TSYGANENKO
C
      IMPLICIT REAL*8 (A-H,O-Z)
      SAVE M,PSI,SPS,CPS
      DATA M,PSI/0,5.D0/
      IF(M.EQ.1.AND.DABS(PS-PSI).LT.1.D-5) GOTO 1
      SPS=DSIN(PS)
      CPS=DCOS(PS)
      PSI=PS
      M=1
  1   P=X**2
      U=Z**2
      V=3.D0*Z*X
      T=Y**2
      Q=30115.D0/DSQRT(P+T+U)**5
      BX=Q*((T+U-2.D0*P)*SPS-V*CPS)
      BY=-3.D0*Y*Q*(X*SPS+Z*CPS)
      BZ=Q*((P+T-2.D0*U)*CPS-V*SPS)
      RETURN
      END


C(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
C
C(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
C=======================================================================================
C
      SUBROUTINE T96_MGNP_D (XN_PD,VEL,XGSM,YGSM,ZGSM,XMGNP,YMGNP,ZMGNP,
     * DIST,ID)

      IMPLICIT REAL*8 (A-H,O-Z)
C
C  DOUBLE-PRECISION VERSION !!!!!!!!   HENCE THE SUFFIX "D" IN THE NAME
C
C
C  FOR ANY POINT OF SPACE WITH GIVEN COORDINATES (XGSM,YGSM,ZGSM), THIS SUBROUTINE DEFINES
C  THE POSITION OF A POINT (XMGNP,YMGNP,ZMGNP) AT THE T96 MODEL MAGNETOPAUSE, HAVING THE
C  SAME VALUE OF THE ELLIPSOIDAL TAU-COORDINATE, AND THE DISTANCE BETWEEN THEM.  THIS IS
C  NOT THE SHORTEST DISTANCE D_MIN TO THE BOUNDARY, BUT DIST ASYMPTOTICALLY TENDS TO D_MIN,
C  AS THE OBSERVATION POINT GETS CLOSER TO THE MAGNETOPAUSE.
C
C  INPUT: XN_PD - EITHER SOLAR WIND PROTON NUMBER DENSITY (PER C.C.) (IF VEL>0)
C                    OR THE SOLAR WIND RAM PRESSURE IN NANOPASCALS   (IF VEL<0)
C         VEL - EITHER SOLAR WIND VELOCITY (KM/SEC)
C                  OR ANY NEGATIVE NUMBER, WHICH INDICATES THAT XN_PD STANDS
C                     FOR THE SOLAR WIND PRESSURE, RATHER THAN FOR THE DENSITY
C
C         XGSM,YGSM,ZGSM - COORDINATES OF THE OBSERVATION POINT IN EARTH RADII
C
C  OUTPUT: XMGNP,YMGNP,ZMGNP - GSM POSITION OF THE BOUNDARY POINT, HAVING THE SAME
C          VALUE OF TAU-COORDINATE AS THE OBSERVATION POINT (XGSM,YGSM,ZGSM)
C          DIST -  THE DISTANCE BETWEEN THE TWO POINTS, IN RE,
C          ID -    POSITION FLAG; ID=+1 (-1) MEANS THAT THE POINT (XGSM,YGSM,ZGSM)
C          LIES INSIDE (OUTSIDE) THE MODEL MAGNETOPAUSE, RESPECTIVELY.
C
C  THE PRESSURE-DEPENDENT MAGNETOPAUSE IS THAT USED IN THE T96_01 MODEL
C  (TSYGANENKO, JGR, V.100, P.5599, 1995; ESA SP-389, P.181, OCT. 1996)
C
c   AUTHOR:  N.A. TSYGANENKO
C   DATE:    AUG.1, 1995, REVISED APRIL 3, 2003.
C
C
C  DEFINE SOLAR WIND DYNAMIC PRESSURE (NANOPASCALS, ASSUMING 4% OF ALPHA-PARTICLES),
C   IF NOT EXPLICITLY SPECIFIED IN THE INPUT:

      IF (VEL.LT.0.D0) THEN
       PD=XN_PD
      ELSE
       PD=1.94D-6*XN_PD*VEL**2
C
      ENDIF
C
C  RATIO OF PD TO THE AVERAGE PRESSURE, ASSUMED EQUAL TO 2 nPa:

      RAT=PD/2.0D0
      RAT16=RAT**0.14D0

C (THE POWER INDEX 0.14 IN THE SCALING FACTOR IS THE BEST-FIT VALUE OBTAINED FROM DATA
C    AND USED IN THE T96_01 VERSION)
C
C  VALUES OF THE MAGNETOPAUSE PARAMETERS FOR  PD = 2 nPa:
C
      A0 =34.586D0
      S00=1.196D0
      X00=3.4397D0
C
C   VALUES OF THE MAGNETOPAUSE PARAMETERS, SCALED BY THE ACTUAL PRESSURE:
C
      A=A0/RAT16
      S0=S00
      X0=X00/RAT16
      XM=X0-A
C
C  (XM IS THE X-COORDINATE OF THE "SEAM" BETWEEN THE ELLIPSOID AND THE CYLINDER)
C
C     (FOR DETAILS ON THE ELLIPSOIDAL COORDINATES, SEE THE PAPER:
C      N.A.TSYGANENKO, SOLUTION OF CHAPMAN-FERRARO PROBLEM FOR AN
C      ELLIPSOIDAL MAGNETOPAUSE, PLANET.SPACE SCI., V.37, P.1037, 1989).
C
       IF (YGSM.NE.0.D0.OR.ZGSM.NE.0.D0) THEN
          PHI=DATAN2(YGSM,ZGSM)
       ELSE
          PHI=0.D0
       ENDIF
C
       RHO=DSQRT(YGSM**2+ZGSM**2)
C
       IF (XGSM.LT.XM) THEN
           XMGNP=XGSM
           RHOMGNP=A*DSQRT(S0**2-1.D0)
           YMGNP=RHOMGNP*DSIN(PHI)
           ZMGNP=RHOMGNP*DCOS(PHI)
           DIST=DSQRT((XGSM-XMGNP)**2+(YGSM-YMGNP)**2+(ZGSM-ZMGNP)**2)
           IF (RHOMGNP.GT.RHO) ID=+1
           IF (RHOMGNP.LE.RHO) ID=-1
           RETURN
       ENDIF
C
          XKSI=(XGSM-X0)/A+1.D0
          XDZT=RHO/A
          SQ1=DSQRT((1.D0+XKSI)**2+XDZT**2)
          SQ2=DSQRT((1.D0-XKSI)**2+XDZT**2)
          SIGMA=0.5D0*(SQ1+SQ2)
          TAU=0.5D0*(SQ1-SQ2)
C
C  NOW CALCULATE (X,Y,Z) FOR THE CLOSEST POINT AT THE MAGNETOPAUSE
C
          XMGNP=X0-A*(1.D0-S0*TAU)
          ARG=(S0**2-1.D0)*(1.D0-TAU**2)
          IF (ARG.LT.0.D0) ARG=0.D0
          RHOMGNP=A*DSQRT(ARG)
          YMGNP=RHOMGNP*DSIN(PHI)
          ZMGNP=RHOMGNP*DCOS(PHI)
C
C  NOW CALCULATE THE DISTANCE BETWEEN THE POINTS {XGSM,YGSM,ZGSM} AND {XMGNP,YMGNP,ZMGNP}:
C   (IN GENERAL, THIS IS NOT THE SHORTEST DISTANCE D_MIN, BUT DIST ASYMPTOTICALLY TENDS
C    TO D_MIN, AS WE ARE GETTING CLOSER TO THE MAGNETOPAUSE):
C
      DIST=DSQRT((XGSM-XMGNP)**2+(YGSM-YMGNP)**2+(ZGSM-ZMGNP)**2)
C
      IF (SIGMA.GT.S0) ID=-1   !  ID=-1 MEANS THAT THE POINT LIES OUTSIDE
      IF (SIGMA.LE.S0) ID=+1   !  ID=+1 MEANS THAT THE POINT LIES INSIDE
C                                           THE MAGNETOSPHERE
      RETURN
      END

c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ c
c
c          ##########################################################################
c          #                                                                        #
c          #                          GEOPACK-2008_dp                               #
c          #                     (MAIN SET OF FORTRAN CODES)                        #
c          #               (double-precision version of 03/21/08)                   #
C          #                (IGRF coefficients updated 01/31/15)                    #
c          ##########################################################################
C
c
c  This collection of subroutines is a result of several upgrades of the original package
c  written by N. A. Tsyganenko in 1978-1979.
c
c  PREFATORY NOTE TO THE VERSION OF FEBRUARY 4, 2008:
c
c  To avoid inappropriate use of obsolete subroutines from earlier versions, a suffix 08 was
c  added to the name of each subroutine in this release.
c
c  A possibility has been added in this version to calculate vector components in the
c  "Geocentric Solar Wind" (GSW) coordinate system, which, to our knowledge, was first
c  introduced by Hones et al., Planet. Space Sci., v.34, p.889, 1986 (aka GSWM, see Appendix,
c  Tsyganenko et al., JGRA, v.103(A4), p.6827, 1998). The GSW system is analogous to the
c  standard GSM, except that its X-axis is antiparallel to the currently observed solar wind
c  flow vector, rather than aligned with the Earth-Sun line. The orientation of axes in the
c  GSW system can be uniquely defined by specifying three components (VGSEX,VGSEY,VGSEZ) of
c  the solar wind velocity, and in the case of a strictly radial anti-sunward flow (VGSEY=
c  VGSEZ=0) the GSW system becomes identical to the standard GSM, which fact was used here
c  to minimize the number of subroutines in the package. To that end, instead of the special
c  case of the GSM coordinates, this version uses a more general GSW system, and three more
c  input parameters are added in the subroutine RECALC_08, the observed values (VGSEX,VGSEY,
c  VGSEZ) of the solar wind velocity. Invoking RECALC_08 with VGSEY=VGSEZ=0 restores the
c  standard (sunward) orientation of the X axis, which allows one to easily convert vectors
c  between GSW and GSM, as well as to/from other standard and commonly used systems. For more
c  details, see the documentation file GEOPACK-2008.DOC.
c
c  Another modification allows users to have more control over the procedure of field line
c  mapping using the subroutine TRACE_08. To that end, three new input parameters were added
c  in that subroutine, allowing one to set (i) an upper limit, DSMAX, on the automatically
c  adjusted step size, (ii) a permissible step error, ERR, and (iii) maximal length, LMAX,
c  of arrays where field line point coordinates are stored. Minor changes were also made in
c  the tracing subroutine, to make it more compact and easier for understanding, and to
c  prevent the algorithm from making uncontrollable large number of multiple loops in some
c  cases with plasmoid-like field structures.
c
C  One more subroutine, named GEODGEO_08, was added to the package, allowing one to convert
c  geodetic coordinates of a point in space (altitude above the Earth's WGS84 ellipsoid and
c  geodetic latitude) to geocentric radial distance and colatitude, and vice versa.
c
C  For a complete list of modifications made earlier in previous versions, see the
c  documentation file GEOPACK-2008.DOC.
c
c----------------------------------------------------------------------------------
c
      SUBROUTINE IGRF_GSW_08 (XGSW,YGSW,ZGSW,HXGSW,HYGSW,HZGSW)
c
C  CALCULATES COMPONENTS OF THE MAIN (INTERNAL) GEOMAGNETIC FIELD IN THE GEOCENTRIC SOLAR-WIND
C  (GSW) COORDINATE SYSTEM, USING IAGA INTERNATIONAL GEOMAGNETIC REFERENCE MODEL COEFFICIENTS
C  (e.g., http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html, revised 22 March, 2005)
c
C  THE GSW SYSTEM IS ESSENTIALLY SIMILAR TO THE STANDARD GSM (THE TWO SYSTEMS BECOME IDENTICAL
C  TO EACH OTHER IN THE CASE OF STRICTLY ANTI-SUNWARD SOLAR WIND FLOW). FOR A DETAILED
C  DEFINITION, SEE INTRODUCTORY COMMENTS FOR THE SUBROUTINE GSWGSE_08 .
C
C  BEFORE THE FIRST CALL OF THIS SUBROUTINE, OR, IF THE DATE/TIME (IYEAR,IDAY,IHOUR,MIN,ISEC),
C  OR THE SOLAR WIND VELOCITY COMPONENTS (VGSEX,VGSEY,VGSEZ) HAVE CHANGED, THE MODEL COEFFICIENTS
c  AND GEO-GSW ROTATION MATRIX ELEMENTS SHOULD BE UPDATED BY CALLING THE SUBROUTINE RECALC_08.
C
C-----INPUT PARAMETERS:
C
C     XGSW,YGSW,ZGSW - CARTESIAN GEOCENTRIC SOLAR-WIND COORDINATES (IN UNITS RE=6371.2 KM)
C
C-----OUTPUT PARAMETERS:
C
C     HXGSW,HYGSW,HZGSW - CARTESIAN GEOCENTRIC SOLAR-WIND COMPONENTS OF THE MAIN GEOMAGNETIC
C                           FIELD IN NANOTESLA
C
C     LAST MODIFICATION:  MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C     THIS VERSION OF THE CODE ACCEPTS DATES FROM 1965 THROUGH 2015.
c
C     AUTHOR: N. A. TSYGANENKO
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /GEOPACK2/ G(105),H(105),REC(105)

      DIMENSION A(14),B(14)

      CALL GEOGSW_08 (XGEO,YGEO,ZGEO,XGSW,YGSW,ZGSW,-1)
      RHO2=XGEO**2+YGEO**2
      R=DSQRT(RHO2+ZGEO**2)
      C=ZGEO/R
      RHO=DSQRT(RHO2)
      S=RHO/R
      IF (S.LT.1.D-10) THEN
        CF=1.D0
        SF=0.D0
      ELSE
        CF=XGEO/RHO
        SF=YGEO/RHO
      ENDIF
C
      PP=1.D0/R
      P=PP
C
C  IN THIS VERSION, THE OPTIMAL VALUE OF THE PARAMETER NM (MAXIMAL ORDER OF THE SPHERICAL
C    HARMONIC EXPANSION) IS NOT USER-PRESCRIBED, BUT CALCULATED INSIDE THE SUBROUTINE, BASED
C      ON THE VALUE OF THE RADIAL DISTANCE R:
C
      IRP3=R+2
      NM=3+30/IRP3
      IF (NM.GT.13) NM=13

      K=NM+1
      DO 150 N=1,K
         P=P*PP
         A(N)=P
150      B(N)=P*N

      P=1.D0
      D=0.D0
      BBR=0.D0
      BBT=0.D0
      BBF=0.D0

      DO 200 M=1,K
         IF(M.EQ.1) GOTO 160
         MM=M-1
         W=X
         X=W*CF+Y*SF
         Y=Y*CF-W*SF
         GOTO 170
160      X=0.D0
         Y=1.D0
170      Q=P
         Z=D
         BI=0.D0
         P2=0.D0
         D2=0.D0
         DO 190 N=M,K
            AN=A(N)
            MN=N*(N-1)/2+M
            E=G(MN)
            HH=H(MN)
            W=E*Y+HH*X
            BBR=BBR+B(N)*W*Q
            BBT=BBT-AN*W*Z
            IF(M.EQ.1) GOTO 180
            QQ=Q
            IF(S.LT.1.D-10) QQ=Z
            BI=BI+AN*(E*X-HH*Y)*QQ
180         XK=REC(MN)
            DP=C*Z-S*Q-XK*D2
            PM=C*Q-XK*P2
            D2=Z
            P2=Q
            Z=DP
190        Q=PM
         D=S*D+C*P
         P=S*P
         IF(M.EQ.1) GOTO 200
         BI=BI*MM
         BBF=BBF+BI
200   CONTINUE
C
      BR=BBR
      BT=BBT
      IF(S.LT.1.D-10) GOTO 210
      BF=BBF/S
      GOTO 211
210   IF(C.LT.0.) BBF=-BBF
      BF=BBF

211   HE=BR*S+BT*C
      HXGEO=HE*CF-BF*SF
      HYGEO=HE*SF+BF*CF
      HZGEO=BR*C-BT*S
C
      CALL GEOGSW_08 (HXGEO,HYGEO,HZGEO,HXGSW,HYGSW,HZGSW,1)
C
      RETURN
      END
C
c==========================================================================================
C
c
      SUBROUTINE IGRF_GEO_08 (R,THETA,PHI,BR,BTHETA,BPHI)
c
C  CALCULATES COMPONENTS OF THE MAIN (INTERNAL) GEOMAGNETIC FIELD IN THE SPHERICAL GEOGRAPHIC
C  (GEOCENTRIC) COORDINATE SYSTEM, USING IAGA INTERNATIONAL GEOMAGNETIC REFERENCE MODEL
C  COEFFICIENTS  (e.g., http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html, revised: 22 March, 2005)
C
C  BEFORE THE FIRST CALL OF THIS SUBROUTINE, OR IF THE DATE (IYEAR AND IDAY) WAS CHANGED,
C  THE MODEL COEFFICIENTS SHOULD BE UPDATED BY CALLING THE SUBROUTINE RECALC_08
C
C-----INPUT PARAMETERS:
C
C   R, THETA, PHI - SPHERICAL GEOGRAPHIC (GEOCENTRIC) COORDINATES:
C   RADIAL DISTANCE R IN UNITS RE=6371.2 KM, COLATITUDE THETA AND LONGITUDE PHI IN RADIANS
C
C-----OUTPUT PARAMETERS:
C
C     BR, BTHETA, BPHI - SPHERICAL COMPONENTS OF THE MAIN GEOMAGNETIC FIELD IN NANOTESLA
C      (POSITIVE BR OUTWARD, BTHETA SOUTHWARD, BPHI EASTWARD)
C
C     LAST MODIFICATION:  MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C     THIS VERSION OF THE  CODE ACCEPTS DATES FROM 1965 THROUGH 2015.
c
C     AUTHOR: N. A. TSYGANENKO
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /GEOPACK2/ G(105),H(105),REC(105)

      DIMENSION A(14),B(14)

      C=DCOS(THETA)
      S=DSIN(THETA)
      CF=DCOS(PHI)
      SF=DSIN(PHI)
C
      PP=1.D0/R
      P=PP
C
C  IN THIS NEW VERSION, THE OPTIMAL VALUE OF THE PARAMETER NM (MAXIMAL ORDER OF THE SPHERICAL
C    HARMONIC EXPANSION) IS NOT USER-PRESCRIBED, BUT CALCULATED INSIDE THE SUBROUTINE, BASED
C      ON THE VALUE OF THE RADIAL DISTANCE R:
C
      IRP3=R+2
      NM=3+30/IRP3
      IF (NM.GT.13) NM=13

      K=NM+1
      DO 150 N=1,K
         P=P*PP
         A(N)=P
150      B(N)=P*N

      P=1.D0
      D=0.D0
      BBR=0.D0
      BBT=0.D0
      BBF=0.D0

      DO 200 M=1,K
         IF(M.EQ.1) GOTO 160
         MM=M-1
         W=X
         X=W*CF+Y*SF
         Y=Y*CF-W*SF
         GOTO 170
160      X=0.D0
         Y=1.D0
170      Q=P
         Z=D
         BI=0.D0
         P2=0.D0
         D2=0.D0
         DO 190 N=M,K
            AN=A(N)
            MN=N*(N-1)/2+M
            E=G(MN)
            HH=H(MN)
            W=E*Y+HH*X
            BBR=BBR+B(N)*W*Q
            BBT=BBT-AN*W*Z
            IF(M.EQ.1) GOTO 180
            QQ=Q
            IF(S.LT.1.D-5) QQ=Z
            BI=BI+AN*(E*X-HH*Y)*QQ
180         XK=REC(MN)
            DP=C*Z-S*Q-XK*D2
            PM=C*Q-XK*P2
            D2=Z
            P2=Q
            Z=DP
190        Q=PM
         D=S*D+C*P
         P=S*P
         IF(M.EQ.1) GOTO 200
         BI=BI*MM
         BBF=BBF+BI
200   CONTINUE
C
      BR=BBR
      BTHETA=BBT
      IF(S.LT.1.D-10) GOTO 210
      BPHI=BBF/S
      RETURN
210   IF(C.LT.0.D0) BBF=-BBF
      BPHI=BBF

      RETURN
      END
C
c==========================================================================================
c
       SUBROUTINE DIP_08 (XGSW,YGSW,ZGSW,BXGSW,BYGSW,BZGSW)
C
C  CALCULATES GSW (GEOCENTRIC SOLAR-WIND) COMPONENTS OF GEODIPOLE FIELD WITH THE DIPOLE MOMENT
C  CORRESPONDING TO THE EPOCH, SPECIFIED BY CALLING SUBROUTINE RECALC_08 (SHOULD BE
C  INVOKED BEFORE THE FIRST USE OF THIS ONE, OR IF THE DATE/TIME, AND/OR THE OBSERVED
C  SOLAR WIND DIRECTION, HAVE CHANGED.
C
C  THE GSW COORDINATE SYSTEM IS ESSENTIALLY SIMILAR TO THE STANDARD GSM (THE TWO SYSTEMS BECOME
C  IDENTICAL TO EACH OTHER IN THE CASE OF STRICTLY RADIAL ANTI-SUNWARD SOLAR WIND FLOW). ITS
C  DETAILED DEFINITION IS GIVEN IN INTRODUCTORY COMMENTS FOR THE SUBROUTINE GSWGSE_08 .

C--INPUT PARAMETERS: XGSW,YGSW,ZGSW - GSW COORDINATES IN RE (1 RE = 6371.2 km)
C
C--OUTPUT PARAMETERS: BXGSW,BYGSW,BZGSW - FIELD COMPONENTS IN GSW SYSTEM, IN NANOTESLA.
C
C  LAST MODIFICATION:   MARCH 21, 2008 (DOUBLE-PRECISION VERSION).
C
C  AUTHOR: N. A. TSYGANENKO
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /GEOPACK1/ AA(10),SPS,CPS,BB(22)
      COMMON /GEOPACK2/ G(105),H(105),REC(105)
C
      DIPMOM=DSQRT(G(2)**2+G(3)**2+H(3)**2)
C
      P=XGSW**2
      U=ZGSW**2
      V=3.D0*ZGSW*XGSW
      T=YGSW**2
      Q=DIPMOM/DSQRT(P+T+U)**5
      BXGSW=Q*((T+U-2.D0*P)*SPS-V*CPS)
      BYGSW=-3.D0*YGSW*Q*(XGSW*SPS+ZGSW*CPS)
      BZGSW=Q*((P+T-2.D0*U)*CPS-V*SPS)
      RETURN
      END

C*******************************************************************
c
      SUBROUTINE SUN_08 (IYEAR,IDAY,IHOUR,MIN,ISEC,GST,SLONG,SRASN,SDEC)
C
C  CALCULATES FOUR QUANTITIES NECESSARY FOR COORDINATE TRANSFORMATIONS
C  WHICH DEPEND ON SUN POSITION (AND, HENCE, ON UNIVERSAL TIME AND SEASON)
C
C-------  INPUT PARAMETERS:
C  IYR,IDAY,IHOUR,MIN,ISEC -  YEAR, DAY, AND UNIVERSAL TIME IN HOURS, MINUTES,
C    AND SECONDS  (IDAY=1 CORRESPONDS TO JANUARY 1).
C
C-------  OUTPUT PARAMETERS:
C  GST - GREENWICH MEAN SIDEREAL TIME, SLONG - LONGITUDE ALONG ECLIPTIC
C  SRASN - RIGHT ASCENSION,  SDEC - DECLINATION  OF THE SUN (RADIANS)
C  ORIGINAL VERSION OF THIS SUBROUTINE HAS BEEN COMPILED FROM:
C  RUSSELL, C.T., COSMIC ELECTRODYNAMICS, 1971, V.2, PP.184-196.
C
C  LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C
C     ORIGINAL VERSION WRITTEN BY:    Gilbert D. Mead
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA RAD/57.295779513D0/
C
      IF(IYEAR.LT.1901.OR.IYEAR.GT.2099) RETURN
      FDAY=DFLOAT(IHOUR*3600+MIN*60+ISEC)/86400.D0
      DJ=365*(IYEAR-1900)+(IYEAR-1901)/4+IDAY-0.5D0+FDAY
      T=DJ/36525.D0
      VL=DMOD(279.696678D0+0.9856473354D0*DJ,360.D0)
      GST=DMOD(279.690983D0+.9856473354D0*DJ+360.D0*FDAY+180.D0,360.D0)/
     * RAD
      G=DMOD(358.475845D0+0.985600267D0*DJ,360.D0)/RAD
      SLONG=(VL+(1.91946D0-0.004789D0*T)*DSIN(G)+0.020094D0
     * *DSIN(2.D0*G))/RAD
      IF(SLONG.GT.6.2831853D0) SLONG=SLONG-6.283185307D0
      IF (SLONG.LT.0.D0) SLONG=SLONG+6.283185307D0
      OBLIQ=(23.45229D0-0.0130125D0*T)/RAD
      SOB=DSIN(OBLIQ)
      SLP=SLONG-9.924D-5
C
C   THE LAST CONSTANT IS A CORRECTION FOR THE ANGULAR ABERRATION DUE TO
C   EARTH'S ORBITAL MOTION
C
      SIND=SOB*DSIN(SLP)
      COSD=DSQRT(1.D0-SIND**2)
      SC=SIND/COSD
      SDEC=DATAN(SC)
      SRASN=3.141592654D0-DATAN2(DCOS(OBLIQ)/SOB*SC,-DCOS(SLP)/COSD)
      RETURN
      END
C
C================================================================================
c
      SUBROUTINE SPHCAR_08 (R,THETA,PHI,X,Y,Z,J)
C
C   CONVERTS SPHERICAL COORDS INTO CARTESIAN ONES AND VICE VERSA
C    (THETA AND PHI IN RADIANS).
C
C                  J>0            J<0
C-----INPUT:   J,R,THETA,PHI     J,X,Y,Z
C----OUTPUT:      X,Y,Z        R,THETA,PHI
C
C  NOTE: AT THE POLES (X=0 AND Y=0) WE ASSUME PHI=0 WHEN CONVERTING
C        FROM CARTESIAN TO SPHERICAL COORDS (I.E., FOR J<0)
C
C   LAST MOFIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C
C   AUTHOR:  N. A. TSYGANENKO
C
      IMPLICIT REAL*8 (A-H,O-Z)

      IF(J.GT.0) GOTO 3
      SQ=X**2+Y**2
      R=DSQRT(SQ+Z**2)
      IF (SQ.NE.0.D0) GOTO 2
      PHI=0.D0
      IF (Z.LT.0.D0) GOTO 1
      THETA=0.D0
      RETURN
  1   THETA=3.141592654D0
      RETURN
  2   SQ=DSQRT(SQ)
      PHI=DATAN2(Y,X)
      THETA=DATAN2(SQ,Z)
      IF (PHI.LT.0.D0) PHI=PHI+6.283185307D0
      RETURN
  3   SQ=R*DSIN(THETA)
      X=SQ*DCOS(PHI)
      Y=SQ*DSIN(PHI)
      Z=R*DCOS(THETA)
      RETURN
      END
C
C===========================================================================
c
      SUBROUTINE BSPCAR_08 (THETA,PHI,BR,BTHETA,BPHI,BX,BY,BZ)
C
C   CALCULATES CARTESIAN FIELD COMPONENTS FROM LOCAL SPHERICAL ONES
C
C-----INPUT:   THETA,PHI - SPHERICAL ANGLES OF THE POINT IN RADIANS
C              BR,BTHETA,BPHI -  LOCAL SPHERICAL COMPONENTS OF THE FIELD
C-----OUTPUT:  BX,BY,BZ - CARTESIAN COMPONENTS OF THE FIELD
C
C   LAST MOFIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C
C   WRITTEN BY:  N. A. TSYGANENKO
C
      IMPLICIT REAL*8 (A-H,O-Z)

      S=DSIN(THETA)
      C=DCOS(THETA)
      SF=DSIN(PHI)
      CF=DCOS(PHI)
      BE=BR*S+BTHETA*C
      BX=BE*CF-BPHI*SF
      BY=BE*SF+BPHI*CF
      BZ=BR*C-BTHETA*S
      RETURN
      END
c
C==============================================================================
C
      SUBROUTINE BCARSP_08 (X,Y,Z,BX,BY,BZ,BR,BTHETA,BPHI)
C
CALCULATES LOCAL SPHERICAL FIELD COMPONENTS FROM THOSE IN CARTESIAN SYSTEM
C
C-----INPUT:   X,Y,Z  - CARTESIAN COMPONENTS OF THE POSITION VECTOR
C              BX,BY,BZ - CARTESIAN COMPONENTS OF THE FIELD VECTOR
C-----OUTPUT:  BR,BTHETA,BPHI - LOCAL SPHERICAL COMPONENTS OF THE FIELD VECTOR
C
C  NOTE: AT THE POLES (THETA=0 OR THETA=PI) WE ASSUME PHI=0,
C        AND HENCE BTHETA=BX, BPHI=BY
C
C   WRITTEN AND ADDED TO THIS PACKAGE:  APRIL 1, 2003
c   LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C   AUTHOR:   N. A. TSYGANENKO
C
      IMPLICIT REAL*8 (A-H,O-Z)

      RHO2=X**2+Y**2
      R=DSQRT(RHO2+Z**2)
      RHO=DSQRT(RHO2)

      IF (RHO.NE.0.D0) THEN
        CPHI=X/RHO
        SPHI=Y/RHO
       ELSE
        CPHI=1.D0
        SPHI=0.D0
      ENDIF

      CT=Z/R
      ST=RHO/R

      BR=(X*BX+Y*BY+Z*BZ)/R
      BTHETA=(BX*CPHI+BY*SPHI)*CT-BZ*ST
      BPHI=BY*CPHI-BX*SPHI

      RETURN
      END
C
c=====================================================================================
C
      SUBROUTINE RECALC_08 (IYEAR,IDAY,IHOUR,MIN,ISEC,VGSEX,VGSEY,VGSEZ)
C
C  1. PREPARES ELEMENTS OF ROTATION MATRICES FOR TRANSFORMATIONS OF VECTORS BETWEEN
C     SEVERAL COORDINATE SYSTEMS, MOST FREQUENTLY USED IN SPACE PHYSICS.
C
C  2. PREPARES COEFFICIENTS USED IN THE CALCULATION OF THE MAIN GEOMAGNETIC FIELD
C      (IGRF MODEL)
C
C  THIS SUBROUTINE SHOULD BE INVOKED BEFORE USING THE FOLLOWING SUBROUTINES:
C  IGRF_GEO_08, IGRF_GSW_08, DIP_08, GEOMAG_08, GEOGSW_08, MAGSW_08, SMGSW_08, GSWGSE_08,
c  GEIGEO_08, TRACE_08, STEP_08, RHAND_08.
C
C  THERE IS NO NEED TO REPEATEDLY INVOKE RECALC_08, IF MULTIPLE CALCULATIONS ARE MADE
C    FOR THE SAME DATE/TIME AND SOLAR WIND FLOW DIRECTION.
C
C-----INPUT PARAMETERS:
C
C     IYEAR   -  YEAR NUMBER (FOUR DIGITS)
C     IDAY  -  DAY OF YEAR (DAY 1 = JAN 1)
C     IHOUR -  HOUR OF DAY (00 TO 23)
C     MIN   -  MINUTE OF HOUR (00 TO 59)
C     ISEC  -  SECONDS OF MINUTE (00 TO 59)
C     VGSEX,VGSEY,VGSEZ - GSE (GEOCENTRIC SOLAR-ECLIPTIC) COMPONENTS OF THE OBSERVED
C                              SOLAR WIND FLOW VELOCITY (IN KM/S)
C
C  IMPORTANT: IF ONLY QUESTIONABLE INFORMATION (OR NO INFORMATION AT ALL) IS AVAILABLE
C             ON THE SOLAR WIND SPEED, OR, IF THE STANDARD GSM AND/OR SM COORDINATES ARE
C             INTENDED TO BE USED, THEN SET VGSEX=-400.0 AND VGSEY=VGSEZ=0. IN THIS CASE,
C             THE GSW COORDINATE SYSTEM BECOMES IDENTICAL TO THE STANDARD GSM.
C
C             IF ONLY SCALAR SPEED V OF THE SOLAR WIND IS KNOWN, THEN SETTING
C             VGSEX=-V, VGSEY=29.78, VGSEZ=0.0 WILL TAKE INTO ACCOUNT THE ~4 degs
C             ABERRATION OF THE MAGNETOSPHERE DUE TO EARTH'S ORBITAL MOTION
C
C             IF ALL THREE GSE COMPONENTS OF THE SOLAR WIND VELOCITY ARE AVAILABLE,
C             PLEASE NOTE THAT IN SOME SOLAR WIND DATABASES THE ABERRATION EFFECT
C             HAS ALREADY BEEN TAKEN INTO ACCOUNT BY SUBTRACTING 29.78 KM/S FROM VYGSE;
C             IN THAT CASE, THE UNABERRATED (OBSERVED) VYGSE VALUES SHOULD BE RESTORED
C             BY ADDING BACK THE 29.78 KM/S CORRECTION. WHETHER OR NOT TO DO THAT, MUST
C             BE EITHER VERIFIED WITH THE DATA ORIGINATOR OR DETERMINED BY AVERAGING
C             VGSEY OVER A SUFFICIENTLY LONG TIME INTERVAL.
C
C-----OUTPUT PARAMETERS:  NONE (ALL OUTPUT QUANTITIES ARE PLACED
C                         INTO THE COMMON BLOCKS /GEOPACK1/ AND /GEOPACK2/)
C
C    OTHER SUBROUTINES CALLED BY THIS ONE: SUN_08
C
C    AUTHOR:  N.A. TSYGANENKO
C    DATE:    DEC.1, 1991
C
C    REVISION OF JANUARY 31, 2015:
C
C     The table of IGRF coefficients was extended to include those for the epoch 2015 (igrf-12)
c         (for details, see http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html)
C
C    REVISION OF NOVEMBER 30, 2010:
C
C     The table of IGRF coefficients was extended to include those for the epoch 2010
c         (for details, see http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html)
C
C    REVISION OF NOVEMBER 15, 2007: ADDED THE POSSIBILITY TO TAKE INTO ACCOUNT THE OBSERVED
C     DEFLECTION OF THE SOLAR WIND FLOW FROM STRICTLY RADIAL DIRECTION. TO THAT END, THREE
C     GSE COMPONENTS OF THE SOLAR WIND VELOCITY WERE ADDED TO THE INPUT PARAMETERS.
C
c    CORRECTION OF MAY 9, 2006:  INTERPOLATION OF THE COEFFICIENTS (BETWEEN
C     LABELS 50 AND 105) IS NOW MADE THROUGH THE LAST ELEMENT OF THE ARRAYS
C     G(105)  AND H(105) (PREVIOUSLY MADE ONLY THROUGH N=66, WHICH IN SOME
C     CASES CAUSED RUNTIME ERRORS)
c
C    REVISION OF MAY 3, 2005:
C     The table of IGRF coefficients was extended to include those for the epoch 2005
c       the maximal order of spherical harmonics was also increased up to 13
c         (for details, see http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html)
c
C    REVISION OF APRIL 3, 2003:
c    The code now includes preparation of the model coefficients for the subroutines
c    IGRF_08 and GEOMAG_08. This eliminates the need for the SAVE statements, used
c    in the old versions, making the codes easier and more compiler-independent.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      SAVE ISW
C
      COMMON /GEOPACK1/ ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI,
     * SPS,CPS,DS3,CGST,SGST,PSI,A11,A21,A31,A12,A22,A32,A13,A23,A33,
     * E11,E21,E31,E12,E22,E32,E13,E23,E33
C
C  THE COMMON BLOCK /GEOPACK1/ CONTAINS ELEMENTS OF THE ROTATION MATRICES AND OTHER
C   PARAMETERS RELATED TO THE COORDINATE TRANSFORMATIONS PERFORMED BY THIS PACKAGE
C
      COMMON /GEOPACK2/ G(105),H(105),REC(105)
C
C  THE COMMON BLOCK /GEOPACK2/ CONTAINS COEFFICIENTS OF THE IGRF FIELD MODEL, CALCULATED
C    FOR A GIVEN YEAR AND DAY FROM THEIR STANDARD EPOCH VALUES. THE ARRAY REC CONTAINS
C    COEFFICIENTS USED IN THE RECURSION RELATIONS FOR LEGENDRE ASSOCIATE POLYNOMIALS.
C
      DIMENSION G65(105),H65(105),G70(105),H70(105),G75(105),H75(105),
     + G80(105),H80(105),G85(105),H85(105),G90(105),H90(105),G95(105),
     + H95(105),G00(105),H00(105),G05(105),H05(105),G10(105),H10(105),
     + G15(105),H15(105),DG15(45),DH15(45)
C
      DATA ISW /0/
c
      DATA G65/0.D0,-30334.D0,-2119.D0,-1662.D0,2997.D0,1594.D0,1297.D0,
     *-2038.D0,1292.D0,856.D0,957.D0,804.D0,479.D0,-390.D0,252.D0,
     *-219.D0,358.D0,254.D0,-31.D0,-157.D0,-62.D0,45.D0,61.D0,8.D0,
     *-228.D0,4.D0,1.D0,-111.D0,75.D0,-57.D0,4.D0,13.D0,-26.D0,-6.D0,
     *13.D0,1.D0,13.D0,5.D0,-4.D0,-14.D0,0.D0,8.D0,-1.D0,11.D0,4.D0,
     *8.D0,10.D0,2.D0,-13.D0,10.D0,-1.D0,-1.D0,5.D0,1.D0,-2.D0,-2.D0,
     *-3.D0,2.D0,-5.D0,-2.D0,4.D0,4.D0,0.D0,2.D0,2.D0,0.D0,39*0.D0/
      DATA H65/0.D0,0.D0,5776.D0,0.D0,-2016.D0,114.D0,0.D0,-404.D0,
     *240.D0,-165.D0,0.D0,148.D0,-269.D0,13.D0,-269.D0,0.D0,19.D0,
     *128.D0,-126.D0,-97.D0,81.D0,0.D0,-11.D0,100.D0,68.D0,-32.D0,-8.D0,
     *-7.D0,0.D0,-61.D0,-27.D0,-2.D0,6.D0,26.D0,-23.D0,-12.D0,0.D0,7.D0,
     *-12.D0,9.D0,-16.D0,4.D0,24.D0,-3.D0,-17.D0,0.D0,-22.D0,15.D0,7.D0,
     *-4.D0,-5.D0,10.D0,10.D0,-4.D0,1.D0,0.D0,2.D0,1.D0,2.D0,6.D0,-4.D0,
     *0.D0,-2.D0,3.D0,0.D0,-6.D0,39*0.D0/
c
      DATA G70/0.D0,-30220.D0,-2068.D0,-1781.D0,3000.D0,1611.D0,1287.D0,
     *-2091.D0,1278.D0,838.D0,952.D0,800.D0,461.D0,-395.D0,234.D0,
     *-216.D0,359.D0,262.D0,-42.D0,-160.D0,-56.D0,43.D0,64.D0,15.D0,
     *-212.D0,2.D0,3.D0,-112.D0,72.D0,-57.D0,1.D0,14.D0,-22.D0,-2.D0,
     *13.D0,-2.D0,14.D0,6.D0,-2.D0,-13.D0,-3.D0,5.D0,0.D0,11.D0,3.D0,
     *8.D0,10.D0,2.D0,-12.D0,10.D0,-1.D0,0.D0,3.D0,1.D0,-1.D0,-3.D0,
     *-3.D0,2.D0,-5.D0,-1.D0,6.D0,4.D0,1.D0,0.D0,3.D0,-1.D0,39*0.D0/
      DATA H70/0.D0,0.D0,5737.D0,0.D0,-2047.D0,25.D0,0.D0,-366.D0,
     *251.D0,-196.D0,0.D0,167.D0,-266.D0,26.D0,-279.D0,0.D0,26.D0,
     *139.D0,-139.D0,-91.D0,83.D0,0.D0,-12.D0,100.D0,72.D0,-37.D0,-6.D0,
     *1.D0,0.D0,-70.D0,-27.D0,-4.D0,8.D0,23.D0,-23.D0,-11.D0,0.D0,7.D0,
     *-15.D0,6.D0,-17.D0,6.D0,21.D0,-6.D0,-16.D0,0.D0,-21.D0,16.D0,6.D0,
     *-4.D0,-5.D0,10.D0,11.D0,-2.D0,1.D0,0.D0,1.D0,1.D0,3.D0,4.D0,-4.D0,
     *0.D0,-1.D0,3.D0,1.D0,-4.D0,39*0.D0/
c
      DATA G75/0.D0,-30100.D0,-2013.D0,-1902.D0,3010.D0,1632.D0,1276.D0,
     *-2144.D0,1260.D0,830.D0,946.D0,791.D0,438.D0,-405.D0,216.D0,
     *-218.D0,356.D0,264.D0,-59.D0,-159.D0,-49.D0,45.D0,66.D0,28.D0,
     *-198.D0,1.D0,6.D0,-111.D0,71.D0,-56.D0,1.D0,16.D0,-14.D0,0.D0,
     *12.D0,-5.D0,14.D0,6.D0,-1.D0,-12.D0,-8.D0,4.D0,0.D0,10.D0,1.D0,
     *7.D0,10.D0,2.D0,-12.D0,10.D0,-1.D0,-1.D0,4.D0,1.D0,-2.D0,-3.D0,
     *-3.D0,2.D0,-5.D0,-2.D0,5.D0,4.D0,1.D0,0.D0,3.D0,-1.D0,39*0.D0/
C
      DATA H75/0.D0,0.D0,5675.D0,0.D0,-2067.D0,-68.D0,0.D0,-333.D0,
     *262.D0,-223.D0,0.D0,191.D0,-265.D0,39.D0,-288.D0,0.D0,31.D0,
     *148.D0,-152.D0,-83.D0,88.D0,0.D0,-13.D0,99.D0,75.D0,-41.D0,-4.D0,
     *11.D0,0.D0,-77.D0,-26.D0,-5.D0,10.D0,22.D0,-23.D0,-12.D0,0.D0,
     *6.D0,-16.D0,4.D0,-19.D0,6.D0,18.D0,-10.D0,-17.D0,0.D0,-21.D0,
     *16.D0,7.D0,-4.D0,-5.D0,10.D0,11.D0,-3.D0,1.D0,0.D0,1.D0,1.D0,3.D0,
     *4.D0,-4.D0,-1.D0,-1.D0,3.D0,1.D0,-5.D0,39*0.D0/
c
      DATA G80/0.D0,-29992.D0,-1956.D0,-1997.D0,3027.D0,1663.D0,1281.D0,
     *-2180.D0,1251.D0,833.D0,938.D0,782.D0,398.D0,-419.D0,199.D0,
     *-218.D0,357.D0,261.D0,-74.D0,-162.D0,-48.D0,48.D0,66.D0,42.D0,
     *-192.D0,4.D0,14.D0,-108.D0,72.D0,-59.D0,2.D0,21.D0,-12.D0,1.D0,
     *11.D0,-2.D0,18.D0,6.D0,0.D0,-11.D0,-7.D0,4.D0,3.D0,6.D0,-1.D0,
     *5.D0,10.D0,1.D0,-12.D0,9.D0,-3.D0,-1.D0,7.D0,2.D0,-5.D0,-4.D0,
     *-4.D0,2.D0,-5.D0,-2.D0,5.D0,3.D0,1.D0,2.D0,3.D0,0.D0,39*0.D0/
C
      DATA H80/0.D0,0.D0,5604.D0,0.D0,-2129.D0,-200.D0,0.D0,-336.D0,
     *271.D0,-252.D0,0.D0,212.D0,-257.D0,53.D0,-297.D0,0.D0,46.D0,
     *150.D0,-151.D0,-78.D0,92.D0,0.D0,-15.D0,93.D0,71.D0,-43.D0,-2.D0,
     *17.D0,0.D0,-82.D0,-27.D0,-5.D0,16.D0,18.D0,-23.D0,-10.D0,0.D0,
     *7.D0,-18.D0,4.D0,-22.D0,9.D0,16.D0,-13.D0,-15.D0,0.D0,-21.D0,
     *16.D0,9.D0,-5.D0,-6.D0,9.D0,10.D0,-6.D0,2.D0,0.D0,1.D0,0.D0,3.D0,
     *6.D0,-4.D0,0.D0,-1.D0,4.D0,0.D0,-6.D0,39*0.D0/
c
      DATA G85/0.D0,-29873.D0,-1905.D0,-2072.D0,3044.D0,1687.D0,1296.D0,
     *-2208.D0,1247.D0,829.D0,936.D0,780.D0,361.D0,-424.D0,170.D0,
     *-214.D0,355.D0,253.D0,-93.D0,-164.D0,-46.D0,53.D0,65.D0,51.D0,
     *-185.D0,4.D0,16.D0,-102.D0,74.D0,-62.D0,3.D0,24.D0,-6.D0,4.D0,
     *10.D0,0.D0,21.D0,6.D0,0.D0,-11.D0,-9.D0,4.D0,4.D0,4.D0,-4.D0,5.D0,
     *10.D0,1.D0,-12.D0,9.D0,-3.D0,-1.D0,7.D0,1.D0,-5.D0,-4.D0,-4.D0,
     *3.D0,-5.D0,-2.D0,5.D0,3.D0,1.D0,2.D0,3.D0,0.D0,39*0.D0/
C
      DATA H85/0.D0,0.D0,5500.D0,0.D0,-2197.D0,-306.D0,0.D0,-310.D0,
     *284.D0,-297.D0,0.D0,232.D0,-249.D0,69.D0,-297.D0,0.D0,47.D0,
     *150.D0,-154.D0,-75.D0,95.D0,0.D0,-16.D0,88.D0,69.D0,-48.D0,-1.D0,
     *21.D0,0.D0,-83.D0,-27.D0,-2.D0,20.D0,17.D0,-23.D0,-7.D0,0.D0,8.D0,
     *-19.D0,5.D0,-23.D0,11.D0,14.D0,-15.D0,-11.D0,0.D0,-21.D0,15.D0,
     *9.D0,-6.D0,-6.D0,9.D0,9.D0,-7.D0,2.D0,0.D0,1.D0,0.D0,3.D0,6.D0,
     *-4.D0,0.D0,-1.D0,4.D0,0.D0,-6.D0,39*0.D0/
c
      DATA G90/0.D0,-29775.D0,-1848.D0,-2131.D0,3059.D0,1686.D0,1314.D0,
     *     -2239.D0,  1248.D0,  802.D0,  939.D0, 780.D0, 325.D0,-423.D0,
     *       141.D0,  -214.D0,  353.D0,  245.D0,-109.D0,-165.D0, -36.D0,
     *        61.D0,    65.D0,   59.D0, -178.D0,   3.D0,  18.D0, -96.D0,
     *        77.D0,   -64.D0,    2.D0,   26.D0,  -1.D0,   5.D0,   9.D0,
     *         0.D0,    23.D0,    5.D0,   -1.D0, -10.D0, -12.D0,   3.D0,
     *         4.D0,     2.D0,   -6.D0,    4.D0,   9.D0,   1.D0, -12.D0,
     *         9.D0,    -4.D0,   -2.D0,    7.D0,   1.D0,  -6.D0,  -3.D0,
     *        -4.D0,     2.D0,   -5.D0,   -2.D0,   4.D0,   3.D0,   1.D0,
     *         3.D0,     3.D0,    0.D0,39*0.D0/

      DATA H90/0.D0,  0.D0,5406.D0,   0.D0,-2279.D0,-373.D0,  0.D0,
     *      -284.D0,293.D0,-352.D0,   0.D0,  247.D0,-240.D0, 84.D0,
     *      -299.D0,  0.D0,  46.D0, 154.D0, -153.D0, -69.D0, 97.D0,
     *         0.D0,-16.D0,  82.D0,  69.D0,  -52.D0,   1.D0, 24.D0,
     *         0.D0,-80.D0, -26.D0,   0.D0,   21.D0,  17.D0,-23.D0,
     *        -4.D0,  0.D0,  10.D0, -19.D0,    6.D0, -22.D0, 12.D0,
     *        12.D0,-16.D0, -10.D0,   0.D0,  -20.D0,  15.D0, 11.D0,
     *        -7.D0, -7.D0,   9.D0,   8.D0,   -7.D0,   2.D0,  0.D0,
     *         2.D0,  1.D0,   3.D0,   6.D0,   -4.D0,   0.D0, -2.D0,
     *         3.D0, -1.D0,  -6.D0,39*0.D0/

      DATA G95/0.D0,-29692.D0,-1784.D0,-2200.D0,3070.D0,1681.D0,1335.D0,
     *     -2267.D0,  1249.D0,  759.D0,  940.D0, 780.D0, 290.D0,-418.D0,
     *       122.D0,  -214.D0,  352.D0,  235.D0,-118.D0,-166.D0, -17.D0,
     *        68.D0,    67.D0,   68.D0, -170.D0,  -1.D0,  19.D0, -93.D0,
     *        77.D0,   -72.D0,    1.D0,   28.D0,   5.D0,   4.D0,   8.D0,
     *        -2.D0,    25.D0,    6.D0,   -6.D0,  -9.D0, -14.D0,   9.D0,
     *         6.D0,    -5.D0,   -7.D0,    4.D0,   9.D0,   3.D0, -10.D0,
     *         8.D0,    -8.D0,   -1.D0,   10.D0,  -2.D0,  -8.D0,  -3.D0,
     *        -6.D0,     2.D0,   -4.D0,   -1.D0,   4.D0,   2.D0,   2.D0,
     *         5.D0,     1.D0,    0.D0,  39*0.D0/

      DATA H95/0.D0,  0.D0,5306.D0,  0.D0,-2366.D0,-413.D0,  0.D0,
     *      -262.D0,302.D0,-427.D0,  0.D0,  262.D0,-236.D0, 97.D0,
     *      -306.D0,  0.D0,  46.D0,165.D0, -143.D0, -55.D0,107.D0,
     *         0.D0,-17.D0,  72.D0, 67.D0,  -58.D0,   1.D0, 36.D0,
     *         0.D0,-69.D0, -25.D0,  4.D0,   24.D0,  17.D0,-24.D0,
     *        -6.D0,  0.D0,  11.D0,-21.D0,    8.D0, -23.D0, 15.D0,
     *        11.D0,-16.D0,  -4.D0,  0.D0,  -20.D0,  15.D0, 12.D0,
     *        -6.D0, -8.D0,   8.D0,  5.D0,   -8.D0,   3.D0,  0.D0,
     *         1.D0,  0.D0,   4.D0,  5.D0,   -5.D0,  -1.D0, -2.D0,
     *         1.D0, -2.D0,  -7.D0,39*0.D0/

      DATA G00/0.D0,-29619.4D0,-1728.2D0,-2267.7D0,3068.4D0,1670.9D0,
     *     1339.6D0,  -2288.D0, 1252.1D0,  714.5D0, 932.3D0, 786.8D0,
     *       250.D0,   -403.D0,  111.3D0, -218.8D0, 351.4D0, 222.3D0,
     *     -130.4D0,  -168.6D0,  -12.9D0,   72.3D0,  68.2D0,  74.2D0,
     *     -160.9D0,    -5.9D0,   16.9D0,  -90.4D0,  79.0D0, -74.0D0,
     *         0.D0,    33.3D0,    9.1D0,    6.9D0,   7.3D0,  -1.2D0,
     *       24.4D0,     6.6D0,   -9.2D0,   -7.9D0, -16.6D0,   9.1D0,
     *        7.0D0,    -7.9D0,    -7.D0,     5.D0,   9.4D0,    3.D0,
     *      - 8.4D0,     6.3D0,   -8.9D0,   -1.5D0,   9.3D0,  -4.3D0,
     *       -8.2D0,    -2.6D0,    -6.D0,    1.7D0,  -3.1D0,  -0.5D0,
     *        3.7D0,      1.D0,     2.D0,    4.2D0,   0.3D0,  -1.1D0,
     *        2.7D0,    -1.7D0,   -1.9D0,    1.5D0,  -0.1D0,   0.1D0,
     *       -0.7D0,     0.7D0,    1.7D0,    0.1D0,   1.2D0,   4.0D0,
     *       -2.2D0,    -0.3D0,    0.2D0,    0.9D0,  -0.2D0,   0.9D0,
     *       -0.5D0,     0.3D0,   -0.3D0,   -0.4D0,  -0.1D0,  -0.2D0,
     *       -0.4D0,    -0.2D0,   -0.9D0,    0.3D0,   0.1D0,  -0.4D0,
     *        1.3D0,    -0.4D0,    0.7D0,   -0.4D0,   0.3D0,  -0.1D0,
     *        0.4D0,      0.D0,    0.1D0/

      DATA H00/0.D0,   0.D0,5186.1D0,   0.D0,-2481.6D0,-458.0D0,   0.D0,
     *     -227.6D0,293.4D0,-491.1D0,   0.D0,  272.6D0,-231.9D0,119.8D0,
     *     -303.8D0,   0.D0,  43.8D0,171.9D0, -133.1D0, -39.3D0,106.3D0,
     *         0.D0,-17.4D0,  63.7D0, 65.1D0,  -61.2D0,   0.7D0, 43.8D0,
     *         0.D0,-64.6D0, -24.2D0,  6.2D0,    24.D0,  14.8D0,-25.4D0,
     *       -5.8D0,  0.0D0,  11.9D0,-21.5D0,    8.5D0, -21.5D0, 15.5D0,
     *        8.9D0,-14.9D0,  -2.1D0,  0.0D0,  -19.7D0,  13.4D0, 12.5D0,
     *       -6.2D0, -8.4D0,   8.4D0,  3.8D0,   -8.2D0,   4.8D0,  0.0D0,
     *        1.7D0,  0.0D0,   4.0D0,  4.9D0,   -5.9D0,  -1.2D0, -2.9D0,
     *        0.2D0, -2.2D0,  -7.4D0,  0.0D0,    0.1D0,   1.3D0, -0.9D0,
     *       -2.6D0,  0.9D0,  -0.7D0, -2.8D0,   -0.9D0,  -1.2D0, -1.9D0,
     *       -0.9D0,  0.0D0,  -0.4D0,  0.3D0,    2.5D0,  -2.6D0,  0.7D0,
     *        0.3D0,  0.0D0,   0.0D0,  0.3D0,   -0.9D0,  -0.4D0,  0.8D0,
     *        0.0D0, -0.9D0,   0.2D0,  1.8D0,   -0.4D0,  -1.0D0, -0.1D0,
     *        0.7D0,  0.3D0,   0.6D0,  0.3D0,   -0.2D0,  -0.5D0, -0.9D0/

      DATA G05/0.D0,  -29554.6D0,-1669.0D0,-2337.2D0,3047.7D0,1657.8D0,
     *     1336.3D0,   -2305.8D0, 1246.4D0,  672.5D0, 920.6D0, 798.0D0,
     *      210.7D0,    -379.9D0,  100.0D0, -227.0D0, 354.4D0, 208.9D0,
     *     -136.5D0,    -168.1D0,  -13.6D0,   73.6D0,  69.6D0,  76.7D0,
     *     -151.3D0,     -14.6D0,   14.6D0,  -86.4D0,  79.9D0, -74.5D0,
     *       -1.7D0,      38.7D0,   12.3D0,    9.4D0,   5.4D0,   1.9D0,
     *       24.8D0,       7.6D0,  -11.7D0,   -6.9D0, -18.1D0,  10.2D0,
     *        9.4D0,     -11.3D0,   -4.9D0,    5.6D0,   9.8D0,   3.6D0,
     *       -6.9D0,       5.0D0,  -10.8D0,   -1.3D0,   8.8D0,  -6.7D0,
     *       -9.2D0,      -2.2D0,   -6.1D0,    1.4D0,  -2.4D0,  -0.2D0,
     *        3.1D0,       0.3D0,    2.1D0,    3.8D0,  -0.2D0,  -2.1D0,
     *        2.9D0,      -1.6D0,   -1.9D0,    1.4D0,  -0.3D0,   0.3D0,
     *       -0.8D0,       0.5D0,    1.8D0,    0.2D0,   1.0D0,   4.0D0,
     *       -2.2D0,      -0.3D0,    0.2D0,    0.9D0,  -0.4D0,   1.0D0,
     *       -0.3D0,       0.5D0,   -0.4D0,   -0.4D0,   0.1D0,  -0.5D0,
     *       -0.1D0,      -0.2D0,   -0.9D0,    0.3D0,   0.3D0,  -0.4D0,
     *        1.2D0,      -0.4D0,    0.8D0,   -0.3D0,   0.4D0,  -0.1D0,
     *        0.4D0,      -0.1D0,   -0.2D0/
C
      DATA H05/0.D0,  0.0D0,5078.0D0,  0.0D0,-2594.5D0,-515.4D0,  0.0D0,
     *     -198.9D0,269.7D0,-524.7D0,  0.0D0,  282.1D0,-225.2D0,145.2D0,
     *     -305.4D0,  0.0D0,  42.7D0,180.3D0, -123.5D0, -19.6D0,103.9D0,
     *        0.0D0,-20.3D0,  54.8D0, 63.6D0,  -63.5D0,   0.2D0, 50.9D0,
     *        0.0D0,-61.1D0, -22.6D0,  6.8D0,   25.4D0,  10.9D0,-26.3D0,
     *       -4.6D0,  0.0D0,  11.2D0,-20.9D0,    9.8D0, -19.7D0, 16.2D0,
     *        7.6D0,-12.8D0,  -0.1D0,  0.0D0,  -20.1D0,  12.7D0, 12.7D0,
     *       -6.7D0, -8.2D0,   8.1D0,  2.9D0,   -7.7D0,   6.0D0,  0.0D0,
     *        2.2D0,  0.1D0,   4.5D0,  4.8D0,   -6.7D0,  -1.0D0, -3.5D0,
     *       -0.9D0, -2.3D0,  -7.9D0,  0.0D0,    0.3D0,   1.4D0, -0.8D0,
     *       -2.3D0,  0.9D0,  -0.6D0, -2.7D0,   -1.1D0,  -1.6D0, -1.9D0,
     *       -1.4D0,  0.0D0,  -0.6D0,  0.2D0,    2.4D0,  -2.6D0,  0.6D0,
     *        0.4D0,  0.0D0,   0.0D0,  0.3D0,   -0.9D0,  -0.3D0,  0.9D0,
     *        0.0D0, -0.8D0,   0.3D0,  1.7D0,   -0.5D0,  -1.1D0,  0.0D0,
     *        0.6D0,  0.2D0,   0.5D0,  0.4D0,   -0.2D0,  -0.6D0, -0.9D0/
C
      DATA G10/0.00D0,-29496.57D0,-1586.42D0,-2396.06D0,3026.34D0,
     *      1668.17D0,  1339.85D0,-2326.54D0, 1232.10D0, 633.73D0, 
     *       912.66D0,   808.97D0,  166.58D0, -356.83D0,  89.40D0,
     *      -230.87D0,   357.29D0,  200.26D0, -141.05D0,-163.17D0,
     *        -8.03D0,    72.78D0,   68.69D0,   75.92D0,-141.40D0, 
     *       -22.83D0,    13.10D0,  -78.09D0,   80.44D0, -75.00D0,
     *        -4.55D0,    45.24D0,   14.00D0,   10.46D0,   1.64D0,
     *         4.92D0,    24.41D0,    8.21D0,  -14.50D0,  -5.59D0,
     *       -19.34D0,    11.61D0,   10.85D0,  -14.05D0,  -3.54D0,
     *         5.50D0,     9.45D0,    3.45D0,   -5.27D0,   3.13D0,
     *       -12.38D0,    -0.76D0,    8.43D0,   -8.42D0, -10.08D0,
     *        -1.94D0,    -6.24D0,    0.89D0,   -1.07D0,  -0.16D0,
     *         2.45D0,    -0.33D0,    2.13D0,    3.09D0,  -1.03D0,
     *        -2.80D0,     3.05D0,   -1.48D0,   -2.03D0,   1.65D0,
     *        -0.51D0,     0.54D0,   -0.79D0,    0.37D0,   1.79D0,
     *         0.12D0,     0.75D0,    3.75D0,   -2.12D0,  -0.21D0,
     *         0.30D0,     1.04D0,   -0.63D0,    0.95D0,  -0.11D0,
     *         0.52D0,    -0.39D0,   -0.37D0,    0.21D0,  -0.77D0,
     *         0.04D0,    -0.09D0,   -0.89D0,    0.31D0,   0.42D0,
     *        -0.45D0,     1.08D0,   -0.31D0,    0.78D0,  -0.18D0,
     *         0.38D0,     0.02D0,    0.42D0,   -0.26D0,  -0.26D0/
C
      DATA H10/0.00D0,  0.00D0,4944.26D0,   0.00D0,-2708.54D0,
     *     -575.73D0,   0.00D0,-160.40D0, 251.75D0, -537.03D0, 0.00D0,
     *      286.48D0,-211.03D0, 164.46D0,-309.72D0,    0.00D0,44.58D0,
     *      189.01D0,-118.06D0,  -0.01D0, 101.04D0,    0.00D0,-20.90D0,
     *       44.18D0,  61.54D0, -66.26D0,   3.02D0,   55.40D0,  0.00D0,
     *      -57.80D0, -21.20D0,   6.54D0,  24.96D0,    7.03D0,-27.61D0,
     *       -3.28D0,   0.00D0,  10.84D0, -20.03D0,   11.83D0,-17.41D0, 
     *       16.71D0,   6.96D0, -10.74D0,   1.64D0,    0.00D0,-20.54D0,
     *       11.51D0,  12.75D0,  -7.14D0,  -7.42D0,    7.97D0,  2.14D0,
     *       -6.08D0,   7.01D0,   0.00D0,   2.73D0,   -0.10D0,  4.71D0,
     *        4.44D0,  -7.22D0,  -0.96D0,  -3.95D0,   -1.99D0, -1.97D0,
     *       -8.31D0,   0.00D0,   0.13D0,   1.67D0,   -0.66D0, -1.76D0,
     *        0.85D0,  -0.39D0,  -2.51D0,  -1.27D0,   -2.11D0, -1.94D0,
     *       -1.86D0,   0.00D0,  -0.87D0,   0.27D0,    2.13D0, -2.49D0,
     *        0.49D0,   0.59D0,   0.00D0,   0.13D0,    0.27D0, -0.86D0,
     *       -0.23D0,   0.87D0,   0.00D0,  -0.87D0,    0.30D0,  1.66D0,
     *       -0.59D0,  -1.14D0,  -0.07D0,   0.54D0,    0.10D0,  0.49D0,
     *        0.44D0,  -0.25D0,  -0.53D0,  -0.79D0/
C
      DATA G15/0.D0,-29442.0D0, -1501.0D0, -2445.1D0,  3012.9D0,
     *     1676.7D0,  1350.7D0, -2352.3D0,  1225.6D0,   582.0D0,907.6D0, 
     *      813.7D0,   120.4D0,  -334.9D0,    70.4D0,  -232.6D0,360.1D0,
     *      192.4D0,  -140.9D0,  -157.5D0,     4.1D0,    70.0D0, 67.7D0,
     *       72.7D0,  -129.9D0,   -28.9D0,    13.2D0,   -70.9D0, 81.6D0,
     *      -76.1D0,    -6.8D0,    51.8D0,    15.0D0,     9.4D0, -2.8D0,
     *        6.8D0,    24.2D0,     8.8D0,   -16.9D0,    -3.2D0,-20.6D0,
     *       13.4D0,    11.7D0,   -15.9D0,    -2.0D0,     5.4D0,  8.8D0, 
     *        3.1D0,    -3.3D0,     0.7D0,   -13.3D0,    -0.1D0,  8.7D0,
     *       -9.1D0,   -10.5D0,    -1.9D0,    -6.3D0,     0.1D0,  0.5D0,
     *       -0.5D0,     1.8D0,    -0.7D0,     2.1D0,     2.4D0, -1.8D0,
     *       -3.6D0,     3.1D0,    -1.5D0,    -2.3D0,     2.0D0, -0.8D0,
     *        0.6D0,    -0.7D0,     0.2D0,     1.7D0,    -0.2D0,  0.4D0,
     *        3.5D0,    -1.9D0,    -0.2D0,     0.4D0,     1.2D0, -0.8D0,
     *        0.9D0,     0.1D0,     0.5D0,    -0.3D0,    -0.4D0,  0.2D0,
     *       -0.9D0,     0.0D0,     0.0D0,    -0.9D0,     0.4D0,  0.5D0,
     *       -0.5D0,     1.0D0,    -0.2D0,     0.8D0,    -0.1D0,  0.3D0,
     *        0.1D0,     0.5D0,    -0.4D0,    -0.3D0/
c
      DATA H15/0.D0,     0.0D0, 4797.1D0,     0.0D0, -2845.6D0,-641.9D0,
     *        0.0D0,  -115.3D0,  244.9D0,  -538.4D0,     0.0D0, 283.3D0,
     *     -188.7D0,   180.9D0, -329.5D0,     0.0D0,    47.3D0, 197.0D0,
     *     -119.3D0,    16.0D0,  100.2D0,     0.0D0,   -20.8D0,  33.2D0,
     *       58.9D0,   -66.7D0,    7.3D0,    62.6D0,     0.0D0, -54.1D0,
     *      -19.5D0,     5.7D0,   24.4D0,     3.4D0,   -27.4D0,  -2.2D0,
     *        0.0D0,    10.1D0,  -18.3D0,    13.3D0,   -14.6D0,  16.2D0,
     *        5.7D0,    -9.1D0,    2.1D0,     0.0D0,   -21.6D0,  10.8D0,
     *       11.8D0,    -6.8D0,   -6.9D0,     7.8D0,     1.0D0,  -4.0D0,
     *        8.4D0,     0.0D0,    3.2D0,    -0.4D0,     4.6D0,   4.4D0,
     *       -7.9D0,    -0.6D0,   -4.2D0,    -2.8D0,    -1.2D0,  -8.7D0,
     *        0.0D0,    -0.1D0,    2.0D0,    -0.7D0,    -1.1D0,   0.8D0,
     *       -0.2D0,    -2.2D0,   -1.4D0,    -2.5D0,    -2.0D0,  -2.4D0,
     *        0.0D0,    -1.1D0,    0.4D0,     1.9D0,    -2.2D0,   0.3D0,
     *        0.7D0,    -0.1D0,    0.3D0,     0.2D0,    -0.9D0,  -0.1D0,
     *        0.7D0,     0.0D0,   -0.9D0,     0.4D0,     1.6D0,  -0.5D0,
     *       -1.2D0,    -0.1D0,    0.4D0,    -0.1D0,     0.4D0,   0.5D0,
     *       -0.3D0,    -0.4D0,   -0.8D0/
c
      DATA DG15/0.0D0, 10.3D0,  18.1D0,  -8.7D0,  -3.3D0,  2.1D0, 3.4D0,
     *         -5.5D0, -0.7D0, -10.1D0,  -0.7D0,   0.2D0, -9.1D0, 4.1D0,
     *         -4.3D0, -0.2D0,   0.5D0,  -1.3D0,  -0.1D0,  1.4D0, 3.9D0,
     *         -0.3D0, -0.1D0,  -0.7D0,   2.1D0,  -1.2D0,  0.3D0, 1.6D0,
     *          0.3D0, -0.2D0,  -0.5D0,   1.3D0,   0.1D0, -0.6D0,-0.8D0,
     *          0.2D0,  0.2D0,   0.0D0,  -0.6D0,   0.5D0, -0.2D0, 0.4D0,
     *          0.1D0, -0.4D0,   0.3D0/
c
      DATA DH15/0.0D0,  0.0D0, -26.6D0,   0.0D0, -27.4D0,-14.1D0, 0.0D0,
     *          8.2D0, -0.4D0,   1.8D0,   0.0D0,  -1.3D0,  5.3D0, 2.9D0,
     *         -5.2D0,  0.0D0,   0.6D0,   1.7D0,  -1.2D0,  3.4D0, 0.0D0,
     *          0.0D0,  0.0D0,  -2.1D0,  -0.7D0,   0.2D0,  0.9D0, 1.0D0,
     *          0.0D0,  0.8D0,   0.4D0,  -0.2D0,  -0.3D0, -0.6D0, 0.1D0,
     *         -0.2D0,  0.0D0,  -0.3D0,   0.3D0,   0.1D0,  0.5D0,-0.2D0,
     *         -0.3D0,  0.3D0,   0.0D0/
C
C
c      IF (VGSEY.EQ.0..AND.VGSEZ.EQ.0..AND.ISW.NE.1) THEN
c      PRINT *, ''
c      PRINT *,
c     *' RECALC_08: RADIAL SOLAR WIND --> GSW SYSTEM IDENTICAL HERE'
c      PRINT *,
c     *' TO STANDARD GSM (I.E., XGSW AXIS COINCIDES WITH EARTH-SUN LINE)'
c      PRINT *, ''
c      ISW=1
c      ENDIF

c      IF ((VGSEY.NE.0.D0.OR.VGSEZ.NE.0.D0).AND.ISW.NE.2) THEN       !  CORRECTED DEC.01, 2010
c      PRINT *, ''
c      PRINT *,
c     *' WARNING: NON-RADIAL SOLAR WIND FLOW SPECIFIED IN RECALC_08;'
c      PRINT *,
c     *' HENCE XGSW AXIS IS ASSUMED ORIENTED ANTIPARALLEL TO V_SW VECTOR'
c      PRINT *, ''
c      ISW=2
c      ENDIF
C
      IY=IYEAR
C
C  WE ARE RESTRICTED BY THE INTERVAL 1965-2020, FOR WHICH EITHER THE IGRF/DGRF COEFFICIENTS OR SECULAR VELOCITIES
c    ARE KNOWN; IF IYEAR IS OUTSIDE THIS INTERVAL, THEN THE SUBROUTINE USES THE
C      NEAREST LIMITING VALUE AND PRINTS A WARNING:
C
      IF(IY.LT.1965) THEN
       IY=1965
       WRITE (*,10) IYEAR,IY
      ENDIF

      IF(IY.GT.2020) THEN
       IY=2020
       WRITE (*,10) IYEAR,IY
      ENDIF
C
C  CALCULATE THE ARRAY REC, CONTAINING COEFFICIENTS FOR THE RECURSION RELATIONS,
C  USED IN THE IGRF SUBROUTINE FOR CALCULATING THE ASSOCIATE LEGENDRE POLYNOMIALS
C  AND THEIR DERIVATIVES:
c
      DO 20 N=1,14
         N2=2*N-1
         N2=N2*(N2-2)
         DO 20 M=1,N
            MN=N*(N-1)/2+M
20    REC(MN)=DFLOAT((N-M)*(N+M-2))/DFLOAT(N2)
C
      IF (IY.LT.1970) GOTO 50          !INTERPOLATE BETWEEN 1965 - 1970
      IF (IY.LT.1975) GOTO 60          !INTERPOLATE BETWEEN 1970 - 1975
      IF (IY.LT.1980) GOTO 70          !INTERPOLATE BETWEEN 1975 - 1980
      IF (IY.LT.1985) GOTO 80          !INTERPOLATE BETWEEN 1980 - 1985
      IF (IY.LT.1990) GOTO 90          !INTERPOLATE BETWEEN 1985 - 1990
      IF (IY.LT.1995) GOTO 100         !INTERPOLATE BETWEEN 1990 - 1995
      IF (IY.LT.2000) GOTO 110         !INTERPOLATE BETWEEN 1995 - 2000
      IF (IY.LT.2005) GOTO 120         !INTERPOLATE BETWEEN 2000 - 2005
      IF (IY.LT.2010) GOTO 130         !INTERPOLATE BETWEEN 2005 - 2010
      IF (IY.LT.2015) GOTO 140         !INTERPOLATE BETWEEN 2010 - 2015
C
C       EXTRAPOLATE BEYOND 2015:
C
      DT=DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-2015.D0
      DO 40 N=1,105
         G(N)=G15(N)
         H(N)=H15(N)
         IF (N.GT.45) GOTO 40
         G(N)=G(N)+DG15(N)*DT
         H(N)=H(N)+DH15(N)*DT
40    CONTINUE
      GOTO 300
C
C       INTERPOLATE BETWEEEN 1965 - 1970:
C
50    F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-1965)/5.D0
      F1=1.D0-F2
      DO 55 N=1,105
         G(N)=G65(N)*F1+G70(N)*F2
55       H(N)=H65(N)*F1+H70(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 1970 - 1975:
C
60    F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-1970)/5.D0
      F1=1.D0-F2
      DO 65 N=1,105
         G(N)=G70(N)*F1+G75(N)*F2
65       H(N)=H70(N)*F1+H75(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 1975 - 1980:
C
70    F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-1975)/5.D0
      F1=1.D0-F2
      DO 75 N=1,105
         G(N)=G75(N)*F1+G80(N)*F2
75       H(N)=H75(N)*F1+H80(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 1980 - 1985:
C
80    F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-1980)/5.D0
      F1=1.D0-F2
      DO 85 N=1,105
         G(N)=G80(N)*F1+G85(N)*F2
85       H(N)=H80(N)*F1+H85(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 1985 - 1990:
C
90    F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-1985)/5.D0
      F1=1.D0-F2
      DO 95 N=1,105
         G(N)=G85(N)*F1+G90(N)*F2
95       H(N)=H85(N)*F1+H90(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 1990 - 1995:
C
100   F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-1990)/5.D0
      F1=1.D0-F2
      DO 105 N=1,105
         G(N)=G90(N)*F1+G95(N)*F2
105      H(N)=H90(N)*F1+H95(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 1995 - 2000:
C
110   F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-1995)/5.D0
      F1=1.D0-F2
      DO 115 N=1,105   !  THE 2000 COEFFICIENTS (G00) GO THROUGH THE ORDER 13, NOT 10
         G(N)=G95(N)*F1+G00(N)*F2
115      H(N)=H95(N)*F1+H00(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 2000 - 2005:
C
120   F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25D0-2000)/5.D0
      F1=1.D0-F2
      DO 125 N=1,105
         G(N)=G00(N)*F1+G05(N)*F2
125      H(N)=H00(N)*F1+H05(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 2005 - 2010:
C
130   F2=(FLOAT(IY)+FLOAT(IDAY-1)/365.25-2005)/5.
      F1=1.-F2
      DO 135 N=1,105
         G(N)=G05(N)*F1+G10(N)*F2
135      H(N)=H05(N)*F1+H10(N)*F2
      GOTO 300
C
C       INTERPOLATE BETWEEN 2010 - 2015:
C
140   F2=(DFLOAT(IY)+DFLOAT(IDAY-1)/365.25-2010)/5.
      F1=1.-F2
      DO 145 N=1,105
         G(N)=G10(N)*F1+G15(N)*F2
145      H(N)=H10(N)*F1+H15(N)*F2
      GOTO 300
C
C   COEFFICIENTS FOR A GIVEN YEAR HAVE BEEN CALCULATED; NOW MULTIPLY
C   THEM BY SCHMIDT NORMALIZATION FACTORS:
C
300   S=1.D0
      DO 150 N=2,14
         MN=N*(N-1)/2+1
         S=S*DFLOAT(2*N-3)/DFLOAT(N-1)
         G(MN)=G(MN)*S
         H(MN)=H(MN)*S
         P=S
         DO 150 M=2,N
            AA=1.D0
            IF (M.EQ.2) AA=2.D0
            P=P*DSQRT(AA*DFLOAT(N-M+1)/DFLOAT(N+M-2))
            MNN=MN+M-1
            G(MNN)=G(MNN)*P
150         H(MNN)=H(MNN)*P

           G_10=-G(2)
           G_11= G(3)
           H_11= H(3)
C
C  NOW CALCULATE GEO COMPONENTS OF THE UNIT VECTOR EzMAG, PARALLEL TO GEODIPOLE AXIS:
C   SIN(TETA0)*COS(LAMBDA0), SIN(TETA0)*SIN(LAMBDA0), AND COS(TETA0)
C         ST0 * CL0                ST0 * SL0                CT0
C
      SQ=G_11**2+H_11**2
      SQQ=DSQRT(SQ)
      SQR=DSQRT(G_10**2+SQ)
      SL0=-H_11/SQQ
      CL0=-G_11/SQQ
      ST0=SQQ/SQR
      CT0=G_10/SQR
      STCL=ST0*CL0
      STSL=ST0*SL0
      CTSL=CT0*SL0
      CTCL=CT0*CL0
C
C  NOW CALCULATE GEI COMPONENTS (S1,S2,S3) OF THE UNIT VECTOR S = EX_GSE
C    POINTING FROM THE EARTH'S CENTER TO SUN
C
      CALL SUN_08 (IY,IDAY,IHOUR,MIN,ISEC,GST,SLONG,SRASN,SDEC)
C
      S1=DCOS(SRASN)*DCOS(SDEC)
      S2=DSIN(SRASN)*DCOS(SDEC)
      S3=DSIN(SDEC)
C
C  NOW CALCULATE GEI COMPONENTS (DZ1,DZ2,DZ3) OF THE UNIT VECTOR EZGSE
C  POINTING NORTHWARD AND ORTHOGONAL TO THE ECLIPTIC PLANE, AS
C  (0,-SIN(OBLIQ),COS(OBLIQ)). FOR THE EPOCH 1978, OBLIQ = 23.44214 DEGS.
C  HERE WE USE A MORE ACCURATE TIME-DEPENDENT VALUE, DETERMINED AS:
C
      DJ=DFLOAT(365*(IY-1900)+(IY-1901)/4 +IDAY)
     * -0.5+DFLOAT(IHOUR*3600+MIN*60+ISEC)/86400.D0
      T=DJ/36525.D0
      OBLIQ=(23.45229D0-0.0130125D0*T)/57.2957795D0
      DZ1=0.D0
      DZ2=-DSIN(OBLIQ)
      DZ3=DCOS(OBLIQ)
C
C  NOW WE OBTAIN GEI COMPONENTS OF THE UNIT VECTOR EYGSE=(DY1,DY2,DY3),
C  COMPLETING THE RIGHT-HANDED SYSTEM. THEY CAN BE FOUND FROM THE VECTOR
C  PRODUCT EZGSE x EXGSE = (DZ1,DZ2,DZ3) x (S1,S2,S3):
C
      DY1=DZ2*S3-DZ3*S2
      DY2=DZ3*S1-DZ1*S3
      DY3=DZ1*S2-DZ2*S1
C
C  NOW LET'S CALCULATE GEI COMPONENTS OF THE UNIT VECTOR X = EXGSW, DIRECTED ANTIPARALLEL
C  TO THE OBSERVED SOLAR WIND FLOW. FIRST, CALCULATE ITS COMPONENTS IN GSE:
C
      V=DSQRT(VGSEX**2+VGSEY**2+VGSEZ**2)
      DX1=-VGSEX/V
      DX2=-VGSEY/V
      DX3=-VGSEZ/V
C
C  THEN IN GEI:
C
      X1=DX1*S1+DX2*DY1+DX3*DZ1
      X2=DX1*S2+DX2*DY2+DX3*DZ2
      X3=DX1*S3+DX2*DY3+DX3*DZ3
C
C  NOW CALCULATE GEI COMPONENTS (DIP1,DIP2,DIP3) OF THE UNIT VECTOR DIP = EZ_SM = EZ_MAG,
C   ALIGNED WITH THE GEODIPOLE AND POINTING NORTHWARD FROM ECLIPTIC PLANE:
C
      CGST=DCOS(GST)
      SGST=DSIN(GST)
C
      DIP1=STCL*CGST-STSL*SGST
      DIP2=STCL*SGST+STSL*CGST
      DIP3=CT0
C
C  THIS ALLOWS US TO CALCULATE GEI COMPONENTS OF THE UNIT VECTOR Y = EYGSW
C   BY TAKING THE VECTOR PRODUCT DIP x X AND NORMALIZING IT TO UNIT LENGTH:
C
      Y1=DIP2*X3-DIP3*X2
      Y2=DIP3*X1-DIP1*X3
      Y3=DIP1*X2-DIP2*X1
      Y=DSQRT(Y1*Y1+Y2*Y2+Y3*Y3)
      Y1=Y1/Y
      Y2=Y2/Y
      Y3=Y3/Y
C
C   AND GEI COMPONENTS OF THE UNIT VECTOR Z = EZGSW = EXGSW x EYGSW = X x Y:
C
      Z1=X2*Y3-X3*Y2
      Z2=X3*Y1-X1*Y3
      Z3=X1*Y2-X2*Y1
C
C   ELEMENTS OF THE MATRIX GSE TO GSW ARE THE SCALAR PRODUCTS:
C
C  E11=(EXGSE,EXGSW)  E12=(EXGSE,EYGSW)  E13=(EXGSE,EZGSW)
C  E21=(EYGSE,EXGSW)  E22=(EYGSE,EYGSW)  E23=(EYGSE,EZGSW)
C  E31=(EZGSE,EXGSW)  E32=(EZGSE,EYGSW)  E33=(EZGSE,EZGSW)
C
      E11= S1*X1 +S2*X2 +S3*X3
      E12= S1*Y1 +S2*Y2 +S3*Y3
      E13= S1*Z1 +S2*Z2 +S3*Z3
      E21=DY1*X1+DY2*X2+DY3*X3
      E22=DY1*Y1+DY2*Y2+DY3*Y3
      E23=DY1*Z1+DY2*Z2+DY3*Z3
      E31=DZ1*X1+DZ2*X2+DZ3*X3
      E32=DZ1*Y1+DZ2*Y2+DZ3*Y3
      E33=DZ1*Z1+DZ2*Z2+DZ3*Z3
C
C   GEODIPOLE TILT ANGLE IN THE GSW SYSTEM: PSI=ARCSIN(DIP,EXGSW)
C
      SPS=DIP1*X1+DIP2*X2+DIP3*X3
      CPS=DSQRT(1.D0-SPS**2)
      PSI=DASIN(SPS)
C
C   ELEMENTS OF THE MATRIX GEO TO GSW ARE THE SCALAR PRODUCTS:
C
C   A11=(EXGEO,EXGSW), A12=(EYGEO,EXGSW), A13=(EZGEO,EXGSW),
C   A21=(EXGEO,EYGSW), A22=(EYGEO,EYGSW), A23=(EZGEO,EYGSW),
C   A31=(EXGEO,EZGSW), A32=(EYGEO,EZGSW), A33=(EZGEO,EZGSW),
C
C   ALL THE UNIT VECTORS IN BRACKETS ARE ALREADY DEFINED IN GEI:
C
C  EXGEO=(CGST,SGST,0), EYGEO=(-SGST,CGST,0), EZGEO=(0,0,1)
C  EXGSW=(X1,X2,X3),  EYGSW=(Y1,Y2,Y3),   EZGSW=(Z1,Z2,Z3)
C                                                           AND  THEREFORE:
C
      A11=X1*CGST+X2*SGST
      A12=-X1*SGST+X2*CGST
      A13=X3
      A21=Y1*CGST+Y2*SGST
      A22=-Y1*SGST+Y2*CGST
      A23=Y3
      A31=Z1*CGST+Z2*SGST
      A32=-Z1*SGST+Z2*CGST
      A33=Z3
C
C  NOW CALCULATE ELEMENTS OF THE MATRIX MAG TO SM (ONE ROTATION ABOUT THE GEODIPOLE AXIS);
C   THEY ARE FOUND AS THE SCALAR PRODUCTS: CFI=GM22=(EYSM,EYMAG)=(EYGSW,EYMAG),
C                                          SFI=GM23=(EYSM,EXMAG)=(EYGSW,EXMAG),
C    DERIVED AS FOLLOWS:
C
C IN GEO, THE VECTORS EXMAG AND EYMAG HAVE THE COMPONENTS (CT0*CL0,CT0*SL0,-ST0)
C  AND (-SL0,CL0,0), RESPECTIVELY.    HENCE, IN GEI THEIR COMPONENTS ARE:
C  EXMAG:    CT0*CL0*COS(GST)-CT0*SL0*SIN(GST)
C            CT0*CL0*SIN(GST)+CT0*SL0*COS(GST)
C            -ST0
C  EYMAG:    -SL0*COS(GST)-CL0*SIN(GST)
C            -SL0*SIN(GST)+CL0*COS(GST)
C             0
C  NOW, NOTE THAT GEI COMPONENTS OF EYSM=EYGSW WERE FOUND ABOVE AS Y1, Y2, AND Y3,
C  AND WE ONLY HAVE TO COMBINE THESE QUANTITIES INTO SCALAR PRODUCTS:
C
      EXMAGX=CT0*(CL0*CGST-SL0*SGST)
      EXMAGY=CT0*(CL0*SGST+SL0*CGST)
      EXMAGZ=-ST0
      EYMAGX=-(SL0*CGST+CL0*SGST)
      EYMAGY=-(SL0*SGST-CL0*CGST)
      CFI=Y1*EYMAGX+Y2*EYMAGY
      SFI=Y1*EXMAGX+Y2*EXMAGY+Y3*EXMAGZ
C
 10   FORMAT(//1X,
     *'**** RECALC_08 WARNS: YEAR IS OUT OF INTERVAL 1965-2020: IYEAR=',
     *I4,/,6X,'CALCULATIONS WILL BE DONE FOR IYEAR=',I4,/)
      RETURN
      END
c
c==================================================================================

      SUBROUTINE GSWGSE_08 (XGSW,YGSW,ZGSW,XGSE,YGSE,ZGSE,J)
C
C  THIS SUBROUTINE TRANSFORMS COMPONENTS OF ANY VECTOR BETWEEN THE STANDARD GSE
C  COORDINATE SYSTEM AND THE GEOCENTRIC SOLAR-WIND (GSW, aka GSWM), DEFINED AS FOLLOWS
C  (HONES ET AL., PLANET.SPACE SCI., V.34, P.889, 1986; TSYGANENKO ET AL., JGRA,
C  V.103(A4), P.6827, 1998):
C
C  IN THE GSW SYSTEM, X AXIS IS ANTIPARALLEL TO THE OBSERVED DIRECTION OF THE SOLAR WIND FLOW.
C  TWO OTHER AXES, Y AND Z, ARE DEFINED IN THE SAME WAY AS FOR THE STANDARD GSM, THAT IS,
C  Z AXIS ORTHOGONAL TO X AXIS, POINTS NORTHWARD, AND LIES IN THE PLANE DEFINED BY THE X-
C  AND GEODIPOLE AXIS. THE Y AXIS COMPLETES THE RIGHT-HANDED SYSTEM.
C
C  THE GSW SYSTEM BECOMES IDENTICAL TO THE STANDARD GSM IN THE CASE OF
C   A STRICTLY RADIAL SOLAR WIND FLOW.
C
C  AUTHOR:  N. A. TSYGANENKO
C  ADDED TO 2008 VERSION OF GEOPACK: JAN 27, 2008.
C  LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C                    J>0                       J<0
C-----INPUT:   J,XGSW,YGSW,ZGSW          J,XGSE,YGSE,ZGSE
C-----OUTPUT:    XGSE,YGSE,ZGSE            XGSW,YGSW,ZGSW
C
C  IMPORTANT THINGS TO REMEMBER:
C
C   (1) BEFORE CALLING GSWGSE_08, BE SURE TO INVOKE SUBROUTINE RECALC_08, IN ORDER
C       TO DEFINE ALL NECESSARY ELEMENTS OF TRANSFORMATION MATRICES
C
C   (2) IN THE ABSENCE OF INFORMATION ON THE SOLAR WIND DIRECTION, E.G., WITH ONLY SCALAR
C       SPEED V KNOWN, THIS SUBROUTINE CAN BE USED TO CONVERT VECTORS TO ABERRATED
C       COORDINATE SYSTEM, TAKING INTO ACCOUNT EARTH'S ORBITAL SPEED OF 29 KM/S.
C       TO DO THAT, SPECIFY THE LAST 3 PARAMETERS IN RECALC_08 AS FOLLOWS:
C       VGSEX=-V, VGSEY=29.0, VGSEZ=0.0.
C
C       IT SHOULD ALSO BE KEPT IN MIND THAT IN SOME SOLAR WIND DATABASES THE ABERRATION
C       EFFECT HAS ALREADY BEEN TAKEN INTO ACCOUNT BY SUBTRACTING 29 KM/S FROM VYGSE;
C       IN THAT CASE, THE ORIGINAL VYGSE VALUES SHOULD BE RESTORED BY ADDING BACK THE
C       29 KM/S CORRECTION. WHETHER OR NOT TO DO THAT, MUST BE VERIFIED WITH THE DATA
C       ORIGINATOR (OR CAN BE DETERMINED BY CALCULATING THE AVERAGE VGSEY OVER
C       A SUFFICIENTLY LONG TIME INTERVAL)
C
C   (3) IF NO INFORMATION IS AVAILABLE ON THE SOLAR WIND SPEED, THEN SET VGSEX=-400.0
c       AND  VGSEY=VGSEZ=0. IN THAT CASE, THE GSW COORDINATE SYSTEM BECOMES
c       IDENTICAL TO THE STANDARD ONE.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /GEOPACK1/ AAA(25),E11,E21,E31,E12,E22,E32,E13,E23,E33
C
C  DIRECT TRANSFORMATION:
C
      IF (J.GT.0) THEN
        XGSE=XGSW*E11+YGSW*E12+ZGSW*E13
        YGSE=XGSW*E21+YGSW*E22+ZGSW*E23
        ZGSE=XGSW*E31+YGSW*E32+ZGSW*E33
      ENDIF
C
C   INVERSE TRANSFORMATION: CARRIED OUT USING THE TRANSPOSED MATRIX:
C
      IF (J.LT.0) THEN
        XGSW=XGSE*E11+YGSE*E21+ZGSE*E31
        YGSW=XGSE*E12+YGSE*E22+ZGSE*E32
        ZGSW=XGSE*E13+YGSE*E23+ZGSE*E33
      ENDIF
C
      RETURN
      END
C
C========================================================================================
C
      SUBROUTINE GEOMAG_08 (XGEO,YGEO,ZGEO,XMAG,YMAG,ZMAG,J)
C
C    CONVERTS GEOGRAPHIC (GEO) TO DIPOLE (MAG) COORDINATES OR VICE VERSA.
C
C                    J>0                       J<0
C-----INPUT:  J,XGEO,YGEO,ZGEO           J,XMAG,YMAG,ZMAG
C-----OUTPUT:    XMAG,YMAG,ZMAG           XGEO,YGEO,ZGEO
C
C  ATTENTION:  SUBROUTINE  RECALC_08  MUST BE INVOKED BEFORE GEOMAG_08 IN TWO CASES:
C     /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES
C     /B/  IF THE VALUES OF IYEAR AND/OR IDAY HAVE BEEN CHANGED
C
C  NO INFORMATION IS REQUIRED HERE ON THE SOLAR WIND VELOCITY, SO ONE
C  CAN SET VGSEX=-400.0, VGSEY=0.0, VGSEZ=0.0 IN RECALC_08.
C
C   LAST MOFIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C
C   AUTHOR:  N. A. TSYGANENKO
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /GEOPACK1/ ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,AB(26)

      IF(J.GT.0) THEN
       XMAG=XGEO*CTCL+YGEO*CTSL-ZGEO*ST0
       YMAG=YGEO*CL0-XGEO*SL0
       ZMAG=XGEO*STCL+YGEO*STSL+ZGEO*CT0
      ELSE
       XGEO=XMAG*CTCL-YMAG*SL0+ZMAG*STCL
       YGEO=XMAG*CTSL+YMAG*CL0+ZMAG*STSL
       ZGEO=ZMAG*CT0-XMAG*ST0
      ENDIF

      RETURN
      END
c
c=========================================================================================
c
      SUBROUTINE GEIGEO_08 (XGEI,YGEI,ZGEI,XGEO,YGEO,ZGEO,J)
C
C   CONVERTS EQUATORIAL INERTIAL (GEI) TO GEOGRAPHICAL (GEO) COORDS
C   OR VICE VERSA.
C                    J>0                J<0
C----INPUT:  J,XGEI,YGEI,ZGEI    J,XGEO,YGEO,ZGEO
C----OUTPUT:   XGEO,YGEO,ZGEO      XGEI,YGEI,ZGEI
C
C  ATTENTION:  SUBROUTINE  RECALC_08  MUST BE INVOKED BEFORE GEIGEO_08 IN TWO CASES:
C     /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES
C     /B/  IF THE CURRENT VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC HAVE BEEN CHANGED
C
C  NO INFORMATION IS REQUIRED HERE ON THE SOLAR WIND VELOCITY, SO ONE
C  CAN SET VGSEX=-400.0, VGSEY=0.0, VGSEZ=0.0 IN RECALC_08.
C
C     LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)

C     AUTHOR:  N. A. TSYGANENKO
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /GEOPACK1/ A(13),CGST,SGST,B(19)
C
      IF(J.GT.0) THEN
       XGEO=XGEI*CGST+YGEI*SGST
       YGEO=YGEI*CGST-XGEI*SGST
       ZGEO=ZGEI
      ELSE
       XGEI=XGEO*CGST-YGEO*SGST
       YGEI=YGEO*CGST+XGEO*SGST
       ZGEI=ZGEO
      ENDIF

      RETURN
      END
C
C=======================================================================================
C
      SUBROUTINE MAGSM_08 (XMAG,YMAG,ZMAG,XSM,YSM,ZSM,J)
C
C  CONVERTS DIPOLE (MAG) TO SOLAR MAGNETIC (SM) COORDINATES OR VICE VERSA
C
C                    J>0              J<0
C----INPUT: J,XMAG,YMAG,ZMAG     J,XSM,YSM,ZSM
C----OUTPUT:    XSM,YSM,ZSM       XMAG,YMAG,ZMAG
C
C  ATTENTION:  SUBROUTINE  RECALC_08  MUST BE INVOKED BEFORE MAGSM_08 IN THREE CASES:
C     /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES, OR
C     /B/  IF THE VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC HAVE CHANGED, AND/OR
C     /C/  IF THE VALUES OF COMPONENTS OF THE SOLAR WIND FLOW VELOCITY HAVE CHANGED
C
C    IMPORTANT NOTE:
C
C        A NON-STANDARD DEFINITION IS IMPLIED HERE FOR THE SOLAR MAGNETIC COORDINATE
C        SYSTEM:  IT IS ASSUMED THAT THE XSM AXIS LIES IN THE PLANE DEFINED BY THE
C        GEODIPOLE AXIS AND THE OBSERVED VECTOR OF THE SOLAR WIND FLOW (RATHER THAN
C        THE EARTH-SUN LINE).  IN ORDER TO CONVERT MAG COORDINATES TO AND FROM THE
C        STANDARD SM COORDINATES, INVOKE RECALC_08 WITH VGSEX=-400.0, VGSEY=0.0, VGSEZ=0.0
C
C     LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C
C     AUTHOR:  N. A. TSYGANENKO
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /GEOPACK1/ A(8),SFI,CFI,B(24)
C
      IF (J.GT.0) THEN
       XSM=XMAG*CFI-YMAG*SFI
       YSM=XMAG*SFI+YMAG*CFI
       ZSM=ZMAG
      ELSE
       XMAG=XSM*CFI+YSM*SFI
       YMAG=YSM*CFI-XSM*SFI
       ZMAG=ZSM
      ENDIF

      RETURN
      END
C
C=====================================================================================
C
       SUBROUTINE SMGSW_08 (XSM,YSM,ZSM,XGSW,YGSW,ZGSW,J)
C
C  CONVERTS SOLAR MAGNETIC (SM) TO GEOCENTRIC SOLAR-WIND (GSW) COORDINATES OR VICE VERSA.
C
C                  J>0                 J<0
C-----INPUT: J,XSM,YSM,ZSM        J,XGSW,YGSW,ZGSW
C----OUTPUT:  XGSW,YGSW,ZGSW       XSM,YSM,ZSM
C
C  ATTENTION:  SUBROUTINE RECALC_08 MUST BE INVOKED BEFORE SMGSW_08 IN THREE CASES:
C     /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES
C     /B/  IF THE VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC HAVE BEEN CHANGED
C     /C/  IF THE VALUES OF COMPONENTS OF THE SOLAR WIND FLOW VELOCITY HAVE CHANGED
C
C    IMPORTANT NOTE:
C
C        A NON-STANDARD DEFINITION IS IMPLIED HERE FOR THE SOLAR MAGNETIC (SM) COORDINATE
C        SYSTEM:  IT IS ASSUMED THAT THE XSM AXIS LIES IN THE PLANE DEFINED BY THE
C        GEODIPOLE AXIS AND THE OBSERVED VECTOR OF THE SOLAR WIND FLOW (RATHER THAN
C        THE EARTH-SUN LINE).  IN ORDER TO CONVERT MAG COORDINATES TO AND FROM THE
C        STANDARD SM COORDINATES, INVOKE RECALC_08 WITH VGSEX=-400.0, VGSEY=0.0, VGSEZ=0.0
C
C     LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C
C     AUTHOR:  N. A. TSYGANENKO
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /GEOPACK1/ A(10),SPS,CPS,B(22)

      IF (J.GT.0) THEN
       XGSW=XSM*CPS+ZSM*SPS
       YGSW=YSM
       ZGSW=ZSM*CPS-XSM*SPS
      ELSE
       XSM=XGSW*CPS-ZGSW*SPS
       YSM=YGSW
       ZSM=XGSW*SPS+ZGSW*CPS
      ENDIF

      RETURN
      END
C
C==========================================================================================
C
      SUBROUTINE GEOGSW_08 (XGEO,YGEO,ZGEO,XGSW,YGSW,ZGSW,J)
C
C CONVERTS GEOGRAPHIC (GEO) TO GEOCENTRIC SOLAR-WIND (GSW) COORDINATES OR VICE VERSA.
C
C                   J>0                   J<0
C----- INPUT:  J,XGEO,YGEO,ZGEO    J,XGSW,YGSW,ZGSW
C---- OUTPUT:    XGSW,YGSW,ZGSW      XGEO,YGEO,ZGEO
C
C  ATTENTION:  SUBROUTINE  RECALC_08  MUST BE INVOKED BEFORE GEOGSW_08 IN THREE CASES:
C     /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES, OR
C     /B/  IF THE VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC  HAVE CHANGED, AND/OR
C     /C/  IF THE VALUES OF COMPONENTS OF THE SOLAR WIND FLOW VELOCITY HAVE CHANGED
C
C  NOTE: THIS SUBROUTINE CONVERTS GEO VECTORS TO AND FROM THE SOLAR-WIND GSW COORDINATE
C        SYSTEM, TAKING INTO ACCOUNT POSSIBLE DEFLECTIONS OF THE SOLAR WIND DIRECTION FROM
C        STRICTLY RADIAL.  BEFORE CONVERTING TO/FROM STANDARD GSM COORDINATES, INVOKE RECALC_08
C        WITH VGSEX=-400.0 and VGSEY=0.0, VGSEZ=0.0
C
C     LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C
C     AUTHOR:  N. A. TSYGANENKO
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /GEOPACK1/ AA(16),A11,A21,A31,A12,A22,A32,A13,A23,A33,B(9)
C
      IF (J.GT.0) THEN
       XGSW=A11*XGEO+A12*YGEO+A13*ZGEO
       YGSW=A21*XGEO+A22*YGEO+A23*ZGEO
       ZGSW=A31*XGEO+A32*YGEO+A33*ZGEO
      ELSE
       XGEO=A11*XGSW+A21*YGSW+A31*ZGSW
       YGEO=A12*XGSW+A22*YGSW+A32*ZGSW
       ZGEO=A13*XGSW+A23*YGSW+A33*ZGSW
      ENDIF

      RETURN
      END
C
C=====================================================================================
C
      SUBROUTINE GEODGEO_08 (H,XMU,R,THETA,J)
C
C  THIS SUBROUTINE (1) CONVERTS VERTICAL LOCAL HEIGHT (ALTITUDE) H AND GEODETIC
C  LATITUDE XMU INTO GEOCENTRIC COORDINATES R AND THETA (GEOCENTRIC RADIAL
C  DISTANCE AND COLATITUDE, RESPECTIVELY; ALSO KNOWN AS ECEF COORDINATES),
C  AS WELL AS (2) PERFORMS THE INVERSE TRANSFORMATION FROM {R,THETA} TO {H,XMU}.
C
C  THE SUBROUTINE USES WORLD GEODETIC SYSTEM WGS84 PARAMETERS FOR THE EARTH'S
C  ELLIPSOID. THE ANGULAR QUANTITIES (GEO COLATITUDE THETA AND GEODETIC LATITUDE
C  XMU) ARE IN RADIANS, AND THE DISTANCES (GEOCENTRIC RADIUS R AND ALTITUDE H
C  ABOVE THE EARTH'S ELLIPSOID) ARE IN KILOMETERS.
C
C  IF J>0, THE TRANSFORMATION IS MADE FROM GEODETIC TO GEOCENTRIC COORDINATES
C   USING SIMPLE DIRECT EQUATIONS.
C  IF J<0, THE INVERSE TRANSFORMATION FROM GEOCENTRIC TO GEODETIC COORDINATES
C   IS MADE BY MEANS OF A FAST ITERATIVE ALGORITHM.
C
c-------------------------------------------------------------------------------
C                   J>0                     |            J<0
c-------------------------------------------|-----------------------------------
C--INPUT:   J        H          XMU         |    J         R          THETA
c         flag  altitude (km)  geodetic     |   flag   geocentric    spherical
c                              latitude     |         distance (km) colatitude
c                              (radians)    |                        (radians)
c-------------------------------------------|-----------------------------------
c                                           |
C----OUTPUT:         R           THETA      |          H              XMU
C                geocentric    spherical    |      altitude (km)    geodetic
C                distance (km) colatitude   |                       latitude
C                              (radians)    |                       (radians)
C-------------------------------------------------------------------------------
C
C   AUTHOR:  N. A. TSYGANENKO
c   DATE:    DEC 5, 2007
C   LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA R_EQ, BETA /6378.137D0, 6.73949674228D-3/
c
c  R_EQ is the semi-major axis of the Earth's ellipsoid, and BETA is its
c  second eccentricity squared
c
      DATA TOL /1.D-6/
c
c   Direct transformation (GEOD=>GEO):
c
      IF (J.GT.0) THEN
       COSXMU=DCOS(XMU)
       SINXMU=DSIN(XMU)
       DEN=DSQRT(COSXMU**2+(SINXMU/(1.0D0+BETA))**2)
       COSLAM=COSXMU/DEN
       SINLAM=SINXMU/(DEN*(1.0D0+BETA))
       RS=R_EQ/DSQRT(1.0D0+BETA*SINLAM**2)
       X=RS*COSLAM+H*COSXMU
       Z=RS*SINLAM+H*SINXMU
       R=DSQRT(X**2+Z**2)
       THETA=DACOS(Z/R)
      ENDIF

c
c   Inverse transformation (GEO=>GEOD):
c
      IF (J.LT.0) THEN
       N=0
       PHI=1.570796327D0-THETA
       PHI1=PHI
  1    SP=DSIN(PHI1)
       ARG=SP*(1.0D0+BETA)/DSQRT(1.0D0+BETA*(2.0D0+BETA)*SP**2)
       XMUS=DASIN(ARG)
       RS=R_EQ/DSQRT(1.0D0+BETA*DSIN(PHI1)**2)
       COSFIMS=DCOS(PHI1-XMUS)
       H=DSQRT((RS*COSFIMS)**2+R**2-RS**2)-RS*COSFIMS
       Z=RS*DSIN(PHI1)+H*DSIN(XMUS)
       X=RS*DCOS(PHI1)+H*DCOS(XMUS)
       RR=DSQRT(X**2+Z**2)
       DPHI=DASIN(Z/RR)-PHI
       PHI1=PHI1-DPHI
       N=N+1
       IF (DABS(DPHI).GT.TOL.AND.N.LT.100) GOTO 1
       XMU=XMUS
      ENDIF

      RETURN
      END
C
C=====================================================================================
C
      SUBROUTINE RHAND_08 (X,Y,Z,R1,R2,R3,IOPT,PARMOD,EXNAME,INNAME)
C
C  CALCULATES THE COMPONENTS OF THE RIGHT HAND SIDE VECTOR IN THE GEOMAGNETIC FIELD
C    LINE EQUATION  (a subsidiary subroutine for the subroutine STEP_08)
C
C     LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C
C     AUTHOR:  N. A. TSYGANENKO
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PARMOD(10)
C
C     EXNAME AND INNAME ARE NAMES OF SUBROUTINES FOR THE EXTERNAL AND INTERNAL
C     PARTS OF THE TOTAL FIELD, E.G., T96_01 AND IGRF_GSW_08
C
      COMMON /GEOPACK1/ A(12),DS3,BB(2),PSI,CC(18)

      CALL EXNAME (IOPT,PARMOD,PSI,X,Y,Z,BXGSW,BYGSW,BZGSW)
      CALL INNAME (X,Y,Z,HXGSW,HYGSW,HZGSW)

      BX=BXGSW+HXGSW
      BY=BYGSW+HYGSW
      BZ=BZGSW+HZGSW
      B=DS3/DSQRT(BX**2+BY**2+BZ**2)
      R1=BX*B
      R2=BY*B
      R3=BZ*B
      RETURN
      END
C
C===================================================================================
C
      SUBROUTINE STEP_08(X,Y,Z,DS,DSMAX,ERRIN,IOPT,PARMOD,EXNAME,INNAME)
C
C   RE-CALCULATES THE INPUT VALUES {X,Y,Z} (IN GSW COORDINATES) FOR ANY POINT ON A FIELD LINE,
C     BY MAKING A STEP ALONG THAT LINE USING RUNGE-KUTTA-MERSON ALGORITHM (G.N. Lance, Numerical
C      methods for high-speed computers, Iliffe & Sons, London 1960.)
C   DS IS A PRESCRIBED VALUE OF THE CURRENT STEP SIZE, DSMAX IS ITS UPPER LIMIT.
C   ERRIN IS A PERMISSIBLE ERROR (ITS OPTIMAL VALUE SPECIFIED IN THE S/R TRACE_08)
C     IF THE ACTUAL ERROR (ERRCUR) AT THE CURRENT STEP IS LARGER THAN ERRIN, THE STEP IS REJECTED,
C       AND THE CALCULATION IS REPEATED ANEW WITH HALVED STEPSIZE DS.
C     IF ERRCUR IS SMALLER THAN ERRIN, THE STEP IS ACCEPTED, AND THE CURRENT VALUE OF DS IS RETAINED
C       FOR THE NEXT STEP.
C     IF ERRCUR IS SMALLER THAN 0.04*ERRIN, THE STEP IS ACCEPTED, AND THE VALUE OF DS FOR THE NEXT STEP
C       IS INCREASED BY THE FACTOR 1.5, BUT NOT LARGER THAN DSMAX.
C   IOPT IS A FLAG, RESERVED FOR SPECIFYNG A VERSION OF THE EXTERNAL FIELD MODEL EXNAME.
C   ARRAY PARMOD(10) CONTAINS INPUT PARAMETERS FOR THE MODEL EXNAME.
C   EXNAME IS THE NAME OF THE SUBROUTINE FOR THE EXTERNAL FIELD MODEL.
C   INNAME IS THE NAME OF THE SUBROUTINE FOR THE INTERNAL FIELD MODEL (EITHER DIP_08 OR IGRF_GSW_08)
C
C   ALL THE ABOVE PARAMETERS ARE INPUT ONES; OUTPUT IS THE RECALCULATED VALUES OF X,Y,Z
C
C     LAST MODIFICATION:    APRIL 21, 2008 (SEE ERRATA AS OF THAT DATE)
C
C     AUTHOR:  N. A. TSYGANENKO
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PARMOD(10)
      COMMON /GEOPACK1/ A(12),DS3,B(21)
      EXTERNAL EXNAME,INNAME

  1   DS3=-DS/3.D0
      CALL RHAND_08 (X,Y,Z,R11,R12,R13,IOPT,PARMOD,EXNAME,INNAME)
      CALL RHAND_08 (X+R11,Y+R12,Z+R13,R21,R22,R23,IOPT,PARMOD,EXNAME,
     * INNAME)
      CALL RHAND_08 (X+.5D0*(R11+R21),Y+.5D0*(R12+R22),Z+.5D0*
     *(R13+R23),R31,R32,R33,IOPT,PARMOD,EXNAME,INNAME)
      CALL RHAND_08 (X+.375D0*(R11+3.D0*R31),Y+.375D0*(R12+3.D0*R32
     *),Z+.375D0*(R13+3.D0*R33),R41,R42,R43,IOPT,PARMOD,EXNAME,INNAME)
      CALL RHAND_08 (X+1.5D0*(R11-3.D0*R31+4.D0*R41),Y+1.5D0*(R12-
     *3.D0*R32+4.D0*R42),Z+1.5D0*(R13-3.D0*R33+4.D0*R43),
     *R51,R52,R53,IOPT,PARMOD,EXNAME,INNAME)
      ERRCUR=DABS(R11-4.5D0*R31+4.D0*R41-.5D0*R51)+DABS(R12-4.5D0*R32
     *+4.D0*R42-.5D0*R52)+DABS(R13-4.5D0*R33+4.D0*R43-.5D0*R53)
C
C  READY FOR MAKING THE STEP, BUT CHECK THE ACCURACY; IF INSUFFICIENT,
C   REPEAT THE STEP WITH HALVED STEPSIZE:
C
      IF (ERRCUR.GT.ERRIN) THEN
      DS=DS*.5D0
      GOTO 1
      ENDIF
C
C  ACCURACY IS ACCEPTABLE, BUT CHECK IF THE STEPSIZE IS NOT TOO LARGE;
C    OTHERWISE REPEAT THE STEP WITH DS=DSMAX
C
      IF (DABS(DS).GT.DSMAX) THEN
      DS=DSIGN(DSMAX,DS)
      GOTO 1
      ENDIF
C
C  MAKING THE STEP:
C
  2   X=X+.5D0*(R11+4.D0*R41+R51)
      Y=Y+.5D0*(R12+4.D0*R42+R52)
      Z=Z+.5D0*(R13+4.D0*R43+R53)
C
C  IF THE ACTUAL ERROR IS TOO SMALL (LESS THAN 4% OF ERRIN) AND DS SMALLER
C   THAN DSMAX/1.5, THEN WE INCREASE THE STEPSIZE FOR THE NEXT STEP BY 50%
C
      IF(ERRCUR.LT.ERRIN*.04D0.AND.DS.LT.DSMAX/1.5D0) DS=DS*1.5D0
      RETURN
      END
C
C==============================================================================
C
      SUBROUTINE TRACE_08 (XI,YI,ZI,DIR,DSMAX,ERR,RLIM,R0,IOPT,PARMOD,
     * EXNAME,INNAME,XF,YF,ZF,XX,YY,ZZ,L,LMAX)
C
C  TRACES A FIELD LINE FROM AN ARBITRARY POINT OF SPACE TO THE EARTH'S
C  SURFACE OR TO A MODEL LIMITING BOUNDARY.
C
C  THIS SUBROUTINE ALLOWS TWO OPTIONS:
C
C  (1) IF INNAME=IGRF_GSW_08, THEN THE IGRF MODEL WILL BE USED FOR CALCULATING
C      CONTRIBUTION FROM EARTH'S INTERNAL SOURCES. IN THIS CASE, SUBROUTINE
C      RECALC_08 MUST BE CALLED BEFORE USING TRACE_08, WITH PROPERLY SPECIFIED DATE,
C      UNIVERSAL TIME, AND SOLAR WIND VELOCITY COMPONENTS, TO CALCULATE IN ADVANCE
C      ALL QUANTITIES NEEDED FOR THE MAIN FIELD MODEL AND FOR TRANSFORMATIONS
C      BETWEEN INVOLVED COORDINATE SYSTEMS.
C
C  (2) IF INNAME=DIP_08, THEN A PURE DIPOLE FIELD WILL BE USED INSTEAD OF THE IGRF MODEL.
C      IN THIS CASE, THE SUBROUTINE RECALC_08 MUST ALSO BE CALLED BEFORE TRACE_08.
C      HERE ONE CAN CHOOSE EITHER TO
C      (a) CALCULATE DIPOLE TILT ANGLE BASED ON DATE, TIME, AND SOLAR WIND DIRECTION,
C   OR (b) EXPLICITLY SPECIFY THAT ANGLE, WITHOUT ANY REFERENCE TO DATE/UT/SOLAR WIND.
C      IN THE LAST CASE (b), THE SINE (SPS) AND COSINE (CPS) OF THE DIPOLE TILT
C      ANGLE MUST BE SPECIFIED IN ADVANCE (BUT AFTER HAVING CALLED RECALC_08) AND FORWARDED
C      IN THE COMMON BLOCK /GEOPACK1/ (IN ITS 11th AND 12th ELEMENTS, RESPECTIVELY).
C      IN THIS CASE THE ROLE OF THE SUBROUTINE RECALC_08 IS REDUCED TO ONLY CALCULATING
C      THE COMPONENTS OF THE EARTH'S DIPOLE MOMENT.
C
C------------- INPUT PARAMETERS:
C
C   XI,YI,ZI - GSW COORDS OF THE FIELD LINE STARTING POINT (IN EARTH RADII, 1 RE = 6371.2 km),
C
C   DIR - SIGN OF THE TRACING DIRECTION: IF DIR=1.0 THEN THE TRACING IS MADE ANTIPARALLEL
C     TO THE TOTAL FIELD VECTOR (E.G., FROM NORTHERN TO SOUTHERN CONJUGATE POINT);
C     IF DIR=-1.0 THEN THE TRACING PROCEEDS IN THE OPPOSITE DIRECTION, THAT IS, PARALLEL TO
C     THE TOTAL FIELD VECTOR.
C
C   DSMAX - UPPER LIMIT ON THE STEPSIZE (SETS A DESIRED MAXIMAL SPACING BETWEEN
C                 THE FIELD LINE POINTS)
C
C   ERR - PERMISSIBLE STEP ERROR. A REASONABLE ESTIMATE PROVIDING A SUFFICIENT ACCURACY FOR MOST
C         APPLICATIONS IS ERR=0.0001. SMALLER/LARGER VALUES WILL RESULT IN LARGER/SMALLER NUMBER
C         OF STEPS AND, HENCE, OF OUTPUT FIELD LINE POINTS. NOTE THAT USING MUCH SMALLER VALUES
C         OF ERR MAY REQUIRE USING A DOUBLE PRECISION VERSION OF THE ENTIRE PACKAGE.
C
C   R0 -  RADIUS OF A SPHERE (IN RE), DEFINING THE INNER BOUNDARY OF THE TRACING REGION
C         (USUALLY, EARTH'S SURFACE OR THE IONOSPHERE, WHERE R0~1.0)
C         IF THE FIELD LINE REACHES THAT SPHERE FROM OUTSIDE, ITS INBOUND TRACING IS
C         TERMINATED AND THE CROSSING POINT COORDINATES XF,YF,ZF  ARE CALCULATED.
C
C   RLIM - RADIUS OF A SPHERE (IN RE), DEFINING THE OUTER BOUNDARY OF THE TRACING REGION;
C         IF THE FIELD LINE REACHES THAT BOUNDARY FROM INSIDE, ITS OUTBOUND TRACING IS
C         TERMINATED AND THE CROSSING POINT COORDINATES XF,YF,ZF ARE CALCULATED.
C
C   IOPT - A MODEL INDEX; CAN BE USED FOR SPECIFYING A VERSION OF THE EXTERNAL FIELD
C       MODEL (E.G., A NUMBER OF THE KP-INDEX INTERVAL). ALTERNATIVELY, ONE CAN USE THE ARRAY
C       PARMOD FOR THAT PURPOSE (SEE BELOW); IN THAT CASE IOPT IS JUST A DUMMY PARAMETER.
C
C   PARMOD -  A 10-ELEMENT ARRAY CONTAINING INPUT PARAMETERS NEEDED FOR A UNIQUE
C      SPECIFICATION OF THE EXTERNAL FIELD MODEL. THE CONCRETE MEANING OF THE COMPONENTS
C      OF PARMOD DEPENDS ON A SPECIFIC VERSION OF THAT MODEL.
C
C   EXNAME - NAME OF A SUBROUTINE PROVIDING COMPONENTS OF THE EXTERNAL MAGNETIC FIELD
C    (E.G., T89, OR T96_01, ETC.).
C   INNAME - NAME OF A SUBROUTINE PROVIDING COMPONENTS OF THE INTERNAL MAGNETIC FIELD
C    (EITHER DIP_08 OR IGRF_GSW_08).
C
C   LMAX - MAXIMAL LENGTH OF THE ARRAYS XX,YY,ZZ, IN WHICH COORDINATES OF THE FIELD
C          LINE POINTS ARE STORED. LMAX SHOULD BE SET EQUAL TO THE ACTUAL LENGTH OF
C          THE ARRAYS, DEFINED IN THE MAIN PROGRAM AS ACTUAL ARGUMENTS OF THIS SUBROUTINE.
C
C-------------- OUTPUT PARAMETERS:
C
C   XF,YF,ZF - GSW COORDINATES OF THE ENDPOINT OF THE TRACED FIELD LINE.
C   XX,YY,ZZ - ARRAYS OF LENGTH LMAX, CONTAINING COORDINATES OF THE FIELD LINE POINTS.
C   L - ACTUAL NUMBER OF FIELD LINE POINTS, GENERATED BY THIS SUBROUTINE.
C
C ----------------------------------------------------------
C
C     LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C
C     AUTHOR:  N. A. TSYGANENKO
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XX(LMAX),YY(LMAX),ZZ(LMAX), PARMOD(10)
      COMMON /GEOPACK1/ AA(12),DD,BB(21)
      EXTERNAL EXNAME,INNAME
C
      L=0
      NREV=0
      DD=DIR
C
C  INITIALIZE THE STEP SIZE AND STARTING PONT:
C
      DS=0.5D0*DIR
      X=XI
      Y=YI
      Z=ZI
c
c  here we call RHAND_08 just to find out the sign of the radial component of the field
c   vector, and to determine the initial direction of the tracing (i.e., either away
c   or towards Earth):
c
      CALL RHAND_08 (X,Y,Z,R1,R2,R3,IOPT,PARMOD,EXNAME,INNAME)
      AD=0.01D0
      IF (X*R1+Y*R2+Z*R3.LT.0.D0) AD=-0.01D0
C
c     |AD|=0.01 and its sign follows the rule:
c (1) if DIR=1 (tracing antiparallel to B vector) then the sign of AD is the same as of Br
c (2) if DIR=-1 (tracing parallel to B vector) then the sign of AD is opposite to that of Br
c     AD is defined in order to initialize the value of RR (radial distance at previous step):

      RR=DSQRT(X**2+Y**2+Z**2)+AD
c
  1   L=L+1
      IF(L.GT.LMAX) GOTO 7
      XX(L)=X
      YY(L)=Y
      ZZ(L)=Z
      RYZ=Y**2+Z**2
      R2=X**2+RYZ
      R=DSQRT(R2)
C
c  check if the line hit the outer tracing boundary; if yes, then terminate
c   the tracing (label 8). The outer boundary is assumed reached, when the line
c   crosses any of the 3 surfaces: (1) a sphere R=RLIM, (2) a cylinder of radius 40Re,
c   coaxial with the XGSW axis, (3) the plane X=20Re:

      IF (R.GT.RLIM.OR.RYZ.GT.1600.D0.OR.X.GT.20.D0) GOTO 8
c
c  check whether or not the inner tracing boundary was crossed from outside,
c  if yes, then calculate the footpoint position by interpolation (go to label 6):
c
      IF (R.LT.R0.AND.RR.GT.R) GOTO 6

c  check if we are moving outward, or R is still larger than 3Re; if yes, proceed further:
c
      IF (R.GE.RR.OR.R.GE.3.D0) GOTO 4
c
c  now we entered inside the sphere R=3: to avoid too large steps (and hence
c  inaccurate interpolated position of the footpoint), enforce the progressively
c  smaller stepsize values as we approach the inner boundary R=R0:
c
      FC=0.2D0
      IF(R-R0.LT.0.05D0) FC=0.05D0
      AL=FC*(R-R0+0.2D0)
      DS=DIR*AL
c
  4   XR=X
      YR=Y
      ZR=Z
c
      DRP=R-RR
      RR=R
c
      CALL STEP_08 (X,Y,Z,DS,DSMAX,ERR,IOPT,PARMOD,EXNAME,INNAME)
c
C  check the total number NREV of changes in the tracing radial direction; (NREV.GT.2) means
c   that the line started making multiple loops, in which case we stop the process:
C
      R=DSQRT(X**2+Y**2+Z**2)
      DR=R-RR
      IF (DRP*DR.LT.0.D0) NREV=NREV+1
      IF (NREV.GT.4) GOTO 8
C
      GOTO 1
c
c  find the footpoint position by interpolating between the current and previous
c   field line points:
c
  6   R1=(R0-R)/(RR-R)
      X=X-(X-XR)*R1
      Y=Y-(Y-YR)*R1
      Z=Z-(Z-ZR)*R1
      GOTO 8
  7   WRITE (*,10)
      L=LMAX
  8   XF=X
      YF=Y
      ZF=Z
C
C  replace the coordinates of the last (L-th) point in the XX,YY,ZZ arrays
C   so that they correspond to the estimated footpoint position {XF,YF,ZF},
c   satisfying:  sqrt(XF**2+YF**2+ZF**2}=R0
C
      XX(L)=XF
      YY(L)=YF
      ZZ(L)=ZF
C
      RETURN
 10   FORMAT(//,1X,'**** COMPUTATIONS IN THE SUBROUTINE TRACE_08 ARE',
     *' TERMINATED: THE NUMBER OF POINTS EXCEEDED LMAX ****'//)
      END
c
C====================================================================================
C
      SUBROUTINE SHUETAL_MGNP_08(XN_PD,VEL,BZIMF,XGSW,YGSW,ZGSW,
     *  XMGNP,YMGNP,ZMGNP,DIST,ID)
C
C  FOR ANY POINT OF SPACE WITH COORDINATES (XGSW,YGSW,ZGSW) AND SPECIFIED CONDITIONS
C  IN THE INCOMING SOLAR WIND, THIS SUBROUTINE:
C
C (1) DETERMINES IF THE POINT (XGSW,YGSW,ZGSW) LIES INSIDE OR OUTSIDE THE
C      MODEL MAGNETOPAUSE OF SHUE ET AL. (JGR-A, V.103, P. 17691, 1998).
C
C (2) CALCULATES THE GSW POSITION OF A POINT {XMGNP,YMGNP,ZMGNP}, LYING AT THE MODEL
C      MAGNETOPAUSE AND ASYMPTOTICALLY TENDING TO THE NEAREST BOUNDARY POINT WITH
C      RESPECT TO THE OBSERVATION POINT {XGSW,YGSW,ZGSW}, AS IT APPROACHES THE MAGNETO-
C      PAUSE.
C
C  INPUT: XN_PD - EITHER SOLAR WIND PROTON NUMBER DENSITY (PER C.C.) (IF VEL>0)
C                    OR THE SOLAR WIND RAM PRESSURE IN NANOPASCALS   (IF VEL<0)
C         BZIMF - IMF BZ IN NANOTESLAS
C
C         VEL - EITHER SOLAR WIND VELOCITY (KM/SEC)
C                  OR ANY NEGATIVE NUMBER, WHICH INDICATES THAT XN_PD STANDS
C                     FOR THE SOLAR WIND PRESSURE, RATHER THAN FOR THE DENSITY
C
C         XGSW,YGSW,ZGSW - GSW POSITION OF THE OBSERVATION POINT IN EARTH RADII
C
C  OUTPUT: XMGNP,YMGNP,ZMGNP - GSW POSITION OF THE BOUNDARY POINT
C          DIST - DISTANCE (IN RE) BETWEEN THE OBSERVATION POINT (XGSW,YGSW,ZGSW)
C                 AND THE MODEL NAGNETOPAUSE
C          ID -  POSITION FLAG:  ID=+1 (-1) MEANS THAT THE OBSERVATION POINT
C          LIES INSIDE (OUTSIDE) OF THE MODEL MAGNETOPAUSE, RESPECTIVELY.
C
C  OTHER SUBROUTINES USED: T96_MGNP_08
C
c          AUTHOR:  N.A. TSYGANENKO,
C          DATE:    APRIL 4, 2003.
C          LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      IF (VEL.LT.0.) THEN
        P=XN_PD
      ELSE
        P=1.94D-6*XN_PD*VEL**2  ! P IS THE SOLAR WIND DYNAMIC PRESSURE (IN nPa)
      ENDIF

c
c  DEFINE THE ANGLE PHI, MEASURED DUSKWARD FROM THE NOON-MIDNIGHT MERIDIAN PLANE;
C  IF THE OBSERVATION POINT LIES ON THE X AXIS, THE ANGLE PHI CANNOT BE UNIQUELY
C  DEFINED, AND WE SET IT AT ZERO:
c
      IF (YGSW.NE.0.D0.OR.ZGSW.NE.0.D0) THEN
         PHI=DATAN2(YGSW,ZGSW)
      ELSE
         PHI=0.D0
      ENDIF
C
C  FIRST, FIND OUT IF THE OBSERVATION POINT LIES INSIDE THE SHUE ET AL BDRY
C  AND SET THE VALUE OF THE ID FLAG:
C
      ID=-1
      R0=(10.22D0+1.29D0*DTANH(.184D0*(BZIMF+8.14D0)))*P**(-.15151515D0)
      ALPHA=(0.58D0-0.007D0*BZIMF)*(1.D0+0.024*DLOG(P))
      R=DSQRT(XGSW**2+YGSW**2+ZGSW**2)
      RM=R0*(2.D0/(1.D0+XGSW/R))**ALPHA
      IF (R.LE.RM) ID=+1
C
C  NOW, FIND THE CORRESPONDING T96 MAGNETOPAUSE POSITION, TO BE USED AS
C  A STARTING APPROXIMATION IN THE SEARCH OF A CORRESPONDING SHUE ET AL.
C  BOUNDARY POINT:
C
      CALL T96_MGNP_08(P,-1.D0,XGSW,YGSW,ZGSW,XMT96,YMT96,ZMT96,DIST,
     *  ID96)
C
      RHO2=YMT96**2+ZMT96**2
      R=DSQRT(RHO2+XMT96**2)
      ST=DSQRT(RHO2)/R
      CT=XMT96/R
C
C  NOW, USE NEWTON'S ITERATIVE METHOD TO FIND THE NEAREST POINT AT THE
C   SHUE ET AL.'S BOUNDARY:
C
      NIT=0

  1   T=DATAN2(ST,CT)
      RM=R0*(2.D0/(1.D0+CT))**ALPHA

      F=R-RM
      GRADF_R=1.D0
      GRADF_T=-ALPHA/R*RM*ST/(1.D0+CT)
      GRADF=DSQRT(GRADF_R**2+GRADF_T**2)

      DR=-F/GRADF**2
      DT= DR/R*GRADF_T

      R=R+DR
      T=T+DT
      ST=DSIN(T)
      CT=DCOS(T)

      DS=DSQRT(DR**2+(R*DT)**2)

      NIT=NIT+1

      IF (NIT.GT.1000) THEN
         PRINT *,
     *' BOUNDARY POINT COULD NOT BE FOUND; ITERATIONS DO NOT CONVERGE'
      ENDIF

      IF (DS.GT.1.D-4) GOTO 1

      XMGNP=R*DCOS(T)
      RHO=  R*DSIN(T)

      YMGNP=RHO*DSIN(PHI)
      ZMGNP=RHO*DCOS(PHI)

      DIST=DSQRT((XGSW-XMGNP)**2+(YGSW-YMGNP)**2+(ZGSW-ZMGNP)**2)

      RETURN
      END
C
C=======================================================================================
C
      SUBROUTINE T96_MGNP_08(XN_PD,VEL,XGSW,YGSW,ZGSW,XMGNP,YMGNP,ZMGNP,
     * DIST,ID)
C
C  FOR ANY POINT OF SPACE WITH GIVEN COORDINATES (XGSW,YGSW,ZGSW), THIS SUBROUTINE DEFINES
C  THE POSITION OF A POINT (XMGNP,YMGNP,ZMGNP) AT THE T96 MODEL MAGNETOPAUSE WITH THE
C  SAME VALUE OF THE ELLIPSOIDAL TAU-COORDINATE, AND THE DISTANCE BETWEEN THEM.  THIS IS
C  NOT THE SHORTEST DISTANCE D_MIN TO THE BOUNDARY, BUT DIST ASYMPTOTICALLY TENDS TO D_MIN,
C  AS THE OBSERVATION POINT GETS CLOSER TO THE MAGNETOPAUSE.
C
C  INPUT: XN_PD - EITHER SOLAR WIND PROTON NUMBER DENSITY (PER C.C.) (IF VEL>0)
C                    OR THE SOLAR WIND RAM PRESSURE IN NANOPASCALS   (IF VEL<0)
C         VEL - EITHER SOLAR WIND VELOCITY (KM/SEC)
C                  OR ANY NEGATIVE NUMBER, WHICH INDICATES THAT XN_PD STANDS
C                     FOR THE SOLAR WIND PRESSURE, RATHER THAN FOR THE DENSITY
C
C         XGSW,YGSW,ZGSW - COORDINATES OF THE OBSERVATION POINT IN EARTH RADII
C
C  OUTPUT: XMGNP,YMGNP,ZMGNP - GSW POSITION OF THE BOUNDARY POINT, HAVING THE SAME
C          VALUE OF TAU-COORDINATE AS THE OBSERVATION POINT (XGSW,YGSW,ZGSW)
C          DIST -  THE DISTANCE BETWEEN THE TWO POINTS, IN RE,
C          ID -    POSITION FLAG; ID=+1 (-1) MEANS THAT THE POINT (XGSW,YGSW,ZGSW)
C          LIES INSIDE (OUTSIDE) THE MODEL MAGNETOPAUSE, RESPECTIVELY.
C
C  THE PRESSURE-DEPENDENT MAGNETOPAUSE IS THAT USED IN THE T96_01 MODEL
C  (TSYGANENKO, JGR, V.100, P.5599, 1995; ESA SP-389, P.181, OCT. 1996)
C
c   AUTHOR:  N.A. TSYGANENKO
C   DATE:    AUG.1, 1995, REVISED APRIL 3, 2003.
C   LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
C
C  DEFINE SOLAR WIND DYNAMIC PRESSURE (NANOPASCALS, ASSUMING 4% OF ALPHA-PARTICLES),
C   IF NOT EXPLICITLY SPECIFIED IN THE INPUT:
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      IF (VEL.LT.0.D0) THEN
       PD=XN_PD
      ELSE
       PD=1.94D-6*XN_PD*VEL**2
C
      ENDIF
C
C  RATIO OF PD TO THE AVERAGE PRESSURE, ASSUMED EQUAL TO 2 nPa:

      RAT=PD/2.0D0
      RAT16=RAT**0.14D0

C (THE POWER INDEX 0.14 IN THE SCALING FACTOR IS THE BEST-FIT VALUE OBTAINED FROM DATA
C    AND USED IN THE T96_01 VERSION)
C
C  VALUES OF THE MAGNETOPAUSE PARAMETERS FOR  PD = 2 nPa:
C
      A0=70.D0
      S00=1.08D0
      X00=5.48D0
C
C   VALUES OF THE MAGNETOPAUSE PARAMETERS, SCALED BY THE ACTUAL PRESSURE:
C
      A=A0/RAT16
      S0=S00
      X0=X00/RAT16
      XM=X0-A
C
C  (XM IS THE X-COORDINATE OF THE "SEAM" BETWEEN THE ELLIPSOID AND THE CYLINDER)
C
C     (FOR DETAILS OF THE ELLIPSOIDAL COORDINATES, SEE THE PAPER:
C      N.A.TSYGANENKO, SOLUTION OF CHAPMAN-FERRARO PROBLEM FOR AN
C      ELLIPSOIDAL MAGNETOPAUSE, PLANET.SPACE SCI., V.37, P.1037, 1989).
C
       IF (YGSW.NE.0.D0.OR.ZGSW.NE.0.D0) THEN
          PHI=DATAN2(YGSW,ZGSW)
       ELSE
          PHI=0.D0
       ENDIF
C
       RHO=DSQRT(YGSW**2+ZGSW**2)
C
       IF (XGSW.LT.XM) THEN
           XMGNP=XGSW
           RHOMGNP=A*DSQRT(S0**2-1.D0)
           YMGNP=RHOMGNP*DSIN(PHI)
           ZMGNP=RHOMGNP*DCOS(PHI)
           DIST=DSQRT((XGSW-XMGNP)**2+(YGSW-YMGNP)**2+(ZGSW-ZMGNP)**2)
           IF (RHOMGNP.GT.RHO) ID=+1
           IF (RHOMGNP.LE.RHO) ID=-1
           RETURN
       ENDIF
C
          XKSI=(XGSW-X0)/A+1.D0
          XDZT=RHO/A
          SQ1=DSQRT((1.D0+XKSI)**2+XDZT**2)
          SQ2=DSQRT((1.D0-XKSI)**2+XDZT**2)
          SIGMA=0.5D0*(SQ1+SQ2)
          TAU=0.5D0*(SQ1-SQ2)
C
C  NOW CALCULATE (X,Y,Z) FOR THE CLOSEST POINT AT THE MAGNETOPAUSE
C
          XMGNP=X0-A*(1.D0-S0*TAU)
          ARG=(S0**2-1.D0)*(1.D0-TAU**2)
          IF (ARG.LT.0.D0) ARG=0.D0
          RHOMGNP=A*DSQRT(ARG)
          YMGNP=RHOMGNP*DSIN(PHI)
          ZMGNP=RHOMGNP*DCOS(PHI)
C
C  NOW CALCULATE THE DISTANCE BETWEEN THE POINTS {XGSW,YGSW,ZGSW} AND {XMGNP,YMGNP,ZMGNP}:
C   (IN GENERAL, THIS IS NOT THE SHORTEST DISTANCE D_MIN, BUT DIST ASYMPTOTICALLY TENDS
C    TO D_MIN, AS WE ARE GETTING CLOSER TO THE MAGNETOPAUSE):
C
      DIST=DSQRT((XGSW-XMGNP)**2+(YGSW-YMGNP)**2+(ZGSW-ZMGNP)**2)
C
      IF (SIGMA.GT.S0) ID=-1   !  ID=-1 MEANS THAT THE POINT LIES OUTSIDE
      IF (SIGMA.LE.S0) ID=+1   !  ID=+1 MEANS THAT THE POINT LIES INSIDE
C                                           THE MAGNETOSPHERE
      RETURN
      END
C
C===================================================================================
C
c
