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
! gfortran ts07d_geopack_example_june2023update.for Geopack-2008_dp.f bendfield.f fac_total.f fac_field.f fac_shield.f conical.f90 deformbirkfield.f dipoleshield.f stretchfield.f tailsheet_shielded.f tailsheet_asym.f tailsheet_sym.f warpfield.f bessjj.f cartharmonic.f cartharmonic_alt.f cylharmonic.f rotate_about_y.f read_ts07.f one_cone.f deformedsheet.f ts07_field.f -o ts07d_geopack_example_june2023update
!======================================================================================================
      P r o g r a m   TS07D
      implicit none 
c----------------------------------------------------------------------
c     Inputs to the model
      REAL*8 XGSM,YGSM,ZGSM      

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

c
      REAL*8 VXGSE,VYGSE,VZGSE,PI
      INTEGER IYEAR,IDAY,IHOUR,MIN,ISEC

c     Time variables, delete when done
      REAL*8 evaltime,start,finish

c--------------------------------------------------------------------------
c//////////////////////////////////////////////////////////////////////////
c     The following reads the static coefficients that are used for the Equatorial Field
c     shielding fields. These only need to be loaded once as they are common for all times.      
      call store_ts07_static(STATICDIR)
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

      call store_ts07_variable_datetime(VARIABLEDIR,
     .     iYear,IDAY,iHour,min,iSec)
      
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

c     Now evaluate the TS07D model, note that we have two variables holding the dipole tilt,
c     PSI which comes from Geopack and TILT which comes from the variables file, which
c     ultimately also was computed using Geopack. So use either one.
      EVALTIME   = 0.d0
      call cpu_time(start)
      call ts07_field(XGSM,YGSM,ZGSM, BXGSM,BYGSM,BZGSM)
c      CALL TS07D_JUNE_2023 (IOPT,PARMOD,PSI,XGSM,YGSM,ZGSM,
c     *     BXGSM,BYGSM,BZGSM)
      PRINT *,' BXGSM,BYGSM,BZGSM=',BXGSM,BYGSM,BZGSM,' nT'      
      call ts07_field_tilt(PSI, XGSM,YGSM,ZGSM,
     *     BXGSM,BYGSM,BZGSM)
      PRINT *,' BXGSM,BYGSM,BZGSM=',BXGSM,BYGSM,BZGSM,' nT'
      call cpu_time(finish)
      evaltime = evaltime+finish-start

      write(*,*) 'Eval time: ', evaltime
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
      SUBROUTINE TS07D_JUNE_2023 (IOPT,PARMOD,PS,X,Y,Z,BX,BY,BZ)     !  A DOUBLE-PRECISION SUBROUTINE
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
C     /TSO/ TSO - An 80x4x5-element array containing the static (time-independent) coefficients that are
C        used to shield the ODD axis-symmetric equatorial expansions
C     /TSE/ TSE - An 80x4x5-element array containing the static (time-independent) coefficients that are
C        used to shield the EVEN axis-symmetric equatorial expansions 
C
C  Outputs:
C     BX,BY,BZ - the evaluated magnetic field vector in the GSM coordinate system in units of nanoTesla (nT)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8    PARMOD(10)

      PARAMETER (NTOT=101)

      COMMON /INPUT/  PDYN
      DIMENSION BXTS(5),BXTO(4,5),BXTE(4,5)
      DIMENSION BYTS(5),BYTO(4,5),BYTE(4,5)
      DIMENSION BZTS(5),BZTO(4,5),BZTE(4,5)
      COMMON /PARAM/ A(NTOT)
      COMMON /TSS/ TSS(80,5)
      COMMON /TSO/ TSO(80,4,5)
      COMMON /TSE/ TSE(80,4,5)
      
      CALL EXTERN (IOPT,A,NTOT, TSS,TSO,TSE,
     *PS,PDYN,X,Y,Z,BXCF,BYCF,BZCF,
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
     *     TSS,TSO,TSE,
     *     tilt,PDYN,X,Y,Z,
     *     BXCF,BYCF,BZCF,
     *     BXTS,BYTS,BZTS,BXTO,BYTO,BZTO,BXTE,BYTE,BZTE,
     *     BXR11a,BYR11a,BZR11a,BXR12a,BYR12a,BZR12a,
     *     BXR21a,BYR21a,BZR21a,BXR21s,BYR21s,BZR21s,
     *     BX,BY,BZ)
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
      INTEGER IOPGEN
      REAL*8 A(NTOT)
      INTEGER NTOT
      REAL*8 TSS(80,Nend)
      REAL*8 TSO(80,Mend,Nend)
      REAL*8 TSE(80,Mend,Nend)
      REAL*8 tilt,PDYN
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
      REAL*8 A0_A /34.586D0/,A0_S0 /1.1960D0/,A0_X0 /3.4397D0/   !   SHUE ET AL. PARAMETERS
      REAL*8 DSIG /0.005D0/, RH2 /-5.2D0/
      REAL*8 X0,AM,S0
C
      REAL*8  BXR22a,BYR22a,BZR22a
      REAL*8  BXR11s,BYR11s,BZR11s,BXR12s,BYR12s,BZR12s
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
      REAL*8 A_dipsh
      REAL*8 A_eq(Nend+2*Mend*Nend)
      REAL*8 A_eq_P(Nend+2*Mend*Nend)
      REAL*8 A_R11a,A_R12a,A_R21a,A_R21s
      REAL*8 D                  ! TAIL SHEET THICKNESS
      REAL*8 RH0,G
      REAL*8 xKappa1,xKappa2    !  SCALING FACTORS FOR BIRKELAND CURRENTS
      REAL*8 TW

      Pdyn_0=2.d0               !   AVERAGE PRESSURE USED FOR NORMALIZATION
      
      XAPPA=(PDYN/2.D0)**0.155   !   0.155 is the value obtained in TS05
      XAPPA3=XAPPA**3

      A_dipsh=A(1)
      
      A_eq=  A(2:46)
      A_eq_P=A(47:91)
      
      A_R11a=A(92)
      A_R12a=A(93)
      A_R21a=A(94)
      A_R21s=A(95)
      
      D=      A(96)
      RH0=    A(97)
      G=      A(98)
      xKappa1=A(99)
      xKappa2=A(100)
      TW=     A(101)       !   THIS PARAMETER CONTROLS THE IMF-INDUCED TWISTING (ADDED 04/21/06)
c
      XX=X*XAPPA  ! pressure scaling has been reinstated here
      YY=Y*XAPPA
      ZZ=Z*XAPPA
      
C
      IF (IOPGEN.LE.1) THEN

C     print *,XX,YY,ZZ
         call dipoleshield(tilt, XX,YY,ZZ, CFX,CFY,CFZ) !  DIPOLE SHIELDING FIELD
         BXCF=CFX  *XAPPA3
         BYCF=CFY  *XAPPA3
         BZCF=CFZ  *XAPPA3
      ELSE
         BXCF=0.D0
         BYCF=0.D0
         BZCF=0.D0
      ENDIF                                              !  DONE

      IF (IOPGEN.EQ.0.OR.IOPGEN.EQ.2) THEN
         CALL deformedsheet (Mend,Nend, D, RH0,G,TW, TSS,TSO,TSE, !  TAIL FIELD (THREE MODES)
     *        tilt,XX,YY,ZZ,      
     *        BXTS,BYTS,BZTS, BXTO,BYTO,BZTO, BXTE,BYTE,BZTE)
      ELSE
         DO n=1,Nend
            
            BXTS(n)=0.d0
            BYTS(n)=0.d0
            BZTS(n)=0.d0
            
         ENDDO

         DO  n=1,Nend
           DO m=1,Mend

              BXTO(m,n)=0.d0
              BYTO(m,n)=0.d0
              BZTO(m,n)=0.d0
              
              BXTE(m,n)=0.d0
              BYTE(m,n)=0.d0
              BZTE(m,n)=0.d0
              
           ENDDO
        ENDDO
      ENDIF
      
      IF (IOPGEN.EQ.0.OR.IOPGEN.EQ.3) THEN
         
         CALL fac_total (xKappa1,xKappa2,
     .        tilt,XX,YY,ZZ,
     .        BXR11a,BYR11a,BZR11a,BXR12a,BYR12a,BZR12a,  ! BIRKELAND FIELD (TWO MODES FOR R1 AND TWO MODES FOR R2)
     *        BXR21a,BYR21a,BZR21a,BXR22a,BYR22a,BZR22a, 
     .        BXR11s,BYR11s,BZR11s,BXR12s,BYR12s,BZR12s,  ! "SYMMETRIC" BIRKELAND FIELD
     .        BXR21s,BYR21s,BZR21s,BXR22s,BYR22s,BZR22s)  ! (TWO MODES FOR R1s AND TWO MODES FOR R2s)
c                                                           (but we actually use from here only R2s modes)
      ELSE
         BXR11a=0.D0
         BYR11a=0.D0
         BZR11a=0.D0
         BXR12a=0.D0
         BYR12a=0.D0
         BZR12a=0.D0
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
      
      RETURN
      END
C
C XXXXXXXXXXXXXXXXXXXXXXXXXXX11/15/05 16:06 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      
