!
!© 2023 The Johns Hopkins University Applied Physics Laboratory LLC
!
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
! ts07d_example_june2023update.for ts07_field.f read_ts07.f fac_total_ts07.f ../deformedsheet.f ../dipoleshield.f ../fac_field.f ../fac_shield.f ../deformbirkfield.f ../stretchfield.f ../one_cone.f ../one_cone_smooth.f ../conical.f90 ../tailsheet_shielded.f ../tailsheet_sym.f ../tailsheet_asym.f ../bendfield.f ../warpfield.f ../cartharmonic.f ../cartharmonic_alt.f ../cylharmonic.f ../rotate_about_y.f ../bessjj.f -o ts07d_example_june2023update
!======================================================================================================
      PROGRAM TS07D
      
      implicit none 
c----------------------------------------------------------------------
c     Inputs to the model
      REAL*8 xGSM,yGSM,zGSM      

c     Output parameters for TS07D model
      REAL*8 bxGSM,byGSM,bzGSM

c     TODO: change the below to paths to your local paths
c     STAICDIR is the directory that contains the static shielding coefficients
      CHARACTER*256 staticDir
      PARAMETER    (staticDir = '/Users/stephgk1/Downloads/TAIL_PAR/')
      
c     VARIABLEDIR is the directory where the variable coefficient directories will be
      CHARACTER*256 variableDir
      PARAMETER(variableDir = '/Users/stephgk1/sst19_fortran/')

c
      INTEGER iYear,iDoY,iHour,min,iSec

c     The following reads the static coefficients that are used for the Equatorial Field
c     shielding fields. These only need to be loaded once as they are common for all times.      
      call store_ts07_static(staticDir)
      PRINT *, '   SHIELDING COEFFICIENTS HAS BEEN READ INTO RAM'

c     Now, initialize the TS07D and IGRF models for the supplied date and time. If multiple times
c     are wanted, the following must be recalled for the new date and time.
      
c     The date and time is specified by the year, day of year (1-366), hour (0-23), minute (0-59),
c     and the second (0-59), note that the second field is not used for accessing the variable
c     parameter and coefficients file.
      iYear = 2005
      iDoY = 253
      iHour = 0
      min = 0
      iSec = 0

      call store_ts07_variable_datetime(variableDir,
     .     iYear,iDoY,iHour,min,iSec)
      
c     Now the TS07D and the IGRF model can be evaluated

c     or just code in the values, e.g. XGSM = 5.0d0, units are in Earth radii (1 RE=6371.2 KM).
      xGSM = -5.0d0
      yGSM = +4.2d0
      zGSM = +1.2d0

c     Now evaluate the TS07D model, note that we have two variables holding the dipole tilt,
c     PSI which comes from Geopack and TILT which comes from the variables file, which
c     ultimately also was computed using Geopack. So use either one.
      call ts07_field(xGSM,yGSM,zGSM, bxGSM,byGSM,bzGSM)

      PRINT *,' '
      PRINT *,'   External field:'
      PRINT *,' BXGSM,BYGSM,BZGSM=',bxGSM,byGSM,bzGSM,' nT'
      PRINT *,' '
c//////////////////////////////////////////////////////////////////////////
c--------------------------------------------------------------------------

c     The program is complete
      END
