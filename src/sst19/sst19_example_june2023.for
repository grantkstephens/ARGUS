!
!© 2023 The Johns Hopkins University Applied Physics Laboratory LLC
!
!===============================================================================
! June 2023: G.K.Stephens (Grant.Stephens@jhuapl.edu)
! This program demonstrates how to setup and evaluate the Stephens et al. (2019)
! empirical geomagnetic field model, termed SST19,
! (https://doi.org/10.1029/2018JA025843).
!
! The model provides the value of the external magnetic field, Bext, as a
! function of time and position in the Geocentric Solar Magnetospheric (GSM)
! system. The model’s morphological structure captures the current systems
! associated with magnetospheric substorms while a data-mining algorithm
! determines its dynamical evolution.
!
! The total time-dependent magnetic field, Btot=Bext+Bint, can be evaluated by
! combining this model with the International Geomagnetic Reference Field
! (IGRF-13) as implemented in Professor Tsyganenko’s Geopack library
! (https://geo.phys.spbu.ru/~tsyganenko/empirical-models/coordinate_systems/geopack).
!
! The source code is programmed in the FORTRAN programming language loosely
! following the FORTRAN 77 standards and is compiled using GFortran
! (https://gcc.gnu.org/wiki/GFortran). This example program can be compiled
! using the following on the command line:
! 
! gfortran sst19_example_june2023.for sst19_field.f read_sst19.f fac_total_sst19.f ../deformedsheet_with_tcs.f ../dipoleshield.f ../fac_field.f ../fac_shield.f ../deformbirkfield.f ../stretchfield.f ../one_cone.f ../one_cone_smooth.f ../conical.f90 ../tailsheet_shielded.f ../tailsheet_sym.f ../tailsheet_asym.f ../bendfield.f ../warpfield.f ../cartharmonic.f ../cartharmonic_alt.f ../cylharmonic.f ../rotate_about_y.f ../bessjj.f -o sst19_example_june2023
!
! This command compiles the program into an executable titled
! sst19_example_june2023. To run the program on a Unix based system, type:
!
! ./sst19_example_june2023
!
!===============================================================================
      PROGRAM SST19
      
      implicit none

c     The input to the model is the position in GSM coordinates in units of
c     Earth radii (1 RE=6371.2 KM) and the output is the vector external
c     magnetic field evaluated at that position in units of nanotesla (nT).
c     However, before the model can be evaluated, two initialization routines
c     must be called.
c      
c     The first, reads and stores a set of static shielding coefficients and the
c     path to a directory containing the files must be provided. This only needs
c     to be called once at the beginning of the program. The second, reads and
c     stores the time-variable coefficients and non-linear parameters that are
c     defined in an ASCII file. This must be called whenever the user wishes to
c     change the time.

c     Inputs to the model
      REAL*8 xGSM,yGSM,zGSM      

c     Output parameters for TS07D model
      REAL*8 bxGSM,byGSM,bzGSM

c     TODO: change the below to paths to your local paths
c     staticDir is the directory that contains the static shielding coefficients
      CHARACTER*256 staticDir
      PARAMETER    (staticDir = '/Users/stephgk1/examples/'//
     .     'staticcoefficients/coeffs_n40_m8')
      
c     variableDir is the directory containing the ASCII files with the
c     time-variable coeffients and non-linear parameters
      CHARACTER*256 variableDir
      PARAMETER(variableDir='/Users/stephgk1/sst19_fortran/src/sst19/')

c     filename for the time-variable ASCII file, set below
      CHARACTER*256 filename
      
c     Date and time fields, set below
      INTEGER iYear,iDoY,iHour,min
      
c     The first initialization routine reads and stores the static coefficients
c     that are used for the equatorial current shielding fields. These only need
c     to be loaded once as they are common for all times.
      call store_sst19_static(staticDir)
      PRINT *, '   shielding coefficients have been stored'

c     The second initialization routine reads and stores the time-variable
c     amplitude coefficients and non-linear parameters. These are saved in
c     human-readable ASCII files. There are two equivalent routines that can be
c     called to perform this. The first, the user provides the filename of the
c     ASCII file directly. The second, the user provides the directory
c     containing the ASCII files and the date and time. The filename standard
c     for these files is YEAR_DOY_HR_MN.par. Either routine must be called each
c     time the user wishes to switch to a new time.

c     Read and store the time-variable coeffs and params by providing the
c     filename.
      filename=trim(variableDir)//"2009_067_11_25.par"
      call store_sst19_variable(filename)
      
c     The date and time is specified by the year, day of year (1-366),
c     hour (0-23), and minute (0-59).
      iYear = 2009
      iDoY = 67
      iHour = 11
      min = 25

c     Alternatively, read and store the time-variable coeffs and params by
c     providing the directory and the date and time. Note, either this routine
c     or the above routine need to be called not both. They are both called
c     here just as an example.
      call store_sst19_variable_datetime(variableDir,
     .     iYear,iDoy,iHour,min)

c     Set the position, units are in Earth radii (1 RE=6371.2 KM).
      xGSM = -5.0d0
      yGSM = +4.2d0
      zGSM = +1.2d0

c     Now evaluate the SST19 model. This routine can be called any number of
c     times.
      call sst19_field(xGSM,yGSM,zGSM, bxGSM,byGSM,bzGSM)
      
      PRINT *,' External field evaluated at: (xGSM,yGSM,zGSM)=',
     .     xGSM,yGSM,zGSM
      PRINT *,' Expected:  (bxGSM,byGSM,bzGSM)='//
     .     '   16.702390096106559       -13.306583526317164'//
     .     '       -36.650180176716844       nT'
      PRINT *,' Evaluated: (bxGSM,byGSM,bzGSM)=',bxGSM,byGSM,bzGSM,' nT'
      PRINT *,' '

c     The program is complete
      END
