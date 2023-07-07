C     May 2019, G.K.Stephens, this subroutine was adabpted to the Saturn
C     magnetic field model.
C      
C     THIS S/R RETURNS THE SHIELDING FIELD FOR THE EARTH'S DIPOLE,
C     REPRESENTED BY  2x3x3=18 "CARTESIAN" HARMONICS, tilted with respect
C     to the z=0 plane (see NB#4, p.74-74)
C     
C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C     The 36 coefficients enter in pairs in the amplitudes of the "cartesian"
c     harmonics (A(1)-A(36).
c     The 14 nonlinear parameters (A(37)-A(50) are the scales Pi,Ri,Qi,and Si
C     entering the arguments of exponents, sines, and cosines in each of the
C     18 "Cartesian" harmonics  PLUS TWO TILT ANGLES FOR THE CARTESIAN HARMONICS
C     (ONE FOR THE PSI=0 MODE AND ANOTHER FOR THE PSI=90 MODE)
C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SUBROUTINE  dipoleshield(sunAngle, x,y,z, bx,by,bz)
c     
c     Computes the dipole shielding magnetic field for Saturn in the KSM
c     coordinate system in units of nT.
c
c     Inputs:
c     real*8 sunAngle: The latitude of the Sun's position relative to
c       Saturn's rotational equatorial plane. This angle is the same as
c       the dipole tilt angle.
c     real*8 x,y,z: The location to evaluate the value of the dipole
c       sheilding magnetic field in the KSM coordinate system in units
c       of Saturn radii.
c     
c     Outputs:
c     real*8 bx,by,bz: The dipole shielding magnetc field evaluated at
c       the supplied location in KSM coordinates in untis of nT.
c      
      IMPLICIT NONE
      
c     Inputs
      REAL*8 sunAngle
      REAL*8 x,y,z
      
c     Outputs
      REAL*8 bx,by,bz

      INTEGER Iend,Kend
      PARAMETER    (Iend = 3)
      PARAMETER    (Kend = 3)
      
c     Internal variables
      REAL*8 pis(1:Iend)
      REAL*8 pisI(1:Iend)
      DATA   pis /9.620648151, 6.082014949, 27.75216226/
      REAL*8 rks(1:Kend)
      REAL*8 rksI(1:Kend)
      DATA   rks /12.44199571, 5.122226936, 6.982039615/
      REAL*8 aPerp(1:Iend,1:Kend)
      REAL*8 aik(1:Iend,1:Kend)
c      DATA   aik /
c     .     -901.2327248, 817.6208321,
c     .     -83.73539535, 336.8781402,
c     .     -311.2947120, 31.94469304,
c     .     125.8739681, -235.4720434,
c     .     21.86305585/
      DATA   aik /
     .     -901.2327248, 336.8781402, 125.8739681,
     .     817.6208321, -311.2947120,-235.4720434,
     .     -83.73539535, 31.94469304, 21.86305585/
      REAL*8 bik(1:Iend,1:Kend)
c      DATA   bik /
c     .     895.8011176, -845.5880889,
c     .     86.58542841, -329.3619944,
c     .     308.6011161, -31.30824526,
c     .     -372.3384278, 286.7594095,
c     .     -27.42344605/
      DATA   bik /
     .     895.8011176, -329.3619944, -372.3384278,
     .     -845.5880889, 308.6011161, 286.7594095,
     .     86.58542841, -31.30824526, -27.42344605/
      
      REAL*8 qis(1:Iend)
      REAL*8 qisI(1:Iend)
      DATA   qis /20.12149582, 6.150973118, 4.663639687/
      REAL*8 sks(1:Kend)
      REAL*8 sksI(1:Kend)
      DATA   sks /15.73319647, 2.303504968, 5.840511214/
      REAL*8 aParr(1:Iend,1:Kend)
      REAL*8 cik(1:Iend,1:Kend)
c      DATA   cik /
c     .     -150.4874688, 1.395023949,
c     .     -56.85224007, -43.48705106,
c     .     1.073551279, 12.21404266,
c     .     5.799964188, -1.044652977,
c     .     3.536082962/
      DATA   cik /
     .     -150.4874688, -43.48705106,  5.799964188,
     .     1.395023949,   1.073551279, -1.044652977,
     .     -56.85224007,  12.21404266, 3.536082962/
      REAL*8 dik(1:Iend,1:Kend)
c      DATA   dik /
c     .     2.669338538, -.5540427503,
c     .     3.681827033, 5.103131905,
c     .     -.6673083508, 4.177465543,
c     .     -.3977802319, .5703560010,
c     .     -3.222069852/
      DATA   dik /
     .     2.669338538,   5.103131905, -.3977802319,
     .     -.5540427503, -.6673083508, .5703560010,
     .     3.681827033,   4.177465543, -3.222069852/
      
      integer i,k
      REAL*8 T1,T2, tilt,cosTilt,sinTilt,sin2Tilt, anglePerp,angleParr
      REAL*8 xr1,yr1,zr1, xr2,yr2,zr2
      REAL*8 bx1,by1,bz1, bxr1,byr1,bzr1
      REAL*8 bx2,by2,bz2, bxr2,byr2,bzr2

c     the dipole tilt angle is the same as the sun latitude angle
      tilt = sunAngle
      
c     invert the coefficients, the coefficients are relatively small,
c     so the inverse is fit instead, so before we evaluate them, we must
c     invert them back
      call invert_a(pis, 1,Iend, pisI)
      call invert_a(rks, 1,Kend, rksI)
      call invert_a(qis, 1,Iend, qisI)
      call invert_a(sks, 1,Kend, sksI)
      
      T1  = 0.08385953499d0
      T2  = 0.3477844929d0
      
      cosTilt = DCOS(tilt)
      sinTilt = DSIN(tilt)
      sin2Tilt = 2.0d0*sinTilt*cosTilt

c     Compute the coefficients that are a function of the dipole tilt
c     angle
      do i=1,Iend
         do k=1,Kend
            aPerp(i,k) = aik(i,k) + bik(i,k)*cosTilt
            aParr(i,k) = cik(i,k)*sinTilt + dik(i,k)*sin2Tilt
         enddo
      enddo
      
      anglePerp = tilt*T1
      angleParr = tilt*T2

c     MAKE THE TERMS IN THE 1ST SUM ("PERPENDICULAR" SYMMETRY):
      bx1 = 0.0d0
      by1 = 0.0d0
      bz1 = 0.0d0
      call rotate_about_y(anglePerp, x,y,z, xr1,yr1,zr1)
      call cartharmonic_alt(.true., .false., Iend,Kend,
     .     pisI,rksI,aPerp,
     .     xr1,yr1,zr1,
     .     bx1,by1,bz1)
      call rotate_about_y(-anglePerp, bx1,by1,bz1, bxr1,byr1,bzr1)
      
      
c     MAKE THE TERMS IN THE 2ND SUM ("PARALLEL" SYMMETRY):
      bx2 = 0.0d0
      by2 = 0.0d0
      bz2 = 0.0d0
      call rotate_about_y(angleParr, x,y,z, xr2,yr2,zr2)
      call cartharmonic(.true., .true., Iend,Kend,
     .     qisI,sksI,aParr,
     .     xr2,yr2,zr2,
     .     bx2,by2,bz2)
      call rotate_about_y(-angleParr, bx2,by2,bz2, bxr2,byr2,bzr2)

c     add the perpendicular and parallel fields together
      bx = bxr1+bxr2
      by = byr1+byr2
      bz = bzr1+bzr2      
       
      RETURN
      END
      
      
      SUBROUTINE invert_a(array, firstindex,lastindex, arrayI)
c     
c     Inverts the elements of an array of doubles.
c
c     Inputs:
c     real*8 array: The array whose elements will be inverted and
c       returned.
c     integer firstindex: The first index of the array.
c     integer lastindex: The last index of the array.
c     
c     Outputs:
c     real*8 arrayI: An array whose elements will be the inverse of the
c       elements of the supplied array.
c      
      
c     Inputs      
      integer i,firstindex,lastindex
      real*8 array(*)

c     Outputs
      real*8 arrayI(*)
      
      do i=firstindex,lastindex
         arrayI(i) = 1.0d0/array(i)
      enddo
      
      return
      END
