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
