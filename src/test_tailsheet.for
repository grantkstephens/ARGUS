c     June 2019, G.K.Stephens: Initial prototype.
      program example_saturn
c
c     This is an example program that demonstrates how to evaluate both
c     the external and internal magnetic field (Dougherty et al. 2018)
c     for Saturn. The coordinate system is KSM (Kronocentric Solar
c     Magnetospheric) and units are in Saturn radii (Rs) and nanoTesla
c     (nT).
c
c     The program can be compiled using the following using gfortran:
c     
c     gfortran example_saturn.f -o example_saturn
c      
c     The expected output of this program is the following:
c     
c     Date:    2013-04-21 22:53:28.39823
c     TDT:        419856875.58223003     
c     sunLat:    0.31991904766201745     
c     (x,y,z):
c       -12.381273999999999    3.8901874300000001    2.1981273432987001
c     (bx,by,bz):
c        6.2702298378118524   -1.2576152407132852    7.0949159196904112
c     (hx,hy,hz):
c       0.74732969552778661   -1.1478438385940311   -9.4198143792161204
c     (r,theta,phi):
c        13.162874607051583    1.4030157792196860    2.8371604325738131
c     (br,btheta,bphi):
c       -7.7000229788157384    1.0538352925030736   0.19132520816042442
c     
      implicit none

      character*(512) pdynfile

      integer year,month,dom, hour,min
      real*8 sec

      real*8 pDyn

      real*8 utc_to_tdt,tdt
      
      real*8 sunLat,sunLong, zPhi,zTheta

      real*8 x,y,z
      
      real*8 bx,by,bz, hx,hy,hz, bTotx,bToty,bTotz
      real*8 r,theta,phi, br,btheta,bphi
      
c     The file containing the dynamic pressure of the solar wind at
c     Saturn. This file is derived from the Tao et al. (2005) model and
c     was obtained from the AMDA archive.
      pdynfile = "Dp.txt"
      pdynfile = trim(pdynfile)

c     Reads and stores the results of the Dynamic pressure file.
c     This call must be made once prior to evaluating the dynamic
c     pressure. Once called, this does not need called again.
      call store_pdyn(pdynfile)

c     The UTC date and time to evaluate the model. The format follows
c     the standard year-month-dayofmonth format, where these fields
c     start from 1, e.g. month=1 corresonds to January, and the time
c     fields start at 0. Year, month, dom, hour, and min are integer,
c     while sec is a double.
      year = 2013
      month = 4
      dom = 21
      hour = 22
      min = 53
      sec = 28.39823d0

c     Evaluates the dynamic pressure at Saturn for the given time.
c     Pressures are linearly interpolated.
      call eval_pdyn(year,month,dom, hour,min,sec,
     .     pdyn)

c     Evalute the Terrestrial Dynamical Time (TDT or Terrestrial Time or
c     TT) for the given UTC date and time.
      tdt = utc_to_tdt(year,month,dom, hour,min,sec)

c     Evaluates the location of the Sun at the supplied TDT time in
c     radians.
      call ksun(tdt, sunLat,sunLong, zPhi,zTheta)

c     The position vector in KSM coordinates in units of Saturn radii.
      x = -12.381274d0
      y = 3.89018743d0
      z = 2.1981273432987d0
      
c     Evaluates the external magnetic field in KSM coordinates in units
c     of nT.
      call saturn_ext(sunLat, pDyn,
     .     x,y,z,
     .     bx,by,bz)

c     Evaluates the internal magnetic field (Dougherty et al. 2018) in
c     KSM coordinates in units of nT.
      call saturn_int(sunLat, x,y,z, hx,hy,hz)

c     The total magnetic field is the sum of the external and internal
c     fields.
      bTotx = bx + hx
      bToty = by + hy
      bTotz = bz + hz

c     Converts the Cartesian total field value to spherical coordinates.
      call cart_field_to_sphere(x,y,z,
     .     bTotx,bToty,bTotz,
     .     r,theta,phi, br,btheta,bphi)

c     Prints out the results.
      write(*,'(A,I4,A,I0.2,A,I0.2,A,I0.2,A,I0.2,A,F8.5)')
     .     "Date:    ",year,"-",month,"-",dom,
     .     " ",hour,":",min,":",sec
      write(*,*) "TDT:     ",TDT
      write(*,*) "sunLat:  ",sunLat
      write(*,*) "(x,y,z):    ",x,y,z
      write(*,*) "(bx,by,bz):    ",bx,by,bz
      write(*,*) "(hx,hy,hz):    ",hx,hy,hz
      write(*,*) "(r,theta,phi):    ",r,theta,phi
      write(*,*) "(br,btheta,bphi):    ",br,btheta,bphi
      
      end program
      
