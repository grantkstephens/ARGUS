C
C-------------------------------------------------------------------------
C
C
      SUBROUTINE fac_shield (A,X_SC,isIeven, tilt,X,Y,Z, BX,BY,BZ)
C
      IMPLICIT NONE

c     Inputs
      REAL*8 A(86)
      REAL*8 X_SC
      LOGICAL isIeven
      REAL*8 tilt
      REAL*8 x,y,z

c     Outputs
      REAL*8 bx,by,bz

c     Internal variables      
      INTEGER Iend,Kend
      PARAMETER    (Iend = 3)
      PARAMETER    (Kend = 3)         
      
      REAL*8 pis(1:Iend)
      REAL*8 pisI(1:Iend)
      REAL*8 rks(1:Iend)
      REAL*8 rksI(1:Iend)
      REAL*8 qis(1:Iend)
      REAL*8 qisI(1:Iend)
      REAL*8 sks(1:Iend)
      REAL*8 sksI(1:Iend)
      
      REAL*8 aPerp(1:Iend,1:Kend)
      REAL*8 aParr(1:Iend,1:Kend)

      integer m,i,k,ic
      REAL*8 cosTilt,sinTilt,sin3Tilt, anglePerp,angleParr
      REAL*8 xr1,yr1,zr1, xr2,yr2,zr2
      REAL*8 bx1,by1,bz1, bxr1,byr1,bzr1
      REAL*8 bx2,by2,bz2, bxr2,byr2,bzr2

      bx=0.d0
      by=0.d0
      bz=0.d0
      
      cosTilt=DCOS(tilt)
      sinTilt=DSIN(tilt)
      sin3Tilt=2.d0*cosTilt

      pis(1:Iend)=A(73:75)
      qis(1:Iend)=A(79:81)
      rks(1:Kend)=A(76:78)
      sks(1:Kend)=A(82:84)

      ic=1
      
      do m=1,2
         do i=1,3
            do k= 1,3
               if (m.eq.1) then
                  aPerp(i,k) = A(ic) + X_SC * A(ic+1) + cosTilt *
     .                 (A(ic+2) + X_SC * A(ic+3))
                  ic=ic+4
               else
                  aParr(i,k) = sinTilt * (A(ic) + X_SC * A(ic+1) +
     .                 sin3Tilt * (A(ic+2) + X_SC * A(ic+3)))
                  ic=ic+4                  
               endif
            enddo
         enddo
      enddo
      
c     invert the coefficients, the coefficients are relatively small,
c     so the inverse is fit instead, so before we evaluate them, we must
c     invert them back
      call invert_a(pis, 1,Iend, pisI)
      call invert_a(rks, 1,Kend, rksI)
      call invert_a(qis, 1,Iend, qisI)
      call invert_a(sks, 1,Kend, sksI)
      
      anglePerp = tilt*A(85)
      angleParr = tilt*A(86)

c     MAKE THE TERMS IN THE 1ST SUM ("PERPENDICULAR" SYMMETRY):
      bx1 = 0.0d0
      by1 = 0.0d0
      bz1 = 0.0d0
      call rotate_about_y(anglePerp, x,y,z, xr1,yr1,zr1)
      call cartharmonic(isIeven, .false., Iend,Kend,
     .     pisI,rksI,aPerp,
     .     xr1,yr1,zr1,
     .     bx1,by1,bz1)
      call rotate_about_y(-anglePerp, bx1,by1,bz1, bxr1,byr1,bzr1)
      
c     MAKE THE TERMS IN THE 2ND SUM ("PARALLEL" SYMMETRY):
      bx2 = 0.0d0
      by2 = 0.0d0
      bz2 = 0.0d0
      call rotate_about_y(angleParr, x,y,z, xr2,yr2,zr2)
      call cartharmonic(isIeven, .true., Iend,Kend,
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
