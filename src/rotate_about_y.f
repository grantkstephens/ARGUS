      SUBROUTINE  rotate_about_y(angle,x,y,z, xr,yr,zr)
c
c     Rotates a vector about the +Y axis by some angle.
c      
c     Inputs:
c     real*8 angle: The angle to ratate about the y axis.
c     real*8 x,y,z: The vector to be rotated.
c
c     Outputs:
c     real*8 xr,yr,zr: The rotated vector.
c     
      IMPLICIT NONE

c     Inputs
      REAL*8 angle,x,y,z
      
c     Outputs
      REAL*8 xr,yr,zr
      
c     Internal variables
      REAL*8 sinT, cosT

      sinT = DSIN(angle)
      cosT = DCOS(angle)

      xr = x*cosT - z*sinT
      yr = y
      zr = x*sinT + z*cosT
      
      END
