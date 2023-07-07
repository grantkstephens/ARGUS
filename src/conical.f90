!
!-----------------------------------------------------------------------
!     Subroutine conical
!-----------------------------------------------------------------------
!
subroutine conical(N,theta0,dTheta, isEven, r,theta,phi, bTheta,bPhi)

! Calculates the theta and phi components of the magnetic field of the conical current
! system described by Tsyganenko (1991). Adapted from Tsyganenko's FIALCOS subroutine.
! See Tsyganenko (1991): https://doi.org/10.1016/0032-0633(91)90058-I
!
!------------input parameters:
!
! N           - mode
! theta0      - polar angle defining the cone (radians)
! dTheta      - angular half width of the current layer (radians)
! isEven      - ??
! r,theta,phi - position in spherical coordinates (radians)
!
!------------output parameters:
! bTheta,bPhi - field components in spherical coordinates, Br=0

  implicit none

  integer*4   N
  real*8      theta0,dTheta
  logical     isEven
  real*8      r,theta,phi
  real*8      bTheta,bPhi
  
  integer*4   M
  real*8      sinTheta,cosTheta
  real*8      sinMPhi,cosMPhi
  real*8      sinPhi,cosPhi
  real*8      sinM1Phi,cosM1Phi
  real*8      T,TG,TG21,CTG
  real*8      thetaNP,thetaNM
  real*8      COSM1,SINM1,TM,DTT,DTT0
  real*8      FC,FC1,RO,TGM,TGM2,TGM2M,TGM2M1,TGP,TGP2,TGP2M
  real*8      CCOS(N),SSIN(N),BTN(N),BPN(N)
  
  sinTheta=sin(theta)
  cosTheta=cos(theta)
  sinPhi=sin(phi)
  cosPhi=cos(phi)

!  sinMphi=sin(N*phi)
!  cosMphi=cos(N*phi)
!  if (isEven .eqv. .true.) then
!     sinMphi=sin(N*phi)
!     cosMphi=cos(N*phi)
!  endif
  
  ! r*sin(theta)
  RO=R*sinTheta
  
  TG=sinTheta/(1.0d0+cosTheta)   !        TAN(THETA/2)
  CTG=sinTheta/(1.0d0-cosTheta)  !        CTG(THETA/2)

  ! the theta boundaries of the cone
  thetaNP=theta0+dTheta
  thetaNM=theta0-dTheta
  
  if (theta .ge. thetaNM) then
     TGP=tan(thetaNP*0.5D0)
     TGM=tan(thetaNM*0.5D0)
     TGM2=TGM*TGM
     TGP2=TGP*TGP
  endif

  ! Initialize recursive values
  cosM1Phi=1.0d0 ! cos(0) = 1
  sinM1Phi=0.0d0 ! sin(0) = 0
  TM=1.0d0
  TGM2M=1.0d0
  TGP2M=1.0d0

  ! Compute all the terms
  do M=1,N

     ! tan(theta/2)^m
     TM=TM*TG

     ! Chebyshev recursive method for computing cos(m*phi) and sin(m*phi) 
     CCOS(M)=cosM1Phi*cosPhi-sinM1Phi*sinPhi
     SSIN(M)=sinM1Phi*cosPhi+cosM1Phi*sinPhi

     ! Store the previous values for the next time through the loop
     cosM1Phi=CCOS(M) ! cos((m-1)*phi)
     sinM1Phi=SSIN(M) ! sin((m-1)*phi)

     FC=1.D0/(TGP-TGM)
     FC1=1.D0/(2*M+1)
     
     ! Evaluate T(theta)
     ! Check 3 cases, (1) theta < thetaNM, (2) thetaNM < theta < thetaNP, (3) theta > thetaNP
     if (theta .lt. thetaNM) then
        T=TM
        DTT=0.5D0*M*TM*(TG+CTG)
     else if (theta .lt. thetaNP) then
        TGM2M=TGM2M*TGM2
        TGM2M1=TGM2M*TGM
        TG21=1.D0+TG*TG
        T=FC*(TM*(TGP-TG)+FC1*(TM*TG-TGM2M1/TM))
        DTT=0.5D0*M*FC*TG21*(TM/TG*(TGP-TG)-FC1*(TM-TGM2M1/(TM*TG)))
     else
        TGP2M=TGP2M*TGP2
        TGM2M=TGM2M*TGM2
        T=FC*FC1*(TGP2M*TGP-TGM2M*TGM)/TM
        DTT=-T*M*0.5D0*(TG+CTG)
     endif

     sinMPhi=SSIN(M)
     cosMPhi=CCOS(M)
     if (isEven .eqv. .true.) then
        sinMphi=CCOS(M)
        cosMphi=-SSIN(M)
     endif
     
     ! Now compute Btheta and Bphi for the mth term, eq. (15)
     BTN(M)=M*T*cosMPhi/RO ! Btheta = m*T*cos(m*phi)/(r*sin(theta))
     BPN(M)=-DTT*sinMPhi/R ! Bphi   = -dT/dtheta*sin(m*phi)/r
     
  enddo
  
  bTheta=BTN(N)
  bPhi  =BPN(N)
  
  return
end subroutine conical
