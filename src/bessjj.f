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

