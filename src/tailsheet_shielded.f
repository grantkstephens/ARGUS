C==========================================================================
C     May 2019, G.K.Stephens, this subroutine was adabpted to the Saturn
C     magnetic field model.
C      
C     July 2017, G.K.Stephens, this routine updated to incorporate Jay Albert's
C     improvements to the Bessel function evaluation.
c
C     CALCULATES GSM COMPONENTS OF THE SHIELDED FIELD OF 45 TAIL MODES WITH UNIT
C     AMPLITUDES,  WITHOUT ANY WARPING OR BENDING.  NONLINEAR PARAMETERS OF THE MODES
C     ARE FORWARDED HERE VIA A COMMON BLOCK /TAIL/.
      SUBROUTINE tailsheet_shielded (Mend,Nend,rho0,D0,
     .     TSS,TSO,TSE,
     .     x,y,z,
     .     bxs,bys,bzs, bxo,byo,bzo, bxe,bye,bze)
c
c     Documentation needed.
c
c     Inputs:
c     integer Mend:
c     ingeger Nend:
c     real*8 rho0:
c     real*8 D0:
c     real*8 x,y,z:
c
c     Outputs:
c     real*8 bxs,bys,bzs:
c     real*8 bxo,byo,bzo:
c     real*8 bxe,bye,bze:
c     
      IMPLICIT NONE

C     Inputs
      INTEGER Mend,Nend
      REAL*8 rho0,D0
      REAL*8 TSS(80,Nend),TSO(80,Mend,Nend),TSE(80,Mend,Nend)
      REAL*8 x,y,z

c     Outputs
      REAL*8 bxs(Nend),bxo(Mend,Nend),bxe(Mend,Nend)
      REAL*8 bys(Nend),byo(Mend,Nend),bye(Mend,Nend)
      REAL*8 bzs(Nend),bzo(Mend,Nend),bze(Mend,Nend)

c     Internal variables
      integer m,n, j
      REAL*8 rho,knRho, AJM(0:Nend),AJMD(0:Nend)
      REAL*8 BXSn,BYSn,BZSn, BXOnm,BYOnm,BZOnm, BXEnm,BYEnm,BZEnm
      REAL*8 HXSn,HYSn,HZSn, HXOnm,HYOnm,HZOnm, HXEnm,HYEnm,HZEnm
      
c      RNOT=20.0     !    Rho_0 - scale parameter along the tail axis

      rho=DSQRT(x*x+y*y)
      
      DO n=1,Nend
         
         knRho = n*rho/rho0
         
C     July 2017, G.K.Stephens, all the Bessel functions are now evaluated first,
C     and passed into the subroutines
         call bessJJ(Nend,knRho, AJM) !!! get all n in one call
         DO j=1,Nend
            AJMD(j)=AJM(j-1)-j*AJM(j)/knRho
         ENDDO
         AJMD(0)=-AJM(1)
         
         CALL tailsheet_sym (n,rho0,D0,AJM,
     .        x,y,z,
     .        BXSn,BYSn,BZSn)
         CALL cylharmonic (.false., TSS(1:80,n), X,Y,Z,
     .        HXSn,HYSn,HZSn)
         
         BXS(n)=BXSn+HXSn
         BYS(n)=BYSn+HYSn
         BZS(n)=BZSn+HZSn
         
         DO m=1,Mend
            
            CALL tailsheet_asym (.false.,m,n,rho0,D0,AJM,AJMD,
     .           x,y,z,
     .           BXOnm,BYOnm,BZOnm)
            CALL cylharmonic (.false.,  TSO(1:80,m,n), X,Y,Z,
     .           HXOnm,HYOnm,HZOnm)
            BXO(m,n)=BXOnm+HXOnm
            BYO(m,n)=BYOnm+HYOnm
            BZO(m,n)=BZOnm+HZOnm
            
            CALL tailsheet_asym (.true.,m,n,rho0,D0,AJM,AJMD,
     .           x,y,z,
     .           BXEnm,BYEnm,BZEnm)
            CALL cylharmonic (.true., TSE(1:80,m,n), X,Y,Z,
     .           HXEnm,HYEnm,HZEnm)
            BXE(m,n)=BXEnm+HXEnm
            BYE(m,n)=BYEnm+HYEnm
            BZE(m,n)=BZEnm+HZEnm
            
         ENDDO
      ENDDO
      
      RETURN
      END

      
      SUBROUTINE tailsheet_shielded_sum (Mend,Nend,rho0,D0,
     .     sym,odd,even,
     .     TSS,TSO,TSE,
     .     x,y,z,
     .     bx,by,bz)

      IMPLICIT NONE
      
C     Inputs
      INTEGER Mend,Nend
      REAL*8 rho0,D0
      REAL*8 sym(Nend),odd(Mend,Nend),even(Mend,Nend)
      REAL*8 TSS(80,Nend),TSO(80,Mend,Nend),TSE(80,Mend,Nend)
      REAL*8 x,y,z

c     Outputs
      REAL*8 bx,by,bz

c     Internal variables
      REAL*8 bxs(Nend),bxo(Mend,Nend),bxe(Mend,Nend)
      REAL*8 bys(Nend),byo(Mend,Nend),bye(Mend,Nend)
      REAL*8 bzs(Nend),bzo(Mend,Nend),bze(Mend,Nend)
      integer m,n

      bx = 0.0d0
      by = 0.0d0
      bz = 0.0d0
      
      call tailsheet_shielded(Mend,Nend,rho0,D0,
     .     TSS,TSO,TSE,
     .     x,y,z,
     .     bxs,bys,bzs, bxo,byo,bzo, bxe,bye,bze)
      
      DO n=1,Nend
         
         bx = bx + sym(n)*BXS(n)
         by = by + sym(n)*BYS(n)
         bz = bz + sym(n)*BZS(n)
         
         DO m=1,Mend

            bx = bx + odd(m,n)*BXO(m,n)
            by = by + odd(m,n)*BYO(m,n)
            bz = bz + odd(m,n)*BZO(m,n)
            
            bx = bx + even(m,n)*BXE(m,n)
            by = by + even(m,n)*BYE(m,n)
            bz = bz + even(m,n)*BZE(m,n)

         ENDDO
      ENDDO

      RETURN
      END
