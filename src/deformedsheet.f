      SUBROUTINE deformedsheet (Mend,Nend, D0, RH0,G,TW,
     .     TSS,TSO,TSE,
     .     tilt, x,y,z,
     *     bxS,byS,bzS, bxO,byO,bzO, bxE,byE,bzE)
C     
C    CALCULATES GSM COMPONENTS OF 104 UNIT-AMPLITUDE TAIL FIELD MODES,
C    TAKING INTO ACCOUNT BOTH EFFECTS OF DIPOLE TILT:
C    WARPING IN Y-Z (DONE BY THE S/R WARPED) AND BENDING IN X-Z (DONE BY THIS SUBROUTINE)
C
      IMPLICIT NONE

C     Inputs
      INTEGER Mend,Nend
      REAL*8  D0,RH0,G,TW
      REAL*8  TSS(80,Nend)
      REAL*8  TSO(80,Mend,Nend)
      REAL*8  TSE(80,Mend,Nend)
      REAL*8  tilt
      REAL*8  x,y,z
      
C     Outputs
      REAL*8 bxS(Nend),bxO(Mend,Nend),bxE(Mend,Nend)
      REAL*8 byS(Nend),byO(Mend,Nend),byE(Mend,Nend)
      REAL*8 bzS(Nend),bzO(Mend,Nend),bzE(Mend,Nend)

C     Internal variables
      INTEGER L,K
      REAL*8  XAS,YAS,ZAS, XASS,YASS,ZASS
C      
      REAL*8 BXASS(Nend),BXASO(Mend,Nend),BXASE(Mend,Nend)
      REAL*8 BYASS(Nend),BYASO(Mend,Nend),BYASE(Mend,Nend)
      REAL*8 BZASS(Nend),BZASO(Mend,Nend),BZASE(Mend,Nend)
C
      REAL*8 BXASSS(Nend),BXASSO(Mend,Nend),BXASSE(Mend,Nend)
      REAL*8 BYASSS(Nend),BYASSO(Mend,Nend),BYASSE(Mend,Nend)
      REAL*8 BZASSS(Nend),BZASSO(Mend,Nend),BZASSE(Mend,Nend)

C     Constants
      REAL*8     rho0
      PARAMETER (rho0=20.0d0)
      
C
C     DEFORM:
C
      call bendpos(tilt, RH0, X,Y,Z, XAS,YAS,ZAS)
      call warppos(tilt, G,TW, XAS,YAS,ZAS, XASS,YASS,ZASS)
      call tailsheet_shielded (Mend,Nend, rho0, D0,
     .     TSS,TSO,TSE,
     .     XASS,YASS,ZASS,
     *     BXASSS,BYASSS,BZASSS,
     *     BXASSO,BYASSO,BZASSO,
     *     BXASSE,BYASSE,BZASSE)
      
C     
C --- New tail structure -------------
      
      DO K=1,Nend
C ------------------------------------------- Deforming symmetric modules
         call warpfield(BXASSS(K),BYASSS(K),BZASSS(K),
     .        BXASS(K),BYASS(K),BZASS(K))
         CALL bendfield(BXASS(K),BYASS(K),BZASS(K),
     .        bxS(K),byS(K),bzS(K))
         
      ENDDO
      
      DO K=1,Nend
         DO L=1,Mend

C -------------------------------------------- Deforming odd modules            
            CALL warpfield(BXASSO(L,K),BYASSO(L,K),BZASSO(L,K),
     .           BXASO(L,K),BYASO(L,K),BZASO(L,K))
            CALL bendfield(BXASO(L,K),BYASO(L,K),BZASO(L,K),
     .           bxO(L,K),byO(L,K),bzO(L,K))
            
C ------------------------------------------- Deforming even modules            
            CALL warpfield(BXASSE(L,K),BYASSE(L,K),BZASSE(L,K),
     .           BXASE(L,K),BYASE(L,K),BZASE(L,K))
            CALL bendfield(BXASE(L,K),BYASE(L,K),BZASE(L,K),
     .           bxE(L,K),byE(L,K),bzE(L,K))
            
         ENDDO
      ENDDO
      
C ------------------------------------
C
      RETURN
      END
