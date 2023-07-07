      SUBROUTINE deformedsheet_with_tcs (Mend,Nend, D,DTCS, RH0,G,TW,
     .     TSS,TSO,TSE,
     .     tilt, x,y,z,
     .     bxS1,byS1,bzS1, bxO1,byO1,bzO1, bxE1,byE1,bzE1,
     .     bxS2,byS2,bzS2, bxO2,byO2,bzO2, bxE2,byE2,bzE2)
C     
C    CALCULATES GSM COMPONENTS OF 104 UNIT-AMPLITUDE TAIL FIELD MODES,
C    TAKING INTO ACCOUNT BOTH EFFECTS OF DIPOLE TILT:
C    WARPING IN Y-Z (DONE BY THE S/R WARPED) AND BENDING IN X-Z (DONE BY THIS SUBROUTINE)
C
      IMPLICIT NONE

C     Inputs
      INTEGER Mend,Nend
      REAL*8  D,DTCS, RH0,G,TW
      REAL*8  TSS(80,Nend)
      REAL*8  TSO(80,Mend,Nend)
      REAL*8  TSE(80,Mend,Nend)
      REAL*8  tilt
      REAL*8  x,y,z
      
C     Outputs
      REAL*8 bxS1(Nend),bxO1(Mend,Nend),bxE1(Mend,Nend)
      REAL*8 byS1(Nend),byO1(Mend,Nend),byE1(Mend,Nend)
      REAL*8 bzS1(Nend),bzO1(Mend,Nend),bzE1(Mend,Nend)
C      
      REAL*8 bxS2(Nend),bxO2(Mend,Nend),bxE2(Mend,Nend)
      REAL*8 byS2(Nend),byO2(Mend,Nend),byE2(Mend,Nend)
      REAL*8 bzS2(Nend),bzO2(Mend,Nend),bzE2(Mend,Nend)
      
C     Internal variables
      INTEGER L,K
      REAL*8  XAS,YAS,ZAS, XASS,YASS,ZASS
C      
      REAL*8 BXASS1(Nend),BXASO1(Mend,Nend),BXASE1(Mend,Nend)
      REAL*8 BYASS1(Nend),BYASO1(Mend,Nend),BYASE1(Mend,Nend)
      REAL*8 BZASS1(Nend),BZASO1(Mend,Nend),BZASE1(Mend,Nend)
C
      REAL*8 BXASSS1(Nend),BXASSO1(Mend,Nend),BXASSE1(Mend,Nend)
      REAL*8 BYASSS1(Nend),BYASSO1(Mend,Nend),BYASSE1(Mend,Nend)
      REAL*8 BZASSS1(Nend),BZASSO1(Mend,Nend),BZASSE1(Mend,Nend)
C
      REAL*8 BXASS2(Nend),BXASO2(Mend,Nend),BXASE2(Mend,Nend)
      REAL*8 BYASS2(Nend),BYASO2(Mend,Nend),BYASE2(Mend,Nend)
      REAL*8 BZASS2(Nend),BZASO2(Mend,Nend),BZASE2(Mend,Nend)
C
      REAL*8 BXASSS2(Nend),BXASSO2(Mend,Nend),BXASSE2(Mend,Nend)
      REAL*8 BYASSS2(Nend),BYASSO2(Mend,Nend),BYASSE2(Mend,Nend)
      REAL*8 BZASSS2(Nend),BZASSO2(Mend,Nend),BZASSE2(Mend,Nend)
      
C     Constants
      REAL*8     rho0
      PARAMETER (rho0=20.0d0)
      
C
C     DEFORM:
C
      call bendpos(tilt, RH0, X,Y,Z, XAS,YAS,ZAS)
      call warppos(tilt, G,TW, XAS,YAS,ZAS, XASS,YASS,ZASS)
      call tailsheet_shielded (Mend,Nend, rho0, D,
     .     TSS,TSO,TSE,
     .     XASS,YASS,ZASS,
     *     BXASSS1,BYASSS1,BZASSS1,
     *     BXASSO1,BYASSO1,BZASSO1,
     *     BXASSE1,BYASSE1,BZASSE1)
      call tailsheet_shielded (Mend,Nend, rho0, DTCS,
     .     TSS,TSO,TSE,
     .     XASS,YASS,ZASS,
     *     BXASSS2,BYASSS2,BZASSS2,
     *     BXASSO2,BYASSO2,BZASSO2,
     *     BXASSE2,BYASSE2,BZASSE2)
      
C     
C --- New tail structure -------------
      
      DO K=1,Nend
C ------------------------------------------- Deforming symmetric modules
         call warpfield(BXASSS1(K),BYASSS1(K),BZASSS1(K),
     .        BXASS1(K),BYASS1(K),BZASS1(K))
         CALL bendfield(BXASS1(K),BYASS1(K),BZASS1(K),
     .        bxS1(K),byS1(K),bzS1(K))
C
         call warpfield(BXASSS2(K),BYASSS2(K),BZASSS2(K),
     .        BXASS2(K),BYASS2(K),BZASS2(K))
         CALL bendfield(BXASS2(K),BYASS2(K),BZASS2(K),
     .        bxS2(K),byS2(K),bzS2(K))         
      ENDDO
      
      DO K=1,Nend
         DO L=1,Mend

C -------------------------------------------- Deforming odd modules            
            CALL warpfield(BXASSO1(L,K),BYASSO1(L,K),BZASSO1(L,K),
     .           BXASO1(L,K),BYASO1(L,K),BZASO1(L,K))
            CALL bendfield(BXASO1(L,K),BYASO1(L,K),BZASO1(L,K),
     .           bxO1(L,K),byO1(L,K),bzO1(L,K))
C
            CALL warpfield(BXASSO2(L,K),BYASSO2(L,K),BZASSO2(L,K),
     .           BXASO2(L,K),BYASO2(L,K),BZASO2(L,K))
            CALL bendfield(BXASO2(L,K),BYASO2(L,K),BZASO2(L,K),
     .           bxO2(L,K),byO2(L,K),bzO2(L,K))            
            
C ------------------------------------------- Deforming even modules            
            CALL warpfield(BXASSE1(L,K),BYASSE1(L,K),BZASSE1(L,K),
     .           BXASE1(L,K),BYASE1(L,K),BZASE1(L,K))
            CALL bendfield(BXASE1(L,K),BYASE1(L,K),BZASE1(L,K),
     .           bxE1(L,K),byE1(L,K),bzE1(L,K))
C
            CALL warpfield(BXASSE2(L,K),BYASSE2(L,K),BZASSE2(L,K),
     .           BXASE2(L,K),BYASE2(L,K),BZASE2(L,K))
            CALL bendfield(BXASE2(L,K),BYASE2(L,K),BZASE2(L,K),
     .           bxE2(L,K),byE2(L,K),bzE2(L,K))            
            
         ENDDO
      ENDDO
      
C ------------------------------------
C
      RETURN
      END
