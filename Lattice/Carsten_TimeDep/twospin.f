      PROGRAM TWOSPIN
      IMPLICIT NONE

      INTEGER NP,NB,I,J,LWORK,INFO
      INCLUDE 'dim.inc'
      PARAMETER (NB=2*NP**2-NP)
      INTEGER IPIV(NB)
      PARAMETER (LWORK = 500)

      DOUBLE PRECISION U0,U1,T,V(NP),BX(NP),BY(NP),BZ(NP),MM(NP),
     &                 E(NB),RWORK(100),N(NP),MX(NP),MY(NP),MZ(NP)
      DOUBLE PRECISION TIME,DT,OMEGA,MTOT(3),BXD(NP)
      DOUBLE COMPLEX M(NB,NB),A(NB,NB),C(NB),R(NB),WORK(LWORK),
     &               ZERO,ONE,IONE
      PARAMETER (ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),IONE=(0.D0,1.D0))

      T = 0.5D0
      WRITE(*,*)'U0?'
      READ(*,*)U0
      U1 = U0/2.D0

      OMEGA = 0.5D0 
      DT = 0.01D0
C**---------------------------------------------------------------------
C**   First, calculate the ground state
C**---------------------------------------------------------------------
      DO I=1,NP
         READ(1,*)V(I),BX(I),BY(I),BZ(I)
      ENDDO    

      CALL MATRIX(M,U0,U1,T,V,BX,BY,BZ)

      CALL ZHEEV( 'V', 'U', NB, M, NB, E, WORK, LWORK, RWORK, INFO )
     
      DO I=1,NB
         C(I) = M(I,1)
C        WRITE(*,*)I,real(E(I))
      ENDDO    
      CALL DENCALC(C,N,MX,MY,MZ,MM,MTOT)

      TIME = 0.D0
      WRITE(10,999)TIME,N(1),N(2),N(3),N(4)
      WRITE(11,999)TIME,MX(1),MX(2),MX(3),MX(4)
      WRITE(12,999)TIME,MY(1),MY(2),MY(3),MY(4)
      WRITE(13,999)TIME,MZ(1),MZ(2),MZ(3),MZ(4)
      WRITE(14,999)TIME,MM(1),MM(2),MM(3),MM(4)
      WRITE(15,999)TIME,MTOT(1),MTOT(2),MTOT(3),
     &             DSQRT(MTOT(1)**2+MTOT(2)**2+MTOT(3)**2)
C**---------------------------------------------------------------------
C**   Now start the time propagation with a sudden switch
C**   (read in the new potential and magnetic field)
C**---------------------------------------------------------------------
      READ(1,*)
      DO I=1,NP
         READ(1,*)V(I),BX(I),BY(I),BZ(I)
      ENDDO   
      
100   CONTINUE
      TIME = TIME + DT

C      DO 16 I=1,NP
C         BXD(I) = BX(I) + 0.01D0*DSIN(OMEGA*(TIME-DT/2.D0))
C16    CONTINUE
C      CALL MATRIX(M,U0,U1,T,V,BXD,BY,BZ)

      CALL MATRIX(M,U0,U1,T,V,BX,BY,BZ)

      DO I=1,NB
         DO J=1,NB
            A(I,J) = -0.5D0*IONE*DT*M(I,J)
            IF (I.EQ.J) A(I,J) = ONE + A(I,J)
         ENDDO
      ENDDO

      DO I=1,NB
         R(I) = ZERO
         DO J=1,NB
            R(I) = R(I) + A(I,J)*C(J)
         ENDDO
      ENDDO   
       
      DO I=1,NB
         DO J=1,NB
            A(I,J) = 0.5D0*IONE*DT*M(I,J)
            IF (I.EQ.J) A(I,J) = ONE + A(I,J)
         ENDDO
      ENDDO

      CALL ZGESV( NB, 1, A, NB, IPIV, R, NB, INFO )
      DO I=1,NB
         C(I) = R(I)
      ENDDO
      CALL DENCALC(C,N,MX,MY,MZ,MM,MTOT)

      WRITE(10,999)TIME,N(1),N(2),N(3),N(4)
      WRITE(11,999)TIME,MX(1),MX(2),MX(3),MX(4)
      WRITE(12,999)TIME,MY(1),MY(2),MY(3),MY(4)
      WRITE(13,999)TIME,MZ(1),MZ(2),MZ(3),MZ(4)
      WRITE(14,999)TIME,MM(1),MM(2),MM(3),MM(4)
      WRITE(15,999)TIME,MTOT(1),MTOT(2),MTOT(3),
     &             DSQRT(MTOT(1)**2+MTOT(2)**2+MTOT(3)**2)
      IF (TIME.LT.199.9999D0) GOTO 100
      
997   FORMAT(F8.3,3F10.5)
998   FORMAT(F8.3,5F10.5)
999   FORMAT(F8.3,4F10.5)
      END
C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE MATRIX(M,U0,U1,T,V,BX,BY,BZ)
      IMPLICIT NONE

      INTEGER NP,NB,NT,NS,I,J,K,L,II,JJ      
      INCLUDE 'dim.inc'
      PARAMETER (NB=2*NP**2-NP)
      DOUBLE PRECISION U0,U1,T
      DOUBLE PRECISION V(NP),BX(NP),BY(NP),BZ(NP)

      DOUBLE COMPLEX M(NB,NB),BP(NP),BM(NP)
      DOUBLE COMPLEX ZERO,ONE,IONE
      PARAMETER (ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),IONE=(0.D0,1.D0))

      DO I=1,NB
         DO J=1,NB
            M(I,J) = ZERO
         ENDDO
      ENDDO

      DO I=1,NP
         BP(I) = (BX(I)*ONE + BY(I)*IONE)/DSQRT(2.D0)
         BM(I) = (BX(I)*ONE - BY(I)*IONE)/DSQRT(2.D0)
      ENDDO   

      NS = NP*(NP+1)/2
      NT = 3*NP*(NP-1)/2

C**------------------------------------------------------------
C**   Singlet block
C**------------------------------------------------------------

C**   First the diagonal elements

      II = -NP
      DO 10 I=1,NP  
         II = II + NP+2-I
         M(II,II) = 2.D0*V(I) + U0
10    CONTINUE

      II = -NP+1
      DO 11 I=1,NP-1    
         II = II + NP+2-I
         M(II,II) = V(I) + V(I+1) + U1
11    CONTINUE

      IF (NP.GT.2) THEN
         JJ=-NP+2
         DO 12 J=1,NP-2
            JJ = JJ + NP+2-J
            DO 13 I=1,NP-J-1
               M(JJ+I-1,JJ+I-1) = V(J) + V(J+I+1) 
13          CONTINUE

12       CONTINUE
      ENDIF

C**   Now the off-diagonal elements

      II = -NP
      DO 15 I=1,NP-1
         II = II + NP+2-I
         M(II,II+1) = -DSQRT(2.D0)*T
         M(II+1,II) = -DSQRT(2.D0)*T
15    CONTINUE
      
      II = -NP+1
      DO 16 I=1,NP-1    
         II = II + NP+2-I
         M(II,II+NP-I) = -DSQRT(2.D0)*T    
         M(II+NP-I,II) = -DSQRT(2.D0)*T   
16    CONTINUE
      
      IF (NP.GT.2) THEN
         JJ=-NP+2
         DO 17 J=1,NP-2
            JJ = JJ + NP+2-J
            DO 18 I=1,NP-J-1
               M(JJ+I-2,JJ+I-1) = -T   
               M(JJ+I-1,JJ+I-2) = -T   

               M(JJ+I-1,JJ+I-1+(NP-J)) = -T                  
               M(JJ+I-1+(NP-J),JJ+I-1) = -T
18          CONTINUE

17       CONTINUE
      ENDIF
C**------------------------------------------------------------
C**   Singlet-Triplet and Triplet-Singlet blocks
C**------------------------------------------------------------

      II = -NP+1
      JJ = NS+1
      DO 20 J=1,NP-1
         II = II + NP-J+2 
         DO 21 I=1,NP-J
             
            M(II+I-1,JJ) = -BP(J) + BP(J+I)
            M(II+I-1,JJ+1) = BZ(J) - BZ(J+I)
            M(II+I-1,JJ+2) = BM(J) - BM(J+I)
            
            M(JJ,II+I-1) = -BM(J) + BM(J+I)
            M(JJ+1,II+I-1) = BZ(J) - BZ(J+I)
            M(JJ+2,II+I-1) = BP(J) - BP(J+I)
            
            JJ=JJ+3  
21       CONTINUE          
20    CONTINUE          

C**------------------------------------------------------------
C**   Triplet block
C**------------------------------------------------------------

      II = NS-1
      DO 40 I=1,NP-1
         DO 41 J=I+1,NP
            II = II + 3
            M(II-1,II-1) = V(I) + V(J) + BZ(I) + BZ(J)
            M(II,II) =  V(I) + V(J)
            M(II+1,II+1) = V(I) + V(J) - BZ(I) - BZ(J)
            
            IF ((J-I).EQ.1) THEN
               M(II-1,II-1) = M(II-1,II-1) + U1
               M(II,II) = M(II,II) + U1
               M(II+1,II+1) = M(II+1,II+1) + U1
            ENDIF
            
            M(II-1,II) = BM(I) + BM(J)
            M(II,II-1) = BP(I) + BP(J)
            M(II,II+1) = BM(I) + BM(J)
            M(II+1,II) = BP(I) + BP(J)
                        
41       CONTINUE          
40    CONTINUE
        
      II = NS-2
      DO 45 I=1,NP-1
         DO 46 J=I+1,NP
            II = II+3
            
            JJ = NS-2
            DO 47 K=1,NP-1
               DO 48 L=K+1,NP
                  JJ = JJ+3
                  IF ((I.EQ.K).AND.(IABS(J-L).EQ.1).OR.
     &                (J.EQ.L).AND.(IABS(I-K).EQ.1)) THEN
                     M(II,JJ) = -T  
                     M(II+1,JJ+1) = -T
                     M(II+2,JJ+2) = -T
                  ENDIF        
48             CONTINUE
47             CONTINUE                   
             
46       CONTINUE
45    CONTINUE             

      RETURN
      END
C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE DENCALC(C,N,MX,MY,MZ,MM,MTOT)
      IMPLICIT NONE
      
      INTEGER NP,NB,NT,NS,I,J,K,II,JJ,KK
      INCLUDE 'dim.inc'
      PARAMETER (NB=2*NP**2-NP)
      DOUBLE COMPLEX C(NB),NUD(NP)
      DOUBLE PRECISION N(NP),MX(NP),MY(NP),MZ(NP),MM(NP),PHI(NP,NP),
     &                 MTOT(3)
           
      DO I=1,NP
         DO J=1,NP
            PHI(I,J) = 0.D0
            IF (I.EQ.J) PHI(I,J) = 1.D0
         ENDDO  
      ENDDO
      
      NS = NP*(NP+1)/2
      NT = 3*NP*(NP-1)/2

      DO 100 I=1,NP
          
C**------------------------------------------------------------
C**   Singlet block contributions to N
C**------------------------------------------------------------
         
         N(I) = 0.D0
         
         II = -NP
         DO 10 J=1,NP
            II = II + NP+2-J
            N(I) = N(I) + PHI(J,I)*CDABS(C(II))**2 
10       CONTINUE  
         
         II = -NP
         DO 11 J=1,NP-1
            II = II + NP+2-J
            DO 12 K=1,NP-J
               JJ = II + K
               N(I) = N(I)+0.5D0*(PHI(J,I)+PHI(J+K,I))*CDABS(C(JJ))**2 
12          CONTINUE                
11       CONTINUE
         
C**------------------------------------------------------------
C**   Triplet block contributions to N
C**------------------------------------------------------------         
         
         JJ = NS-2
         DO 13 J=1,NP-1
            DO 14 K=1,NP-J
               JJ = JJ+3
               N(I) = N(I)+0.5D0*(PHI(J,I)+PHI(J+K,I))*(CDABS(C(JJ))**2 
     &                    + CDABS(C(JJ+1))**2 + CDABS(C(JJ+2))**2 )
14          CONTINUE                
13       CONTINUE
         
C**------------------------------------------------------------
C**   Triplet block contributions to MZ
C**------------------------------------------------------------ 
         
         MZ(I) = 0.D0
         
         JJ = NS-2
         DO 15 J=1,NP-1
            DO 16 K=1,NP-J
               JJ = JJ+3
               MZ(I) = MZ(I)+0.5D0*(PHI(J,I)+PHI(J+K,I))
     &                      *(CDABS(C(JJ))**2 - CDABS(C(JJ+2))**2 )
16          CONTINUE                
15       CONTINUE
         
C**------------------------------------------------------------
C**   Singlet-Triplet block contributions to MZ
C**------------------------------------------------------------         
               
         II = -NP
         KK = NS-1
         DO 20 J=1,NP-1
            II = II + NP+2-J
            DO 21 K=1,NP-J
               JJ = II+K
               KK = KK+3
               MZ(I) = MZ(I)+0.5D0*(PHI(J,I)-PHI(J+K,I))
     &           *( C(JJ)*DCONJG(C(KK)) + DCONJG(C(JJ))*C(KK) )        
21          CONTINUE                
20       CONTINUE

C**------------------------------------------------------------
C**   Singlet-Triplet block contributions to NUD
C**------------------------------------------------------------

         NUD(I) = (0.D0,0.D0)
         
         II = -NP
         KK = NS-1
         DO 25 J=1,NP-1
            II = II + NP+2-J
            DO 26 K=1,NP-J
               JJ = II+K
               KK = KK+3
               NUD(I) = NUD(I)+0.25D0*DSQRT(2.D0)*(PHI(J,I)-PHI(J+K,I))
     &           *( C(JJ)*DCONJG(C(KK+1)) - DCONJG(C(JJ))*C(KK-1) )        
26          CONTINUE                
25       CONTINUE
         
C**------------------------------------------------------------
C**   Triplet block contributions to NUD
C**------------------------------------------------------------

         KK = NS-1
         DO 27 J=1,NP-1
            II = II + NP+2-J
            DO 28 K=1,NP-J
               KK = KK+3
               NUD(I) = NUD(I)+0.25D0*DSQRT(2.D0)*(PHI(J,I)+PHI(J+K,I))
     &           *( C(KK)*DCONJG(C(KK+1)) + DCONJG(C(KK))*C(KK-1) )        
28          CONTINUE                
27       CONTINUE

         MX(I) = 2.D0*DREAL(NUD(I)) 
         MY(I) = -2.D0*DIMAG(NUD(I)) 
      
         N(I) = 2.D0*N(I)
         MX(I) = 2.D0*MX(I)
         MY(I) = 2.D0*MY(I)
         MZ(I) = 2.D0*MZ(I) 
         MM(I) = DSQRT(MX(I)**2+MY(I)**2+MZ(I)**2) 
100   CONTINUE         
    
      MTOT(1) = 0.D0
      MTOT(2) = 0.D0
      MTOT(3) = 0.D0
      DO I=1,NP
         MTOT(1) = MTOT(1) + MX(I)
         MTOT(2) = MTOT(2) + MY(I)
         MTOT(3) = MTOT(3) + MZ(I)
      ENDDO    
      
      RETURN
      END
