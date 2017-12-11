      PROGRAM MATRIX
      USE ROGEN
      use potgen3
      IMPLICIT NONE

      INTEGER N,nocc,LWORK,LWORKC,LC,INFO,I,IP,J,JP,K,L,A,AP,R,RP
      INTEGER JJ,NEX,PERT,M,LCWORK
C      PARAMETER (N=200, NOCC=8, NUN=N-NOCC, LC=NOCC*NUN, 
      PARAMETER (N=200, nocc=7,LC=2*NOCC*NUN,LCWORK=10*LC,
     &           LWORK=10*N, LWORKC=10*LC, PERT=NOCC*(NOCC-1)/2 )
      DOUBLE PRECISION HMAT(N,N),W(N),WORK(LWORK),WC(LC),
     &                 WR(LC),WI(LC),WORKC(LWORKC),CWORK(LCWORK)
      CHARACTER*1 JOBZ, UPLO

      DOUBLE PRECISION DELTA,RR,V(N),RHO1(N),DUM1,DUM2
      DOUBLE PRECISION OMEGA_AI(LC),CMAT(LC,LC),KER(N,N)
      DOUBLe PRECISION N1(N),PHM(N,N),TDM(N,N)
      DOUBLE PRECISION TEST(NOCC,NOCC),IDENTITY(NOCC,NOCC)
      DOUBLE PRECISION ANGLES(PERT),x,DUM
      DOUBLE PRECISION RMAT(NOCC,NOCC),AMAT(LC/2,LC/2),
     &                   KMAT(LC/2,LC/2),VL(LC,LC),VR(LC,LC)

      DOUBLE PRECISION ujb(lc/2,lc),vjb(lc/2,lc),xi(n,nocc)
      DOUBLE PRECISION rmatdagger(nocc,nocc),dummy(n,n)

      CHARACTER(LEN=30) :: matrixformt,pformat

      write(matrixformt,'(a,i3,a)') '(', NOCC ,'E12.3)'

      write(pformat,'(a,i0, a)') '(', 1,'E13.4)'

      JOBZ = 'V'
      UPLO = 'U'

C*********************************************************************
C**   Define the grid spacing DELTA and set the potential
C*********************************************************************
      DELTA = 10.D0/(N-1)

      call potential(v,n)

      DO I=1,N
        WRITE(1,*)I,V(I)
      END DO

C*********************************************************************
C**   Set up the Hamiltonian matrix and solve the Schrodinger equation.
C**   The eigenstates will be contained in HMAT.
C*********************************************************************
      DO 5 I=1,N
      DO 5 J=1,N
         dummy(i,j)=0.d0
5        HMAT(I,J)=0.D0

      DO 6 I=1,N
6        HMAT(I,I) = 1.D0/DELTA**2 + V(I)

      DO 7 I=1,N-1
         HMAT(I+1,I) = -0.5D0/DELTA**2
7        HMAT(I,I+1) = -0.5D0/DELTA**2

      CALL DSYEV( JOBZ, UPLO, N, HMAT, N, W, WORK, LWORK, INFO )
      
      DO 9 I=1,7
      DO 9 J=1,N
9        WRITE(9+I,*)J,HMAT(J,I)

      do i=1,n
         do j=1,n
            dummy(i,j) = hmat(i,j)
         end do
      end do

      DO J=1,N
         WRITE(1000,matrixformt) HMAT(J,:)
      END DO

!      DO I=1,PERT
!         ANGLES(I) = 0.D0
!      END DO

      DO I=1,PERT
         ANGLES(I) = DBLE(I)/20.D0
      END DO

      DO I=1,NOCC
         DO J=1,NOCC
            IF (I.EQ.J) THEN
               TEST(J,I)=1.D0
               IDENTITY(J,I)=1.D0
            ELSE
               TEST(J,I)=0.D0
               IDENTITY(J,I)=0.D0
            END IF
         END DO
      END DO

      A = 0

      OPEN(999,file='occupied.txt')

!      DO I=1,PERT
!         READ(999,*) ANGLES(i)
!      END DO

      close(999)

      write(*,*) angles

!      DO I=1,NOCC-1
!         DO J=I+1,NOCC

!            TEST=IDENTITY
!            A = A + 1

!               TEST(I,I)= DCOS(ANGLES(A))
!               TEST(J,J)= DCOS(ANGLES(A))
!               TEST(J,I)= DSIN(ANGLES(A))
!               TEST(I,J)= -DSIN(ANGLES(A))

!               DO M=1,Nocc
!                  DO L=1,NOCC
!                     X=0.D0
!                     DO K=1,N
!                        X = X + HMAT(K,M)*TEST(L,M)
!                     END DO
!                     HMAT(K,L)=X
!                  END DO
!               END DO

!         END DO
!      END DO

      DO J=1,N
         DO I=1,NOCC
            WRITE(99+I,*)J, HMAT(J,I)
         END DO
      END DO

C*********************************************************************
C**   Set up the Casida equation. The index i of occupied orbitals
C**   runs from i=1,...,NOCC, and the index a runs from a=NOCC+1,...,N.
C**   
C**   We use the RPA.
C*********************************************************************
      DO 10 I=1,NOCC
      DO 10 J=1,NUN
         A = NOCC+J
         OMEGA_AI(J + NUN*(I-1)) = W(A) - W(I)
         WRITE(79,*)J + NUN*(I-1),OMEGA_AI(J + NUN*(I-1))
10    CONTINUE

      DO 20 I=1,NOCC
         WRITE(*,*)I,'   of',NOCC
      DO 20 J=1,NUN
         A = NOCC+J
         DO 21 IP=1,NOCC
         DO 21 JP=1,NUN
            AP = NOCC+JP

            DO 25 L=1,N
               DUM1 = HMAT(L,I)*HMAT(L,A)/DELTA
            DO 25 K=1,N
               IF (L.EQ.K) THEN
                  KER(L,K) = 0.D0
               ELSE
                  KER(L,K) = DUM1*HMAT(K,IP)*HMAT(K,AP)/ABS(L-K)
               ENDIF
25          CONTINUE

            DUM1 = 0.D0
            DO 30 L=1,N
               DUM2 = 0.D0
               DO 31 K=1,N
                  DUM2 = DUM2 + KER(L,K)
31             CONTINUE
               DUM1 = DUM1 + DUM2
30          CONTINUE

      KMAT(J+NUN*(I-1),JP+NUN*(IP-1)) = 2.D0*DUM1*DELTA**2

C*********************************************************************
C**    The factor 2.D0 is because we need A+B later
C*********************************************************************
21    CONTINUE
20    CONTINUE

      DO 40 I=1,NOCC
      DO 40 J=1,NUN
         A = NOCC+J
         DO 41 IP=1,NOCC
         DO 41 JP=1,NUN
            AP = NOCC+JP

            IF (I.EQ.IP.AND.A.EQ.AP) THEN
               KMAT(J+NUN*(I-1),JP+NUN*(IP-1)) =
     &         KMAT(J+NUN*(I-1),JP+NUN*(IP-1)) + OMEGA_AI(J+NUN*(I-1))
            ENDIF

            KMAT(J+NUN*(I-1),JP+NUN*(IP-1)) =
     &           DSQRT(OMEGA_AI(J+NUN*(I-1))) *
     &           KMAT(J+NUN*(I-1),JP+NUN*(IP-1)) *
     &           DSQRT(OMEGA_AI(JP+NUN*(IP-1)))

41       CONTINUE
40    CONTINUE

C*********************************************************************
C**   Diagonalize Casida equation
C*********************************************************************

      CALL DSYEV( JOBZ, UPLO, LC, KMAT, LC, WC, WORKC, LWORKC, INFO )

      DO 50 I=1,LC
         WRITE(31,*)I,DSQRT(WC(I))
50    CONTINUE

      DO 59 I=1,5
      DO 59 J=1,LC
59       WRITE(20+I,*)J,KMAT(J,I)


      DO 54 I=1,LC
54       WRITE(31,*)I-NOCC*NUN,WR(I)!,WI(I)

!      DO 59 I=1,5   
!      DO 59 J=1,LC
!59       WRITE(200+I,*)J,KMAT(J,I)

!^^^^^^^^^^^^^^^^^^      Correct to this point   ^^^^^^^^^^^^^^^^^^^^^

      call romat(nocc,pert,angles,rmat)

      do i=1,nocc
         do j=1,nocc
            rmatdagger(i,j) = rmat(j,i)
         end do
      end do

      do i=1,n
         do j=1,nocc
            xi(i,j) = 0.d0
         end do
      end do

      DO I=1,N
         DO J=1,NOCC
            DO K=1,NOCC
               xi(i,j) = xi(i,j) + rmatdagger(j,k)*hmat(i,k)
            end do
         end do
      end do

      do i=1,n
         do j=1,nocc
            hmat(i,j) = xi(i,j)
         end do
      end do

      DO J=1,N
         DO I=1,NOCC
            WRITE(99+I,*)J, HMAT(J,I)
         END DO
      END DO

C*********************************************************************
C**   Calculate the density response of the NEX-th excitation
C*********************************************************************
      DO 100 NEX=1,5

!      write(*,matrixformt) rmat

      do i=1,nocc
         do j=1,nun
            x = 0.d0
            do k=1,nocc
               x = x + rmat(i,k)*kmat(j+nun*(k-1),nex)
            end do
            ujb(i+nocc*(j-1),nex) = x
         end do
      end do

      do j=1,lc/2
         write(250+nex,*) ujb(j,nex)
      end do

      do i=1,nocc
         do j=1,nun
            x = 0.d0
            do k=1,nocc
               x = x + rmat(i,k)*kmat(lc/2+j+nun*(k-1),nex)/2.d0
            end do
            vjb(i+nocc*(j-1),nex) = x
         end do
      end do

      do j=1,lc/2
         write(350+nex,*) vjb(j,nex)
      end do
 
      DO 60 K=1,N
         N1(K) = 0.D0
         DO 61 I=1,NOCC
         DO 61 J=1,NUN        
            A = NOCC+J
            N1(K) = N1(K) + xi(K,I)*hmat(K,A)
!     &              *(VR(J+NUN*(I-1),NEX)+VR(LC/2+J+NUN*(I-1),NEX))
     &               *(ujb(j+nun*(i-1),nex)+vjb(j+nun*(i-1),nex))
     &               *dsqrt(OMEGA_AI(J + NUN*(I-1)))

61       CONTINUE
         WRITE(41+10*NEX,*)K,N1(K)
60    CONTINUE
C*********************************************************************
C**   Calculate the TDM and the PHM of the n-th excitation
C*********************************************************************

      DO I=1,N
         DO J=1,N
            PHM(I,J) = 0.D0
            TDM(I,J) = 0.D0
         END DO
      END DO


!      DO NEX=1,5
         DO R=1,N
         DO RP=1,N
            DO I=1,NOCC
            DO J=1,NUN
               A=NOCC+J

!               TDM(R,RP) = TDM(R,RP) + HMAT(R,I)*HMAT(RP,A) *
!     &                        vr(J+NUN*(I-1),NEX) +
!     &              HMAT(RP,I) * HMAT(R,A) * vr(LC/2+J+NUN*(I-1),NEX)

            END DO
            END DO
         END DO
         END DO
!      END DO
       DO R = 1,N
          DO RP=1,N
             write(92,*) RP,R,TDM(R,RP)
          END DO
       END DO


C      DO 75 K=1,N
C      DO 76 L=1,N
C76       WRITE(42+10*NEX,*)PHM(L,K),TDM(L,K)
C75       WRITE(42+10*NEX,900)

      DO 80 K=1,N
         WRITE(460+10*NEX,901) (PHM(K,L), L=1,N,1)
         WRITE(570+10*NEX,901) (TDM(K,L), L=1,N,1)
80    CONTINUE

100   CONTINUE

900   FORMAT()
901   FORMAT(200E12.4)
      END


