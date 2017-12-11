      PROGRAM MATRIX
      USE FBMethod
      !use boys
      USE ROGEN
      use potgen3
!      use potgen2
      IMPLICIT NONE

      INTEGER N,LWORK,LWORKC,LC,INFO,I,IP,J,JP,K,L,A,AP,R,RP
      INTEGER JJ,NEX,M,LCWORK,iter,itmax,counter
C      PARAMETER (N=200, NOCC=8, NUN=N-NOCC, LC=NOCC*NUN, 
      PARAMETER (N=200, LC=2*NOCC*NUN,LCWORK=10*LC,
     &           LWORK=10*N, LWORKC=10*LC, itmax=250000 )
      DOUBLE PRECISION HMAT(N,N),W(N),WORK(LWORK),WC(LC),
     &                 WR(LC),WI(LC),WORKC(LWORKC),CWORK(LCWORK)
      CHARACTER*1 JOBZ, UPLO

      double precision,PARAMETER :: thresh=2.5d-3
      DOUBLE PRECISION RR,V(N),RHO1(N),DUM1,DUM2
      DOUBLE PRECISION OMEGA_AI(LC),CMAT(LC,LC),KER(N,N)
      DOUBLE PRECISION N1(N),PHM(N,N),TDM(N,N)
      DOUBLE PRECISION d(nocc,nocc)!,dboys(nocc,nocc),d(nocc,nocc)
      DOUBLE PRECISION TEST(NOCC,NOCC),IDENTITY(NOCC,NOCC)
      DOUBLE PRECISION ANGLES(PERT),x,DUM,HRtMat(n,nocc)
      DOUBLE PRECISION RMAT(NOCC,NOCC),AMAT(LC/2,LC/2),
     &                   KMAT(LC/2,LC/2),VL(LC,LC),VR(LC,LC)

      DOUBLE PRECISION ujb(lc/2,lc),vjb(lc/2,lc),xi(n,nocc)
      DOUBLE PRECISION rmatdagger(nocc,nocc),dummy(lc/2,lc/2)
      DOUBLE PRECISION zvec(lc/2),jacobiev(nocc)

      CHARACTER(LEN=30) :: matrixfmt,pformat,kform1,kform2

      write(matrixfmt,'(a,i3,a)') '(', NOCC ,'E12.3)'
      write(kform1,'(a,i3,a)') '(', lc/2 ,'E12.3)'
      write(kform2,'(a,i3,a)') '(', lc ,'E12.3)'
      write(pformat,'(a,i0, a)') '(', 1,'E13.4)'

      JOBZ = 'V'
      UPLO = 'U'

C*********************************************************************
C**   Define the grid spacing DELTA and set the potential
C*********************************************************************

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
5        HMAT(I,J)=0.D0

      DO 6 I=1,N
6        HMAT(I,I) = 1.D0/DELTA**2 + V(I)

      DO 7 I=1,N-1
         HMAT(I+1,I) = -0.5D0/DELTA**2
7        HMAT(I,I+1) = -0.5D0/DELTA**2

      CALL DSYEV( JOBZ, UPLO, N, HMAT, N, W, WORK, LWORK, INFO )
      
      DO 9 I=1,nocc
      DO 9 J=1,N
9        WRITE(9+I,*)J,HMAT(J,I)

      DO I=1,PERT
         ANGLES(I) = 0.d0
      END DO

      call romat(nocc,pert,angles,rmat)

!      write(*,matrixfmt) rmat

      do i=1,nocc
         do j=1,nocc
            rmatdagger(i,j) = rmat(j,i)
         end do
      end do

!      write(*,*) '********'

      do i=1,n
         do j=1,nocc
            xi(i,j) = hmat(i,j)
         end do
      end do

      hrtmat = matmul(xi,rmat)

      do i=1,nocc
!         x = 0.d0
         do j=1,n
            write(999+i,*)j,hrtmat(j,i)
!            x = x + hrtmat(j,i)**2
         end do
!         write(*,*) '*************'
!         write(*,*) x
      end do

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

      DO J=1,lc/2
         DO I=1,lc/2
            amat(j,i) = 0.d0
            kmat(j,i) = 0.d0
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
         WRITE(78,*)J + NUN*(I-1),OMEGA_AI(J + NUN*(I-1))
10    CONTINUE

      DO 20 I=1,NOCC
         WRITE(*,*)I,'   of',NOCC
      DO 20 J=1,NUN        
         A = NOCC+J
         DO 21 IP=1,NOCC
         DO 21 JP=1,NUN
            AP = NOCC+JP
         
            DO 25 L=1,N
               DUM1 = hmat(l,i)*HMAT(L,A)/DELTA
            DO 25 K=1,N
               IF (L.EQ.K) THEN 
                  KER(L,K) = 0.D0
               ELSE
                  KER(L,K) = DUM1*hmat(k,IP)*HMAT(K,AP)/ABS(L-K)
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

      KMAT(J+NUN*(I-1),JP+NUN*(IP-1)) = DUM1*DELTA**2

21    CONTINUE
20    CONTINUE

      DO I=1,NOCC
      DO J=1,NUN
         A=NOCC+J
         DO IP=1,NOCC
         DO JP=1,NUN
            AP=NOCC+JP

            IF (I.EQ.IP .AND. A.EQ.AP) THEN
               AMAT(J+NUN*(I-1),JP+NUN*(IP-1)) =
     &            KMAT(J+NUN*(I-1),JP+NUN*(IP-1)) +
     &               OMEGA_AI(J + NUN*(I-1))
            ELSE
               AMAT(J+NUN*(I-1),JP+NUN*(IP-1)) =
     &            KMAT(J+NUN*(I-1),JP+NUN*(IP-1))
!     &            0.d0
            END IF

         END DO
         END DO
      END DO
      END DO

      DO I=1,LC/2
         DO J=1,LC/2
            CMAT(I,J)=-AMAT(I,J)
            CMAT(LC/2+I,J)= KMAT(I,J)!*0.d0
            CMAT(I,LC/2+J)= -KMAT(I,J)!*0.d0
            CMAT(LC/2+I,LC/2+J)= AMAT(I,J)
         END DO
      END DO
    
C*********************************************************************
C**   Diagonalize Casida equation
C*********************************************************************
      CALL DGEEV( 'N', 'V', LC, CMAT, LC, WR, WI, VL, LC, VR,
     $                  LC, CWORK, LCWORK, INFO )

*********************************************************************
C     Sort the eigenvalues in ascending order
C*********************************************************************
      DO 51 I=1,LC-1
         DO 52 J=I+1,LC
            IF (WR(I).GE.WR(J)) THEN
               DUM = WR(I)
               WR(I) = WR(J)
               WR(J) = DUM

               DUM = WI(I)
               WI(I) = WI(J)
               WI(J) = DUM

               DO 53 JJ=1,LC
                  DUM = VR(JJ,I)
                  VR(JJ,I) = VR(JJ,J)
                  VR(JJ,J) = DUM

                  DUM = VL(JJ,I)
                  VL(JJ,I) = VL(JJ,J)
                  VL(JJ,J) = DUM
53             CONTINUE

            ENDIF
52       CONTINUE
51    CONTINUE

      DO 54 I=1,LC
54       WRITE(32,*)I-NOCC*NUN,WR(I)!,WI(I)

      do i=lc/2+1,lc/2+5
         x = 0.d0
         do j=1,lc/2
            x = x + (vr(j+lc/2,i)**2 - vr(j,i)**2)
         end do
        
         do j=1,lc
            vr(j,i) = vr(j,i)/dsqrt(x)
         end do        
!         write(*,*) x
      end do

      counter = 0
      x = 1.d0
      do iter=1,itmax

         counter = counter + 1

         call boys(hrtmat,d)

         call jacobi(d,nocc,nocc,jacobiev,rmatdagger,info)
   
         if (counter.eq.1) then
            do i=1,n
               do j=1,nocc
                  xi(i,j) = 0.d0
               end do
            end do
 
            do i=1,nocc
               do j=1,nocc
                  if (i.eq.j) then
                     rmat(i,j) = 1.d0
                  else
                     rmat(i,j) = 0.d0
                  endif
               end do
            end do

            xi=matmul(hrtmat,rmatdagger)

            test=matmul(rmat,rmatdagger)

         else

            xi=matmul(hrtmat,rmatdagger)

            rmat=matmul(rmatdagger,test)

         end if   

!         do i=1,n
!            do j=1,nocc
!               hrtmat(i,j) = xi(i,j)
!            end do
!         end do

         do i=1,nocc
            d(i,i) = 0.d0
         end do

         x = maxval(d)

         if (dabs(x).lt.thresh) exit
          
         do i=1,nocc
            do j=1,nocc
               test(i,j) = rmat(i,j)
            end do
         end do

      end do   

      do i=1,nocc
         do j=1,nocc
            rmatdagger(i,j) = rmat(j,i)
         end do
      end do

      test = matmul(rmat,rmatdagger)

      write(*,matrixfmt) test

      write(*,*) counter
!!      write(*,*) '**** rdagger *****'
!!      write(*,matrixfmt) rmatdagger  
!      write(*,*) '**** d matrix *****'
!      write(*,matrixfmt) d 

      xi = matmul(hrtmat,rmatdagger)

      DO J=1,N
         DO I=1,NOCC
            WRITE(99+I,*)J, xi(J,I)
         END DO
      END DO

C*********************************************************************
C**   Calculate the density response of the NEX-th excitation
C*********************************************************************
      DO 100 NEX=lc/2+1,lc/2+5

!      write(*,*) '******'
!      write(*,*) rmat(3,4)
!      write(*,matrixfmt) rmat

      do i=1,nocc
         do j=1,nun
            x = 0.d0
            do k=1,nocc
               x = x + rmatdagger(k,i)*vr(j+nun*(k-1),nex)
            end do
            ujb(j+nun*(i-1),nex) = x
!            ujb(i+nocc*(j-1),nex) = x
         end do
      end do

      do j=1,lc/2
         write(250+(nex-lc/2),*) ujb(j,nex)
      end do

      do i=1,nocc
         do j=1,nun
            x = 0.d0
            do k=1,nocc
               x = x + rmatdagger(k,i)*vr(lc/2+j+nun*(k-1),nex)
            end do
            vjb(j+nun*(i-1),nex) = x
!            vjb(i+nocc*(j-1),nex) = x
         end do
      end do

      do j=1,lc/2
         write(350+(nex-lc/2),*) vjb(j,nex)
      end do

      DO 60 K=1,N
         N1(K) = 0.D0
         DO 61 I=1,NOCC
         DO 61 J=1,NUN        
            A = NOCC+J
            N1(K) = N1(K) + xi(K,I)*HMAT(K,A)
!     &              *(VR(J+NUN*(I-1),NEX)+VR(LC/2+J+NUN*(I-1),NEX))
     &               *(ujb(j+nun*(i-1),nex)+vjb(j+nun*(i-1),nex))

61       CONTINUE
         WRITE(43+10*(NEX-lc/2),*)K,N1(K)
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


         DO R=1,N
         DO RP=1,N
            DO I=1,NOCC
            DO J=1,NUN
               A=NOCC+J

               TDM(R,RP) = TDM(R,RP) + HMAT(R,I)*HMAT(RP,A) *
     &                        vr(J+NUN*(I-1),NEX) +
     &              HMAT(RP,I) * HMAT(R,A) * vr(LC/2+J+NUN*(I-1),NEX)

               PHM(r,rp) = PHM(r,rp) + xi(r,i)**2 * 
     &                   (xi(rp,i)*hmat(rp,a) * ujb(j+nun*(i-1),nex) +
     &                   xi(rp,i) * hmat(rp,a) * vjb(j+nun*(i-1),nex)) 
               

            END DO
            END DO
         END DO
         END DO

C      DO 75 K=1,N
C      DO 76 L=1,N
C76       WRITE(42+10*NEX,*)PHM(L,K),TDM(L,K)
C75       WRITE(42+10*NEX,900)

      DO 80 K=1,N
         WRITE(461+10*(NEX-lc/2),901) (PHM(L,K), L=1,N,1)
         WRITE(460+10*(NEX-lc/2),901) (TDM(K,L), L=1,N,1)
80    CONTINUE

100   CONTINUE

900   FORMAT()
901   FORMAT(200E12.4)
      END
