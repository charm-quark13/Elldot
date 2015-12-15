      PROGRAM MATTEST

      DOUBLE COMPLEX M(4,4), N(4,4), MATRIX(4,4)

      DOUBLE COMPLEX IONE,ONE,ZERO
      INTEGER I,J,K,L

      PARAMETER(IONE = (0.D0,1.D0),ONE=(1.D0,0.D0),ZERO=(0.D0,0.D0))

      DO 10 I=1,4
      DO 10 K=1,4
         M(I,K) = ZERO
         N(I,K) = ZERO
         MATRIX(I,K) = ZERO
10    CONTINUE

      DO 20 I=1,4

         M(I,I) = IONE

20    CONTINUE

      N(1,1) = (2.D0,1.D0)
      N(1,4) = (3.D0,0.D0)
      N(2,2) = (2.D0,1.D0)
      N(3,3) = IONE
      N(4,1) = ONE
      N(4,4) = (1.D0,1.D0)


      DO 30 I=1,4
      DO 30 K=1,4

         MATRIX(I,K) = M(I,K)*CONJG(M(I,K))

30    CONTINUE

        WRITE(*,*) MATRIX

        WRITE(*,*)'____________________________________________________'

        WRITE(*,*) REAL(MATRIX(1,1)),REAL(MATRIX(2,2)),
     &              REAL(MATRIX(3,3)),REAL(MATRIX(4,4))
      DO 40 I=1,4
      DO 40 K=1,4

         MATRIX(I,K) = ZERO

40    CONTINUE

      DO 50 I=1,4
      DO 50 K=1,4

         MATRIX(I,K) = N(I,K)*CONJG(N(I,K))

50    CONTINUE 

         WRITE(*,*) '__________________________________________________'
         WRITE(*,*) MATRIX

      END
