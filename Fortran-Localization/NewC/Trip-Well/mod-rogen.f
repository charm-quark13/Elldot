      Module ROGEN

      IMPLICIT NONE

      contains


      SUBROUTINE ROMAT(OCC,PERT,XIE,rmat)

      IMPLICIT NONE

      INTEGER :: OCC

      INTEGER :: I,N,J,K,COUNTER,M,L,PERT

      DOUBLE PRECISION :: X

      DOUBLE PRECISION RMAT(OCC,OCC)

      DOUBLE PRECISION V(OCC,OCC),IDEN(OCC,OCC),TEST(OCC,OCC)

      DOUBLE PRECISION XIE(PERT),VECTOR(PERT)

      DO J=1,OCC
         DO K=1,OCC
            V(K,J)=0.D0
            IDEN(K,J)=0.D0
            TEST(K,J)=0.D0
         END DO
      END DO

      DO K=1,OCC
         V(K,K)=1.D0
         IDEN(K,K)=1.D0
      END DO

      RMAT=IDEN

      test=0.d0
      COUNTER = 0
      X = 0.D0
               DO I=1,OCC-1
                  DO J=I+1,OCC
                     V=IDEN

                     COUNTER = COUNTER + 1

                     V(I,I)= DCOS(XIE(COUNTER))
                     V(J,J)= DCOS(XIE(COUNTER))
                     V(J,I)= DSIN(XIE(COUNTER))
                     V(I,J)= -DSIN(XIE(COUNTER))

                     !write(*,*) v
                     !write(*,*) '   '
                     !write(*,*) rmat

                     DO K=1,OCC
                        DO L=1,OCC
                           x=0.d0
                           DO M=1,OCC
                              x = x + V(M,L)*rmat(K,M)
                           END DO
                           test(K,L) = X
                        END DO
                     END DO
                  rmat=test
                  END DO
               END DO

         rmat = test

      END SUBROUTINE

      END MODULE

