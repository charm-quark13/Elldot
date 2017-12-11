      Module POTGEN2

      IMPLICIT NONE

      contains


      SUBROUTINE POTENTIAL(V,GRID)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: GRID

      DOUBLE PRECISION, INTENT(OUT) :: V(GRID)

      INTEGER :: I

      DO I=1,GRID
         V(I) = 0.D0

         IF (I.GE.1.AND.I.LE.40) THEN
            V(I) = -6.D0

         !ELSEIF (I.GE.56.AND.I.LE.75) THEN
            !V(I) = 2.0D0

         !ELSEIF (I.GE.81.AND.I.LE.100) THEN
            !V(I) = -2.0D0

         ELSEIF (I.GE.41.AND.I.LE.154) THEN
            V(I) = -3.D0

         ELSEIF (I.GE.155.AND.I.LE.200) THEN
            V(I) = -6.2D0

         ELSE
            V(I) = 0.D0

         ENDIF

      END DO

      END SUBROUTINE

      END MODULE

