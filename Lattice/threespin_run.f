      PROGRAM THREESPIN
      IMPLICIT NONE

      INTEGER I,J,LWORK,MODE,INFO,outer,it,time
      PARAMETER (LWORK = 500)

      DOUBLE PRECISION U0,C,CP,T,TP,V(3),BX(3),BY(3),BZ(3), dt
      DOUBLE PRECISION E(20),RWORK(100),N(3),MX(3),MY(3),MZ(3),ediff
      DOUBLE COMPLEX M(20,20),WORK(LWORK), psin(20), psinp1(20)
     &               , m0(20,20)

      ! T = 0.5D0
!      WRITE(*,*)'C?'
!      READ(*,*)C
!      WRITE(*,*)'Chain or triangle? (1,2)'
!      READ(*,*)MODE
!      IF (MODE.EQ.1) THEN
!         TP = 0.D0
!         CP = 0.D0
!      ELSEIF (MODE.EQ.2) THEN
!         TP = T
!         CP = C
!      ELSE
!         WRITE(*,*)'Invalid choice'
!         STOP
!      ENDIF
      dt = .01d0

      U0 = 1.d0

      do i=1,3
        v(i) = 0.d0
        bx(i) = 0.d0
        by(i) = 0.d0
        bz(i) = 0.d0
      end do

      ! bz(1) = 1.d-11

      open(26, file='full_energy_surfaces.txt')

      do outer = 1, 100

        T = dble(outer) * 0.01D0

        TP = T

        MODE = 2

!        DO 5 I=1,3
!           READ(1,*)V(I),BX(I),BY(I),BZ(I)
!5       CONTINUE

        DO 101 J=1,251
          C = dble(J - 1) * .2d0

          CP = C

          CALL MATRIX(M,U0,C,CP,T,TP,V,BX,BY,BZ)

          CALL ZHEEV( 'V', 'U', 20, M, 20, E, WORK, LWORK, RWORK, INFO )

          do it = 1, 20
            psin = m(it, 1)
          end do

C**----------------------------------------------------------------------
C**   calculate the densities:
C**----------------------------------------------------------------------
          CALL DENCALC(M,N,MX,MY,MZ)

    !       IF (outer.EQ.1.and.j.eq.1) THEN
    !         write(20,*)'# N'
    !         write(20,*)'# U0 ,   site 1   ,   site2   ,  site3  ',
    !  &                 ',  t/U0  ,  C/U0'

    !         write(21,*)'# MX'
    !         write(21,*)'# U0 ,   site 1   ,   site2   ,  site3  ',
    !  &                 ',  t/U0  ,  C/U0'

    !         write(22,*)'# MY'
    !         write(22,*)'# U0 ,   site 1   ,   site2   ,  site3  ',
    !  &                 ',  t/U0  ,  C/U0'

    !         write(23,*)'# MZ'
    !         write(23,*)'# U0 ,   site 1   ,   site2   ,  site3  ',
    !  &                 ',  t/U0  ,  C/U0'

    !       END IF

    !       WRITE(20,*)real(U0),real(N(1)),real(N(2)),real(N(3)), t/U0,
    !  &                 C/U0
    !       WRITE(21,*)real(U0),real(MX(1)),real(MX(2)),real(MX(3)), t/U0,
    !  &                 C/U0
    !       WRITE(22,*)real(U0),real(MY(1)),real(MY(2)),real(MY(3)), t/U0,
    !  &                 C/U0
    !       WRITE(23,*)real(U0),real(MZ(1)),real(MZ(2)),real(MZ(3)), t/U0,
    !  &                 C/U0

          ! ediff = e(2) - (e(3))

           ! WRITE(25,*)real(U0),real(E(1))
          if (dabs(ediff).lt.1d-4) then

            write(*,*) 'entering time_prop'

            ! v = 0.d0
            ! bx = 0.d0
            ! by = 0.d0
            ! bz = 0.d0

            CALL MATRIX(M,U0,C,CP,T,TP,V,BX,BY,BZ)

            t = 0.d0
            M0 = M
            do time = 1, 500
              t = t + time * dt
              call time_prop(M0, M, psin, psinp1, dt)

              write(*,*) M
              write(*,*) '___________________________________'
              write(*,*) M0
              call exit(-1)
              call DENCALC(M, N, MX, MY, MZ)
              write(*,*) time
              write(*,*) n
              write(*,*) mx
              write(*,*) my
              write(*,*) mz
              M0 = M
            end do

          end if

          ! write(*,*) outer

          ! write(*,*) j

101     CONTINUE

      end do

      close(26)

      END

C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE time_prop(M0, t_matrix, psin, psinp1, dt)
      IMPLICIT NONE

      INTEGER I, J, ipiv(20), INFO
      DOUBLE PRECISION U0,C,CP,T,TP,dt
      DOUBLE PRECISION V(3),BX(3),BY(3),BZ(3)

      DOUBLE COMPLEX M0(20,20), t_matrix(20, 20)
     &               , psin(20), psinp1(20), rhs(20)
      DOUBLE COMPLEX ZERO,ONE,IONE
      PARAMETER (ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),IONE=(0.D0,1.D0))


!   Right-side of the crank-nicolson propagation.

      t_matrix = zero

      do i = 1, 20
        do j = 1, 20
          t_matrix(i, j) = -0.5d0 * ione * M0(i,j) * dt 
          if (i.eq.j) then
            t_matrix(i, j) = one + t_matrix(i, j)
          end if
        end do
      end do

!   Multiplying wavefunction by Hamiltonian propagation. 

      DO I=1, 20
        rhs(i) = zero
        DO J=1,20
            rhs(i) = Rhs(i) + t_matrix(I,J) * psin(J)
        ENDDO
      ENDDO

!   Left-side of the crank-nicolson propagation.

      do i = 1, 20
        do j = 1, 20
          t_matrix(i, j) = 0.5d0 * ione * M0(i,j) * dt 
          if (i.eq.j) then
            t_matrix(i, j) = one + t_matrix(i, j)
          end if
        end do
      end do

!   Call linear equation solver for Psi^(n+1)(t)
      call zgesv(20, 1, t_matrix, 20, IPIV, rhs, 20, INFO)

      do i = 1, 20
        psinp1(i) = rhs(i)
      end do
    
      end subroutine

C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE MATRIX(M,U0,C,CP,T,TP,V,BX,BY,BZ)
      IMPLICIT NONE

      INTEGER I,J
      DOUBLE PRECISION U0,C,CP,T,TP
      DOUBLE PRECISION V(3),BX(3),BY(3),BZ(3)

      DOUBLE COMPLEX M(20,20)
      DOUBLE COMPLEX ZERO,ONE,IONE
      PARAMETER (ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),IONE=(0.D0,1.D0))

      DO 1 I=1,20
      DO 1 J=1,20
1        M(I,J) = ZERO

      M(1,1) = V(1) + V(2) + V(3) + BZ(1) - BZ(2) - BZ(3)
      M(2,2) = V(1) + V(2) + V(3) - BZ(1) + BZ(2) - BZ(3)
      M(3,3) = V(1) + V(2) + V(3) - BZ(1) - BZ(2) + BZ(3)
      M(4,4) = V(1) + V(2) + V(3) - BZ(1) + BZ(2) + BZ(3)
      M(5,5) = V(1) + V(2) + V(3) + BZ(1) - BZ(2) + BZ(3)
      M(6,6) = V(1) + V(2) + V(3) + BZ(1) + BZ(2) - BZ(3)

      M(7,7) = 2.D0*V(1) + V(2) - BZ(2) + U0
      M(8,8) = 2.D0*V(1) + V(3) - BZ(3) + U0
      M(9,9) = V(1) + 2.D0*V(2) - BZ(1) + U0
      M(10,10) = 2.D0*V(2)+V(3) - BZ(3) + U0
      M(11,11) = V(1)+2.D0*V(3) - BZ(1) + U0
      M(12,12) = V(2)+2.D0*V(3) - BZ(2) + U0
      M(13,13) = 2.D0*V(1) + V(2) + BZ(2) + U0
      M(14,14) = 2.D0*V(1) + V(3) + BZ(3) + U0
      M(15,15) = V(1) + 2.D0*V(2) + BZ(1) + U0
      M(16,16) = 2.D0*V(2)+V(3) + BZ(3) + U0
      M(17,17) = V(1)+2.D0*V(3) + BZ(1) + U0
      M(18,18) = V(2)+2.D0*V(3) + BZ(2) + U0

      M(19,19) = V(1) + V(2) + V(3) - BZ(1) - BZ(2) - BZ(3)
      M(20,20) = V(1) + V(2) + V(3) + BZ(1) + BZ(2) + BZ(3)

      M(1,7) =  TP*ONE + IONE*CP
      M(1,8) =  T*ONE - IONE*C
      M(1,10) = T*ONE - IONE*C
      M(1,12) = TP*ONE + IONE*CP

      M(2,8) =  -T*ONE - IONE*C
      M(2,9) =   T*ONE - IONE*C
      M(2,10) = -T*ONE - IONE*C
      M(2,11) = -T*ONE + IONE*C

      M(3,7) =  -TP*ONE + IONE*CP
      M(3,9) =  -T*ONE - IONE*C
      M(3,11) =  T*ONE + IONE*C
      M(3,12) = -TP*ONE + IONE*CP

      M(4,13) = -TP*ONE + IONE*CP
      M(4,14) = -T*ONE - IONE*C
      M(4,16) = -T*ONE - IONE*C
      M(4,18) = -TP*ONE + IONE*CP

      M(5,14) =  T*ONE - IONE*C
      M(5,15) = -T*ONE - IONE*C
      M(5,16) =  T*ONE - IONE*C
      M(5,17) =  T*ONE + IONE*C

      M(6,13) =  TP*ONE + IONE*CP
      M(6,15) =  T*ONE - IONE*C
      M(6,17) = -T*ONE + IONE*C
      M(6,18) =  TP*ONE + IONE*CP

      M(7,8) =    T*ONE + IONE*C
      M(7,9) =   -T*ONE + IONE*C
      M(8,11) =  -TP*ONE - IONE*CP
      M(9,10) =  -TP*ONE + IONE*CP
      M(10,12) = -T*ONE + IONE*C
      M(11,12) = -T*ONE - IONE*C

      M(13,14) =  T*ONE - IONE*C
      M(13,15) = -T*ONE - IONE*C
      M(14,17) = -TP*ONE + IONE*CP
      M(15,16) = -TP*ONE - IONE*CP
      M(16,18) = -T*ONE - IONE*C
      M(17,18) = -T*ONE + IONE*C

C**   Now the transverse magnetic field parts:

      M(1,5) = ONE*BX(3) + IONE*BY(3)
      M(1,6) = ONE*BX(2) + IONE*BY(2)
      M(1,19) = ONE*BX(1) - IONE*BY(1)
      M(2,4) = ONE*BX(3) + IONE*BY(3)
      M(2,6) = ONE*BX(1) + IONE*BY(1)
      M(2,19) = ONE*BX(2) - IONE*BY(2)
      M(3,4) = ONE*BX(2) + IONE*BY(2)
      M(3,5) = ONE*BX(1) + IONE*BY(1)
      M(3,19) = ONE*BX(3) - IONE*BY(3)
      M(4,20) = ONE*BX(1) + IONE*BY(1)
      M(5,20) = ONE*BX(2) + IONE*BY(2)
      M(6,20) = ONE*BX(3) + IONE*BY(3)
      M(7,13) = ONE*BX(2) + IONE*BY(2)
      M(8,14) = ONE*BX(3) + IONE*BY(3)
      M(9,15) = ONE*BX(1) + IONE*BY(1)
      M(10,16) = ONE*BX(3) + IONE*BY(3)
      M(11,17) = ONE*BX(1) + IONE*BY(1)
      M(12,18) = ONE*BX(2) + IONE*BY(2)

      DO 2 I=1,19
      DO 2 J=I+1,20
         M(J,I) = DCONJG(M(I,J))
2     CONTINUE

      RETURN
      END
C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE DENCALC(M,N,MX,MY,MZ)
      IMPLICIT NONE

      INTEGER I
      DOUBLE COMPLEX M(20,20),C(20),NUD(3)
      DOUBLE PRECISION N(3),MX(3),MY(3),MZ(3),
     &                 PHI1(3),PHI2(3),PHI3(3)

      DO 1 I=1,20
         C(I) = M(I,1)
1     CONTINUE

      PHI1(1) = 1.D0
      PHI1(2) = 0.D0
      PHI1(3) = 0.D0

      PHI2(1) = 0.D0
      PHI2(2) = 1.D0
      PHI2(3) = 0.D0

      PHI3(1) = 0.D0
      PHI3(2) = 0.D0
      PHI3(3) = 1.D0

      DO 10 I=1,3

         N(I) = ( CDABS(C(1))**2 + CDABS(C(2))**2 + CDABS(C(3))**2
     &          + CDABS(C(4))**2 + CDABS(C(5))**2 + CDABS(C(6))**2
     &          + CDABS(C(19))**2 + CDABS(C(20))**2)
     &         *(PHI1(I)+PHI2(I)+PHI3(I))
     &    + (CDABS(C(7))**2 + CDABS(C(13))**2)*(2.D0*PHI1(I) + PHI2(I))
     &    + (CDABS(C(8))**2 + CDABS(C(14))**2)*(2.D0*PHI1(I) + PHI3(I))
     &    + (CDABS(C(9))**2 + CDABS(C(15))**2)*(PHI1(I) + 2.D0*PHI2(I))
     &    + (CDABS(C(10))**2 + CDABS(C(16))**2)*(2.D0*PHI2(I) + PHI3(I))
     &    + (CDABS(C(11))**2 + CDABS(C(17))**2)*(PHI1(I) + 2.D0*PHI3(I))
     &    + (CDABS(C(12))**2 + CDABS(C(18))**2)*(PHI2(I) + 2.D0*PHI3(I))

         MZ(I) = CDABS(C(1))**2*(PHI1(I)-PHI2(I)-PHI3(I))
     &         + CDABS(C(2))**2*(PHI2(I)-PHI1(I)-PHI3(I))
     &         + CDABS(C(3))**2*(PHI3(I)-PHI1(I)-PHI2(I))
     &         + CDABS(C(4))**2*(PHI2(I)+PHI3(I)-PHI1(I))
     &         + CDABS(C(5))**2*(PHI1(I)+PHI3(I)-PHI2(I))
     &         + CDABS(C(6))**2*(PHI1(I)+PHI2(I)-PHI3(I))
     &         + (CDABS(C(13))**2-CDABS(C(7))**2)*PHI2(I)
     &         + (CDABS(C(14))**2-CDABS(C(8))**2)*PHI3(I)
     &         + (CDABS(C(15))**2-CDABS(C(9))**2)*PHI1(I)
     &         + (CDABS(C(16))**2-CDABS(C(10))**2)*PHI3(I)
     &         + (CDABS(C(17))**2-CDABS(C(11))**2)*PHI1(I)
     &         + (CDABS(C(18))**2-CDABS(C(12))**2)*PHI2(I)
     &         + (CDABS(C(20))**2-CDABS(C(19))**2)
     &          *(PHI1(I)+PHI2(I)+PHI3(I))

         NUD(I) = C(1)*DCONJG(C(19))*PHI1(I)
     &          + C(2)*DCONJG(C(19))*PHI2(I)
     &          + C(3)*DCONJG(C(19))*PHI3(I)
     &          + C(4)*DCONJG(C(2))*PHI3(I)
     &          + C(4)*DCONJG(C(3))*PHI2(I)
     &          + C(5)*DCONJG(C(3))*PHI1(I)
     &          + C(5)*DCONJG(C(1))*PHI3(I)
     &          + C(6)*DCONJG(C(1))*PHI2(I)
     &          + C(6)*DCONJG(C(2))*PHI1(I)
     &          + C(13)*DCONJG(C(7))*PHI2(I)
     &          + C(14)*DCONJG(C(8))*PHI3(I)
     &          + C(15)*DCONJG(C(9))*PHI1(I)
     &          + C(16)*DCONJG(C(10))*PHI3(I)
     &          + C(17)*DCONJG(C(11))*PHI1(I)
     &          + C(18)*DCONJG(C(12))*PHI2(I)
     &          + C(20)*DCONJG(C(4))*PHI1(I)
     &          + C(20)*DCONJG(C(5))*PHI2(I)
     &          + C(20)*DCONJG(C(6))*PHI3(I)

         MX(I) = 2.D0*DREAL(NUD(I))
         MY(I) = -2.D0*DIMAG(NUD(I))

10    CONTINUE


      RETURN
      END
