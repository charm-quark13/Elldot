PROGRAM THREESPIN
implicit none

integer :: I,J,MODE,INFO,outer,it,time
integer, PARAMETER :: LWORK = 500

real(8) :: U0,C,CP,T,TP,V(3),BX(3),BY(3),BZ(3),t2,dt, &
           v0(3),BX0(3),BY0(3),BZ0(3)
real(8) :: E(20),RWORK(100),N(3),MX(3),MY(3),MZ(3),ediff

complex(8) :: M(20,20),WORK(LWORK), psin(20), psinp1(20) &
               , m0(20,20), temp_mat(20,20), m_org(20,20)
complex(8), PARAMETER :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),IONE=(0.D0,1.D0)

character(50) :: wmat, wsite

write(wmat,'(a)') '(5f15.9)'

dt = .01d0
t2 = 0.d0
U0 = 1.d0

open(1, file='planar_field.txt')

do i=1,3
    read(1,*) v(i),bx(i),by(i),bz(i)
end do

close(1)

v0 = v
bx0 = bx
by0 = BY
bz0 = BZ

T = .6d0
TP = T

C = 0.5d0
CP = C

open(1, file='constants.txt')
read(1, *) c, t
close(1)

MODE = 2

!        DO 5 I=1,3
!           READ(1,*)V(I),BX(I),BY(I),BZ(I)
!5       CONTINUE



CALL MATRIX(M,U0,C,CP,T,TP,V,BX,BY,BZ)

m_org = m

CALL ZHEEV( 'V', 'U', 20, M, 20, E, WORK, LWORK, RWORK, INFO )

do it = 1, 20
    psin(it) = m(it, 1)
end do

! C**----------------------------------------------------------------------
! C**   calculate the densities:
! C**----------------------------------------------------------------------
CALL DENCALC(M,N,MX,MY,MZ)

do i = 1, 3
    write(wsite, '(a, i1, a)') 'site-', i, '_tddensity.txt'
    open(20+i, file=wsite)
end do

do i = 1, 3
    write(20+i, wmat) t2, dble(n(i)), dble(mx(i)), dble(my(i)), dble(mz(i))
end do



! do i = 1, 20
!     write(*,*) m(i, 1) - psin(i)
! end do

temp_mat = zero

! do i = 1, 3
!     v(i) = 0.0d0
!     ! bx(i) = 0.d0
!     bx(i) = 0.d0
!     by(i) = 0.d0
!     ! bz(i) = 1.d-12
!     bz(i) = 1.d-11
! end do

! v(1) = -1.d0

CALL MATRIX(M,U0,C,CP,T,TP,V,BX,BY,BZ)

do time = 1, 1000
    t2 = t2 + dt

    ! do i=1, 20
    !     write(*,*) m(i,1) - m0(i,1)
    ! end do

    if (time.eq.1) then
        do i = 1, 3
            v(i) = 0.0d0
            ! bx(i) = 0.d0
            bx(i) = 0.d0
            by(i) = 0.d0
            ! bz(i) = 1.d-12
            bz(i) = 1.d-5
        end do
        CALL MATRIX(M,U0,C,CP,T,TP,V,BX,BY,BZ)
    elseif (time.ge.2.and.time.le.10) then
        do i = 1, 3
            v(i) = 0.0d0
            ! bx(i) = 0.d0
            bx(i) = 0.d0
            by(i) = 0.d0
            ! bz(i) = 1.d-12
            bz(i) = 0.d0
        end do
        CALL MATRIX(M,U0,C,CP,T,TP,V,BX,BY,BZ)
    end if

    ! CALL MATRIX(m,U0,C,CP,T,TP,V,BX0,BY0,BZ0)
    call time_prop(m0, M, psin, psinp1, dt)

    do j=1, 20
        do i=1, 20
            temp_mat(i,1) = psinp1(i)   
        !     if (i.eq.1) then
        !         temp_mat(i,j) = psinp1(i)
        !     else
        !         temp_mat(i,j) = zero
        !     end if
        end do
    end do

    ! do i = 1, 20
    !     write(*,*) psin(i) - psinp1(i)
    ! end do

    call DENCALC(temp_mat, N, MX, MY, MZ)

    ! if (time.eq.1) call exit(-1)

    do i = 1, 3
        write(20+i, wmat) t2, dble(n(i)), dble(mx(i)), dble(my(i)), dble(mz(i))
        ! write(20+i,*) '##############################'
    end do
    ! write(*,*) t2, (n(1) + n(2) + n(3))
    ! CALL MATRIX(M,U0,C,CP,T,TP,V,BX,BY,BZ)

    ! CALL ZHEEV( 'V', 'U', 20, M, 20, E, WORK, LWORK, RWORK, INFO )

    ! do i = 1, 3
    !     call MATRIX(M_org,U0,C,CP,T,TP,V0,BX0,BY0,BZ0)
    !     CALL ZHEEV( 'V', 'U', 20, M_org, 20, E, WORK, LWORK, RWORK, INFO )
    !     call dencalc(m_org, n, mx, my, mz)
    !     write(20+i, wmat) t2, dble(n(i)), dble(mx(i)), dble(my(i)), dble(mz(i))
    !     write(20+i,*) '##############################'
    ! end do

    psin = psinp1

    ! M0 = M
end do

do i=1, 3
    close(20+i)
end do

! write(*,*) outer

! write(*,*) j



END PROGRAM 

!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------

SUBROUTINE time_prop(m0, t_matrix, psin, psinp1, dt)
IMPLICIT NONE

integer :: I, J, ipiv(20), INFO
real(8) :: U0,C,CP,T,TP,dt
real(8) :: V(3),BX(3),BY(3),BZ(3)

complex(8) :: M0(20,20), t_matrix(20, 20) &
              , psin(20), psinp1(20), rhs(20)
complex(8), PARAMETER :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),IONE=(0.D0,1.D0)


!   Right-side of the crank-nicolson propagation.
! write(*,*) '**********************'
! do i=1, 20
!     do j=1, 20
        
!         write(*,*) t_matrix(i,j)
        
!     end do
! end do
! write(*,*) '**********************'

do i = 1, 20
    do j = 1, 20
        m0(i, j) = -0.5d0 * ione * t_matrix(i,j) * dt 
        if (i.eq.j) then
            m0(i, j) = one + m0(i, j)
        end if
    end do
end do

!   Multiplying wavefunction by Hamiltonian propagation. 

DO I=1, 20
    rhs(i) = zero
    DO J=1,20
        rhs(i) = Rhs(i) + m0(I,J) * psin(J)
    ENDDO
ENDDO

!   Left-side of the crank-nicolson propagation.

do i = 1, 20
    do j = 1, 20
        m0(i, j) = 0.5d0 * ione * t_matrix(i,j) * dt 
        if (i.eq.j) then
            m0(i, j) = one + m0(i, j)
        end if
    end do
end do

! write(*,*) '**********************'
! do i=1, 20
       
!         write(*,*) psin(i), rhs(i)        
! end do

!   Call linear equation solver for Psi^(n+1)(t)
call zgesv(20, 1, m0, 20, IPIV, rhs, 20, INFO)

! write(*,*) '**********************'
! do i=1, 20
       
!         write(*,*) rhs(i)
        
! end do
! call exit(-1)


do i = 1, 20
    psinp1(i) = rhs(i)
end do

end subroutine

!************************************************************************
!************************************************************************
!************************************************************************

SUBROUTINE MATRIX(M,U0,C,CP,T,TP,V,BX,BY,BZ)
IMPLICIT NONE

integer :: I,J
real(8) :: U0,C,CP,T,TP
real(8) :: V(3),BX(3),BY(3),BZ(3)

complex(8) :: M(20,20)
complex(8), PARAMETER :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),IONE=(0.D0,1.D0)

DO I=1,20
    DO J=1,20
        M(I,J) = ZERO
    end do
end do
    

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

!**   Now the transverse magnetic field parts:

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

DO I=1,19
    DO J=I+1,20
        M(J,I) = DCONJG(M(I,J))
    end do
end do

RETURN
END subroutine

!************************************************************************
!************************************************************************
!************************************************************************

SUBROUTINE DENCALC(M,N,MX,MY,MZ)
IMPLICIT NONE

integer :: I
complex(8) :: M(20,20),C(20),NUD(3)
real(8) :: N(3),MX(3),MY(3),MZ(3), &
           PHI1(3),PHI2(3),PHI3(3)

DO I=1,20
    C(I) = M(I,1)
end do

PHI1(1) = 1.D0
PHI1(2) = 0.D0
PHI1(3) = 0.D0

PHI2(1) = 0.D0
PHI2(2) = 1.D0
PHI2(3) = 0.D0

PHI3(1) = 0.D0
PHI3(2) = 0.D0
PHI3(3) = 1.D0

DO I=1,3

    N(I) = ( CDABS(C(1))**2 + CDABS(C(2))**2 + CDABS(C(3))**2 &
              + CDABS(C(4))**2 + CDABS(C(5))**2 + CDABS(C(6))**2 &
              + CDABS(C(19))**2 + CDABS(C(20))**2) &
             *(PHI1(I)+PHI2(I)+PHI3(I)) &
        + (CDABS(C(7))**2 + CDABS(C(13))**2)*(2.D0*PHI1(I) + PHI2(I)) &
        + (CDABS(C(8))**2 + CDABS(C(14))**2)*(2.D0*PHI1(I) + PHI3(I)) &
        + (CDABS(C(9))**2 + CDABS(C(15))**2)*(PHI1(I) + 2.D0*PHI2(I)) &
        + (CDABS(C(10))**2 + CDABS(C(16))**2)*(2.D0*PHI2(I) + PHI3(I)) &
        + (CDABS(C(11))**2 + CDABS(C(17))**2)*(PHI1(I) + 2.D0*PHI3(I)) &
        + (CDABS(C(12))**2 + CDABS(C(18))**2)*(PHI2(I) + 2.D0*PHI3(I))

    MZ(I) = CDABS(C(1))**2*(PHI1(I)-PHI2(I)-PHI3(I)) &
             + CDABS(C(2))**2*(PHI2(I)-PHI1(I)-PHI3(I)) &
             + CDABS(C(3))**2*(PHI3(I)-PHI1(I)-PHI2(I)) &
             + CDABS(C(4))**2*(PHI2(I)+PHI3(I)-PHI1(I)) &
             + CDABS(C(5))**2*(PHI1(I)+PHI3(I)-PHI2(I)) &
             + CDABS(C(6))**2*(PHI1(I)+PHI2(I)-PHI3(I)) &
             + (CDABS(C(13))**2-CDABS(C(7))**2)*PHI2(I) &
             + (CDABS(C(14))**2-CDABS(C(8))**2)*PHI3(I) &
             + (CDABS(C(15))**2-CDABS(C(9))**2)*PHI1(I) &
             + (CDABS(C(16))**2-CDABS(C(10))**2)*PHI3(I) &
             + (CDABS(C(17))**2-CDABS(C(11))**2)*PHI1(I) &
             + (CDABS(C(18))**2-CDABS(C(12))**2)*PHI2(I) &
             + (CDABS(C(20))**2-CDABS(C(19))**2) &
              *(PHI1(I)+PHI2(I)+PHI3(I))

    NUD(I) = C(1)*DCONJG(C(19))*PHI1(I) &
              + C(2)*DCONJG(C(19))*PHI2(I) &
              + C(3)*DCONJG(C(19))*PHI3(I) &
              + C(4)*DCONJG(C(2))*PHI3(I) &
              + C(4)*DCONJG(C(3))*PHI2(I) &
              + C(5)*DCONJG(C(3))*PHI1(I) &
              + C(5)*DCONJG(C(1))*PHI3(I) &
              + C(6)*DCONJG(C(1))*PHI2(I) &
              + C(6)*DCONJG(C(2))*PHI1(I) &
              + C(13)*DCONJG(C(7))*PHI2(I) &
              + C(14)*DCONJG(C(8))*PHI3(I) &
              + C(15)*DCONJG(C(9))*PHI1(I) &
              + C(16)*DCONJG(C(10))*PHI3(I) &
              + C(17)*DCONJG(C(11))*PHI1(I) &
              + C(18)*DCONJG(C(12))*PHI2(I) &
              + C(20)*DCONJG(C(4))*PHI1(I) &
              + C(20)*DCONJG(C(5))*PHI2(I) &
              + C(20)*DCONJG(C(6))*PHI3(I)

    MX(I) = 2.D0*DREAL(NUD(I))
    MY(I) = -2.D0*DIMAG(NUD(I))

end do


RETURN
END subroutine 
