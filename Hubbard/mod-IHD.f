      Module interactingHD

      implicit none

      real(8) :: cn(6)

      contains

      subroutine interHam(v,U,ham)
      implicit none

      real(8),intent(in) :: v(dim*2),U
      real(8),intent(out) :: ham(6,6)
      integer,parameter :: lwork=136

      integer :: i,j,k,jj,info
      real(8) :: wn(dim),work(lwork),vr(dim,dim),vl(dim,dim)
      real(8) :: x,dum
      real(8) :: vec(3),Ba(2),Bb(2)

      ham = 0.d0

      do i=1,3
        ham(i,i+1) = -dsqrt(2.d0)*t
        ham(i+1,i) = ham(i,i+1)
      end do

      do i=1,2
        Ba(i) = v(3) + v(5)*(-1.d0)**mod(i+1,2)
        Bb(i) = v(4) + v(6)*(-1.d0)**mod(i+1,2)
      end do

      ham(1,1) = 2.d0*v(1) + U
      ham(2,2) = V(1) + v(2)
      ham(3,3) = 2.d0*v(2) + U

      k = 0
      do i=1,3,2
        k = k + 1
        vec(i) = (-1.d0)**(k)*Ba(k) + (-1.d0)**(k+1)*Bb(k)
      end do

!      vec(1) = -1.d0*Ba(1) + Bb(1)
      vec(2) = v(7)-v(8)
!      vec(3) = Ba(2) - Bb(2)

      do j=1,3
        ham(2,j) = vec(j)
        ham(j+3,2) = vec(j)
      end do

      do i=4,6
        do j=4,6
          if (i.eq.j) then
            do k=1,2
              ham(i,j) = ham(i,j) + v(k)
            end do
            ham(i,j) = ham(i,j)
     &                      + (v(7)+v(8))*(-1.d0)**(mod(i+1,3))
          end if

          ham(i,i+1) = Ba(2) + Bb(2)
          ham(i+1,i) = Ba(1) + Bb(1)
        end do
      end do

      ham(5,5) = ham(5,5) - 1.d0

      call dgeev('n','v',6,ham,6,cn,wn,VL,6,VR,6
     &                                   ,WORK,LWORK,INFO)

*********************************************************************
***     Sort the eigenvalues in ascending order
*********************************************************************
      do I=1,6
        do J=i+1,6
          if (cn(I).GE.cn(J)) THEN
             DUM = cn(I)
             cn(I) = cn(J)
             cn(J) = DUM

             DUM = Wn(I)
             Wn(I) = Wn(J)
             Wn(J) = DUM

            do JJ=1,dim
              DUM = VR(JJ,I)
              VR(JJ,I) = VR(JJ,J)
              VR(JJ,J) = DUM

              DUM = VL(JJ,I)
              VL(JJ,I) = VL(JJ,J)
              VL(JJ,J) = DUM
            end do
          end if
        end do
      end do
      ham = vr

      end subroutine interHam

      end module
