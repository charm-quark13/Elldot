      program dcomp
      implicit none

      integer :: i,j
      real(8) :: rmat(8,8),gmat(8,8),fmat(8,8)
      real(8) :: hamgen(4,4),hamreal(4,4),eng(4),enr(4)
      character(20) :: dmat,hamilton,vector

      write(dmat,'(a,i1,a)') '(', 8,'f14.10)'
      write(hamilton,'(a, i3, a)') '(', 4, 'f13.7)'
      write(vector,'(a, i3, a)') '(', 1, 'f16.10)'

      gmat = 0.d0
      rmat = 0.d0
      fmat = 0.d0

      hamgen = 0.d0
      hamreal = 0.d0

!      do i=1,8
!        do j=1,8
!          if (dabs(fmat(i,j)).gt.1.d-9) then!1.d-8 ) then
!            write(*,*) 'fmat(',i,',',j,')',' is greater than zero.'
!            write(*,*) 'the difference is:',fmat(i,j)
!            write(*,*) 'the generalized value is:',gmat(i,j)
!            write(*,*) 'the real version value is:',rmat(i,j)
!          end if
!        end do
!      end do

      open(30,file='hmat-gen.txt')
      read(30,hamilton) hamgen
      close(30)

      open(40,file='hmat-real.txt')
      read(40,hamilton) hamreal
      close(40)

      write(*,*) 'hamiltonian comparison'

      write(*,hamilton) transpose(hamgen-hamreal)
      write(*,*)
      write(*,hamilton) transpose(hamgen)
      write(*,*)
      write(*,hamilton) transpose(hamreal)
      write(*,*)
!      do i=1,4
!        do j=1,4
!          if (dabs(hamgen(i,j)-hamreal(i,j)).ge.1.d-9) then
!            write(*,*) 'matrix value',i,j,'not equal.'
!            write(*,*) 'general hamiltonian value:',hamgen(i,j)
!            write(*,*) 'real ver hamiltonian value:',hamreal(i,j)
!          end if
!        end do
!      end do

      open(50,file='real-energy.txt')
      read(50,vector) enr
      close(50)

      open(60,file='gen-energy.txt')
      read(60,vector) eng
      close(60)

      do i=1,4
        write(*,*) dabs(eng(i))-dabs(enr(i))
        write(*,*)
      end do
      write(*,*) '^^^ energy comparisons ^^^'
      end program
