      program cgtest
      use test
      implicit none

      integer, parameter :: num=8
      integer :: i,iter
      real(8) :: pot(8),dv(8),fret,ftol

      write(vector,'(a, i3, a)') '(', 1, 'f16.10)'

      do i=1,num
        pot(i) = 0.d0
        pot(i) = 1.d0*(-1.d0)**i
      end do

      do i=1,num
        dd(i) = mod(i,2)
      end do

      call frprmn(pot,num,ftol,iter,fret)

      write(*,*) 'through frprmn'

      write(*,vector) pot
      !write(*,*) ftol, iter, fret

      end program
