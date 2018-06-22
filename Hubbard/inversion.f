      Program Inversion
C  Here I use a module to allow for the passing of global variables as well
C  as portability of the code to future programs.
      USE CGMethod
      Implicit none

      !integer, parameter :: n=4!,lwork=136
      !real(kind=8), parameter :: t=0.5d0
      real(kind=8), parameter :: dab=2.d0,dbc=-3.d0,dcd=1.5d0

      integer :: i,j,k,iter!,info
      real(kind=8) :: x(dim),fret,ftol,v(dim),h0(dim,dim)!,work(lwork)
      real(8) :: check(dim)

      write(matrix,'(a, i3, a)') '(', dim, 'f8.3)'
      write(vector,'(a, i3, a)') '(', 1, 'f16.10)'

      v(1)=0.d0
      v(2)=-dab
      v(3)=v(2)-dbc
      v(4)=v(3)-dcd

      x = v

C  Subroutine that builds out the hamiltonian by solving the linear system
C  via the DSYEV lapack routine.
      call hbuild(v,h0)

      !write(*,matrix) v
      !write(*,*) '********'
C  Calculating the difference between each point potential.
!      do i=1,n-1
!        diffv(i) = v(i)-v(i+1)
!      end do

C  Calculating the target density, which is equal to the density from solving
C  the Schrodinger Eqn in the hbuild subroutine.
      ntarget = 0.d0
      do i = 1,dim
        do j=1,2
          ntarget(i) = ntarget(i) + 2.d0*h0(i,j)**2
        end do
      end do

C  Setting an arbitrary potential difference for points 2, 3, and 4.
      do i=2,dim
        v(1) = 0.d0
        v(i) = i*.5d0
      end do

!      write(*,matrix) diffv

C  Here the subroutine frprmn is the numerical recipes subroutine that uses the
C  conjugate gradient method.
      call frprmn(v,dim,ftol,iter,fret)
!      call frprmn(diffv,n-1,ftol,iter,fret)


C  Recalculating the potential after using the CG subroutine.
!      do i=2,n
!        v(1) = 0.d0
!        v(i) = v(i-1)-diffv(i-1)
!      end do

      call dfunc(v,check)

      write(*,*) '*********** dS/dV *************'
      write(*,matrix) check
      write(*,*) '*********** V_final *************'
      write(*,matrix) v
      write(*,*) '********** V_original ***********'
      write(*,*) iter
      write(*,matrix) x

      end

!*********************************************************************
!***  Begin functions and Subroutines                              ***
!*********************************************************************
