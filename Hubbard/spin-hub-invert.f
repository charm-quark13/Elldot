      Program Inversion
C  Here I use a module to allow for the passing of global variables as well
C  as portability of the code to future programs.
      USE ComplexCGMethod
      Implicit none

      integer :: i,j,k,counter

      complex(8) :: dens(dim*2),ntarget(dim*2)
      complex (8) :: V(dim*2),h0(dim,dim),phi(spin,sites)
      complex (8) :: phi(dim)

      write(matrix,'(a, i3, a)') '(', dim, 'f8.3)'

      Bx = zero
      By = zero
      Bz(1) = one
      Bz(2) = cmplx(0.5, kind=8)*one

***************************************************************************
***   Setting the initial potentials which will generate our target density.
***************************************************************************
      do i=1,sites
        v(1) = zero
        if (i.ne.1) then
          v(i) = cmplx(i,kind=8)*(-1.d0)**(i)
        end if
      end do

      do i=1,sites
        v(sites+i) = Bx(i)
        v(sites*2+i) = By(i)
        v(sites*3+i) = Bz(i)
      end do

***************************************************************************
***   Solving the initial Schrodinger equation for our system via hbuild.
***************************************************************************
      call hbuild(v,h0)

      ntarget = 0.d0

      phi = zero

***************************************************************************
***   Creating the target spin-density vector which will be used in the
***   optimization subroutine to find our Kohn-Sham potential.
***************************************************************************
      call wfmatrix(phi,h0)

      call densvec(ntarget,phi)

      do i=1,sites
        v(1) = zero
        if (i.ne.1) then
          v(i) = cmplx(-.75,kind=8)
        end if
      end do

      Bx = zero
      By = zero
      Bz(1) = .3d0*one
      Bz(2) = cmplx(0.25, kind=8)*one

      do i=1,sites
        v(sites+i) = Bx(i)
        v(sites*2+i) = By(i)
        v(sites*3+i) = Bz(i)
      end do

***************************************************************************
***   Calling Numerical Recipes' conjugate gradient method optimiztion
***   subroutine.
***************************************************************************
      call frprmn(v,dim*sites,ftol,iter,fret)

      write(*,matrix) v

      end
