      Program Inversion
C  Here I use a module to allow for the passing of global variables as well
C  as portability of the code to future programs.
      USE testCG
!      USE interactingHD

      Implicit none

      integer :: i,iter

      real(8) :: x,ftol,fret
      real (8) :: V(dim*2),h0(dim,dim),vi(dim*2)
      real(8) :: Bx(sites), By(sites), Bz(sites)
      real(8) :: htest(6,6),u

      character(20) :: intprint

      write(matrix,'(a, i3, a)') '(', dim, 'f13.7)'
      write(matprint,'(a, i3, a)') '(', sites, 'e16.6)'
      write(vector,'(a, i3, a)') '(', 1, 'f16.10)'
      write(intprint,'(a, i3, a)') '(', 6, 'e16.6)'

      call Pauli(sig)

      vi = 0.d0

      Bx(1) = .05d0
      Bx(2) = -.02d0
!      Bx = 0.d0
      By = 0.d0
!      Bz(1) = .06d0
      Bz(2) = -.15d0

***************************************************************************
***   Setting the initial potentials which will generate our target density.
***************************************************************************
      do i=1,sites
        v(1) = 0.d0
        if (i.ne.1) then
          v(i) = -1.5d0
        end if
      end do

      do i=1,sites
        v(sites+i) = Bx(i)
        v(sites*2+i) = By(i)
        v(sites*3+i) = Bz(i)
      end do

      v(7) = 0.d0

      vi=v

      write(*,vector) v
      write(*,*) 'V0', '^^^^^^^^^^^^^^'

***************************************************************************
***   Solving the initial Schrodinger equation for our system via hbuild.
***************************************************************************
      call hbuild(v,h0)

      U = 1.d0

      call interHam(v,U,htest)

      write(*,intprint) transpose(htest)
      call exit(-1)
      ntarget = 0.d0

!      phi = 0.d0

***************************************************************************
***   Creating the target spin-density vector which will be used in the
***   optimization subroutine to find our Kohn-Sham potential.
***************************************************************************

      call densvec(ntarget,h0)


***   Using normalized density to our advantage to reduce the dimensionality of
***   the problem, as we can simply calculate the density difference between
***   sites. This means, we can choose to set either v(1) or v(2) to zero.

      write(*,vector) ntarget
      write(*,*) '^^^^^^^^^^^ ntarget ^^^^^^^^^^^^'
!      do iter = 1,100
        x = -3.d0
        do i=1,sites
          v(1) = 0.d0
          if (i.ne.1) then
            v(i) = .5d0
          end if
        end do

        Bx(1) = -.125d0
        Bx(2) = .03d0
!        Bx = 0.d0
        By = 0.d0

        Bz(1) = 0.d0
        Bz(2) = x

        do i=1,sites
          v(sites+i) = Bx(i)
          v(sites*2+i) = By(i)
          v(sites*3+i) = Bz(i)
        end do

***************************************************************************
***   Calling Numerical Recipes' conjugate gradient method optimiztion
***   subroutine.
***************************************************************************

      write(*,vector) v
      write(*,*) '^^^^^ v_test ^^^^^'
      !call hbuild(v,hmat)

      call frprmn(v,dim*2,ftol,iter,fret)

      write(*,*) '*********************************'
      write(*,*) iter, ftol, fret

      write(*,vector) v
      write(*,*) 'v_final', '^^^^^^^^^^^^^^^^^^^^^'

      write(*,vector) vi
      write(*,*) 'v_target ^^^^^^^^^^^^^^^^^^^^^^^'

***   Final check to ensure ntarget is equivalent to the density found by our
***   conjugate gradient method.
      call hbuild(v,hmat)
      call densvec(v,hmat)
      v = ntarget - v

      write(*,vector) v

      end
