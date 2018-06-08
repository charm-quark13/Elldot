      Program Inversion
C  Here I use a module to allow for the passing of global variables as well
C  as portability of the code to future programs.
      USE RealCGMethod
      Implicit none

      integer :: i,iter, itol

      real(8) :: x,ftol,fret
      real (8) :: V(dim*2),h0(dim,dim),phi(spin,sites),dSdV(dim*2)
      real(8) :: Bx(sites), By(sites), Bz(sites)
!      character(len=20):: Pmat
      write(matrix,'(a, i3, a)') '(', dim, 'f13.7)'
      write(matprint,'(a, i3, a)') '(', sites, 'e16.6)'
      write(vector,'(a, i3, a)') '(', 1, 'f13.8)'
!      write(Pmat,'(a, i3, a)') '(', 8, 'f13.7)'

      call Pauli(sig)

!      write(*,Pmat) transpose(sig)

!      call exit(-1)

      Bx = 0.d0
      By = 0.d0
      Bz(1) = 1.d0
      Bz(2) = -.5d0

***************************************************************************
***   Setting the initial potentials which will generate our target density.
***************************************************************************
      do i=1,sites
        v(1) = 0.d0
        if (i.ne.1) then
          v(i) = i*(-1.d0)**(i)
        end if
      end do

      do i=1,sites
        v(sites+i) = Bx(i)
        v(sites*2+i) = By(i)
        v(sites*3+i) = Bz(i)
      end do

      dsdv = v

      write(*,vector) v
      write(*,*) 'V0', '^^^^^^^^^^^^^^'

***************************************************************************
***   Solving the initial Schrodinger equation for our system via hbuild.
***************************************************************************
      call hbuild(v,h0)

      ntarget = 0.d0

!      phi = 0.d0

***************************************************************************
***   Creating the target spin-density vector which will be used in the
***   optimization subroutine to find our Kohn-Sham potential.
***************************************************************************
      call wfmatrix(phi,h0)

***   h0 and phi are being created and positions assigned correctly (6/4/18 EP)

      call densvec(ntarget,phi)

***   densvec and test matrix within densvec is being correctly calculated
***   and values assigned to correct positions (6/4/18 EP)

***   Using normalized density to our advantage to reduce the dimensionality of
***   the problem. This allows us to pin Bz(b) to zero. We can do the same with
***   the density, leading us to just the differene, or n(b) = 1 - n(a).
      !ntarget(2) = 1-ntarget(1)
      !ntarget(sites*4) = 0.d0
      !do i=sites*2+1,sites*3
      !  ntarget(i) = 0.d0
      !end do

      write(*,vector) ntarget

!      do iter = 1,100
        x = -1.d0
        do i=1,sites
          v(1) = 0.d0
          if (i.ne.1) then
            v(i) = .5d0
          end if
        end do

        Bx = 0.d0
        By = 0.d0
        Bz(1) = 1.d0
        Bz(2) = -x

        do i=1,sites
          v(sites+i) = Bx(i)
          v(sites*2+i) = By(i)
          v(sites*3+i) = Bz(i)
        end do
!
!        ftol = func(v)
!
!        write(*,*) '************'
!        write(*,*) ftol
!        write(*,*) '***',x,'***'
!
!        call dfunc(v,dSdV)
!
!        write(2,*) '***',x,'***'
!        write(2,vector) dSdV
!        write(2,*) '***********'
!      end do

***************************************************************************
***   Calling Numerical Recipes' conjugate gradient method optimiztion
***   subroutine.
***************************************************************************
      call frprmn(v,dim*2,ftol,iter,fret)

      write(*,*) '*********************************'
      write(*,*) iter, ftol, fret

      write(*,*) '**************************'
      write(*,matrix) transpose(h0)
      write(*,*) '*************  h0     ************'

      write(*,vector) v
      write(*,*) 'V_final', '^^^^^^^^^^^^^^^^^^^^^'

      write(*,vector) dsdv
      write(*,*) 'V_initial ^^^^^^^^^^^^^^^^^^^^^^^'

      end
