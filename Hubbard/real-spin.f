      Program Inversion
C  Here I use a module to allow for the passing of global variables as well
C  as portability of the code to future programs.
      USE testCG
!      USE interactingHD

      Implicit none

      integer :: i,k,iter,it

      real(8) :: ftol,fret,x
      real(8) :: V(dim*2),h0(dim,dim),intn(dim*2),dens(dim*2)
      real(8) :: vstart(dim*2),vhxc(2),bxc(6)
      real(8) :: Bx(sites), By(sites), Bz(sites)
      real(8) :: htest(6,6),u,tau1(3),tau2(3),txc1,txc2

      write(matrix,'(a, i3, a)') '(', dim, 'f13.7)'
      write(matprint,'(a, i3, a)') '(', sites, 'e16.6)'
      write(vector,'(a, i3, a)') '(', 1, 'f16.10)'
      write(intprint,'(a, i3, a)') '(', 6, 'e16.6)'

      call Pauli(sig)

      intn = 0.d0
      v = 0.d0

      Bx(1) = .4d0
      Bx(2) = .0d0
!      Bx = 0.d0
      By = 0.d0
      Bz(1) = .0d0
      Bz(2) = .4d0

***************************************************************************
***   Setting the initial potentials which will generate our target density.
***************************************************************************
      do i=1,sites
        v(1) = 2.5d0
        if (i.ne.1) then
          v(i) = -2.5d0
        end if
      end do

      do i=1,sites
        v(sites+i) = Bx(i)
        v(sites*2+i) = By(i)
        v(sites*3+i) = Bz(i)
      end do

!      v(7) = 0.d0

      !vi=v

      vstart = v

***************************************************************************
***   Solving the initial Schrodinger equation for our system via hbuild.
***************************************************************************
      do it=10,10

        do i=1,sites
          v(1) = 2.5d0
          if (i.ne.1) then
            v(i) = -2.5d0
          end if
        end do

        do i=1,sites
          v(sites+i) = Bx(i)
          v(sites*2+i) = By(i)
          v(sites*3+i) = Bz(i)
        end do

        U = it/10.d0

        ntarget = 0.d0

        call interHam(v,U,htest)
        call intdens(ntarget,htest)

        write(*,vector) ntarget
        write(*,*) '^^^^^ int_n ^^^^^'

        h0 = 0.d0

***************************************************************************
***   Creating the target spin-density vector which will be used in the
***   optimization subroutine to find our Kohn-Sham potential.
***************************************************************************


***   Using normalized density to our advantage to reduce the dimensionality of
***   the problem, as we can simply calculate the density difference between
***   sites. This means, we can choose to set either v(1) or v(2) to zero.

***************************************************************************
***   Calling Numerical Recipes' conjugate gradient method optimiztion
***   subroutine.
***************************************************************************


        call frprmn(v,dim*2,ftol,iter,fret)

        write(*,*) '*********************************'
        write(*,*) iter, ftol, fret


***   Final check to ensure ntarget is equivalent to the density found by our
***   conjugate gradient method.
        call hbuild(v,hmat)
        call densvec(dens,hmat)

        write(*,vector) v
        write(*,*) '^^^^^ v_KS ^^^^^'
        write(*,vector) vstart
        write(*,*) '^^^^^^ v_start ^^^^^^'

        vhxc = 0.d0
        do i=1,2
          vhxc(i) = v(i) - vstart(i)
        end do

        bxc = 0.d0
        do i=1,6
          k=i+2
          bxc(i) = v(k) - vstart(k)
        end do

        write(*,vector) bxc
        write(*,*) '^^^^^^ Bxc ^^^^^^'

        tau1(1) = dens(5)*bxc(5) - dens(7)*bxc(3)
        tau1(2) = -(dens(3)*bxc(5)-dens(7)*bxc(1))
        tau1(3) = dens(3)*bxc(3) - dens(5)*bxc(1)

        tau2(1) = dens(6)*bxc(6) - dens(8)*bxc(4)
        tau2(2) = -(dens(4)*bxc(6)-dens(8)*bxc(2))
        tau2(3) = dens(4)*bxc(4) - dens(6)*bxc(2)

!        write(*,vector) v
        txc1 = 0.d0
        txc2 = 0.d0
        do i=1,3
          txc1 = txc1 + tau1(i)**2
          txc2 = txc2 + tau2(i)**2
        end do

        txc1 = dsqrt(txc1)
        txc2 = dsqrt(txc2)

        open(100,file='xc-torque-04.txt',access='append')

        if (it.eq.1) then
          write(100,*) 'U','xc torque 1', 'xc torque 2'
        end if
        write(*,*) u,txc1,-txc2
!        write(*,vector) tau1
!        write(*,*) '******'
!        write(*,vector) tau2
        write(*,*) '*** end of loop ***'
      end do

      end
