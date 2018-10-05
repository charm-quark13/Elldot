      Module generalCG

      IMPLICIT NONE

      integer,parameter :: spin=2, sites = 3, occ=2
      integer, parameter :: dim=spin*sites, intd=2*sites**2-sites
      real(kind=8), parameter :: t=0.5d0

      integer :: number
      real(8) :: ntarget(dim*2),En(dim)
      complex(8) :: hmat(dim,dim)
!      real(8) :: tBx(sites), tBy(sites), tBz(sites)
      complex(8) :: sig(spin*4,spin*4)
      complex(8), parameter :: Zero = (0.d0,0.d0), One = (1.d0,0.d0)
      complex(8), parameter :: IOne = (0.d0,1.d0)

      character(20) :: intprint,dmat
      character (len=30) :: matrix, matprint, vector

      contains

C  Beginning of the function being minimized through the CG method.
C  In this case, func = \int (n{v(x)} - ntarget{v(x)})**2 dx
!      function func(diffv) result(integral)
      function func(v) result(integral)
      implicit none
      integer :: i,j
      real(8) :: integral,v(dim*2),dens(dim*2)

***   We can set one scalar potential to zero, as the difference between the
***   relative difference between the two potentials is all that matters.
***   Additionally, because of the nature of the Hubbard Dimer, we can set
***   one magnetic potential to zero, as the magnetizations of the two sites
***   must satisfy: |vec(m_1)| = |vec(m_2)|.

!      v(1) = 0.d0
      v(7) = 0.d0

C  Recalculating the eigenvectors through exact diagonalization of the Hamiltonian.
      call hbuild(v,hmat)

      dens = 0.d0

      call densvec(dens,hmat)

      if (number.eq.0) then
        open(10,file='hmat-gen.txt')
        write(10,matrix) dreal(hmat)
        close(10)

        open(50,file='gen-energy.txt')
        write(50,vector) En
        close(50)
      end if

!      do j=3,4
!        hmat(:,j) = hmat(:,j)*(-1.d0)
!      end do

      do i=1,dim
        do j=1,dim
          if (dabs(dreal(hmat(i,j))).lt.1d-10) then
            hmat(i,j) = dcmplx(0.d0,dimag(hmat(i,j)))
          elseif (dabs(dimag(hmat(i,j))).lt.1d-10) then
            hmat(i,j) = dcmplx(dreal(hmat(i,j)),0.d0)
          end if
        end do
      end do

***   Using normalized density to our advantage to reduce the dimensionality of
***   the problem. This allows us to pin Bz(b) to zero. We can do the same with
***   the density, leading us to just the differene, or n(b) = 1 - n(a).

C  Calculating the difference integral between the target and current densities.
      integral = 0.d0
      do i =1,dim*2
!        integral = integral + dreal((dens(i)-ntarget(i))
!     &                          *(conjg(dens(i)-ntarget(i))))

        if (i.eq.5.or.i.eq.6) then
          integral = integral + dreal((dens(i)*ione-ntarget(i)*ione)*
     &                    (dens(i)*(-ione)-ntarget(i)*(-ione)))
        else
          integral = integral + (dens(i)-ntarget(i))**2
        end if
      end do

      end function
*************************************************************************
C  Derivative subroutine of func.
C  Derivative is of the form:
C  dF(n{v(x)}) = 2 * \int (n{v(x)} - ntarget{v(x)})*(dn/dv) dx
      !subroutine dfunc(diffv,diffn)
      subroutine dfunc(v,dSdV)

      implicit none
      integer :: j,k,counter,num
      real(8) :: v(dim*2),dens(dim*2),dSdV(dim*2)
      complex(8) :: x
      complex(8) :: dn(dim*2)!, wfmat(spin,sites)
      complex(8) :: vec(8,(occ+sites)*sites)

      call vmat(vec)
      call densvec(dens,hmat)

      dSdV = 0.d0
      dn = zero

      number = number + 1

***   Beginning of the loop to calculate the integral of dS/dV.
***   The subroutine dnvec generates the gradient vector of the density or
***   magnetization with respect to each potential on each site.

      counter = 0
      do num=1,dim
        do k=1,sites
          call dnvec(num,k,dn,vec)
          x = zero
          do j=1,dim*2
            if (j.eq.5.or.j.eq.6) then
              x = x + (dens(j)*ione-ntarget(j)*ione)*dn(j)
            else
              x = x + dcmplx((dens(j)-ntarget(j))*dn(j))
            end if
!            write(*,*) x
          end do
!          write(*,*)
          counter = counter + 1

          dSdV(counter) = 2.d0*(dreal(x)+dimag(x))
        end do
      end do

***   The derivative, dS/dV, will be all zeros when a minimum state of dens
***   has been reached (i.e.- when dens = ntarget).

      write(*,vector) dSdV
      write(*,*) '^^^^^^^^^^^^ dSdV ^^^^^^^^^^^^'

!      if (number.eq.1) then
!        call exit(-1)
!      end if

      end subroutine

***************************************************************************
      subroutine vmat(mat)
      implicit none
      complex(8),intent(out) :: mat(8,(occ+sites)*sites)
      integer :: j,k,s,num,alpha,counter
      complex(8) :: x

***   vmat generates a matrix containing the values of the derivatives of
***   psi_up and psi_dn for each occupied state with respect to each potential.
***   Please see the LaTeX write up for exact ordering of the matrix. The values
***   are used in the dnvec calculations to calculate the magnetization
***   derivatives.

      counter = 0
      do k=1,sites
        x = zero
        do j = 1,sites
          do alpha=1,2
            counter = counter + 1
            do s=1,2
              do num=1,dim
                mat(dim*(s-1)+num,counter) =
     &                      derivative(alpha,j,k,num,s)
              end do
            end do
          end do
        end do
      end do

      end subroutine vmat

***************************************************************************
      function derivative(alpha,j,k,num,sigma) result(dp)
      implicit none

      integer :: l,m,beta
      integer, intent(in) :: j,k,alpha,sigma,num
      !real(8), intent(in) :: v(dim*2)
      real(8) :: denom
      complex(8) :: numer
      complex(8):: dp, x

***************************************************************************
***   The eigenvector is stored as (a_up,b_up,a_dn,b_dn). Then, l is a
***   pointer that directs the calculation to the correct spinor/position
***   value.
***   j is the position label, where l points to spin up or spin down.
***   m is then a compound variable that encompasses spatial and
***   spin components.
***************************************************************************
      dp = zero

      if (sigma.eq.1) then
        l = 0
      else
        l = sites
      end if

      m = l + j

      do beta = 1, dim
        denom = En(alpha) - En(beta)

***   Carrying out perturbation theory multiplication, using potential
***   matrix and spinors gives:
***   ( {[V0 + V3]Phi_aup + [V1-iV2]Phi*_bup + ... } * Phi_j spinor
***
***   Taking the derivative with respect to each potential, gives a set of
***   equations for potential derivative. (i.e. - Va = V(1) so j<=sites,
***   which means it would use the first equation. Likewise,
***   Bxb=V(sites + 2), putting it in the second equation.)

        if (num.eq.1) then
          numer = hmat(m,beta)*(dconjg(hmat(k,beta))*hmat(k,alpha)+
     &                   dconjg(hmat(k+sites,beta))*hmat(k+sites,alpha))

        elseif (num.eq.2) then
          numer = hmat(m,beta)
     &                *(dconjg(hmat(k,beta))*hmat(k+sites,alpha)+
     &                        dconjg(hmat(k+sites,beta))*hmat(k,alpha))

        elseif (num.eq.3) then
!          numer = 0.d0
          numer = hmat(m,beta)
     &           *(-dconjg(hmat(k,beta))*hmat(k+sites,alpha)+
     &               dconjg(hmat(k+sites,beta))*hmat(k,alpha))*ione

        else
          numer = hmat(m,beta)
     &            *(dconjg(hmat(k,beta))*hmat(k,alpha) -
     &                dconjg(hmat(k+sites,beta))*hmat(k+sites,alpha))

        end if

        x = numer/denom

***   From nondegenerate perturbation theory, if psi_alpha = psi_beta, then we
***   must get zero for that summation value. The if statement fixes this.
***   The final line is the explicit summation over each value.

        if (alpha.eq.beta) then
          x = zero
        end if

        dp = dp + x

!        if (number.eq.1) then
!          write(*,*)
!     &         'n:', num-1, 'alpha:',alpha,'j:',j,'k:',k,'spin:',sigma
!          write(*,*) 'numerator:',dreal(numer)
!          write(*,*) 'denom:',denom
!          write(*,*)
!        end if

      end do

!      if (number.eq.1) then
!        write(*,vector) En
!        write(*,*) '^^^ En ^^^'
!      end if

      end function

***************************************************************************

      subroutine dnvec(num,k,dn,vec)
      implicit none

      integer :: i,j,k,num,alpha,counter,dum
      complex(8) :: x,y
      complex(8),intent(in) :: vec(8,(occ+sites)*sites)
      complex(8),intent(out) :: dn(dim*2)

**************************************************************************
***   Subroutine calculating the magnetizations and densities of the system at
***   each site with an occupation number of occ.
***   For the subroutine, num refers to each potential (V_0, V_1, ...), k is the
***   site where the potential is being applied, j is the site at which the
***   magnetization/density is being calculated, and alpha are the occupied
***   orbitals.
**************************************************************************

***   The equation for i is a pointer that sets the columns that will be looped
***   over in the vec matrix used to calculate the magnetizations/densities.
***   For explicit representation of the vec matrix, see the LaTeX document.

      i = (k-1)*(occ+sites)

***   Solving the first magnetization (m_0)
***   See the LaTeX for the explicit representation of the vec matrix.
***   The x value for each loop is the first part of each magnetization/density
***   equation. For m_0, it is:
***   x = [(dpsi_up/dV)(psi_up)* + (psi_up)(d(psi_up)*/dV)]
***   y = [(dpsi_dn/dV)(psi_dn)* + (psi_dn)(d(psi_dn)*/dV)]
***   d(m_0)/dV = x + y

      do j=1,sites
        x = zero
        y = zero
        do alpha = 1, 2
          x = x + vec(num,i+(j-1)*occ+alpha)*dconjg(hmat(j,alpha)) +
     &           hmat(j,alpha)*dconjg(vec(num,i+(j-1)*occ+alpha))
          y = y + vec(num+4,i+(j-1)*occ+alpha)
     &                      *dconjg(hmat(j+sites,alpha)) +
     &           hmat(j+sites,alpha)
     &                      *dconjg(vec(num+4,i+(j-1)*occ+alpha))

        end do
        dn(j) = x + y
      end do

!      if (num.eq.1.and.j.eq.1) then
!        write(*,*) 'dn(1)           :',dn(1)
!        write(*,*) 'dn(2)           :',dn(2)
!      end if

***************************************************************************
***   Solving the second magnetization (m_1)
***************************************************************************

      counter = sites
      do j=1,sites
        x = zero
        y = zero
        do alpha = 1, 2
          x = x + vec(num,i+(j-1)*occ+alpha)
     &                     *dconjg(hmat(j+sites,alpha)) +
     &           hmat(j,alpha)*dconjg(vec(num+4,i+(j-1)*occ+alpha))
          y = y + vec(num+4,i+(j-1)*occ+alpha)*dconjg(hmat(j,alpha)) +
     &           hmat(j+sites,alpha)*dconjg(vec(num,i+(j-1)*occ+alpha))
        end do
        dn(counter + j) = x + y
      end do

***************************************************************************
***   Solving the third magnetization (m_2)
***************************************************************************

      counter = sites*2
      do j=1,sites
        x = zero
        y = zero
        do alpha = 1, 2
          x = x + vec(num,i+(j-1)*occ+alpha)
     &                    *dconjg(hmat(j+sites,alpha)) +
     &           hmat(j,alpha)*dconjg(vec(num+4,i+(j-1)*occ+alpha))
          y = y + vec(num+4,i+(j-1)*occ+alpha)*dconjg(hmat(j,alpha)) +
     &           hmat(j+sites,alpha)*dconjg(vec(num,i+(j-1)*occ+alpha))
        end do
        dn(counter + j) = (x - y)*ione

!        if ((counter+j).eq.5.or.(counter+j).eq.6) then
!          write(*,*) vec(num,i+(j-1)*occ+alpha),j,k
!          write(*,*) 'num:',num,'derivative:',i+(j-1)*occ+alpha
!        end if

      end do

***************************************************************************
***   Solving the fourth magnetization (m_3)
***************************************************************************

      counter = sites*3
      do j=1,sites
        x = zero
        y = zero
        do alpha = 1, 2
          x = x + vec(num,i+(j-1)*occ+alpha)*dconjg(hmat(j,alpha)) +
     &           hmat(j,alpha)*dconjg(vec(num,i+(j-1)*occ+alpha))
          y = y + vec(num+4,i+(j-1)*occ+alpha)
     &                    *dconjg(hmat(j+sites,alpha)) +
     &           hmat(j+sites,alpha)
     &                    *dconjg(vec(num+4,i+(j-1)*occ+alpha))
        end do
        dn(counter + j) = x - y
      end do

      end subroutine dnvec

****************************************************************************
***   Density vector calculation subroutine
****************************************************************************
***   Calculates the spin density matrix in vectorized form. Each site is
***   listed as a 4-vector stacked on top of one another.
***   (i.e. - n(1) = |a_up|^2, n(3) = a_up x a_down*, n(5) = |b_up|^2)
****************************************************************************

      subroutine densvec(realn0,h0)
      implicit none

      real(8), intent(out) :: realn0(dim*2)
      complex(8) :: n0(dim*2)
      complex(8), intent(in) :: h0(dim,dim)
      integer :: i,j,k,m
      complex(8) :: x
      complex(8) :: test(spin,spin)
      complex(8) :: sigma(spin,spin), dum(spin,spin)

      test = zero

      n0 = zero

      do i=1,sites
        call spinmat(i,occ,h0,test)

**************************************************************************
***   Beginning of calculation of n, mx, my, mz. The index m runs from 1-
***   to-4 for the density and then each magnetization direction.
**************************************************************************
        !write(*,matprint) test
        !write(*,*) '^^^^^^^^ test ^^^^^^^^'
        do m=0,3
          do j=1,spin
            do k=1,spin
              sigma(j,k) = sig(m*spin+j,m*spin+k)
            end do
          end do

          dum = matmul(sigma,test)

          x = zero
          do j=1,spin
            x = x + dum(j,j)
          end do

          n0(i+m*sites) = x

        end do
      end do

      realn0 = dreal(n0)

      !write(*,*) n0(5)
      !write(*,*) n0(6)

      end subroutine densvec

****************************************************************************
      subroutine wfmatrix(nocc,psi,h0)
      implicit none

      integer, intent(in) :: nocc
      complex(8), intent(in) :: h0(dim,dim)
      complex(8), intent(out) :: psi(spin,sites)

      integer :: counter, i, j

****************************************************************************
***   Wavefunction matrix population subroutine.
****************************************************************************
***   Matrix values are wave amplitudes populated according to
***   spin and occupation site in row and column, respectively.
***   (i.e.- phi_up(a) = m_11, phi_dn(c) = m_23)
****************************************************************************
      counter = 0
      do i=1,spin
        do j=1,sites
          counter = counter + 1
          psi(i,j) = h0(counter, nocc)
        end do
      end do

      end subroutine wfmatrix

****************************************************************************

      subroutine spinmat(loc,nocc,ham,mat)
      implicit none

      integer, intent(in) :: loc,nocc
      complex(8), intent(in) :: ham(dim,dim)
      complex(8), intent(inout) :: mat(spin,spin)

      integer :: j,k,m
      complex(8) :: psi(spin,sites)

***************************************************************************
***   Generating the spin density matrix for each lattice site for all occupied
***   states.
***************************************************************************

      mat = zero
      do m=1,nocc
        call wfmatrix(m,psi,ham)
        do j=1,spin
          do k=1,spin
            mat(j,k) = mat(j,k) + psi(j,loc)*conjg(psi(k,loc))
            if (dabs(dreal(mat(j,k))).lt.1d-10) then
              mat(j,k) = complex(0.d0,dimag(mat(j,k)))
            elseif (dabs(dimag(mat(j,k))).lt.1d-10) then
              mat(j,k) = complex(dreal(mat(j,k)),0.d0)
            end if
          end do
        end do
      end do

      end subroutine spinmat

****************************************************************************

      subroutine intdens(n0,inmat)
      implicit none

      real(8),intent(out) :: n0(dim*2)
      complex(8),intent(in) :: inmat(6,6)

      integer :: i,j
      real(8) :: x

      n0 = 0.d0

      j=0
      do i=1,3,2
        j=j+1
        x = 0.d0
        x = x + dreal(inmat(i,1)*conjg(inmat(i,1)))
        x = x + .5d0*dreal(inmat(2,1)*conjg(inmat(2,1))
     &                 + inmat(4,1)*conjg(inmat(4,1))
     &                      + inmat(5,1)*conjg(inmat(5,1))
     &                         + inmat(6,1)*conjg(inmat(6,1)))

!          x = x + inmat(k,i)**2
!          x = x + .5d0*(inmat(k,2)**2 + inmat(k,4)**2
!     &                     + inmat(k,5)**2 + inmat(k,6)**2)

!          write(*,*) x/3.d0

        n0(j) = n0(j) + x*2.d0!/3.d0
      end do

      do i=1,2
        j=j+1
        x = 0.d0
        x = x + dreal(conjg(inmat(6,1))*
     &                  (inmat(5,1)+(-1.d0)**(i+1)*inmat(2,1)))
     &        + dreal(inmat(4,1)*
     &          (conjg(inmat(5,1))+(-1.d0)**(i)*conjg(inmat(2,1))))
!        x = x + inmat(6,1)*(inmat(5,1)+(-1.d0)**(i+1)*inmat(2,1))
!        x = x + inmat(4,1)*(inmat(5,1)+(-1.d0)**(i)*inmat(2,1))

!        x = x + inmat(k,6)*(inmat(k,5)+(-1.d0)**(i+1)*inmat(k,2))
!        x = x + inmat(k,4)*(inmat(k,5)+(-1.d0)**(i)*inmat(k,2))

        n0(j) = n0(j) + x*2.d0/dsqrt(2.d0)!/(3.d0*dsqrt(2.d0))
      end do

      do i=1,2
        j=j+1
        x = 0.d0
        x = x + dimag(conjg(inmat(6,1))*
     &                  (inmat(5,1)+(-1.d0)**(i+1)*inmat(2,1)))
     &        + dimag(inmat(4,1)*
     &          (conjg(inmat(5,1))+(-1.d0)**(i)*conjg(inmat(2,1))))
!        x = x + inmat(6,1)*(inmat(5,1)+(-1.d0)**(i+1)*inmat(2,1))
!        x = x + inmat(4,1)*(inmat(5,1)+(-1.d0)**(i)*inmat(2,1))

!          x = x + inmat(k,6)*(inmat(k,5)+(-1.d0)**(i+1)*inmat(k,2))
!          x = x + inmat(k,4)*(inmat(k,5)+(-1.d0)**(i)*inmat(k,2))

        n0(j) = n0(j) + x*2.d0/dsqrt(2.d0)
      end do

      do i=1,2
        j=j+1
        x = 0.d0
        x = x + dreal(inmat(4,1)*conjg(inmat(4,1)) -
     &                   inmat(6,1)*conjg(inmat(6,1)))
     &        + dreal(inmat(2,1)*conjg(inmat(5,1))
     &                      + conjg(inmat(2,1))*inmat(5,1))
     &                                             *(-1.d0)**(i+1)

!        x = x + inmat(4,k)**2 - inmat(6,k)**2
!        x = x + (inmat(2,k)*inmat(5,k)+inmat(2,k)*inmat(5,k))
!     &                                             *(-1.d0)**(i+1)

!          x = x + inmat(k,4)**2 - inmat(k,6)**2
!          x = x + (inmat(k,2)*inmat(k,5)+inmat(k,2)*inmat(k,5))
!     &                                             *(-1.d0)**(i+1)

        n0(j) = n0(j) + x!/(3.d0*2.d0)
      end do

      end subroutine intdens

***************************************************************************

      subroutine hbuild(v,mat)
      implicit none

      integer, parameter :: lwork=100
      integer :: i,j,jj,info

      real(8) :: rwork(100),v(2*dim)
      complex(8) :: mat(dim,dim),test(sites,sites),work(lwork)
!      complex(8) :: vr(dim,dim),vl(dim,dim),dum,dumv

***************************************************************************
***   Solving the Schrodinger Eqn for our system through direct diagonalization
***   of the Hamiltonian matrix.
***************************************************************************

      mat = zero

C  Setting up noninteracting Hubbard-like system Hamiltonian with t as the
C  hopping constant between lattice sites.
      do i=1,sites
        do j=1,sites
          if (i.eq.j) Then
            test(i,j) = dcmplx(v(i) + v((dim*2-sites) + i),0.d0)
            write(*,*) (dim-1)*2+i,dim*2-sites+i
          elseif (abs(i-j).eq.1) Then
            test(i,j) = dcmplx(-t,0.d0)
          else
            test(i,j) = zero
          endif
        end do
      end do

      do i=1,sites
        do j=1,sites
          mat(i,j) = test(i,j)
        end do
      end do

      do i=1,sites
        do j=1,sites
          if (i.eq.j) Then
            test(i,j) = dcmplx(v(i) - v((dim*2-sites) + i),0.d0)
          elseif (abs(i-j).eq.1) Then
            test(i,j) = dcmplx(-t,0.d0)
          else
            test(i,j) = zero
          endif
        end do
      end do

      do i=1,sites
        do j=1,sites
          mat(dim/2+i,dim/2+j) = test(i,j)
        end do
      end do

      do i=1,sites
        mat(i,dim/2+i) = dcmplx(v(sites+i),-v(sites*2+i))
        mat(dim/2+i,i) = dcmplx(v(sites+i),v(sites*2+i))
      end do

      do i=1,dim
        do j=1,dim
          if (dabs(dreal(mat(i,j))).lt.1d-10) then
            mat(i,j) = dcmplx(0.d0,dimag(mat(i,j)))
          elseif(dabs(dimag(mat(i,j))).lt.1d-10) then
            mat(i,j) = dcmplx(dreal(mat(i,j)),0.d0)
          end if
        end do
      end do
C  Eigenvalue solver for a complex, non-symmetric matrix.

      call ZHEEV('v','l', dim, mat, dim, en, work, lwork, rwork, info)

!      call ZGEEV('n','v', dim, mat, dim, En, vl, dim, vr, dim,
!     &                           WORK, LWORK, RWORK, INFO )
*********************************************************************
***     Sort the eigenvalues in ascending order
*********************************************************************

      if (number.eq.0) then
        write(*,*)
        write(*,matrix) dreal(transpose(mat))
        write(*,*)
      end if

      end subroutine hbuild

***************************************************************************

      subroutine interHam(v,U,ham)
      implicit none

      real(8), intent(in) :: U, v(dim*2)
      complex(8),intent(out) :: ham(intd,intd)
      integer,parameter :: lwork=600

      integer :: i,j,k,jj,info
      real(8) :: rwork(100),cn(intd)
      complex(8) :: work(lwork)
!      complex(8) :: dum,vl(6,6),vr(6,6),cn(6)
      complex(8) :: vec(3),Bp(sites),Bm(sites)

      ham = zero

      do i=1,2
        ham(i,i+1) = -dsqrt(2.d0)*t
        ham(i+1,i) = -dsqrt(2.d0)*t
      end do

      do i=1,sites
        Bp(i) = (v(sites+i)+ione*v(sites*2+i))/dsqrt(2.d0)
        Bm(i) = (v(sites+i)-ione*v(sites*2+i))/dsqrt(2.d0)
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
        ham(2,j+3) = vec(j)
      end do

      k=0
      do i=1,3,2
        k=k+1
        vec(i) = Ba(mod(k,2)+1)*(-1.d0)**(k)
     &                   + Bb(mod(k,2)+1)*(-1.d0)**(k+1)
      end do

      do j=1,3
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

        end do
      end do

      do i=4,5
        ham(i,i+1) = Ba(2) + Bb(2)
        ham(i+1,i) = Ba(1) + Bb(1)
      end do

      ham(5,5) = ham(5,5) - (v(7)+v(8))

!      call ZGEEV('n','v',6, ham, 6, cn, vl, 6, vr, 6,
!     &                           WORK, LWORK, RWORK, INFO )

      call ZHEEV('v','l', 6, ham, 6, cn, work, lwork, rwork, info)
*********************************************************************
***     Sort the eigenvalues in ascending order
*********************************************************************
!      do I=1,6
!        do J=i+1,6
!          IF (dreal(cn(i)).GE.dreal(cn(J))) THEN
!             DUM = cn(I)
!             cn(I) = cn(J)
!             cn(J) = DUM
!
!            do JJ=1,6
!              DUM = vr(JJ,I)
!              vr(JJ,I) = vr(JJ,J)
!              vr(JJ,J) = DUM
!            end do
!          end if
!        end do
!      end do
!
!      ham = vr

      end subroutine interHam

***************************************************************************
      subroutine Pauli(sigma)
      implicit none

      integer :: i,j
      complex(8), intent(out) :: sigma(spin*4,spin*4)
***************************************************************************
***   Subroutine Pauli generates the well-known Pauli spin matrices or the
***   matrices forming the basis of the SU(2) space. The matrices are stored
***   block-diagonally, with sigma_0 starting at [1,1]. All off block-diagonal
***   elements are set to zero.
***************************************************************************
      sigma = zero
      do i=1,spin
        do j=1,spin
          if (i.eq.j) then
            sigma(i,j) = one
            sigma(i+3*spin,j+3*spin) = (-one)**(i-1)
          else
            sigma(i+spin,j+spin) = one
            sigma(i+2*spin,j+2*spin) = ione*(-1.d0)**(i)
          end if
        end do
      end do

      end subroutine Pauli

****************************************************************************
****  Numerical Recipes subroutines and functions                       ****
****************************************************************************
        SUBROUTINE frprmn(p,n,ftol,iter,fret)
        INTEGER iter,n,NMAX,ITMAX
        REAL (8) :: fret,ftol,p(n),EPS
        PARAMETER (NMAX=50,ITMAX=100000,EPS=1.d-12)
cU    USES dfunc,func,linmin
        INTEGER its,j
        REAL (8) :: dgg,fp,gam,gg,g(NMAX),h(NMAX),xi(NMAX)

        fp=func(p)
        call dfunc(p,xi)

        do 11 j=1,n
          g(j)=-xi(j)
          h(j)=g(j)
          xi(j)=h(j)
11    continue
        do 14 its=1,ITMAX
          iter=its
          call linmin(p,xi,n,fret)
          if(2.d0*dabs(fret-fp).le.ftol*(dabs(fret)+dabs(fp)+EPS))return
          fp=func(p)
          call dfunc(p,xi)
          gg=0.d0
          dgg=0.d0
          do 12 j=1,n
            gg=gg+g(j)**2
C         dgg=dgg+xi(j)**2
            dgg=dgg+(xi(j)+g(j))*xi(j)
12      continue
          if(gg.eq.0.d0)return
          gam=dgg/gg
          do 13 j=1,n
            g(j)=-xi(j)
            h(j)=g(j)+gam*h(j)
            xi(j)=h(j)
13      continue
14    continue
        pause 'frprmn maximum iterations exceeded'
        return
      END subroutine
C  (C) Copr. 1986-92 Numerical Recipes Software 6?6>)AY.
****************************************************************************
      FUNCTION f1dim(x)
      INTEGER NMAX
      REAL (8) :: x, f1dim
      PARAMETER (NMAX=50)
C     USES func
      INTEGER j,ncom
      REAL (8) :: pcom(NMAX),xicom(NMAX),xt(NMAX)
      COMMON /f1com/ pcom,xicom,ncom
      do 11 j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
11    continue
      f1dim=func(xt)
      return
      END function
C  (C) Copr. 1986-92 Numerical Recipes Software 6?6>)AY.
*****************************************************************************
        SUBROUTINE linmin(p,xi,n,fret)
        INTEGER n,NMAX
        REAL (8) :: fret,p(n),xi(n),TOL
        PARAMETER (NMAX=50,TOL=1.d-4)
CU    USES brent,f1dim,mnbrak
        INTEGER j,ncom
        REAL (8) :: ax,bx,fa,fb,fx,xmin,xx,pcom(NMAX),xicom(NMAX)
        COMMON /f1com/ pcom,xicom,ncom
        !real (8), EXTERNAL :: f1dim
        ncom=n
        do 11 j=1,n
          pcom(j)=p(j)
          xicom(j)=xi(j)
11    continue
        ax=0.d0
        xx=1.d0
        call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
        fret=brent(ax,xx,bx,f1dim,TOL,xmin)
        do 12 j=1,n
          xi(j)=xmin*xi(j)
          p(j)=p(j)+xi(j)
12    continue
        return
        END subroutine
C  (C) Copr. 1986-92 Numerical Recipes Software 6?6>)AY.
****************************************************************************

        SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
        REAL (8) :: ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,TINY
        EXTERNAL func
        PARAMETER (GOLD=1.618034d0, GLIMIT=100.d0, TINY=1.d-20)
        REAL (8) :: dum,fu,q,r,u,ulim
        fa=func(ax)
        fb=func(bx)
        if(fb.gt.fa)then
          dum=ax
          ax=bx
          bx=dum
          dum=fb
          fb=fa
          fa=dum
        endif
        cx=bx+GOLD*(bx-ax)
        fc=func(cx)
1     if(fb.ge.fc)then
          r=(bx-ax)*(fb-fc)
          q=(bx-cx)*(fb-fa)
          u=bx-((bx-cx)*q-(bx-ax)*r)
     &                      /(2.d0*sign(max(dabs(q-r),TINY),q-r))
          ulim=bx+GLIMIT*(cx-bx)
          if((bx-u)*(u-cx).gt.0.d0)then
            fu=func(u)
            if(fu.lt.fc)then
              ax=bx
              fa=fb
              bx=u
              fb=fu
              return
            else if(fu.gt.fb)then
              cx=u
              fc=fu
              return
            endif
            u=cx+GOLD*(cx-bx)
            fu=func(u)
          else if((cx-u)*(u-ulim).gt.0.d0)then
            fu=func(u)
            if(fu.lt.fc)then
              bx=cx
              cx=u
              u=cx+GOLD*(cx-bx)
              fb=fc
              fc=fu
              fu=func(u)
            endif
          else if((u-ulim)*(ulim-cx).ge.0.d0)then
            u=ulim
            fu=func(u)
          else
            u=cx+GOLD*(cx-bx)
            fu=func(u)
          endif
          ax=bx
          bx=cx
          cx=u
          fa=fb
          fb=fc
          fc=fu
          goto 1
        endif
        return
        END subroutine
****************************************************************************

        FUNCTION brent(ax,bx,cx,f,tol,xmin)
        INTEGER ITMAX
        REAL (8) :: brent,ax,bx,cx,tol,xmin,f,CGOLD,ZEPS
        EXTERNAL f
        PARAMETER (ITMAX=100,CGOLD=.3819660d0,ZEPS=1.0d-10)
        INTEGER iter
        REAL (8) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
        a=min(ax,cx)
        b=max(ax,cx)
        v=bx
        w=v
        x=v
        e=0.
        fx=f(x)
        fv=fx
        fw=fx
        do 11 iter=1,ITMAX
          xm=0.5d0*(a+b)
          tol1=tol*dabs(x)+ZEPS
          tol2=2.d0*tol1
          if(dabs(x-xm).le.(tol2-.5d0*(b-a))) goto 3
          if(dabs(e).gt.tol1) then
            r=(x-w)*(fx-fv)
            q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r
            q=2.d0*(q-r)
            if(q.gt.0.d0) p=-p
            q=dabs(q)
            etemp=e
            e=d
      if(dabs(p).ge.dabs(.5d0*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x))
     *goto 1
            d=p/q
            u=x+d
            if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
            goto 2
          endif
1       if(x.ge.xm) then
            e=a-x
          else
            e=b-x
          endif
          d=CGOLD*e
2       if(dabs(d).ge.tol1) then
            u=x+d
          else
            u=x+sign(tol1,d)
          endif
          fu=f(u)
          if(fu.le.fx) then
            if(u.ge.x) then
              a=x
            else
              b=x
            endif
            v=w
            fv=fw
            w=x
            fw=fx
            x=u
            fx=fu
          else
            if(u.lt.x) then
              a=u
            else
              b=u
            endif
            if(fu.le.fw .or. w.eq.x) then
              v=w
              fv=fw
              w=u
              fw=fu
            else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
              v=u
              fv=fu
            endif
          endif
11    continue
        pause 'brent exceed maximum iterations'
3     xmin=x
        brent=fx
        return
      END function
C  (C) Copr. 1986-92 Numerical Recipes Software 6?6>)AY.

*****************************************************************************

      FUNCTION dbrent(ax,bx,cx,f,df,tol,xmin)
      INTEGER ITMAX
      REAL(8) dbrent,ax,bx,cx,tol,xmin,df,f,ZEPS
      EXTERNAL df,f
      PARAMETER (ITMAX=100,ZEPS=1.0d-10)
      INTEGER iter
      REAL(8) a,b,d,d1,d2,du,dv,dw,dx,e
      REAL(8) fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm
      LOGICAL ok1,ok2
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.d0
      fx=f(x)
      fv=fx
      fw=fx
      dx=df(x)
      dv=dx
      dw=dx
      do 11 iter=1,ITMAX
        xm=0.5d0*(a+b)
        tol1=tol*dabs(x)+ZEPS
        tol2=2.d0*tol1
        if(dabs(x-xm).le.(tol2-.5d0*(b-a))) goto 3
        if(dabs(e).gt.tol1) then
          d1=2.d0*(b-a)
          d2=d1
          if(dw.ne.dx) d1=(w-x)*dx/(dx-dw)
          if(dv.ne.dx) d2=(v-x)*dx/(dx-dv)
          u1=x+d1
          u2=x+d2
          ok1=((a-u1)*(u1-b).gt.0.d0).and.(dx*d1.le.0.d0)
          ok2=((a-u2)*(u2-b).gt.0.d0).and.(dx*d2.le.0.d0)
          olde=e
          e=d
          if(.not.(ok1.or.ok2))then
            goto 1
          else if (ok1.and.ok2)then
            if(dabs(d1).lt.dabs(d2))then
              d=d1
            else
              d=d2
            endif
          else if (ok1)then
            d=d1
          else
            d=d2
          endif
          if(dabs(d).gt.dabs(0.5d0*olde))goto 1
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(dx.ge.0.d0) then
          e=a-x
        else
          e=b-x
        endif
        d=0.5d0*e
2       if(dabs(d).ge.tol1) then
          u=x+d
          fu=f(u)
        else
          u=x+sign(tol1,d)
          fu=f(u)
          if(fu.gt.fx)goto 3
        endif
        du=df(u)
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          dv=dw
          w=x
          fw=fx
          dw=dx
          x=u
          fx=fu
          dx=du
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            dv=dw
            w=u
            fw=fu
            dw=du
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
            dv=du
          endif
        endif
11    continue
      pause 'dbrent exceeded maximum iterations'
3     xmin=x
      dbrent=fx
      return
      END Function
C  (C) Copr. 1986-92 Numerical Recipes Software 6?6>)AY.


*************************************************************************
      FUNCTION df1dim(x)
      INTEGER NMAX
      REAL(8) df1dim,x
      PARAMETER (NMAX=50)
CU    USES dfunc
      INTEGER j,ncom
      REAL(8) df(NMAX),pcom(NMAX),xicom(NMAX),xt(NMAX)
      COMMON /f1com/ pcom,xicom,ncom
      do 11 j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
11    continue
      call dfunc(xt,df)
      df1dim=0.d0
      do 12 j=1,ncom
        df1dim=df1dim+df(j)*xicom(j)
12    continue
      return
      END Function
C  (C) Copr. 1986-92 Numerical Recipes Software 6?6>)AY.


      End Module
