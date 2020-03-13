      Module coplanarCG

      IMPLICIT NONE

      integer,parameter :: spin=2, sites=4, occ=2
      integer, parameter :: dim=spin*sites, intd=2*sites**2-sites
      real(kind=8), parameter :: t=0.5d0

      integer :: number,rest,xc,corr,flag
      real(8) :: ntarget(dim*2),En(dim)
      complex(8) :: hmat(dim,dim)
!      real(8) :: tBx(sites), tBy(sites), tBz(sites)
      complex(8) :: sig(spin*4,spin*4)
      complex(8), parameter :: Zero = (0.d0,0.d0), One = (1.d0,0.d0)
      complex(8), parameter :: IOne = (0.d0,1.d0)

      character(20) :: intprint,dmat
      character (len=30) :: matrix, tprint, vector

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

      v(2) = 0.d0

!      if (flag.eq.1) then
!        v(sites+2) = 0.25d0
!      end if
C  Recalculating the eigenvectors through exact diagonalization of the Hamiltonian.
      call hbuild(v,hmat)
      dens = 0.d0
      call densvec(dens,hmat)

      do i=1,dim
        do j=1,dim
          if (dabs(dreal(hmat(i,j))).lt.1d-10) then
            hmat(i,j) = dcmplx(0.d0,dimag(hmat(i,j)))
          elseif (dabs(dimag(hmat(i,j))).lt.1d-10) then
            hmat(i,j) = dcmplx(dreal(hmat(i,j)),0.d0)
          end if
        end do
      end do

C  Calculating the difference integral between the target and current densities.
      integral = 0.d0
      do i =1,dim*2
!        integral = integral + dreal((dens(i)-ntarget(i))
!     &                          *(conjg(dens(i)-ntarget(i))))

        if (i.gt.2*sites.or.i.le.3*sites) then
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
      complex(8) :: vec(8,occ*sites**2)

      call vmat(vec)
      call densvec(dens,hmat)

      dSdV = 0.d0
      dn = zero

      number = number + 1

***   Beginning of the loop to calculate the integral of dS/dV.
***   The subroutine dnvec generates the gradient vector of the density or
***   magnetization with respect to each potential on each site.

      counter = 0
      do num=1,4
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

!      write(*,vector) dSdV
!      write(*,*) '^^^^^^^^^^^^ dSdV ^^^^^^^^^^^^'

!      if (number.eq.1) then
!        call exit(-1)
!      end if

      end subroutine

***************************************************************************
      subroutine vmat(mat)
      implicit none
      complex(8),intent(out) :: mat(8,(occ*sites**2))
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
              do num=1,4
                mat(4*(s-1)+num,counter) =
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

      end do

      end function

***************************************************************************

      subroutine dnvec(num,k,dn,vec)
      implicit none

      integer :: i,j,k,num,alpha,counter,dum
      complex(8) :: x,y
      complex(8),intent(in) :: vec(8,occ*sites**2)
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

      i = (k-1)*(occ*sites)

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
      complex(8),intent(in) :: inmat(intd,intd)

      integer :: i,ii,j,jj,k,kk,nsng,ntrp,x,y,z
      real(8) :: r
      real(8) :: mx(sites),my(sites),mz(sites),n(sites)
      complex(8) :: coef(intd),phi(sites,sites),nud(sites)

      n0 = 0.d0
      n = 0.d0
      mx = 0.d0
      my = 0.d0
      mz = 0.d0
      nud = zero

      x = sites
      y = sites*2
      z = sites*3

      do i=1,intd
        coef(i) = inmat(i,1)
      end do

      DO I=1,sites
        DO J=1,sites
          IF (I.EQ.J) then
            PHI(I,J) = 1.D0
          else
            phi(i,j) = 0.d0
          end if
        end do
      end do

      NSng = sites*(sites+1)/2
      NTrp = 3*sites*(sites-1)/2

      DO I=1,sites
C**------------------------------------------------------------
C**   Singlet block contributions to N
C**------------------------------------------------------------
        N(I) = 0.D0

        II = -sites
        DO J=1,sites
          II = II + sites+2-J
          N(I) = N(I) + dreal(PHI(J,I)*CDABS(Coef(II))**2 )
        end do

        II = -sites
        DO J=1,sites-1
          II = II + sites+2-J
          DO K=1,sites-J
            JJ = II + K
            N(I) = N(I)+0.5D0*
     &             dreal((PHI(J,I)+PHI(J+K,I))*CDABS(Coef(JJ))**2)
          end do
        end do
C**------------------------------------------------------------
C**   Triplet block contributions to N
C**------------------------------------------------------------

        JJ = NSng-2
        DO J=1,sites-1
          DO K=1,sites-J
            JJ = JJ+3
            N(I) = N(I)+0.5D0* dreal(
     &                  (PHI(J,I)+PHI(J+K,I))*(CDABS(Coef(JJ))**2
     &                 + CDABS(Coef(JJ+1))**2+CDABS(Coef(JJ+2))**2)
     &                          )
          end do
        end do
C**------------------------------------------------------------
C**   Triplet block contributions to MZ
C**------------------------------------------------------------

        JJ = nsng-2
        DO J=1,sites-1
          DO K=1,sites-J
            JJ = JJ+3
            MZ(I) = MZ(I)+0.5D0*dreal(
     &                      (PHI(J,I)+PHI(J+K,I))
     &            *(CDABS(Coef(JJ))**2 - CDABS(Coef(JJ+2))**2)
     &                            )
          end do
        end do
C**------------------------------------------------------------
C**   Singlet-Triplet block contributions to MZ
C**------------------------------------------------------------

        II = -sites
        KK = NSng-1
        DO J=1,sites-1
          II = II + sites+2-J
          DO K=1,sites-J
            JJ = II+K
            KK = KK+3
            MZ(I) = MZ(I)+0.5D0* dreal(
     &              (PHI(J,I)-PHI(J+K,I))
     &        *(Coef(JJ)*DCONJG(Coef(KK))
     &                      + DCONJG(Coef(JJ))*Coef(KK)) )
          end do
        end do
C**------------------------------------------------------------
C**   Singlet-Triplet block contributions to NUD
C**------------------------------------------------------------

        NUD(I) = (0.D0,0.D0)

        II = -sites
        KK = NSng-1
        DO J=1,sites-1
          II = II + sites+2-J
          DO K=1,sites-J
            JJ = II+K
            KK = KK+3
            NUD(I) = NUD(I) +
     &                  0.25D0*DSQRT(2.D0)*(PHI(J,I)-PHI(J+K,I))
     &         *(Coef(JJ)*DCONJG(Coef(KK+1))
     &                  - DCONJG(Coef(JJ))*Coef(KK-1))
          end do
        end do
C**------------------------------------------------------------
C**   Triplet block contributions to NUD
C**------------------------------------------------------------

        KK = NSng-1
        DO J=1,sites-1
          II = II + sites+2-J
          DO K=1,sites-J
            KK = KK+3
            NUD(I) = NUD(I) +
     &              0.25D0*DSQRT(2.D0)*(PHI(J,I)+PHI(J+K,I))
     &         *(Coef(KK)*DCONJG(Coef(KK+1))
     &                + DCONJG(Coef(KK))*Coef(KK-1) )
          end do
        end do

        MX(I) = 2.D0*DREAL(NUD(I))
        MY(I) = -2.D0*DIMAG(NUD(I))

        N(I) = 2.D0*N(I)
        MX(I) = 2.D0*MX(I)
        MY(I) = 2.D0*MY(I)
        MZ(I) = 2.D0*MZ(I)

      end do

      do i = 1, sites
        n0(i) = n(i)
        n0(x+i) = mx(i)
        n0(y+i) = my(i)
        n0(z+i) = mz(i)
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

      end subroutine hbuild

****************************************************************************

      subroutine interHam(v,U0,U1,ham,cn)
      implicit none

      real(8), intent(in) :: v(dim*2),U0,U1
      complex(8),intent(out) :: ham(intd,intd)
      integer,parameter :: lwork=600

      integer :: Ntrp,Nsng,i,j,k,l,ii,jj,info
      real(8) :: rwork(100),cn(intd)
      complex(8) :: work(lwork)
!      complex(8) :: dum,vl(6,6),vr(6,6),cn(6)
      complex(8) :: Bp(sites),Bm(sites)

      ham = ZERO

      DO I=1,sites
        BP(I) = (v(sites+i)*ONE + v(sites*2+I)*IONE)/DSQRT(2.D0)
        BM(I) = (v(sites+i)*ONE - v(sites*2+I)*IONE)/DSQRT(2.D0)
      end do

      NSng = sites*(sites+1)/2
      NTrp = 3*sites*(sites-1)/2

C**------------------------------------------------------------
C**   Singlet block
C**------------------------------------------------------------

C**   First the diagonal elements

      II = -sites
      DO I=1,sites
         II = II + sites+2-I
         ham(II,II) = 2.D0*V(I) + U0
      end do

      II = -sites+1
      DO I=1,sites-1
         II = II + sites+2-I
         ham(II,II) = V(I) + V(I+1) + U1
      end do

      IF (sites.GT.2) THEN
        JJ=-sites+2
        DO J=1,sites-2
          JJ = JJ + sites+2-J
          DO I=1,sites-J-1
            ham(JJ+I-1,JJ+I-1) = V(J) + V(J+I+1)
          end do
        end do
      ENDIF

C**   Now the off-diagonal elements

      II = -sites
      DO I=1,sites-1
        II = II + sites+2-I
        ham(II,II+1) = -DSQRT(2.D0)*T
        ham(II+1,II) = -DSQRT(2.D0)*T
      end do

      II = -sites+1
      DO I=1,sites-1
        II = II + sites+2-I
        ham(II,II+sites-I) = -DSQRT(2.D0)*T
        ham(II+sites-I,II) = -DSQRT(2.D0)*T
      end do

      IF (sites.GT.2) THEN
        JJ=-sites+2
        DO J=1,sites-2
          JJ = JJ + sites+2-J
          DO I=1,sites-J-1
            ham(JJ+I-2,JJ+I-1) = -T
            ham(JJ+I-1,JJ+I-2) = -T

            ham(JJ+I-1,JJ+I-1+(sites-J)) = -T
            ham(JJ+I-1+(sites-J),JJ+I-1) = -T
          end do
        end do
      ENDIF

C**------------------------------------------------------------
C**   Singlet-Triplet and Triplet-Singlet blocks
C**------------------------------------------------------------

      II = -sites+1
      JJ = NSng+1
      DO J=1,sites-1
         II = II + sites-J+2
         DO I=1,sites-J

            ham(II+I-1,JJ) = -BP(J) + BP(J+I)
            ham(II+I-1,JJ+1) = v(sites*3+j) - v(sites*3+J+I)
            ham(II+I-1,JJ+2) = BM(J) - BM(J+I)

            ham(JJ,II+I-1) = -BM(J) + BM(J+I)
            ham(JJ+1,II+I-1) = v(sites*3+j) - v(sites*3+J+I)
            ham(JJ+2,II+I-1) = BP(J) - BP(J+I)

            JJ=JJ+3
         end do
      end do
C**------------------------------------------------------------
C**   Triplet block
C**------------------------------------------------------------

      II = Nsng-1
      DO I=1,sites-1
         DO J=I+1,sites
            II = II + 3
            ham(II-1,II-1) = V(I) + V(J)
     &                          + v(sites*3+I) + v(sites*3+J)
            ham(II,II) =  V(I) + V(J)
            ham(II+1,II+1) = V(I) + V(J)
     &                          - v(sites*3+I) - v(sites*3+J)

            IF ((J-I).EQ.1) THEN
               ham(II-1,II-1) = ham(II-1,II-1) + U1
               ham(II,II) = ham(II,II) + U1
               ham(II+1,II+1) = ham(II+1,II+1) + U1
            ENDIF

            ham(II-1,II) = BM(I) + BM(J)
            ham(II,II-1) = BP(I) + BP(J)
            ham(II,II+1) = BM(I) + BM(J)
            ham(II+1,II) = BP(I) + BP(J)
         end do
      end do

      II = NSng-2
      DO I=1,sites-1
         DO J=I+1,sites
            II = II+3

            JJ = NSng-2
            DO K=1,sites-1
               DO L=K+1,sites
                  JJ = JJ+3
                  IF ((I.EQ.K).AND.(IABS(J-L).EQ.1).OR.
     &                (J.EQ.L).AND.(IABS(I-K).EQ.1)) THEN
                     ham(II,JJ) = -T
                     ham(II+1,JJ+1) = -T
                     ham(II+2,JJ+2) = -T
                  ENDIF
               end do
            end do

         end do
      end do

      call ZHEEV('v','l', intd, ham, intd, cn, work, lwork, rwork, info)

      end subroutine
***************************************************************************

      subroutine btrans(den,v,v0,mm,btx,bty,btz)
      implicit none

      real(8),intent(in) :: den(dim*2),v(dim*2),v0(dim*2)
      real(8),intent(out) :: mm(sites),btx(sites),bty(sites),
     &                       btz(sites)

      integer :: I,x,y,z

      real(8) :: bxcx(sites),bxcy(sites),bxcz(sites),
     &           mx(sites),my(sites),mz(sites),torq(3),
     &           tx(sites),ty(sites),tz(sites)

      x = sites
      y = 2*sites
      z = 3*sites

      torq = 0.d0

      bxcx = 0.d0
      bxcy = 0.d0
      bxcz = 0.d0

      do i = 1,sites
        bxcx(i) = v(x+i) - v0(x+i)
        bxcy(i) = v(y+i) - v0(y+i)
        bxcz(i) = v(z+i) - v0(z+i)
        mx(i) = den(x+i)
        my(i) = den(y+i)
        mz(i) = den(z+i)
      end do

      do i = 1,sites
        mm(i) = dsqrt(mx(i)**2 + my(i)**2 + mz(i)**2)
      end do

      do i = 1, sites
        tx(i) = my(i)*BXCZ(I) - mz(i)*BXCY(I)
        ty(i) = mz(i)*BXCX(I) - mx(i)*BXCZ(I)
        tz(i) = mx(i)*BXCY(I) - my(i)*BXCX(I)
      end do

      do i = 1, sites
        btx(i) = (ty(i)*mz(i) - tz(i)*my(i))/mm(i)**2
        bty(i) = (tz(i)*mx(i) - tx(i)*mz(i))/mm(i)**2
        btz(i) = (tx(i)*my(i) - ty(i)*mx(i))/mm(i)**2
      end do

      end subroutine

***********************************************************************

      subroutine blong(den,v,v0,mm,blx,bly,blz)
      implicit none

      real(8),intent(in) :: den(dim*2),v(dim*2),v0(dim*2)
      real(8),intent(out) :: mm(sites),blx(sites),bly(sites),
     &                       blz(sites)

      integer :: I,x,y,z

      real(8) :: bxcx(sites),bxcy(sites),bxcz(sites),
     &           mx(sites),my(sites),mz(sites),mdotb(sites)

      x = sites
      y = 2*sites
      z = 3*sites

      bxcx = 0.d0
      bxcy = 0.d0
      bxcz = 0.d0

      do i = 1,sites
        bxcx(i) = v(x+i) - v0(x+i)
        bxcy(i) = v(y+i) - v0(y+i)
        bxcz(i) = v(z+i) - v0(z+i)
        mx(i) = den(x+i)
        my(i) = den(y+i)
        mz(i) = den(z+i)
      end do

      do i = 1,sites
        mm(i) = dsqrt(mx(i)**2 + my(i)**2 + mz(i)**2)
      end do

      do i = 1, sites
        mdotb(i) = mx(i)*bxcx(i) + my(i)*BXCY(I) + mz(i)*bxcz(i)
      end do

      do i = 1, sites
        blx(i) = mdotb(i)*mx(i)/mm(i)**2
        bly(i) = mdotb(i)*my(i)/mm(i)**2
        blz(i) = mdotb(i)*mz(i)/mm(i)**2
      end do

      end subroutine

***************************************************************************

      subroutine bpar(den,v,v0,mm,blx,bly,blz)
      implicit none

      real(8),intent(in) :: den(dim*2),v(dim*2),v0(dim*2)
      real(8),intent(out) :: mm(sites),blx(sites),bly(sites),
     &                       blz(sites)

      integer :: I,x,y,z

      real(8) :: bxcx(sites),bxcy(sites),bxcz(sites),
     &           mx(sites),my(sites),mz(sites),bmag(sites)

      x = sites
      y = 2*sites
      z = 3*sites

      bxcx = 0.d0
      bxcy = 0.d0
      bxcz = 0.d0

      do i = 1,sites
        bxcx(i) = v(x+i) - v0(x+i)
        bxcy(i) = v(y+i) - v0(y+i)
        bxcz(i) = v(z+i) - v0(z+i)
        mx(i) = den(x+i)
        my(i) = den(y+i)
        mz(i) = den(z+i)
      end do

      do i = 1,sites
        mm(i) = dsqrt(mx(i)**2 + my(i)**2 + mz(i)**2)
      end do

      do i = 1, sites
        bmag(i) = dsqrt(bxcx(i)**2 + BXCY(I)**2 + bxcz(i)**2)
      end do

      do i = 1, sites
        blx(i) = bmag(i)*mx(i)/mm(i)**2
        bly(i) = bmag(i)*my(i)/mm(i)**2
        blz(i) = bmag(i)*mz(i)/mm(i)**2
      end do

      end subroutine

***************************************************************************

      subroutine tcalc(den,v,v0,torq,tx,ty,tz)
      implicit none

      real(8),intent(in) :: den(dim*2),v(dim*2),v0(dim*2)
      real(8),intent(out) :: torq(3),tx(sites),ty(sites),
     &                       tz(sites)

      integer :: I,x,y,z

      real(8) :: bxcx(sites),bxcy(sites),bxcz(sites),
     &           mx(sites),my(sites),mz(sites)

      x = sites
      y = 2*sites
      z = 3*sites

      torq = 0.d0

      bxcx = 0.d0
      bxcy = 0.d0
      bxcz = 0.d0

      do i = 1,sites
        bxcx(i) = v(x+i) - v0(x+i)
        bxcy(i) = v(y+i) - v0(y+i)
        bxcz(i) = v(z+i) - v0(z+i)
        mx(i) = den(x+i)
        my(i) = den(y+i)
        mz(i) = den(z+i)
      end do

      do i = 1, sites
        tx(i) = (my(i)*BXCZ(I) - mz(i)*BXCY(I))
        ty(i) = (mz(i)*BXCX(I) - mx(i)*BXCZ(I))
        tz(i) = (mx(i)*BXCY(I) - my(i)*BXCX(I))
      end do

      do i = 1, sites
        torq(1) = torq(1) + tx(i)
        torq(2) = torq(2) + ty(i)
        torq(3) = torq(3) + tz(i)
      end do

      END subroutine

***************************************************************************

      subroutine tortc(den, v0, textx, texty, textz)
      implicit none

      real(8),intent(in) :: den(dim*2), v0(dim*2)
      real(8),intent(out) :: textx(sites),texty(sites),
     &                       textz(sites)

      integer :: I,x,y,z

      real(8) :: mx(sites), my(sites), mz(sites)

      x = sites
      y = 2*sites
      z = 3*sites

      do i = 1,sites
        mx(i) = den(x+i)
        my(i) = den(y+i)
        mz(i) = den(z+i)
      end do

      do i = 1, sites
        textx(i) = (my(i)*v0(z+i) - mz(i)*v0(y+i))
        texty(i) = (mz(i)*v0(x+i) - mx(i)*v0(z+i))
        textz(i) = (mx(i)*v0(y+i) - my(i)*v0(x+i))
      end do

      END subroutine

***************************************************************************

      subroutine metric(bmet,btx,bty,btz,v,v0)
      implicit none

      real(8),intent(in) :: btx(sites),bty(sites),btz(sites),
     &                      v(dim*2),v0(dim*2)
      real(8),intent(out) :: bmet

      real(8) :: bxcx(sites),bxcy(sites),bxcz(sites),s1,s2,
     &           btm(sites),bxcm(sites)
      integer :: i,x,y,z

      x = sites
      y = 2*sites
      z = 3*sites

      do i = 1, sites
        bxcx(i) = v(x+i) - v0(x+i)
        bxcy(i) = v(y+i) - v0(y+i)
        bxcz(i) = v(z+i) - v0(z+i)
      end do

      do i = 1,sites
        btm(i) = dsqrt(btx(i)**2+bty(i)**2+btz(i)**2)
        bxcm(i) = dsqrt(bxcx(i)**2+bxcy(i)**2+bxcz(i)**2)
      end do

      s1=0.d0
      s2=0.d0
      do i = 1, sites
        s1 = s1 + btm(i)
        s2 = s2 + bxcm(i)
      end do

      bmet = s1/s2

      end subroutine

*************************************************************

      subroutine CalcAngle(tta, exden, ksden)
      implicit none

      real(8), intent(in) :: exden(dim*2), ksden(dim*2)
      real(8), intent(out) :: tta(sites)

      integer :: i, x, y, z
      real(8) :: exmm(sites), ksmm(sites), xmx(sites),
     &           xmy(sites), xmz(sites), kmx(sites), kmy(sites),
     &           kmz(sites), mdot(sites)

      x = sites
      y = sites*2
      z = sites*3

      do i = 1, sites
        xmx(i) = exden(x + i)
        xmy(i) = exden(y + i)
        xmz(i) = exden(z + i)
        kmx(i) = ksden(x + i)
        kmy(i) = ksden(y + i)
        kmz(i) = ksden(z + i)
      end do

      do i = 1, sites
        exmm(i) = dsqrt(xmx(i)**2 + xmy(i)**2 + xmz(i)**2)
        ksmm(i) = dsqrt(kmx(i)**2 + kmy(i)**2 + kmz(i)**2)
        mdot(i) = xmx(i)*kmx(i) + xmy(i)*kmy(i) + xmz(i)*kmz(i)
      end do

      do i = 1, sites
        tta(i) = ACOS(mdot(i)/(exmm(i)*ksmm(i)))
      end do

      end subroutine

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

      subroutine StepMin(steps, dx, v, ftol, iter, fret, fretlo)
      implicit none

      integer, intent(in) :: steps
      integer, intent(out) :: iter
      real(8), intent(in) :: dx
      real(8), intent(out) :: ftol
      real(8), intent(inout) :: v(dim*2), fret, fretlo

      integer :: i, j, pass
      real(8) :: vstart(dim*2), vprev(dim*2)

      do pass = 1, steps
        write(*,*) 'Step #:', pass
        do i = 1, dim*2
          v(i) = vstart(i) + steps*dx
          call frprmn(v, dim*2, ftol, iter, fret)
          if (fret<fretlo) then
            fretlo = fret
          end if
          if (fretlo<1d-14) return
        end do

        vprev = vstart
        v = vstart
        do i = 1, dim*2-1
          do j = i+1, dim*2
            v = vprev
            v(i) = v(i) + steps*dx
            v(j) = v(j) + steps*dx
            call frprmn(v, dim*2, ftol, iter, fret)
            if (fretlo<1d-14) return
          end do
          vprev(i) = vstart(i) + steps*dx
        end do
      end do

      end subroutine

****************************************************************************
****  Numerical Recipes subroutines and functions                       ****
****************************************************************************
        SUBROUTINE frprmn(p,n,ftol,iter,fret)
        INTEGER iter,n,NMAX,ITMAX
        REAL (8) :: fret,ftol,p(n),EPS
        PARAMETER (NMAX=50,ITMAX=1000000,EPS=1.d-10)
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

        if (iter.eq.itmax) then
          write(*,*) 'maximum iterations reached'
          write(*,*) iter
          return
        end if

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
