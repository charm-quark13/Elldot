      program dertest
      implicit none

      integer,parameter :: occ=2,sites=2,dim=2*sites
      complex(8),parameter:: zero=(0.d0,0.d0),one=(1.d0,0.d0)
      complex(8),parameter:: ione=(0.d0,1.d0)

      real(8) :: en(dim),rmat(dim,dim),imat(dim,dim)
      complex(8) :: mat(8,(occ+sites)*sites)
      integer :: i,j,k,s,num,alpha,counter
      complex(8) :: x,hmat(dim,dim)

      character(20) :: intprint,dmat
      character (len=30) :: matrix, matprint, vector

      write(matrix,'(a, i3, a)') '(', dim, 'f13.7)'
      write(matprint,'(a, i3, a)') '(', sites, 'e16.6)'
      write(vector,'(a, i3, a)') '(', 1, 'f16.10)'
      write(intprint,'(a, i3, a)') '(', 6, 'e16.6)'
      write(dmat,'(a,i3,a)') '(', dim*2,'f14.10)'

***   vmat generates a matrix containing the values of the derivatives of
***   psi_up and psi_dn for each occupied state with respect to each potential.
***   Please see the LaTeX write up for exact ordering of the matrix. The values
***   are used in the dnvec calculations to calculate the magnetization
***   derivatives.

      open(20,file='renergy.txt')
      read(20,vector) en
      close(20)

      write(*,vector) en

      open(120,file='hamre.txt')
      read(120,matrix) rmat
      close(120)

      open(220,file='hamim.txt')
      read(220,matrix) imat
      close(220)

      hmat = rmat + ione*imat

      do i=1,dim
        do j=1,dim
          write(*,*) hmat(i,j)
        end do
        write(*,*) '***'
      end do

      call exit(-1)
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

      contains

***************************************************************************
      function derivative(alpha,j,k,num,sigma) result(dp)
      implicit none

      integer :: l,m,beta
      integer, intent(in) :: j,k,alpha,sigma,num
      !real(8), intent(in) :: v(dim*2)
      complex(8) :: numer, denom
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

************************************************************************

      end program
