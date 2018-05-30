      Program gradient
      implicit none

      integer, parameter :: n=4
      real(kind=8), parameter :: dab=2.d0,dbc=-3.d0,dcd=1.5d0
      integer :: i,j,k,alpa,bta,counter,iter
      real (kind=8) :: ntarget(n),En(n),dSdV(n),derivative
      real (kind=8) :: diffn(n),InnerP(n,2),diffv(n-1)
      real(kind=8) :: v(n),xi(2),x,integral(n),dens(n),hmat(n,n)
      character (len=15) :: matrix!, matrix2
      write(matrix,'(a, i3, a)') '(', n, 'f15.12)'
!      write(matrix2,'(a, i3, a)') '(', 2, 'f8.3)'

C  Dfunc is fed the potential difference between the lattice points, so
C  v must be calculated to recalculate the eigenvectors of the system.
      v(1) = 0.d0
      v(2) = v(1) - dab
      v(3) = v(2) - dbc
      v(4) = v(3) - dcd
C  Do loop that calculates the perturbation step of the derivative subroutine.
C  Specifically, this routine calculates dn/dv, using perturbation theory:
C  \Sum < phi_m | vk | phi_i > * (E_i - E_m)^(-1) | phi_m >
      call hbuild(v,En,hmat,n)

      !write(*,matrix) transpose(hmat)

      ntarget = 0.d0
      do i=1,n
        do j=1,2
          ntarget(i) = ntarget(i) + 2.d0*(hmat(i,j)**2)
        end do
      end do

!      do i=2,n
!        v(1) = 0.d0
!        v(i) = dble(i)*1d-1
!      end do

      call hbuild(v,En,hmat,n)
      dens = 0.d0
      do i=1,n
        do j=1,2
          dens(i) = dens(i) + 2.d0*(hmat(i,j)**2)
        end do
      end do

      do j=1,n
        x = 0.d0
        do i=1,n
          do alpa = 1,2
            do bta = 1,n
!              write(*,*) derivative(v,hmat,En,alpa,bta,i,j)
              x = x + 4.d0*(dens(i)-ntarget(i))
     &                *derivative(v,hmat,En,alpa,bta,i,j)
            end do
          end do
        end do

        dSdV(j) = x
      end do

      write(*,matrix) dSdV
!      write(*,matrix) transpose(hmat)
!      write(*,*) '**************'
!      write(*,matrix) En

      end

***************************************************************************
      function derivative(v,hmat,En,alpha,beta,i,j) result(dn)
      implicit none

      integer,parameter :: n=4
      integer :: i,j,alpha,beta
      real(8) :: num, denom, dn, hmat(n,n), En(n), v(n)

      if (alpha.eq.beta) then
        dn = 0.d0
      else
        num = hmat(i,alpha)*hmat(i,beta)*hmat(j,alpha)*hmat(j,beta)
        denom = En(beta) - En(alpha)
        dn = num/denom
      end if

      end function

***************************************************************************

        subroutine hbuild(v,En,mat,n)
        implicit none

        integer, parameter :: lwork=136
        integer :: i,j,info,n
        real (kind=8), parameter :: t=0.5d0
        real (kind=8) :: v(n), x, mat(n,n)
        real(kind=8) :: work(lwork),En(n)

C  Setting up noninteracting Hubbard-like system Hamiltonian with t as the
C  hopping constant between lattice sites.
        do i = 1,n
          do j=1,n
            if (i.eq.j) Then
              mat(i,j) = v(i)
            elseif (abs(i-j).eq.1) Then
              mat(i,j) = -t
            else
              mat(i,j) = 0.d0
            endif
          end do
        end do

C  Eigenvalue solver for a symmetric matrix.
        call DSYEV('v','u',n,mat,n,En,work,lwork,info)

C  Normalization step to normalize each orbital to 1.
        do j = 1,n
          x = 0.d0
            do i = 1, n
              x = x + mat(i,j)**2
            end do
            do i = 1, n
              mat(i,j) = mat(i,j)/dsqrt(x)
            end do
        end do

        end subroutine hbuild

**************************************************************************
        Subroutine perturb(point,InnerP,v,en,hmat,n)

        implicit none

        integer :: i,j,k,point,n
        real (kind = 8) :: InnerP(n,2),dV(n,n),v(n)
        real (kind = 8) :: x,en(n),hmat(n,n)
        character (len=15) :: matrix,matrix2

        write(matrix,'(a, i3, a)') '(', n, 'f8.3)'
      write(matrix2,'(a, i3, a)') '(', 2, 'f8.3)'


!      do i=1,n
!         v(i) = dble(i)
!      end do

        do i=1,n
          do j=1,n
             if (i.eq.j) then
                dV(i,j) = v(i)
!                dV(j,i) = v(i)
             else
                dV(i,j) = 0.d0
             end if
          end do
        end do

!        write(*,*) 'iter inside perturb loop:'
!        write(*,*) point

C  Step calculating the inner product of each eigenvector on the perturbation
C  potential. The calculated values are stored in a 4 x 2 matrix, where the
C  row and column denote the bra and ket of the orbitals in the calculation,
C  respectively.

      do k=1,2
        x = 0.d0
        do j=1,n
          InnerP(j,k)=0.d0
          if (j.eq.k) then
             x = x + 0.d0
          else
             do i=1,n
                x = x + hmat(i,j)*dv(i,point)*hmat(i,k)/(En(k)-En(j))
             end do
          end if

          if (j.eq.k) then
             InnerP(j,k) = 0.d0
          else
             InnerP(j,k) = x
          end if
        end do
      end do

      write(*,matrix2) transpose(InnerP)

        end subroutine
