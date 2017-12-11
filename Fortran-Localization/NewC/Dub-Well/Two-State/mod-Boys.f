      MODULE Boys
      IMPLICIT NONE

!      INTEGER, PARAMETER :: NOCC=7,GRID=200

      INTEGER, PARAMETER :: NOCC=2,GRID=200
      INTEGER, PARAMETER :: NUN = 5
      INTEGER, PARAMETER :: PERT=NOCC*(NOCC-1)/2

      double precision, parameter :: pi=dacos(-1.d0)
      double precision, parameter :: DELTA=10.0D0/(grid-1)

      contains


********************************************************************************

      SUBROUTINE localize(phi)

      double precision :: a(nocc,nocc),b(nocc,nocc),d(nocc,nocc)
      double precision :: dmax,ak,bk,angle,dum1,dum2,x,y,z
      double precision :: reference(pert,3)
      double precision :: phi(grid,nocc),original(grid,nocc)
      double precision, parameter :: thresh=0.1d-3
      integer :: i,j,k,l,iter,counter
      integer, parameter :: itmax = 10000

      counter = 0

      do i=1,grid
         do j=1,nocc
            original(i,j) = phi(i,j)
         end do
      end do

      do iter = 1,itmax

      counter = counter + 1

      if ((iter.gt.3).and.(maxval(d).lt.thresh)) then
         exit
      else
   
         do j=1,nocc
            do k=1,nocc

!***********************************************************************
!          Here Kernels A, A2, and B correspond to the integrals   
!         calculated in Edmiston-Ruedenberg. Specifically, eqns 17
!        and 18, where they are calculated according to eqns 5 & 5' 
!***********************************************************************
               x = 0.d0
               do l=1,grid
                  dum1 = 0.d0
                  call integral(l,phi(l,j),phi(l,k),dum1)
                  x = x + dum1
               end do
           
               y = 0.d0
               do l=1,grid
                  dum1 = 0.d0
                  call integral(l,phi(l,j),phi(l,j),dum1)
                  y = y + dum1
               end do

               z = 0.d0
               do l=1,grid
                  dum1 = 0.d0
                  call integral(l,phi(l,k),phi(l,k),dum1)
                  z = z + dum1
               end do

               a(j,k) = x**2 - (y - z)**2/4.d0

               b(j,k) = x * (y-z)

               if (j.eq.k) then 
                  d(j,k) = 0.d0
               end if

               d(j,k) = a(j,k) + 
     &                 dsqrt(a(j,k)**2 + b(j,k)**2)

               if (j.eq.k) then
                  d(j,k) = 0.d0
               end if

      !         write(*,*) d(j,k)
            end do
         end do

         do j=1,nocc
            do k=1,nocc
               if (d(j,k).eq.maxval(d)) then
                  dum1 = dacos(-a(j,k)/dsqrt(a(j,k)**2+b(j,k)**2))/4.d0
                  call rotate(phi,dum1,j,k)
               end if
         !      do l = 1,grid
         !         x = original(l,j) - phi(l,j)
         !      end do
            end do
         end do
         !write(*,*) dum1
      end if
      end do

      do i=1,nocc
         do l=1,grid
            write(99+i,*)l,phi(l,j)
         end do
      end do

      write(*,*) counter

      End subroutine

*****************************************************************************

      SUBROUTINE rotate(phi,alpha,i,j)

      DOUBLE PRECISION :: alpha,phi(grid,nocc)
      DOUBLE PRECISION :: lambda_i(grid),lambda_j(grid)
 
      INTEGER :: i,j,r

      do r=1,grid
         lambda_i(r) = dcos(alpha)*phi(r,i) + dsin(alpha)*phi(r,j)
         lambda_j(r) = -dsin(alpha)*phi(r,i) + dcos(alpha)*phi(r,j)
      end do

      do r=1,grid
         phi(r,i) = lambda_i(r)
         phi(r,j) = lambda_j(r)
      end do

      end subroutine
****************************************************************************

      SUBROUTINE Integral(l,phi1,phi2,dum1) 

      INTEGER :: l
      DOUBLE PRECISION :: dum1
      DOUBLE PRECISION :: phi1,phi2

      DUM1 = dum1 + phi1*phi2 * delta*L 

      END SUBROUTINE

***************************************************************************

      SUBROUTINE jacobi(a,n,np,d,v,nrot)

      INTEGER n,np,nrot,NMAX

      double precision a(np,np),d(np),v(np,np)

      PARAMETER (NMAX=500)

      INTEGER i,ip,iq,j

      double precision c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)

      do 12 ip=1,n

        do 11 iq=1,n

          v(ip,iq)=0.d0

11      continue

        v(ip,ip)=1.d0

12    continue

      do 13 ip=1,n

        b(ip)=a(ip,ip)

        d(ip)=b(ip)

        z(ip)=0.d0

13    continue

      nrot=0

      do 24 i=1,50

        sm=0.d0

        do 15 ip=1,n-1

          do 14 iq=ip+1,n

            sm=sm+dabs(a(ip,iq))

14        continue

15      continue

        if(sm.eq.0.d0)return

        if(i.lt.4)then

!          tresh=0.2d0*sm/10.d10**2

           tresh=0.d0

        else

          tresh=0.d0

        endif

        do 22 ip=1,n-1

          do 21 iq=ip+1,n

            g=100.d0*dabs(a(ip,iq))

            if((i.gt.4).and.(dabs(d(ip))+g.eq.dabs(d(ip)))
     &                .and.(dabs(d(iq))+g.eq.dabs(d(iq))))then

              a(ip,iq)=0.d0

            else if(dabs(a(ip,iq)).gt.tresh)then

              h=d(iq)-d(ip)

              if(dabs(h)+g.eq.dabs(h))then

                t=a(ip,iq)/h

              else

                theta=0.5d0*h/a(ip,iq)

                t=1.d0/(dabs(theta)+dsqrt(1.d0+theta**2))

                if(theta.lt.0.d0)t=-t

              endif

              c=1.d0/dsqrt(1.d0+t**2)

              s=t*c

              tau=s/(1.d0+c)

              h=t*a(ip,iq)

              z(ip)=z(ip)-h

              z(iq)=z(iq)+h

              d(ip)=d(ip)-h

              d(iq)=d(iq)+h

              a(ip,iq)=0.d0

              do 16 j=1,ip-1

                g=a(j,ip)

                h=a(j,iq)

                a(j,ip)=g-s*(h+g*tau)

                a(j,iq)=h+s*(g-h*tau)

16            continue

              do 17 j=ip+1,iq-1

                g=a(ip,j)

                h=a(j,iq)

                a(ip,j)=g-s*(h+g*tau)

                a(j,iq)=h+s*(g-h*tau)

17            continue

              do 18 j=iq+1,n

                g=a(ip,j)

                h=a(iq,j)

                a(ip,j)=g-s*(h+g*tau)

                a(iq,j)=h+s*(g-h*tau)

18            continue

              do 19 j=1,n

                g=v(j,ip)

                h=v(j,iq)

                v(j,ip)=g-s*(h+g*tau)

                v(j,iq)=h+s*(g-h*tau)

19            continue

              nrot=nrot+1

            endif

21        continue

22      continue

        do 23 ip=1,n

          b(ip)=b(ip)+z(ip)

          d(ip)=b(ip)

          z(ip)=0.d0

23      continue

24    continue

      pause 'too many iterations in jacobi'

      return

      END subroutine

*****************************************************************************

      END MODULE



















