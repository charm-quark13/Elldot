      Module CGMethod

      IMPLICIT NONE


      integer,parameter :: spin=2, sites = 2
      integer, parameter :: dim=spin*sites!,lwork=136
      real(kind=8), parameter :: t=0.5d0
      !real(kind=8), parameter :: dab=2.d0,dbc=-3.d0,dcd=1.5d0

      !complex(8) :: hmat(dim,dim)
      real(kind=8) :: ntarget(dim*2),d0(dim),hmat(dim,dim)
      real(kind=8) :: density,dphi(dim,dim),En(dim)!,omega(n-1)

      character (len=30) :: matrix

      contains

C  Beginning of the function being minimized through the CG method.
C  In this case, func = \int (n{v(x)} - ntarget{v(x)})**2 dx
!      function func(diffv) result(integral)
      function func(v) result(integral)
      implicit none
      integer :: i,j
      real(kind=8) :: dens(dim), integral, v(dim)!, diffv(n-1)

!      write(*,matrix) v

C  Func is fed the potential difference between the lattice points, so
C  v must be calculated to recalculate the eigenvectors of the system.
!      do i=2,n
        v(1) = 0.d0
!        v(i) = v(i-1) - diffv(i-1)
!      end do

C  Recalculating the eigenvectors through the Schrodiner Eqn.
      call hbuild(v,hmat)

!      write(*,matrix) v

      dens = 0.d0

C  Recalculating the density with the new orbitals/eigenvectors.
      do i=1,dim
        do j=1,2
        dens(i) = dens(i) + 2.d0*hmat(i,j)**2
        end do
      end do

C  Calculating the difference integral between the target and current densities.
      integral = 0.d0
      do i =1,dim
        integral = integral + (dens(i)-ntarget(i))**2
      end do

      end function
*************************************************************************
C  Derivative subroutine of func.
C  Derivative is of the form:
C  dF(n{v(x)}) = 2 * \int (n{v(x)} - ntarget{v(x)})*(dn/dv) dx
      !subroutine dfunc(diffv,diffn)
      subroutine dfunc(v,dSdV)

      implicit none
      integer :: i,j,alpha,beta
      real(8) :: v(dim),x,dens(dim),dSdV(dim)

      dens = 0.d0
      do i=1,dim
        do j=1,2
        dens(i) = dens(i) + 2.d0*hmat(i,j)**2
        end do
      end do

      do j=1,dim
        x = 0.d0
        do i=1,dim
          do alpha = 1,2
            do beta = 1,dim
!              write(*,*) derivative(v,hmat,En,alpa,bta,i,j)
              x = x + 4.d0*(dens(i)-ntarget(i))
     &                *derivative(v,alpha,beta,i,j)
            end do
          end do
        end do

        dSdV(j) = x
      end do

      write(*,matrix) dSdV

      end subroutine

***************************************************************************
      function derivative(v,alpha,beta,i,j) result(dn)
      implicit none

      integer :: i,j,alpha,beta
      real(8) :: num, denom, dn, v(dim)

      if (alpha.eq.beta) then
        dn = 0.d0
      else
        num = hmat(i,alpha)*hmat(i,beta)*hmat(j,alpha)*hmat(j,beta)
        denom = En(beta) - En(alpha)
        dn = num/denom
      end if

      end function

***************************************************************************
        subroutine hbuild(v,mat)
        implicit none

        integer, parameter :: lwork=136
        integer :: i,j,info
        real (kind=8) :: v(dim), x, mat(dim,dim)
        real(kind=8) :: work(lwork)

C  Setting up noninteracting Hubbard-like system Hamiltonian with t as the
C  hopping constant between lattice sites.
        do i = 1,dim
          do j=1,dim
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
        call DSYEV('v','u',dim,mat,dim,En,work,lwork,info)

C  Normalization step to normalize each orbital to 1.
        do j = 1,dim
          x = 0.d0
            do i = 1, dim
              x = x + mat(i,j)**2
            end do
            do i = 1, dim
              mat(i,j) = mat(i,j)/dsqrt(x)
            end do
        end do

        end subroutine hbuild

****************************************************************************
****  Numerical Recipes subroutines and functions                       ****
****************************************************************************
        SUBROUTINE frprmn(p,n,ftol,iter,fret)
        INTEGER iter,n,NMAX,ITMAX
        REAL (8) :: fret,ftol,p(n),EPS
        PARAMETER (NMAX=50,ITMAX=200,EPS=1.d-10)
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
          if(2.d0*abs(fret-fp).le.ftol*(abs(fret)+abs(fp)+EPS))return
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
        END
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
          u=bx-((bx-cx)*q-(bx-ax)*r)/(2.d0*sign(max(abs(q-r),TINY),q-r))
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
2       if(abs(d).ge.tol1) then
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
      End Module
