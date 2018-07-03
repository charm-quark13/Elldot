      Module test

      IMPLICIT NONE

      integer,parameter :: spin=2, sites = 2!, nmax = 100
      integer, parameter :: dim=spin*sites!, ncom=dim*2!,lwork=136
      real(8), parameter :: t=0.5d0

      real(8) :: ntarget(dim*2),dd(8),hmat(dim,dim)
      real(8) :: density,dphi(dim,dim),En(dim)!,omega(n-1)
      character (30) :: matrix, vector

      contains

C  Beginning of the function being minimized through the CG method.
C  In this case, func = \int (n{v(x)} - ntarget{v(x)})**2 dx
!      function func(diffv) result(integral)
      function func(v)
      implicit none
      integer :: i
      real(8) :: func, x, v(8)!, diffv(n-1)

      !write(*,vector) v
      !write(*,*) '^^^^^^ pot ^^^^^^'
      !write(*,vector) dens
      !write(*,*) '^^^^^^^^ dens ^^^^^^^^'
      x = 0.d0
      do i=1,8
        x = x + (v(i)-dd(i))**2
      end do

      func = x

      end function
*************************************************************************
C  Derivative subroutine of func.
C  Derivative is of the form:
C  dF(n{v(x)}) = 2 * \int (n{v(x)} - ntarget{v(x)})*(dn/dv) dx
      !subroutine dfunc(diffv,diffn)
      subroutine dfunc(p,df)
      implicit none
      integer :: i
      real(8) :: p(8),df(8)

      df = 0.d0
      do i=1,8
        df(i) = 2.d0*(p(i)-dd(i))
      end do

      write(*,vector) df
      write(*,*) '^^^^^^^^^^^ dSdV ^^^^^^^^^^^^^'
      end subroutine

****************************************************************************
****  Numerical Recipes subroutines and functions                       ****
****************************************************************************
        SUBROUTINE frprmn(p,n,ftol,iter,fret)
        INTEGER iter,n,NMAX,ITMAX
        REAL (8) :: fret,ftol,p(n),EPS
        PARAMETER (NMAX=8,ITMAX=200,EPS=1.d-10)
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
      PARAMETER (NMAX=8)
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
      PARAMETER (NMAX=8,TOL=1.d-4)
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
      call mnbrak(ax,xx,bx,fa,fx,fb)
      !write(*,*) ax, xx, bx, 'ax, xx, bx'
      fret=brent(ax,xx,bx,TOL,xmin)
      !write(*,*) xmin, 'xmin'
      do 12 j=1,n
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
12    continue

      do j=1,n
        !write(*,*) xi(j)
      end do
      !write(*,*) '^^^^^^^^^^^^ xi ^^^^^^^^^^^'
      return
      END subroutine
C  (C) Copr. 1986-92 Numerical Recipes Software 6?6>)AY.
****************************************************************************

        SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc)
        REAL (8) :: ax,bx,cx,fa,fb,fc,GOLD,GLIMIT,TINY
        PARAMETER (GOLD=1.618034d0, GLIMIT=100.d0, TINY=1.d-20)
        REAL (8) :: dum,fu,q,r,u,ulim
        fa=f1dim(ax)
        fb=f1dim(bx)
        if(fb.gt.fa)then
          dum=ax
          ax=bx
          bx=dum
          dum=fb
          fb=fa
          fa=dum
        endif
        cx=bx+GOLD*(bx-ax)
        fc=f1dim(cx)
1     if(fb.ge.fc)then
          r=(bx-ax)*(fb-fc)
          q=(bx-cx)*(fb-fa)
          u=bx-((bx-cx)*q-(bx-ax)*r)/(2.d0*sign(max(abs(q-r),TINY),q-r))
          ulim=bx+GLIMIT*(cx-bx)
          if((bx-u)*(u-cx).gt.0.d0)then
            fu=f1dim(u)
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
            fu=f1dim(u)
          else if((cx-u)*(u-ulim).gt.0.d0)then
            fu=f1dim(u)
            if(fu.lt.fc)then
              bx=cx
              cx=u
              u=cx+GOLD*(cx-bx)
              fb=fc
              fc=fu
              fu=f1dim(u)
            endif
          else if((u-ulim)*(ulim-cx).ge.0.d0)then
            u=ulim
            fu=f1dim(u)
          else
            u=cx+GOLD*(cx-bx)
            fu=f1dim(u)
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

        FUNCTION brent(ax,bx,cx,tol,xmin)
        INTEGER ITMAX
        REAL (8) :: brent,ax,bx,cx,tol,xmin,CGOLD,ZEPS
        PARAMETER (ITMAX=1000,CGOLD=.3819660d0,ZEPS=1.0d-10)
        INTEGER iter
        REAL (8) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
        a=min(ax,cx)
        b=max(ax,cx)
        v=bx
        w=v
        x=v
        e=0.
        fx=f1dim(x)
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
          fu=f1dim(u)
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
