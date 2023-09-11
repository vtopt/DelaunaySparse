c  This is a minimal version of the bqpd code by R. Fletcher as needed by
c  the PROJECT subroutine from DELAUNAYSPARSE.
c
c  This subset of code has been "flattened" into a single file.
c
c  For the original code, contact Sven Leyffer at Argonne National Lab.


C  --- Begin file: auxil.f

cut here >>>>>>>>>>>>>>>>>

c*********** auxiliary routines requiring additional storage ************

c  Copyright, University of Dundee (R.Fletcher), June 1996
c  Current version dated 18/02/99

      subroutine major_s(p,n,k,kmax,a,la,aa,ll,an,bn,rh,ak,bk,ck,lv,
     *  ws,lws,jmin,jmax,w,ls,sgs,pp,smx,ep)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),aa(*),ll(*),an(*),bn(*),rh(*),
     *  ak(*),bk(*),ck(*),lv(*),ws(*),lws(*),w(*),ls(*)
      common/noutc/nout
      common/iprintc/iprint
c     write(nout,*)'major_s  p =',p
      call tfbsub(n,a,la,p,an,an,aa,ll,ep,.true.)
c  compute  s = Y.e_p - Z.M^(-1).Zt.G.Y.e_p  in an
      do i=1,n
        w(i)=an(i)
      enddo
c     write(nout,*)'s =',(w(i),i=1,n)
      if(kmax.gt.0)then
        call linf(n,an,smx,pp)
        call gdotx(n,an,ws,lws,bn)
c       write(nout,*)'G.Y.e_p =',(bn(i),i=1,n)
        sgs=scpr(0.D0,an,bn,n)
        if(k.gt.0)then
c         write(nout,*)'Y.e_p =',(an(i),i=1,n)
          call ztaq(0,n,k,a,la,bn,bk,w,lv,aa,ll)
c  compute  R^(-T).Zt.G.Y.e_p  in ak (ak and bk also used in extend)
c         write(nout,*)'Z#t.G.Y.e_p =',(bk(i),i=1,k),sgs
          do i=1,k
            ak(i)=bk(i)
          enddo
          call rtsol(k,kk,kmax,rh,ak)
c         write(nout,*)'R^(-T).Zt.G.Y.e_p =',(ak(i),i=1,k+1)
          do i=1,k
            ck(i)=ak(i)
          enddo
          call rsol(k,kk,kmax,rh,ck)
c         write(nout,*)'ck =',(ck(i),i=1,k)
          call zprod(p,n,k,a,la,w,ck,lv,w,ls,aa,ll)
          bk(k+1)=sgs
          sgs=-scpr(-sgs,ak,ak,k)
        endif
c       write(nout,*)'sgs =',sgs
      else
        sgs=0.D0
      endif
c     write(nout,*)'s =',(w(i),i=1,n)
c  form At.s
      do j=jmin,jmax
        i=abs(ls(j))
        if(i.gt.n)then
          w(i)=aiscpr(n,a,la,i-n,w,0.D0)
        endif
      enddo
      return
      end

      subroutine minor_s(n,k,kmax,a,la,aa,ll,an,rh,ak,rg,rs,lv,
     *  jmin,jmax,w,ls,sg,sgs)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),aa(*),ll(*),an(*),rh(*),ak(*),rg(*),rs(*),
     *  lv(*),w(*),ls(*)
      common/noutc/nout
      common/minorc/c
c     write(nout,*)'minor_s'
      if(sgs.gt.0.D0)then
        do i=1,k
          rs(i)=rg(i)
        enddo
        call rtsol(k,kk,kmax,rh,rs)
        sgs=scpr(0.D0,rs,rs,k)
        sg=sgs
        if(sgs.eq.0.D0)return
        call rsol(k,kk,kmax,rh,rs)
      else
        do i=1,k
          rs(i)=ak(i)
        enddo
        call rtsol(k,kk,kmax,rh,rs)
        sgs=scpr(0.D0,rs,rs,k)
        c=min(1.D0-sgs,0.D0)
        sgs=sgs*c
        call rsol(k,kk,kmax,rh,rs)
        sg=scpr(0.D0,rs,rg,k)
        if(sg.lt.0.D0)then
          do i=1,k
            rs(i)=-rs(i)
          enddo
          sg=-sg
          c=-c
        endif
      endif
c     write(nout,*)'rs =',(rs(i),i=1,k)
c     write(nout,*)'sg,sgs,c',sg,sgs,c
      call zprod(0,n,k,a,la,an,rs,lv,w,ls,aa,ll)
c     write(nout,*)'s =',(an(i),i=1,n)
c  form At.s
      do j=jmin,jmax
        i=abs(ls(j))
        if(i.gt.n)then
          w(i)=aiscpr(n,a,la,i-n,an,0.D0)
        else
          w(i)=an(i)
        endif
      enddo
      return
      end

      subroutine extend(k,kmax,p,rh,ak,bk,lv,sgs)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension rh(*),ak(*),bk(*),lv(*)
      common/noutc/nout
      common/epsc/eps,tol,emin
c  extend factors of M
      kp=k+1
      lv(kp)=p
      ij=kp
      if(sgs.gt.0.D0)then
c       write(nout,*)'ak =',(ak(i),i=1,k),'   sgs =',sgs
        do i=1,k
          rh(ij)=ak(i)
          ij=ij+kmax-i
        enddo
        rh(ij)=sqrt(sgs)
      else
        u=bk(kp)
c       if(u.gt.1.D0)then
c         b=sqrt(u)
c         c=1.D0/b
c       else
c         b=1.D0
c         c=1.D0/(1.D0+sqrt(1.D0-u))
c       endif
c       a=(c*u-1.D0/c)*5.D-1
        bsq=max(tol,scpr(0.D0,bk,bk,k)/rh(1)**2+abs(u))
        b=sqrt(bsq)
        c=1.D0/(b+sqrt(bsq-u))
        a=c*u-b
        do i=1,k
          bk(i)=c*bk(i)
          ak(i)=bk(i)
          rh(ij)=0.D0
          ij=ij+kmax-i
        enddo
        ak(kp)=a
        bk(kp)=b
        ij=1
        do i=1,k
          call angle(rh(ij),bk(i),cos,sin)
          ij=ij+1
          call rot(kp-i,rh(ij),bk(i+1),cos,sin)
          ij=ij+kmax-i
        enddo
        rh(ij)=bk(kp)
      endif
      k=kp
c     write(nout,*)'extend: V-list =',(lv(i),i=1,k)
c     write(nout,*)'extended reduced Hessian factor'
c     ij=0
c     do i=1,k
c       write(nout,*)(rh(ij+j),j=1,k-i+1)
c       ij=ij+kmax-i+1
c     enddo
c     if(sgs.le.0.D0)write(nout,*)'negative part:  a =',(ak(i),i=1,k)
      return
      end

      subroutine rgup(k,ak,rg,alpha,sgs)
      implicit double precision (a-h,o-z)
      dimension ak(*),rg(*)
      common/noutc/nout
      common/minorc/c
      if(sgs.gt.0.D0)then
        c=1.D0-alpha
        do i=1,k
          rg(i)=rg(i)*c
        enddo
      else
        call mysaxpy(-alpha*c,ak,rg,k)
      endif
c     write(nout,*)'reduced gradient =',(rg(i),i=1,k)
      return
      end

      subroutine zprod(p,n,k,a,la,an,ak,lv,w,ls,aa,ll)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),an(*),ak(*),lv(*),w(*),ls(*),aa(*),ll(*)
      common/noutc/nout
      do i=1,n-k
        w(abs(ls(i)))=0.D0
      enddo
      if(p.gt.0)w(p)=1.D0
      do i=1,k
        w(lv(i))=-ak(i)
      enddo
      call tfbsub(n,a,la,0,w,an,aa,ll,ep,.false.)
      return
      end

      subroutine qinv(q,qv,k,ak,lv)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension ak(*),lv(*)
c  look for q in V-list
      do i=1,k
        if(q.eq.lv(i))then
          do j=1,k
            ak(j)=0.D0
          enddo
          qv=i
          ak(qv)=1.D0
          return
        endif
      enddo
      qv=0
      return
      end

      subroutine ztaq(q,n,k,a,la,an,ak,w,lv,aa,ll)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),an(*),ak(*),w(*),lv(*),aa(*),ll(*)
      common/noutc/nout
      call fbsub(n,1,k,a,la,q,an,w,lv,aa,ll,.false.)
      do i=1,k
        ak(i)=w(lv(i))
      enddo
c     write(nout,*)'Zt.aq =',(ak(i),i=1,k)
      return
      end

      subroutine reduce(k,kmax,p,pv,rh,rg,lv,ak,bk,ck,sgs,rg_up)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension rh(*),rg(*),lv(*),ak(*),bk(*),ck(*)
      common/noutc/nout
      logical rg_up
      external r_shift
c     write(nout,*)'reduce:  p =',p,'   pv =',pv
c  remove column pv from V and update reduced Hessian factors
c     write(nout,*)'ak =',(ak(i),i=1,k)
c     write(nout,*)'bk =',(bk(i),i=1,k)
      ij=pv
      do i=1,pv
        ck(i)=rh(ij)
        call r_shift(rh(ij),k-pv,1)
        ij=ij+kmax-i
      enddo
      kk=ij-kmax+pv
      kmx1=kmax+1
      ij=ij+1
      za=ak(pv)
      zb=bk(pv)
      zr=rg(pv)
      do i=pv+1,k
        ck(i)=rh(ij)
        call r_shift(rh(ij),k-i,1)
        im=i-1
        ak(im)=ak(i)
        bk(im)=bk(i)
        rg(im)=rg(i)
        lv(im)=lv(i)
        ij=ij+kmx1-i
      enddo
c     write(nout,*)'ck =',(ck(i),i=1,k)
      k=k-1
      if(k.eq.0)return
c     write(nout,*)'V list =',(lv(i),i=1,k)
c  return to triangular form
      if(p.ne.0)then
        call brots(k,kmax,pv,kk,rh,ck)
        call mysaxpy(-ck(1)/zb,bk,rh,k)
        call frots(k,k,kmax,rh,ck)
        if(sgs.le.0.D0)call mysaxpy(-za/zb,bk,ak,k)
      else
        call frots(k-pv+1,k-pv+1,kmax-pv+1,rh(kk),ck(pv))
      endif
      if(rg_up)call mysaxpy(-zr/zb,bk,rg,k)
      if(sgs.le.0.D0)then
c  check if negative part can be removed
        do i=1,k
          bk(i)=ak(i)
        enddo
        call rtsol(k,kk,kmax,rh,bk)
        sgs=-scpr(-1.D0,bk,bk,k)
        if(sgs.gt.0.D0)then
          call brots(k,kmax,k,kk,rh,bk)
          sgs=sqrt(sgs)
          do i=1,k
            rh(i)=rh(i)*sgs
          enddo
          call frots(k,k-1,kmax,rh,bk)
        endif
      endif
c     write(nout,*)'reduced Hessian factor'
c     ij=0
c     do i=1,k
c       write(nout,*)(rh(ij+j),j=1,k-i+1)
c       ij=ij+kmx1-i
c     enddo
c     if(sgs.le.0.D0)write(nout,*)'negative part:  a =',(ak(i),i=1,k)
c     write(nout,*)'reduced gradient',(rg(i),i=1,k)
      return
      end

      subroutine revise(k,kmax,rh,rg,lv,q,ak,bk,ck,sgs)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension rh(*),rg(*),lv(*),ak(*),bk(*),ck(*)
      common/noutc/nout
c     write(nout,*)'revise'
c     write(nout,*)'bk =',(bk(i),i=1,k)
      smx=bk(k)
      lv(k)=q
      km=k-1
      ij=k
      do i=1,k
        ck(i)=rh(ij)/smx
        rh(ij)=0.D0
        ij=ij+kmax-i
      enddo
      ij=ij-kmax+k
c     write(nout,*)'ak =',(ak(i),i=1,k)
c     write(nout,*)'ck =',(ck(i),i=1,k)
      call brots(km,kmax,k,ij,rh,ck)
      call mysaxpy(-ck(1),bk,rh,km)
      rh(k)=ck(1)
      call frots(k,km,kmax,rh,ck)
      if(sgs.le.0.D0)then
        ak(k)=ak(k)/smx
        call mysaxpy(-ak(k),bk,ak,km)
      endif
      rg(k)=rg(k)/smx
      call mysaxpy(-rg(k),bk,rg,km)
c     write(nout,*)'V list =',(lv(i),i=1,k)
c     write(nout,*)'revised reduced Hessian factor'
c     ij=0
c     do i=1,k
c       write(nout,*)(rh(ij+j),j=1,k-i+1)
c       ij=ij+kmax-i+1
c     enddo
c     if(sgs.le.0.D0)write(nout,*)'negative part:  a =',(ak(i),i=1,k)
c     write(nout,*)'reduced gradient',(rg(i),i=1,k)
      return
      end

      subroutine hot_start(n,k,kmax,a,la,x,an,bn,rh,ak,lv,w,ls,aa,ll,
     *  ws,lws,hstep)
      implicit double precision (a-h,o-z)
      dimension a(*),la(*),x(*),an(*),bn(*),rh(*),ak(*),lv(*),
     *  w(*),ls(*),aa(*),ll(*),ws(*),lws(*)
      common/noutc/nout
      common/iprintc/iprint
c     dimension Z(100,20)
c     write(nout,*)'Zt matrix'
c     do i=1,k
c       call tfbsub(n,a,la,lv(i),w,Z(1,i),aa,ll,ep,.false.)
c       write(nout,*)lv(i),':',(Z(j,i),j=1,n)
c     enddo
      call gdotx(n,x,ws,lws,an)
      call saipy(1.D0,a,la,0,an,n)
c     write(nout,*)'g =',(an(i),i=1,n)
      call ztaq(0,n,k,a,la,an,ak,w,lv,aa,ll)
c     write(nout,*)'Zt.g =',(ak(i),i=1,k)
      call rtsol(k,kk,kmax,rh,ak)
      call rsol(k,kk,kmax,rh,ak)
c     write(nout,*)'M^(-1).Zt.g =',(ak(i),i=1,k)
      call zprod(0,n,k,a,la,an,ak,lv,w,ls,aa,ll)
c     write(nout,*)'Z.M^(-1).Zt.g =',(an(i),i=1,n)
      do j=1,n-k
        i=abs(ls(j))
        if(i.le.n)an(i)=0.D0
      enddo
      hstep=sqrt(scpr(0.D0,an,an,n))
      if(iprint.ge.1)write(nout,*)'norm of horizontal step =',hstep
      do i=1,n
        bn(i)=x(i)
        x(i)=x(i)+an(i)
      enddo
c     write(nout,*)'x =',(x(i),i=1,n)
      return
      end

      subroutine checkrh(n,a,la,k,kmax,rh,rg,ak,lv,aa,ll,x,ws,lws,
     *  sgs,tol)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),rh(*),rg(*),ak(*),lv(*),aa(*),ll(*),
     *  x(*),ws(*),lws(*)
      dimension Z(500,20),E(20,20),v(20),w(500)
      common/noutc/nout
      if(n.gt.500.or.k.gt.20)then
        write(6,*)'extend local storage for checkrh'
        stop
      endif
c     write(nout,*)'checkrh'
c     write(nout,*)'Zt matrix'
      do i=1,k
        call tfbsub(n,a,la,lv(i),w,Z(1,i),aa,ll,ep,.false.)
c       write(nout,*)lv(i),':',(Z(j,i),j=1,n)
      enddo
      call gdotx(n,x,ws,lws,w)
      call saipy(1.D0,a,la,0,w,n)
c     write(nout,*)'g =',(w(i),i=1,n)
      do i=1,k
        v(i)=scpr(0.D0,w,Z(1,i),n)
      enddo
c     write(nout,*)'reduced gradient vector Zt.g',(v(i),i=1,k)
      emax=0.D0
      do i=1,k
        emax=max(emax,abs(rg(i)-v(i)))
      enddo
      if(emax.gt.tol)
     *  write(nout,*)'max error in reduced gradient =',emax
c     write(nout,*)'Zt.G.Z matrix'
      do i=1,k
        call gdotx(n,Z(1,i),ws,lws,w)
        do j=1,k
          E(i,j)=scpr(0.D0,w,Z(1,j),n)
        enddo
c       write(nout,*)(E(i,j),j=1,k)
      enddo
      ij=1
      do i=1,k
        do j=i,k
          w(j)=rh(ij+j-i)
        enddo
        do j=i,k
          do jj=i,k
            E(j,jj)=E(j,jj)-w(j)*w(jj)
          enddo
        enddo
        ij=ij+kmax-i+1
      enddo
      if(sgs.lt.0.D0)then
        do j=1,k
          do jj=1,k
            E(j,jj)=E(j,jj)+ak(j)*ak(jj)
          enddo
        enddo
      endif
c     write(nout,*)'error  Zt.G.Z - Rt.R'
      emax=0.D0
      do i=1,k
        do j=1,k
          emax=max(emax,abs(E(i,j)))
        enddo
c       write(nout,*)(E(i,j),j=1,k)
      enddo
      if(emax.gt.tol)write(nout,*)'max error in Rt.R =',emax
      return
      end

C  --- Begin file: bqpd.f

cut here >>>>>>>>>>>>>>>>

c  Copyright, University of Dundee (R.Fletcher), June 1996
c  Current version dated  5 October 2011

      subroutine bqpd(n,m,k,kmax,a,la,x,bl,bu,f,fmin,g,r,w,e,ls,
     *  alp,lp,mlp,peq,ws,lws,m0de,ifail,info,iprint,nout)
      implicit double precision (a-h,r-z), integer (i-q)

c  This routine finds a KT point for the bounded QP problem

c       minimize    f(x) = ct.x + xt.G.x/2

c       subject to  l <= [I : A]t.x <= u                  (t = transpose)

c  where x and c are n-vectors, G is a symmetric n*n matrix, and A is an
c  n*m matrix. If G is also positive semi-definite then the KT point is a
c  global solution, else usually a local solution. The method may also be
c  used efficiently to solve an LP problem (G=0). A recursive form of an
c  active set method is used, using Wolfe's method to resolve degeneracy.
c  Matrix information is made available and processed by calls to external
c  subroutines. Details of these are given in an auxiliary file named either
c  'denseL.f' or 'sparseL.f'.

c  parameter list  (variables in a line starting with C must be set on entry)
c  **************

C  n     number of variables
C  m     number of general constraints (columns of A)
c  k     dimension of the 'reduced space' obtained by eliminating the active
c        constraints (only to be set if mode>=2). The number of constraints in
c        the active set is n-k
C  kmax  maximum value of k (set kmax=0 iff the problem is an LP problem)
C  a(*)  storage of reals associated with c and A. This storage may be provided
C        in either dense or sparse format. Refer to either denseA.f or sparseA.f
C        for information on how to set a(*) and la(*).
C  la(*) storage of integers associated with c and A
C  x(n)  contains the vector of variables. Initially an estimate of the solution
C        must be set, replaced by the solution (if it exists) on exit.
C  bl(n+m)  vector of lower bounds for variables and general constraints
C  bu(n+m)  vector of upper bounds (use numbers less than about 1.e30, and
C        where possible supply realistic bounds on the x variables)
c  f     returns the value of f(x) when x is a feasible solution
c        Otherwise f stores the sum of constraint infeasibilities
C  fmin  set a strict lower bound on f(x) (used to identify an unbounded QP)
c  g(n)  returns the gradient vector of f(x) when x is feasible
c  r(n+m) workspace: stores residual vectors at different levels
c        The sign convention is such that residuals are nonnegative at
c        a solution (except for multipliers of equality constraints)
c  w(n+m) workspace: stores denominators for ratio tests
c  e(n+m) stores steepest-edge normalization coefficients: if mode>2 then
c        information in this vector from a previous call should not be changed.
c        (In mode = 3 or 5 these values provide approximate coefficients)
c  ls(n+m) stores indices of the active constraints in locations 1:n-k and of
c        the inactive constraints in locations n-k+1:n+m. The simple bounds
c        on the variables are indexed by 1:n and the general constraints by
c        n+1:n+m. The sign of ls(j) indicates whether the lower bound (+) or
c        the upper bound (-) of constraint ls(j) is currently significant.
c        Within the set of active constraints, locations 1:peq store the
c        indices of active equality constraints, and indices peq+1:lp(1)-1
c        are indices of any pseudo-bounds (active constraints which are not
c        at their bound). If mode>=2, the first n-k elements of ls must be set
c        on entry
c  alp(mlp) workspace associated with recursion
c  lp(mlp)  list of pointers to recursion information in ls
C  mlp   maximum number of levels of recursion allowed (typically mlp=20 would
c        usually be adequate but mlp=m is an upper bound
c  peq   pointer to the end of equality constraint indices in ls
c  ws(*) real workspace for gdotx (see below), bqpd and denseL.f (or sparseL.f).
c          bqpd requires kmax*(kmax+9)/2+2*n+m locations. Refer to denseL.f
c          (or sparseL.f) for the number of real locations required by these
c          routines. Set the total number in mxws (see "Common" below).
c  lws(*) integer workspace for gdotx, bqpd and denseL.f (or sparseL.f). bqpd
c          requires kmax locations. Refer to denseL.f (or sparseL.f) for the
c          number of integer locations required by these routines and set the
c          total number in mxlws (see "Common" below).
c        The storage maps for ws and lws are set by the routine stmap below
C  m0de  mode of operation (larger numbers imply extra information):
C          0 = cold start (no other information available, takes simple
C                bounds for the initial active set)
C          1 = as 0 but includes all equality constraints in initial active set
C          2 = user sets n-k active constraint indices in ls(j), j=1,..,n-k.
c                For a general constraint the sign of ls(j) indicates which
c                bound to use. For a simple bound the current value of x is used
C          3 = takes active set and other information from a previous call.
C                Steepest edge weights are approximated using previous values.
C          4 = as 3 but it is also assumed that A (but not G) is unchanged so
C                that factors of the basis matrix stored in ws and lws are valid
C                (changes in the vectors c, l, u and the matrix G are allowed)
C          5 = as 3 but it is assumed that the reduced Hessian matrix Zt.G.Z is
C                unchanged - use with care as this mode allows A (and hence Z)
C                to change
C          6 = as for 4 but it is also assumed that the reduced Hessian matrix
C                is unchanged, e.g. when both A and G are unchanged
C        A local copy (mode) of m0de is made and may be changed by bqpd
c  ifail   outcome of the process
c              0 = solution obtained
c              1 = unbounded problem detected (f(x)<=fmin would occur)
c              2 = bl(i) > bu(i) for some i
c              3 = infeasible problem detected in Phase 1
c              4 = incorrect setting of m, n, kmax, mlp, mode or tol
c              5 = not enough space in lp
c              6 = not enough space for reduced Hessian matrix (increase kmax)
c              7 = not enough space for sparse factors (sparse code only)
c              8 = maximum number of unsuccessful restarts taken
c             >8 = possible use by later sparse matrix codes
c  info  a vector giving information about the progress of the routines for
c        manipulating the basis matrix (see denseL.f or sparseL.f for details)
c        (info(1) always counts the number of pivots taken)
C  iprint  switch for diagnostic printing (0 = off, 1 = summary,
C                 2 = scalar information, 3 = verbose)
C  nout  channel number for output

c  Common
c  ******
c  User information about the lengths of ws and lws is supplied to bqpd in
c    common/wsc/kk,ll,kkk,lll,mxws,mxlws
c  mxws and mxlws must be set to the total length of ws and lws.  kk and ll
c  refer to the length of ws and lws that is used by gdotx. kkk and lll
c  are the numbers of locations used by bqpd and are set by bqpd.

c  User subroutine
c  ***************
c  The user must provide a subroutine to calculate the vector v := G.x from a
c  given vector x. The header of the routine is
c            subroutine gdotx(n,x,ws,lws,v)
c            dimension x(*),ws(*),lws(*),v(*)
c  In the case that G=0 (i.e. kmax=0) the subroutine is not called by bqpd and
c  a dummy routine may be supplied. Otherwise the user may use the parameters
c  ws and lws (see above) for passing real or integer arrays relating to G.
c  Locations ws(1),...,ws(kk) are available to the user for storing any real
c  information to be used by gdotx. Likewise locations lws(1),...,lws(ll) are
c  available for storing any integer information. Default values are kk=ll=0.
c  Any other setting is made by changing  common/wsc/kk,ll,kkk,lll,mxws,mxlws

c  Tolerances and accuracy
c  ***********************
c  bqpd uses tolerance and accuracy information stored in
c     common/epsc/eps,tol,emin
c     common/repc/sgnf,nrep,npiv,nres
c     common/refactorc/nup,nfreq
c  eps must be set to the machine precision (unit round-off) and tol is a
c  tolerance which gives the hoped for relative accuracy in the solution.
c  The tolerance strategy in the code assumes that the problem is well-scaled
c  and a pre-scaling routine cscale is supplied in denseA.f or sparseA.f.
c  The parameter sgnf is used to measure the maximum allowable relative error
c  in two numbers that would be equal in exact arithmetic.
c     Default values are set in block data but can be reset by the user.
c  Suitable values for both single and double precision calculation are given
c  on line 1227. IT IS IMPORTANT THAT THESE MATCH THE PRECISION BEING USED.
c     The code allows one or more refinement steps after the
c  calculation has terminated, to improve the accuracy of the solution,
c  and a fixed number nrep of such repeats is allowed. However the code
c  terminates without further repeats if no more than npiv pivots are taken.
c    When poor agreement is detected in numbers that would be equal in exact
c  arithmetic, and also in other cases of breakdown, the code is restarted,
c  usually in mode 2. The maximum number of unsuccessful restarts allowed
c  is set in nres.
c    The basis matrix may be refactorised on occasions, for example to prevent
c  build-up of round-off in the factors or (when using sparse.f) to improve
c  the sparsity of the factors. The maximum interval between refactorizations
c  is set in nfreq.

      parameter (ncyc=4,ncyc2=8)
      dimension lcyc(ncyc2)
      dimension a(*),x(*),bl(*),bu(*),g(*),r(*),w(*),e(*),alp(*),ws(*),
     *  la(*),ls(*),lp(*),lws(*),info(*)
      character*32 spaces
      common/bqpdc/irh1,na,na1,nb,nb1,ka1,kb1,kc1,irg1,lu1,lv,lv1,ll1
      common/epsc/eps,tol,emin
      common/vstepc/vstep
      common/repc/sgnf,nrep,npiv,nres
      common/wsc/kk,ll,kkk,lll,mxws,mxlws
      common/refactorc/nup,nfreq
      common/alphac/alpha,rp,pj,qqj,qqj1
      logical psb,degen,eqty
      spaces='         '
      mode=m0de
      iter=0
      info(1)=0
      if(m.lt.0.or.n.le.0.or.mlp.lt.2.or.mode.lt.0.or.mode.gt.6
     *  .or.kmax.lt.0)then
        ifail=4
        return
      endif
      small=max(1.D1*tol,sqrt(eps))
      smallish=max(sqrt(small),1.D1*small)
c     smallish=max(eps/tol,1.D1*small)
      tolp1=tol+1.D0
      nm=n+m
      nmi=nm
      if(iprint.ge.3)then
        write(nout,1000)'lower bounds',(bl(i),i=1,nm)
        write(nout,1000)'upper bounds',(bu(i),i=1,nm)
      endif
      irep=0
      ires=0
      mres=0
      bestf=1.D37
      phase=1
      do i=1,nm
        if(bl(i).gt.bu(i))then
          ifail=2
          return
        endif
      enddo
      if(mode.le.2)then
        do i=1,n
          x(i)=min(bu(i),max(bl(i),x(i)))
        enddo
        call stmap(n,nm,kmax)
        if(mode.eq.0)then
          nk=0
        elseif(mode.eq.1)then
c  collect near-equality c/s
          nk=0
          do i=1,nm
            if(bu(i)-bl(i).le.tol)then
              nk=nk+1
              ls(nk)=i
            endif
          enddo
c         write(nout,*)'number of eqty c/s =',nk
        else
          nk=n-k
        endif
      else
        nk=n-k
      endif
    3 continue
      lev=1
      if(mode.le.3.or.mode.eq.5)then
c  set up factors of basis matrix and permutation vectors
        ifail=mode
        call start_up(n,nm,nmi,a,la,nk,e,ls,ws(lu1),lws(ll1),mode,ifail)
        if(ifail.gt.0)return
        if(mode.le.3)nk=n
      endif
c  collect active near-equality c/s
      peq=0
      do j=1,nk
        i=abs(ls(j))
        if(bu(i)-bl(i).le.tol)then
          if(i.le.n)x(i)=min(bu(i),max(bl(i),x(i)))
          peq=peq+1
          call iexch(ls(peq),ls(j))
        endif
      enddo
      k=n-nk
    5 continue
c  refinement step loop
      mpiv=iter+npiv
      if(mode.eq.0)goto8
      call warm_start(n,nm,nk,a,la,x,bl,bu,r,ls,lws(lv1),ws(lu1),
     *  lws(ll1),ws(na1),vstep)
      if(vstep.gt.tol)mpiv=0
      if(mode.eq.4)then
        do i=lv1,lv+k
          q=lws(i)
          do j=nk+1,nm
            if(abs(ls(j)).eq.q)then
              nk=nk+1
              call iexch(ls(nk),ls(j))
              if(bu(q)-bl(q).le.tol)then
                x(q)=min(bu(q),max(bl(q),x(q)))
                call iexch(ls(nk),ls(lp(1)))
                peq=peq+1
                call iexch(ls(peq),ls(lp(1)))
                lp(1)=lp(1)+1
              endif
              goto6
            endif
          enddo
    6     continue
        enddo
        nk=n
        k=0
      endif
      if(k.gt.0)then
        call hot_start(n,k,kmax,a,la,x,ws(na1),ws(nb1),ws(irh1),
     *    ws(ka1),lws(lv1),w,ls,ws(lu1),lws(ll1),ws,lws,hstep)
        if(hstep.gt.tol)mpiv=0
      endif
c     write(nout,*)'x =',(x(i),i=1,n)
    8 continue
c  calculate residuals of inactive c/s
      call residuals(n,nk+1,nm,a,la,x,bl,bu,r,ls,qj,ninf)
      if(ninf.eq.0.or.k.eq.0)goto9
c  try dual asm iterations
      if(ninf.gt.k)then
        do i=1,n
          x(i)=ws(nb+i)
        enddo
        mode=4
        goto5
      endif
      q=abs(ls(qj))
c     print *,'q,ninf,k',q,ninf,k
c     print *,'V list',(lws(i),i=lv1,lv+k)
      call qinv(q,qqv,k,ws(kb1),lws(lv1))
      if(qqv.eq.0)then
        call ztaq(q,n,k,a,la,w,ws(kb1),w,lws(lv1),ws(lu1),lws(ll1))
        call linf(k,ws(kb1),zmx,pv)
        if(zmx.le.tol)then
          mode=4
          goto5
        endif
      endif
      do i=0,k-1
        ws(kc1+i)=ws(kb1+i)
      enddo
      call rtsol(k,kk_,kmax,ws(irh1),ws(kc1))
      sgs=scpr(0.D0,ws(kc1),ws(kc1),k)
      call rsol(k,kk_,kmax,ws(irh1),ws(kc1))
      call zprod(0,n,k,a,la,ws(na1),ws(kc1),lws(lv1),w,ls,ws(lu1),
     *  lws(ll1))
c     write(nout,*)'dual step =',(ws(i),i=na1,na+n)
      if(qqv.eq.0)then
        p=lws(lv+pv)
        if(iprint.ge.2)write(nout,*)'replace',p,' by',q
        call pivot(p,q,n,nmi,a,la,e,ws(lu1),lws(ll1),ifail,info)
        if(ifail.ge.1)then
          if(ifail.ge.2)return
          if(iprint.ge.1)write(nout,*)'dual step fails'
          mode=2
          goto3
        endif
      else
        p=0
        pv=qqv
      endif
c     write(nout,*)'r(q), sgs =',r(q),sgs
      alpha=r(q)/sign(sgs,dble(ls(qj)))
      do i=1,n
        ws(nb+i)=x(i)
        x(i)=x(i)+alpha*ws(na+i)
      enddo
c     write(nout,*)'x =',(x(i),i=1,n)
      call reduce(k,kmax,p,pv,ws(irh1),ws(irg1),lws(lv1),ws(ka1),
     *  ws(kb1),ws(kc1),sgs,.false.)
      nk=nk+1
      call iexch(ls(qj),ls(nk))
      if(bu(q)-bl(q).le.tol)then
        ws(nb+q)=bl(q)
        peq=peq+1
        call iexch(ls(peq),ls(nk))
      endif
c     call  checkrh(n,a,la,k,kmax,ws(irh1),ws(irg1),ws(ka1),lws(lv1),
c    *  ws(lu1),lws(ll1),x,ws,lws,sgs,tol)
      goto8
    9 continue
c  remove redundant equations
      if(mode.eq.1)then
        do j=nm,n+1,-1
          i=abs(ls(j))
          if(r(i).eq.0.D0.and.bu(i)-bl(i).le.tol)then
            call iexch(ls(j),ls(nm))
            nm=nm-1
          endif
        enddo
        if(iprint.ge.1)write(nout,*)nmi-nm,'  redundant equations'
      endif
c  collect pseudo-bounds
      infb=0
      lp(1)=peq+1
      do j=lp(1),nk
        i=abs(ls(j))
        if(i.le.n)then
          if(abs(x(i)-bl(i)).le.tol)then
            x(i)=bl(i)
            ls(j)=i
          elseif(abs(x(i)-bu(i)).le.tol)then
            x(i)=bu(i)
            ls(j)=-i
          else
            call iexch(ls(j),ls(lp(1)))
            if(x(i).lt.bl(i).or.x(i).gt.bu(i))infb=infb+1
            lp(1)=lp(1)+1
          endif
        endif
      enddo
      ninf=ninf+infb
      gtol=tol
c  phase 1 loop
   10 continue
      if(ninf.eq.0)then
        if(iprint.ge.1)write(nout,*)'FEASIBILITY OBTAINED at level 1'
        call setfg2(n,kmax,a,la,x,f,g,ws,lws)
        fold=f
        psb1=lp(1)-peq
        peq1=peq
        fb=0.D0
        qqq=0
        if(iprint.ge.1)write(nout,'(''pivots ='',I5,
     *    ''  level = 1    f ='',E16.8)')info(1),f
      else
        call setfg1(n,nm,a,la,x,bl,bu,f,fb,g,peq,lp(1),nk+1,ls,r)
        if(iprint.ge.1)write(nout,'(''pivots ='',I5,
     *    ''  level = 1    f ='',E16.8,''   ninf ='',I4)')
     *    info(1),f,ninf
      endif
c     write(nout,*)'g =',(g(i),i=1,n)
      gnorm=sqrt(scpr(0.D0,g,g,n))
      if(iprint.ge.2)write(nout,*)'gnorm =',gnorm
      gtol=max(tol*gnorm,gtol)
      call newg
      do i=1,ncyc2
        lcyc(i)=-i
      enddo
      goto16
c  start of major iteration
   15 continue
      if(iprint.ge.1)then
        if(ninf.eq.0)then
          if(k.gt.0)then
            write(nout,'(''pivots ='',I5,
     *        ''  level = 1    f ='',E16.8,''   k ='',I4)')info(1),f,k
          else
            write(nout,'(''pivots ='',I5,
     *        ''  level = 1    f ='',E16.8)')info(1),f
          endif
        else
          write(nout,'(''pivots ='',I5,''  level = 1    f ='',
     *      E16.8,''   ninf ='',I4)')info(1),f,ninf
        endif
      endif
      if(ninf.eq.0.and.kmax.gt.0)then
        call gdotx(n,x,ws,lws,g)
        call saipy(1.D0,a,la,0,g,n)
        gnorm=max(gnorm,sqrt(scpr(0.D0,g,g,n)))
        if(iprint.ge.2)write(nout,*)'gnorm =',gnorm
        if(fold-f.le.max(tol,abs(fold)*smallish).and.
     *    psb1.eq.lp(1)-peq.and.peq1.eq.peq)gtol=2.D0*gtol
        fold=f
        psb1=lp(1)-peq
        peq1=peq
        gtol=max(tol*gnorm,gtol)
        call newg
      endif
c  cycle detection
      do i=1,ncyc
        if(lcyc(i).ne.lcyc(i+ncyc))goto16
      enddo
      gtol=1.D1*gtol
   16 continue
c     if(info(1).ge.117)iprint=3
c     if(info(1).ge.50)stop
c  calculate multipliers
      call fbsub(n,1,nk,a,la,0,g,w,ls,ws(lu1),lws(ll1),.true.)
      call signst(1,nk,.true.,r,w,e,ls,gtol)
c  opposite bound or reset multiplier loop
   20 continue
c  make psb multipliers negative sign
      do j=peq+1,lp(1)-1
        i=abs(ls(j))
        if(r(i).gt.0.D0)then
          r(i)=-r(i)
          ls(j)=-ls(j)
        endif
      enddo
      if(iprint.ge.3)then
        write(nout,1001)'costs vector and indices',
     *    (ls(j),r(abs(ls(j))),j=1,nk)
        write(nout,1000)'steepest edge coefficients',
     *    (e(abs(ls(j))),j=1,nk)
        if(lp(1).gt.1)write(nout,*)'active equality c/s =',peq,
     *    '   pseudo-bounds =',lp(1)-peq-1
      endif
c     call check(n,nk,nmi,kmax,g,a,la,x,bl,bu,r,ls,ws(nb1),f,
c    *  ws,lws,ninf,peq,lp(1),nk+1,1,p,rp)
      if(ninf.eq.0.and.qqq.gt.0)r(qqq)=max(0.D0,r(qqq))
      rp=0.D0
      call optest(peq+1,nk,r,e,ls,rp,pj)
   25 continue
      if(rp.eq.0.D0)then
        if(ninf.gt.0)then
          if(iprint.ge.1)write(nout,'(''pivots ='',I5,
     *      ''  level = 1    f ='',E16.8,''   ninf ='',I4)')
     *      info(1),f,ninf
        else
          if(iprint.ge.1)write(nout,'(''pivots ='',I5,
     *      ''  level = 1    f ='',E16.8)')info(1),f
        endif
        if(iprint.ge.2)write(nout,*)'OPTIMAL at level 1'
        if(iprint.ge.3)then
c         write(nout,1000)'x variables',(x(i),i=1,n)
          write(nout,1001)'residual vector and indices',
     *      (ls(j),r(abs(ls(j))),j=nk+1,nm)
        endif
        irep=irep+1
        if(irep.le.nrep.and.iter.gt.mpiv)then
          if(iprint.ge.1)write(nout,*)'refinement step #',irep
          mode=6
          goto5
        endif
c     call checkout(n,a,la,ws(lu1),lws(ll1),lws(ll1+n),lws(ll1+n+500))       
c     call checkout(n,a,la,lws(ll1),lws(ll1+n),lws(ll1+2*n),
c    *  lws(ll1+3*n+nmi),lws(ll1+4*n+nmi),lws(ll1+5*n+nmi),
c    *  lws(ll1+6*n+nmi),ws(lu1+5*n),mxws-kk-kkk-5*n,ws(lu1),
c    *  ws(lu1+3*n))
        if(iprint.ge.1)write(nout,*)'total number of restarts =',mres
        if(ninf.gt.0)then
          ifail=3
          return
        endif
c  tidy up x
        do j=lp(1),nk
          i=abs(ls(j))
          if(i.le.n)then
            if(ls(j).ge.0)then
              x(i)=bl(i)
            else
              x(i)=bu(i)
            endif
          endif
        enddo
        do i=1,n
          x(i)=max(min(x(i),bu(i)),bl(i))
        enddo
        do j=nk+1,nm
          i=abs(ls(j))
          if(r(i).eq.0.D0.and.i.le.n)then
            if(ls(j).ge.0)then
              x(i)=bl(i)
            else
              x(i)=bu(i)
            endif
          endif
        enddo
        ifail=0
        return
      endif
      degen=.false.
      psb=pj.lt.lp(1)
      if(psb)lp(1)=lp(1)-1
      call iexch(ls(pj),ls(lp(1)))
   30 continue
      pj=lp(1)
      p=abs(ls(pj))
      qqq=0
c  calculate denominators
      call major_s(p,n,k,kmax,a,la,ws(lu1),lws(ll1),ws(na1),ws(nb1),
     *  ws(irh1),ws(ka1),ws(kb1),ws(kc1),lws(lv1),ws,lws,nk+1,nm,
     *  w,ls,sgs,pp,smx,e(p))
      call signs(nk+1,nm,ls(pj).lt.0,w,ls,tol)
      if(degen)then
c  reset r and w from level 2
        do j=lp(2),nm
          i=abs(ls(j))
          w(i)=min(w(i),0.D0)
          r(i)=0.D0
        enddo
c       call check(n,nk,nmi,kmax,g,a,la,x,bl,bu,r,ls,ws(nb1),f,
c    *    ws,lws,ninf,peq,lp(1),nk+1,1,p,rp)
      endif
c  check consistency of r(p)
      rp=scpr(0.D0,g,ws(na1),n)
      if(ls(pj).lt.0)rp=-rp
      if(abs(rp-r(p)).ge.-sgnf*r(p))then
        call refine_s(p,nk,n,a,la,ls,ws(na1),ws(nb1),ws(lu1),lws(ll1))
c       write(nout,*)'norm_ds =',sqrt(scpr(0.D0,ws(nb1),ws(nb1),n))
        rp=scpr(0.D0,g,ws(na1),n)
        if(ls(pj).lt.0)rp=-rp
        if(iprint.ge.1)
     *    write(nout,*)'new and old values of r(p)',rp,r(p),e(p)
        if(abs(rp).le.gtol*e(p))then
c         if(iprint.ge.1)write(nout,*)'refined r(p) is negligible',rp
          if(psb)lp(1)=pj+1
          r(p)=0.D0
          goto20
        endif
c       if(abs(rp-r(p)).ge.-sgnf*r(p))goto98
        if(abs(rp-r(p)).ge.-sgnf*r(p))then
c         if(iprint.ge.1)write(nout,*)'refined r(p) is negligible',rp
          if(psb)lp(1)=pj+1
          r(p)=0.D0
          goto20
        endif
        r(p)=rp
      endif
      rp=r(p)
      if(ninf.eq.0.and.sgs.gt.eps)then
        ff=f-5.D-1*(rp/sgs)*rp
        if(ff.ge.f*tolp1)then
          if(iprint.ge.1)write(nout,*)
     *      'minimizing step fails to reduce f enough: truncate r(p)'
          if(psb)lp(1)=pj+1
          r(p)=0.D0
          goto20
        endif
      endif
   35 continue
      if(iprint.ge.3)then
        write(nout,1000)'x variables',(x(i),i=1,n)
        write(nout,1001)'residual vector and indices',
     *    (ls(j),r(abs(ls(j))),j=nk+1,nm)
        write(nout,1000)'denominators',(w(abs(ls(j))),j=nk+1,nm)
      endif
      sml=smallish
   40 continue
c  line search
      if(psb)then
        if(ls(pj).gt.0)then
          if(x(p).gt.bu(p))then
            amax=1.D37
          elseif(x(p).lt.bl(p))then
            amax=bl(p)-x(p)
          else
            amax=bu(p)-x(p)
          endif
        else
          if(x(p).lt.bl(p))then
            amax=1.D37
          elseif(x(p).gt.bu(p))then
            amax=x(p)-bu(p)
          else
            amax=x(p)-bl(p)
          endif
        endif
      else
        amax=bu(p)-bl(p)
      endif
      qqj=pj
      if(ninf.eq.0)then
        if(sgs*amax.ge.-rp)then
          amax=-rp/sgs
          qqj=0
        endif
      endif
c  level 1 ratio test
      call ratio(nk+1,nm,r,w,bl,bu,ls,amax,alpha,qqj,qqj1)
c  test update of f
      if(ninf.eq.0)then
        ff=f+alpha*(rp+5.D-1*alpha*sgs)
        if(ff.le.fmin)then
          irep=irep+1
          if(irep.le.nrep.and.iter.gt.mpiv)then
            mode=4
            if(sgs.gt.0.D0)mode=6
            if(iprint.ge.1)write(nout,*)
     *        'unbounded solution identified: refinement step #',irep
            goto5
          else
            ifail=1
c  tidy up x
            do i=1,n
              x(i)=max(min(x(i),bu(i)),bl(i))
            enddo
            do j=nk+1,nm
              i=abs(ls(j))
              if(r(i).eq.0.D0.and.i.le.n)then
                if(ls(j).ge.0)then
                  x(i)=bl(i)
                else
                  x(i)=bu(i)
                endif
              endif
            enddo
            return
          endif
        endif
      else
        ff=f+alpha*rp
        if(qqj1.eq.0.and.ff.le.-small)then
          if(iprint.ge.1)write(nout,*)'negative phase 1 objective'
          if(psb)lp(1)=pj+1
          r(p)=0.D0
          goto20
        endif
      endif
      call ishift(lcyc,ncyc2-1,1)
      lcyc(ncyc2)=p
      if(qqj.eq.0)then
        if(iprint.ge.2)write(nout,*)'line minimum found:  alpha =',
     *      alpha,'   p =',p
        goto50
      endif
      qq=abs(ls(qqj))
      if(iprint.ge.2)then
        write(nout,*)'alpha =',alpha,'   p =',p,'   qq =',qq
        write(nout,*)'r(p),r(qq),w(qq) =',r(p),r(qq),w(qq)
      endif
c  test significance of w(qq)
      if(k.eq.0.and.qqj.ne.pj.and.abs(w(qq)).le.sml)then
        call refine_s(p,nk,n,a,la,ls,ws(na1),ws(nb1),ws(lu1),lws(ll1))
c       write(nout,*)'norm_ds =',sqrt(scpr(0.D0,ws(nb1),ws(nb1),n))
        wqq=w(qq)
        sml=small
        call reset_w(nk+1,nm,ls(pj).lt.0,n,a,la,w,ls,ws(na1),sml)
        if(degen)then
          do j=lp(2),nm
            i=abs(ls(j))
            w(i)=min(w(i),0.D0)
          enddo
        endif
        if(iprint.ge.1)
     *    write(nout,*)'new and old values of w(qq)',w(qq),wqq
        if(iprint.ge.3)write(nout,1000)
     *    'refined denominators',(w(abs(ls(j))),j=nk+1,nm)
        if(wqq*w(qq).le.0.D0)then
          w(qq)=0.D0
          goto40
        endif
        if(nup.gt.0.and.abs(wqq-w(qq)).ge.sgnf*abs(w(qq)))then
          call refactor(n,nm,a,la,ws(lu1),lws(ll1),ifail)
          if(iprint.gt.1)write(nout,*)'refactor: ifail =',ifail
          if(ifail.eq.1)goto98
          if(ifail.gt.0)return
        endif
c       if(w(qq).eq.0.D0)goto40
c       if(abs(wqq-w(qq)).ge.sgnf*abs(w(qq)))goto98
      endif
      eqty=bu(qq)-bl(qq).le.tol
      if(ff.ge.f.and.qqj.ne.pj.and..not.(psb.or.eqty))then
c  potential degeneracy block at level 1
        if(k.gt.0)then
          call qinv(qq,qqv,k,ws(kb1),lws(lv1))
          if(qqv.gt.0)goto50
          call ztaq(qq,n,k,a,la,w,ws(kb1),w,lws(lv1),ws(lu1),lws(ll1))
          call linf(k,ws(kb1),zmx,pv)
          if(zmx.gt.tol)goto50
        endif
        if(abs(r(qq)).gt.tol.or.(r(qq).ge.0.D0.and.w(qq).lt.0.D0))then
          if(iprint.gt.1)write(nout,*)'* truncate r(p)'
          if(psb)lp(1)=pj+1
          r(p)=0.D0
          goto20
        endif
        plev=nm+1
        if(k.eq.0)then
          do j=nm,nk+1,-1
            i=abs(ls(j))
            if(r(i).eq.0.D0)then
              plev=plev-1
              call iexch(ls(j),ls(plev))
              if(bu(i)-bl(i).gt.tol)r(i)=1.D0
            endif
          enddo
        else
          if(degen)then
            plev=lp(2)-1
          else
            plev=nm
          endif
          call iexch(ls(qqj),ls(plev))
          i=abs(ls(plev))
          if(bu(i)-bl(i).gt.tol)r(i)=1.D0
        endif
        if(mlp.lt.3)then
          ifail=5
          return
        endif
        lp(2)=plev
        lev=2
        alp(1)=f
        f=0.D0
        if(iprint.ge.2)write(nout,*)'degeneracy: increase level to 2'
        qj=pj
        q=p
        if(iprint.ge.1)write(nout,'(''pivots ='',I5,''     level = 2'',
     *    ''    f ='',E16.8)')info(1),f
        goto86
      endif
      df=f-ff
      if(k.gt.0.and.qqj.ne.pj)goto50
   45 continue
      if(.not.eqty)qqq=qq
      if(qqj.ne.pj)then
        if(iprint.ge.2)write(nout,*)'replace',p,' by',qq
        call pivot(p,qq,n,nmi,a,la,e,ws(lu1),lws(ll1),ifail,info)
        if(ifail.ge.1)then
          if(ifail.ge.2)return
c         call iexch(ls(pj),ls(qqj))
          if(iprint.ge.1)write(nout,*)'near singularity in pivot (1)'
          goto98
        endif
      endif
c  update f, x, and r
      call setrp(psb,p.le.n,alpha,x(p),bl(p),bu(p),g(p),rpu,ls(pj),
     *  infb,fb,tol)
      if(alpha.gt.0.D0)then
        iter=iter+1
        call update(n,nk+1,nm,a,la,g,r,w,x,bl,bu,ls,alpha,n_inf,fff,tol)
        fff=fff+fb
        n_inf=n_inf+infb
        if(rpu.lt.0.D0)then
          fff=fff-rpu
          n_inf=n_inf+1
        endif
        if(ninf.gt.0)then
          if(n_inf.lt.ninf)then
            ff=fff
            gnorm=sqrt(scpr(0.D0,g,g,n))
            if(iprint.ge.2)write(nout,*)'gnorm =',gnorm
            gtol=max(tol*gnorm,gtol)
            call newg
          else
            ff=min(ff,fff)
          endif
          if(ff.gt.0.D0.and.ff.lt.smallish.and.f.ge.smallish)then
            call iexch(ls(pj),ls(qqj))
            mode=6
            goto99
          endif
        endif
      else
        n_inf=ninf
      endif
      f=ff
      if(qqj.eq.pj)then
c  opposite bound comes active
        if(ninf.eq.0)then
          if(kmax.gt.0)goto15
          if(iprint.ge.1)write(nout,'(''pivots ='',I5,
     *      ''  level = 1    f ='',E16.8)')info(1),f
        elseif(f.le.0.D0)then
          ninf=0
          do j=nk+1,nm
            i=abs(ls(j))
            r(i)=max(r(i),0.D0)
          enddo
          goto10
        elseif(n_inf.lt.ninf)then
          ninf=n_inf
          goto15
        else
          if(iprint.ge.1)write(nout,'(''pivots ='',I5,
     *      ''  level = 1    f ='',E16.8,''   ninf ='',I4)')
     *      info(1),f,ninf
        endif
        r(p)=-rp
        goto20
      endif
      r(p)=rpu
      call iexch(ls(pj),ls(qqj))
      if(eqty)then
        peq=peq+1
        call iexch(ls(peq),ls(pj))
        lp(1)=pj+1
      endif
      if(ninf.eq.0)then
        goto15
      elseif(f.le.0.D0)then
        ninf=0
        do j=nk+1,nm
          i=abs(ls(j))
          r(i)=max(r(i),0.D0)
        enddo
        goto10
      elseif(n_inf.lt.ninf)then
        ninf=n_inf
        goto15
      endif
      goto15
c  code for QP-like iterations
   50 continue
      iter=iter+1
c  set reduced gradient
      do i=irg1,irg1+k-1
        ws(i)=0.D0
      enddo
      if(ls(pj).ge.0)then
        ws(irg1+k)=rp+alpha*sgs
      else
        ws(irg1+k)=-rp-alpha*sgs
      endif
c  update f, x, r
      f=ff
      call setrp(psb,p.le.n,alpha,x(p),bl(p),bu(p),g(p),r(p),ls(pj),
     *  infb,fb,tol)
      if(alpha.gt.0.D0)
     *  call update(n,nk+1,nm,a,la,g,r,w,x,bl,bu,ls,alpha,n_inf,ff,tol)
      if(k.ge.kmax)then
        ifail=6
        return
      endif
c  free constraint p
      call iexch(ls(pj),ls(nk))
      nk=nk-1
      call extend(k,kmax,p,ws(irh1),ws(ka1),ws(kb1),lws(lv1),sgs)
      if(p.gt.n.or.smx.gt.4.D0)then
        call ztaq(pp,n,k,a,la,w,ws(kb1),w,lws(lv1),ws(lu1),lws(ll1))
        if(iprint.ge.2)write(nout,*)'for stability, replace',p,' by',pp
        call pivot(p,pp,n,nmi,a,la,e,ws(lu1),lws(ll1),ifail,info)
        if(ifail.ge.1)then
          if(ifail.ge.2)return
          if(iprint.ge.1)write(nout,*)'near singularity in pivot (2)'
          goto98
        endif
        call revise(k,kmax,ws(irh1),ws(irg1),lws(lv1),pp,ws(ka1),
     *    ws(kb1),ws(kc1),sgs)
      endif
c     call  checkrh(n,a,la,k,kmax,ws(irh1),ws(irg1),ws(ka1),lws(lv1),
c    *  ws(lu1),lws(ll1),x,ws,lws,sgs,tol)
      if(qqj.eq.0)goto15
c     write(nout,*)'reduced gradient =',(ws(i),i=irg1,irg1+k-1)
c  minor iteration loop
   60 continue
      call qinv(qq,qqv,k,ws(kb1),lws(lv1))
      if(qqv.eq.0)then
        call ztaq(qq,n,k,a,la,w,ws(kb1),w,lws(lv1),ws(lu1),lws(ll1))
        call linf(k,ws(kb1),zmx,pv)
        if(zmx.eq.0.D0)then
          if(iprint.ge.1)write(nout,*)'no available pivot in V-list'
          goto98
        endif
        p=lws(lv+pv)
        if(iprint.ge.2)write(nout,*)'replace',p,' by',qq
        call pivot(p,qq,n,nmi,a,la,e,ws(lu1),lws(ll1),ifail,info)
        if(ifail.ge.1)then
          if(ifail.ge.2)return
          if(iprint.ge.1)write(nout,*)'near singularity in pivot (3)'
          goto98
        endif
      else
        p=0
        pv=qqv
        zmx=1.D0
      endif
      call reduce(k,kmax,p,pv,ws(irh1),ws(irg1),lws(lv1),ws(ka1),
     *  ws(kb1),ws(kc1),sgs,.true.)
c     write(nout,*)'reduced gradient =',(ws(i),i=irg1,irg1+k-1)
      nk=nk+1
      call iexch(ls(nk),ls(qqj))
      if(eqty)then
        call iexch(ls(lp(1)),ls(nk))
        peq=peq+1
        call iexch(ls(peq),ls(lp(1)))
        lp(1)=lp(1)+1
      endif
      if(k.eq.0)goto15
      if(iprint.ge.2)write(nout,*)'minor iteration loop'
      if(iprint.ge.1)write(nout,'(''pivots ='',I5,''  level = 1    f =''
     *  ,E16.8,''   k ='',I4,''   minor'')')info(1),f,k
c     call  checkrh(n,a,la,k,kmax,ws(irh1),ws(irg1),ws(ka1),lws(lv1),
c    *  ws(lu1),lws(ll1),x,ws,lws,sgs,tol)
      call minor_s(n,k,kmax,a,la,ws(lu1),lws(ll1),ws(na1),ws(irh1),
     *  ws(ka1),ws(irg1),ws(kc1),lws(lv1),nk+1,nm,w,ls,rp,sgs)
      if(rp.eq.0.D0.and.sgs.ge.0.D0)goto15
      call signs(nk+1,nm,.false.,w,ls,tol)
      if(iprint.ge.3)then
        write(nout,1000)'x variables',(x(i),i=1,n)
        write(nout,1001)'residual vector and indices',
     *    (ls(j),r(abs(ls(j))),j=nk+1,nm)
        write(nout,1000)'denominators',(w(abs(ls(j))),j=nk+1,nm)
      endif
      if(sgs.gt.0.D0)then
        amax=1.D0
      else
        t=2.D0*(f-fmin)
        amax=t/(sqrt(rp**2-t*sgs)+rp)
      endif
      alpha=amax
      qqj=0
      call ratio(nk+1,nm,r,w,bl,bu,ls,amax,alpha,qqj,qqj1)
      f=min(f-alpha*(rp-5.D-1*sgs*alpha),f)
      if(qqj.eq.0)then
        if(sgs.le.0.D0.or.f.le.fmin)then
          irep=irep+1
          if(irep.le.nrep.and.iter.gt.mpiv)then
            mode=4
            if(sgs.gt.0.D0)mode=6
            if(iprint.ge.1)write(nout,*)
     *        'unbounded solution identified: refinement step #',irep
            goto5
          else
            ifail=1
            return
          endif
        else
          if(iprint.ge.2)
     *      write(nout,*)'line minimum obtained:  alpha =',alpha
          qqj=0
        endif
      endif
      call update(n,nk+1,nm,a,la,g,r,w,x,bl,bu,ls,alpha,n_inf,ff,tol)
      if(qqj.eq.0)goto15
      call rgup(k,ws(ka1),ws(irg1),alpha,sgs)
      qq=abs(ls(qqj))
      eqty=bu(qq)-bl(qq).le.tol
      if(iprint.ge.2)write(nout,*)'alpha =',alpha,'   amax =',amax,
     *  '   qq =',qq
      goto60
c  recursive code for resolving degeneracy (Wolfe's method)
   80 continue
c     if(info(1).ge.117)iprint=3
c     if(info(1).ge.120)stop
c  calculate multipliers
      call fbsub(n,1,nk,a,la,0,g,w,ls,ws(lu1),lws(ll1),.true.)
      call signst(1,nk,.true.,r,w,e,ls,gtol)
c  opposite bound or reset multiplier loop
   82 continue
      if(iprint.ge.3)then
        write(nout,1001)'costs vector and indices',
     *    (ls(j),r(abs(ls(j))),j=1,nk)
        write(nout,1000)'steepest edge coefficients',
     *    (e(abs(ls(j))),j=1,nk)
        if(lp(1).gt.1)write(nout,*)'active equality c/s =',peq,
     *    '   pseudo-bounds =',lp(1)-peq-1
      endif
c  psb contribution to optimality test
      do j=peq+1,lp(1)-1
        i=abs(ls(j))
        if(r(i).gt.0.D0)then
          r(i)=-r(i)
          ls(j)=-ls(j)
        endif
      enddo
      call optest(peq+1,nk,r,e,ls,rp,pj)
c  84 continue
      if(rp.eq.0.D0.or.pj.lt.lp(1))then
        if(iprint.ge.2)write(nout,*)'return to level 1'
        lev=1
        f=alp(1)
        do j=lp(2),nm
          r(abs(ls(j)))=0.D0
        enddo
        lev=1
        if(rp.eq.0.D0)goto25
        goto20
      endif
      call iexch(ls(pj),ls(lp(1)))
      pj=lp(1)
      p=abs(ls(pj))
      call major_s(p,n,k,0,a,la,ws(lu1),lws(ll1),ws(na1),ws(nb1),
     *  ws(irh1),ws(ka1),ws(kb1),ws(kc1),lws(lv1),ws,lws,lp(lev),nm,
     *  w,ls,sgs,pp,smx,e(p))
      call signs(lp(lev),nm,ls(pj).lt.0,w,ls,tol)
c  check consistency of r(p)
      rp=scpr(0.D0,g,ws(na1),n)
      if(ls(pj).lt.0)rp=-rp
      if(abs(rp-r(p)).ge.-sgnf*r(p))then
        call refine_s(p,nk,n,a,la,ls,ws(na1),ws(nb1),ws(lu1),lws(ll1))
c       write(nout,*)'norm_ds =',sqrt(scpr(0.D0,ws(nb1),ws(nb1),n))
        rp=scpr(0.D0,g,ws(na1),n)
        if(ls(pj).lt.0)rp=-rp
        if(iprint.ge.1)
     *    write(nout,*)'new and old values of r(p)',rp,r(p),e(p)
        if(-rp.le.gtol*e(p))then
c         if(iprint.ge.1)write(nout,*)'refined r(p) is negligible',rp
          r(p)=0.D0
          goto82
        endif
c       if(abs(rp-r(p)).ge.-sgnf*r(p))goto98
        if(abs(rp-r(p)).ge.-sgnf*r(p))then
          if(iprint.ge.1)write(nout,*)'refined r(p) is negligible',rp
          r(p)=0.D0
          goto82
        endif
        r(p)=rp
      endif
      rp=r(p)
   86 continue
      if(iprint.ge.3)then
        write(nout,1001)'residual vector and indices',
     *    (ls(j),r(abs(ls(j))),j=lp(lev),nm)
        write(nout,1000)'denominators',(w(abs(ls(j))),j=lp(lev),nm)
      endif
      sml=smallish
   88 continue
c  ratio test at higher levels
      alpha=1.D37
      qqj=0
      do 90 j=lp(lev),nm
        i=abs(ls(j))
        wi=w(i)
        if(wi.le.0.D0)goto90
        if(r(i).lt.0.D0)goto90
        z=r(i)/wi
        if(z.ge.alpha)goto90
        alpha=z
        qqj=j
   90 continue
      if(qqj.eq.0)then
        do j=lp(lev),nm
          r(abs(ls(j)))=0.D0
        enddo
        jmin=max(nk+1,lp(lev-1))
        do j=jmin,lp(lev)-1
          i=abs(ls(j))
          if(i.le.n)then
            w(i)=ws(na+i)
          else
            w(i)=aiscpr(n,a,la,i-n,ws(na1),0.D0)
          endif
        enddo
        call signs(jmin,lp(lev)-1,ls(pj).lt.0,w,ls,tol)
        lev=lev-1
        f=alp(lev)
        if(iprint.ge.2)write(nout,*)'UNBOUNDED:   p =',p,
     *    '   return to level',lev
        if(lev.gt.1)goto86
        degen=.true.
        if(kmax.ne.0.and.ninf.eq.0)goto30
        goto35
      endif
      qq=abs(ls(qqj))
c  test significance of w(qq)
      if(abs(w(qq)).le.sml)then
        call refine_s(p,nk,n,a,la,ls,ws(na1),ws(nb1),ws(lu1),lws(ll1))
c       write(nout,*)'norm_ds =',sqrt(scpr(0.D0,ws(nb1),ws(nb1),n))
        wqq=w(qq)
        sml=small
        call reset_w(lp(lev),nm,ls(pj).lt.0,n,a,la,w,ls,ws(na1),sml)
        if(iprint.ge.1)
     *    write(nout,*)'new and old values of w(qq)',w(qq),wqq,qq
        if(iprint.ge.3)write(nout,1000)
     *    'refined denominators',(w(abs(ls(j))),j=lp(lev),nm)
c       if(w(qq).eq.0.D0)goto88
        if(wqq*w(qq).le.0.D0)then
          w(qq)=0.D0
          goto88
        endif
        if(nup.gt.0.and.abs(wqq-w(qq)).ge.sgnf*abs(w(qq)))then
          call refactor(n,nm,a,la,ws(lu1),lws(ll1),ifail)
          if(iprint.gt.1)write(nout,*)'refactor: ifail =',ifail
          if(ifail.eq.1)goto98
          if(ifail.gt.0)return
        endif
c       if(abs(wqq-w(qq)).ge.sgnf*abs(w(qq)))goto98
      endif
      ff=f+alpha*rp
      if(iprint.ge.2)then
        write(nout,*)'alpha =',alpha,'   p =',p,'   qq =',qq
        write(nout,*)'r(p),r(qq),w(qq) =',r(p),r(qq),w(qq)
      endif
c  test for equality c/s
      if(bu(qq)-bl(qq).le.tol)then
        do j=lp(2),nm
          r(abs(ls(j)))=0.D0
        enddo
        lev=1
        f=alp(1)
        ff=f
        alpha=0.D0
        eqty=.true.
        if(iprint.ge.2)write(nout,*)'EQTY:   p =',p,'   qq =',qq,
     *    '   return to level 1'
        goto45
      endif
      if(ff.ge.f)then
c  potential degeneracy block at level lev
        if(r(qq).gt.tol)then
          if(iprint.gt.1)write(nout,*)'** truncate r(p)'
          r(p)=0.D0
          goto82
        endif
        if(lev+2.gt.mlp)then
          ifail=5
          return
        endif
        r(qq)=0.D0
        plev=nm+1
        do j=nm,lp(lev),-1
          i=abs(ls(j))
          if(r(i).eq.0.D0)then
            plev=plev-1
            call iexch(ls(j),ls(plev))
            if(bu(i)-bl(i).gt.tol)r(i)=1.D0
          endif
        enddo
        alp(lev)=f
        f=0.D0
        lev=lev+1
        lp(lev)=plev
        if(iprint.ge.2)write(nout,*)
     *    'degeneracy: increase level to ',lev       
        if(iprint.ge.1)write(nout,'(''pivots ='',I5,A,''level ='',
     *    I2,''    f ='',E16.8)')info(1),spaces(:3*lev-1),lev,f
        goto86
      endif
      iter=iter+1
      if(iprint.ge.2)write(nout,*)'replace',p,' by',qq
      call pivot(p,qq,n,nmi,a,la,e,ws(lu1),lws(ll1),ifail,info)
      if(ifail.ge.1)then
        if(ifail.ge.2)return
c       call iexch(ls(pj),ls(qqj))
        if(iprint.ge.1)write(nout,*)'near singularity in pivot (4)'
        goto98
      endif
c  update r and f
      do j=lp(lev),nm
        i=abs(ls(j))
        ri=r(i)-alpha*w(i)
        if(abs(ri).le.tol)ri=0.D0
        r(i)=ri
      enddo
      f=ff
c  exchange a constraint
      r(p)=alpha
      if(r(p).le.tol)r(p)=0.D0
      call iexch(ls(pj),ls(qqj))
      if(iprint.ge.1)write(nout,'(''pivots ='',I5,A,''level ='',
     *    I2,''    f ='',E16.8)')info(1),spaces(:3*lev-1),lev,f
      goto80
c  restart sequence
   98 continue
      k=0
      nk=n
      mode=2
      mres=mres+1
      if(iprint.ge.1)write(nout,*)'major restart #',mres
   99 continue
      if(lev.gt.1)f=alp(1)
      if(ninf.eq.0)then
        if(phase.eq.1.or.f.le.bestf-smallish*gnorm)then
          bestf=f
          phase=2
          ires=0
        else
          ires=ires+1
        endif
      elseif(phase.eq.1.and.f.lt.bestf-smallish*gnorm)then
        bestf=f
        ires=0
      else
        ires=ires+1
      endif
      if(ires.le.nres)goto3
      ifail=8
      return
 1000 format(a/(e16.5,4e16.5))
 1001 format(a/(i4,1x,e11.5,4(i4,1x,e11.5)))
c1000 format(a/(e18.8,3e19.8))
c1001 format(a/(i3,1x,e14.8,3(i4,1x,e14.8)))
      end

      block data defaults
      implicit double precision (a-h,o-z)
      common/epsc/eps,tol,emin
      common/repc/sgnf,nrep,npiv,nres
      common/refactorc/nup,nfreq
      common/wsc/kk,ll,kkk,lll,mxws,mxlws
c  single length tolerances
c     data  eps,   tol, emin, sgnf, nrep, npiv, nres, nfreq, kk, ll
c    *   /6.D-8, 1.D-6, 0.D0, 1.D-1,  2,    3,    3,   100,   0,  0/
c  double length tolerances
      data  eps,    tol,   emin, sgnf, nrep, npiv, nres, nfreq, kk, ll 
     * /1111.D-19, 1.D-12, 0.D0, 1.D-4,  2,    3,   2,   500,   0,  0/
      end

      subroutine ratio(jmin,jmax,r,w,bl,bu,ls,amax,alpha,qqj,qqj1)
c  two-sided ratio test
      implicit double precision (a-h,r-z), integer (i-q)
      common/noutc/nout
      dimension r(*),w(*),bl(*),bu(*),ls(*)
c     write(nout,*)'qqj,alpha,amax',qqj,alpha,amax
      alpha1=1.D37
      qqj1=0
      do 1 j=jmin,jmax
        i=abs(ls(j))
        wi=w(i)
        if(wi.eq.0.D0)goto1
        ri=r(i)
        if(wi.gt.0.D0)then
          if(ri.lt.0.D0)goto1
          z=ri/wi
        else
          if(ri.lt.0.D0)then
            z=ri/wi
            if(z.lt.alpha1)then
              alpha1=z
              qqj1=j
            endif
          endif
          z=(bl(i)-bu(i)+ri)/wi
        endif
        if(z.ge.amax)goto1
        amax=z
        qqj=j
    1 continue
c     write(nout,*)'qqj,qqj1,alpha1,amax',qqj,qqj1,alpha1,amax
      if(qqj1.gt.0.and.alpha1.le.amax)then
        do 2 j=jmin,jmax
          i=abs(ls(j))
          wi=w(i)
          if(wi.ge.0.D0)goto2
          ri=r(i)
          if(ri.lt.0.D0)then
            z=ri/wi
            if(z.gt.alpha1.and.z.le.amax)then
              alpha1=z
              qqj1=j
            endif
          endif
    2   continue
        alpha=alpha1
        qqj=qqj1
      else
        qqj1=0
        alpha=amax
      endif
c     write(nout,*)'alpha,amax',alpha,amax
      return
      end

      subroutine update(n,jmin,nm,a,la,g,r,w,x,bl,bu,ls,
     *  alpha,n_inf,f,tol)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),g(*),r(*),w(*),x(*),bl(*),bu(*),ls(*)
      n_inf=0
      f=0.D0
      do j=jmin,nm
        i=abs(ls(j))
        y=alpha*w(i)
        if(y.ne.0.D0)then
          if(i.le.n)then
            if(ls(j).gt.0)then
              x(i)=x(i)-y
            else
              x(i)=x(i)+y
            endif
          endif
          ri=r(i)-y
          s=sign(1.D0,dble(ls(j)))
          if(y.lt.0.D0)then
            ro=bu(i)-bl(i)-ri
            if(ro.lt.ri)then
              ri=max(ro,0.D0)
              w(i)=-w(i)
              ls(j)=-ls(j)
            endif
          endif
          if(abs(ri).le.tol)ri=0.D0
          if(r(i).lt.0.D0)then
            if(ri.ge.0.D0)then
c  remove contribution to gradient
              if(i.gt.n)then
                call saipy(s,a,la,i-n,g,n)
              else
                g(i)=g(i)+s
              endif
            else
              n_inf=n_inf+1
              f=f-ri
            endif
          endif
          r(i)=ri
        elseif(r(i).lt.0.D0)then
          n_inf=n_inf+1
          f=f-r(i)
        endif
      enddo
      return
      end

      subroutine optest(jmin,jmax,r,e,ls,rp,pj)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension r(*),e(*),ls(*)
      rp=0.D0
      do 1 j=jmin,jmax
        i=abs(ls(j))
        if(e(i).eq.0.D0)print *,'e(i).eq.0.D0: i =',i
        ri=r(i)/e(i)
        if(ri.ge.rp)goto1
        rp=ri
        pj=j
    1 continue
      return
      end

      subroutine signs(jmin,jmax,plus,w,ls,tol)
      implicit double precision (a-h,o-z)
      dimension w(*),ls(*)
      logical plus
c  change signs as necessary
      do 1 j=jmin,jmax
        i=abs(ls(j))
        wi=w(i)
        if(wi.eq.0.D0)goto1
        if(abs(wi).le.tol)then
          w(i)=0.D0
          goto1
        endif
        if(plus.neqv.ls(j).ge.0)then
          w(i)=-wi
        endif
    1 continue
      return
      end

      subroutine signst(jmin,jmax,plus,r,w,e,ls,tol)
      implicit double precision (a-h,o-z)
      dimension r(*),w(*),e(*),ls(*)
      logical plus
c  transfer with sign change as necessary
      do j=jmin,jmax
        i=abs(ls(j))
        if(abs(w(i)).le.e(i)*tol)then
          r(i)=0.D0
        elseif(plus.eqv.ls(j).ge.0)then
          r(i)=w(i)
        else
          r(i)=-w(i)
        endif
      enddo
      return
      end

      subroutine stmap(n,nm,kmax)
c  set storage map for workspace in bqpd and auxiliary routines
      implicit double precision (a-h,r-z), integer (i-q)
      common/wsc/kk,ll,kkk,lll,mxws,mxlws
      common/bqpdc/irh1,na,na1,nb,nb1,ka1,kb1,kc1,irg1,lu1,lv,lv1,ll1
c  real storage (ws)
      kkk=kmax*(kmax+9)/2+nm+n
c  (number of real locations required by bqpd)
c  slot for user workspace for gdotx
      irh1=kk+1
c  slot of length kmax*(kmax+1)/2 for reduced Hessian matrix
      na1=irh1+kmax*(kmax+1)/2
      na=na1-1
c  scratch slots of length n+m, n and four of length kmax follow
      nb1=na1+nm
      nb=nb1-1
      ka1=nb1+n
      kb1=ka1+kmax
      kc1=kb1+kmax
      irg1=kc1+kmax
      lu1=irg1+kmax
c  remaining space for use by denseL.f or sparseL.f
c  integer storage (lws)
      lll=kmax
c  (number of integer locations required by bqpd)
c  slot for user workspace for gdotx
      lv=ll
      lv1=ll+1
c  slot for V-list
      ll1=lv1+kmax
c  remaining space for use by denseL.f or sparseL.f
      return
      end

      subroutine check(n,nk,nm,kmax,g,a,la,x,bl,bu,r,ls,an,f,
     *  ws,lws,ninf,peq,lp1,nk1,lev,p,alp2)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension g(*),a(*),la(*),x(*),bl(*),bu(*),r(*),ls(*),
     *  an(*),ws(*),lws(*)
      common/noutc/nout
      common/epsc/eps,tol,emin
      if(lev.eq.2)then
        do i=1,n
          an(i)=g(i)
        enddo
        e=alp2*sign(1.D0,dble(p))
        i=abs(p)
        if(i.le.n)then
          an(i)=an(i)-e
        else
          call saipy(-e,a,la,i-n,an,n)
        endif
        goto1
      endif
      j=nm*(nm+1)/2
      do i=1,nm
        j=j-abs(ls(i))
      enddo
      if(j.ne.0)write(nout,*)'indexing error'
      e=0.
      do j=nk+1,nm
        i=abs(ls(j))
        if(i.le.n)then
          s=x(i)
        else
          s=aiscpr(n,a,la,i-n,x,0.D0)
        endif
        if(ls(j).gt.0)then
          s=r(i)-s+bl(i)
        else
          s=r(i)+s-bu(i)
        endif
        if(abs(s).le.tol*max(1.D0,abs(r(i))))s=0.D0
        if(abs(s).gt.e)then
          e=abs(s)
          ie=i
        endif
      enddo
      if(e.gt.tol)write(nout,*)'residual error at level 1 = ',e,ie
      if(ninf.eq.0)then
        call setfg2(n,kmax,a,la,x,ff,an,ws,lws)
      else
        call setfg1(n,nm,a,la,x,bl,bu,ff,fb,an,peq,lp1,nk1,ls,r)
      endif
      gnm=sqrt(scpr(0.D0,an,an,n))
      if(lev.eq.1)then
        e=abs(ff-f)
        if(e.gt.tol*max(1.D0,abs(f)))write(nout,*)'function error = ',e,
     *    '   f(x) =',ff
      endif
    1 continue
      e=0.D0
      do j=1,nk
c       write(nout,*)'an =',(an(i),i=1,n)
        i=abs(ls(j))
        s=sign(1.D0,dble(ls(j)))
        if(i.le.n)then
          an(i)=an(i)-s*r(i)
          if(ls(j).gt.0)then
            s=x(i)-bl(i)
          else
            s=bu(i)-x(i)
          endif
        else
          call saipy(-s*r(i),a,la,i-n,an,n)
          if(ls(j).gt.0)then
            s=aiscpr(n,a,la,i-n,x,-bl(i))
          else
            s=-aiscpr(n,a,la,i-n,x,-bu(i))
          endif
        endif
        if(abs(s).gt.e)then
          if(j.gt.peq.and.j.lt.lp1)s=0.
          e=abs(s)
          ie=i
        endif
      enddo
      if(e.gt.tol)write(nout,*)'residual error at level 2 = ',e,ie
      e=0.D0
      do i=1,n
        if(abs(an(i)).gt.e)then
          e=abs(an(i))
          ie=i
        endif
      enddo
      if(e.gt.gnm*tol)write(nout,*)'KT condition error = ',e,ie,gnm
c     if(e.gt.gnm*tol)write(nout,*)'KT cond_n errors = ',(an(i),i=1,n)
      return
      end

      subroutine warm_start(n,nm,nk,a,la,x,bl,bu,b,ls,lv,aa,ll,an,vstep)
      implicit double precision (a-h,o-z)
      dimension a(*),la(*),x(*),bl(*),bu(*),b(*),ls(*),lv(*),
     *  aa(*),ll(*),an(*)
      DOUBLE PRECISION daiscpr
      common/epsc/eps,tol,emin
      common/noutc/nout
      common/iprintc/iprint
      do j=1,nk
        i=abs(ls(j))
        if(i.le.n)then
          b(i)=0.D0
          ri=x(i)-bl(i)
          ro=bu(i)-x(i)
          if(ri.le.ro)then
            ls(j)=i
            if(abs(ri).le.tol)x(i)=bl(i)
          else
            ls(j)=-i
            if(abs(ro).le.tol)x(i)=bu(i)
          endif
        endif
      enddo
c     write(nout,*)'ls(1:nk) =',(ls(i),i=1,nk)
      do j=1,nk
        i=abs(ls(j))
        if(i.gt.n)then
          if(ls(j).ge.0)then
            b(i)=daiscpr(n,a,la,i-n,x,-bl(i))
          else
            b(i)=daiscpr(n,a,la,i-n,x,-bu(i))
          endif
        endif
      enddo
      do j=1,n-nk
        b(lv(j))=0.D0
      enddo
c     write(nout,1000)'x =',(x(i),i=1,n)
c     write(nout,1000)'r =',(b(abs(ls(j))),j=1,nk)
      call tfbsub(n,a,la,0,b,an,aa,ll,ep,.false.)
c     write(nout,1000)'d =',(an(i),i=1,n)
      call linf(n,an,vstep,i)
      if(iprint.ge.1)
     *  write(nout,*)'infinity norm of vertical step =',vstep
      do j=1,nk
        i=abs(ls(j))
        if(i.le.n)an(i)=0.D0
      enddo
      do i=1,n
        x(i)=x(i)-an(i)
      enddo
c     write(nout,*)'x =',(x(i),i=1,n)
 1000 format(a/(e16.5,4e16.5))
c1001 format(a/(i4,1x,e11.5,4(i4,1x,e11.5)))
      return
      end

      subroutine residuals(n,jmin,nm,a,la,x,bl,bu,r,ls,qj,ninf)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),x(*),bl(*),bu(*),r(*),ls(*)
      common/epsc/eps,tol,emin
      DOUBLE PRECISION daiscpr,z
      rq=0.D0
      ninf=0
      do j=jmin,nm
        i=abs(ls(j))
        if(i.gt.n)then
          z=daiscpr(n,a,la,i-n,x,0.D0)
        else
          z=dble(x(i))
        endif
        ri=z-dble(bl(i))
        ro=dble(bu(i))-z
        if(ri.le.ro)then
          ls(j)=i
        else
          ri=ro
          ls(j)=-i
        endif
        if(abs(ri).le.tol)ri=0.D0
        if(ri.lt.0.D0)then
          ninf=ninf+1
          if(ri.lt.rq)then
            qj=j
            rq=ri
          endif
        endif
        r(i)=ri
      enddo
      return
      end

      subroutine setrp(psb,plen,alpha,xp,blp,bup,gp,rpu,lspj,infb,fb,
     *  tol)
      implicit double precision (a-h,o-z)
      logical psb,plen
      if(psb)then
        t=sign(alpha,dble(lspj))
        if(xp.lt.blp)then
          fb=fb+xp-blp
          xp=xp+t
          lspj=abs(lspj)
          rpu=xp-blp
          infb=infb-1
          if(fb.le.tol)fb=0.D0
          if(rpu.ge.-tol)then
            xp=blp
            gp=gp+1.D0
          endif
        elseif(xp.gt.bup)then
c         fp=fp+bup-xp
          fb=fb+bup-xp
          xp=xp+t
          lspj=-abs(lspj)
          rpu=bup-xp
          infb=infb-1
          if(fb.le.tol)fb=0.D0
          if(rpu.ge.-tol)then
            xp=bup
            gp=gp-1.D0
          endif
        else
          xp=xp+t
          rpu=bup-xp
          t=xp-blp
          if(t.le.rpu)then
            rpu=max(t,0.D0)
            lspj=abs(lspj)
          else
            rpu=max(rpu,0.D0)
            lspj=-abs(lspj)
          endif
        endif
      else
        if(plen)xp=xp+sign(alpha,dble(lspj))
        rpu=max(bup-blp-alpha,0.D0)
        if(alpha.le.rpu)then
          rpu=alpha
          lspj=lspj
        else
          lspj=-lspj
        endif
      endif
      if(abs(rpu).le.tol)rpu=0.D0
      return
      end

      subroutine refine_s(p,nk,n,a,la,ls,s,ds,aa,ll)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),ls(*),s(*),ds(*),aa(*),ll(*)
      common/noutc/nout
      DOUBLE PRECISION daiscpr
c     write(nout,*)'refine_s'
      do j=1,nk
        i=abs(ls(j))
        if(i.gt.n)then
          s(i)=daiscpr(n,a,la,i-n,s,0.D0)
        endif
      enddo
      if(p.le.n)then
        s(p)=0.D0
      else
        s(p)=daiscpr(n,a,la,p-n,s,-1.D0)
      endif
      call tfbsub(n,a,la,0,s,ds,aa,ll,ep,.false.)
      do i=1,n
        s(i)=s(i)-ds(i)
      enddo
      if(p.le.n)s(p)=1.D0
      return
      end

      subroutine reset_w(jmin,jmax,plus,n,a,la,w,ls,s,sml)
      implicit double precision (a-h,o-z)
      dimension a(*),la(*),w(*),ls(*),s(*)
      logical plus
      do 1 j=jmin,jmax
        i=abs(ls(j))
        wi=w(i)
        if(wi.eq.0.D0)goto1
        if(i.le.n)then
          z=s(i)
        else
          z=aiscpr(n,a,la,i-n,s,0.D0)
        endif
        if(abs(z).le.sml)then
          w(i)=0.D0
          goto1
        endif
        if(plus.neqv.ls(j).ge.0)then
          w(i)=-z
        else
          w(i)=z
        endif
    1 continue
      return
      end

      subroutine setfg1(n,nm,a,la,x,bl,bu,f,fb,g,peq,lp1,nk1,ls,r)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),x(*),bl(*),bu(*),g(*),ls(*),r(*)
      do i=1,n
        g(i)=0.D0
      enddo
      fb=0.D0
      do j=peq+1,lp1-1
        i=abs(ls(j))
        if(x(i).lt.bl(i))then
          fb=fb+bl(i)-x(i)
          g(i)=g(i)-1.D0
        elseif(x(i).gt.bu(i))then
          fb=fb+x(i)-bu(i)
          g(i)=g(i)+1.D0
        endif
      enddo
      f=fb
      do j=nk1,nm
        i=abs(ls(j))
        if(r(i).lt.0.D0)then
          f=f-r(i)
          if(i.le.n)then
            g(i)=g(i)-sign(1.D0,dble(ls(j)))
          else
            sign_=sign(1.D0,dble(ls(j)))
            call saipy(-sign_,a,la,i-n,g,n)
          endif
        endif
      enddo
      return
      end

      subroutine setfg2(n,kmax,a,la,x,f,g,ws,lws)
      implicit double precision (a-h,o-z)
      dimension a(*),la(*),x(*),g(*),ws(*),lws(*)
      common/noutc/nout
      if(kmax.eq.0)then
        do i=1,n
          g(i)=0.D0
        enddo
        call saipy(1.D0,a,la,0,g,n)
        f=scpr(0.D0,x,g,n)
      else
        call gdotx(n,x,ws,lws,g)
        call saipy(1.D0,a,la,0,g,n)
        f=5.D-1*scpr(aiscpr(n,a,la,0,x,0.D0),g,x,n)
      endif
      return
      end


C  --- Begin file: sparseA.f

c  Copyright (C) 1996 Roger Fletcher

c  Current version dated 9 December 2010

c  THE ACCOMPANYING PROGRAM IS PROVIDED UNDER THE TERMS OF THE ECLIPSE PUBLIC
c  LICENSE ("AGREEMENT"). ANY USE, REPRODUCTION OR DISTRIBUTION OF THE PROGRAM
c  CONSTITUTES RECIPIENT'S ACCEPTANCE OF THIS AGREEMENT

c  ******************************************
c  Specification of A in sparse matrix format
c  ******************************************

c  The matrix A contains gradients of the linear terms in the objective
c  function (column 0) and the general constraints (columns 1:m).
c  No explicit reference to simple bound constraints is required in A.
c  The information is set in the parameters a(*) (double precision real) and
c  la(*) (integer).

c  In this sparse format, these vectors have dimension  a(1:maxa)  and
c   la(0:maxla-1),  where maxa is at least nnza (the number of nonzero elements
c  in A), and maxla is at least  nnza+m+3.  la(0) and the last m+2 elementss
c  in la are pointers.

c  The vectors a(.) and la(.) must be set as follows:

c  a(j) and la(j) for j=1,nnza are set to the values and row indices (resp.)
c  of all the nonzero elements of A. Entries for each column are grouped
c  together in increasing column order. Within each column group, it is
c  not necessary to have the row indices in increasing order.

c  la(0) is a pointer which points to the start of the pointer information in
c  la. la(0) must be set to nnza+1 (or a larger value if it is desired to
c  allow for future increases to nnza).

c  The last m+2 elements of la(.) contain pointers to the first elements in
c  each of the column groupings. Thus la(la(0)+i)) for i=0,m is set to the
c  location in a(.) containing the first nonzero element for column i of A.
c  Also la(la(0)+m+1)) is set to nnza+1 (the first unused location in a(.)).

c  Copyright, University of Dundee (R.Fletcher), June 1996
c  Current version dated 31/01/07

      subroutine saipy(s,a,la,i,y,n)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),y(*)
c  saxpy with column i of A
      if(s.eq.0.D0)return
      jp=la(0)+i
      do j=la(jp),la(jp+1)-1
        ir=la(j)
        y(ir)=y(ir)+s*a(j)
      enddo
    1 format(A,15I2)
    2 format(A,5E15.7)
    3 format(A/(20I4))
    4 format(A/(5E15.7))

      return
      end

      subroutine msaipy(s,a,la,i,y,n)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),y(*)
c  saxpy with modulus of column i of A
      jp=la(0)+i
      do j=la(jp),la(jp+1)-1
        ir=la(j)
        y(ir)=y(ir)+s*abs(a(j))
      enddo
      return
      end

c     subroutine daipy(s,a,la,i,y,n)
c     DOUBLE PRECISION a(*),y(*),d
c     dimension la(0:*)
c     if(s.eq.0.D0)return
c     d=dble(s)
c     jp=la(0)+i
c     do j=la(jp),la(jp+1)-1
c       ir=la(j)
c       y(ir)=y(ir)+d*dble(a(j))
c     enddo
c     return
c     end

      subroutine isaipy(s,a,la,i,y,n,lr,li)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),y(*),lr(*),li(*)
c  indirectly addressed saxpy with column i of A
      if(s.eq.0.D0)return
      jp=la(0)+i
      do j=la(jp),la(jp+1)-1
        ir=li(la(j))
        y(ir)=y(ir)+s*a(j)
      enddo
      return
      end

c the old isaipy was what might be called isaipy2

      subroutine isaipy1(s,a,la,i,y,n,lr,li,m1)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),y(*),lr(*),li(*)
c  indirectly addressed saxpy with column i of A_1
      if(s.eq.0.D0)return
      jp=la(0)+i
      do j=la(jp),la(jp+1)-1
        ir=li(la(j))
        if(ir.le.m1)y(ir)=y(ir)+s*a(j)
      enddo
      return
      end

c     subroutine ssaipy(s,a,la,i,y,n)
c     implicit double precision (a-h,o-z)
c     dimension a(*),la(0:*),y(*)
c  saxpy with squares of column i of A
c     if(s.eq.0.D0)return
c     jp=la(0)+i
c     do j=la(jp),la(jp+1)-1
c       ir=la(j)
c       y(ir)=y(ir)+s*(a(j))**2
c     enddo
c     return
c     end

      function aiscpr(n,a,la,i,x,b)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),x(*)
c  scalar product with column i of A
      aiscpr=b
      jp=la(0)+i
      do j=la(jp),la(jp+1)-1
        ir=la(j)
        aiscpr=aiscpr+x(ir)*a(j)
      enddo
      return
      end

      function daiscpr(n,a,la,i,x,b)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),x(*)
      DOUBLE PRECISION daiscpr
      daiscpr=dble(b)
      jp=la(0)+i
      do j=la(jp),la(jp+1)-1
        ir=la(j)
        daiscpr=daiscpr+dble(x(ir))*dble(a(j))
      enddo
      return
      end

      function aiscpri(n,a,la,i,x,b,lr,li)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),x(*),lr(*),li(*)
c  indirectly addressed scalar product with column i of A
      aiscpri=b
      jp=la(0)+i
      do j=la(jp),la(jp+1)-1
        ir=li(la(j))
        aiscpri=aiscpri+x(ir)*a(j)
      enddo
      return
      end

      function daiscpri(n,a,la,i,x,b,lr,li)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),x(*),lr(*),li(*)
      DOUBLE PRECISION daiscpri
      daiscpri=dble(b)
      jp=la(0)+i
      do j=la(jp),la(jp+1)-1
        ir=li(la(j))
        daiscpri=daiscpri+dble(x(ir))*dble(a(j))
      enddo
      return
      end

c the old aiscpri was what might be called aiscpri2

      function aiscpri1(n,a,la,i,x,b,lr,li,m1)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),x(*),lr(*),li(*)
c  indirectly addressed scalar product with column i of A_1
      aiscpri1=b
      jp=la(0)+i
      do j=la(jp),la(jp+1)-1
        ir=li(la(j))
        if(ir.le.m1)aiscpri1=aiscpri1+x(ir)*a(j)
      enddo
      return
      end

      function ailen(n,a,la,i)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*)
c  L2 length of column i of A
      ailen=0.D0
      jp=la(0)+i
      do j=la(jp),la(jp+1)-1
        ailen=ailen+a(j)**2
      enddo
      ailen=sqrt(ailen)
      return
      end

      subroutine iscatter(a,la,i,li,an,n)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),li(*),an(*)
c  indirect scatter into previously zeroed vector an
      jp=la(0)+i
      do j=la(jp),la(jp+1)-1
        an(li(la(j)))=a(j)
      enddo
      return
      end

      subroutine iunscatter(a,la,i,li,an,n)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),li(*),an(*)
c  undo effect of iscatter
      jp=la(0)+i
      do j=la(jp),la(jp+1)-1
        an(li(la(j)))=0.D0
      enddo
      return
      end

      function aij(i,j,a,la)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*)
c  get element A(i,j)
      jp=la(0)+j
      do ij=la(jp),la(jp+1)-1
        ir=la(ij)
        if(ir.eq.i)then
          aij=a(ij)
          return
        endif
      enddo
      aij=0.D0
      return
      end

      subroutine setaij(aij,i,j,a,la)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*)
c  set element A(i,j)
      jp=la(0)+j
      do jj=la(jp+1)-1,la(jp),-1
        ir=la(jj)
        if(ir.eq.i)then
          a(jj)=aij
          return
        endif
      enddo
      if(aij.eq.0.D0)return
      print *,'malfunction: no slot for A(i,j)'
      stop
      end

      subroutine cscale(n,m,a,la,x,bl,bu,s,menu,ifail)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),x(*),bl(*),bu(*),s(*)

c     Constraint scaling procedure for use prior to calling bqpd when using
c     sparseA.f

c     Parameters are set as for bqpd, except for s, menu and ifail

c     The user must set the parameter menu to control how the
c     x-variables are scaled (or equivalently how constraints i = 1:n
c     are scaled), as follows

c     menu = 1 indicates that a unit scaling applies to the x-variables

c     menu = 2 the user provides estimates s(i)>0 of the magnitude of
c              x(i) for i = 1:n. In this case the elements  x(i), bl(i), bu(i)
c              are divided by s(i) for i = 1:n.

c     In all cases, cscale goes on to scale the general constraints, in
c     such a way that the normal vector of each nontrivial constraint in
c     the scaled problem has an l_2 norm of unity. This scaling is also
c     applied to the right hand sides  bl(i), bu(i) for i = n+1:n+m.
c     The scaled data overwrites the original data.

c     cscale also scales the constant vector of the quadratic function,
c     which is found in a(1:n). However if a non-unit x-variable scaling
c     is used, it is necessary for the user to scale the Hessian matrix
c     G appropriately. This can be done by passing the x-variable scale
c     factors s(i) i = 1:n into the subroutine gdotx using the
c     parameter ws, and multiplying G(i,j) by s(i)*s(j) (possibly
c     implicitly).

c     cscale sets ifail = 1 to indicate that some s(i)< = 0,
c             and ifail = 2 to indicate an incorrect setting of menu.
c       Otherwise ifail = 0.

      integer pjp

      ifail=2
      if(menu.lt.1.or.menu.gt.2)return
      pjp=la(0)
c     z=1.D0/log(2.D0)
      if(menu.eq.1)then
        do j=1,n
          s(j)=1.D0
        enddo
      else
        ifail=1
        do j=1,n
          if(s(j).le.0.D0)return
        enddo
c       if(menu.eq.2)then
c         do j=1,n
c           s(j)=2.D0**nint(log(s(j))*z)
c         enddo
c       endif
        do j=1,n
          if(s(j).ne.1.D0)then
            x(j)=x(j)/s(j)
            bl(j)=bl(j)/s(j)
            bu(j)=bu(j)/s(j)
          endif
        enddo
        do j=1,la(pjp+1)-1
          a(j)=a(j)*s(la(j))
        enddo
      endif
      do i=1,m
        t=0.D0
        do j=la(pjp+i),la(pjp+i+1)-1
          a(j)=s(la(j))*a(j)
          t=t+a(j)**2
        enddo
        t=sqrt(t)
        if(t.eq.0.D0)then
          s(n+i)=1.D0
        else
c         t=2.D0**nint(log(t)*z)
          s(n+i)=t
          do j=la(pjp+i),la(pjp+i+1)-1
            a(j)=a(j)/t 
          enddo
          bl(n+i)=bl(n+i)/t
          bu(n+i)=bu(n+i)/t
        endif
      enddo
      ifail=0
      return
      end

      subroutine modify(n,m,sigma,s,a,la,maxa,iws)
      implicit double precision (a-h,o-z)
      dimension s(*),a(*),la(0:*),iws(*)

c  Modifies the sparse data structure to add an extra variable and duplicate
c  the general constraints, to enable scaled L-infinity QPs to be solved.
c  Scale factors given in s(1:m) and the coefficient of the objective function in sigma
c  For unscaled problems set s=ones and sigma=1.
c  Needs m+1 locations of integer workspace in iws(*)

      n1=n+1
      n1m=n1+m
      m1=m+1
      la0=la(0)
      nextra=la(la0+m1)-la(la0+1)+m1+m
      ij=la(la0+m1)+nextra
c     print 1,'la0,nextra,ij',la0,nextra,ij
      if(ij-1.gt.maxa)then
        print *,'not enough space:  reset maxa to at least ',ij-1
        stop
      endif
      do i=1,m1
        iws(i)=la(la0+i)
      enddo
      la0=ij
      la(la0+m1+m)=ij
c  set lower bounds
      do i=m,1,-1
        ij=ij-1
c       a(ij)=-1.D0
        a(ij)=-s(i)
        la(ij)=n1
        do j=iws(i+1)-1,iws(i),-1
          ij=ij-1
          a(ij)=-a(j)
          la(ij)=la(j)
        enddo
        la(la0+m+i)=ij
      enddo
c  set upper bounds
      do i=m,1,-1
        ij=ij-1
c       a(ij)=-1.D0
        a(ij)=-s(i)
        la(ij)=n1
        do j=iws(i+1)-1,iws(i),-1
          ij=ij-1
          a(ij)=a(j)
          la(ij)=la(j)
        enddo
        la(la0+i)=ij
      enddo
      ij=ij-1
      la(ij)=n1
      a(ij)=sigma
      la(la0)=1
      la(0)=la0
c     print 3,'pointers =',(la(i),i=la0,la0+m1+m)
c     print 4,'a =',(a(i),i=1,la(la0+m1+m)-1)
c     print 3,'la =',(la(i),i=1,la(la0+m1+m)-1)
    1 format(A,15I4)
    2 format(A,5E15.7)
    3 format(A/(20I4))
    4 format(A/(5E15.7))
      return
      end

      subroutine restore(n,m,a,la)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*)

c  restores the changes made by subroutine modify

      la0=la(0)
      do i=1,m
        do j=la(la0+i),la(la0+i+1)-2
          la(j-i)=la(j)
          a(j-i)=a(j)
        enddo
        la(la0+i)=la(la0+i)-i
      enddo
      la(la0+m+1)=la(la0+m+1)-m-1
c     print 3,'pointers =',(la(i),i=la0,la0+m+1)
c     print 4,'a =',(a(i),i=1,la(la0+m+1)-1)
c     print 3,'la =',(la(i),i=1,la(la0+m+1)-1)
    3 format(A/(20I4))
    4 format(A/(5E15.7))
      return
      end

      subroutine extend_la(n,m,la,lax)
      implicit double precision (a-h,o-z)
      dimension la(0:*),lax(0:*)

c  Modifies the sparse data structure to add an extra variable and duplicate
c  the general constraints, to enable scaled L-infinity QPs to be solved.
c  The gradient vector is assumed to have a single non-zero entry (n+1).

      n1=n+1
      la0=la(0)
      lax0=2*(la(la0+m+1)-la(la0+1)+m+1)
      lax(0)=lax0
      lax(1)=n1
      lax(lax0)=1
      ijx=2
      jp=lax0+1
      do k=1,2
        do j=1,m
          lax(jp)=ijx
          do ij=la(la0+j),la(la0+j+1)-1
            lax(ijx)=la(ij)
            ijx=ijx+1
          enddo
          lax(ijx)=n1
          ijx=ijx+1
          jp=jp+1
        enddo
      enddo
      lax(jp)=ijx
c     print 3,'lax pointers =',(lax(i),i=lax0,lax0+m+m+1)
c     print 3,'lax =',(lax(i),i=1,lax(lax0+m+m+1)-1)
    3 format(A/(20I4))
      return
      end

      subroutine extend_a(n,m,a,la,ax,lax,s,sigma)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),ax(*),lax(0:*),s(*)

c  Extends the sparse data values to enable scaled L-infinity QPs to be solved.
c  Scale factors given in s(1:m) and the coefficient of the objective function in sigma
c  For unscaled problems set s=ones and sigma=1.

      n1=n+1
      la0=la(0)
      ax(1)=sigma
      ijx=2
c  set lower bounds
      do j=1,m
        do ij=la(la0+j),la(la0+j+1)-1
          ax(ijx)=-a(ij)
          ijx=ijx+1
        enddo
        ax(ijx)=-s(j)
        ijx=ijx+1
      enddo
c  set upper bounds
      do j=1,m
        do ij=la(la0+j),la(la0+j+1)-1
          ax(ijx)=a(ij)
          ijx=ijx+1
        enddo
        ax(ijx)=-s(j)
        ijx=ijx+1
      enddo
c     print 4,'ax =',(ax(i),i=1,lax(lax(0)+m+m+1)-1)
    4 format(A/(5E15.7))
      return
      end


C  --- Begin file: sparseL.f

cut here >>>>>>>>>>>>>>>>>

c***************** sparse matrix routines for manipulating L *******************

c           ***************************************************
c           Basis matrix routines for bqpd with sparse matrices
c           ***************************************************

c  These routines form and update L-Implicit-U factors LPB=U of a matrix B
c  whose columns are the normal vectors of the active constraints. In this
c  method only the unit lower triangular matrix L and the diagonal of U (in
c  addition to the row permutation P) is stored. B is represented in block form

c    | I  A_2 |    where the last m1 columns (A_2 and A_1) come from the
c    | 0  A_1 |    general constraint normals (columns of the matrix A in bqpd)

c  and the remaining unit columns come from simple bounds. The matrix A must be
c  specified in sparse format and the user is referred to the file  sparseA.f.

c  The data structure used for L is that of a profile or skyline scheme, in
c  which the nontrivial rows of L are stored as dense row spikes. The use of
c  a Tarjan+spk1 ordering algorithm to control the length of these spikes has
c  proved quite effective. The factors are updated by a variant of the
c  Fletcher-Matthews method, which has proved very reliable in practice.
c  However the B matrix is re-factored every 30 updates to control growth in
c  the total spike length.

c  Workspace
c  *********
c  The user needs to supply storage for the rows of L, although the amount
c  required is unknown a-priori.
c  sparse.f requires
c     5*n+nprof          locations of real workspace, and
c     9*n+m              locations of integer workspace
c  where nprof is the space required for storing the row spikes of the L matrix.
c  Storage for sparseL.f is situated at the end of the workspace arrays ws
c  and lws in bqpd.
c  Allow as much space for nprof as you can afford: the routine will report if
c  there is not enough. So far 10^6 locations has proved adequate for problems
c  of up to 5000 variables.

c  In addition the current version of bqpd.f requires
c     kmax*(kmax+9)/2+2*n+m   locations of real workspace in ws
c     kmax                    locations of integer workspace in lws
c  The user is also allowed to reserve storage in ws and lws, for use in the
c  user-supplied routine gdotx. This storage is situated at the start of the
c  arrays ws and lws. The user specifies the amount required by
c  setting the parameters kk and ll in the common block
c     common/wsc/kk,ll,kkk,lll,mxws,mxlws
c  The user MUST also set mxws and mxlws to be (respectively) the total amount
c  of real and integer workspace for the arrays ws and lws.

c  Other information
c  *****************

c  The methodology behind the L-Implicit-U factors and the row spike storage
c  scheme for L is described in the references
c    Fletcher R., Dense Factors of Sparse Matrices, in "Approximation Theory
c    and Optimization. Tributes to M.J.D. Powell", (M.D. Buhmann and A. Iserles,
c    eds), Cambridge University Press (1997), pp. 145-166.
c  and
c    Fletcher R., Block Triangular Orderings and Factors for Sparse Matrices
c    in LP, in "Numerical analysis 1997" (D.F. Griffiths, D.J. Higham and
c    G.A. Watson, eds.), Pitman Research Notes in Mathematics 380, (1998),
c    Longman, Harlow, pp. 91-110.

c  The file contains routines for solving systems with B or its transpose
c  which might be of use in association with bqpd. These routines are
c  documented below.

c  Steepest edge coefficients e(i) are also updated in these routines

c  Copyright, University of Dundee (R.Fletcher), January 1998
c  Current version dated 16/04/02

      subroutine start_up(n,nm,nmi,a,la,nk,e,ls,aa,ll,mode,ifail)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(0:*),e(*),ls(*),aa(*),ll(*)
      common/noutc/nout
      common/wsc/kk,ll_,kkk,lll,mxws,mxlws
      common/epsc/eps,tol,emin
      common/sparsec/ns,ns1,nt,nt1,nu,nu1,nx,nx1,np,np1,nprof,
     *  lc,lc1,li,li1,lm,lm1,lp,lp1,lq,lq1,lr,lr1,ls_,ls1,lt,lt1
      common/factorc/m1,m2,mp,mq,lastr,irow
      common/refactorc/nup,nfreq
      nfreq=min(30,nfreq)
      nup=0
      ns=kk+kkk+5*n
      nt=ll_+lll+8*n+nmi
      nprof=mxws-ns
      if(nprof.le.0.or.nt.gt.mxlws)then
        write(nout,*)'not enough real (ws) or integer (lws) workspace'
        write(nout,*)'you give values for mxws and mxlws as',mxws,mxlws
        write(nout,*)'minimum values for mxws and mxlws are',ns,nt
        ifail=7
        return
      endif
    3 format(A/(20I5))
    4 format(A/(5E15.7))
c  set storage map for sparse factors
      ns=n
      ns1=ns+1
      nt=ns+n
      nt1=nt+1
      nu=nt+n
      nu1=nu+1
      nx=nu+n
      nx1=nx+1
      np=nx+n
      np1=np+1
      lc=n
      lc1=lc+1
      li=lc+n
      li1=li+1
      lm=li+nmi
      lm1=lm+1
      lp=lm+n
      lp1=lp+1
      lq=lp+n
      lq1=lq+1
      lr=lq+n
      lr1=lr+1
      ls_=lr+n
      ls1=ls_+1
      lt=ls_+n
      lt1=lt+1
      m=nm-n
      mp=-1
      mq=-1
c     write(nout,*)'ls',(ls(ij),ij=1,nk)
      if(mode.eq.3)then
        call re_order(n,nm,a,la(1),la(la(0)),ll,ll(lc1),ll(li1),
     *    ll(lm1),ll(lp1),ll(lq1),ll(lr1),ll(ls1),ll(lt1),aa(np1),
     *    nprof,ifail)
        if(ifail.ge.1)then
c         write(nout,*)'failure in re_order (1)'
          if(ifail.eq.7)return
          mode=2
          goto1
        endif
        call re_factor(n,nm,a,la,ll,ll(lc1),ll(li1),
     *    ll(lm1),ll(lp1),ll(lq1),ll(lr1),ll(ls1),ll(lt1),aa(np1),
     *    nprof,aa,ifail)
        if(ifail.eq.7)return
        call check_L(n,aa,ll(lp1),ifail)
        if(ifail.eq.1)then
          mode=2
          goto1
        endif
        if(nk.eq.n)return
c  reset ls from e
        do j=1,nk
          i=-ls(j)
          if(i.gt.0)e(i)=-e(i)
        enddo
        j=0
        nk=nmi
        do i=1,nmi
          if(e(i).ne.0.D0)then
            j=j+1
            if(e(i).gt.0.D0)then
              ls(j)=i
            else
              ls(j)=-i
              e(i)=-e(i)
            endif
          else
            ls(nk)=i
            nk=nk-1
          endif
        enddo
        if(j.ne.n)then
          write(nout,*)'malfunction in reset sequence in start_up'
          stop
        endif
        return
      endif
    1 continue
      if(emin.eq.0.D0)then
c  set a lower bound on e(i): setting emin=0.D0 will force emin to be recalculated: do this only if mode<3
        emin=1.D0
        do i=1,nmi-n
          emin=max(emin,ailen(n,a,la,i))
        enddo
        emin=1.D0/emin
      endif
      do i=1,n
        ll(i)=i
        ll(li+i)=i
        e(i)=1.D0
      enddo
      do i=n+1,nmi
        ll(li+i)=0
        e(i)=0.D0
      enddo
      nu_=0
      if(mode.ne.0)then
c  shift designated bounds to end and order the resulting rows and columns
        do j=1,nk
          i=abs(ls(j))
          if(i.le.n)then
            nn=n-nu_
            nu_=nu_+1
            call iexch(ls(nu_),ls(j))
            ii=ll(li+i)
            ll(ii)=ll(nn)
            ll(li+ll(ii))=ii
            ll(nn)=i
            ll(li+i)=nn
          endif
        enddo
        call order(n,nu_,nk,la,ll,ls,ll(li1),ll(lp1),ll(lq1),ll(lr1),
     *    aa(np1),nprof,ifail)
        if(ifail.gt.0)return
      endif
      call factor(n,nmi,nu_,nk,a,la,e,ls,aa(ns1),aa(nt1),aa(nu1),
     *  aa(nx1),ll,ll(lc1),ll(li1),ll(lm1),ll(lp1),ll(lq1),ll(lr1),
     *  ll(ls1),aa(np1),nprof,aa,ifail)
      if(ifail.gt.0)return
c     write(nout,*)'steepest edge coefficients',(e(ij),ij=1,nm)
c     emax=0.D0
c     do i=1,nm
c       if(e(i).gt.0.D0)then
c         call eptsol(n,a,la,i,a,aa(ns1),aa(nt1),aa,aa(np1),
c    *      ll,ll(lc1),ll(li1),ll(lp1),ll(lq1))
c         ei=xlen(0.D0,aa(ns1),n)
c         ei=sqrt(scpr(0.D0,aa(ns1),aa(ns1),n))
c         emax=max(emax,abs(ei-e(i)))
c       endif
c     enddo
c     if(emax.ge.tol)
c    *  write(nout,*)'error in steepest edge coefficients =',emax
      return
      end

      subroutine refactor(n,nm,a,la,aa,ll,ifail)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),aa(*),ll(*)
      common/sparsec/ns,ns1,nt,nt1,nu,nu1,nx,nx1,np,np1,nprof,
     *  lc,lc1,li,li1,lm,lm1,lp,lp1,lq,lq1,lr,lr1,ls,ls1,lt,lt1
      common/factorc/m1,m2,mp,mq,lastr,irow
      common/noutc/nout
c     write(nout,*)'refactor'
      m=nm-n
      call re_order(n,nm,a,la(1),la(la(0)),ll,ll(lc1),ll(li1),
     *  ll(lm1),ll(lp1),ll(lq1),ll(lr1),ll(ls1),ll(lt1),aa(np1),
     *  nprof,ifail)
      if(ifail.ge.1)then
c       write(nout,*)'failure in re_order (2)'
        return
      endif
      call re_factor(n,nm,a,la,ll,ll(lc1),ll(li1),ll(lm1),
     *  ll(lp1),ll(lq1),ll(lr1),ll(ls1),ll(lt1),aa(np1),
     *  nprof,aa,ifail)
      if(ifail.eq.7)return
      call check_L(n,aa,ll(lp1),ifail)
      return
      end

      subroutine pivot(p,q,n,nm,a,la,e,aa,ll,ifail,info)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(0:*),e(*),aa(*),ll(*),info(*)
      common/noutc/nout
      common/iprintc/iprint
      common/sparsec/ns,ns1,nt,nt1,nu,nu1,nx,nx1,np,np1,nprof,
     *  lc,lc1,li,li1,lm,lm1,lp,lp1,lq,lq1,lr,lr1,ls,ls1,lt,lt1
      common/factorc/m1,m2,mp,mq,lastr,irow
      common/mxm1c/mxm1
      common/refactorc/nup,nfreq
      common/epsc/eps,tol,emin
c     write(nout,*)'pivot: p,q =',p,q
      ifail=0
      if(p.ne.mp)then
        call eptsol(n,a,la,p,a,aa(ns1),aa(nt1),aa,aa(np1),
     *    ll,ll(lc1),ll(li1),ll(lp1),ll(lq1))
        if(p.gt.n)then
          e(p)=xlen(0.D0,aa(ns1+m2),m1)
        else
          e(p)=xlen(1.D0,aa(ns1+m2),m1)
        endif
        epp=e(p)
        mp=p
      endif
      if(q.ne.mq)then
        call aqsol(n,a,la,q,a,aa(nt1),aa(nx1),aa,aa(np1),
     *    ll,ll(lc1),ll(li1),ll(lp1),ll(lq1))
        mq=q
      endif
c  update steepest edge coefficients
      tp=aa(nt+ll(li+p))
      if(tp.eq.0.D0)tp=eps
      ep=e(p)
      eq=2.D0/ep
c     do i=1,m2-1
c       aa(nu+i)=0.D0
c     enddo
c     do i=m2,n
      do i=1,n
        aa(nu+i)=eq*aa(ns+i)
      enddo
      call aqsol(n,a,la,-1,a,aa(nu1),aa(nx1),aa,aa(np1),
     *  ll,ll(lc1),ll(li1),ll(lp1),ll(lq1))
c     write(nout,*)'row perm',(ll(ij),ij=1,n)
c     write(nout,*)'column perm',(ll(lc+ij),ij=m2+1,n)
c     write(nout,*)'s =',(aa(ns+ij),ij=1,n)
c     write(nout,*)'t =',(aa(nt+ij),ij=1,n)
c     write(nout,*)'u =',(aa(nu+ij),ij=1,n)
      e(p)=0.D0
      eq=ep/tp
      do i=1,nm
        if(e(i).gt.0.D0)then
          j=ll(li+i)
          ei=e(i)
          wi=aa(nt+j)*eq
          awi=abs(wi)
          if(ei.ge.awi)then
            wi=wi/ei
            e(i)=max(emin,ei*sqrt(max(0.D0,1.D0+wi*(wi-aa(nu+j)/ei))))
          else
            wi=ei/wi
            e(i)=max(emin,awi*sqrt(max(0.D0,1.D0+wi*(wi-aa(nu+j)/ei))))
          endif
        endif
      enddo
      e(q)=max(emin,abs(eq))
      info(1)=info(1)+1
      if(nup.ge.nfreq)then
c     if(nup.ge.30)then
c  refactorize L
        ip=ll(li+p)
        if(p.gt.n)then
          m2=m2+1
          qq=ll(lc+m2)
          ll(lc+ip)=qq
          ll(li+qq)=ip
          ll(li+p)=0
        else
          ll(ip)=ll(m2)
          ll(li+ll(ip))=ip
          ll(m2)=p
          ll(li+p)=m2
        endif
        if(q.gt.n)then
          ll(lc+m2)=q
          ll(li+q)=m2
          m2=m2-1
        else
          iq=ll(li+q)
          ll(iq)=ll(m2)
          ll(li+ll(iq))=iq
          ll(m2)=q
          ll(li+q)=m2
        endif
        m1=n-m2
        call re_order(n,nm,a,la(1),la(la(0)),ll,ll(lc1),ll(li1),
     *    ll(lm1),ll(lp1),ll(lq1),ll(lr1),ll(ls1),ll(lt1),aa(np1),
     *    nprof,ifail)
        if(ifail.ge.1)then
c         write(nout,*)'failure in re_order (3)'
          return
        endif
        call re_factor(n,nm,a,la,ll,ll(lc1),ll(li1),
     *    ll(lm1),ll(lp1),ll(lq1),ll(lr1),ll(ls1),ll(lt1),aa(np1),
     *    nprof,aa,ifail)
      else
c  update L
        call update_L(p,q,n,nm,a,la,ll,ll(lc1),ll(li1),ll(lm1),ll(lp1),
     *    ll(lq1),ll(lr1),ll(ls1),aa(np1),nprof,aa,aa(ns1),ifail)
      endif
      if(ifail.eq.7)return
      mp=-1
      mq=-1
      call check_L(n,aa,ll(lp1),ifail)
c     write(nout,*)'steepest edge coefficients',(e(ij),ij=1,nm)
c     emax=0.D0
c     do i=1,nm
c       if(e(i).gt.0.D0)then
c         call eptsol(n,a,la,i,a,aa(ns1),aa(nt1),aa,aa(np1),
c    *      ll,ll(lc1),ll(li1),ll(lp1),ll(lq1))
c         ei=xlen(0.D0,aa(ns1),n)
c         ei=sqrt(scpr(0.D0,aa(ns1),aa(ns1),n))
c         emax=max(emax,abs(ei-e(i)))
c       endif
c     enddo
c     if(emax.ge.tol)
c    *  write(nout,*)'error in steepest edge coefficients =',emax
      return
      end

      subroutine fbsub(n,jmin,jmax,a,la,q,b,x,ls,aa,ll,save)
      implicit double precision (a-h,r-z), integer (i-q)
      logical save
      dimension a(*),la(*),b(*),x(*),ls(*),aa(*),ll(*)

c  solves a system  B.x=b

c  Parameter list
c  **************
c   n   number of variables (as for bqpd)
c   jmin,jmax  (see description of ls below)
c   a,la   specification of QP problem data (as for bqpd)
c   q   an integer which, if in the range 1:n+m, specifies that the rhs vector
c       b is to be column q of the matrix A of general constraint normals.
c       In this case the parameter b is not referenced by fbsub.
c       If q=0 then b is taken as the vector given in the parameter b.
c   b(n)  must be set to the r.h.s. vector b (but only if q=0)
c   x(n+m)  contains the required part of the solution x, set according to the
c       index number of that component (in the range 1:n for a simple bound and
c       n+1:n+m for a general constraint)
c   ls(*)  an index vector, listing the components of x that are required.
c       Only the absolute value of the elements of ls are used (this allows
c       the possibility of using of the contents of the ls parameter of bqpd).
c       Elements of x in the range abs(ls(j)), j=jmin:jmax are set by fbsub.
c       These contortions allow bqpd to be independent of the basis matrix code.
c   aa(*)  real storage used by the basis matrix code (supply the vector
c       ws(lu1) with ws as in the call of bqpd and lu1 as in common/bqpdc/...)
c   ll(*)  integer storage used by the basis matrix code (supply the vector
c       lws(ll1) with lws as in the call of bqpd and ll1 as in common/bqpdc/...)
c   save   indicates if fbsub is to save its copy of the solution for possible
c       future use. We suggest that the user only sets save = .false.

      common/noutc/nout
      common/sparsec/ns,ns1,nt,nt1,nu,nu1,nx,nx1,np,np1,nprof,
     *  lc,lc1,li,li1,lm,lm1,lp,lp1,lq,lq1,lr,lr1,ls_,ls1,lt,lt1
      common/factorc/m1,m2,mp,mq,lastr,irow
c     write(nout,*)'fbsub  q =',q
      if(save)then
        if(q.ne.mq)then
          call aqsol(n,a,la,q,b,aa(nt1),aa(nx1),aa,aa(np1),
     *      ll,ll(lc1),ll(li1),ll(lp1),ll(lq1))
          mq=q
        endif
        do j=jmin,jmax
          i=abs(ls(j))
          x(i)=aa(nt+ll(li+i))
        enddo
      else
        call aqsol(n,a,la,q,b,aa(nu1),aa(nx1),aa,aa(np1),
     *    ll,ll(lc1),ll(li1),ll(lp1),ll(lq1))
        do j=jmin,jmax
          i=abs(ls(j))
          x(i)=aa(nu+ll(li+i))
        enddo
c       print *,'x =',(x(i),i=1,18)
      endif
      return
      end

      subroutine ztg(n,k,rg,lv,aa,ll)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension rg(*),lv(*),aa(*),ll(*)
      common/sparsec/ns,ns1,nt,nt1,nu,nu1,nx,nx1,np,np1,nprof,
     *  lc,lc1,li,li1,lm,lm1,lp,lp1,lq,lq1,lr,lr1,ls_,ls1,lt,lt1
c     print *,'aa =',(aa(nu+i),i=1,18)
      do j=1,k
        rg(j)=aa(nu+ll(li+lv(j)))
      enddo
      return
      end

      subroutine tfbsub(n,a,la,p,b,x,aa,ll,ep,save)
      implicit double precision (a-h,r-z), integer (i-q)
      logical save
      dimension a(*),la(*),b(*),x(*),aa(*),ll(*)

c  solves a system  Bt.x=b

c  Parameter list
c  **************
c   n   number of variables (as for bqpd)
c   a,la   specification of QP problem data (as for bqpd)
c   p    an integer which, if in the range 1:n+m, specifies that the rhs vector
c        b is a unit vector appropriate to the position of p in the current
c        ordering. In this case b is not referenced by tfbsub.
c   b(n+m)  If p=0, this must be set to the r.h.s. vector b. Only the components
c        of b need be set, according to the index number of each component (in
c        the range 1:n for a simple bound and n+1:n+m for a general constraint)
c   x(n)  contains the solution x (in natural ordering)
c   aa(*)  real storage used by the basis matrix code (supply the vector
c       ws(lu1) with ws as in the call of bqpd and lu1 as in common/bqpdc/...)
c   ll(*)  integer storage used by the basis matrix code (supply the vector
c       lws(ll1) with lws as in the call of bqpd and ll1 as in common/bqpdc/...)
c   ep  if p.ne.0 and save is true, ep contains the l_2 length of x on exit
c   save  indicates if tfbsub is to save its copy of the solution for possible
c       future use. We suggest that the user only sets save = .false.

      common/noutc/nout
      common/sparsec/ns,ns1,nt,nt1,nu,nu1,nx,nx1,np,np1,nprof,
     *  lc,lc1,li,li1,lm,lm1,lp,lp1,lq,lq1,lr,lr1,ls,ls1,lt,lt1
      common/factorc/m1,m2,mp,mq,lastr,irow
c     write(nout,*)'tfbsub  p =',p
      if(save)then
        if(p.ne.mp)then
          call eptsol(n,a,la,p,b,aa(ns1),aa(nt1),aa,aa(np1),
     *      ll,ll(lc1),ll(li1),ll(lp1),ll(lq1))
          mp=p
        endif
        do i=1,n
          x(ll(i))=aa(ns+i)
        enddo
        if(p.gt.n)then
          ep=xlen(0.D0,aa(ns1+m2),m1)
        elseif(p.gt.0)then
          ep=xlen(1.D0,aa(ns1+m2),m1)
        endif
      else
        call eptsol(n,a,la,p,b,aa(nu1),aa(nt1),aa,aa(np1),
     *    ll,ll(lc1),ll(li1),ll(lp1),ll(lq1))
        do i=1,n
          x(ll(i))=aa(nu+i)
        enddo
      endif
c     write(nout,*)'x =',(x(i),i=1,n)
      return
      end

      subroutine newg
      common/factorc/m1,m2,mp,mq,lastr,irow
      mq=-1
      return
      end

c******** The following routines are internal to sparseL.f **************

      subroutine check_L(n,d,p,ifail)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension d(*),p(*)
      common/noutc/nout
      common/factorc/m1,nu,mp,mq,lastr,irow
      common/epsc/eps,tol,emin
c     write(nout,*)'check_L'
      ifail=1
c     dmin=1.D37
      do k=nu+1,n
c       dmin=min(dmin,abs(d(k)))
        if(abs(d(k)).le.tol)return
      enddo
c     write(nout,*)'dmin =',dmin
c     len=0
c     do i=1,n
c       len=len+p(i)
c     enddo
c     write(nout,*)m1*(m1+1)/2,len+m1
c     write(nout,*)'m1 =',m1,'   file length =',len,'   total =',len+m1
      ifail=0
      return
      end

      subroutine aqsol(n,a,la,q,b,tn,xn,d,ws,lr,lc,li,pp,qq)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),b(*),tn(*),xn(*),d(*),ws(*),
     *  lr(*),lc(*),li(*),pp(*),qq(*)
      common/noutc/nout
      common/factorc/m1,m2,mp,mq,lastr,irow
c     write(nout,*)'aqsol  q =',q
      if(q.gt.0)then
        do i=1,n
          tn(i)=0.D0
        enddo
        if(q.le.n)then
          tn(li(q))=1.D0
        else
          call iscatter(a,la,q-n,li,tn,n)
        endif
      elseif(q.eq.0)then
        do i=1,n
          tn(li(i))=b(i)
        enddo
      endif
c     write(nout,*)'tn =',(tn(i),i=1,n)
      do i=n,m2+1,-1
        ir=lr(i)
        pri=pp(ir)
        if(pri.eq.0)then
          xn(i)=tn(i)/d(i)
        else
          xn(i)=(scpr(tn(i),ws(qq(ir)+1),tn(i-pri),pri))/d(i)
        endif
        call isaipy(-xn(i),a,la,lc(i)-n,tn,n,lr,li)
      enddo
      do i=m2+1,n
        tn(i)=xn(i)
      enddo
c     write(nout,*)'tn =',(tn(i),i=1,n)
      return
      end

      subroutine eptsol(n,a,la,p,b,sn,tn,d,ws,lr,lc,li,pp,qq)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),b(*),sn(*),tn(*),d(*),ws(*),
     *  lr(*),lc(*),li(*),pp(*),qq(*)
      common/noutc/nout
      common/iprintc/iprint
      common/epsc/eps,tol,emin
      common/factorc/m1,m2,mp,mq,lastr,irow
c     write(nout,*)'eptsol  p =',p
      if(p.eq.0)then
        do i=1,m2
          sn(i)=b(lr(i))
        enddo
        do i=m2+1,n
          sn(i)=0.D0
        enddo
        do i=m2+1,n
          j=lc(i)
          sn(i)=-aiscpri(n,a,la,j-n,sn,-b(j),lr,li)/d(i)
          ir=lr(i)
          pri=pp(ir)
          if(pri.gt.0)call mysaxpy(sn(i),ws(qq(ir)+1),sn(i-pri),pri)
        enddo
      else
        do i=1,n
          sn(i)=0.D0
        enddo
        pr=li(p)
        if(p.le.n)then
          if(pr.gt.m2)goto1
          sn(pr)=1.D0
          do i=m2+1,n
            sn(i)=-aiscpri(n,a,la,lc(i)-n,sn,0.D0,lr,li)/d(i)
            ir=lr(i)
            pri=pp(ir)
            if(pri.gt.0)call mysaxpy(sn(i),ws(qq(ir)+1),sn(i-pri),pri)
          enddo
        else
          if(pr.le.m2)goto1
          do i=m2+1,n
            bi=0.D0
            if(i.eq.pr)bi=-1.D0
            sn(i)=-aiscpri(n,a,la,lc(i)-n,sn,bi,lr,li)/d(i)
            ir=lr(i)
            pri=pp(ir)
            if(pri.gt.0)call mysaxpy(sn(i),ws(qq(ir)+1),sn(i-pri),pri)
          enddo
        endif
      endif
c     write(nout,*)'sn =',(sn(i),i=1,n)
      return
    1 continue
      write(nout,*)'malfunction detected in eptsol: p =',p
      stop
      end

      subroutine order(n,nu,nc,la,lr,ls,li,p,q,r,ws,mxws,ifail)
      implicit integer (c-t)
      double precision ws
      dimension la(0:*),lr(*),ls(*),li(*),p(*),q(*),r(*),ws(*)
      common/noutc/nout
c     character star(1000,80)
c     write(nout,*)'order'
c  spk1 ordering on full matrix
      ifail=0
      if(nu.eq.n)return
c  set row and column counts and row-wise data structure
      nn=n-nu
      ii=mxws/nn
      do j=1,nn
        rowj=lr(j)
        p(rowj)=(j-1)*ii
        r(rowj)=0
      enddo
      do j=nn+1,n
        r(lr(j))=0
      enddo
    1 continue
      do i=nu+1,nc
        coli=abs(ls(i))
        li(coli)=0
        jp=la(0)+coli-n
        do j=la(jp),la(jp+1)-1
          rowj=la(j)
          if(li(rowj).le.nn)then
            li(coli)=li(coli)+1
            r(rowj)=r(rowj)+1
            ij=p(rowj)+r(rowj)
            if(ij.gt.mxws)then
              ij=mxws
              ifail=1
            endif
            ws(ij)=dble(coli)
          endif
        enddo
      enddo
c  check for no overlaps
      qrj=0
      do j=1,nn
        rowj=lr(j)
        if(p(rowj).lt.qrj)ifail=1
        qrj=p(rowj)+r(rowj)
        q(rowj)=qrj
        p(rowj)=p(rowj)+1
      enddo
      if(ifail.eq.1.or.qrj.gt.mxws)then
        qrj=0
        do j=1,nn
          rowj=lr(j)
          p(rowj)=qrj
          qrj=qrj+r(rowj)
          r(rowj)=0
        enddo
        if(qrj.gt.mxws)then
          write(nout,*)'not enough space for ws in order:  mxws =',mxws
          ifail=7
          return
        endif
        ifail=0
        goto1
      endif
      ifirstc=nu+1
      ifirstr=1
    2 continue
c  move zero-column-count columns to lhs and find minimum column count
      mcc=n
      do i=ifirstc,nc
        coli=abs(ls(i))
        if(li(coli).eq.0)then
          call iexch(ls(i),ls(ifirstc))
          li(coli)=ifirstr-1
          ifirstc=ifirstc+1
        else
          mcc=min(mcc,li(coli))
        endif
      enddo
c     write(nout,*)'ifirstc,ifirstr,mcc',ifirstc,ifirstr,mcc
c     write(nout,*)'lr =',(lr(j),j=1,n)
c     write(nout,*)'ls =',(ls(i),i=nu+1,nc)
c     write(nout,*)'row counts =',(r(lr(j)),j=1,n)
c     write(nout,*)'column counts =',(li(abs(ls(i))),i=nu+1,nc)
      if(ifirstc.gt.nc)goto4
c  apply tie-break rule
      tie=0
      do i=ifirstc,nc
        coli=abs(ls(i))
        if(li(coli).eq.mcc)then
          ti=0
          jp=la(0)+coli-n
          do j=la(jp),la(jp+1)-1
            rowj=la(j)
            if(li(rowj).ge.ifirstr)ti=ti+r(rowj)
          enddo
          if(ti.gt.tie)then
            tie=ti
            mccc=coli
          endif
        endif
      enddo
c     write(nout,*)'tie,mccc',tie,mccc
c  permute rows of m-c-c column to top and update column counts
      jp=la(0)+mccc-n
      do j=la(jp),la(jp+1)-1
        rowj=la(j)
        jr=li(rowj)
        if(jr.lt.ifirstr)goto3
        if(jr.gt.nn)goto3
        lr(jr)=lr(ifirstr)
        li(lr(jr))=jr
        lr(ifirstr)=rowj
        li(rowj)=ifirstr
        ifirstr=ifirstr+1
        do i=p(rowj),q(rowj)
          coli=int(ws(i))
          li(coli)=li(coli)-1
        enddo
    3   continue
      enddo
      goto2
    4 continue
c  print star diagram
c     if(nc-nu.gt.80.or.n.gt.1000)stop
c     write(nout,*)'spk1 ordering'
c     ij=li(abs(ls(nc)))
c     do i=1,ij
c       do j=1,nc-nu
c         star(i,j)=' '
c       enddo
c     enddo
c     do j=1,nc-nu
c       jp=la(0)+abs(ls(nu+j))-n
c       do i=la(jp),la(jp+1)-1
c         star(li(la(i)),j)='*'
c       enddo
c     enddo
c     do i=1,ij
c       write(nout,*)(star(i,j),j=1,nc-nu)
c     enddo
c     write(nout,*)'lr =',(lr(i),i=1,n)
c     write(nout,*)'ls =',(ls(i),i=nu+1,nc)
c     write(nout,*)'lower profile =',(li(abs(ls(i))),i=nu+1,nc)
      return
      end

      subroutine factor(n,nm,nu,nc,a,la,e,ls,sn,tn,un,xn,lr,lc,li,
     *  mao,p,q,r,s,ws,mxws,d,ifail)
      implicit double precision (a-h,r-z), integer (i-q)
      integer coli,r,s,rowi,rowp,tl,tu
      dimension a(*),la(0:*),e(*),ls(*),sn(*),tn(*),un(*),xn(*),
     *  lr(*),lc(*),li(*),mao(*),p(*),q(*),r(*),s(*),ws(*),d(*)
c     character star(1000,80)
      common/factorc/m1,m2,mp,mq,lastr,irow
      common/iprintc/iprint
      common/refactorc/nup,nfreq
      common/epsc/eps,tol,emin
      common/noutc/nout
      parameter (thresh=1.D-1)
c  factorize LPA=U when A is rectangular
c    p(row) stores the number of stored elements of a natural row
c    q(row) stores the base address in ws of a natural row
c    r(row) stores the previous row stored in ws (or 0 if the first row in ws)
c    s(row) stores the next row stored in ws (or 0 if the last row in ws)
c    li(n+*) stores the lower profile of the sparse matrix
c    irow stores the natural row number of the initial row stored in ws
c    lastr stores the natural row number of the previous row put into ws
c     write(nout,*)'factor'
      nup=0
      lastr=0
      irow=0
      do i=1,n
        p(i)=0
      enddo
      m1=0
      tl=1
      do ii=nu+1,nc
        coli=abs(ls(ii))
c       write(nout,*)'coli =',coli
        tu=li(coli)
        do i=1,n
          tn(i)=0.D0
        enddo
        call iscatter(a,la,coli-n,li,tn,n)
        do i=m1,1,-1
          rowi=lr(i)
          pri=p(rowi)
          if(pri.eq.0)then
            xn(i)=tn(i)/d(i)
          else
            xn(i)=(scpr(tn(i),ws(q(rowi)+1),tn(i-pri),pri))/d(i)
          endif
          call isaipy(-xn(i),a,la,lc(i)-n,tn,n,lr,li)
        enddo
        do i=1,m1
          tn(i)=xn(i)
        enddo
        m1p=m1+1
c       write(nout,*)'lr =',(lr(i),i=1,n)
c       write(nout,*)'tn =',(tn(i),i=1,tu)
c  threshold pivot selection
        call linf(tu-m1,tn(m1p),z,iz)
        if(z.le.tol)then
          li(coli)=0
          goto2
        endif
        zz=max(tol,z*thresh)
        do i=tl,tu
          q(lr(i))=m1p
        enddo
c       write(nout,*)'q =',(q(lr(i)),i=m1p,tu)
        iz=iz+m1
        if(iz.lt.tl)then
          z=0.D0
          qri=m1p
          do j=m1p,tu
            tnj=abs(tn(j))
            if(tnj.ge.zz)then
              qrj=q(lr(j))
              if(qrj.eq.qri)then
                if(tnj.gt.z)then
                  z=tnj
                  iz=j
                endif
              elseif(qrj.gt.qri)then
                z=tnj
                iz=j
                qri=qrj
              endif
            endif
          enddo
        endif
        tl=tu+1
c       write(nout,*)'zz,z,iz,m1,qri',zz,z,iz,m1,qri
        if(iz.gt.m1p)then
          call rexch(tn(m1p),tn(iz))
          call iexch(lr(m1p),lr(iz))
          li(lr(m1p))=m1p
          li(lr(iz))=iz
        endif
        rowp=lr(m1p)
c  reset q values
        qrp=q(rowp)
        do i=m1p+1,tu
          if(abs(tn(i)).gt.tol)then
            rowi=lr(i)
            if(qrp.lt.q(rowi))q(rowi)=qrp
          endif
        enddo
        tnp=tn(m1p)
        do i=1,n
          sn(i)=0.D0
        enddo
        sn(m1p)=1.D0
        do i=1,m1
          sn(i)=-aiscpri(n,a,la,lc(i)-n,sn,0.D0,lr,li)/d(i)
          rowi=lr(i)
          pri=p(rowi)
          if(pri.gt.0)call mysaxpy(sn(i),ws(q(rowi)+1),sn(i-pri),pri)
        enddo
c       write(nout,*)'sn =',(sn(i),i=1,m1)
c  update steepest edge coefficients
        ep=e(rowp)
        e(rowp)=0.D0
        eq=2.D0/ep
        do i=1,n
          un(i)=eq*sn(i)
        enddo
        do i=m1,1,-1
          rowi=lr(i)
          pri=p(rowi)
          if(pri.eq.0)then
            xn(i)=un(i)/d(i)
          else
            xn(i)=(scpr(un(i),ws(q(rowi)+1),un(i-pri),pri))/d(i)
          endif
          call isaipy(-xn(i),a,la,lc(i)-n,un,n,lr,li)
        enddo
        do i=1,m1
          un(i)=xn(i)
        enddo
c       write(nout,*)'un =',(un(i),i=1,n)
        eq=ep/tnp
        do i=1,nm
          if(e(i).gt.0.D0)then
            j=li(i)
            ei=e(i)
            wi=tn(j)*eq
            awi=abs(wi)
            if(ei.ge.awi)then
              wi=wi/ei
              e(i)=max(emin,ei*sqrt(max(0.D0,1.D0+wi*(wi-un(j)/ei))))
            else
              wi=ei/wi
              e(i)=max(emin,awi*sqrt(max(0.D0,1.D0+wi*(wi-un(j)/ei))))
            endif
          endif
        enddo
        e(coli)=max(emin,abs(eq))
        do j=qrp,m1
          if(abs(sn(j)).gt.tol)goto1
        enddo
        j=m1p
    1   continue
        pri=m1p-j
        if(pri.gt.0)then
          call newslot(rowp,pri,lastr,irow,p,q,r,s,ws,mxws,i,ifail)
          if(ifail.gt.0)return
          p(rowp)=pri
          i=q(rowp)
          do j=j,m1
            i=i+1
            ws(i)=sn(j)
          enddo
        endif
        m1=m1p
        ls(m1)=ls(ii)
        lc(m1)=coli
        li(coli)=m1
        d(m1)=tnp
    2   continue
      enddo
c  complete ls and reorder lr, lc and d
      do i=m1+1,n
        ls(i)=lr(i)
      enddo
      j=n
      do i=1,nm
        if(e(i).eq.0.D0)then
          j=j+1
          ls(j)=i
        endif
      enddo
      m2=n-m1
      do i=n,m2+1,-1
        lc(i)=lc(i-m2)
        li(lc(i))=i
        lr(i)=lr(i-m2)
        li(lr(i))=i
        d(i)=d(i-m2)
      enddo
      do i=1,m2
        lr(i)=ls(m1+i)
        li(lr(i))=i
      enddo
c  reset mao
      ilast=n
      ii=ilast
      do i=ilast,m2+1,-1
        mao(i)=ilast
        ii=min(ii,i-p(lr(i)))
        if(ii.eq.i)ilast=i-1
      enddo
c     write(nout,*)'PAQ factors:  m1 =',m1
c     write(nout,*)'d =',(d(ij),ij=m2+1,n)
c     do j=m2+1,n
c       rowp=lr(j)
c       if(p(rowp).ne.0)then
c         write(nout,*)'L(',rowp,')',
c    *      (ws(k),k=q(rowp)+1,q(rowp)+p(rowp))
c       endif
c     enddo
c  print star diagram
c     write(nout,*)'factored ordering:  m1 =',m1
c     if(m1.gt.80.or.n.gt.1000)stop
c     do i=1,n
c       do j=1,m1
c         star(i,j)=' '
c       enddo
c     enddo
c     do j=1,m1
c       jp=la(0)+lc(m2+j)-n
c       do i=la(jp),la(jp+1)-1
c         star(li(la(i)),j)='*'
c       enddo
c     enddo
c     do i=m2+1,n
c       write(nout,*)(star(i,j),j=1,m1)
c     enddo
c     write(nout,*)'ls =',(ls(j),j=1,n)
c     write(nout,*)'s.e. coeffs =',(e(i),i=1,nm)
c     write(nout,*)'lr =',(lr(j),j=1,n)
c     write(nout,*)'lc =',(lc(j),j=m2+1,n)
c     write(nout,*)'mao =',(mao(j),j=m2+1,n)
c     call checkout(n,a,la,lr,lc,li,p,q,r,s,ws,mxws,d)
      return
      end

      subroutine re_order(n,nm,a,la,point,lr,lc,li,mao,p,q,r,s,
     *  t,ws,mxws,ifail)
      implicit double precision (a-h,u-z), integer (i-t)
      dimension a(*),la(*),point(0:*),lr(*),lc(*),li(*),mao(*),
     *  p(*),q(*),r(*),s(*),t(*),ws(*)
      common/factorc/m1,nu,mp,mq,lastr,irow
      common/noutc/nout
      logical backtrack
c     character star(1000,80)
c  print star diagram
c     if(n-nu.gt.80.or.n.gt.1000)stop
c     write(nout,*)'initial ordering'
c     do i=1,n
c       do j=1,n-nu
c         star(i,j)=' '
c       enddo
c     enddo
c     do j=1,n-nu
c       ilp=lc(nu+j)-n
c       do i=point(ilp),point(ilp+1)-1
c         star(li(la(i)),j)='*'
c       enddo
c     enddo
c     do i=nu+1,n
c       write(nout,*)(star(i,j),j=1,n-nu)
c     enddo
c     write(nout,*)'re_order'
      if(nu.eq.n)then
        ifail=0
        return
      endif
      m=nm-n
c  transversal search
      do iq=nu+1,n
        backtrack=.false.
        istack=nu
        inode=iq
        nodec=lc(inode)
        nodec_n=nodec-n
        lap=point(nodec_n+1)-point(nodec_n)
c       write(nout,*)'column node =',nodec,'  look-ahead rows =',
c    *    (la(j),j=point(nodec_n),point(nodec_n)+lap-1)
c  look-ahead loop
    1   continue
          lap=lap-1
          nextr=la(point(nodec_n)+lap)
          inext=li(nextr)
          if(inext.ge.iq)goto4
          if(lap.gt.0)goto1
          li(nodec)=0
    2   continue
c  reassignment depth first search
        t(inode)=point(nodec_n+1)-point(nodec_n)
c       write(nout,*)'column node =',nodec,'  unfathomed rows =',
c    *    (la(j),j=point(nodec_n),point(nodec_n)+t(inode)-1)
    3   continue
c  examine successor nodes
        if(t(inode).eq.0)then
          if(istack.eq.nu)then
            ifail=1
c           ifail=iq
c           write(nout,*)'exit: ifail =',iq
            return
          endif
          istack=istack-1
          backtrack=.true.
          if(istack.eq.nu)then
            inode=iq
          else
            inode=mao(istack)
          endif
c         write(nout,*)'backtrack to node at address =',inode
          nodec=lc(inode)
          nodec_n=nodec-n
c         write(nout,*)'column node =',nodec,'  unfathomed rows =',
c    *      (la(j),j=point(nodec_n),point(nodec_n)+t(inode)-1)
          goto3
        endif
        t(inode)=t(inode)-1
        nextr=la(point(nodec_n)+t(inode))
        inext=li(nextr)
        if(inext.le.nu)goto3
        if(t(inext).ge.0)goto3
c  extend depth first search
c       write(nout,*)'nextr,inext',nextr,inext
        inode=inext
c       write(nout,*)'put node address on stack'
        istack=istack+1
        mao(istack)=inode
c       write(nout,*)'stack =',(mao(j),j=nu+1,istack)
        nodec=lc(inode)
        nodec_n=nodec-n
        lap=li(nodec)
        if(lap.eq.0)goto2
c       write(nout,*)'column node =',nodec,'  look-ahead rows =',
c    *    (la(j),j=point(nodec_n),point(nodec_n)+lap-1)
        goto1
    4   continue
c       write(nout,*)'new assignment found in row',nextr
c       write(nout,*)'istack,inext,nextr',istack,inext,nextr
c       if(istack.gt.nu)write(nout,*)'stack =',(mao(j),j=nu+1,istack)
        li(nodec)=lap
c  perform row permutation
        lr(inext)=lr(iq)
        li(lr(inext))=inext
        inode=iq
        do i=nu+1,istack
          inext=mao(i)
          lr(inode)=lr(inext)
          li(lr(inode))=inode
          inode=inext
        enddo
        lr(inode)=nextr
        li(nextr)=inode
c       write(nout,*)'lr =',(lr(j),j=nu+1,n)
c       write(nout,*)'look-ahead lengths =',(li(lc(j)),j=nu+1,iq)
        t(iq)=-1
        if(backtrack.or.istack.gt.nu+1)then
          do i=nu+1,iq-1
            t(i)=-1
          enddo
        endif
        do i=1,n
          if(li(i).gt.n)then
            write(nout,*)'iq =',iq
            stop
          endif
        enddo
      enddo
c     write(nout,*)'transversal found'
c     write(nout,*)'lr =',(lr(j),j=1,n)
c     write(nout,*)'lc =',(lc(j),j=nu+1,n)
c  print star diagram
c     if(n-nu.gt.80.or.n.gt.1000)stop
c     write(nout,*)'transversal ordering'
c     do i=1,n
c       do j=1,n-nu
c         star(i,j)=' '
c       enddo
c     enddo
c     do j=1,n-nu
c       ilp=lc(nu+j)-n
c       do i=point(ilp),point(ilp+1)-1
c         star(li(la(i)),j)='*'
c       enddo
c     enddo
c     do i=nu+1,n
c       write(nout,*)(star(i,j),j=1,n-nu)
c     enddo

c  tarjan ordering
      do i=1,n
        q(i)=0
        r(i)=0
      enddo
c  reset li and pair off columns with rows
      do i=nu+1,n
        nodec=lc(i)
        li(nodec)=i
        t(lr(i))=nodec
        s(i)=0
      enddo
      do i=nu+1,n
        noder=lr(i)
        nodec=t(noder)
        lc(noder)=point(nodec-n+1)-point(nodec-n)
        li(nodec)=-1
      enddo
      ifath=nu
      istack=n+1
c  tarjan loop
   10 continue
        istack=istack-1
        inode=istack
        noder=lr(inode)
        if(lc(noder).eq.0)then
          write(nout,*)'malfunction: zero length'
          stop
        endif
        nodec=t(noder)
   11   continue
        li(nodec)=lc(noder)
        mao(inode)=istack
c       write(nout,*)'put new node',noder,' on stack'
c       write(nout,*)'active part of lr =',(lr(j),j=ifath+1,n)
c       write(nout,*)'ifath,istack =',ifath,istack
c       write(nout,*)'column node =',nodec,'  unfathomed rows =',
c    *    (la(j),j=point(nodec-n),point(nodec-n)+li(nodec)-1)
   12   continue
          if(li(nodec).eq.0)then
c           write(nout,*)'backtrack to previous nodes'
   13       continue
              if(inode.eq.n)goto14
              inext=inode+1
              nextr=lr(inext)
              if(mao(inode).lt.mao(inext))goto14
              inode=inext
              noder=nextr
              nodec=t(noder)
              if(li(nodec).eq.0)goto13
c           write(nout,*)'stack =',(lr(j),j=istack,n)
c           write(nout,*)'lengths =',(li(t(lr(j))),j=istack,n)
c           write(nout,*)'column node =',nodec,'  unfathomed rows =',
c    *        (la(j),j=point(nodec-n),point(nodec-n)+li(nodec)-1)
            goto12
          endif
c  examine successors of current node
          li(nodec)=li(nodec)-1
          nextr=la(point(nodec-n)+li(nodec))
          inext=li(nextr)
          if(inext.le.ifath)goto12
          q(nextr)=q(nextr)+1
          nextc=t(nextr)
c         write(nout,*)'nextc,nextr,inext',nextc,nextr,inext
          if(li(nextc).ge.0)then
            mx=mao(inext)
            if(mao(inode).ge.mx)goto12
            do j=istack,n
              if(mao(j).eq.mx)goto12
              mao(j)=mx
            enddo
            write(nout,*)'malfunction'
            stop
          endif
          nodec=nextc
          noder=nextr
          istack=istack-1
          inode=istack
          lr(inext)=lr(inode)
          li(lr(inext))=inext
          lr(inode)=noder
          li(noder)=inode
          goto11
   14   continue
c       write(nout,*)'strong component identified'
c       write(nout,*)'active part of lr =',(lr(j),j=ifath+1,n)
c       write(nout,*)'ifath,istack,inode =',ifath,istack,inode,n
c  shift forward strong component
        inext=istack-1
        ir=inode-inext
        do j=istack,inode
          mao(j)=lr(j)
        enddo
        do j=inext+ir,ifath+1+ir,-1
          lr(j)=lr(j-ir)
          li(lr(j))=j
        enddo
        mx=ifath+ir
        iq=inext-ifath
        ifath=ifath+1
        do j=ifath,mx
          lr(j)=mao(j+iq)
          li(lr(j))=j
          mao(j)=mx
        enddo
        istack=inode+1
        ifath=mx
c       write(nout,*)'active part of lr =',(lr(j),j=ifath+1,n)
c       write(nout,*)'ifath,istack =',ifath,istack
        if(istack.le.n)then
          inode=istack
          noder=lr(inode)
          nodec=t(noder)
          nodec_n=nodec-n
c         write(nout,*)'column node =',nodec,'  unfathomed rows =',
c    *      (la(j),j=point(nodec-n),point(nodec-n)+li(nodec)-1)
          goto12
        endif
      if(ifath.lt.n)goto10
c  end of tarjan process
c  reset lc and li
      do i=nu+1,n
        lc(i)=t(lr(i))
        li(lc(i))=i
      enddo
c     write(nout,*)'mao =',(mao(j),j=nu+1,n)
c     write(nout,*)'q =',(q(j),j=1,n)
c     write(nout,*)'lr =',(lr(j),j=1,n)
c     write(nout,*)'lc =',(lc(j),j=nu+1,n)
c     write(nout,*)'li =',(li(j),j=1,n+m)
c  print star diagram
c     if(n-nu.gt.80.or.n.gt.1000)stop
c     write(nout,*)'tarjan ordering'
c     do i=1,n
c       do j=1,n-nu
c         star(i,j)=' '
c       enddo
c     enddo
c     do j=1,n-nu
c       ilp=lc(nu+j)-n
c       do i=point(ilp),point(ilp+1)-1
c         star(li(la(i)),j)='*'
c       enddo
c     enddo
c     do i=nu+1,n
c       write(nout,*)(star(i,j),j=1,n-nu)
c     enddo
c  set up pointers for row-wise sparse structure
      p(1)=1
      do i=1,n-1
        p(i+1)=p(i)+q(i)
        q(i)=p(i)-1
      enddo
      if(p(n)+q(n).gt.mxws)then
        ifail=7
        return
      endif
      q(n)=p(n)-1
      i=nu+1
   20 continue
      if(i.eq.mao(i))then
        t(i)=i
      else
c  spk1 ordering on tarjan block
c  set row and column counts
        do inode=i,mao(i)
          nodec=lc(inode)
          do j=point(nodec-n),point(nodec-n+1)-1
            noder=la(j)
            if(li(noder).ge.i)then
              q(noder)=q(noder)+1
              ws(q(noder))=dble(nodec)
              s(inode)=s(inode)+1
            endif
          enddo
        enddo
c       print *,'r-c counts: i =',i,'   mao(i) =',mao(i)
c       print *,'q =',(q(j),j=i,mao(i))
c       print *,'s =',(s(j),j=i,mao(i))
c  find minimum-column-count column
        mcc=n
        do inode=i,mao(i)
          noder=lr(inode)
          r(noder)=q(noder)-p(noder)+1
          mcc=min(mcc,s(inode))
        enddo
c     write(nout,*)'i,mao(i),mcc',i,mao(i),mcc
c     write(nout,*)'p =',(p(lr(j)),j=i,mao(i))
c     write(nout,*)'q =',(q(lr(j)),j=i,mao(i))
c     write(nout,*)'r =',(r(lr(j)),j=i,mao(i))
c     write(nout,*)'s =',(s(j),j=i,mao(i))
c  check for fully dense block
        if(mcc.gt.mao(i)-i)then
          do inode=i,mao(i)
            t(inode)=mao(i)
          enddo
          goto22
        endif
c  determine spk1 ordering
        ifirstr=i
        ifirstc=i
   21   continue
c  apply tie-break rule
        tie=0
        do inode=ifirstc,mao(i)
          if(s(inode).eq.mcc)then
            nodec=lc(inode)-n
            ti=0
            do j=point(nodec),point(nodec+1)-1
              noder=la(j)
              if(li(noder).ge.ifirstr)ti=ti+r(noder)
            enddo
            if(ti.gt.tie)then
              tie=ti
              mccc=nodec
            endif
          endif
        enddo
c       write(nout,*)'tie,mccc',tie,mccc+n
c  permute rows of m-c-c column to top and update column counts
        do j=point(mccc),point(mccc+1)-1
          noder=la(j)
          ir=li(noder)
          if(ir.ge.ifirstr)then
            lr(ir)=lr(ifirstr)
            li(lr(ir))=ir
            lr(ifirstr)=noder
            li(noder)=ifirstr
            ifirstr=ifirstr+1
            do ir=p(noder),q(noder)
              inode=li(int(ws(ir)))
              s(inode)=s(inode)-1
            enddo
          endif
        enddo
c       write(nout,*)'s =',(s(ij),ij=i,mao(i))
c       write(nout,*)'lr =',(lr(ij),ij=i,mao(i))
c  move zero-column-count columns to lhs and find minimum column count
        mcc=n
        do inode=ifirstc,mao(i)
          if(s(inode).eq.0)then
            nodec=lc(inode)
            lc(inode)=lc(ifirstc)
            li(lc(inode))=inode
            lc(ifirstc)=nodec
            li(nodec)=ifirstc
            s(inode)=s(ifirstc)
            t(ifirstc)=ifirstr-1
            ifirstc=ifirstc+1
          else
            mcc=min(mcc,s(inode))
          endif
        enddo
c       write(nout,*)'lc =',(lc(ij),ij=i,mao(i))
c       write(nout,*)'ifirstc,mcc',ifirstc,mcc
        if(ifirstc.lt.mao(i))goto21
      endif
   22 continue
      i=mao(i)+1
      if(i.le.n)goto20
c  print star diagram
c     if(n-nu.gt.80.or.n.gt.1000)stop
c     write(nout,*)'tarjan + spk1 ordering'
c     do i=1,n
c       do j=1,n-nu
c         star(i,j)=' '
c       enddo
c     enddo
c     do j=1,n-nu
c       ilp=lc(nu+j)-n
c       do i=point(ilp),point(ilp+1)-1
c         star(li(la(i)),j)='*'
c       enddo
c     enddo
c     do i=nu+1,n
c       write(nout,*)(star(i,j),j=1,n-nu)
c     enddo
c     write(nout,*)'lr =',(lr(j),j=nu+1,n)
c     write(nout,*)'lc =',(lc(j),j=nu+1,n)
c     write(nout,*)'lower profile =',(t(j),j=nu+1,n)
      ifail=0
      return
      end

      subroutine re_factor(n,nm,a,la,lr,lc,li,mao,p,q,r,s,
     *  t,ws,mxws,d,ifail)
      implicit double precision (a-h,u-z), integer (i-t)
      dimension a(*),la(0:*),lr(*),lc(*),li(*),mao(*),
     *  p(*),q(*),r(*),s(*),t(*),d(*),ws(*)
c     character star(1000,80)
      common/factorc/m1,nu,mp,mq,lastr,irow
      common/iprintc/iprint
      common/refactorc/nup,nfreq
      common/epsc/eps,tol,emin
      common/noutc/nout
      double precision thresh,tol
      parameter (thresh=1.D-1)
c  factorize LPA=U
c    p(row) stores the number of stored elements of a natural row
c    q(row) stores the base address in ws of a natural row
c    r(row) stores the previous row stored in ws (or 0 if the first row in ws)
c    s(row) stores the next row stored in ws (or 0 if the last row in ws)
c    t(*) stores the lower profile of the sparse matrix
c    irow stores the natural row number of the initial row stored in ws
c    lastr stores the natural row number of the previous row put into ws
c     write(nout,*)'re_factor'
      nup=0
      m=nm-n
      lastr=0
      irow=0
      do i=1,n
        p(i)=0
      enddo
      if(m1.eq.0)return
      i=nu+1
    1 continue
      if(i.eq.mao(i))then
        d(i)=aij(lr(i),lc(i)-n,a,la)
        if(d(i).eq.0.D0)d(i)=eps
c       write(nout,*)'row,col,d(i) =',lr(i),lc(i),d(i)
      else
c       write(nout,*)'lc =',(lc(j),j=i,mao(i))
        do inode=i,mao(i)-1
          nodec=lc(inode)-n
          im=inode-1
c  form L.a_q
          z=0.
c         write(nout,*)'inode,t(inode)',inode,t(inode)
          do j=inode,t(inode)
            rowj=lr(j)
            prj=p(rowj)
            if(prj.gt.0)then
              d(j)=aiscpri2(n,a,la,rowj,nodec,ws(q(rowj)+1),1.D0,im,
     *          prj,li)
            else
              d(j)=aij(rowj,nodec,a,la)
            endif
            z=max(z,abs(d(j)))
          enddo
c         write(nout,*)'d =',(d(ij),ij=inode,t(inode))
c  threshold pivot selection
          zz=z*thresh
          z=0.D0
          pri=n
          do j=inode,t(inode)
            dj=abs(d(j))
            if(dj.ge.zz)then
              prj=p(lr(j))
              if(prj.eq.pri)then
                if(dj.gt.z)then
                  z=dj
                  iz=j
                endif
              elseif(prj.lt.pri)then
                z=dj
                iz=j
                pri=prj
              endif
            endif
          enddo
c       write(nout,*)'zz,z,iz,pri',zz,z,iz,pri
          if(iz.gt.inode)then
c  pivot interchange
            call rexch(d(inode),d(iz))
            call iexch(lr(inode),lr(iz))
            li(lr(iz))=iz
            li(lr(inode))=inode
          endif
          if(d(inode).eq.0.D0)d(inode)=eps
c  update L
          qri=q(lr(inode))
          zz=-d(inode)
          do j=inode+1,t(inode)
            z=d(j)/zz
            rowj=lr(j)
            prj=p(rowj)
            qrj=q(rowj)
c  find space available in-situ in ws
            if(prj.eq.0)then
              len=0
            elseif(s(rowj).eq.0)then
              len=mxws-qrj
            else
              len=q(s(rowj))-qrj
            endif
            if(abs(z).le.tol)then
c  special case of a zero multiplier
              if(prj.eq.0)goto2
              len_=prj+1
              if(len_.gt.len)then
                call newslot(rowj,len_,lastr,irow,p,q,r,s,ws,mxws,qrj,
     *            ifail)
                if(ifail.gt.0)return
                qrj_=q(rowj)
                do k=1,prj
                  ws(qrj_+k)=ws(qrj+k)
                enddo
                ws(qrj_+len_)=z
              else
                ws(qrj+len_)=z
              endif
              p(rowj)=len_
              goto2
            endif
            len_=max(pri,prj)+1
            if(len_.gt.len.or.pri.gt.prj)then
c  create a new slot and use saxpyz ...
              call newslot(rowj,len_,lastr,irow,p,q,r,s,ws,mxws,qrj,
     *          ifail)
              if(ifail.gt.0)return
              qrj_=q(rowj)
              len=prj-pri
              if(len.ge.0)then
                do k=1,len
                  ws(qrj_+k)=ws(qrj+k)
                enddo
                len=len+1
                call saxpyz(z,ws(qri+1),ws(qrj+len),ws(qrj_+len),
     *            len_-len)
              else
                len=-len
                do k=1,len
                  ws(qrj_+k)=z*ws(qri+k)
                enddo
                len=len+1
                call saxpyz(z,ws(qri+len),ws(qrj+1),ws(qrj_+len),
     *            len_-len)
              endif
              ws(qrj_+len_)=z
            else
c  ... else saxpy in-situ
              if(pri.gt.0)
     *          call mysaxpy(z,ws(qri+1),ws(qrj+prj-pri+1),pri)
              ws(qrj+len_)=z
            endif
            p(rowj)=len_
c           do rj=1,n
c             if(p(rj).ne.0)then
c               write(nout,*)'storage for row',rj,'  p,q,r,s =',
c    *            p(rj),q(rj),r(rj),s(rj)
c             endif
c           enddo
    2       continue
          enddo
c         write(nout,*)'lr =',(lr(j),j=i,mao(i))
c         do j=i,mao(i)
c           rowj=lr(j)
c           if(p(rowj).ne.0)then
c             write(nout,*)'L(',rowj,')',
c    *          (ws(k),k=q(rowj)+1,q(rowj)+p(rowj))
c           endif
c         enddo
        enddo
        inode=mao(i)
        noder=lr(inode)
        pri=p(noder)
        if(pri.gt.0)then
         d(inode)=aiscpri2(n,a,la,noder,lc(inode)-n,ws(q(noder)+1),
     *     1.D0,inode-1,pri,li)
        else
          d(inode)=aij(noder,lc(inode)-n,a,la)
        endif
        if(d(inode).eq.0.D0)d(inode)=eps
      endif
      i=mao(i)+1
      if(i.le.n)goto1
c     write(nout,*)'PAQ factors:  nu =',nu
c     write(nout,*)'column perm =',(lc(j),j=nu+1,n)
c     write(nout,*)'row perm =',(lr(j),j=nu+1,n)
c     write(nout,*)'d =',(d(ij),ij=nu+1,n)
c     do j=nu+1,n
c       rowj=lr(j)
c       if(p(rowj).ne.0)then
c         write(nout,*)'L(',rowj,')',
c    *      (ws(k),k=q(rowj)+1,q(rowj)+p(rowj))
c       endif
c     enddo
c     call checkout(n,a,la,lr,lc,li,p,q,r,s,ws,mxws,d)
c  print star diagram
c     if(m1.gt.80.or.n.gt.1000)stop
c     write(nout,*)'factored tarjan + spk1 ordering:  nu =',nu
c     do i=1,n
c       do j=1,m1
c         star(i,j)=' '
c       enddo
c     enddo
c     do j=1,m1
c       jp=la(0)+lc(nu+j)-n
c       do i=la(jp),la(jp+1)-1
c         star(li(la(i)),j)='*'
c       enddo
c     enddo
c     do i=nu+1,n
c       write(nout,*)(star(i,j),j=1,m1)
c     enddo
c     write(nout,*)'lr =',(lr(j),j=nu+1,n)
c     write(nout,*)'lc =',(lc(j),j=nu+1,n)
      mp=-1
      mq=-1
      ifail=0
      return
      end

      function aiscpri2(n,a,la,rowi,coli,ws,di,im,pri,li)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),ws(*),li(*)
      integer rowi,coli,rowj,pri
      aiscpri2=0.D0
      jp=la(0)+coli
      do j=la(jp),la(jp+1)-1
        rowj=la(j)
        if(rowj.eq.rowi)then
          aiscpri2=aiscpri2+di*a(j)
        else
          ir=li(rowj)-im
          if(ir.gt.0)goto1
          ir=ir+pri
          if(ir.gt.0)aiscpri2=aiscpri2+ws(ir)*a(j)
        endif
    1   continue
      enddo
      return
      end

      subroutine update_L(pp,qq,n,nm,a,la,lr,lc,li,mao,p,q,r,s,
     *  ws,mxws,d,sn,ifail)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(0:*),lr(*),lc(*),li(*),mao(*),
     *  p(*),q(*),r(*),s(*),ws(*),d(*),sn(*)
c     character star(1000,80)
      double precision l11,l21
      integer r,s,rowim,rowi,rowj,rrj
      common/factorc/m1,nu,mp,mq,lastr,irow
      common/refactorc/nup,nfreq
      common/iprintc/iprint
      common/epsc/eps,tol,emin
      common/noutc/nout
      parameter (thresh=1.D-1,growth=1.D1)
c     write(nout,*)'update_L:  p,q =',pp,qq
      nup=nup+1
      if(qq.gt.n)then
        ilast=nu
        jp=la(0)+qq-n
        do j=la(jp),la(jp+1)-1
          ip=li(la(j))
          if(ip.gt.nu)ilast=max(ilast,mao(ip))
        enddo
        qqq=qq
      else
c  row flma procedure to remove row qq (includes qq amongst the unit vectors)
        iq=li(qq)
        if(iq.le.nu)goto99
        ilast=mao(iq)
        l11=1.D0
        u11=d(iq)
        ss=-sn(iq)
        nu=nu+1
        do i=iq,nu+1,-1
          lr(i)=lr(i-1)
          li(lr(i))=i
          sn(i)=sn(i-1)
          d(i)=d(i-1)
        enddo
        lr(nu)=qq
        li(qq)=nu
c  update mao
        do j=iq-1,nu,-1
          if(mao(j).lt.ilast)goto5
        enddo
        j=nu-1
    5   continue
        do j=j,nu,-1
          mao(j+1)=mao(j)+1
        enddo
        prq=p(qq)
        if(prq.gt.0)qrq=q(qq)
        do i=iq+1,ilast
          im=i-1
          rowi=lr(i)
          pri=p(rowi)
          u22=d(i)
          if(prq.gt.0)then
            u12=aiscpri2(n,a,la,qq,lc(i)-n,ws(qrq+1),l11,im,prq,li)
          else
            u12=l11*aij(qq,lc(i)-n,a,la)
          endif
          if(abs(u12).le.tol)u12=0.D0
          if(pri.gt.0)then
            qri=q(rowi)
            is=im-iq
            ii=pri-is
            if(ii.le.0)then
              l21=0.
            else
              l21=ws(qri+ii)
              if(abs(l21).le.tol)l21=0.D0
              if(ii.eq.1)then
                call trim_(rowi,pri,qri,q,ws)
                if(pri.eq.0)call erase(rowi,lastr,irow,r,s)
                if(s(rowi).eq.0)then
                  qr_=mxws
                else
                  qr_=q(s(rowi))
                endif
                if(qri+pri.ge.qr_)then
                  call r_shift(ws(qri),pri,1)
                  qri=qri-1
                  q(rowi)=qri
                endif
              else
                pri=pri-1
                call r_shift(ws(qri+ii),is,1)
              endif
            endif
            p(rowi)=pri
          else
            l21=0.D0
          endif
          rr=-l21/l11
          del=rr*u12+u22
          test=abs(rr)*max(abs(u11),abs(u22))
c         write(nout,*)'l11,l21,u11,u12,u22,del,test',
c    *      l11,l21,u11,u12,u22,del,test
          is=pri-prq
          if(is.lt.0)test=test*growth
          if(u12.eq.0.D0.and.is.gt.0)test=test*thresh
c           write(nout,*)'rowi,pri,qri =',rowi,pri,qri
c           write(nout,*)'rowq,prq,qrq =',qq,prq,qrq
c           write(nout,*)'j,p(j),q(j),r(j),s(j)   irow =',irow
c           do j=1,n
c             if(p(j).ne.0)write(nout,*)j,p(j),q(j),r(j),s(j)
c           enddo
c           write(nout,*)'rowq =',(ws(qrq+ij),ij=1,prq)
c           write(nout,*)'rowi =',(ws(qri+ij),ij=1,pri)
          if(abs(del).le.test)then
c  no-perm operation for row flma
c           write(nout,*)'no-perm operation for row flma'
            if(is.gt.0)then
              pr_=prq
              prq=pri+1
              call newslot(qq,prq,lastr,irow,p,q,r,s,ws,mxws,qr_,ifail)
              if(ifail.gt.0)return
              qrq=q(qq)
              qri=q(rowi)
              call r_shift(ws(qrq+1),pri,qri-qrq)
              call mysaxpy(rr,ws(qr_+1),ws(qri+is+1),pr_)
            else
              if(prq.eq.0)then
                call erase(rowi,lastr,irow,r,s)
                p(rowi)=0
                call newslot(qq,1,lastr,irow,p,q,r,s,ws,mxws,qr_,ifail)
                if(ifail.gt.0)return
                prq=1
                qrq=q(qq)
              else
                is=-is
                do j=1,is
                  ws(qrq+j)=rr*ws(qrq+j)
                enddo
                if(pri.gt.0)then
                  call saxpyx(rr,ws(qrq+is+1),ws(qri+1),pri)
                else
                  call newslot(rowi,1,lastr,irow,p,q,r,s,ws,mxws,qr_,
     *              ifail)
                  if(ifail.gt.0)return
                  qri=q(rowi)
                  qrq=q(qq)
                endif
                if(abs(ws(qrq+1)).le.tol)call trim_(qq,prq,qrq,q,ws)
c  rename qq as rowi and vice-versa
                if(qri.lt.qrq)then
                  if(s(rowi).eq.qq)then
                    r(qq)=r(rowi)
                    r(rowi)=qq
                    s(rowi)=s(qq)
                    s(qq)=rowi
                  else
                    call iexch(r(qq),r(rowi))
                    call iexch(s(qq),s(rowi))
                    r(s(qq))=qq
                    s(r(rowi))=rowi
                  endif
                  if(r(qq).gt.0)then
                    s(r(qq))=qq
                  else
                    irow=qq
                  endif
                  if(s(rowi).gt.0)r(s(rowi))=rowi
                else
                  if(s(qq).eq.rowi)then
                    r(rowi)=r(qq)
                    r(qq)=rowi
                    s(qq)=s(rowi)
                    s(rowi)=qq
                  else
                    call iexch(r(rowi),r(qq))
                    call iexch(s(rowi),s(qq))
                    r(s(rowi))=rowi
                    s(r(qq))=qq
                  endif 
                  if(r(rowi).gt.0)then
                    s(r(rowi))=rowi
                  else
                    irow=rowi
                  endif
                  if(s(qq).gt.0)r(s(qq))=qq
                endif
                call iexch(pri,prq)
                call iexch(qri,qrq)
                call iexch(q(rowi),q(qq))
                if(pri.eq.0)call erase(rowi,lastr,irow,r,s)
                prq=prq+1
              endif
            endif
            p(rowi)=pri
            p(qq)=prq
            ws(qrq+prq)=1.D0
            d(i)=rr*u11
            u11=u22
            l11=l21
          else
c  perm operation for row flma
c           write(nout,*)'perm operation for row flma'
            if(rr.ne.0.D0)then
              if(is.ge.0)then
                if(prq.gt.0)then
                  call mysaxpy(rr,ws(qrq+1),ws(qri+is+1),prq)
                  if(abs(ws(qri+1)).le.tol)call trim_(rowi,pri,qri,q,ws)
                  if(pri.eq.0)call erase(rowi,lastr,irow,r,s)
                endif
                is=pri-prq
              else
                pr_=pri
                pri=prq
                call newslot(rowi,pri,lastr,irow,p,q,r,s,ws,mxws,qr_,
     *            ifail)
                if(ifail.gt.0)return
                qrq=q(qq)
                qri=q(rowi)
                is=-is
                do j=1,is
                  ws(qri+j)=rr*ws(qrq+j)
                enddo
                call saxpyz(rr,ws(qrq+is+1),ws(qr_+1),ws(qri+is+1),pr_)
                is=0
              endif
            endif
            p(rowi)=pri
            if(u12.ne.0.D0)then
              u12=-u12/del
              if(is.gt.0)then
                pr_=prq
                prq=pri+1
                call newslot(qq,prq,lastr,irow,p,q,r,s,ws,mxws,qr_,
     *            ifail)
                if(ifail.gt.0)return
                qrq=q(qq)
                qri=q(rowi)
                do j=1,is
                  ws(qrq+j)=u12*ws(qri+j)
                enddo
                call saxpyz(u12,ws(qri+is+1),ws(qr_+1),ws(qrq+is+1),pr_)
                ws(qrq+prq)=u12
                goto7
              else
                if(pri.gt.0)then
                  is=-is
                  call mysaxpy(u12,ws(qri+1),ws(qrq+is+1),pri)
                  if(abs(ws(qrq+1)).le.tol)then
                    call trim_(qq,prq,qrq,q,ws)
                    if(prq.eq.0)call erase(qq,lastr,irow,r,s)
                    p(qq)=prq
                  endif
                endif
              endif
            endif
            if(prq.gt.0.or.u12.ne.0.D0)then
              if(prq.eq.0)then
                len=0
              elseif(s(qq).eq.0)then
                len=mxws-qrq
              else
                len=q(s(qq))-qrq
              endif
              if(len.eq.prq)then
                call newslot(qq,prq+1,lastr,irow,p,q,r,s,ws,mxws,qr_,
     *            ifail)     
                if(ifail.gt.0)return
                qrq=q(qq)
                qri=q(rowi)
                call r_shift(ws(qrq+1),prq,qr_-qrq)
              endif
              prq=prq+1
              ws(qrq+prq)=u12
            endif
   7        continue
            p(rowi)=pri
            p(qq)=prq
            d(i)=del
            u11=u11*u22/del
            call iexch(lc(i),lc(im))
          endif
c           write(nout,*)'rowi,pri,qri =',rowi,pri,qri
c           write(nout,*)'rowq,prq,qrq =',qq,prq,qrq
c           write(nout,*)'j,p(j),q(j),r(j),s(j)   irow =',irow
c           do j=1,n
c             if(p(j).ne.0)write(nout,*)j,p(j),q(j),r(j),s(j)
c           enddo
c           write(nout,*)'rowq* =',(ws(qrq+ij),ij=1,prq)
c           write(nout,*)'rowi* =',(ws(qri+ij),ij=1,pri)
        enddo
        if(prq.gt.0)then
c         write(nout,*)'ss,l11,ilast,n,prq',ss,l11,ilast,n,prq
c         write(nout,*)'sn =',(sn(ij),ij=nu+1,n)
          call mysaxpy(ss/l11,ws(qrq+1),sn(ilast-prq+1),prq)
          call erase(qq,lastr,irow,r,s)
          p(qq)=0
        endif
        qqq=lc(ilast)
        do i=ilast,nu+1,-1
          lc(i)=lc(i-1)
          li(lc(i))=i
        enddo
c       if(pp.le.n)then
c         ip=li(pp)
c         write(nout,*)'check sn'
c         do i=nu+1,ilast
c           nodec=lc(i)
c           u12=aiscpri2(n,a,la,pp,lc(i)-n,sn(nu+1),1.D0,ilast,
c             ilast-nu,li)
c           if(abs(u12).gt.tol)write(nout,*)'error,nodec =',u12,nodec
c         enddo
c       endif
c       write(nout,*)'intermediate PAQ factors:  new q =',qqq
c       write(nout,*)'lr =',(lr(j),j=nu+1,n)
c       write(nout,*)'lc =',(lc(j),j=nu+1,n)
c       write(nout,*)'d =',(d(ij),ij=nu+1,n)
c       do j=nu+1,n
c         rowj=lr(j)
c         if(p(rowj).ne.0)then
c           write(nout,*)'L(',rowj,')',
c    *        (ws(k),k=q(rowj)+1,q(rowj)+p(rowj))
c         endif
c       enddo
c       call checkout(n,a,la,lr,lc,li,p,q,r,s,ws,mxws,d)
      endif
      ip=li(pp)
      if(pp.gt.n)then
        li(pp)=0
        if(pp.eq.qqq)goto30
        if(ip.le.nu)goto99
        iout=ip
        rowim=lr(ip)
        prim=p(rowim)
        if(prim.gt.0)qrim=q(rowim)
      else
        if(ip.gt.nu.or.p(pp).gt.0)goto99
        lr(ip)=lr(nu)
        li(lr(ip))=ip
c  check for growth in sn
c       write(nout,*)'sn =',(sn(i),i=nu+1,n)
        iout=ilast
        i=nu+1
        if(i.gt.ilast)goto13
   11   continue
          do j=i,mao(i)
            if(abs(sn(j)).gt.growth)then
              iout=i-1
              goto13
            endif
          enddo
          i=mao(i)+1
          if(i.le.ilast)goto11
   13   continue
        do j=nu+1,iout
          if(abs(sn(j)).gt.tol)goto14
        enddo
        j=iout+1
   14   continue
        rowim=pp
        prim=iout-j+1
        if(prim.gt.0)then
          call newslot(pp,prim,lastr,irow,p,q,r,s,ws,mxws,qr_,ifail)
          if(ifail.gt.0)return
          p(pp)=prim
          qrim=q(pp)
          ii=qrim
          do j=j,iout
            ii=ii+1
            ws(ii)=sn(j)
          enddo
        endif
        do i=nu,iout-1
          lr(i)=lr(i+1)
          li(lr(i))=i
          lc(i)=lc(i+1)
          li(lc(i))=i
          d(i)=d(i+1)
        enddo
        lr(iout)=pp
        li(pp)=iout
c       write(nout,*)'lr =',(lr(ij),ij=nu,iout)
c       write(nout,*)'lc =',(lc(ij),ij=nu,iout-1)
c       if(prim.gt.0)write(nout,*)'L(',pp,') =',(ws(qrim+j),j=1,prim)
        nu=nu-1
      endif
c     write(nout,*)'iout,ilast,rowim,prim =',iout,ilast,rowim,prim
c  column flma operations to restore L to triangular form
      iswap=0
      do i=iout+1,ilast
        im=i-1
        lc(im)=lc(i)
        li(lc(im))=im
        rowi=lr(i)
        pri=p(rowi)
c       if(pri.gt.0)write(nout,*)'L(',rowi,') =',(ws(q(rowi)+j),j=1,pri)
        u22=d(i)
        if(prim.gt.0)then
          u12=aiscpri2(n,a,la,rowim,lc(i)-n,ws(qrim+1),1.D0,im-1,prim,
     *      li)
          if(abs(u12).le.tol)u12=0.D0
        else
          u12=aij(rowim,lc(i)-n,a,la)
        endif
        if(pri.gt.0)then
c         write(nout,*)'pri,iswap',pri,iswap
          qri=q(rowi)
          ii=pri-iswap
          if(ii.le.0)then
            l21=0.D0
          else
            l21=ws(qri+ii)
            if(abs(l21).le.tol)l21=0.D0
            if(ii.eq.1)then
              call trim_(rowi,pri,qri,q,ws)
              if(pri.eq.0)call erase(rowi,lastr,irow,r,s)
              if(s(rowi).eq.0)then
                qr_=mxws
              else
                qr_=q(s(rowi))
              endif
              if(qri+pri.ge.qr_)then
                call r_shift(ws(qri),pri,1)
                qri=qri-1
                q(rowi)=qri
              endif
            else
              pri=pri-1
              call r_shift(ws(qri+ii),iswap,1)
            endif
            p(rowi)=pri
c           write(nout,*)'rowi =',(ws(qri+ij),ij=1,pri)
          endif
        else
          l21=0.D0
        endif
        del=u22-l21*u12
        test=abs(u12)*max(1.D0,abs(l21))
c       write(nout,*)'l21,u12,u22,del,test',l21,u12,u22,del,test
        is=pri-prim
        if(is.gt.0)test=growth*test
        if(l21.eq.0.D0.and.is.lt.0)test=thresh*test
c         write(nout,*)'rowim,prim,qrim =',rowim,prim,qrim
c         write(nout,*)'rowi,pri,qri =',rowi,pri,qri
c         write(nout,*)'j,p(j),q(j),r(j),s(j)   irow =',irow
c         do j=1,n
c           if(p(j).ne.0)write(nout,*)j,p(j),q(j),r(j),s(j)
c         enddo
c         write(nout,*)'rowim =',(ws(qrim+ij),ij=1,prim)
c         write(nout,*)'rowi =',(ws(qri+ij),ij=1,pri)
        if(abs(del).le.test)then
c  no-perm operation for column flma
c         write(nout,*)'no-perm operation for column flma'
          rr=-u22/u12
          l21=l21+rr
          if(abs(l21).le.tol)l21=0.D0
          if(is.ge.0)then
            if(prim.gt.0)then
              call mysaxpy(rr,ws(qrim+1),ws(qri+is+1),prim)
              if(abs(ws(qri+1)).le.tol)call trim_(rowi,pri,qri,q,ws)
              if(pri.eq.0)then
                call erase(rowi,lastr,irow,r,s)
                p(rowi)=0
              endif
            endif
            if(pri.gt.0.or.l21.ne.0.D0)then
              if(pri.eq.0)then
                len=0
              elseif(s(rowi).eq.0)then
                len=mxws-qri
              else
                len=q(s(rowi))-qri
              endif
              if(len.eq.pri)then
                call newslot(rowi,pri+1,lastr,irow,p,q,r,s,ws,mxws,qr_,
     *            ifail)
                if(ifail.gt.0)return
                qrim=q(rowim)
                qri=q(rowi)
                call r_shift(ws(qri+1),pri,qr_-qri)
              endif
              pri=pri+1
              ws(qri+pri)=l21
            endif
          else
            pr_=pri
            pri=prim+1
            call newslot(rowi,pri,lastr,irow,p,q,r,s,ws,mxws,qr_,ifail)      
            if(ifail.gt.0)return
            qrim=q(rowim)
            qri=q(rowi)
            is=-is
            do j=1,is
              ws(qri+j)=rr*ws(qrim+j)
            enddo
            call saxpyz(rr,ws(qrim+is+1),ws(qr_+1),ws(qri+is+1),pr_)
            ws(qri+pri)=l21
          endif
c           write(nout,*)'rowim,prim,qrim =',rowim,prim,qrim
c           write(nout,*)'rowi,pri,qri =',rowi,pri,qri
c           write(nout,*)'j,p(j),q(j),r(j),s(j)   irow =',irow
c           do j=1,n
c             if(p(j).ne.0)write(nout,*)j,p(j),q(j),r(j),s(j)
c           enddo
c           write(nout,*)'rowim* =',(ws(qrim+ij),ij=1,prim)
c           write(nout,*)'rowi* =',(ws(q(rowi)+ij),ij=1,p(rowi))
          p(rowi)=pri
          rowim=rowi
          prim=pri
          qrim=qri
          d(im)=u12
c  perform accumulated cyclic permutation in subsequent rows
          if(iswap.gt.0)then
            do j=i+1,ilast
              rowj=lr(j)
              prj=p(rowj)
              is=prj-j+i
              if(is.gt.0)then
                qrj=q(rowj)
                if(is.gt.iswap)then
                  ii=is-iswap
                  l21=ws(qrj+ii)
                  call r_shift(ws(qrj+ii),iswap,1)
                  ws(qrj+is)=l21
                  if(abs(ws(qrj+1)).le.tol)call trim_(rowj,prj,qrj,q,ws)
                  if(prj.eq.0)call erase(rowj,lastr,irow,r,s)
                else
                  prj=prj+1
                  rrj=r(rowj)
                  if(rrj.eq.0)then
                    len=qrj
                  else
                    len=qrj-q(rrj)-p(rrj)
                  endif
                  if(len.gt.0)then
                    call r_shift(ws(qrj),is,1)
                    ws(qrj+is)=0.D0
                    qrj=qrj-1
                    q(rowj)=qrj
                  else
                    call newslot(rowj,prj,lastr,irow,p,q,r,s,ws,mxws,
     *                qr_,ifail) 
                    if(ifail.gt.0)return
                    qrj=q(rowj)
                    qrim=q(rowim)
                    call r_shift(ws(qrj+1),is,qr_-qrj)
                    ws(qrj+is+1)=0.D0
                    call r_shift(ws(qrj+is+2),j-i,qr_-qrj-1)
                  endif
                endif
                p(rowj)=prj
c               write(nout,*)'L(',rowj,')* =',(ws(qrj+ij),ij=1,prj)
              endif
            enddo
          endif
          iswap=0
        else
c  perm operation for column flma
c         write(nout,*)'perm operation for column flma'
          rr=-l21
          if(rr.ne.0.D0)then
            if(is.ge.0)then
              if(prim.gt.0)then
                call mysaxpy(rr,ws(qrim+1),ws(qri+is+1),prim)
                if(abs(ws(qri+1)).le.tol)call trim_(rowi,pri,qri,q,ws)
                if(pri.eq.0)call erase(rowi,lastr,irow,r,s)
              endif
              is=pri-prim
            else
              pr_=pri
              pri=prim
              call newslot(rowi,pri,lastr,irow,p,q,r,s,ws,mxws,qr_,
     *          ifail)
              if(ifail.gt.0)return
              qrim=q(rowim)
              qri=q(rowi)
              is=-is
              do j=1,is
                ws(qri+j)=rr*ws(qrim+j)
              enddo
              call saxpyz(rr,ws(qrim+is+1),ws(qr_+1),ws(qri+is+1),pr_)
              is=0
            endif
          endif
          p(rowi)=pri
          if(u12.ne.0.D0)then
            u12=-u12/del
            if(is.gt.0)then
              pr_=prim
              prim=pri+1
              call newslot(rowim,prim,lastr,irow,p,q,r,s,ws,mxws,qr_,
     *          ifail)
              if(ifail.gt.0)return
              qrim=q(rowim)
              qri=q(rowi)
              do j=1,is
                ws(qrim+j)=u12*ws(qri+j)
              enddo
              call saxpyz(u12,ws(qri+is+1),ws(qr_+1),ws(qrim+is+1),pr_)
              ws(qrim+prim)=u12
              goto27
            else
              if(pri.gt.0)then
                is=-is
                call mysaxpy(u12,ws(qri+1),ws(qrim+is+1),pri)
                if(abs(ws(qrim+1)).le.tol)then
                  call trim_(rowim,prim,qrim,q,ws)
                  if(prim.eq.0)call erase(rowim,lastr,irow,r,s)
                  p(rowim)=prim
                endif
              endif
            endif
          endif
          if(prim.gt.0.or.u12.ne.0.D0)then
            if(prim.eq.0)then
              len=0
            elseif(s(rowim).eq.0)then
              len=mxws-qrim
            else
              len=q(s(rowim))-qrim
            endif
            if(len.eq.prim)then
              call newslot(rowim,prim+1,lastr,irow,p,q,r,s,ws,mxws,qr_,
     *          ifail)
              if(ifail.gt.0)return
              qrim=q(rowim)
              qri=q(rowi)
              call r_shift(ws(qrim+1),prim,qr_-qrim)
            endif
            prim=prim+1
            ws(qrim+prim)=u12
          endif
   27     continue
          p(rowim)=prim
          p(rowi)=pri
c           write(nout,*)'rowim,prim,qrim =',rowim,prim,qrim
c           write(nout,*)'rowi,pri,qri =',rowi,pri,qri
c           write(nout,*)'j,p(j),q(j),r(j),s(j)   irow =',irow
c           do j=1,n
c             if(p(j).ne.0)write(nout,*)j,p(j),q(j),r(j),s(j)
c           enddo
c           write(nout,*)'rowim* =',(ws(qrim+ij),ij=1,prim)
c           write(nout,*)'rowi* =',(ws(q(rowi)+ij),ij=1,p(rowi))
          d(im)=del
          call iexch(lr(i),lr(i-1))
          call iexch(li(lr(i)),li(lr(i-1)))
          iswap=iswap+1
        endif
      enddo
      lc(ilast)=qqq
      li(qqq)=ilast
c     write(nout,*)'rowim* =',(ws(qrim+ij),ij=1,prim)
c     write(nout,*)'ilast,prim,qrim',ilast,prim,qrim
      if(prim.gt.0)then
       d(ilast)=aiscpri2(n,a,la,rowim,qqq-n,ws(qrim+1),1.D0,ilast-1,
     *    prim,li)
      else
        d(ilast)=aij(rowim,qqq-n,a,la)
      endif
c  reset mao
      iout=ilast
      do i=ilast,nu+1,-1
        mao(i)=ilast
        iout=min(iout,i-p(lr(i)))
        if(iout.eq.i)ilast=i-1
      enddo
   30 continue
      m1=n-nu
c     write(nout,*)'PAQ factors:  nu =',nu
c     write(nout,*)'d =',(d(ij),ij=nu+1,n)
c     do j=nu+1,n
c       rowj=lr(j)
c       if(p(rowj).ne.0)then
c         write(nout,*)'L(',rowj,')',
c    *      (ws(k),k=q(rowj)+1,q(rowj)+p(rowj))
c       endif
c     enddo
c     call checkout(n,a,la,lr,lc,li,p,q,r,s,ws,mxws,d)
c  print star diagram
c     if(m1.gt.80.or.n.gt.1000)stop
c     write(nout,*)'updated ordering:  nu =',nu
c     do i=1,n
c       do j=1,m1
c         star(i,j)=' '
c       enddo
c     enddo
c     do j=1,m1
c       jp=la(0)+lc(nu+j)-n
c       do i=la(jp),la(jp+1)-1
c         star(li(la(i)),j)='*'
c       enddo
c     enddo
c     do i=nu+1,n
c       write(nout,*)(star(i,j),j=1,m1)
c     enddo
c     write(nout,*)'lr =',(lr(j),j=nu+1,n)
c     write(nout,*)'lc =',(lc(j),j=nu+1,n)
c     write(nout,*)'mao =',(mao(j),j=nu+1,n)
      return
   99 continue
      write(nout,*)'malfunction in update_L:  p,q =',pp,qq
      stop
      end

      subroutine newslot(row,len,lastr,irow,p,q,r,s,ws,mxws,qr_,
     *  ifail)
      implicit double precision (a-h,u-z), integer (i-t)
      parameter (igap=10)
      dimension p(*),q(*),r(*),s(*),ws(*)
      common/noutc/nout
c     write(nout,*)'newslot: row =',row,'   len =',len
c     write(nout,*)'irow,lastr,mxws =',irow,lastr,mxws
      ifail=0
      if(lastr.eq.0)then
        if(mxws.lt.len)then
          write(nout,*)'insufficient space available for profile'
          ifail=7
        else
          irow=row
          q(row)=0
          r(row)=0
          s(row)=0
          lastr=row
        endif
        return
      endif
      igp=igap
    1 continue
      len_=len+igp
      thisr=lastr
    2 continue
      qrow=q(thisr)+p(thisr)
      nextr=s(thisr)
c     write(nout,*)'thisr,nextr,qrow,p(thisr),len_',
c    *  thisr,nextr,qrow,p(thisr),len_
      if(nextr.ne.0)then
        if(q(nextr).ge.qrow+len_)then
c  free slot after this row
          goto4
        else
          thisr=nextr
          if(thisr.ne.lastr)goto2
        endif
      else
        if(mxws-qrow.ge.len_)then
c  free slot at end of ws
          goto4
        elseif(q(irow).ge.len_)then
c  free slot at beginning of ws
          qrow=0
          thisr=0
          nextr=irow
          irow=row
          igp=0
          goto4
        endif
        thisr=irow
        if(thisr.ne.lastr)goto2
      endif
c  no free space: try minimum value of len
      if(igp.gt.0)then
        igp=0
        goto1
      endif
c  compress ws
      thisr=irow
      qrow=0
    3 continue
      call r_shift(ws(qrow+1),p(thisr),q(thisr)-qrow)
      q(thisr)=qrow
      qrow=qrow+p(thisr)
      if(s(thisr).ne.0)then
        thisr=s(thisr)
        goto3
      endif
      if(mxws.lt.qrow+len_)then
        write(nout,*)'insufficient space available for profile'
        write(nout,*)'mxws,qrow,len_',mxws,qrow,len_
        ifail=7
        return
      endif
c  insert at end of compressed file
      nextr=0
    4 continue
      qr_=q(row)
      q(row)=qrow+igp
      if(p(row).gt.0)then
        if(r(row).eq.thisr.or.s(row).eq.nextr)return
c  insert after row thisr and take out old row
        call erase(row,lastr,irow,r,s)
      endif
      lastr=row
      r(row)=thisr
      if(thisr.gt.0)s(thisr)=row
      s(row)=nextr
      if(nextr.gt.0)r(nextr)=row
      i=0
      return
      end

      subroutine erase(row,lastr,irow,r,s)
c  remove slot for row from the data file
      implicit integer (i-s)
      dimension r(*),s(*)
      common/noutc/nout
c     write(nout,*)'erase: row,irow,lastr =',row,irow,lastr
      if(r(row).eq.0)then
        if(s(row).eq.0)then
          irow=0
          lastr=0
          return
        endif
        irow=s(row)
        r(irow)=0
      elseif(s(row).eq.0)then
        s(r(row))=0
      else
        s(r(row))=s(row)
        r(s(row))=r(row)
      endif
      if(row.eq.lastr)lastr=irow
      return
      end

      subroutine trim_(rowi,pri,qri,q,ws)
c  trim leading zeros off slot for row i
      implicit double precision (a-h,s-z), integer (i-r)
      dimension q(*),ws(*)
      common/epsc/eps,tol,emin
    1 continue
      qri=qri+1
      pri=pri-1
      if(pri.eq.0)return
      if(abs(ws(qri+1)).le.tol)goto1
      q(rowi)=qri
      return
      end

      subroutine checkout(n,a,la,lr,lc,li,p,q,r,s,ws,mxws,d)
      implicit double precision (a-h,r-z), integer (i-q)
      integer r,s,rowj,thisr
      dimension a(*),la(*),lr(*),lc(*),li(*),p(*),q(*),r(*),s(*),ws(*),
     *  d(*)
      common/factorc/m1,nu,mp,mq,lastr,irow
      common/noutc/nout
      common/epsc/eps,tol,emin
c  check indexing
      do j=1,nu
        if(p(lr(j)).ne.0)then
          write(nout,*)'p(lr(j)).ne.0'
          goto11
        endif
      enddo
      np=0
      do i=nu+1,n
        if(p(lr(i)).gt.0)np=np+1
      enddo
      if(irow.gt.0)then
        if(r(irow).ne.0)then
          write(nout,*)'r(irow).ne.0'
          goto11
        endif
        thisr=irow
    1   continue
        if(p(thisr).le.0)then
          write(nout,*)'p(thisr).le.0'
          goto11
        endif
        np=np-1
        nextr=s(thisr)
        if(nextr.eq.0)then
          if(q(thisr)+p(thisr).gt.mxws)then
            write(nout,*)'q(thisr)+p(thisr).gt.mxws'
            goto11
          endif
        else
          if(r(nextr).ne.thisr)then
            write(nout,*)'r(nextr).ne.thisr'
            goto11
          endif
          if(nextr.ne.s(thisr))then
            write(nout,*)'nextr.ne.s(thisr)'
            goto11
          endif
          if(q(thisr)+p(thisr).gt.q(nextr))then
            write(nout,*)'q(thisr)+p(thisr).gt.q(nextr)'
            goto11
          endif
          thisr=nextr
          goto1
        endif
      endif
      if(np.ne.0)then
        write(nout,*)'np.ne.0'
        goto11
      endif
      last=0
      emax=0.D0
      length=0
      do inode=nu+1,n
        nodec=lc(inode)
c  form L.a_q
        rowj=lr(inode)
        prj=p(rowj)
        length=length+prj
        if(prj.lt.0)then
          write(nout,*)'prj.lt.0'
          goto11
        elseif(prj.eq.0)then
          e=abs(aij(rowj,nodec-n,a,la)-d(inode))
        else
          e=abs(d(inode)-aiscpri2(n,a,la,rowj,nodec-n,ws(q(rowj)+1),
     *      1.D0,inode-1,prj,li))
        endif
c       if(e.gt.tol)write(nout,*)'error =',e,
c    *    '  inode,nodec,rowj =',inode,nodec,rowj
        emax=max(emax,e)
        do j=inode+1,n
          rowj=lr(j)
          prj=p(rowj)
          if(prj.gt.0)then
            e=abs(aiscpri2(n,a,la,rowj,nodec-n,ws(q(rowj)+1),1.D0,j-1,
     *         prj,li))
          else
            e=abs(aij(rowj,nodec-n,a,la))
          endif
c         if(e.gt.tol)write(nout,*)'error =',e,
c    *      '  inode,nodec,j,rowj =',inode,nodec,j,rowj
          emax=max(emax,e)
        enddo
      enddo
      write(nout,*)'checkout:  m1 =',m1,'  file length =',length
      if(emax.gt.tol)write(nout,*)'error =',emax
      return
   11 continue
      write(nout,*)'thisr,nextr =',thisr,nextr
      write(nout,*)'i,p(i),q(i),r(i),s(i):  irow =',irow
      do i=1,n
        if(p(i).ne.0)write(nout,*)i,p(i),q(i),r(i),s(i)
      enddo
      stop
      end


C  --- Begin file: util.f

c  Copyright (C) 1996 Roger Fletcher

c  Current version dated 26 May 2011

c  THE ACCOMPANYING PROGRAM IS PROVIDED UNDER THE TERMS OF THE ECLIPSE PUBLIC
c  LICENSE ("AGREEMENT"). ANY USE, REPRODUCTION OR DISTRIBUTION OF THE PROGRAM
c  CONSTITUTES RECIPIENT'S ACCEPTANCE OF THIS AGREEMENT

c*********************** dense matrix utilities ************************

      subroutine rsol(n,nn,nmax,R,b)
      implicit double precision (a-h,o-z)
      dimension R(*),b(*)
c  solves Rx=b where R is nxn upper triangular. Solution overwrites b.
c  R is a single suffix array: the first nmax elements contain the first row
c  of R in positions 1:n, the next nmax-1 elements contain the second row of R,
c  and so on. nn indexes the element R(n,n) (where nn=n*(3-n)/2+(n-1)*nmax)
      n1=nmax+1
      ii=nn
      b(n)=b(n)/R(nn)
      do i=n-1,1,-1
        ii=ii-n1+i
        b(i)=-scpr(-b(i),R(ii+1),b(i+1),n-i)/R(ii)
      enddo
      return
      end

      subroutine rtsol(n,nn,nmax,R,b)
      implicit double precision (a-h,o-z)
      dimension R(*),b(*)
c  solves Rt.x=b with same conventions as above
c  nn is not required on entry but is set on exit
      n2=nmax+2
      nn=1
      b(1)=b(1)/R(1)
      do i=2,n
        i1=i-1
        call mysaxpy(-b(i1),R(nn+1),b(i),n-i1)
        nn=nn+n2-i
        b(i)=b(i)/R(nn)
      enddo
      return
      end

      subroutine Qprod(n,nmax,Q,x,b)
      implicit double precision (a-h,o-z)
      dimension Q(*),x(*),b(*)
c  forms b=M.x where Q is nxn, stored by columns, with stride nmax
      do i=1,n
        b(i)=0.D0
      enddo
      i1=1
      do i=1,n
        call mysaxpy(x(i),Q(i1),b,n)
        i1=i1+nmax
      enddo
      return
      end

      subroutine Qtprod(n,nmax,Q,x,b)
      implicit double precision (a-h,o-z)
      dimension Q(*),x(*),b(*)
c  forms b=M'.x where Q is nxn, stored by columns, with stride nmax
      i1=1
      do i=1,n
        b(i)=scpr(0.D0,Q(i1),x,n)
        i1=i1+nmax
      enddo
      return
      end

      subroutine brots(n,nmax,k,kk,R,v)
      implicit double precision (a-h,o-z)
      dimension R(*),v(*)
      ipip=kk
      do i=k-1,1,-1
        ip=i+1
        ipi=ipip-nmax+i
        ii=ipi-1
        call angle(v(i),v(ip),cos,sin)
        call rot(n-i,R(ipi),R(ipip),cos,sin)
        v(ip)=sin*R(ii)
        R(ii)=cos*R(ii)
        ipip=ii
      enddo
      return
      end

      subroutine frots(nr,nc,nmax,R,v)
      implicit double precision (a-h,o-z)
      dimension R(*),v(*)
c nr is either nc or nc+1
      ii=1
      do i=1,nc
        ip=i+1
        ipi=ii+1
        ipip=ipi+nmax-i
        call angle(R(ii),v(ip),cos,sin)
        call rot(nr-i,R(ipi),R(ipip),cos,sin)
        ii=ipip
      enddo
      return
      end

      subroutine angle(a,b,cos,sin)
      implicit double precision (a-h,o-z)
      z=sqrt(a**2+b**2)
      if(z.eq.0.D0)then
        cos=1.D0
        sin=0.D0
        return
      endif
      cos=a/z
      sin=b/z
      a=z
      b=0.D0
      return
      end

      subroutine rot(n,a,b,cos,sin)
      implicit double precision (a-h,o-z)
      dimension a(*),b(*)
      if(sin.eq.0.D0)then
        if(cos.gt.0.D0)then
          do i=1,n
            b(i)=-b(i)
          enddo
        else
          do i=1,n
            a(i)=-a(i)
          enddo
        endif
      elseif(cos.eq.0.D0)then
        if(sin.ge.0.D0)then
          do i=1,n
            z=a(i)
            a(i)=b(i)
            b(i)=z
          enddo
        else
          do i=1,n
            z=a(i)
            a(i)=-b(i)
            b(i)=-z
          enddo
        endif
      else
        do i=1,n
          z=a(i)
          a(i)=cos*z+sin*b(i)
          b(i)=sin*z-cos*b(i)
        enddo
      endif
      return
      end

      subroutine mysaxpy(a,x,y,n)
      implicit double precision (a-h,o-z)
      dimension x(*),y(*)
      if(a.eq.0.D0)return
      do i=1,n
        y(i)=y(i)+a*x(i)
      enddo
      return
      end

      subroutine saxpys(a,x,is,y,n)
      implicit double precision (a-h,o-z)
c  saxpy with stride
      dimension x(*),y(*)
      if(a.eq.0.D0)return
      ix=1
      do i=1,n
        y(i)=y(i)+a*x(ix)
        ix=ix+is
      enddo
      return
      end

      subroutine saxpyx(a,x,y,n)
      implicit double precision (a-h,o-z)
c  saxpy with result in x
      dimension x(*),y(*)
      if(a.eq.0.D0)then
        do i=1,n
          x(i)=y(i)
        enddo
      else
        do i=1,n
          x(i)=y(i)+a*x(i)
        enddo
      endif
      return
      end

      subroutine saxpyz(a,x,y,z,n)
      implicit double precision (a-h,o-z)
c  saxpy with result in z
      dimension x(*),y(*),z(*)
      if(a.eq.0.D0)then
        do i=1,n
          z(i)=y(i)
        enddo
      else
        do i=1,n
          z(i)=y(i)+a*x(i)
        enddo
      endif
      return
      end

      subroutine saxpyi(a,x,y,n)
      implicit double precision (a-h,o-z)
c  saxpy with interchange of x and y
      dimension x(*),y(*)
      if(a.eq.0.D0)then
        do i=1,n
          call rexch(x(i),y(i))
        enddo
      else
        do i=1,n
          z=y(i)
          y(i)=x(i)+a*y(i)
          x(i)=z
        enddo
      endif
      return
      end

      function scpr(a,x,y,n)
      implicit double precision (a-h,o-z)
      dimension x(*),y(*)
      scpr=a
      do i=1,n
        scpr=scpr+x(i)*y(i)
      enddo
      return
      end

c     function xlen(a,x,n)
c     implicit double precision (a-h,o-z)
c     dimension x(*)
c  finds the l_2 length of [a:x] where a is either 0.D0 or 1.D0
c  if overflow occurs the function is calculated in a less efficient way.
c  Users who cannot trap overflow should either use this method of calculation,
c  or use the alternative routine "xlen" below which is not quite so well
c  protected against overflow.
c     external  ieee_handler, abort
c     integer   ieee_flags, ieeer, ieee_handler
c     external  ieee_flags
c     character out*16
c     out = ''
c     ieeer = ieee_flags ( 'clearall','all','',out )
c     ieeer=ieee_handler('clear','overflow',abort)
c  this call of ieee_handler assumes that 
c         ieeer=ieee_handler('set','overflow',abort)
c  has been set in the driver. If not this call of ieee_handler and that below
c  should be removed
c     xlen=a
c     do i=1,n
c       xlen=xlen+x(i)**2
c     enddo
c     xlen=sqrt(xlen)
c     ieeer=ieee_flags ( 'get','exception','',out )
c     if(out.eq.'overflow')then
c       call linf(n,x,xmx,i)
c       xmx=max(xmx,1.D0) %this is needed if normalization is always used
c       xlen=(a/xmx)**2
c       do i=1,n
c         xlen=xlen+(x(i)/xmx)**2
c       enddo
c       xlen=xmx*sqrt(xlen)
c       ieeer=ieee_flags ( 'clear','overflow','',out )
c     endif
c     ieeer=ieee_handler('set','overflow',abort)
c     return
c     end
      
      function xlen(a,x,n)
      implicit double precision (a-h,o-z)
      dimension x(*)
      xlen=a
      do i=1,n
        xlen=xlen+x(i)**2
      enddo
      xlen=sqrt(xlen)
      return
      end
      
      subroutine linf(n,x,z,iz)
      implicit double precision (a-h,o-z)
      dimension x(*)
      z=0.D0
      do i=1,n
        a=abs(x(i))
        if(a.gt.z)then
          z=a
          iz=i
        endif
      enddo
      return
      end

      subroutine r_shift(r,n,k)
      implicit double precision (a-h,o-z)
      dimension r(*)
      if(k.gt.0)then
        do i=1,n
          r(i)=r(i+k)
        enddo
      elseif(k.lt.0)then
        do i=n,1,-1
          r(i)=r(i+k)
        enddo
      endif
      return
      end

      subroutine ishift(l,n,k)
      implicit double precision (a-h,o-z)
      dimension l(*)
      if(k.gt.0)then
        do i=1,n
          l(i)=l(i+k)
        enddo
      elseif(k.lt.0)then
        do i=n,1,-1
          l(i)=l(i+k)
        enddo
      endif
      return
      end

      subroutine rexch(a,b)
      double precision a,b,z
      z=a
      a=b
      b=z
      return
      end

      subroutine vexch(a,b,n)
      double precision a,b,z
      dimension a(*),b(*)
      do i=1,n
        z=a(i)
        a(i)=b(i)
        b(i)=z
      enddo
      return
      end

      subroutine iexch(i,j)
      k=i
      i=j
      j=k
      return
      end
