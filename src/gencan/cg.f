      subroutine cgf(nind,ind,n,x,m,lambda,rho,g,delta,l,u,eps,epsnqmp,
     +maxitnqmp,maxit,nearlyq,gtype,htvtype,trtype,iprint,ncomp,s,iter,
     +rbdtype,rbdind,inform,w,y,r,d,sprev,theta,sterel,steabs,epsrel,
     +epsabs,infrel,infabs)

      implicit none

C     SCALAR ARGUMENTS
      logical nearlyq
      integer gtype,htvtype,inform,iprint,iter,m,maxit,maxitnqmp,n,
     +        ncomp,nind,trtype,rbdind,rbdtype
      double precision delta,eps,epsnqmp,epsabs,epsrel,infrel,infabs,
     +        steabs,sterel,theta

C     ARRAY ARGUMENTS
      integer ind(nind)
      double precision d(n),g(n),l(n),lambda(m),r(n),rho(m),s(n),
     +        sprev(n),u(n),w(n),x(n),y(n)

C     This subroutine implements the Conjugate Gradients method for 
C     minimizing the quadratic approximation q(s) of f(x) at x, where
C
C     q(s) = 1/2 s^T H s + g^T s,
C
C        H = \nabla^2 f(x),
C
C        g = \nabla f(x),
C
C     subject to || s || <= delta and l <= x + s <= u.
C
C     In the constraint ''|| s || <= delta'', the norm will be the
C     Euclidian norm if the input parameter trtype is equal to 0, and
C     it will be the Sup norm if trtype is equal to 1.
C
C     The method returns an approximation s to the solution such that 
C     ||H s + g||_2 <= eps * ||g||_2; or converges to the boundary of 
C     ||s||_2 <= delta and l <= x + s <= u; or finds a point s and a 
C     direction d such that q(s + alpha d) = q(s) for any alpha, i.e., 
C     d^T H d = g^T d = 0.
C
C     On Entry
C
C     nind     integer
C              number of free variables (this is thee dimension in 
C              which this subroutine will work)
C
C     ind      integer ind(n)
C              array which contains, in the first nind positions, the
C              identifiers of the free variables
C
C     n        integer
C              dimension of the full space
C
C     x        double precision x(n)
C              point at which f function is being approximated by the
C              quadratic model
C
C              The first nind positions of x contains the free variables 
C              x_ind(1), x_ind(2), ..., x_ind(nind).
C
C     m        integer
C     lambda   double precision lambda(m)
C     rho      double precision rho(m)
C              These three parameters are not used nor modified by 
C              GENCAN and they are passed as arguments to the user-
C              defined subroutines evalal and evalnal to compute the 
C              objective function and its gradient, respectively. 
C              Clearly, in an Augmented Lagrangian context, if GENCAN is 
C              being used to solve the bound-constrainted subproblems, m 
C              would be the number of constraints, lambda the Lagrange 
C              multipliers approximation and rho the penalty parameters
C
C     g        double precision g(n)
C              linear coefficient of the quadratic function
C
C              This is \nabla f(x) and it also contains in the first 
C              nind positions the components g_ind(1), g_ind(2), ..., 
C              g_ind(nind).
C    
C              IMPORTANT: the linear algebra of this subroutine lies in 
C              a space of dimension nind. The value of the full 
C              dimension n, the non-free variables (which are at the end 
C              of array x) and its gradient components (which are at the 
C              and of array g) are, at this moment, being used to 
C              approximate the Hessian times vector products by 
C              incremental quotients.
C
C     delta    double precision
C              trust region radius (||s||_2 <= delta) 
C    
C     l        double precision l(n)
C              lower bounds on x + s. It components are ordered in the 
C              same way as x and g.
C
C     u        double precision u(n)
C              upper bounds on x + s. It components are ordered in the 
C              same way as x, g and l.
C
C     eps      double precision
C              tolerance for the stopping criterion ||H s + g||_2 < eps 
C              * ||g||_2
C
C     epsnqmp  double precision
C              See below
C
C     maxitnqmp integer
C              This and the previous one parameter are used for a 
C              stopping criterion of the conjugate gradient 
C              subalgorithm. If the progress in the quadratic model is 
C              less or equal than a fraction of the best progress 
C              ( epsnqmp * bestprog ) during maxitnqmp consecutive 
C              iterations then CG is stopped by not enough progress of 
C              the quadratic model.
C
C              RECOMMENDED: epsnqmp = 1.0d-4, maxitnqmp = 5
C
C     maxit    integer
C              maximum number of iterations allowed
C
C     nearlyq  logical
C              if function f is (nearly) quadratic, use the option 
C              nearlyq = TRUE. Otherwise, keep the default option.
C
C              if, in an iteration of CG we find a direction d such that
C              d^T H d <= 0 then we take the following decision:
C
C              (i) if nearlyq = TRUE then take direction d and try to go 
C              to the boundary choosing the best point among the two 
C              point at the boundary and the current point. 
C
C              (ii) if nearlyq = FALSE then we stop at the current 
C              point.
C
C              RECOMMENDED: nearlyq = FALSE
C
C     gtype    integer
C              type of gradient calculation
C              gtype = 0 means user suplied evalg subroutine,
C              gtype = 1 means central difference approximation.
C
C              RECOMMENDED: gtype = 0
C
C              (provided you have the evalg subroutine)
C
C     htvtype  integer
C              type of Hessian times vector product calculation
C              htvtype = 0 means user supplied evalhd subroutine,
C              htvtype = 1 means incremental quotients approximation.
C
C              RECOMMENDED: htvtype = 1
C
C              (you take some risk using this option but, unless you 
C              have a good evalhd subroutine, incremental quotients is a
C              very cheap option)
C
C     trtype   integer
C              type of trust-region radius
C              trtype = 0 means 2-norm trust-region
C              trtype = 1 means infinite-norm trust-region
C
C              RECOMMENDED: trtype = 0
C
C     iprint   integer
C              Commands printing. Nothing is printed if iprint is 
C              smaller than 2. If iprint is greater than or equal to 
C              2, GENCAN iterations information is printed. If iprint 
C              is greater than or equal to 3, line searches and 
C              Conjugate Gradients information is printed.
C
C              RECOMMENDED: iprint = 2
C
C              CONSTRAINTS: allowed values are just 2 or 3.
C
C     ncomp    integer
C              This constant is just for printing. In a detailed 
C              printing option, ncomp component of some vectors will be 
C              printed
C
C              RECOMMENDED: ncomp = 5
C
C              CONSTRAINTS: ncomp >= 0
C
C     w        double precision w(n)
C     y        double precision y(n)
C     r        double precision r(n)
C     d        double precision d(n)
C     sprev    double precision sprev(n)
C              working vectors
C
C     theta    double precision
C              constant for the angle condition, i.e., at iteration k we 
C              need a direction d_k such that <gk,dk> <= - theta 
C              ||gk||_2 ||dk||_2, where gk is \nabla f(xk)
C
C              RECOMMENDED: theta = 10^{-6}
C
C     sterel   double precision
C     steabs   double precision
C              this constants mean a ``relative small number'' and ``an 
C              absolute small number'' for the increments in finite
C              difference approximations of derivatives
C
C              RECOMMENDED: epsrel = 10^{-7}, epsabs = 10^{-10} 
C
C     epsrel   double precision
C     epsabs   double precision
C     infrel   double precision
C     infabs   double precision 
C              this constants mean a ``relative small number'', ``an 
C              absolute small number'', and ``infinite or a very big 
C              number''. Basically, a quantity A is considered 
C              negligible with respect to another quantity B if |A| < 
C              max ( epsrel * |B|, epsabs ) 
C
C              RECOMMENDED: epsrel = 10^{-10}, epsabs = 10^{-20},
C                           infrel = 10^{+20}, infabs = 10^{+99}
C
C     On Return
C
C     s        double precision s(n)
C              final estimation of the solution
C
C     iter     integer
C              number of Conjugate Gradient iterations performed
C
C     inform   integer
C              termination parameter:
C
C              0 = convergence with ||H s + g||_2 <= eps * ||g||_2;
C
C              1 = convergence to the boundary of ||s||_2 <= delta;
C
C              2 = convergence to the boundary of l - x <= s <= u - x;
C
C              3 = stopping with s = sk  such that <gk,sk> <= -t heta 
C                  ||gk||_2 ||sk||_2 and <gk,s_{k+1}> > - theta 
C                  ||gk||_2 ||s_{k+1}||_2;
C
C              4 = not enough progress of the quadratic model during
C                  maxitnqmp iterations, i.e., during maxitnqmp 
C                  iterations | q - qprev | <= max ( epsrel * | q |, 
C                  epsabs );
C
C              6 = very similar consecutive iterates, for two 
C                  consecutive iterates x and y, for all i | x(i) - 
C                  y(i) | <= max ( epsrel * | x(i) |, epsabs );
C
C              7 = stopping with d such that d^T H d = 0 and g^T d = 0;
C
C              8 = too many iterations;
C
C            < 0 = error in evalhd subroutine.

C     LOCAL SCALARS
      character * 5 rbdtypea
      logical samep,detail_trace
      integer i,itnqmp,rbdnegaind,rbdnegatype,rbdposaind,rbdposatype,
     +        detail_call_target,detail_status,cg_call_id
      integer, save :: cg_call_counter = 0
      double precision aa,alpha,amax,amax1,amax1n,amaxn,amax2,amax2n,
     +        amax2nx,amax2x,bb,bestprog,beta,cc,currprog,dd,dnorm2,dtr,
     +        dts,dtw,gnorm2,gts,norm2s,q,qamax,qamaxn,qprev,rnorm2,
     +        rnorm2prev,snorm2,snorm2prev
      character(len=32) detail_value

C     ==================================================================
C     Initialization
C     ==================================================================

      cg_call_counter = cg_call_counter + 1
      cg_call_id = cg_call_counter
      detail_call_target = 0
      detail_value = ''
      call get_environment_variable('PACKMOL_GENCAN_CG_COMPARE_CALL',
     + detail_value, status=detail_status)
      if ( detail_status .eq. 0 ) then
          read(detail_value,*,err=10) detail_call_target
      end if
 10   continue
      detail_trace = detail_call_target .gt. 0 .and.
     + detail_call_target .eq. cg_call_id

      gnorm2   = norm2s(nind,g)

      iter     =      0
      itnqmp   =      0
      qprev    = infabs
      bestprog =  0.0d0

      do i = 1,nind
          s(i) = 0.0d0
          r(i) =  g(i)
      end do

      q        =  0.0d0
      gts      =  0.0d0
      snorm2   =  0.0d0
      rnorm2   = gnorm2

      if ( detail_trace ) then
          write(*,*) '[gencan-fortran-cg] entry call=', cg_call_id,
     +     ' nind=', nind, ' delta=', delta, ' eps=', eps,
     +     ' gnorm2=', gnorm2
      end if

C     ==================================================================
C     Print initial information
C     ==================================================================

      if ( iprint .ge. 4 ) then
          write(*, 980) maxit,eps
          if ( trtype .eq. 0 ) then
              write(*, 981) delta
          else if ( trtype .eq. 1 ) then
              write(*, 982) delta
          else
              write(*, 983)
          end if
          write(*, 984) iter,rnorm2,sqrt(snorm2),q

          write(10,980) maxit,eps
          if ( trtype .eq. 0 ) then
              write(10,981) delta
          else if ( trtype .eq. 1 ) then
              write(10,982) delta
          else
              write(10,983)
          end if
          write(10,984) iter,rnorm2,sqrt(snorm2),q

      end if

C     ==================================================================
C     Main loop
C     ==================================================================

 100  continue

C     ==================================================================
C     Test stopping criteria
C     ==================================================================

C     if ||r||_2 = ||H s + g||_2 <= eps * ||g||_2 then stop

      if ( rnorm2 .le. 1.0d-16 .or.
     +     ( ( rnorm2 .le. eps ** 2 * gnorm2 .or.
     +       ( rnorm2 .le. 1.0d-10 .and. iter .ne. 0 ) ) 
     +       .and. iter .ge. 4 ) ) then

          inform = 0

          if ( iprint .ge. 4 ) then
              write(*, 990) inform
              write(10,990) inform
          end if
  
          go to 500

      end if

C     if the maximum number of iterations was achieved then stop

      if ( iter .ge. max(4, maxit) ) then

          inform = 8

          if ( iprint .ge. 4 ) then
              write(*, 998) inform
              write(10,998) inform
          end if
  
          go to 500

      end if

C     ==================================================================
C     Compute direction
C     ==================================================================

      if ( iter .eq. 0 ) then

          do i = 1,nind
              d(i) = - r(i)
          end do

          dnorm2 =   rnorm2
          dtr    = - rnorm2

      else

          beta = rnorm2 / rnorm2prev

          do i = 1,nind
              d(i) = - r(i) + beta * d(i)
          end do

          dnorm2 = rnorm2 - 2.0d0 * beta * ( dtr + alpha * dtw ) + 
     +             beta ** 2 * dnorm2
          dtr    = - rnorm2 + beta * ( dtr + alpha * dtw )

      end if

C     Force d to be a descent direction of q(s), i.e.,
C     <\nabla q(s), d> = <H s + g, d> = <r, d> \le 0.

      if ( dtr .gt. 0.0d0 ) then

          do i = 1,nind
              d(i) = - d(i)
          end do
          dtr = - dtr

      end if

C     ==================================================================
C     Compute d^T H d
C     ==================================================================

C     w = A d

      if ( htvtype .eq. 0 ) then
          call calchd(nind,ind,x,d,g,n,x,m,lambda,rho,w,y,sterel,steabs,
     +    inform)

      else if ( htvtype .eq. 1 ) then
          call calchddiff(nind,ind,x,d,g,n,x,m,lambda,rho,gtype,w,y,
     +    sterel,steabs,inform)
      end if

      if ( inform .lt. 0 ) then

          if ( iprint .ge. 4 ) then
              write(*, 1000) inform
              write(10,1000) inform
          end if

          return
      end if

C     Compute d^T w and ||w||^2

      dtw = 0.0d0
      do i = 1,nind
          dtw = dtw + d(i) * w(i)
      end do 

      if ( detail_trace ) then
          write(*,*) '[gencan-fortran-cg] step-pre call=', cg_call_id,
     +     ' iter=', iter, ' dtr=', dtr, ' dtw=', dtw,
     +     ' rnorm2=', rnorm2, ' dnorm2=', dnorm2
      end if

C     ==================================================================
C     Compute maximum step
C     ==================================================================

C     amax1 > 0 and amax1n < 0 are the values of alpha such that 
C     ||s + alpha * d||_2 or ||s + alpha * d||_\infty = delta 

      dts = 0.0d0
      do i = 1,nind
          dts = dts + d(i) * s(i)
      end do

C     Euclidian-norm trust radius

      if ( trtype .eq. 0 ) then

          aa = dnorm2
          bb = 2.0d0 * dts
          cc = snorm2 - delta ** 2
          dd = sqrt( bb ** 2 - 4.0d0 * aa * cc )

          amax1  = ( - bb + dd ) / ( 2.0d0 * aa )
          amax1n = ( - bb - dd ) / ( 2.0d0 * aa )

C     Sup-norm trust radius

      else if ( trtype .eq. 1 ) then

          amax1  =  infabs
          amax1n = -infabs

          do i = 1,nind
              if ( d(i) .gt. 0.0d0 ) then
                  amax1  = min( amax1,  (   delta - s(i) ) / d(i) )
                  amax1n = max( amax1n, ( - delta - s(i) ) / d(i) )
              else if ( d(i) .lt. 0.0d0 ) then
                  amax1  = min( amax1,  ( - delta - s(i) ) / d(i) )
                  amax1n = max( amax1n, (   delta - s(i) ) / d(i) )
              end if
          end do

      end if

C     amax2 > 0 and amax2n < 0 are the maximum and the minimum values of 
C     alpha such that l - x <= s + alpha * d <= u - x, respectively

      amax2  =   infabs
      amax2n = - infabs

      do i = 1,nind
          if ( d(i) .gt. 0.0d0 ) then
C             if (u(i).lt.infrel) then
                  amax2x = ( u(i) - x(i) - s(i) ) / d(i)
                  if ( amax2x .lt. amax2 ) then
                      amax2       = amax2x
                      rbdposaind  = i
                      rbdposatype = 2
                  end if
C             end if
C             if (l(i).gt.-infrel) then    
                  amax2nx = ( l(i) - x(i) - s(i) ) / d(i)
                  if ( amax2nx .gt. amax2n ) then
                      amax2n      = amax2nx
                      rbdnegaind  = i
                      rbdnegatype = 1
                  end if
C             end if
          else if ( d(i) .lt. 0.0d0 ) then
C             if (l(i).gt.-infrel) then    
                  amax2x = ( l(i) - x(i) - s(i) ) / d(i)
                  if ( amax2x .lt. amax2 ) then
                      amax2       = amax2x
                      rbdposaind  = i
                      rbdposatype = 1
                  end if
C             end if
C             if (u(i).lt.infrel) then
                  amax2nx = ( u(i) - x(i) - s(i) ) / d(i)
                  if ( amax2nx .gt. amax2n ) then
                      amax2n      = amax2nx
                      rbdnegaind  = i
                      rbdnegatype = 2
                  end if
C             end if
          end if
      end do

C     Compute amax as the minimum among amax1 and amax2, and amaxn as 
C     the minimum among amax1n and amax2n. Moreover change amaxn by 
C     - amaxn to have amax and amaxn as maximum steps along d direction 
C     (and not -d in the case of amaxn)

      amax  = min( amax1 , amax2  )
      amaxn = max( amax1n, amax2n )

      if ( detail_trace ) then
          write(*,*) '[gencan-fortran-cg] bounds call=', cg_call_id,
     +     ' iter=', iter, ' amax1=', amax1, ' amax1n=', amax1n,
     +     ' amax2=', amax2, ' amax2n=', amax2n
          write(*,*) '[gencan-fortran-cg] bounds-ext call=', cg_call_id,
     +     ' iter=', iter, ' amax=', amax, ' amaxn=', amaxn
      end if

C     ==================================================================
C     Compute the step (and the quadratic functional value at the new 
C     point)
C     ==================================================================

      qprev = q

C     If d^T H d > 0 then take the conjugate gradients step

      if ( dtw .gt. 0.0d0 ) then

          alpha = min( amax, rnorm2 / dtw )

          q = q + 0.5d0 * alpha ** 2 * dtw + alpha * dtr

C     If d^T H d <= 0 and function f is nearly quadratic then take the 
C     point with the minimum functional value (q) among the current one 
C     and the ones which are at the boundary, i.e., the best one between 
C     q(s), q(s + amax*d) and q(s + amaxn*d).

      else

          qamax = q + 0.5d0 * amax ** 2 * dtw + amax * dtr

C         If we are at iteration zero then take the maximum positive
C         step in the minus gradient direction

          if ( iter .eq. 0 ) then

              alpha = amax
              q     = qamax

C         If we are not in the first iteration then if function f is 
C         nearly quadratic and q(s + amax * d) or q(s + amaxn * d) is 
C         smaller than q(s), go to the best point in the boundary

          else 

              qamaxn = q + 0.5d0 * amaxn ** 2 * dtw + amaxn * dtr

              if ( nearlyq .and. 
     +             ( qamax .lt. q .or. qamaxn .lt. q ) ) then

                  if ( qamax .lt. qamaxn ) then
                      alpha = amax
                      q     = qamax
                  else
                      alpha = amaxn
                      q     = qamaxn
                  end if

C         Else, stop at the current point

              else

                  inform = 7

                  if ( detail_trace ) then
                      write(*,*) '[gencan-fortran-cg] neg-curv-stop'
                      write(*,*) ' call=', cg_call_id, ' iter=', iter,
     +                 ' q=', q, ' qamax=', qamax, ' qamaxn=', qamaxn
                      write(*,*) ' dtr=', dtr, ' dtw=', dtw
                  end if

                  if ( iprint .ge. 4 ) then
                      write(*, 997) inform
                      write(10,997) inform
                  end if
  
                  go to 500

              end if

          end if
      end if

      if ( detail_trace ) then
          write(*,*) '[gencan-fortran-cg] step call=', cg_call_id,
     +     ' iter=', iter, ' alpha=', alpha, ' q=', q
      end if

C     ==================================================================
C     Compute new s
C     ==================================================================

      do i = 1,nind
          sprev(i) = s(i)
          s(i)     = s(i) + alpha * d(i)
      end do

      snorm2prev = snorm2
      snorm2     = snorm2 + alpha ** 2 * dnorm2 + 2.0d0 * alpha * dts

C     ==================================================================
C     Compute the residual r = H s + g
C     ==================================================================

      rnorm2prev = rnorm2

      do i = 1,nind
          r(i) = r(i) + alpha * w(i)
      end do

      rnorm2 = norm2s(nind,r)

C     ==================================================================
C     Increment number of iterations
C     ==================================================================

      iter = iter + 1

      if ( detail_trace ) then
          write(*,*) '[gencan-fortran-cg] post call=', cg_call_id,
     +     ' iter=', iter, ' snorm2=', snorm2, ' rnorm2=', rnorm2
      end if

C     ==================================================================
C     Print information of this iteration
C     ==================================================================

      if ( iprint .ge. 4 ) then
          write(*, 984) iter,sqrt(rnorm2),sqrt(snorm2),q
          write(10,984) iter,sqrt(rnorm2),sqrt(snorm2),q
      end if

C     ==================================================================
C     Test other stopping criteria
C     ==================================================================

C     Test angle condition

      gts = 0.0d0
      do i = 1,nind
          gts = gts + g(i) * s(i)
      end do

      if ( gts .gt. 0.0d0 .or. 
     +     gts ** 2 .lt. theta ** 2 * gnorm2 * snorm2 ) then

          do i = 1,nind
              s(i) = sprev(i)
          end do

          snorm2 = snorm2prev

          q = qprev

          inform = 3

          if ( iprint .ge. 4 ) then
              write(*, 993) inform
              write(10,993) inform
          end if

          go to 500

      end if

C     If we are in the boundary of the box also stop

      if ( alpha .eq. amax2 .or. alpha .eq. amax2n ) then

          if ( alpha .eq. amax2 ) then
              rbdind  = rbdposaind
              rbdtype = rbdposatype
          else ! if (alpha.eq.amax2n) then
              rbdind  = rbdnegaind
              rbdtype = rbdnegatype
          end if

          if ( rbdtype .eq. 1 ) then
              rbdtypea = 'lower'
          else ! if (rbdtype.eq.2) then
              rbdtypea = 'upper'
          end if
 
          inform = 2

          if ( iprint .ge. 4 ) then
              write(*, 992) inform,ind(rbdind),rbdtypea
              write(10,992) inform,ind(rbdind),rbdtypea
          end if
  
          go to 500

      end if

C     If we are in the boundary of the trust region then stop

      if ( alpha .eq. amax1 .or. alpha .eq. amax1n ) then

          inform = 1

          if ( iprint .ge. 4 ) then
              write(*, 991) inform
              write(10,991) inform
          end if
  
          go to 500

      end if

C     If two consecutive iterates are much close then stop

      samep = .true.
      do i = 1,nind
         if ( abs( alpha * d(i) ) .gt. 
     +        max( epsrel * abs( s(i) ), epsabs ) ) then
              samep = .false.
          end if
      end do

      if ( samep ) then

          inform = 6

          if ( iprint .ge. 4 ) then
              write(*, 996) inform
              write(10,996) inform
          end if
  
          go to 500

      end if

C     Test whether we performed many iterations without good progress of
C     the quadratic model

C     if (abs( q - qprev ) .le. max( epsrel * abs( qprev ), epsabs ) ) 
C    +then

C         itnqmp = itnqmp + 1

C         if ( itnqmp .ge. maxitnqmp ) then

C             inform = 4

C             if ( iprint .ge. 4 ) then
C                 write(*,994)  inform,itnqmp
C                 write(10,994) inform,itnqmp
C             end if
  
C             go to 500

C         endif

C     else
C         itnqmp= 0
C     endif

C     Test whether we performed many iterations without good progress of
C     the quadratic model 

      currprog = qprev - q
      bestprog = max( currprog, bestprog )

      if ( currprog .le. epsnqmp * bestprog ) then

          itnqmp = itnqmp + 1

          if ( itnqmp .ge. maxitnqmp ) then
              inform = 4

              if ( iprint .ge. 4 ) then
                  write(*, 994) inform,itnqmp,epsnqmp,bestprog
                  write(10,994) inform,itnqmp,epsnqmp,bestprog
              end if

              go to 500
          endif

      else
          itnqmp = 0
      endif

C     ==================================================================
C     Iterate
C     ==================================================================

      go to 100

C     ==================================================================
C     End of main loop
C     ==================================================================

C     ==================================================================
C     Return
C     ==================================================================

 500  continue

      if ( detail_trace ) then
          write(*,*) '[gencan-fortran-cg] exit call=', cg_call_id,
     +     ' inform=', inform, ' iter=', iter, ' snorm2=', snorm2
      end if

C     Print final information

      if ( iprint .ge. 4 ) then
          write(*, 985) min0(nind,ncomp),(s(i),i=1,min0(nind,ncomp))
          write(10,985) min0(nind,ncomp),(s(i),i=1,min0(nind,ncomp))
      end if

      return

C     Non-executable statements

 980  format(/,6x,'Conjugate gradients (maxit= ',I7,' acc= ',1PD11.4,
     *')')
 981  format(6x,'Using Euclidian trust region (delta= ',1PD11.4,
     *')')
 982  format(6x,'Using sup-norm trust region (delta= ',1PD11.4,')')
 983  format(6x,'Unknown trust-region type')
 984  format(6x,'CG iter= ',I5,' rnorm: ',1PD11.4,' snorm= ',1PD11.4,
     *' q= ',1PD11.4)
 985  format(/,6x,'Truncated Newton direction (first ',I6, 
     *' components): ',/,1(6x,6(1PD11.4,1x)))
 990  format(6x,'Flag of CG = ',I3,' (Convergence with small residual)')
 991  format(6x,'Flag of CG = ',I3,
     *' (Convergence to the trust region boundary)')
 992  format(6x,'Flag of CG = ',I3,
     *' (Convergence to the boundary of the box constraints,',/,6x,
     *'taking step >= 1, variable ',I6,' will reaches its ',A5,
     *' bound)')
 993  format(6x,'Flag of CG = ',I3,
     *' (The next CG iterate will not satisfy the angle condition)')
 994  format(6x,'Flag of CG = ',I3,
     *' (Not enough progress in the quadratic model. This means',/,6x,
     *'that the progress of the last ',I7,' iterations was smaller ', 
     *'than ',/,6x,1PD11.4,' times the best progress (',1PD11.4,')')
 996  format(6x,'Flag of CG = ',I3,
     *' (Very near consecutive iterates)')
 997  format(6x,'Flag of CG= ',I3,
     *' (d such that d^T H d = 0 and g^T d = 0 was found)')
 998  format(6x,'Flag of CG = ',I3,' (Too many GC iterations)')
 1000 format(6x,'Flag of CG = ',I3,' Fatal Error')

      end

C     *****************************************************************
C     *****************************************************************
