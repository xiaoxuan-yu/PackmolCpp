      subroutine gencf(n,x,l,u,m,lambda,rho,epsgpen,epsgpsn,maxitnfp,
     +epsnfp,maxitngp,fmin,maxit,maxfc,udelta0,ucgmaxit,cgscre,cggpnf,
     +cgepsi,cgepsf,epsnqmp,maxitnqmp,nearlyq,nint,next,mininterp,
     +maxextrap,gtype,htvtype,trtype,iprint,ncomp,f,g,gpeucn2,gpsupn,
     +iter,fcnt,gcnt,cgcnt,spgiter,spgfcnt,tniter,tnfcnt,tnstpcnt,
     +tnintcnt,tnexgcnt,tnexbcnt,tnintfe,tnexgfe,tnexbfe,inform,s,y,d,
     +ind,lastgpns,w,eta,delmin,lspgma,lspgmi,theta,gamma,beta,sigma1,
     +sigma2,sterel,steabs,epsrel,epsabs,infrel,infabs)

      use compute_data, only : fdist, frest
      implicit none

C     SCALAR ARGUMENTS
      logical nearlyq
      integer cgcnt,cgscre,fcnt,gcnt,gtype,htvtype,inform,iprint,iter,m,
     +        maxextrap,maxfc,maxit,maxitnfp,maxitngp,maxitnqmp,
     +        mininterp,n,ncomp,spgfcnt,spgiter,tnexbcnt,tnexbfe,
     +        tnexgcnt,tnexgfe,tnfcnt,tnintcnt,tnintfe,tniter,tnstpcnt,
     +        trtype,ucgmaxit    
      double precision beta,cgepsf,cgepsi,cggpnf,delmin,epsabs,epsgpen,
     +        epsgpsn,epsnfp,epsnqmp,epsrel,eta,f,fmin,gamma,gpeucn2,
     +        gpsupn,infabs,infrel,lspgma,lspgmi,next,nint,sigma1,
     +        sigma2,steabs,sterel,theta,udelta0

C     ARRAY ARGUMENTS
      integer ind(n)
      double precision d(n),g(n),l(n),lambda(m),lastgpns(0:maxitngp-1),
     +        rho(m),s(n),u(n),w(5*n),x(n),y(n)

C     Solves the box-constrained minimization problem
C
C                         Minimize f(x)
C
C                         subject to 
C 
C                                  l <= x <= u
C     
C     using a method described in 
C
C     E. G. Birgin and J. M. Martinez, ''Large-scale active-set box-
C     constrained optimization method with spectral projected 
C     gradients'', Computational Optimization and Applications 23, pp. 
C     101-125, 2002.  
C
C     Subroutine evalal must be supplied by the user to evaluate the 
C     objective function. The prototype of evalal subroutine must be
C
C           subroutine evalal(n,x,m,lambda,rho,f,flag)
C
C     C     On Entry:
C     C
C     C     n     integer
C     C           number of variables
C     C
C     C     x     double precision x(n)
C     C           current point
C     C
C     C     m     integer
C     C           number of constraints (equalities plus inequalities)
C     C
C     C     lambda double precision lambda(m)
C     C           current estimation of the Lagrange multipliers
C     C
C     C     rho   double precision rho(m)
C     C           penalty parameters
C     C
C     C     NOTE: arguments m, lambda and rho are useful when GENCAN is 
C     C     being used for solving the box-constrained subproblems of an 
C     C     Augmented Lagrangian framework. When GENCAN is being used 
C     C     stand-alone for solving a bound-constrained problem, these 
C     C     arguments are dummy arguments and must be ignored.
C     C
C     C     On Return
C     C
C     C     f     double precision
C     C           objective function value at x
C     C
C     C     flag  integer
C     C           0 means ''no errors''
C     C           any other value means ''there was an error in the 
C     C           objective function calculation''.
C     C
C     C     SCALAR ARGUMENTS
C           integer flag,m,n
C           double precision f
C     
C     C     ARRAY ARGUMENTS
C           double precision lambda(m),rho(m),x(n)
C
C     C     ''Here it should be the body of evalal subroutine that saves 
C     C     in f the objective function value at x. Moreover, it sets 
C     C     flag equal to 0 if the calculation was successfully done and 
C     C     sets flag equal to any other value different from 0 if the 
C     C     objective function is not well defined at the current point 
C     C     x.''
C  
C           end
C
C     Subroutine evalnal to calculate the gradient of the objective 
C     function may be supplied by the user or not, depending on the 
C     value of gtype argument (gtype equal to 0 means that the evalnal 
C     subroutine will be supplied by the user and gtype equal to 1 means 
C     that an internal GENCAN subroutine will be used to estimate the 
C     gradient vector by central finite differences). In any case, a 
C     subroutine named evalnal with the following prototype must 
C     present.
C
C           subroutine evalnal(n,x,m,lambda,rho,g,flag)
C
C     C     On Entry:
C     
C     C     n     integer
C     C           number of variables
C     C
C     C     x     double precision x(n)
C     C           current point
C     C
C     C     m     integer
C     C           number of constraints (equalities plus inequalities)
C     C
C     C     lambda double precision lambda(m)
C     C           current estimation of the Lagrange multipliers
C     C
C     C     rho   double precision rho(m)
C     C           penalty parameters
C     C
C     C     NOTE: arguments m, lambda and rho are useful when GENCAN is 
C     C     being used for solving the box-constrained subproblems of an 
C     C     Augmented Lagrangian framework. When GENCAN is being used 
C     C     stand-alone for solving a bound-constrained problem, these 
C     C     arguments are dummy arguments and must be ignored.
C     C
C     C     On Return
C     C
C     C     g     double precision g(n)
C     C           gradient of the objective function at x
C     C
C     C     flag  integer
C     C           0 means ''no errors'',
C     C           any other value means ''there was an error in the 
C     C           gradient calculation''.
C     C
C     C     SCALAR ARGUMENTS
C           integer flag,m,n
C     
C     C     ARRAY ARGUMENTS
C           double precision g(n),lambda(m),rho(m),x(n)
C
C     C     ''Here it should be the body of evalnal subroutine that 
C     C     saves in g the gradient vector of the objective function at 
C     C     x. Moreover, it sets flag equal to 0 if the calculation was 
C     C     successfully done and sets flag equal to any other value 
C     C     different from 0 if the gradient vector is not well defined 
C     C     at the current point x. If GENCAN gtype argument was setted 
C     C     to 1, i.e., the finite difference approximation provided by 
C     C     GENCAN will be used, then this subroutine must even be 
C     C     present for compilation purpose but it will never be 
C     C     called.''
C  
C      end
C
C     Subroutine evalhd to calculate of the Hessian of the objective 
C     function times a given vector may be supplied by the user or not, 
C     depending on the value of htvtype argument (htvtype equal to 0 
C     means that the evalhd subroutine will be supplied by the user and 
C     htvtype equal to 1 means tha an internal GENCAN subroutine will be
C     used to estimate the product by incremental quotients). In any 
C     case, a subroutine named evalhd with the following prototype must 
C     present.
C
C           subroutine evalhd(nind,ind,n,x,m,lambda,rho,d,hd,flag)
C
C     C     On Entry:
C     C
C     C     nind  integer
C     C           number of component of the Hessian-vector product that
C     C           must be computed
C     C
C     C     ind   integer ind(nind)
C     C           the component that must be computed are ind(1)-th ... 
C     C           ind(nind)-th
C     C
C     C     n     integer
C     C           number of variables
C     C
C     C     x     double precision x(n)
C     C           current point
C     C
C     C     m     integer
C     C           number of constraints (equalities plus inequalities)
C     C
C     C     lambda double precision lambda(m)
C     C           current estimation of the Lagrange multipliers
C     C
C     C     rho   double precision rho(m)
C     C           penalty parameters
C     C
C     C     NOTE: arguments m, lambda and rho are useful when GENCAN is 
C     C     being used for solving the box-constrained subproblems of an 
C     C     Augmented Lagrangian framework. When GENCAN is being used 
C     C     stand-alone for solving a bound-constrained problem, these 
C     C     arguments are dummy arguments and must be ignored.
C     C
C     C     d     double precision d(n)
C     C           vector of the Hessian-vector product
C     C
C     C     On Return
C     C
C     C     hd    double precision g(n)
C     C           Hessian-vector product
C     C
C     C     flag  integer
C     C           0 means ''no errors'',
C     C           any other value means ''there was an error in the 
C     C           product calculation''. Just as an example, as it has
C     C           no sense that an error occurs in a matrix-vector
C     C           product, the error could happen in the Hessian
C     C           calculation. But the possible errors will depend
C     C           on the way this Hessian-vector product is computed
C     C           or approximated.
C
C     C     SCALAR ARGUMENTS
C           integer flag,m,n,nind
C     
C     C     ARRAY ARGUMENTS
C           integer ind(nind)
C           double precision d(n),hd(n),lambda(m),rho(m),x(n)
C     
C     C     ''Here it should be the body of evalhd subroutine that saves 
C     C     in hd the product of the Hessian of the objective function 
C     C     times vector d. Moreover, it sets flag equal to 0 if the 
C     C     calculation was successfully done and sets flag equal to any 
C     C     other value different from 0 if the Hessian matrix is not 
C     C     well defined at the current point x. If GENCAN htvtype 
C     C     argument was setted to 1, i.e., the incremental quotients 
C     C     approximation provided by GENCAN will be used, then this 
C     C     subroutine must even be present for compilation purposes 
C     C     but it will never be called.''
C  
C           end
C
C     In evalhd subroutine, the information about the matrix H must be 
C     passed by means of common declarations. This subroutine must be 
C     coded by the user, taking into account that only nind components 
C     of d are nonnull and that ind is the set of indices of those 
C     components. In other words, the user must write evalhd in such a 
C     way that hd is the vector whose i-th entry is
C 
C               hd(i) = \Sum_{j=1}^{nind} H_{i,ind(j)} d_ind(j)
C
C     Moreover, the only components of hd that must be computed are 
C     those which correspond to the indices ind(1),...,ind(nind). 
C     However, observe that it must be assumed that, in d, the whole 
C     dense vector is present, with its n components, even the null 
C     ones. So, if the user decides to code evalhd without taking into 
C     account the presence of ind and nind, it can be easily done. A 
C     final observation: probably, if nind is close to n, it is not 
C     worthwhile to use ind, due to the cost of accessing the correct 
C     indices. 
C
C     Example: Assume that H is dense. The main steps of evalhd could 
C     be:
C
C          do i = 1,nind
C              indi     = ind(i)
C              hd(indi) = 0.0d0
C              do j = 1,nind
C                  indj     = ind(j)
C                  hd(indi) = hd(indi) + H(indi,indj) * d(indj)
C              end do
C          end do
C
C
C     Description of the GENCAN arguments:
C
C     On Entry
C
C     n        integer 
C              number of variables
C
C     x        double precision x(n)
C              initial estimation of the solution
C
C     l        double precision l(n)
C              lower bounds on the variables
C
C     u        double precision u(n)
C              upper bounds on the variables
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
C     epsgpen  double precision
C              epsgpen means EPSilon for the Projected Gradient Euclidian
C              Norm. It is a small positive number for declaring 
C              convergence when the Euclidian norm of the continuous 
C              projected gradient is less than or equal to epsgpen
C
C              RECOMMENDED: epsgpen = 1.0d-05
C
C              CONSTRAINTS: epsgpen >= 0.0
C
C     epsgpsn  double precision
C              epsgpsn means EPSilon for the Projected Gradient Sup Norm.
C              It is a small positive number for declaring convergence 
C              when the sup norm of the continuous projected gradient is 
C              less than or equal to epsgpsn
C
C              RECOMMENDED: epsgpsn = 1.0d-05
C
C              CONSTRAINTS: epsgpsn >= 0.0
C
C     maxitnfp integer
C              maxitnfp means MAXimum of ITerations with No Function 
C              Progress. See below for more details.
C
C     epsnfp   double precision
C              epsnfp means EPSilon for No Function Progress. It is a
C              small positive number for declaring ''lack of progress in 
C              the objective function value'' if f(x_k) - f(x_{k+1}) <= 
C              epsnfp * max{ f(x_j) - f(x_{j+1}, j < k } during maxitnfp 
C              consecutive iterations. This stopping criterion may be 
C              inhibited setting maxitnfp equal to maxit.
C
C              RECOMMENDED: maxitnfp = 5 and epsnfp = 1.0d-02
C
C              CONSTRAINTS: maxitnfp >= 1 and epsnfp >= 0.0
C
C     maxitngp integer
C              maxitngp means MAXimum of ITerations with No Gradient
C              Progress. If the order of the Euclidian norm of the 
C              continuous projected gradient did not change during 
C              maxitngp consecutive iterations then the execution stops. 
C
C              RECOMMENDED: maxitngp = 10
C
C              CONSTRAINTS: maxitngp >= 1
C
C     fmin     double precision
C              function value for the stopping criteria f <= fmin
C
C              There is a stopping criterion that stops GENCAN if a 
C              point with a functional value smaller than fmin is found. 
C              The idea behind this stopping criterion is to stop the 
C              method if the objective function is not bounded from 
C              below.
C
C              RECOMMENDED: fmin = - infabs
C
C              CONSTRAINTS: there are no constraints for this argument
C
C     maxit    integer
C              maximum number of allowed iterations
C
C              RECOMMENDED: maxit = 1000
C
C              CONSTRAINTS: maxit >= 0
C
C     maxfc    integer
C              maximum allowed number of functional evaluations
C
C              RECOMMENDED: maxfc = 5 * maxit
C
C              CONSTRAINTS: maxfc >= 1
C
C     udelta0  double precision
C              initial ''trust-radius'' for Conjugate Gradients. The 
C              default value max( delmin, 0.1 * max( 1, ||x|| ) ) is 
C              used if the user sets udelta0 <= 0. 
C
C              RECOMMENDED: udelta0 = - 1.0
C
C              CONSTRAINTS: there are no constraints for this argument
C
C     ucgmaxit integer
C              maximum allowed number of iterations for each run of the 
C              Conjugate Gradient subalgorithm
C
C              The default values for this argument is max( 1, 10 * 
C              log( nind ) ), where nind is the number of free 
C              variables, and it will be used if the user sets ucgmaxit 
C              to any non-positive value. 
C
C          RECOMMENDED: ucgmaxit = - 1
C
C          CONSTRAINTS: there are no constraints for this argument
C
C     cgscre   integer
C              See below
C
C     cggpnf   double precision
C              cgscre means conjugate gradient stopping criterion 
C              relation, and cggpnf means Conjugate Gradients projected 
C              gradient final norm. Both are related to a stopping 
C              criterion of Conjugate Gradients. This stopping criterion 
C              depends on the norm of the residual of the linear system. 
C              The norm of the residual should be less or equal than a 
C              ''small'' quantity which decreases as we are 
C              approximating the solution of the minimization problem 
C              (near the solution, better the truncated-Newton direction 
C              we aim). Then, the log of the required accuracy requested 
C              to Conjugate Gradient has a linear dependence on the log 
C              of the norm of the continuous projected gradient. This 
C              linear relation uses the squared Euclidian norm of the 
C              projected gradient if cgscre is equal to 1 and uses the 
C              sup-norm if cgscre is equal to 2. In addition, the 
C              precision required to CG is equal to cgepsi (conjugate 
C              gradient initial epsilon) at x0 and cgepsf (conjugate 
C              gradient final epsilon) when the Euclidian- or sup-norm 
C              of the projected gradient is equal to cggpnf (conjugate 
C              gradients projected gradient final norm) which is an 
C              estimation of the value of the Euclidian- or sup-norm of 
C              the projected gradient at the solution.
C
C              RECOMMENDED: cgscre = 1, cggpnf = epsgpen; or
C                           cgscre = 2, cggpnf = epsgpsn.
C
C              CONSTRAINTS:  allowed values for cgscre are just 1 or 2
C                            cggpnf >= 0.0
C
C     cgepsi   double precision
C              See below
C
C     cgepsf   double precision
C              small positive numbers for declaring convergence of the 
C              Conjugate Gradients subalgorithm when ||r||_2 < cgeps * 
C              ||rhs||_2, where r is the residual and rhs is the right 
C              hand side of the linear system, i.e., CG stops when the 
C              relative error of the solution is smaller than cgeps. 
C
C              cgeps varies from cgepsi to cgepsf in a way that depends 
C              on cgscre as follows:
C
C              i) CASE cgscre = 1: log10(cgeps^2) depends linearly on 
C              log10(||g_P(x)||_2^2) which varies from ||g_P(x_0)||_2^2 
C              to epsgpen^2
C
C              ii)  CASE cgscre = 2: log10(cgeps) depends linearly on 
C              log10(||g_P(x)||_inf) which varies from ||g_P(x_0)||_inf 
C              to epsgpsn
C
C              RECOMMENDED: cgepsi = 1.0d-01, cgepsf = 1.0d-05
C
C              CONSTRAINTS: cgepsi >= cgepsf >= 0.0
C
C     epsnqmp  double precision
C              See below
C
C     maxitnqmp integer
C              This and the previous argument are used for a stopping 
C              criterion of the Conjugate Gradients subalgorithm. If the 
C              progress in the quadratic model is smaller than fraction 
C              of the best progress ( epsnqmp * bestprog ) during 
C              maxitnqmp consecutive iterations then CG is stopped 
C              declaring ''not enough progress of the quadratic model''.
C
C              RECOMMENDED: epsnqmp = 1.0d-04, maxitnqmp = 5
C
C              CONSTRAINTS: epsnqmp >= 0.0, maxitnqmp >= 1.
C
C     nearlyq  logical
C              If the objective function is (nearly) quadratic, use the 
C              option nearlyq = TRUE. Otherwise, keep the default 
C              option.
C
C              If, in an iteration of CG we find a direction d such that 
C              d^T H d <= 0 then we take the following decision:
C
C              (i) If nearlyq = TRUE then we take direction d and try to 
C              go to the boundary choosing the best point among the two 
C              points at the boundary and the current point. 
C
C              (ii) If nearlyq = FALSE then we stop at the current point.
C
C              Moreover, if the objective function is quadratic more 
c              effort is due in computing the Truncated Newton direction.
C
C              RECOMMENDED: nearlyq = FALSE
C
C              CONSTRAINTS: allowed values are just TRUE or FALSE.
C
C     nint     double precision
C              Constant for the interpolation. See the description of 
C              sigma1 and sigma2 above. Sometimes, in a line search, we 
C              take the new trial step as the previous one divided by 
C              nint
C
C              RECOMMENDED: nint = 2.0
C
C              CONSTRAINTS: nint > 1.0.
C
C     next     double precision
C              Constant for the extrapolation. When extrapolating we 
C              try alpha_new = alpha * next
C
C              RECOMMENDED: next = 2.0
C
C              CONSTRAINTS: next > 1.0
C
C     mininterp integer
C              Constant for testing if, after having made at least 
C              mininterp interpolations, the steplength is too small. In
C              that case, failure of the line search is declared (may be 
C              the direction is not a descent direction due to an error 
C              in the gradient calculations). Use mininterp greater 
C              than or equal to maxfc for inhibit this stopping 
C              criterion
C
C              RECOMMENDED: mininterp = 4 
C
C              CONSTRAINTS: mininterp >= 1
C
C     maxextrap integer
C              Constant to limit the number of extrapolations in the 
C              Truncated Newton direction.
C
C              RECOMMENDED: maxextrap = 100 
C
C              CONSTRAINTS: maxextrap >= 0
C
C     gtype    integer
C              gtype indicates in which way the gradient of the 
C              objective function will be computed. If the user have 
C              been implemented the user-supplied evalnal subroutine to 
C              compute the gradient of the objective function then 
C              gtype argument must be set to 0 (ZERO) and the user-
C              supplied evalnal subroutine will be called by GENCAN any 
C              time the gradient would be required.
C
C                    subroutine evalnal(n,x,m,lambda,rho,g,flag)
C
C              C     On Entry:
C     
C              C     n     integer,
C              C           number of variables,
C              C
C              C     x     double precision x(n),
C              C           current point,
C              C
C              C     m     integer,
C              C           number of constraints (equalities plus 
C              C           inequalities),
C              C
C              C     lambda double precision lambda(m),
C              C           current estimation of the Lagrange 
C              C           multipliers,
C              C
C              C     rho   double precision rho(m)
C              C           penalty parameters,
C              C
C              C     NOTE: arguments m, lambda and rho are useful when 
C              C     GENCAN is being used for solving the box-
C              C     constrained subproblems of an Augmented Lagrangian 
C              C     framework. When GENCAN is being used stand-alone 
C              C     for solving a bound-constrained problem, these 
C              C     arguments are dummy arguments.
C              C
C              C     On Return
C              C
C              C     g     double precision g(n),
C              C           gradient of the objective function at x,
C              C
C              C     flag  integer
C              C           0 means ''no errors'',
C              C           1 means ''some error occurs in the gradient 
C              C             evaluation''.
C              C
C              C     SCALAR ARGUMENTS
C                    integer flag,m,n
C     
C              C     ARRAY ARGUMENTS
C                    double precision g(n),lambda(m),rho(m),x(n)
C
C              C     ''Here it should be the body of evalnal subroutine 
C              C     that saves in g the gradient vector of the 
C              C     objective at x. Moreover, it sets flag equal to 0 
C              C     if the calculation was successfully done and sets 
C              C     flag equal to any other value different from 0 if 
C              C     the gradient vector is not well defined at the 
C              C     current point x. If GENCAN gtype argument was 
C              C     setted to 1, i.e., the finite difference 
C              C     approximation provided by GENCAN will be used, then 
C              C     this subroutine must even be present for 
C              C     compilation purposes but it will never be called.''
C  
C               end
C
C              If, on the other hand, the user is not able to provide 
C              evalnal subroutine, gtype argument must be set to 1 
C              (ONE). In this case, every time GENCAN needs to compute 
C              the gradient of the objective function, an internal 
C              subroutine that approximates it by finite-differences 
C              will be used (be aware that it maybe very time 
C              consuming). Moreover, note that the evalnal subroutine 
C              must still be present (with an empty body).
C
C              RECOMMENDED: gtype = 0 (provided you have the evalg 
C                           subroutine)
C
C              CONSTRAINTS: allowed values are just 0 or 1.
C
C     htvtype  integer
C              htvtype indicates in which way the product of the Hessian 
C              of the objective function times an arbitrary vector will be 
C              computed. If the user has not been implemented the user-
C              supplied evalhd subroutine to do this task then htvtype 
C              argument must be set to 1 (ONE). In this case an internal 
C              subroutine that approximates this product by incremental 
C              quotients will be used. Note that, even in this case, 
C              evalhd subroutine must be present (with an empty body). 
C              This is the default option and the empty-body subroutine 
C              follows:
C
C              subroutine evalhd(nind,ind,n,x,m,lambda,rho,d,hd,flag)
C
C              C     SCALAR ARGUMENTS
C                    integer nind,n,m,flag
C
C              C     ARRAY ARGUMENTS
C                    integer ind(nind)
C                    double precision x(n),lambda(m),rho(m),d(n),hd(n) 
C
C                    flag = - 1
C
C                    end
C
C              If, on the other hand, the user prefers to implement his/
C              her own evalhd subroutine then htvtype argument must be 
C              set to 0 (ZERO). In this case, the product of the Hessian 
C              times vector d (input argument of evalhd subroutine) must 
C              be saved in vector hd (output argument of evalhd 
C              subroutine). The other arguments description as well as 
C              some hints on how to implement your own evalhd subroutine 
C              can be found in the GENCAN arguments description.
C
C              RECOMMENDED: htvtype = 1
C
C              (you take some risk using this option but, unless you 
C              have a good evalhd subroutine, incremental quotients is a 
C              very cheap option)
C
C              CONSTRAINTS: allowed values are just 0 or 1.
C
C     trtype   integer
C              Type of Conjugate Gradients ''trust-radius''. trtype 
C              equal to 0 means Euclidian-norm trust-radius and trtype 
C              equal to 1 means sup-norm trust radius
C
C              RECOMMENDED: trtype = 0
C
C              CONSTRAINTS: allowed values are just 0 or 1.
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
C     s        double precision s(n)
C     y        double precision y(n)
C     d        double precision d(n)
C     ind      integer ind(n)
C     lastgpns double precision lastgpns(maxitngp)
C     w        double precision w(5*n)
C              working vectors
C
C     eta      double precision
C              Constant for deciding abandon the current face or not. We 
C              abandon the current face if the norm of the internal 
C              gradient (here, internal components of the continuous 
C              projected gradient) is smaller than ( 1 - eta ) times the 
C              norm of the continuous projected gradient. Using eta = 
C              0.9 is a rather conservative strategy in the sense that 
C              internal iterations are preferred over SPG iterations. 
C
C              RECOMMENDED: eta = 0.9
C
C              CONSTRAINTS: 0.0 < eta < 1.0
C
C     delmin   double precision
C              Smaller Conjugate Gradients ''trust radius'' to compute 
C              the Truncated Newton direction
C
C              RECOMMENDED: delmin = 0.1
C
C              CONSTRAINTS: delmin > 0.0
C
C     lspgmi   double precision
C              See below
C
C     lspgma   double precision
C              The spectral steplength, called lamspg, is projected onto 
C              the box [lspgmi,lspgma] 
C
C              RECOMMENDED: lspgmi = 1.0d-10 and lspgma = 1.0d+10
C 
C              CONSTRAINTS: lspgma >= lspgmi > 0.0
C
C     theta    double precision
C              Constant for the angle condition, i.e., at iteration k we 
C              need a direction dk such that <gk,dk> <= - theta 
C              ||gk||_2 ||dk||_2, where gk is \nabla f(xk)
C
C              RECOMMENDED: theta = 10^{-6}
C
C              CONSTRAINTS: 0.0 < theta < 1.0
C
C     gamma    double precision
C              Constant for the Armijo criterion
C              f(x + alpha d) <= f(x) + gamma * alpha * <g,d>
C
C              RECOMMENDED: gamma = 1.0d-04
C
C              CONSTRAINTS: 0.0 < gamma < 0.5.
C
C     beta     double precision
C              Constant for the beta condition <dk, g(xk + dk)>  < beta 
C              * <dk,gk>. If (xk + dk) satisfies the Armijo condition 
C              but does not satisfy the beta condition then the point is 
C              accepted, but if it satisfied the Armijo condition and 
C              also satisfies the beta condition then we know that there 
C              is the possibility for a successful extrapolation
C
C              RECOMMENDED: beta = 0.5
C
C              CONSTRAINTS: 0.0 < beta < 1.0.
C
C     sigma1   double precision
C              See below
C
C     sigma2   double precision
C              Constant for the safeguarded interpolation. If alpha_new 
C              is not inside the interval [sigma1, sigma * alpha] then 
C              we take alpha_new = alpha / nint
C
C              RECOMMENDED: sigma1 = 0.1 and sigma2 = 0.9
C
C              CONSTRAINTS: 0 < sigma1 < sigma2 < 1.
C
C     sterel   double precision
C              See below
C
C     steabs   double precision
C              This constants mean a ''relative small number'' and ''an 
C              absolute small number'' for the increments in finite 
C              difference approximations of derivatives
C
C              RECOMMENDED: epsrel = 1.0d-07 and epsabs = 1.0d-10 
C
C              CONSTRAINTS: sterel >= steabs > 0
C
C     epsrel   double precision
C              See below
C
C     epsabs   double precision
C              See below
C
C     infrel   double precision
C              See below
C
C     infabs   double precision
C              This four constants mean a ''relative small number'', 
C              ''an absolute small number'', ''a relative large number'' 
C              and ''an absolute large number''. Basically, a quantity A 
C              is considered negligible with respect to another quantity 
C              B if |A| < max ( epsrel * |B|, epsabs ) 
C
C              RECOMMENDED: epsrel = 1.0d-10, epsabs = 1.0d-20, 
C                           infrel = 1.0d+20, infabs = 1.0d+99
C
C              CONSTRAINTS: epsrel >= epsabs >= 0.0
C                           infabs >= infrel >= 0.0
C
C     On Return
C
C     x        double precision x(n)
C              Final estimation to the solution
C
C     f        double precision
C              Function value at the final estimation 
C
C     g        double precision g(n)
C              Gradient at the final estimation
C
C     gpeucn2  double precision
C              Squared Euclidian norm of the continuous projected 
C              gradient at the final estimation
C
C     gpsupn   double precision
C              the same as before but with sup-norm
C
C     iter     integer
C              number of iterations
C
C     fcnt     integer
C              number of function evaluations   
C
C     gcnt     integer
C              number of gradient evaluations   
C
C     cgcnt    integer
C              number of Conjugate Gradients iterations   
C
C     spgiter  integer
C              number of Spectral Projected Gradient iterations
C
C     spgfcnt  integer
C              number of functional evaluations along Spectral Projected
C              Gradient directions
C
C     tniter   integer
C              number of Truncated-Newton iterations
C
C     tnfcnt   integer
C              number of functional evaluations along Truncated-Newton
C              directions
C
C     tnintcnt integer
C              number of times a backtracking in a Truncated-Newton
C              direction was needed
C
C     tnexgcnt integer
C              number of times an extrapolation in a Truncated-Newton
C              direction successfully decreased the objective funtional
C              value
C
C     tnexbcnt integer
C              number of times an extrapolation was aborted in the first
C              extrapolated point by an increase in the objective 
C              functional value
C
C     tnstpcnt integer
C              number of times the Newton point was accepted (without
C              interpolations nor extrapolations)
C
C     tnintfe  integer
C              number of functional evaluations used in interpolations 
C              along Truncated-Newton directions
C
C     tnexgfe  integer
C              number of functional evaluations used in successful 
C              extrapolations along Truncated-Newton directions
C
C     tnexbfe  integer
C              number of functional evaluations used in unsuccessful 
C              extrapolations along Truncated-Newton directions
C
C     inform   integer
C              This output parameter tells what happened in this 
C              subroutine, according to the following conventions:
C 
C              0 = convergence with small Euclidian norm of the 
C                  continuous projected gradient (smaller than epsgpen);
C
C              1 = convergence with small sup-norm of the continuous 
C                  projected gradient (smaller than epsgpsn);
C
C              2 = the algorithm stopped by ''lack of progress'', that 
C                  means that f(xk) - f(x_{k+1}) <= epsnfp * 
C                  max{ f(x_j) - f(x_{j+1}, j < k } during maxitnfp 
C                  consecutive iterations. If desired, set maxitnfp 
C                  equal to maxit to inhibit this stopping criterion.
C
C              3 = the algorithm stopped because the order of the 
C                  Euclidian norm of the continuous projected gradient 
C                  did not change during maxitngp consecutive 
C                  iterations. Probably, we are asking for an 
C                  exaggerated small norm of continuous projected 
C                  gradient for declaring convergence. If desired, set
C                  maxitngp equal to maxit to inhibit this stopping 
C                  criterion.
C
C              4 = the algorithm stopped because the functional value 
c                  is very small (smaller than fmin). If desired, set 
C                  fmin equal to minus infabs to inhibit this stopping 
C                  criterion.
C
C              6 = too small step in a line search. After having made at 
C                  least mininterp interpolations, the steplength 
C                  becames small. ''small steplength'' means that we are 
C                  at point x with direction d and step alpha, and 
C
C                  alpha * ||d||_infty < max( epsabs, epsrel * 
C                  ||x||_infty ). 
C 
C                  In that case failure of the line search is declared 
C                  (may be the direction is not a descent direction due 
C                  to an error in the gradient calculations). If 
C                  desired, set mininterp equal to maxfc to inhibit this 
C                  stopping criterion.
C
C              7 = it was achieved the maximum allowed number of 
C                  iterations (maxit);
C
C              8 = it was achieved the maximum allowed number of 
C                  function evaluations (maxfc);
C
C            < 0 = error in evalal, evalnal or evalhd subroutines.

C     LOCAL SCALARS
      character * 3 ittype
      integer cgiter,cgmaxit,fcntprev,i,infotmp,itnfp,nind,nprint,
     +        rbdind,rbdtype,tnexbprev,tnexgprev,tnintprev
      double precision acgeps,amax,amaxx,bestprog,bcgeps,cgeps,currprog,
     +        delta,epsgpen2,fprev,gieucn2,gpeucn20,gpi,gpnmax,gpsupn0,
     +        kappa,lamspg,ometa2,sts,sty,xnorm
      logical debug_trace,packmolprecision

C     ==================================================================
C     Initialization
C     ==================================================================

C     Set some initial values:

C     counters,
      iter     =  0
      fcnt     =  0
      gcnt     =  0
      cgcnt    =  0

      spgiter  =  0
      spgfcnt  =  0

      tniter   =  0
      tnfcnt   =  0

      tnstpcnt =  0
      tnintcnt =  0
      tnexgcnt =  0
      tnexbcnt =  0

      tnintfe  =  0
      tnexgfe  =  0
      tnexbfe  =  0

C     just for printing,
      nprint   = min0( n, ncomp )

C     for testing convergence,
      epsgpen2 = epsgpen ** 2

C     for testing whether to abandon the current face or not,
C     (ometa2 means '(one minus eta) squared')
      ometa2   = ( 1.0d0 - eta ) ** 2

C     for testing progress in f, and
      fprev    = infabs
      bestprog =  0.0d0
      itnfp    =      0

C     for testing progress in the projected gradient norm.
      do i = 0,maxitngp - 1
          lastgpns(i) = infabs
      end do

C     Print problem information

      if( iprint .ge. 3 ) then
          write(*, 977) n
          write(*, 978) nprint,(l(i),i=1,nprint)
          write(*, 979) nprint,(u(i),i=1,nprint)
          write(*, 980) nprint,(x(i),i=1,nprint)

          write(10,977) n
          write(10,978) nprint,(l(i),i=1,nprint)
          write(10,979) nprint,(u(i),i=1,nprint)
          write(10,980) nprint,(x(i),i=1,nprint)
      end if

C     Project initial guess. If the initial guess is infeasible, 
C     projection puts it into the box.

      do i = 1,n
          x(i) = max( l(i), min( x(i), u(i) ) )
      end do

C     Compute x Euclidian norm

      xnorm = 0.0d0
      do i = 1,n
          xnorm = xnorm + x(i) ** 2
      end do
      xnorm = sqrt( xnorm )

C     Compute function and gradient at the initial point

      call evalal(n,x,m,lambda,rho,f,inform)

c LM: Added packmolprecision function test, for Packmol

      if ( packmolprecision(n,x) ) then
        if(iprint.gt.0) then
          write(*,780)
780       format('  Current point is a solution.') 
        end if
        return
      end if

      fcnt = fcnt + 1

      if ( inform .lt. 0 ) then

          if ( iprint .ge. 3 ) then
              write(*, 1000) inform
              write(10,1000) inform
          end if

          return
      end if

      if ( gtype .eq. 0 ) then
          call evalnal(n,x,m,lambda,rho,g,inform)
      else ! if ( gtype .eq. 1 ) then
          call evalnaldiff(n,x,m,lambda,rho,g,sterel,steabs,inform)
      end if
      gcnt = gcnt + 1

      if ( inform .lt. 0 ) then

          if ( iprint .ge. 3 ) then
              write(*, 1000) inform
              write(10,1000) inform
          end if

          return
      end if

C     Compute continuous-project-gradient Euclidian and Sup norms,
C     internal gradient Euclidian norm, and store in nind the number of
C     free variables and in array ind their identifiers.

      nind    = 0
      gpsupn  = 0.0d0
      gpeucn2 = 0.0d0
      gieucn2 = 0.0d0
      do i = 1,n
          gpi     = min( u(i), max( l(i), x(i) - g(i) ) ) - x(i)
          gpsupn  = max( gpsupn, abs( gpi ) )
          gpeucn2 = gpeucn2 + gpi ** 2
          if ( x(i) .gt. l(i) .and. x(i) .lt. u(i) ) then
              gieucn2   = gieucn2 + gpi ** 2
              nind      = nind + 1
              ind(nind) = i
          end if
      end do

C     Compute a linear relation between gpeucn2 and cgeps2, i.e.,
C     scalars a and b such that 
c
C         a * log10(||g_P(x_0)||_2^2) + b = log10(cgeps_0^2) and
c
C         a * log10(||g_P(x_f)||_2^2) + b = log10(cgeps_f^2),
c
C     where cgeps_0 and cgeps_f are provided. Note that if 
C     cgeps_0 is equal to cgeps_f then cgeps will be always 
C     equal to cgeps_0 and cgeps_f.

C     We introduce now a linear relation between gpsupn and cgeps also.

c LM: changed to avoid error with gpsupn=0
      call gp_ieee_signal1(gpsupn, acgeps, bcgeps, cgepsf,
     +      cgepsi, cggpnf)
c      if ( gpsupn .gt. 0.d0 ) then
c         acgeps = log10( cgepsf / cgepsi ) / log10( cggpnf / gpsupn )
c         bcgeps = log10( cgepsi ) - acgeps * log10( gpsupn )
c      else
c         acgeps = 0.0d0
c         bcgeps = cgepsf
c      end if

c      if ( cgscre .eq. 1 ) then
c          acgeps = 2.0d0 * log10( cgepsf / cgepsi ) / 
c     +                     log10( cggpnf ** 2 / gpeucn2 )
c          bcgeps = 2.0d0 * log10( cgepsi ) - acgeps * log10( gpeucn2 )
c      else ! if ( cgscre .eq. 2 ) then
c          acgeps = log10( cgepsf / cgepsi ) / log10( cggpnf / gpsupn )
c          bcgeps = log10( cgepsi ) - acgeps * log10( gpsupn )
c      end if 

C     And it will be used for the linear relation of cgmaxit

      gpsupn0  = gpsupn
      gpeucn20 = gpeucn2

      call packmol_gencan_debug_enabled_f(debug_trace)
      if ( debug_trace ) then
          call packmol_gencan_debug_trace_iter('entry', f, epsgpsn,
     +    fdist, frest, gpsupn, gpeucn2, gieucn2, xnorm,
     +    gieucn2 .le. ometa2 * gpeucn2)
      end if

C     Print initial information

c LM: progress bar for packmol
      call printbar(iprint, iter, maxit)

      if( iprint .ge. 3 ) then
          write(*, 981) iter
          write(*, 985) nprint,(x(i),i=1,nprint)
          write(*, 986) nprint,(g(i),i=1,nprint)
          write(*, 987) nprint,(min(u(i),max(l(i),x(i)-g(i)))-x(i),i=1,
     +    nprint)
          write(*, 988) min0(nprint,nind),nind,(ind(i),i=1,min0(nprint,
     +    nind))
          write(*, 1002) f,sqrt(gpeucn2),sqrt(gieucn2),gpsupn,nind,n,
     +    spgiter,tniter,fcnt,gcnt,cgcnt

          write(10,981) iter
          write(10,985) nprint,(x(i),i=1,nprint)
          write(10,986) nprint,(g(i),i=1,nprint)
          write(10,987) nprint,(min(u(i),max(l(i),x(i)-g(i)))-x(i),i=1,
     +    nprint)
          write(10,988) min0(nprint,nind),nind,(ind(i),i=1,min0(nprint,
     +    nind))
          write(10,1002) f,sqrt(gpeucn2),sqrt(gieucn2),gpsupn,nind,n,
     +    spgiter,tniter,fcnt,gcnt,cgcnt
      end if

C     ==================================================================
C     Main loop
C     ==================================================================
      
 100  continue

      if ( debug_trace ) then
          call packmol_gencan_debug_trace_iter('loop', f, epsgpsn,
     +    fdist, frest, gpsupn, gpeucn2, gieucn2, xnorm,
     +    gieucn2 .le. ometa2 * gpeucn2)
      end if

C     ==================================================================
C     Test stopping criteria
C     ==================================================================

c LM: Added packmolprecision function test, for Packmol

      if ( packmolprecision(n,x) ) then
        goto 500
      end if

C     Test whether the continuous-projected-gradient Euclidian norm
C     is small enough to declare convergence

      if ( gpeucn2 .le. epsgpen2 ) then
          inform = 0

          if ( iprint .ge. 3 ) then
              write(*, 990) inform,epsgpen
              write(10,990) inform,epsgpen
          end if

          go to 500
      end if

C     Test whether the continuous-projected-gradient Sup norm
C     is small enough to declare convergence

      if ( gpsupn .le. epsgpsn ) then
          inform = 1

          if ( iprint .ge. 3 ) then
              write(*, 991) inform,epsgpsn
              write(10,991) inform,epsgpsn
          end if

          go to 500
      end if

C     Test whether we performed many iterations without good progress
C     of the functional value

      currprog = fprev - f
      bestprog = max( currprog, bestprog )

      if ( currprog .le. epsnfp * bestprog ) then

          itnfp = itnfp + 1

          if ( itnfp .ge. maxitnfp ) then
              inform = 2

              if ( iprint .ge. 3 ) then
                  write(*, 992) inform,epsnfp,maxitnfp
                  write(10,992) inform,epsnfp,maxitnfp
              end if

              go to 500
          endif

      else
          itnfp = 0
      endif

C     Test whether we have performed many iterations without good 
C     reduction of the euclidian-norm of the projected gradient

      gpnmax = 0.0d0
      do i = 0,maxitngp - 1
          gpnmax = max( gpnmax, lastgpns(i) )
      end do

      lastgpns(mod( iter, maxitngp )) = gpeucn2

      if ( gpeucn2 .ge. gpnmax ) then

          inform = 3

          if ( iprint .ge. 3 ) then
              write(*, 993) inform,maxitngp
              write(10,993) inform,maxitngp
          end if

          go to 500

      endif

C     Test whether the functional value is very small

      if ( f .le. fmin ) then

          inform = 4

          if ( iprint .ge. 3 ) then
              write(*, 994) inform,fmin
              write(10,994) inform,fmin
          end if

          go to 500

      end if

C     Test whether the number of iterations is exhausted

      if ( iter .ge. maxit ) then

          inform = 7

          if ( iprint .ge. 3 ) then
              write(*, 997) inform,maxit
              write(10,997) inform,maxit
          end if

          go to 500

      end if

C     Test whether the number of functional evaluations is exhausted

      if ( fcnt .ge. maxfc ) then

          inform = 8

          if ( iprint .ge. 3 ) then
              write(*, 998) inform,maxfc
              write(10,998) inform,maxfc
          end if

          go to 500

      end if

C     ==================================================================
C     The stopping criteria were not satisfied, a new iteration will be 
C     made
C     ==================================================================

      iter = iter + 1

C     ==================================================================
C     Save current values, f, x and g
C     ==================================================================

      fprev = f

      do i = 1,n
          s(i) = x(i)
          y(i) = g(i)
      end do

C     ==================================================================
C     Compute new iterate
C     ==================================================================

C     We abandon the current face if the norm of the internal gradient
C     (here, internal components of the continuous projected gradient)
C     is smaller than (1-eta) times the norm of the continuous 
C     projected gradient. Using eta=0.9 is a rather conservative 
C     strategy in the sense that internal iterations are preferred over 
C     SPG iterations. Replace eta = 0.9 by other tolerance in (0,1) if 
C     you find it convenient. 

      if ( gieucn2 .le. ometa2 * gpeucn2 ) then

C         ==============================================================
C         Some constraints should be abandoned. Compute the new iterate 
C         using an SPG iteration
C         ==============================================================

          ittype  = 'SPG'
          spgiter = spgiter + 1

C         Compute spectral steplength

          if ( iter .eq. 1 .or. sty .le. 0.0d0 ) then
              lamspg = max( 1.0d0, xnorm ) / sqrt( gpeucn2 )
          else
              lamspg = sts / sty
          end if
          lamspg = min( lspgma, max( lspgmi, lamspg ) )

C         Perform a line search with safeguarded quadratic interpolation 
C         along the direction of the spectral continuous projected 
C         gradient

          fcntprev = fcnt

          call spgls(n,x,m,lambda,rho,f,g,l,u,lamspg,nint,mininterp,
     +    fmin,maxfc,iprint,fcnt,inform,w(1),w(n+1),gamma,sigma1,sigma2,
     +    sterel,steabs,epsrel,epsabs,infrel,infabs) 

          spgfcnt = spgfcnt + ( fcnt - fcntprev ) 

          if ( inform .lt. 0 ) then

              if ( iprint .ge. 3 ) then
                  write(*, 1000) inform
                  write(10,1000) inform
              end if

              return
          end if

C         Compute the gradient at the new iterate

          if ( gtype .eq. 0 ) then
              call evalnal(n,x,m,lambda,rho,g,inform)
          else ! if ( gtype .eq. 1 ) then
              call evalnaldiff(n,x,m,lambda,rho,g,sterel,steabs,inform)
          end if
          gcnt = gcnt + 1
 
          if ( inform .lt. 0 ) then

              if ( iprint .ge. 3 ) then
                  write(*, 1000) inform
                  write(10,1000) inform
              end if

              return
          end if

      else

C         ==============================================================
C         The new iterate will belong to the closure of the current face
C         ==============================================================

          ittype = 'TN '
          tniter = tniter + 1

C         Compute trust-region radius

          if ( iter .eq. 1 ) then
              if( udelta0 .le. 0.0d0 ) then
                  delta = max( delmin, 0.1d0 * max( 1.0d0, xnorm ) )
              else
                  delta = udelta0
              end if
          else
              delta = max( delmin, 10.0d0 * sqrt( sts ) )
          end if

C         Shrink the point, its gradient and the bounds

          call shrink(nind,ind,n,x)
          call shrink(nind,ind,n,g)
          call shrink(nind,ind,n,l)
          call shrink(nind,ind,n,u)

C         Compute the descent direction solving the newtonian system by 
C         conjugate gradients

C         Set conjugate gradient stopping criteria. Default values are 
C         taken if you set ucgeps < 0 and ucgmaxit < 0, respectively. 
C         Otherwise, the parameters cgeps and cgmaxit will be the ones 
C         set by the user.

        call gp_ieee_signal2(
     +    cgmaxit, nind, nearlyq, ucgmaxit, cgscre,
     +    kappa, gpeucn2, gpeucn20, epsgpen2, epsgpsn,
     +    cgeps, acgeps, bcgeps, cgepsf, cgepsi, gpsupn, gpsupn0)

          if ( debug_trace ) then
              call packmol_gencan_debug_trace_tn(iter, delta, cgeps,
     +        cgmaxit, kappa, gpsupn, gpeucn2, xnorm)
          end if

c          if( ucgmaxit .le. 0 ) then
c              if ( nearlyq ) then
c                  cgmaxit = nind
c              else
c                  if ( cgscre .eq. 1 ) then
c                      kappa = log10( gpeucn2 / gpeucn20 )/
c     +                        log10( epsgpen2 / gpeucn20 )
c                  else ! if ( cgscre .eq. 2 ) then
c                      kappa= log10( gpsupn / gpsupn0 ) / 
c     +                       log10( epsgpsn / gpsupn0 )
c                  end if
c                  kappa = max( 0.0d0, min( 1.0d0, kappa ) )
c                  cgmaxit = int(
c     +            ( 1.0d0 - kappa ) * max( 1.0d0, 10.0d0 * 
c     +            log10( dfloat( nind ) ) ) + kappa * dfloat( nind ) )
cc L. Martinez added to accelerate the iterations near the solution 
c                  cgmaxit = min(20,cgmaxit)
c              end if
cc              cgmaxit = 2 * nind
c          else
c              cgmaxit = ucgmaxit
c          end if
c
c          if ( cgscre .eq. 1 ) then
c              cgeps = sqrt( 10.0d0 ** ( acgeps * log10( gpeucn2 ) + 
c     +                bcgeps ) )
c          else ! if ( cgscre .eq. 2 ) then
c              cgeps = 10.0d0 ** ( acgeps * log10( gpsupn ) + bcgeps )
c          end if
c          cgeps = max( cgepsf, min( cgepsi, cgeps ) )

C         Call conjugate gradients

          call cg(nind,ind,n,x,m,lambda,rho,g,delta,l,u,cgeps,epsnqmp,
     +    maxitnqmp,cgmaxit,nearlyq,gtype,htvtype,trtype,iprint,ncomp,d,
     +    cgiter,rbdtype,rbdind,inform,w(1),w(n+1),w(2*n+1),w(3*n+1),
     +    w(4*n+1),theta,sterel,steabs,epsrel,epsabs,infrel,infabs)

          cgcnt = cgcnt + cgiter

          if ( debug_trace ) then
              call packmol_gencan_debug_trace_tn_result('tn-cg',
     +        iter, inform, cgiter, rbdtype, rbdind, f)
          end if

          if ( inform .lt. 0 ) then

              if ( iprint .ge. 3 ) then
                  write(*, 1000) inform
                  write(10,1000) inform
              end if

              return

          end if

C         Compute maximum step

          if ( inform .eq. 2 ) then
              amax = 1.0d0
          else
              amax = infabs
              do i = 1,nind
                  if ( d(i) .gt. 0.0d0 ) then
                      amaxx = ( u(i) - x(i) ) / d(i)
                      if ( amaxx .lt. amax ) then
                          amax    = amaxx
                          rbdind  = i
                          rbdtype = 2
                      end if
                  else if ( d(i) .lt. 0.0d0 ) then
                      amaxx = ( l(i) - x(i) ) / d(i)
                      if ( amaxx .lt. amax ) then
                          amax    = amaxx
                          rbdind  = i
                          rbdtype = 1
                      end if
                  end if
               end do
          end if

C         Perform the line search

          tnintprev = tnintcnt
          tnexgprev = tnexgcnt
          tnexbprev = tnexbcnt

          fcntprev = fcnt

          call tnls(nind,ind,n,x,m,lambda,rho,l,u,f,g,d,amax,rbdtype,
     +    rbdind,nint,next,mininterp,maxextrap,fmin,maxfc,gtype,iprint,
     +    fcnt,gcnt,tnintcnt,tnexgcnt,tnexbcnt,inform,w(1),w(n+1),
     +    w(2*n+1),gamma,beta,sigma1,sigma2,sterel,steabs,epsrel,epsabs,
     +    infrel,infabs)

          if ( debug_trace ) then
              call packmol_gencan_debug_trace_tn_result('tn-line',
     +        iter, inform, cgiter, rbdtype, rbdind, f)
          end if

          if ( inform .lt. 0 ) then

              if ( iprint .ge. 3 ) then
                  write(*, 1000) inform
                  write(10,1000) inform
              end if

              return

          end if

          if ( tnintcnt .gt. tnintprev ) then
              tnintfe = tnintfe + ( fcnt - fcntprev )
          else if ( tnexgcnt .gt. tnexgprev ) then
              tnexgfe = tnexgfe + ( fcnt - fcntprev )
          else if ( tnexbcnt .gt. tnexbprev ) then
              tnexbfe = tnexbfe + ( fcnt - fcntprev )
          else
              tnstpcnt = tnstpcnt + 1
          end if

          tnfcnt = tnfcnt + ( fcnt - fcntprev )

C         Expand the point, its gradient and the bounds

          call expand(nind,ind,n,x)
          call expand(nind,ind,n,g)
          call expand(nind,ind,n,l)
          call expand(nind,ind,n,u)

C         If the line search (interpolation) in the Truncated Newton
C         direction stopped due to a very small step (inform = 6), we 
C         will discard this iteration and force a SPG iteration

C         Note that tnls subroutine was coded in such a way that in case
C         of inform = 6 termination the subroutine discards all what was 
C         done and returns with the same point it started

          if ( inform .eq. 6 ) then

              if ( iprint .ge. 3 ) then
                  write(*,*)  
                  write(*,*)  
     +            '     The previous TN iteration was discarded due to',
     +            '     a termination for very small step in the line ',
     +            '     search. A SPG iteration will be forced now.   '

                  write(10,*)  
                  write(10,*)  
     +            '     The previous TN iteration was discarded due to',
     +            '     a termination for very small step in the line ',
     +            '     search. A SPG iteration will be forced now.   '
              end if

              ittype  = 'SPG'
              spgiter = spgiter + 1

C             Compute spectral steplength

              if ( iter .eq. 1 .or. sty .le. 0.0d0 ) then
                  lamspg = max( 1.0d0, xnorm ) / sqrt( gpeucn2 )
              else
                  lamspg = sts / sty
              end if
              lamspg = min( lspgma, max( lspgmi, lamspg ) )

C             Perform a line search with safeguarded quadratic 
C             interpolation along the direction of the spectral 
C             continuous projected gradient

              fcntprev = fcnt

              call spgls(n,x,m,lambda,rho,f,g,l,u,lamspg,nint,mininterp,
     +        fmin,maxfc,iprint,fcnt,inform,w(1),w(n+1),gamma,sigma1,
     +        sigma2,sterel,steabs,epsrel,epsabs,infrel,infabs) 

              spgfcnt = spgfcnt + ( fcnt - fcntprev )

              if ( inform .lt. 0 ) then

                  if ( iprint .ge. 3 ) then
                      write(*, 1000) inform
                      write(10,1000) inform
                  end if

                  return
              end if

C             Compute the gradient at the new iterate

              infotmp = inform

              if ( gtype .eq. 0 ) then
                  call evalnal(n,x,m,lambda,rho,g,inform)
              else ! if ( gtype .eq. 1 ) then
                  call evalnaldiff(n,x,m,lambda,rho,g,sterel,steabs,
     +            inform)
              end if
              gcnt = gcnt + 1
 
              if ( inform .lt. 0 ) then

                  if ( iprint .ge. 3 ) then
                      write(*, 1000) inform
                      write(10,1000) inform
                  end if

                  return
              end if

              inform = infotmp

          end if

      end if

C     ==================================================================
C     Prepare for the next iteration 
C     ==================================================================

C     This adjustment/projection is ''por lo que las putas pudiera''

      do i = 1,n
          if ( x(i) .le. l(i) + max( epsrel * abs( l(i) ), epsabs ) ) 
     +    then
              x(i) = l(i)
          else if (x(i). ge. u(i) - max( epsrel * abs( u(i) ), epsabs )) 
     +    then  
              x(i) = u(i)
          end if
      end do

C     Compute x Euclidian norm

      xnorm = 0.0d0
      do i = 1,n
          xnorm = xnorm + x(i) ** 2
      end do
      xnorm = sqrt( xnorm )

C     Compute s = x_{k+1} - x_k, y = g_{k+1} - g_k, <s,s> and <s,y>

      sts = 0.0d0
      sty = 0.0d0
      do i = 1,n
          s(i) = x(i) - s(i)
          y(i) = g(i) - y(i)
          sts  = sts + s(i) ** 2
          sty  = sty + s(i) * y(i)
      end do

C     Compute continuous-project-gradient Euclidian and Sup norms,
C     internal gradient Euclidian norm, and store in nind the number of
C     free variables and in array ind their identifiers.

      nind    = 0
      gpsupn  = 0.0d0
      gpeucn2 = 0.0d0
      gieucn2 = 0.0d0
      do i = 1,n
          gpi     = min( u(i), max( l(i), x(i) - g(i) ) ) - x(i)
          gpsupn  = max( gpsupn, abs( gpi ) )
          gpeucn2 = gpeucn2 + gpi ** 2
          if ( x(i) .gt. l(i) .and. x(i) .lt. u(i) ) then
              gieucn2   = gieucn2 + gpi ** 2
              nind      = nind + 1
              ind(nind) = i
          end if
      end do

C     Print information of this iteration
c LM: progress bar for packmol
      call printbar(iprint, iter, maxit)

      if ( iprint .ge. 3 ) then 
          write(*, 983) iter,ittype
          write(*, 985) nprint,(x(i),i=1,nprint)
          write(*, 986) nprint,(g(i),i=1,nprint)
          write(*, 987) nprint,(min(u(i),max(l(i),x(i)-g(i)))-x(i),i=1,
     +    nprint)
          write(*, 988) min0(nprint,nind),nind,(ind(i),i=1,min0(nprint,
     +    nind))
          write(*, 1002) f,sqrt(gpeucn2),sqrt(gieucn2),gpsupn,nind,n,
     +    spgiter,tniter,fcnt,gcnt,cgcnt

          write(10,983) iter,ittype
          write(10,985) nprint,(x(i),i=1,nprint)
          write(10,986) nprint,(g(i),i=1,nprint)
          write(10,987) nprint,(min(u(i),max(l(i),x(i)-g(i)))-x(i),i=1,
     +    nprint)
          write(10,988) min0(nprint,nind),nind,(ind(i),i=1,min0(nprint,
     +    nind))
          write(10,1002) f,sqrt(gpeucn2),sqrt(gieucn2),gpsupn,nind,n,
     +    spgiter,tniter,fcnt,gcnt,cgcnt
      end if

C     ==================================================================
C     Test some stopping criteria that may occur inside the line 
C     searches 
C     ==================================================================

      if ( inform .eq. 6 ) then

          if ( iprint .ge. 3 ) then
              write(*, 996) inform,mininterp,epsrel,epsabs
              write(10,996) inform,mininterp,epsrel,epsabs
          end if

          go to 500

      end if

C     ==================================================================
C     Iterate 
C     ==================================================================

      go to 100

C     ==================================================================
C     End of main loop
C     ==================================================================

C     ==================================================================
C     Report output status and return
C     ==================================================================

 500  continue

C     Print final information

      if ( iprint .ge. 3 ) then
          write(*, 982) iter
          write(*, 985) nprint,(x(i),i=1,nprint)
          write(*, 986) nprint,(g(i),i=1,nprint)
          write(*, 987) nprint,(min(u(i),max(l(i),x(i)-g(i)))-x(i),i=1,
     +    nprint)
          write(*, 988) min0(nprint,nind),nind,(ind(i),i=1,min0(nprint,
     +    nind))
          write(*, 1002) f,sqrt(gpeucn2),sqrt(gieucn2),gpsupn,nind,n,
     +    spgiter,tniter,fcnt,gcnt,cgcnt

          write(10,982) iter
          write(10,985) nprint,(x(i),i=1,nprint)
          write(10,986) nprint,(g(i),i=1,nprint)
          write(10,987) nprint,(min(u(i),max(l(i),x(i)-g(i)))-x(i),i=1,
     +    nprint)
          write(10,988) min0(nprint,nind),nind,(ind(i),i=1,min0(nprint,
     +    nind))
          write(10,1002) f,sqrt(gpeucn2),sqrt(gieucn2),gpsupn,nind,n,
     +    spgiter,tniter,fcnt,gcnt,cgcnt
      end if

      return 

C     Non-executable statements

 977  format(/1X, 'Entry to GENCAN. Number of variables: ',I7)
 978  format(/1X,'Lower bounds (first ',I6, ' components): ',
     */,6(1X,1PD11.4))
 979  format(/1X,'Upper bounds (first ',I6, ' components): ',
     */,6(1X,1PD11.4))
 980  format(/1X,'Initial point (first ',I6, ' components): ',
     */,6(1X,1PD11.4))
 981  format(/1X,'GENCAN iteration: ',I6, ' (Initial point)')
 982  format(/1X,'GENCAN iteration: ',I6, ' (Final point)')
 983  format(/,1X,'GENCAN iteration: ',I6,
     *' (This point was obtained using a ',A3,' iteration)')
 985  format(1X,'Current point (first ',I6, ' components): ',
     */,6(1X,1PD11.4))
 986  format(1X,'Current gradient (first ',I6, ' components): ',
     */,6(1X,1PD11.4))
 987  format(1X,'Current continuous projected gradient (first ',I6, 
     *' components): ',/,6(1X,1PD11.4))
 988  format(1X,'Current free variables (first ',I6,
     *', total number ',I6,'): ',/,10(1X,I6))
 990  format(/1X,'Flag of GENCAN = ',I3,
     *' (convergence with Euclidian-norm of the projected gradient',
     */,1X,'smaller than ',1PD11.4,')')
 991  format(/1X,'Flag of GENCAN = ',I3,
     *' (convergence with sup-norm of the projected gradient',
     */,1X,'smaller than ',1PD11.4,')')
 992  format(/1X,'Flag of GENCAN= ',I3,
     *' (The algorithm stopped by lack of enough progress. This means',
     */,1X,'that  f(x_k) - f(x_{k+1}) .le. ',1PD11.4,
     *' * max [ f(x_j)-f(x_{j+1}, j < k ]',/,1X,'during ',I7,
     *' consecutive iterations')
 993  format(/1X,'Flag of GENCAN = ',I3,
     *' (The algorithm stopped because the order of the',
     */,1X,'Euclidian-norm of the continuous projected gradient did',
     *' not change during ',/,1X,I7,' consecutive iterations.',
     *' Probably, an exaggerated small norm of the',/,1X,'continuous',
     *' projected gradient is required for declaring convergence')
 994  format(/1X,'Flag of GENCAN = ',I3,
     *' (The algorithm stopped because the functional value is',
     */,1X,'smaller than ',1PD11.4)
 996  format(/1X,'Flag of GENCAN = ',I3,
     *' (Too small step in a line search. After having made at ',
     */,1X,'least ',I7,' interpolations, the steplength becames small.',
     *' Small means that',/,1X,'we were at point x with direction d',
     *' and took a step  alpha such that',/,1X,'alpha * |d_i| .lt.',
     *' max [',1PD11.4,' * |x_i|,',1PD11.4,' ] for all i)')
 997  format(/1X,'Flag of GENCAN = ',I3,
     *' (It was exceeded the maximum allowed number of iterations',
     */,1X,'(maxit=',I7,')')
 998  format(/1X,'Flag of GENCAN = ',I3,
     *' (It was exceeded the maximum allowed number of functional',
     */,1X,'evaluations (maxfc=',I7,')')
 1002 format(1X,'Functional value: ', 1PD11.4,
     */,1X,'Euclidian-norm of the continuous projected gradient: ',
     *1PD11.4,
     */,1X,'Euclidian-norm of the internal projection of gp: ',1PD11.4,
     */,1X,'Sup-norm of the continuous projected gradient: ',1PD11.4,
     */,1X,'Free variables at this point: ',I7,
     *' (over a total of ',I7,')',
     */,1X,'SPG iterations: ',I7,
     */,1X,'TN iterations: ',I7,
     */,1X,'Functional evaluations: ',I7,
     */,1X,'Gradient evaluations: ',I7,
     */,1X,'Conjugate gradient iterations: ',I7)
C1003 format(6X,I6,T22,D17.6,T43,D17.6)
C1003 format(6X,'Iter = ',I6,' f = ',1PD11.4,' gpsupn = ',1PD11.4)
 1000 format(/1X,'Flag of GENCAN = ',I3,' Fatal Error')

      end

C     ******************************************************************
C     ******************************************************************
