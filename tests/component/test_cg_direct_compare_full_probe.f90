program test_cg_direct_compare_full_probe
   use iso_c_binding, only : c_bool, c_double, c_int
   use compute_data, only : init1
   use test_runtime_setup
   use usegencan, only : g, l, u, maxit
   implicit none

   interface
      subroutine packmol_cg_cpp_direct_c(nind, ind, n, x, m, lambda, rho, g, delta, l, u, eps, epsnqmp, &
         maxitnqmp, maxit, nearlyq, gtype, htvtype, trtype, iprint, ncomp, s, iter, rbdtype, rbdind, inform, &
         w, y, r, d, sprev, theta, sterel, steabs, epsrel, epsabs, infrel, infabs, used_cpp) bind(C, name="packmol_cg_cpp_direct_c")
         import :: c_bool, c_double, c_int
         integer(c_int), intent(in) :: nind, n, m, maxitnqmp, maxit, gtype, htvtype, trtype, iprint, ncomp
         logical(c_bool), intent(in) :: nearlyq
         integer(c_int), intent(in) :: ind(*)
         real(c_double), intent(inout) :: x(*), g(*)
         real(c_double), intent(in) :: lambda(*), rho(*), delta, l(*), u(*), eps, epsnqmp
         real(c_double), intent(inout) :: s(*), w(*), y(*), r(*), d(*), sprev(*)
         integer(c_int), intent(out) :: iter, rbdtype, rbdind, inform, used_cpp
         real(c_double), intent(in) :: theta, sterel, steabs, epsrel, epsabs, infrel, infabs
      end subroutine packmol_cg_cpp_direct_c
   end interface

   integer, parameter :: m = 1
   integer :: n, nind, i, cgmaxit, maxitnqmp, cgscre
   integer(c_int) :: iter_cpp, rbdtype_cpp, rbdind_cpp, inform_cpp, used_cpp
   integer :: iter_fortran, rbdtype_fortran, rbdind_fortran, inform_fortran
   integer, allocatable :: ind(:)
   double precision :: lambda(m), rho(m), delta, cgeps, diff_s_max, diff_s_rms
   double precision :: gpsupn, gpeucn2, gieucn2, xnorm, gpi, ometa2, fx
   double precision :: acgeps, bcgeps, cgepsi, cgepsf, cggpnf, epsgpsn, epsgpen2, epsnqmp, kappa
   double precision, parameter :: theta = 1.d-6, sterel = 1.d-7, steabs = 1.d-10
   double precision, parameter :: epsrel = 1.d-10, epsabs = 1.d-20, infrel = 1.d20, infabs = 1.d99
   double precision, allocatable :: x_full(:), x_fortran(:), x_cpp(:), g_fortran(:), g_cpp(:)
   double precision, allocatable :: l_fortran(:), l_cpp(:), u_fortran(:), u_cpp(:)
   double precision, allocatable :: s_fortran(:), s_cpp(:), w_fortran(:), w_cpp(:), y_fortran(:), y_cpp(:)
   double precision, allocatable :: r_fortran(:), r_cpp(:), d_fortran(:), d_cpp(:), sprev_fortran(:), sprev_cpp(:)
   character(len=1000) :: input_path
   logical :: nearlyq
   external :: initial, computef, computeg, shrink, cg, gp_ieee_signal1, gp_ieee_signal2

   call get_command_argument(1, input_path)
   call load_runtime_problem(trim(input_path), n, x_full)
   call initial(n, x_full)

   init1 = .false.
   call computef(n, x_full, fx)
   allocate(ind(n), x_fortran(n), x_cpp(n), g_fortran(n), g_cpp(n))
   allocate(l_fortran(n), l_cpp(n), u_fortran(n), u_cpp(n))
   allocate(s_fortran(n), s_cpp(n), w_fortran(n), w_cpp(n), y_fortran(n), y_cpp(n))
   allocate(r_fortran(n), r_cpp(n), d_fortran(n), d_cpp(n), sprev_fortran(n), sprev_cpp(n))

   call computeg(n, x_full, g_fortran)
   g_cpp = g_fortran
   gpsupn = 0.d0
   gpeucn2 = 0.d0
   gieucn2 = 0.d0
   nind = 0
   do i = 1, n
      gpi = min(u(i), max(l(i), x_full(i) - g_fortran(i))) - x_full(i)
      gpsupn = max(gpsupn, abs(gpi))
      gpeucn2 = gpeucn2 + gpi * gpi
      if (x_full(i) .gt. l(i) .and. x_full(i) .lt. u(i)) then
         gieucn2 = gieucn2 + gpi * gpi
         nind = nind + 1
         ind(nind) = i
      end if
   end do

   xnorm = sqrt(sum(x_full(1:n) * x_full(1:n)))
   ometa2 = (1.d0 - 0.9d0) ** 2
   if (gieucn2 .gt. ometa2 * gpeucn2) then
      write(*,'(A,1X,ES24.16,1X,A,1X,ES24.16,1X,A)') 'gieucn2', gieucn2, 'gpeucn2', gpeucn2, 'tn_expected false'
      stop 1
   end if

   delta = max(1.d-2, 0.1d0 * max(1.d0, xnorm))
   epsgpsn = 1.d-6
   cggpnf = max(1.d-4, epsgpsn)
   cgscre = 2
   cgepsi = 1.d-1
   cgepsf = 1.d-5
   epsnqmp = 1.d-4
   maxitnqmp = 5
   nearlyq = .false.
   epsgpen2 = 0.d0
   call gp_ieee_signal1(gpsupn, acgeps, bcgeps, cgepsf, cgepsi, cggpnf)
   call gp_ieee_signal2(cgmaxit, nind, nearlyq, -1, cgscre, kappa, gpeucn2, gpeucn2, epsgpen2, epsgpsn, &
      cgeps, acgeps, bcgeps, cgepsf, cgepsi, gpsupn, gpsupn)

   x_fortran = x_full
   x_cpp = x_full
   l_fortran = l
   l_cpp = l
   u_fortran = u
   u_cpp = u
   call shrink(nind, ind, n, x_fortran)
   call shrink(nind, ind, n, x_cpp)
   call shrink(nind, ind, n, g_fortran)
   call shrink(nind, ind, n, g_cpp)
   call shrink(nind, ind, n, l_fortran)
   call shrink(nind, ind, n, l_cpp)
   call shrink(nind, ind, n, u_fortran)
   call shrink(nind, ind, n, u_cpp)

   lambda = 0.d0
   rho = 0.d0
   s_fortran = 0.d0
   s_cpp = 0.d0
   w_fortran = 0.d0
   w_cpp = 0.d0
   y_fortran = 0.d0
   y_cpp = 0.d0
   r_fortran = 0.d0
   r_cpp = 0.d0
   d_fortran = 0.d0
   d_cpp = 0.d0
   sprev_fortran = 0.d0
   sprev_cpp = 0.d0
   iter_fortran = 0

   call cg(nind, ind, n, x_fortran, m, lambda, rho, g_fortran, delta, l_fortran, u_fortran, cgeps, epsnqmp, &
      maxitnqmp, cgmaxit, nearlyq, 0, 1, 1, 0, n, s_fortran, iter_fortran, rbdtype_fortran, rbdind_fortran, &
      inform_fortran, w_fortran, y_fortran, r_fortran, d_fortran, sprev_fortran, theta, sterel, steabs, &
      epsrel, epsabs, infrel, infabs)

   call packmol_cg_cpp_direct_c(nind, ind, n, x_cpp, m, lambda, rho, g_cpp, delta, l_cpp, u_cpp, cgeps, epsnqmp, &
      maxitnqmp, cgmaxit, .false._c_bool, 0, 1, 1, 0, n, s_cpp, iter_cpp, rbdtype_cpp, rbdind_cpp, inform_cpp, &
      w_cpp, y_cpp, r_cpp, d_cpp, sprev_cpp, theta, sterel, steabs, epsrel, epsabs, infrel, infabs, used_cpp)

   diff_s_max = maxval(abs(s_fortran(1:nind) - s_cpp(1:nind)))
   diff_s_rms = sqrt(sum((s_fortran(1:nind) - s_cpp(1:nind))**2) / dble(max(1, nind)))

   write(*,'(A,1X,I0,1X,A,1X,I0,1X,A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16)') &
      'n', n, 'nind', nind, 'xnorm', xnorm, 'delta', delta, 'gpsupn', gpsupn
   write(*,'(A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,I0)') &
      'gpeucn2', gpeucn2, 'gieucn2', gieucn2, 'cgeps', cgeps, 'cgmaxit', cgmaxit
   write(*,'(A,1X,I0,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,I0)') &
      'fortran_inform', inform_fortran, 'fortran_iter', iter_fortran, 'fortran_rbdtype', rbdtype_fortran, 'fortran_rbdind', rbdind_fortran
   write(*,'(A,1X,I0,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,I0)') &
      'cpp_inform', inform_cpp, 'cpp_iter', iter_cpp, 'cpp_rbdtype', rbdtype_cpp, 'cpp_rbdind', rbdind_cpp, 'used_cpp', used_cpp
   write(*,'(A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16)') &
      'fortran_snorm2', sum(s_fortran(1:nind) * s_fortran(1:nind)), 'cpp_snorm2', sum(s_cpp(1:nind) * s_cpp(1:nind)), &
      'diff_s_max', diff_s_max, 'diff_s_rms', diff_s_rms
end program test_cg_direct_compare_full_probe
