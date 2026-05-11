program test_cg_direct_compare_probe
   use iso_c_binding, only : c_bool, c_double, c_int
   use compute_data, only : init1
   use test_runtime_setup
   use usegencan, only : g, l, u
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
   integer :: n, nind, itype, i
   integer(c_int) :: iter_cpp, rbdtype_cpp, rbdind_cpp, inform_cpp, used_cpp
   integer :: iter_fortran, rbdtype_fortran, rbdind_fortran, inform_fortran
   integer, allocatable :: ind(:)
   double precision :: lambda(m), rho(m), xnorm, delta, cgeps, diff_s_max, diff_s_rms
   double precision :: gnorm2, gpsupn
   double precision, allocatable :: x_fortran(:), x_cpp(:), g_fortran(:), g_cpp(:)
   double precision, allocatable :: s_fortran(:), s_cpp(:), w_fortran(:), w_cpp(:), y_fortran(:), y_cpp(:)
   double precision, allocatable :: r_fortran(:), r_cpp(:), d_fortran(:), d_cpp(:), sprev_fortran(:), sprev_cpp(:)
   character(len=1000) :: input_path
   logical :: nearlyq
   external :: initial, swaptype, computeg, cg

   call get_command_argument(1, input_path)
   call load_runtime_problem(trim(input_path), n, x_fortran)
   call initial(n, x_fortran)

   itype = 1
   init1 = .true.
   call swaptype(n, x_fortran, itype, 0)
   call swaptype(n, x_fortran, itype, 1)

   nind = n
   allocate(ind(nind), x_cpp(n), g_fortran(n), g_cpp(n))
   allocate(s_fortran(n), s_cpp(n), w_fortran(n), w_cpp(n), y_fortran(n), y_cpp(n))
   allocate(r_fortran(n), r_cpp(n), d_fortran(n), d_cpp(n), sprev_fortran(n), sprev_cpp(n))

   do i = 1, nind
      ind(i) = i
   end do

   x_cpp = x_fortran
   lambda = 0.d0
   rho = 0.d0
   call computeg(n, x_fortran, g_fortran)
   g_cpp = g_fortran

   xnorm = sqrt(sum(x_fortran(1:n) * x_fortran(1:n)))
   delta = max(1.d-2, 0.1d0 * max(1.d0, xnorm))
   cgeps = 1.d-1
   gnorm2 = sum(g_fortran(1:n) * g_fortran(1:n))
   gpsupn = maxval(abs(g_fortran(1:n)))

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

   nearlyq = .false.
   call cg(nind, ind, n, x_fortran, m, lambda, rho, g_fortran, delta, l, u, cgeps, 1.d-4, 5, 20, nearlyq, 0, 1, 1, 0, n, &
      s_fortran, iter_fortran, rbdtype_fortran, rbdind_fortran, inform_fortran, w_fortran, y_fortran, r_fortran, d_fortran, &
      sprev_fortran, 1.d-6, 1.d-7, 1.d-10, 1.d-10, 1.d-20, 1.d20, 1.d99)

   call packmol_cg_cpp_direct_c(nind, ind, n, x_cpp, m, lambda, rho, g_cpp, delta, l, u, cgeps, 1.d-4, 5, 20, .false._c_bool, &
      0, 1, 1, 0, n, s_cpp, iter_cpp, rbdtype_cpp, rbdind_cpp, inform_cpp, w_cpp, y_cpp, r_cpp, d_cpp, sprev_cpp, 1.d-6, 1.d-7, &
      1.d-10, 1.d-10, 1.d-20, 1.d20, 1.d99, used_cpp)

   diff_s_max = maxval(abs(s_fortran(1:nind) - s_cpp(1:nind)))
   diff_s_rms = sqrt(sum((s_fortran(1:nind) - s_cpp(1:nind))**2) / dble(max(1, nind)))

   write(*,'(A,1X,I0,1X,A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16)') &
      'n', n, 'xnorm', xnorm, 'delta', delta, 'gpsupn', gpsupn
   write(*,'(A,1X,ES24.16,1X,A,1X,I0)') 'gnorm2', gnorm2, 'used_cpp', used_cpp
   write(*,'(A,1X,I0,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,I0)') &
      'fortran_inform', inform_fortran, 'fortran_iter', iter_fortran, 'fortran_rbdtype', rbdtype_fortran, 'fortran_rbdind', rbdind_fortran
   write(*,'(A,1X,I0,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,I0)') &
      'cpp_inform', inform_cpp, 'cpp_iter', iter_cpp, 'cpp_rbdtype', rbdtype_cpp, 'cpp_rbdind', rbdind_cpp
   write(*,'(A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16)') &
      'fortran_snorm2', sum(s_fortran(1:nind) * s_fortran(1:nind)), 'cpp_snorm2', sum(s_cpp(1:nind) * s_cpp(1:nind)), &
      'diff_s_max', diff_s_max, 'diff_s_rms', diff_s_rms

   call swaptype(n, x_fortran, itype, 3)
   init1 = .false.
end program test_cg_direct_compare_probe
