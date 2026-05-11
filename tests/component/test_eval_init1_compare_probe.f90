program test_eval_init1_compare_probe
   use iso_c_binding, only : c_double, c_int
   use compute_data, only : init1, fdist, frest
   use test_runtime_setup
   implicit none

   interface
      subroutine packmol_evalal_cpp_direct_c(n, x, m, lambda, rho, f, flag, used_cpp) bind(C, name="packmol_evalal_cpp_direct_c")
         import :: c_double, c_int
         integer(c_int), intent(in) :: n, m
         real(c_double), intent(inout) :: x(*)
         real(c_double), intent(in) :: lambda(*), rho(*)
         real(c_double), intent(out) :: f
         integer(c_int), intent(out) :: flag, used_cpp
      end subroutine packmol_evalal_cpp_direct_c

      subroutine packmol_evalnal_cpp_direct_c(n, x, m, lambda, rho, g, flag, used_cpp) bind(C, name="packmol_evalnal_cpp_direct_c")
         import :: c_double, c_int
         integer(c_int), intent(in) :: n, m
         real(c_double), intent(inout) :: x(*)
         real(c_double), intent(in) :: lambda(*), rho(*)
         real(c_double), intent(out) :: g(*)
         integer(c_int), intent(out) :: flag, used_cpp
      end subroutine packmol_evalnal_cpp_direct_c
   end interface

   integer :: n, itype, i, nind
   integer(c_int) :: inform_cpp, used_cpp
   double precision :: f_fortran, f_cpp, diff_g_max, diff_g_rms
   double precision :: fortran_fdist, fortran_frest, cpp_fdist, cpp_frest
   double precision :: gpsupn, gpeucn2, gieucn2, gpi, xnorm, ometa2
   double precision :: lambda(1), rho(1)
   double precision, allocatable :: x(:), x_cpp(:), g_fortran(:), g_cpp(:)
   character(len=1000) :: input_path
   external :: initial, computef, computeg, swaptype

   call get_command_argument(1, input_path)
   call load_runtime_problem(trim(input_path), n, x)
   call initial(n, x)

   allocate(x_cpp(size(x)), g_fortran(size(x)), g_cpp(size(x)))
   x_cpp = x
   lambda = 0.d0
   rho = 0.d0

   itype = 1
   init1 = .true.
   call swaptype(n, x, itype, 0)
   call swaptype(n, x, itype, 1)
   call swaptype(n, x_cpp, itype, 0)
   call swaptype(n, x_cpp, itype, 1)

   call computef(n, x, f_fortran)
   fortran_fdist = fdist
   fortran_frest = frest
   call packmol_evalal_cpp_direct_c(n, x_cpp, 1_c_int, lambda, rho, f_cpp, inform_cpp, used_cpp)
   cpp_fdist = fdist
   cpp_frest = frest

   call computeg(n, x, g_fortran(1:n))
   call packmol_evalnal_cpp_direct_c(n, x_cpp, 1_c_int, lambda, rho, g_cpp(1:n), inform_cpp, used_cpp)
   diff_g_max = maxval(abs(g_fortran(1:n) - g_cpp(1:n)))
   diff_g_rms = sqrt(sum((g_fortran(1:n) - g_cpp(1:n))**2) / dble(max(1, n)))
   gpsupn = 0.d0
   gpeucn2 = 0.d0
   gieucn2 = 0.d0
   nind = 0
   do i = 1, n
      gpi = -g_fortran(i)
      gpsupn = max(gpsupn, abs(gpi))
      gpeucn2 = gpeucn2 + gpi * gpi
      gieucn2 = gieucn2 + gpi * gpi
      nind = nind + 1
   end do
   xnorm = sqrt(sum(x(1:n) * x(1:n)))
   ometa2 = (1.d0 - 0.9d0) ** 2

   write(*,'(A,1X,I0,1X,A,1X,I0,1X,A,1X,ES24.16,1X,A,1X,ES24.16)') &
      'n_init1', n, 'used_cpp', used_cpp, 'f_fortran', f_fortran, 'f_cpp', f_cpp
   write(*,'(A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16)') &
      'fdist_fortran', fortran_fdist, 'fdist_cpp', cpp_fdist, 'frest_fortran', fortran_frest, 'frest_cpp', cpp_frest
   write(*,'(A,1X,ES24.16,1X,A,1X,ES24.16)') &
      'diff_g_max', diff_g_max, 'diff_g_rms', diff_g_rms
   write(*,'(A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,I0)') &
      'gpsupn', gpsupn, 'gpeucn2', gpeucn2, 'gieucn2', gieucn2, 'xnorm', xnorm, 'nind', nind
   write(*,'(A,1X,I0)') 'spg_test', merge(1, 0, gieucn2 <= ometa2 * gpeucn2)

   call swaptype(n, x, itype, 3)
   call swaptype(n, x_cpp, itype, 3)
   init1 = .false.
end program test_eval_init1_compare_probe
