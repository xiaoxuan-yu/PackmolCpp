program test_eval_reduced_compare_probe
   use iso_c_binding, only : c_double, c_int
   use test_runtime_setup
   implicit none

   interface
      subroutine packmol_calcf_cpp_reduced_direct_c(nind, ind, x, n, xc, m, lambda, rho, f, inform) bind(C, name="packmol_calcf_cpp_reduced_direct_c")
         import :: c_double, c_int
         integer(c_int), intent(in) :: nind, n, m
         integer(c_int), intent(in) :: ind(*)
         real(c_double), intent(inout) :: x(*)
         real(c_double), intent(in) :: xc(*), lambda(*), rho(*)
         real(c_double), intent(out) :: f
         integer(c_int), intent(out) :: inform
      end subroutine packmol_calcf_cpp_reduced_direct_c

      subroutine packmol_calcg_cpp_reduced_direct_c(nind, ind, x, n, xc, m, lambda, rho, gtype, g, sterel, steabs, inform) bind(C, name="packmol_calcg_cpp_reduced_direct_c")
         import :: c_double, c_int
         integer(c_int), intent(in) :: nind, n, m, gtype
         integer(c_int), intent(in) :: ind(*)
         real(c_double), intent(inout) :: x(*)
         real(c_double), intent(in) :: xc(*), lambda(*), rho(*), sterel, steabs
         real(c_double), intent(out) :: g(*)
         integer(c_int), intent(out) :: inform
      end subroutine packmol_calcg_cpp_reduced_direct_c
   end interface

   integer :: n, nind, i
   integer(c_int) :: inform_cpp
   integer, allocatable :: ind(:)
   double precision :: lambda(1), rho(1), f_fortran, f_cpp
   double precision :: diff_x_max, diff_g_max, diff_g_rms
   double precision, allocatable :: xc(:), xred_fortran(:), xred_cpp(:), g_fortran(:), g_cpp(:)
   character(len=1000) :: input_path
   external :: initial, calcf, calcg

   call get_command_argument(1, input_path)
   call load_runtime_problem(trim(input_path), n, xc)
   call initial(n, xc)

   nind = min(max(2, n / 2), n)
   allocate(ind(nind), xred_fortran(n), xred_cpp(n), g_fortran(n), g_cpp(n))
   do i = 1, nind
      ind(i) = 2 * i - 1
      if (ind(i) > n) ind(i) = i
   end do

   xred_fortran = -777.d0
   xred_cpp = -777.d0
   do i = 1, nind
      xred_fortran(i) = xc(ind(i))
      xred_cpp(i) = xc(ind(i))
   end do
   if (nind >= 1) xred_fortran(1) = xred_fortran(1) + 0.137d0
   if (nind >= 2) xred_fortran(2) = xred_fortran(2) - 0.219d0
   if (nind >= 3) xred_fortran(3) = xred_fortran(3) + 0.083d0
   if (nind >= 4) xred_fortran(4) = xred_fortran(4) - 0.051d0
   xred_cpp(1:nind) = xred_fortran(1:nind)

   lambda = 0.d0
   rho = 0.d0

   call calcf(nind, ind, xred_fortran, n, xc, 1, lambda, rho, f_fortran, i)
   call packmol_calcf_cpp_reduced_direct_c(nind, ind, xred_cpp, n, xc, 1_c_int, lambda, rho, f_cpp, inform_cpp)
   diff_x_max = maxval(abs(xred_fortran(1:nind) - xred_cpp(1:nind)))

   call calcg(nind, ind, xred_fortran, n, xc, 1, lambda, rho, g_fortran, i)
   call packmol_calcg_cpp_reduced_direct_c(nind, ind, xred_cpp, n, xc, 1_c_int, lambda, rho, 0_c_int, g_cpp, 1.d-6, 1.d-8, inform_cpp)
   diff_g_max = maxval(abs(g_fortran(1:nind) - g_cpp(1:nind)))
   diff_g_rms = sqrt(sum((g_fortran(1:nind) - g_cpp(1:nind))**2) / dble(nind))

   write(*,'(A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16)') &
      'f_fortran', f_fortran, 'f_cpp', f_cpp, 'diff_x_max', diff_x_max
   write(*,'(A,1X,ES24.16,1X,A,1X,ES24.16)') 'diff_g_max', diff_g_max, 'diff_g_rms', diff_g_rms
end program test_eval_reduced_compare_probe
