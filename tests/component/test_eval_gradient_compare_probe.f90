program test_eval_gradient_compare_probe
   use iso_c_binding, only : c_double, c_int
   use test_runtime_setup
   implicit none

   interface
      subroutine packmol_evalnal_cpp_direct_c(n, x, m, lambda, rho, g, flag, used_cpp) bind(C, name="packmol_evalnal_cpp_direct_c")
         import :: c_double, c_int
         integer(c_int), intent(in) :: n, m
         real(c_double), intent(inout) :: x(*)
         real(c_double), intent(in) :: lambda(*), rho(*)
         real(c_double), intent(out) :: g(*)
         integer(c_int), intent(out) :: flag, used_cpp
      end subroutine packmol_evalnal_cpp_direct_c
   end interface

   integer :: n
   integer(c_int) :: inform_cpp, used_cpp
   double precision :: lambda(1), rho(1)
   double precision, allocatable :: x(:), x_cpp(:), x_perturbed(:), x_cpp_perturbed(:)
   double precision, allocatable :: g_fortran(:), g_cpp(:), g_fortran_perturbed(:), g_cpp_perturbed(:)
   double precision :: diff_max, diff_rms, diff_max_perturbed, diff_rms_perturbed
   character(len=1000) :: input_path
   external :: initial, computeg

   call get_command_argument(1, input_path)
   call load_runtime_problem(trim(input_path), n, x)
   call initial(n, x)

   allocate(x_cpp(n), x_perturbed(n), x_cpp_perturbed(n))
   allocate(g_fortran(n), g_cpp(n), g_fortran_perturbed(n), g_cpp_perturbed(n))
   x_cpp = x
   x_perturbed = x
   x_cpp_perturbed = x
   if (n >= 6) then
      x_perturbed(1:6) = x_perturbed(1:6) + (/ 0.137d0, -0.219d0, 0.083d0, 0.11d0, -0.07d0, 0.05d0 /)
      x_cpp_perturbed(1:6) = x_perturbed(1:6)
   end if
   if (n >= 12) then
      x_perturbed(7:12) = x_perturbed(7:12) + (/ -0.091d0, 0.064d0, -0.147d0, -0.03d0, 0.09d0, -0.12d0 /)
      x_cpp_perturbed(7:12) = x_perturbed(7:12)
   end if
   lambda = 0.d0
   rho = 0.d0

   call computeg(n, x, g_fortran)
   call packmol_evalnal_cpp_direct_c(n, x_cpp, 1_c_int, lambda, rho, g_cpp, inform_cpp, used_cpp)
   diff_max = maxval(abs(g_fortran - g_cpp))
   diff_rms = sqrt(sum((g_fortran - g_cpp)**2) / dble(n))

   call computeg(n, x_perturbed, g_fortran_perturbed)
   call packmol_evalnal_cpp_direct_c(n, x_cpp_perturbed, 1_c_int, lambda, rho, g_cpp_perturbed, inform_cpp, used_cpp)
   diff_max_perturbed = maxval(abs(g_fortran_perturbed - g_cpp_perturbed))
   diff_rms_perturbed = sqrt(sum((g_fortran_perturbed - g_cpp_perturbed)**2) / dble(n))

   write(*,'(A,1X,I0,1X,A,1X,I0,1X,A,1X,ES24.16,1X,A,1X,ES24.16)') &
      'used_cpp', used_cpp, 'inform_cpp', inform_cpp, 'diff_max', diff_max, 'diff_rms', diff_rms
   write(*,'(A,1X,ES24.16,1X,A,1X,ES24.16)') &
      'perturbed_diff_max', diff_max_perturbed, 'perturbed_diff_rms', diff_rms_perturbed
end program test_eval_gradient_compare_probe
