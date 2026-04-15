program test_eval_objective_compare_probe
   use iso_c_binding, only : c_double, c_int
   use compute_data
   use input
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
   end interface

   integer :: n, inform_fortran
   integer(c_int) :: inform_cpp, used_cpp
   double precision :: f_fortran, f_cpp, f_fortran_perturbed, f_cpp_perturbed
   double precision :: fortran_fdist, fortran_frest, fortran_pair_sum, fortran_constraint_sum
   double precision :: cpp_fdist, cpp_frest, cpp_pair_sum, cpp_constraint_sum
   double precision :: fortran_fdist_perturbed, fortran_frest_perturbed
   double precision :: cpp_fdist_perturbed, cpp_frest_perturbed
   double precision :: fortran_pair_sum_perturbed, fortran_constraint_sum_perturbed
   double precision :: cpp_pair_sum_perturbed, cpp_constraint_sum_perturbed
   integer :: fortran_pair_count, fortran_constraint_count, cpp_pair_count, cpp_constraint_count
   integer :: fortran_pair_count_perturbed, fortran_constraint_count_perturbed
   integer :: cpp_pair_count_perturbed, cpp_constraint_count_perturbed
   character(len=1000) :: input_path
   double precision, allocatable :: x(:), x_cpp(:), x_perturbed(:), x_cpp_perturbed(:)
   real(c_double) :: lambda(1), rho(1)
   external :: initial, computef

   call get_command_argument(1, input_path)
   call load_runtime_problem(trim(input_path), n, x)
   call initial(n, x)

   allocate(x_cpp(n), x_perturbed(n), x_cpp_perturbed(n))
   x_cpp = x
   x_perturbed = x
   x_cpp_perturbed = x
   lambda = 0.d0
   rho = 0.d0

   if (n >= 6) then
      x_perturbed(1:6) = x_perturbed(1:6) + (/ 0.137d0, -0.219d0, 0.083d0, 0.11d0, -0.07d0, 0.05d0 /)
      x_cpp_perturbed(1:6) = x_perturbed(1:6)
   end if
   if (n >= 12) then
      x_perturbed(7:12) = x_perturbed(7:12) + (/ -0.091d0, 0.064d0, -0.147d0, -0.03d0, 0.09d0, -0.12d0 /)
      x_cpp_perturbed(7:12) = x_perturbed(7:12)
   end if

   call computef(n, x, f_fortran)
   inform_fortran = 0
   fortran_fdist = fdist
   fortran_frest = frest
   fortran_pair_sum = pair_penalty_sum
   fortran_constraint_sum = constraint_penalty_sum
   fortran_pair_count = pair_penalty_count
   fortran_constraint_count = constraint_penalty_count

   call packmol_evalal_cpp_direct_c(n, x_cpp, 1_c_int, lambda, rho, f_cpp, inform_cpp, used_cpp)
   cpp_fdist = fdist
   cpp_frest = frest
   cpp_pair_sum = pair_penalty_sum
   cpp_constraint_sum = constraint_penalty_sum
   cpp_pair_count = pair_penalty_count
   cpp_constraint_count = constraint_penalty_count

    call computef(n, x_perturbed, f_fortran_perturbed)
   fortran_fdist_perturbed = fdist
   fortran_frest_perturbed = frest
   fortran_pair_sum_perturbed = pair_penalty_sum
   fortran_constraint_sum_perturbed = constraint_penalty_sum
   fortran_pair_count_perturbed = pair_penalty_count
   fortran_constraint_count_perturbed = constraint_penalty_count

   call packmol_evalal_cpp_direct_c(n, x_cpp_perturbed, 1_c_int, lambda, rho, f_cpp_perturbed, inform_cpp, used_cpp)
   cpp_fdist_perturbed = fdist
   cpp_frest_perturbed = frest
   cpp_pair_sum_perturbed = pair_penalty_sum
   cpp_constraint_sum_perturbed = constraint_penalty_sum
   cpp_pair_count_perturbed = pair_penalty_count
   cpp_constraint_count_perturbed = constraint_penalty_count

   write(*,'(A,1X,I0,1X,A,1X,I0,1X,A,1X,ES24.16,1X,A,1X,ES24.16)') &
      'used_cpp', used_cpp, 'inform_cpp', inform_cpp, 'f_fortran', f_fortran, 'f_cpp', f_cpp
   write(*,'(A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16)') &
      'fdist_fortran', fortran_fdist, 'fdist_cpp', cpp_fdist, 'frest_fortran', fortran_frest, 'frest_cpp', cpp_frest
   write(*,'(A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16)') &
      'pair_sum_fortran', fortran_pair_sum, 'pair_sum_cpp', cpp_pair_sum, 'constraint_sum_fortran', fortran_constraint_sum, 'constraint_sum_cpp', cpp_constraint_sum
   write(*,'(A,1X,I0,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,I0)') &
      'pair_count_fortran', fortran_pair_count, 'pair_count_cpp', cpp_pair_count, &
      'constraint_count_fortran', fortran_constraint_count, 'constraint_count_cpp', cpp_constraint_count
   write(*,'(A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16)') &
      'perturbed_f_fortran', f_fortran_perturbed, 'perturbed_f_cpp', f_cpp_perturbed, &
      'perturbed_fdist_fortran', fortran_fdist_perturbed, 'perturbed_fdist_cpp', cpp_fdist_perturbed
   write(*,'(A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16)') &
      'perturbed_frest_fortran', fortran_frest_perturbed, 'perturbed_frest_cpp', cpp_frest_perturbed, &
      'perturbed_pair_sum_fortran', fortran_pair_sum_perturbed, 'perturbed_pair_sum_cpp', cpp_pair_sum_perturbed
   write(*,'(A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,I0,1X,A,1X,I0)') &
      'perturbed_constraint_sum_fortran', fortran_constraint_sum_perturbed, &
      'perturbed_constraint_sum_cpp', cpp_constraint_sum_perturbed, &
      'perturbed_pair_count_fortran', fortran_pair_count_perturbed, 'perturbed_pair_count_cpp', cpp_pair_count_perturbed
   write(*,'(A,1X,I0,1X,A,1X,I0)') &
      'perturbed_constraint_count_fortran', fortran_constraint_count_perturbed, &
      'perturbed_constraint_count_cpp', cpp_constraint_count_perturbed
end program test_eval_objective_compare_probe
