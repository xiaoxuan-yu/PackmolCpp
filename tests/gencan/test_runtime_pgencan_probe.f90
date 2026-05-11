program test_runtime_pgencan_probe
   use iso_c_binding, only : c_int
   use compute_data
   use input
   use test_runtime_setup
   use usegencan, only : maxit
   implicit none

   interface
      integer(c_int) function packmol_gencan_impl_mode_c() bind(C, name="packmol_gencan_impl_mode_c")
         import :: c_int
      end function packmol_gencan_impl_mode_c
   end interface

   integer :: n, mode_id, i
   double precision :: fx, fx_before, fx_check, dx_max, dx_rms
   double precision, allocatable :: x(:), x0(:)
   character(len=1000) :: input_path
   external :: initial, pgencan, computef

   call get_command_argument(1, input_path)
   call load_runtime_problem(trim(input_path), n, x)
   allocate(x0(n))
   call initial(n, x)
   x0 = x
   call computef(n, x, fx_before)
   call pgencan(n, x, fx)

   call computef(n, x, fx_check)
   dx_max = 0.d0
   dx_rms = 0.d0
   do i = 1, n
      dx_max = max(dx_max, abs(x(i) - x0(i)))
      dx_rms = dx_rms + (x(i) - x0(i))**2
   end do
   dx_rms = sqrt(dx_rms / dble(max(1, n)))

   mode_id = packmol_gencan_impl_mode_c()
   write(*,'(A,1X,I0,1X,A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,ES24.16)') &
      'mode', mode_id, 'fx_before', fx_before, 'fx', fx, 'fx_check', fx_check, 'fdist', fdist, 'frest', frest
   write(*,'(A,1X,ES24.16,1X,A,1X,ES24.16,1X,A,1X,I0)') &
      'dx_max', dx_max, 'dx_rms', dx_rms, 'maxit', maxit
end program test_runtime_pgencan_probe
