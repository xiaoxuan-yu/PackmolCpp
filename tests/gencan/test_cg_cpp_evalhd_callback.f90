program test_cg_cpp_evalhd_callback
   use iso_c_binding, only : c_int
   use compute_data
   use test_assertions
   implicit none

   interface
      integer(c_int) function packmol_gencan_impl_mode_c() bind(C, name="packmol_gencan_impl_mode_c")
         import :: c_int
      end function packmol_gencan_impl_mode_c

      subroutine packmol_gencan_reset_evalhd_call_count_c() bind(C, name="packmol_gencan_reset_evalhd_call_count_c")
      end subroutine packmol_gencan_reset_evalhd_call_count_c

      integer(c_int) function packmol_gencan_evalhd_call_count_c() bind(C, name="packmol_gencan_evalhd_call_count_c")
         import :: c_int
      end function packmol_gencan_evalhd_call_count_c
   end interface

   integer, parameter :: n = 2, m = 1
   integer :: ind(n), iter, rbdtype, rbdind, inform
   double precision :: x(n), g(n), l(n), u(n), s(n), w(n), y(n), r(n), d(n), sprev(n)
   double precision :: lambda(m), rho(m)
   external :: cg

   ntotmol = 0
   ntype = 0
   init1 = .true.
   move = .false.

   ind = (/ 1, 2 /)
   x = (/ 0.d0, 0.d0 /)
   g = (/ 1.d0, 0.d0 /)
   l = -1.d20
   u = 1.d20
   s = 0.d0
   w = 0.d0
   y = 0.d0
   r = 0.d0
   d = 0.d0
   sprev = 0.d0
   lambda = 0.d0
   rho = 0.d0

   call assert_equal_int(packmol_gencan_impl_mode_c(), 1, "test must run in cpp impl mode")
   call packmol_gencan_reset_evalhd_call_count_c()

   call cg(n, ind, n, x, m, lambda, rho, g, 1.d0, l, u, 1.d-4, 1.d-4, 5, 5, .false., 0, 0, 1, 0, n, s, iter, &
      rbdtype, rbdind, inform, w, y, r, d, sprev, 1.d-6, 1.d-7, 1.d-10, 1.d-10, 1.d-20, 1.d20, 1.d99)

   call assert_true(packmol_gencan_evalhd_call_count_c() > 0, "cpp cg with htvtype=0 should call evalhd wrapper")
   call assert_true(inform >= 0, "cpp cg with htvtype=0 should not fail in evalhd path")
end program test_cg_cpp_evalhd_callback
