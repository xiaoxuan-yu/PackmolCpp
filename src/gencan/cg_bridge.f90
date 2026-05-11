subroutine cg(nind, ind, n, x, m, lambda, rho, g, delta, l, u, eps, epsnqmp, &
   maxitnqmp, maxit, nearlyq, gtype, htvtype, trtype, iprint, ncomp, s, iter, &
   rbdtype, rbdind, inform, w, y, r, d, sprev, theta, sterel, steabs, epsrel, &
   epsabs, infrel, infabs)
   use iso_c_binding, only : c_int, c_double, c_bool
   implicit none

   integer, intent(in) :: nind, n, m, maxitnqmp, maxit, gtype, htvtype, trtype, iprint, ncomp
   logical, intent(in) :: nearlyq
   logical(c_bool) :: nearlyq_c
   integer, intent(inout) :: iter
   integer, intent(out) :: rbdtype, rbdind, inform
   integer, intent(in) :: ind(nind)
   double precision, intent(inout) :: x(n), g(n)
   double precision, intent(in) :: lambda(m), rho(m), delta, l(n), u(n), eps, epsnqmp
   double precision, intent(inout) :: s(n), w(n), y(n), r(n), d(n), sprev(n)
   double precision, intent(in) :: theta, sterel, steabs, epsrel, epsabs, infrel, infabs

   interface
      subroutine packmol_gencan_cg_bridge(nind, ind, n, x, m, lambda, rho, g, delta, l, u, eps, epsnqmp, &
         maxitnqmp, maxit, nearlyq, gtype, htvtype, trtype, iprint, ncomp, s, iter, rbdtype, rbdind, inform, &
         w, y, r, d, sprev, theta, sterel, steabs, epsrel, epsabs, infrel, infabs) bind(C, name="packmol_gencan_cg_bridge")
         import :: c_int, c_double, c_bool
         integer(c_int), intent(in) :: nind, n, m, maxitnqmp, maxit, gtype, htvtype, trtype, iprint, ncomp
         logical(c_bool), intent(in) :: nearlyq
         integer(c_int), intent(inout) :: iter
         integer(c_int), intent(out) :: rbdtype, rbdind, inform
         integer(c_int), intent(in) :: ind(*)
         real(c_double), intent(inout) :: x(*), g(*)
         real(c_double), intent(in) :: lambda(*), rho(*), delta, l(*), u(*), eps, epsnqmp
         real(c_double), intent(inout) :: s(*), w(*), y(*), r(*), d(*), sprev(*)
         real(c_double), intent(in) :: theta, sterel, steabs, epsrel, epsabs, infrel, infabs
      end subroutine packmol_gencan_cg_bridge
   end interface

   nearlyq_c = nearlyq

   call packmol_gencan_cg_bridge(nind, ind, n, x, m, lambda, rho, g, delta, l, u, eps, epsnqmp, &
      maxitnqmp, maxit, nearlyq_c, gtype, htvtype, trtype, iprint, ncomp, s, iter, rbdtype, rbdind, inform, &
      w, y, r, d, sprev, theta, sterel, steabs, epsrel, epsabs, infrel, infabs)
end subroutine cg

subroutine packmol_cg_fortran_c(nind, ind, n, x, m, lambda, rho, g, delta, l, u, eps, epsnqmp, &
   maxitnqmp, maxit, nearlyq, gtype, htvtype, trtype, iprint, ncomp, s, iter, rbdtype, rbdind, inform, &
   w, y, r, d, sprev, theta, sterel, steabs, epsrel, epsabs, infrel, infabs) bind(C, name="packmol_cg_fortran_c")
   use iso_c_binding, only : c_int, c_double, c_bool
   implicit none

   integer(c_int), intent(in) :: nind, n, m, maxitnqmp, maxit, gtype, htvtype, trtype, iprint, ncomp
   logical(c_bool), intent(in) :: nearlyq
   logical :: nearlyq_f
   integer(c_int), intent(inout) :: iter
   integer(c_int), intent(out) :: rbdtype, rbdind, inform
   integer(c_int), intent(in) :: ind(*)
   real(c_double), intent(inout) :: x(*), g(*)
   real(c_double), intent(in) :: lambda(*), rho(*), delta, l(*), u(*), eps, epsnqmp
   real(c_double), intent(inout) :: s(*), w(*), y(*), r(*), d(*), sprev(*)
   real(c_double), intent(in) :: theta, sterel, steabs, epsrel, epsabs, infrel, infabs

   external :: cgf

   nearlyq_f = nearlyq

   call cgf(nind, ind, n, x, m, lambda, rho, g, delta, l, u, eps, epsnqmp, maxitnqmp, maxit, nearlyq_f, &
      gtype, htvtype, trtype, iprint, ncomp, s, iter, rbdtype, rbdind, inform, w, y, r, d, sprev, theta, &
      sterel, steabs, epsrel, epsabs, infrel, infabs)
end subroutine packmol_cg_fortran_c
