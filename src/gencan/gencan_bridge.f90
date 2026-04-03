subroutine easygencan(n, x, l, u, m, lambda, rho, epsgpsn, maxit, maxfc, trtype, iprint, ncomp, &
   f, g, gpsupn, iter, fcnt, gcnt, cgcnt, inform, wi, wd, delmin)
   use iso_c_binding, only : c_int, c_double
   use compute_data, only : init1
   implicit none

   integer, intent(in) :: n, m, maxit, maxfc, trtype, iprint, ncomp
   integer, intent(out) :: iter, fcnt, gcnt, cgcnt, inform
   integer, intent(inout) :: wi(n)
   double precision, intent(inout) :: x(n), f, g(n), wd(8 * n), delmin
   double precision, intent(in) :: l(n), u(n), lambda(m), rho(m), epsgpsn
   double precision, intent(out) :: gpsupn
   integer(c_int) :: init1_flag

   interface
      subroutine packmol_gencan_set_init1_phase_c(flag) bind(C, name="packmol_gencan_set_init1_phase_c")
         import :: c_int
         integer(c_int), intent(in) :: flag
      end subroutine packmol_gencan_set_init1_phase_c

      subroutine packmol_gencan_easy_bridge(n, x, l, u, m, lambda, rho, epsgpsn, maxit, maxfc, trtype, &
         iprint, ncomp, f, g, gpsupn, iter, fcnt, gcnt, cgcnt, inform, wi, wd, delmin) &
         bind(C, name="packmol_gencan_easy_bridge")
         import :: c_int, c_double
         integer(c_int), intent(in) :: n, m, maxit, maxfc, trtype, iprint, ncomp
         integer(c_int), intent(out) :: iter, fcnt, gcnt, cgcnt, inform
         integer(c_int), intent(inout) :: wi(*)
         real(c_double), intent(inout) :: x(*), f, g(*), wd(*), delmin
         real(c_double), intent(in) :: l(*), u(*), lambda(*), rho(*), epsgpsn
         real(c_double), intent(out) :: gpsupn
      end subroutine packmol_gencan_easy_bridge
   end interface

   if (init1) then
      init1_flag = 1_c_int
   else
      init1_flag = 0_c_int
   end if
   call packmol_gencan_set_init1_phase_c(init1_flag)
   call packmol_gencan_easy_bridge(n, x, l, u, m, lambda, rho, epsgpsn, maxit, maxfc, trtype, &
      iprint, ncomp, f, g, gpsupn, iter, fcnt, gcnt, cgcnt, inform, wi, wd, delmin)
end subroutine easygencan

subroutine gencan(n, x, l, u, m, lambda, rho, epsgpen, epsgpsn, maxitnfp, epsnfp, maxitngp, fmin, maxit, &
   maxfc, udelta0, ucgmaxit, cgscre, cggpnf, cgepsi, cgepsf, epsnqmp, maxitnqmp, nearlyq, nint, next, &
   mininterp, maxextrap, gtype, htvtype, trtype, iprint, ncomp, f, g, gpeucn2, gpsupn, iter, fcnt, gcnt, &
   cgcnt, spgiter, spgfcnt, tniter, tnfcnt, tnstpcnt, tnintcnt, tnexgcnt, tnexbcnt, tnintfe, tnexgfe, &
   tnexbfe, inform, s, y, d, ind, lastgpns, w, eta, delmin, lspgma, lspgmi, theta, gamma, beta, sigma1, &
   sigma2, sterel, steabs, epsrel, epsabs, infrel, infabs)
   use iso_c_binding, only : c_int, c_double, c_bool
   use compute_data, only : init1
   implicit none

   logical, intent(in) :: nearlyq
   logical(c_bool) :: nearlyq_c
   integer(c_int) :: init1_flag
   integer, intent(in) :: n, m, maxitnfp, maxitngp, maxit, maxfc, ucgmaxit, cgscre, maxitnqmp, mininterp
   integer, intent(in) :: maxextrap, gtype, htvtype, trtype, iprint, ncomp
   integer, intent(inout) :: iter, fcnt, gcnt, cgcnt, spgiter, spgfcnt, tniter, tnfcnt, tnstpcnt, tnintcnt
   integer, intent(inout) :: tnexgcnt, tnexbcnt, tnintfe, tnexgfe, tnexbfe
   integer, intent(out) :: inform
   integer, intent(inout) :: ind(n)
   double precision, intent(inout) :: x(n), f, g(n), gpeucn2, gpsupn, s(n), y(n), d(n), w(5 * n)
   double precision, intent(inout) :: lastgpns(0:maxitngp-1)
   double precision, intent(in) :: l(n), u(n), lambda(m), rho(m), epsgpen, epsgpsn, epsnfp, fmin, udelta0
   double precision, intent(in) :: cggpnf, cgepsi, cgepsf, epsnqmp, nint, next, eta, delmin, lspgma, lspgmi
   double precision, intent(in) :: theta, gamma, beta, sigma1, sigma2, sterel, steabs, epsrel, epsabs
   double precision, intent(in) :: infrel, infabs

   interface
      subroutine packmol_gencan_set_init1_phase_c(flag) bind(C, name="packmol_gencan_set_init1_phase_c")
         import :: c_int
         integer(c_int), intent(in) :: flag
      end subroutine packmol_gencan_set_init1_phase_c

      subroutine packmol_gencan_gencan_bridge(n, x, l, u, m, lambda, rho, epsgpen, epsgpsn, maxitnfp, epsnfp, &
         maxitngp, fmin, maxit, maxfc, udelta0, ucgmaxit, cgscre, cggpnf, cgepsi, cgepsf, epsnqmp, maxitnqmp, &
         nearlyq, nint, next, mininterp, maxextrap, gtype, htvtype, trtype, iprint, ncomp, f, g, gpeucn2, &
         gpsupn, iter, fcnt, gcnt, cgcnt, spgiter, spgfcnt, tniter, tnfcnt, tnstpcnt, tnintcnt, tnexgcnt, &
         tnexbcnt, tnintfe, tnexgfe, tnexbfe, inform, s, y, d, ind, lastgpns, w, eta, delmin, lspgma, lspgmi, &
         theta, gamma, beta, sigma1, sigma2, sterel, steabs, epsrel, epsabs, infrel, infabs) &
         bind(C, name="packmol_gencan_gencan_bridge")
         import :: c_int, c_double, c_bool
         integer(c_int), intent(in) :: n, m, maxitnfp, maxitngp, maxit, maxfc, ucgmaxit, cgscre, maxitnqmp
         integer(c_int), intent(in) :: mininterp, maxextrap, gtype, htvtype, trtype, iprint, ncomp
         logical(c_bool), intent(in) :: nearlyq
         integer(c_int), intent(inout) :: iter, fcnt, gcnt, cgcnt, spgiter, spgfcnt, tniter, tnfcnt, tnstpcnt
         integer(c_int), intent(inout) :: tnintcnt, tnexgcnt, tnexbcnt, tnintfe, tnexgfe, tnexbfe
         integer(c_int), intent(out) :: inform
         integer(c_int), intent(inout) :: ind(*)
         real(c_double), intent(inout) :: x(*), f, g(*), gpeucn2, gpsupn, s(*), y(*), d(*), lastgpns(*), w(*)
         real(c_double), intent(in) :: l(*), u(*), lambda(*), rho(*), epsgpen, epsgpsn, epsnfp, fmin, udelta0
         real(c_double), intent(in) :: cggpnf, cgepsi, cgepsf, epsnqmp, nint, next, eta, delmin, lspgma, lspgmi
         real(c_double), intent(in) :: theta, gamma, beta, sigma1, sigma2, sterel, steabs, epsrel, epsabs
         real(c_double), intent(in) :: infrel, infabs
      end subroutine packmol_gencan_gencan_bridge
   end interface

   nearlyq_c = nearlyq
   if (init1) then
      init1_flag = 1_c_int
   else
      init1_flag = 0_c_int
   end if
   call packmol_gencan_set_init1_phase_c(init1_flag)

   call packmol_gencan_gencan_bridge(n, x, l, u, m, lambda, rho, epsgpen, epsgpsn, maxitnfp, epsnfp, maxitngp, &
      fmin, maxit, maxfc, udelta0, ucgmaxit, cgscre, cggpnf, cgepsi, cgepsf, epsnqmp, maxitnqmp, nearlyq_c, &
      nint, next, mininterp, maxextrap, gtype, htvtype, trtype, iprint, ncomp, f, g, gpeucn2, gpsupn, iter, &
      fcnt, gcnt, cgcnt, spgiter, spgfcnt, tniter, tnfcnt, tnstpcnt, tnintcnt, tnexgcnt, tnexbcnt, tnintfe, &
      tnexgfe, tnexbfe, inform, s, y, d, ind, lastgpns, w, eta, delmin, lspgma, lspgmi, theta, gamma, beta, &
      sigma1, sigma2, sterel, steabs, epsrel, epsabs, infrel, infabs)
end subroutine gencan

subroutine packmol_easyg_fortran_c(n, x, l, u, m, lambda, rho, epsgpsn, maxit, maxfc, trtype, iprint, ncomp, &
   f, g, gpsupn, iter, fcnt, gcnt, cgcnt, inform, wi, wd, delmin) bind(C, name="packmol_easyg_fortran_c")
   use iso_c_binding, only : c_int, c_double
   implicit none

   integer(c_int), intent(in) :: n, m, maxit, maxfc, trtype, iprint, ncomp
   integer(c_int), intent(out) :: iter, fcnt, gcnt, cgcnt, inform
   integer(c_int), intent(inout) :: wi(*)
   real(c_double), intent(inout) :: x(*), f, g(*), wd(*), delmin
   real(c_double), intent(in) :: l(*), u(*), lambda(*), rho(*), epsgpsn

   real(c_double), intent(out) :: gpsupn
   external :: easygcf

   call easygcf(n, x, l, u, m, lambda, rho, epsgpsn, maxit, maxfc, trtype, iprint, ncomp, f, g, gpsupn, iter, &
      fcnt, gcnt, cgcnt, inform, wi, wd, delmin)
end subroutine packmol_easyg_fortran_c

subroutine packmol_gencan_fortran_c(n, x, l, u, m, lambda, rho, epsgpen, epsgpsn, maxitnfp, epsnfp, maxitngp, &
   fmin, maxit, maxfc, udelta0, ucgmaxit, cgscre, cggpnf, cgepsi, cgepsf, epsnqmp, maxitnqmp, nearlyq, nint, &
   next, mininterp, maxextrap, gtype, htvtype, trtype, iprint, ncomp, f, g, gpeucn2, gpsupn, iter, fcnt, gcnt, &
   cgcnt, spgiter, spgfcnt, tniter, tnfcnt, tnstpcnt, tnintcnt, tnexgcnt, tnexbcnt, tnintfe, tnexgfe, tnexbfe, &
   inform, s, y, d, ind, lastgpns, w, eta, delmin, lspgma, lspgmi, theta, gamma, beta, sigma1, sigma2, sterel, &
   steabs, epsrel, epsabs, infrel, infabs) bind(C, name="packmol_gencan_fortran_c")
   use iso_c_binding, only : c_int, c_double, c_bool
   implicit none

   logical(c_bool), intent(in) :: nearlyq
   logical :: nearlyq_f
   integer(c_int), intent(in) :: n, m, maxitnfp, maxitngp, maxit, maxfc, ucgmaxit, cgscre, maxitnqmp
   integer(c_int), intent(in) :: mininterp, maxextrap, gtype, htvtype, trtype, iprint, ncomp
   integer(c_int), intent(inout) :: iter, fcnt, gcnt, cgcnt, spgiter, spgfcnt, tniter, tnfcnt, tnstpcnt
   integer(c_int), intent(inout) :: tnintcnt, tnexgcnt, tnexbcnt, tnintfe, tnexgfe, tnexbfe
   integer(c_int), intent(out) :: inform
   integer(c_int), intent(inout) :: ind(*)
   real(c_double), intent(inout) :: x(*), f, g(*), gpeucn2, gpsupn, s(*), y(*), d(*), lastgpns(*), w(*)
   real(c_double), intent(in) :: l(*), u(*), lambda(*), rho(*), epsgpen, epsgpsn, epsnfp, fmin, udelta0
   real(c_double), intent(in) :: cggpnf, cgepsi, cgepsf, epsnqmp, nint, next, eta, delmin, lspgma, lspgmi
   real(c_double), intent(in) :: theta, gamma, beta, sigma1, sigma2, sterel, steabs, epsrel, epsabs
   real(c_double), intent(in) :: infrel, infabs

   external :: gencf

   nearlyq_f = nearlyq

   call gencf(n, x, l, u, m, lambda, rho, epsgpen, epsgpsn, maxitnfp, epsnfp, maxitngp, fmin, maxit, maxfc, &
      udelta0, ucgmaxit, cgscre, cggpnf, cgepsi, cgepsf, epsnqmp, maxitnqmp, nearlyq_f, nint, next, mininterp, &
      maxextrap, gtype, htvtype, trtype, iprint, ncomp, f, g, gpeucn2, gpsupn, iter, fcnt, gcnt, cgcnt, spgiter, &
      spgfcnt, tniter, tnfcnt, tnstpcnt, tnintcnt, tnexgcnt, tnexbcnt, tnintfe, tnexgfe, tnexbfe, inform, s, y, &
      d, ind, lastgpns, w, eta, delmin, lspgma, lspgmi, theta, gamma, beta, sigma1, sigma2, sterel, steabs, &
      epsrel, epsabs, infrel, infabs)
end subroutine packmol_gencan_fortran_c
