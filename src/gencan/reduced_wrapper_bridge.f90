subroutine packmol_calchddiff_fortran_c(nind, ind, n, x, d, g, m, lambda, rho, gtype, hd, xtmp, sterel, steabs, inform) &
   bind(C, name="packmol_calchddiff_fortran_c")
   use iso_c_binding, only : c_int, c_double
   implicit none

   integer(c_int), intent(in) :: nind, n, m, gtype
   integer(c_int), intent(out) :: inform
   integer(c_int), intent(in) :: ind(*)
   real(c_double), intent(inout) :: x(*), d(*), g(*), hd(*), xtmp(*)
   real(c_double), intent(in) :: lambda(*), rho(*), sterel, steabs

   external :: calchddiff

   call calchddiff(nind, ind, x, d, g, n, x, m, lambda, rho, gtype, hd, xtmp, sterel, steabs, inform)
end subroutine packmol_calchddiff_fortran_c
