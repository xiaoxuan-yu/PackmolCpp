module eval_state_bridge

   use iso_c_binding, only : c_bool, c_double, c_int, c_loc, c_null_ptr, c_ptr
   implicit none

contains

subroutine packmol_eval_meta_state_c( &
   maxrest_out, mrperatom_out, ntype_out, ntotmol_out, ntotat_out, nfixedat_out, ntype_with_fixed_out, lcellfirst_out, &
   ncells_out, scale_out, scale2_out, pbc_min_out, pbc_length_out, cell_length_out, &
   init1_out, move_out, fix_out, using_pbc_out) bind(C, name="packmol_eval_meta_state_c")
   use compute_data, only : cell_length, init1, lcellfirst, move, ncells, nfixedat, ntotat, ntotmol, ntype, &
                            scale, scale2
   use input, only : fix, ntype_with_fixed
   use pbc, only : pbc_length, pbc_min, using_pbc
   use sizes, only : maxrest, mrperatom
   implicit none

   integer(c_int), intent(out) :: maxrest_out, mrperatom_out, ntype_out, ntotmol_out, ntotat_out
   integer(c_int), intent(out) :: nfixedat_out, ntype_with_fixed_out, lcellfirst_out
   integer(c_int), intent(out) :: ncells_out(3)
   real(c_double), intent(out) :: scale_out, scale2_out
   real(c_double), intent(out) :: pbc_min_out(3), pbc_length_out(3), cell_length_out(3)
   logical(c_bool), intent(out) :: init1_out, move_out, fix_out, using_pbc_out

   maxrest_out = maxrest
   mrperatom_out = mrperatom
   ntype_out = ntype
   ntotmol_out = ntotmol
   ntotat_out = ntotat
   nfixedat_out = nfixedat
   ntype_with_fixed_out = ntype_with_fixed
   lcellfirst_out = lcellfirst
   ncells_out = ncells
   scale_out = scale
   scale2_out = scale2
   pbc_min_out = pbc_min
   pbc_length_out = pbc_length
   cell_length_out = cell_length
   init1_out = init1
   move_out = move
   fix_out = fix
   using_pbc_out = using_pbc
end subroutine packmol_eval_meta_state_c

subroutine packmol_eval_numeric_pointers_c( &
   nmols_ptr, natoms_ptr, idfirst_ptr, nratom_ptr, iratom_ptr, ityperest_ptr, ibmol_ptr, ibtype_ptr, &
   latomnext_ptr, latomfirst_ptr, lcellnext_ptr, xcart_ptr, coor_ptr, restpars_ptr, radius_ptr, &
   radius_ini_ptr, fscale_ptr, short_radius_ptr, short_radius_scale_ptr, gxcar_ptr, fdist_atom_ptr, &
   frest_atom_ptr) bind(C, name="packmol_eval_numeric_pointers_c")
   use compute_data, only : coor, fdist_atom, frest_atom, fscale, gxcar, ibmol, ibtype, idfirst, iratom, &
                            ityperest, latomfirst, latomnext, lcellnext, natoms, nmols, nratom, radius, &
                            radius_ini, restpars, short_radius, short_radius_scale, xcart
   implicit none

   type(c_ptr), intent(out) :: nmols_ptr, natoms_ptr, idfirst_ptr, nratom_ptr, iratom_ptr, ityperest_ptr
   type(c_ptr), intent(out) :: ibmol_ptr, ibtype_ptr, latomnext_ptr, latomfirst_ptr, lcellnext_ptr
   type(c_ptr), intent(out) :: xcart_ptr, coor_ptr, restpars_ptr, radius_ptr, radius_ini_ptr, fscale_ptr
   type(c_ptr), intent(out) :: short_radius_ptr, short_radius_scale_ptr, gxcar_ptr
   type(c_ptr), intent(out) :: fdist_atom_ptr, frest_atom_ptr

   nmols_ptr = c_null_ptr
   if (allocated(nmols)) nmols_ptr = c_loc(nmols(1))
   natoms_ptr = c_null_ptr
   if (allocated(natoms)) natoms_ptr = c_loc(natoms(1))
   idfirst_ptr = c_null_ptr
   if (allocated(idfirst)) idfirst_ptr = c_loc(idfirst(1))
   nratom_ptr = c_null_ptr
   if (allocated(nratom)) nratom_ptr = c_loc(nratom(1))
   iratom_ptr = c_null_ptr
   if (allocated(iratom)) iratom_ptr = c_loc(iratom(1,1))
   ityperest_ptr = c_null_ptr
   if (allocated(ityperest)) ityperest_ptr = c_loc(ityperest(1))
   ibmol_ptr = c_null_ptr
   if (allocated(ibmol)) ibmol_ptr = c_loc(ibmol(1))
   ibtype_ptr = c_null_ptr
   if (allocated(ibtype)) ibtype_ptr = c_loc(ibtype(1))
   latomnext_ptr = c_null_ptr
   if (allocated(latomnext)) latomnext_ptr = c_loc(latomnext(1))
   latomfirst_ptr = c_null_ptr
   if (allocated(latomfirst)) latomfirst_ptr = c_loc(latomfirst(1,1,1))
   lcellnext_ptr = c_null_ptr
   if (allocated(lcellnext)) lcellnext_ptr = c_loc(lcellnext(1))
   xcart_ptr = c_null_ptr
   if (allocated(xcart)) xcart_ptr = c_loc(xcart(1,1))
   coor_ptr = c_null_ptr
   if (allocated(coor)) coor_ptr = c_loc(coor(1,1))
   restpars_ptr = c_null_ptr
   if (allocated(restpars)) restpars_ptr = c_loc(restpars(1,1))
   radius_ptr = c_null_ptr
   if (allocated(radius)) radius_ptr = c_loc(radius(1))
   radius_ini_ptr = c_null_ptr
   if (allocated(radius_ini)) radius_ini_ptr = c_loc(radius_ini(1))
   fscale_ptr = c_null_ptr
   if (allocated(fscale)) fscale_ptr = c_loc(fscale(1))
   short_radius_ptr = c_null_ptr
   if (allocated(short_radius)) short_radius_ptr = c_loc(short_radius(1))
   short_radius_scale_ptr = c_null_ptr
   if (allocated(short_radius_scale)) short_radius_scale_ptr = c_loc(short_radius_scale(1))
   gxcar_ptr = c_null_ptr
   if (allocated(gxcar)) gxcar_ptr = c_loc(gxcar(1,1))
   fdist_atom_ptr = c_null_ptr
   if (allocated(fdist_atom)) fdist_atom_ptr = c_loc(fdist_atom(1))
   frest_atom_ptr = c_null_ptr
   if (allocated(frest_atom)) frest_atom_ptr = c_loc(frest_atom(1))
end subroutine packmol_eval_numeric_pointers_c

subroutine packmol_eval_copy_type_flags_c(comptype_out) bind(C, name="packmol_eval_copy_type_flags_c")
   use compute_data, only : comptype, ntype
   implicit none

   integer(c_int), intent(out) :: comptype_out(*)
   integer :: i

   do i = 1, ntype
      if (comptype(i)) then
         comptype_out(i) = 1_c_int
      else
         comptype_out(i) = 0_c_int
      end if
   end do
end subroutine packmol_eval_copy_type_flags_c

subroutine packmol_eval_copy_atom_flags_c(fixedatom_out, use_short_radius_out) bind(C, name="packmol_eval_copy_atom_flags_c")
   use compute_data, only : fixedatom, ntotat, use_short_radius
   implicit none

   integer(c_int), intent(out) :: fixedatom_out(*), use_short_radius_out(*)
   integer :: i

   do i = 1, ntotat
      if (fixedatom(i)) then
         fixedatom_out(i) = 1_c_int
      else
         fixedatom_out(i) = 0_c_int
      end if
      if (use_short_radius(i)) then
         use_short_radius_out(i) = 1_c_int
      else
         use_short_radius_out(i) = 0_c_int
      end if
   end do
end subroutine packmol_eval_copy_atom_flags_c

subroutine packmol_eval_set_lcellfirst_c(lcellfirst_in) bind(C, name="packmol_eval_set_lcellfirst_c")
   use compute_data, only : lcellfirst
   implicit none

   integer(c_int), intent(in) :: lcellfirst_in

   lcellfirst = lcellfirst_in
end subroutine packmol_eval_set_lcellfirst_c

end module eval_state_bridge
