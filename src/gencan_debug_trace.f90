subroutine packmol_gencan_debug_enabled_f(enabled)
   implicit none
   logical, intent(out) :: enabled
   character(len=32) :: value
   integer :: status

   value = ''
   call get_environment_variable('PACKMOL_GENCAN_DEBUG', value, status=status)
   enabled = .false.
   if (status == 0) then
      if (trim(value) == '1' .or. trim(value) == 'true' .or. trim(value) == 'TRUE' .or. &
          trim(value) == 'on' .or. trim(value) == 'ON' .or. trim(value) == 'yes' .or. trim(value) == 'YES') then
         enabled = .true.
      end if
   end if
end subroutine packmol_gencan_debug_enabled_f

subroutine packmol_gencan_debug_trace_iter(tag, f, eps, fdist, frest, gpsupn, gpeucn2, gieucn2, xnorm, spg_test)
   implicit none
   character(len=*), intent(in) :: tag
   double precision, intent(in) :: f, eps, fdist, frest, gpsupn, gpeucn2, gieucn2, xnorm
   logical, intent(in) :: spg_test

   write(*,'(A,1X,A,1X,A,1PE24.16,1X,A,1PE24.16,1X,A,1PE24.16,1X,A,1PE24.16,1X,A,1PE24.16,1X,A,1PE24.16,1X,A,1PE24.16,1X,A,1PE24.16,1X,A,1X,I0)') &
      '[gencan-fortran]', trim(tag), 'f=', f, 'precision=', eps, 'fdist=', fdist, 'frest=', frest, &
      'gpsupn=', gpsupn, 'gpeucn2=', gpeucn2, 'gieucn2=', gieucn2, 'xnorm=', xnorm, 'spg_test=', merge(1,0,spg_test)
end subroutine packmol_gencan_debug_trace_iter

subroutine packmol_gencan_debug_trace_tn(iter, delta, cgeps, cgmaxit, kappa, gpsupn, gpeucn2, xnorm)
   implicit none
   integer, intent(in) :: iter, cgmaxit
   double precision, intent(in) :: delta, cgeps, kappa, gpsupn, gpeucn2, xnorm

   write(*,'(A,1X,A,1X,I0,1X,A,1PE24.16,1X,A,1PE24.16,1X,A,1X,I0,1X,A,1PE24.16,1X,A,1PE24.16,1X,A,1PE24.16,1X,A,1PE24.16)') &
      '[gencan-fortran]', 'tn-setup iter=', iter, 'delta=', delta, 'cgeps=', cgeps, 'cgmaxit=', cgmaxit, &
      'kappa=', kappa, 'gpsupn=', gpsupn, 'gpeucn2=', gpeucn2, 'xnorm=', xnorm
end subroutine packmol_gencan_debug_trace_tn

subroutine packmol_gencan_debug_trace_tn_result(tag, iter, inform, cgiter, rbdtype, rbdind, f)
   implicit none
   character(len=*), intent(in) :: tag
   integer, intent(in) :: iter, inform, cgiter, rbdtype, rbdind
   double precision, intent(in) :: f

   write(*,'(A,1X,A,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,I0,1X,A,1PE24.16)') &
      '[gencan-fortran]', trim(tag), 'iter=', iter, 'inform=', inform, 'cg_iter=', cgiter, &
      'rbdtype=', rbdtype, 'rbdind=', rbdind, 'f=', f
end subroutine packmol_gencan_debug_trace_tn_result

subroutine packmol_gencan_debug_trace_initial(tag, itype, loop_idx, fx, frest)
   implicit none
   character(len=*), intent(in) :: tag
   integer, intent(in) :: itype, loop_idx
   double precision, intent(in) :: fx, frest

   write(*,'(A,1X,A,1X,A,1X,I0,1X,A,1X,I0,1X,A,1PE24.16,1X,A,1PE24.16)') &
      '[gencan-initial]', trim(tag), 'itype=', itype, 'loop=', loop_idx, 'fx=', fx, 'frest=', frest
end subroutine packmol_gencan_debug_trace_initial
