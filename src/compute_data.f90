!
!  Written by Leandro Martínez, 2009-2011.
!  Copyright (c) 2009-2018, Leandro Martínez, Jose Mario Martinez,
!  Ernesto G. Birgin.
!
module compute_data

   use sizes
   use iso_c_binding, only : c_double, c_int

   integer :: ntotmol, ntype, nfixedat, ntotat
   integer :: ncells(3)

   integer, allocatable, target :: nmols(:) ! (ntype)
   integer, allocatable, target :: natoms(:) ! (ntype)
   integer, allocatable, target :: idfirst(:) ! (ntype)
   integer, allocatable, target :: nratom(:) ! (ntotat)
   integer, allocatable, target :: iratom(:,:) ! (ntotat,mrperatom)
   integer, allocatable, target :: ityperest(:) ! (maxrest)
   integer, allocatable, target :: ibmol(:) ! (ntotat)
   integer, allocatable, target :: ibtype(:) ! (ntotat)

   double precision :: scale, scale2
   real(c_double), bind(C, name="packmol_fdist") :: fdist
   real(c_double), bind(C, name="packmol_frest") :: frest
   real(c_double), bind(C, name="packmol_pair_penalty_sum") :: pair_penalty_sum
   real(c_double), bind(C, name="packmol_constraint_penalty_sum") :: constraint_penalty_sum
   integer(c_int), bind(C, name="packmol_pair_penalty_count") :: pair_penalty_count
   integer(c_int), bind(C, name="packmol_constraint_penalty_count") :: constraint_penalty_count
   double precision :: sizemin(3), sizemax(3)
   double precision :: cell_length(3), system_length(3)
   double precision :: radmax

   double precision, allocatable, target :: xcart(:,:) ! (ntotat,3)
   double precision, allocatable, target :: coor(:,:) ! (ntotat,3)
   double precision, allocatable, target :: restpars(:,:) ! (maxrest,9)
   double precision, allocatable :: rot_bound(:,:,:) ! (ntype,3,2)
   double precision, allocatable, target :: radius(:), radius_ini(:), fscale(:) ! (ntotat)
   double precision, allocatable, target :: short_radius(:), short_radius_scale(:) ! ntotat
   double precision, allocatable, target :: gxcar(:,:) ! (ntotat,3)

   double precision, allocatable, target :: fdist_atom(:), frest_atom(:) ! (ntotat)
   double precision, allocatable :: dmax(:) ! (ntype)

   logical, allocatable :: constrain_rot(:,:) ! (ntype,3)
   logical, allocatable :: comptype(:) ! (ntype)
   logical, allocatable :: fixedatom(:) ! (ntotat)
   logical, allocatable :: use_short_radius(:) ! ntotat
   logical :: init1, move

   ! For linked lists
   integer, allocatable, target :: latomnext(:) ! (ntotat)
   integer, allocatable, target :: latomfirst(:,:,:) !  (ncells(1),ncells(2),ncells3))
   integer, allocatable :: latomfix(:,:,:) ! (ncells(1),ncells(2),ncells(3))

   ! For movebad
   double precision, allocatable :: fmol(:), radiuswork(:) ! (ntotat)

   ! For restmol
   double precision, allocatable :: xmol(:) ! (nn)
   logical, allocatable :: compsafe(:) ! (ntype)

   ! For cells with atoms linked lists
   integer :: lcellfirst
   integer, allocatable, target :: lcellnext(:) ! (ncells(1)*ncells(2)*ncells(3))
   logical, allocatable :: empty_cell(:,:,:) ! (ncells(1),ncells(2),ncells(3))

end module compute_data
