!=========================================
subroutine set_ti (sibpt)
!=========================================
!Sets time invariant parameters for sib from
!  TI tables (biovart3 & biovar4).
!
!kdcorbin, 02/11

use sibtype
use sib_const_module
use sib_io_module, only: drvr_type

implicit none

type(sib_t), intent(inout) :: sibpt

integer biometype,cinttype
integer pt

!------------------------
biometype=sibpt%param%biome
if (biometype >= 20) then
   if (sibpt%diag%phen_switch == 1) then
      biometype=temp_biome_crop
   else 
      biometype=temp_biome_bare
   endif
endif

do pt=1,phys  !phystype loop
   if (sibpt%param%phystype(pt) == 3) then
      ! C3 Table
      sibpt%param%z2 = biovart3(biometype,1)
      sibpt%param%z1 = biovart3(biometype,2)
      sibpt%param%chil = biovart3(biometype,4)
      sibpt%param%phc = biovart3(biometype,7)
      sibpt%param%tran(1,1) = biovart3(biometype,8)
      sibpt%param%tran(2,1) = biovart3(biometype,9)
      sibpt%param%tran(1,2) = biovart3(biometype,10)
      sibpt%param%tran(2,2) = biovart3(biometype,11)
      sibpt%param%ref(1,1) = biovart3(biometype,12)
      sibpt%param%ref(2,1) = biovart3(biometype,13)
      sibpt%param%ref(1,2) = biovart3(biometype,14)
      sibpt%param%ref(2,2) = biovart3(biometype,15)
      sibpt%param%vmax0(1) = biovart3(biometype,16)
      sibpt%param%effcon(1) = biovart3(biometype,17)
      sibpt%param%gradm(1) = biovart3(biometype,18)
      sibpt%param%binter(1) = biovart3(biometype,19)
      sibpt%param%atheta(1) = biovart3(biometype,20)
      sibpt%param%btheta(1) = biovart3(biometype,21)
      sibpt%param%trda(1) = biovart3(biometype,22)
      sibpt%param%trdm(1) = biovart3(biometype,23)
      sibpt%param%trop(1) = biovart3(biometype,24)
      sibpt%param%respcp(1) = biovart3(biometype,25)
      sibpt%param%slti(1) = biovart3(biometype,26)
      sibpt%param%hltii(1) = biovart3(biometype,27)
      sibpt%param%shti(1) = biovart3(biometype,28)
      sibpt%param%hhti(1) = biovart3(biometype,29)
   elseif (sibpt%param%phystype(pt) == 4) then
      ! C4 Table
      sibpt%param%z2 = biovart4(biometype,1)
      sibpt%param%z1 = biovart4(biometype,2)
      sibpt%param%chil = biovart4(biometype,4)
      sibpt%param%phc = biovart4(biometype,7)
      sibpt%param%tran(1,1) = biovart4(biometype,8)
      sibpt%param%tran(2,1) = biovart4(biometype,9)
      sibpt%param%tran(1,2) = biovart4(biometype,10)
      sibpt%param%tran(2,2) = biovart4(biometype,11)
      sibpt%param%ref(1,1) = biovart4(biometype,12)
      sibpt%param%ref(2,1) = biovart4(biometype,13)
      sibpt%param%ref(1,2) = biovart4(biometype,14)
      sibpt%param%ref(2,2) = biovart4(biometype,15)
      sibpt%param%vmax0(2) = biovart4(biometype,16)
      sibpt%param%effcon(2) = biovart4(biometype,17)
      sibpt%param%gradm(2) = biovart4(biometype,18)
      sibpt%param%binter(2) = biovart4(biometype,19)
      sibpt%param%atheta(2) = biovart4(biometype,20)
      sibpt%param%btheta(2) = biovart4(biometype,21)
      sibpt%param%trda(2) = biovart4(biometype,22)
      sibpt%param%trdm(2) = biovart4(biometype,23)
      sibpt%param%trop(2) = biovart4(biometype,24)
      sibpt%param%respcp(2) = biovart4(biometype,25)
      sibpt%param%slti(2) = biovart4(biometype,26)
      sibpt%param%hltii(2) = biovart4(biometype,27)
      sibpt%param%shti(2) = biovart4(biometype,28)
      sibpt%param%hhti(2) = biovart4(biometype,29)
   elseif (sibpt%param%phystype(pt) .ne. 0) then
      stop 'WE DO NOT HAVE PHYSIOLOGY OTHER THAN C3/C4 YET'
   endif
enddo

!Now change from basic variables to crop values, but only
!  for certain parameters
biometype=sibpt%param%biome-19

if (biometype > 0) then
     cinttype=crop_cint(biometype)
      !set physfrac - kdcorbin, 02/11
      if ( drvr_type .ne. 'single' ) then
          sibpt%param%physfrac(1) = crop_physfrac1(biometype)
          sibpt%param%physfrac(2) = crop_physfrac2(biometype)
          sibpt%param%physfrac1 = sibpt%param%physfrac
          sibpt%param%physfrac2 = sibpt%param%physfrac
      endif

      if (sibpt%diag%phen_switch > 0) then
          sibpt%param%soref(1) = crop_soref1(biometype)
          sibpt%param%soref(2) = crop_soref2(biometype)
          sibpt%param%vcover = crop_vcover(biometype)
          sibpt%param%z2 = crop_z2(biometype)
          sibpt%param%vmax0(cinttype) = crop_vmax0a(biometype)
          sibpt%param%effcon(cinttype) = crop_effcon(biometype)
          sibpt%param%gradm(cinttype) = crop_gradm(biometype)
          sibpt%param%binter(cinttype) = crop_binter(biometype)
          sibpt%param%atheta(cinttype) = crop_atheta(biometype)
          sibpt%param%btheta(cinttype) = crop_btheta(biometype)
          sibpt%param%respcp(cinttype) = crop_respcp(biometype)
          sibpt%param%slti(cinttype) = crop_slti(biometype)
          sibpt%param%hltii(cinttype) = crop_hltii(biometype)
          sibpt%param%shti(cinttype) = crop_shti(biometype)
          sibpt%param%hhti(cinttype) = crop_hhti(biometype)
      else
          sibpt%param%soref(1) = bare_soref1
          sibpt%param%soref(2) = bare_soref2
          sibpt%param%vcover = bare_vcover
      endif
endif

end subroutine set_ti

