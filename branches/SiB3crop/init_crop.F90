!Subroutine to initialize variables for the crop model
!kdcorbin, 01/11

subroutine init_crop(sib,time)

use sibtype
use timetype
use sib_bc_module
use sib_const_module
use sib_io_module, only: drvr_type

implicit none

!input
type(sib_t), dimension(subcount), intent(inout) :: sib
type(time_struct),intent(inout) :: time

!local variables
integer :: i

!time dependant, output variables
!kdcorbin, 02/11 - added scatp, scatg, park
type time_dep_var
     real (kind=real_kind) :: fPAR    ! Canopy absorbed fraction of PAR
     real (kind=real_kind) :: LAI     ! Leaf-area index
     real (kind=real_kind) :: Green   ! Canopy greeness fraction of LAI
     real (kind=real_kind) :: zo      ! Canopy roughness scheme
     real (kind=real_kind) :: zp_disp ! Zero plane displacement
     real (kind=real_kind) :: RbC     ! RB Coefficient (c1)
     real (kind=real_kind) :: RdC     ! RC Coefficient (c2)
     real (kind=real_kind) :: scatp   ! Canopy transmittance + reflectance of PAR
     real (kind=real_kind) :: scatg   ! Ground transmittance + reflectance of PAR
     real (kind=real_kind) :: park     ! Mean canopy absorption optical depth wrt PAR
     real (kind=real_kind) :: gmudmu  ! Time-mean leaf projection
end type time_dep_var

type(time_dep_var) TimeVar

type(aero_var),dimension(50,50) :: tempaerovar
real(kind=real_kind),dimension(2,2) :: temptran,tempref
integer(kind=int_kind) :: biomeclass
real(kind=real_kind) :: crop_doy

do i=1,subcount

   if (sib(i)%param%biome >= 20.) then

       !kdcorbin, 09/09 - change temporary biome class depending on 
       !  if crops exist or if fallow
       if (sib(i)%diag%phen_switch == 0) then
           biomeclass = temp_biome_bare
       else
           biomeclass = temp_biome_crop
       endif

       tempaerovar = aerovar(:,:,biomeclass)
 
       temptran(1,1) = sib(i)%param%tran(1,1)
       temptran(1,2) = sib(i)%param%tran(1,2)
       temptran(2,1) = sib(i)%param%tran(2,1)
       temptran(2,2) = sib(i)%param%tran(2,2)

       tempref(1,1) = sib(i)%param%ref(1,1)
       tempref(1,2) = sib(i)%param%ref(1,2)
       tempref(2,1) = sib(i)%param%ref(2,1)
       tempref(2,2) = sib(i)%param%ref(2,2)

      !kdcorbin, 02/11
       crop_doy = time%doy

       if (sib(i)%diag%phen_switch == 0) then  !FALLOW
           call mapper(          &
                latsib(i),        &
                crop_doy,     &
                min_ndvi_crop,      &
                min_ndvi_crop,      &
                min_fvcov_crop,    &
                sib(i)%param%chil,   &
                temptran,         &
                tempref,          &
                morphtab(biomeclass),      &
                tempaerovar,      &
                laigrid,          &
                fvcovergrid,      &
                timevar)

                !kdcorbin, 07/09 
                sib(i)%diag%phen_lai = min_lai_crop

        elseif(sib(i)%diag%phen_switch == 1) then !CROP
              
            call leaf_weight(sib,time,i)

            sib(i)%diag%phen_lai = sib(i)%diag%leafwt_c*2*0.02
            timevar%lai = sib(i)%diag%phen_lai

           !kdcorbin, 02/11 - changed time from pmonth
            call phen_mapper(   &
                 latsib(i),  &
                 crop_doy,  &
                 sib(i)%param%vcover,  &
                 sib(i)%param%chil,    &
                 temptran,   &
                 tempref,    &
                 morphtab(biomeclass),   &
                 tempaerovar,  &
                 laigrid,   &
                 fvcovergrid,  &
                 timevar)

           sib(i)%param%aparc1 = timevar%fpar
           sib(i)%param%zlt1 = timevar%lai
           sib(i)%param%green1 = timevar%green
           sib(i)%param%z0d1 = timevar%zo
           sib(i)%param%zp_disp1 = timevar%zp_disp
           sib(i)%param%rbc1 = timevar%rbc
           sib(i)%param%rdc1 = timevar%rdc
           sib(i)%param%gmudmu1 = timevar%gmudmu

       else

           print*,'phen_switch=',sib%diag%phen_switch
           stop'INCORRECT VALUE FOR PHENLOLOGY SWITCH'

       endif

           sib(i)%param%aparc2 = timevar%fpar
           sib(i)%param%zlt2 = timevar%lai
           sib(i)%param%green2 = timevar%green
           sib(i)%param%z0d2 = timevar%zo
           sib(i)%param%zp_disp2 = timevar%zp_disp
           sib(i)%param%rbc2 = timevar%rbc
           sib(i)%param%rdc2 = timevar%rdc
           sib(i)%param%gmudmu2 = timevar%gmudmu
  
           sib(i)%param%scatp = timevar%scatp
           sib(i)%param%scatg = timevar%scatg
           sib(i)%param%park = timevar%park

   endif  !crop variable biome class

   !Set initial parameters for crops - kdcorbin, 02/11
   !!!if (drvr_type .ne. 'single') then    
        call set_ti(sib(i))
   !!!endif

enddo  !subcount

end



