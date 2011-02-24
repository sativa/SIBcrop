subroutine calculate_td (sib, midmon, curNDVI)

! Calls mapper for each sib point.  Mapper calculates the time dependent
!   boundary condition data.  This data is then assigned to appropriate 
!   variables.

!     CREATED BY:
!        Owen Leonard   August 10, 2001
!     MODIFICATIONS:
!     SUBROUTINES CALLED:
!     FUNCTIONS CALLED:

use sibtype
use sib_const_module
use sib_bc_module
use kinds

implicit none

! parameters
type(sib_t), dimension(subcount), intent(inout) :: sib
real(kind=real_kind), intent (in) :: midmon     ! middle of month
!ndvi variables (local)
real(kind=real_kind), dimension(subcount), intent(in) :: curndvi

! local variables
integer(kind=int_kind) :: i,k   ! index variable
real(kind=real_kind) :: temptran (2,2)
real(kind=real_kind) :: tempref (2,2)

! kdcorbin, 02/11 - added scatp, scatg, and park
type time_dep_var
    real(kind=real_kind) :: fpar       ! canopy absorbed fraction of par
    real(kind=real_kind) :: lai        ! leaf-area index
    real(kind=real_kind) :: green      ! canopy greeness fraction of lai
    real(kind=real_kind) :: zo         ! canopy roughness coeff 
    real(kind=real_kind) :: zp_disp    ! zero plane displacement
    real(kind=real_kind) :: rbc        ! rb coefficient (c1)
    real(kind=real_kind) :: rdc        ! rc coefficient (c2)
    real (kind=real_kind) :: scatp   ! Canopy transmittance + reflectance of PAR
    real (kind=real_kind) :: scatg   ! Ground transmittance + reflectance of PAR
    real (kind=real_kind) :: park     ! Mean canopy absorption optical depth wrt PAR
    real(kind=real_kind) :: gmudmu     ! time-mean leaf projection
end type time_dep_var
type(time_dep_var) :: timevar

type(aero_var) :: tempaerovar(50,50)



    do i = 1, subcount 

        k = int(sib(i)%param%biome)

!itb_crop...
        if(sib(i)%param%biome >= 20.0) k = 12
!itb_crop...

        temptran = sib(i)%param%tran(:,:)
        tempref = sib(i)%param%ref(:,:)
        tempaerovar = aerovar(:,:,k)

        call mapper(          &
            latsib(subset(i)),&
            midmon,           &
            prevndvi(i),      &
            curndvi(i),       &
            sib(i)%param%vcover, &
            sib(i)%param%chil,   &
            temptran,         &
            tempref,          &
            morphtab(k),      &
            tempaerovar,      &
            laigrid,          &
            fvcovergrid,      &
            timevar)

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

    enddo      

end subroutine calculate_td
