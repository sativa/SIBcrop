!==================SUBROUTINE CROP_ACCUM=======================================
subroutine crop_accum(sib)

use kinds
use sibtype
use physical_parameters, only: tice    

implicit none

integer(kind=int_kind) :: i0
real(kind=dbl_kind)    :: temp_accum

!----------------------------------------------------------------------

type(sib_t), intent(inout) :: sib

!----------------------------------------------------------------------  

   temp_accum = 0.0_dbl_kind

   do i0 = 1, sib%diag%tb_indx
    
      temp_accum = temp_accum + sib%diag%tb_temp(i0)

   enddo

   sib%diag%ta_bar = temp_accum / float(sib%diag%tb_indx)


!itb_crop...GROWING DEGREE DAY; defined as a day with mean
!itb_crop...temperature (CAS) above 20C/293K

   if(sib%diag%ta_bar > 20.0_dbl_kind + tice) then

     sib%diag%gdd = sib%diag%gdd + sib%diag%ta_bar      &
                                 - 20.0_dbl_kind + tice

   endif


!itb_crop...reset things, now that accumulation has 
!itb_crop...taken place

   sib%diag%tb_indx = 0

   print*,'subroutine crop_accum'

end subroutine crop_accum
