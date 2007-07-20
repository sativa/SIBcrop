!==================SUBROUTINE CROP_ACCUM=======================================
subroutine crop_accum(sib)

use kinds
use sibtype
use timetype
use physical_parameters, only: tice    

implicit none

integer(kind=int_kind) :: i0,n,yr_index,year,doy,ndf60,pd
real(kind=dbl_kind)    :: temp_accum,ta_bar
real, allocatable	   :: tempf(:)

!----------------------------------------------------------------------

type(sib_t), intent(inout) :: sib

!----------------------------------------------------------------------  

   temp_accum = 0.0_dbl_kind

   do i0 = 1, sib%diag%tb_indx
    
      temp_accum = temp_accum + sib%diag%tb_temp(i0)
!print*,i0,sib%diag%tb_temp(i0),temp_accum
!pause
!      print*,i0,sib%diag%tb_temp(i0),temp_accum
 

   enddo

   sib%diag%ta_bar = temp_accum / float(sib%diag%tb_indx)

!conversion of avg. daily temperature (ta_bar) from Kelvin to Fahrenheit- EL
allocate(tempf(365))
tempf(sib%diag%doy)=((sib%diag%ta_bar-273.15)*1.8)+32.0	


!itb_crop...GROWING DEGREE DAY; defined as a day with mean
!itb_crop...temperature (CAS) above 20C/293K

!Calling for different phenology schemes based on the year and the crop- EL
	if(mod(sib%diag%year,2)==1) then  
!		call soy_phen
!	else
		call corn_phen	
	endif
!print*,sib%diag%year,tempf(sib%diag%doy),sib%diag%doy

contains

!------------------------------------------------------------------------------
subroutine corn_phen
!------------------------------------------------------------------------------

!--------------------------
!Calculate the planting date
!---------------------------

if (tempf(sib%diag%doy)<60.0) then
    ndf60=0			!ndf60= no. of days withe avg. temperature above 60F
elseif (tempf(sib%diag%doy)>=60.0) then
    ndf60=ndf60+1
endif

if (ndf60==10) then

    pd=sib%diag%doy


    pd=sib%diag%doy

endif

print*,pd


!----------------------------
!Calculate growing degree days
!-----------------------------


	if (pd>0 							.and.		& !this line was added to avoid gdd calculation before real planting date, since pd is printed out as 0 before the real planting date based on the above ndf60==10 criterion
		sib%diag%doy >= pd            .and.          &
        tempf(sib%diag%doy)>50.0      .and.          &
        tempf(sib%diag%doy)<86.0)      then

    	sib%diag%gdd=sib%diag%gdd + tempf(sib%diag%doy)- 50.0_dbl_kind


   	sib%diag%gdd=sib%diag%gdd + tempf(sib%diag%doy)- 50.0_dbl_kind
	
	endif

!   if(sib%diag%ta_bar > 20.0_dbl_kind + tice) then

!     sib%diag%gdd = sib%diag%gdd + sib%diag%ta_bar      &
!                                 - 20.0_dbl_kind + tice
!   endif
print*,pd,tempf(sib%diag%doy),sib%diag%gdd





!itb_crop...reset things, now that accumulation has 
!itb_crop...taken place



   sib%diag%tb_indx = 0

 print*,'subroutine crop_accum' 

end subroutine corn_phen

 

end subroutine crop_accum
