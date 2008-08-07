!================ SUBROUTINE LEAF WEIGHT ===============


!itb_crop...Short program to determine final leaf weight,
!itb_crop...based on crop, time since planting, and 
!itb_crop...cumulative dry weight

   subroutine leaf_weight(sib,time)

use kinds
use sibtype
use timetype

!------------------------------------------------------------------
type(sib_t), intent(inout) :: sib
type(time_struct), intent(in) :: time
!------------------------------------------------------------------


!itb_crop...CORN

   if(sib%param%biome == 20.0) then

      if( sib%diag%gdd == 0.0 ) then

           sib%diag%leafwt_c =  0.01 
      
      elseif (sib%diag%gdd > 0.0001 .and. sib%diag%gdd < 2500) then

	       sib%diag%leafwt_c = sib%diag%cum_drywt(2)
	
      elseif (sib%diag%gdd >= 2650.0 .and. sib%diag%gdd < 2900.0) then

          sib%diag%leafwt_c = 0.95 * sib%diag%cum_drywt(2) -   &
          (0.95-0.1) * sib%diag%cum_drywt(2) *      &
                              ((sib%diag%gdd - 2650.0) / 250.0)
       
      elseif( sib%diag%gdd >= 2900.0 .or. time%doy>=sib%diag%pd+175) then

          sib%diag%leafwt_c=sib%diag%cum_drywt(2)*0.01

      endif

!itb_crop...SOY
   elseif(sib%param%biome == 21.0) then

      if (sib%diag%pd > 0                     .AND.      &
         time%doy    >= (sib%diag%pd+10)     .AND.      &
         time%doy    <  (sib%diag%pd+75))    then

	       sib%diag%leafwt_c = sib%diag%cum_drywt(2)
	
      elseif (sib%diag%pd >  0                    .AND.      &
             time%doy    >= (sib%diag%pd+75)     .AND.  &
              time%doy   <  (sib%diag%pd+100))   then

	       sib%diag%leafwt_c = sib%diag%cum_drywt(2) * 0.85 
  
            elseif (sib%diag%pd >  0                    .AND.      &
             time%doy    >= (sib%diag%pd+100)    .AND.      &
             time%doy    <  (sib%diag%pd+140))   then

               x=(time%doy - (sib%diag%pd+100)) / 40.0

	           sib%diag%leafwt_c = sib%diag%cum_drywt(2)*0.85 -   &
                 (0.85-0.1) * sib%diag%cum_drywt(2) * x 


              elseif (sib%diag%pd >  0                .AND.      &
             time%doy    >= sib%diag%pd+121) then
                 sib%diag%leafwt_c = sib%diag%cum_drywt(2)*0.01


	 else

	       sib%diag%leafwt_c = sib%diag%cum_drywt(2) * 0.01

       

       endif

  

   elseif(sib%param%biome == 22.0) then

      if( sib%diag%gdd == 0.0 ) then

           sib%diag%leafwt_c =  0.01 
      
      elseif (sib%diag%gdd > 0.0001 .and. sib%diag%gdd < 2230) then

	       sib%diag%leafwt_c = sib%diag%cum_drywt(2)
	
      elseif (sib%diag%gdd >= 2230.0 .and. sib%diag%gdd < 2450.0) then

          sib%diag%leafwt_c = 0.95 * sib%diag%cum_drywt(2) -   &
          (0.95-0.1) * sib%diag%cum_drywt(2) *      &
                              ((sib%diag%gdd - 2230.0) / 220.0)
       
      elseif( sib%diag%gdd >= 2450.0) then

          sib%diag%leafwt_c=sib%diag%cum_drywt(2)*0.01

      endif
 else

      print*,'INCORRECT CROP TYPE SELECTED'
      print*,'Biome Number=',sib%param%biome
      stop'SUBROUTINE LEAF_WEIGHT'
   endif


   end subroutine leaf_weight

