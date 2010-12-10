!================ SUBROUTINE LEAF WEIGHT ===============


!itb_crop...Short program to determine final leaf weight,
!itb_crop...based on crop, time since planting, and 
!itb_crop...cumulative dry weight

   subroutine leaf_weight(sib,time)

use kinds
use sibtype
use timetype

implicit none

real(kind=dbl_kind) :: x

!------------------------------------------------------------------
type(sib_t), intent(inout) :: sib
type(time_struct), intent(in) :: time
!------------------------------------------------------------------


!itb_crop...CORN

   if(sib%param%biome == 20.0) then

 if( sib%diag%gdd == 0.0 ) then

           sib%diag%leafwt_c =  0.01 
      
      elseif (sib%diag%gdd > 0.0001 .and. sib%diag%gdd < 2650) then

	       sib%diag%leafwt_c = sib%diag%cum_drywt(2)
	
      elseif (sib%diag%gdd >= 2650.0 .and. sib%diag%gdd < 2730.0) then

          sib%diag%leafwt_c = 0.95 * sib%diag%cum_drywt(2) -   &
          (0.95-0.5) * sib%diag%cum_drywt(2) *      &
                              ((sib%diag%gdd - 2650.0) / 80.0)

!EL...allowing for carbon remobilization from senescing leaves  to growing products
  
          sib%diag%cum_drywt(4)=sib%diag%cum_drywt(4)+((sib%diag%cum_drywt(2)- sib%diag%leafwt_c))*0.03

       elseif (sib%diag%gdd >= 2730.0 .and. sib%diag%gdd < 2900.0) then 
         
          sib%diag%leafwt_c = 0.5 * sib%diag%cum_drywt(2) -   &
          (0.5-0.01) * sib%diag%cum_drywt(2) *      &
                              ((sib%diag%gdd - 2730.0) / 170.0)


 !EL...introducing harvest towards  the end of field drying     
      elseif( sib%diag%gdd >= 2900.0) then

          sib%diag%leafwt_c=0.0001
 
      endif


!itb_crop...SOY
   elseif(sib%param%biome == 21.0) then

	 if (sib%diag%pd > 0                     .AND.      &
         time%doy    >= (sib%diag%pd+10)     .AND.      &
         time%doy    <  (sib%diag%pd+60))    then

	       sib%diag%leafwt_c = sib%diag%cum_drywt(2)
	
         elseif (sib%diag%pd >  0                    .AND.      &
             time%doy    >= (sib%diag%pd+60)     .AND.  &
              time%doy   <  (sib%diag%pd+90))   then


              x=(time%doy - (sib%diag%pd+60)) / 30.0

	           sib%diag%leafwt_c = sib%diag%cum_drywt(2)*1.0 -   &
                 (1.0-0.85) * sib%diag%cum_drywt(2) * x 
  
     
	 elseif (sib%diag%pd >  0                    .AND.      &
             time%doy    >= (sib%diag%pd+90)    .AND.      &
             time%doy    <  (sib%diag%pd+121))   then

 

               x=(time%doy - (sib%diag%pd+90)) / 31.0

	           sib%diag%leafwt_c = sib%diag%cum_drywt(2)*0.85 -   &
                 (0.85-0.01) * sib%diag%cum_drywt(2) * x 


          elseif (sib%diag%pd >  0                .AND.      &
             time%doy    >= sib%diag%pd+121) then
                 sib%diag%leafwt_c = 0.0001


	 else
                 sib%diag%leafwt_c = sib%diag%cum_drywt(2) * 0.01

        endif


   elseif(sib%param%biome == 22.0 .OR.sib%param%biome == 23.0 ) then


     if( sib%diag%gdd == 0.0 ) then

           sib%diag%leafwt_c =  0.01 

      
      elseif (sib%diag%gdd > 0.0 .and. sib%diag%gdd <= 2269.0) then
	       sib%diag%leafwt_c = sib%diag%cum_drywt(2)


      elseif (sib%diag%gdd >= 2269.0 .and. sib%diag%gdd < 2440.0) then

          sib%diag%leafwt_c = sib%diag%cum_drywt(2) -   &
          (1.0-0.01) * sib%diag%cum_drywt(2) *      &
                              ((sib%diag%gdd - 2269.0) / 171.0)    
       
 !EL...introducing harvest towards  the end of field drying     

      elseif(sib%diag%gdd >2440.0) then

          sib%diag%leafwt_c=0.0001

      endif
 else

      print*,'INCORRECT CROP TYPE SELECTED'
      print*,'Biome Number=',sib%param%biome
      stop'SUBROUTINE LEAF_WEIGHT'
   endif


   end subroutine leaf_weight

