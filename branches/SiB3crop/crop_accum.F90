
!==================SUBROUTINE CROP_ACCUM=======================================

subroutine crop_accum(sib,time,timevar)

use kinds
use sibtype
use timetype
use sib_bc_module
use physical_parameters, only: tice    
use sib_const_module, only:     latsib

implicit none

integer(kind=int_kind) :: i0,n,i,j
real(kind=dbl_kind)    :: temp_accum,assim_accum,allocwt_accum,drate,   &
                          dgrowth,max_wmain,assimd_new,x

!itb_crop
real(kind=dbl_kind)  :: coeff  ! catch-all variable for coefficients 
                               !   used multiple times



! begin time dependant, output variables
type time_dep_var
   real (kind=real_kind) :: fPAR    ! Canopy absorbed fraction of PAR
   real (kind=real_kind) :: LAI     ! Leaf-area index
   real (kind=real_kind) :: Green   ! Canopy greeness fraction of LAI
   real (kind=real_kind) :: zo      ! Canopy roughness coeff 
   real (kind=real_kind) :: zp_disp ! Zero plane displacement
   real (kind=real_kind) :: RbC     ! RB Coefficient (c1)
   real (kind=real_kind) :: RdC     ! RC Coefficient (c2)
   real (kind=real_kind) :: gmudmu  ! Time-mean leaf projection
end type time_dep_var

type(time_dep_var) TimeVar
type(aero_var),dimension(50,50) :: tempaerovar

real(kind=real_kind),dimension(2,2) :: temptran,tempref

integer(kind=int_kind) :: temp_biome

integer(kind=int_kind) :: gdd_flag = 0

!------------------------------------------------------------------
type(sib_t), intent(inout) :: sib
type(time_struct), intent(in) :: time
!----------------------------------------------------------------------  

   temp_accum = 0.0_dbl_kind
   	
!print*,'in crop_accum: leafwt_c = ',sib%diag%leafwt_c

   do i0 = 1, sib%diag%tb_indx
    
      temp_accum = temp_accum + sib%diag%tb_temp(i0) 

   enddo



   sib%diag%ta_bar = temp_accum / float(sib%diag%tb_indx)


!EL-conversion of avg. daily temperature (ta_bar) from Kelvin to Fahrenheit

sib%diag%tempf=((sib%diag%ta_bar-273.15)*1.8)+32.0	



!EL-conversion of avg. daily temperature (ta_bar) from Kelvin to Celcius

sib%diag%tempc=sib%diag%ta_bar - tice  !tice=273K


!EL-Calling for different phenology schemes based on the biome type

!itb_crop...BIOMES 20-29 reserved for crops
!itb_crop...
!itb_crop...20 - corn
!itb_crop...21 - soy

   if(sib%param%biome >= 20.0) then
        temp_biome = 12
   else
        temp_biome = int(sib%param%biome)
   endif

print*,'crop accum: BIOME=',sib%param%biome

	if(sib%param%biome == 20) then  
		call corn_phen	
	elseif(sib%param%biome == 21 ) then
		call soy_phen
	endif


!itb_crop...need to calculate veg params daily...


   if(sib%diag%phen_switch == 1 ) then


      tempaerovar = aerovar(:,:,temp_biome)

      timevar%lai = sib%diag%phen_lai

      temptran(1,1) = sib%param%tran(1,1)
      temptran(1,2) = sib%param%tran(1,2)
      temptran(2,1) = sib%param%tran(2,1)
      temptran(2,2) = sib%param%tran(2,2)

      tempref(1,1) = sib%param%ref(1,1)
      tempref(1,2) = sib%param%ref(1,2)
      tempref(2,1) = sib%param%ref(2,1)
      tempref(2,2) = sib%param%ref(2,2)

!print*,'phen_mapper:',timevar%lai,sib%diag%leafwt_c
        	
      call phen_mapper(                              &
         latsib(1),                             &
         time%mid_month(time%pmonth),           &
         sib%param%vcover,                  &
         sib%param%chil,                     &
         temptran,                              &
         tempref,                               & 
         morphtab(temp_biome),          &
         tempaerovar,                           &
         laigrid,                               &
         fvcovergrid,                           &
         timevar)		
			

	     sib%param%aparc1 = timevar%fpar
         sib%param%zlt1 = timevar%lai
         sib%param%green1 = timevar%green
         sib%param%z0d1 = timevar%zo
         sib%param%zp_disp1 = timevar%zp_disp
         sib%param%rbc1 = timevar%rbc
         sib%param%rdc1 = timevar%rdc
         sib%param%gmudmu1 = timevar%gmudmu

	     sib%param%aparc2 = timevar%fpar
         sib%param%zlt2 = timevar%lai
         sib%param%green2 = timevar%green
         sib%param%z0d2 = timevar%zo
         sib%param%zp_disp2 = timevar%zp_disp
         sib%param%rbc2 = timevar%rbc
         sib%param%rdc2 = timevar%rdc
         sib%param%gmudmu2 = timevar%gmudmu

         call bc_interp(sib, time)
			
   endif

!at the end of each day tb_index is set to zero
   sib%diag%tb_indx = 0   

!print'(a,2g12.5)','CROP ACCUM:',sib%diag%phen_LAI,sib%diag%leafwt_c

contains


















!------------------------------------------------------------------------------
subroutine corn_phen
!------------------------------------------------------------------------------




!itb...local vars

   real(kind=dbl_kind) :: temp1

!--------------------------
!Calculate the planting date
!---------------------------

!EL...sib%diag%ndf_opt= no. of days with avg. temperature above 57F

    if (sib%diag%tempf<57.0) then

       sib%diag%ndf_opt=0			!sib%diag%ndf_opt= no. of days withe
                                    !     avg. temperature above 57F
    
    elseif (sib%diag%tempf>=57.0) then
  
       sib%diag%ndf_opt=sib%diag%ndf_opt+1

   endif

   if (sib%diag%ndf_opt==7)  then

      sib%diag%pd7_est=time%doy
      sib%diag%pdindx7=sib%diag%pdindx7+(time%doy)

   endif

!EL...pdindx7 was added to avoid any second date with ndf_opt=7 
!EL...becoming a planting date..

   if (sib%diag%ndf_opt == 7 .AND.      &
       sib%diag%pdindx7 < (sib%diag%pd7_est+7) ) then

      sib%diag%pd7 = time%doy
      sib%diag%pd  = sib%diag%pd7


   endif
 


   if (sib%diag%pd7 > 0 .AND.    &
       (time%doy == (sib%diag%pd7+1)  .OR. time%doy == (sib%diag%pd7+1)    &
   .OR. time%doy == (sib%diag%pd7+3)  .OR. time%doy == (sib%diag%pd7+4)    &
   .OR. time%doy == (sib%diag%pd7+5)  .OR. time%doy == (sib%diag%pd7+6)    &
   .OR. time%doy == (sib%diag%pd7+7)) .AND. sib%diag%tempf < 53.0) then

      sib%diag%gdd = 0.0
      sib%diag%pd  = sib%diag%pd7+21

   endif
  


!----------------------------
!Calculate growing degree days
!-----------------------------

!itb_crop...gdd flag to determine initial LAI on day that
!itb_crop...seeds emerge from ground

   if(sib%diag%gdd > 100.0_dbl_kind) gdd_flag = 1


!EL...added to avoid gdd calculation before real planting date, 
!EL...since pd is printed out as 0 before the real planting 
!EL...date based on the above ndf_opt criterion

	if (sib%diag%pd > 0                  .AND.         & 
		time%doy  >=  sib%diag%pd        .AND.         &
        sib%diag%tempf > 50.0   .AND.         &
        sib%diag%tempf < 86.0)      then

  	     sib%diag%gdd=sib%diag%gdd + sib%diag%tempf- 50.0_dbl_kind
	
	endif


!itb_crop...harvest: reset
   
	if (time%doy>sib%diag%pd+175) then

             sib%diag%gdd=0.0001

	endif
      


!------------------------------
!	Reading and summing assimn
!-------------------------------
   assim_accum=0.0_dbl_kind

   do i0 = 1, sib%diag%tb_indx

!EL...multiplied by the no. secs per each timestep (i.e. tbsib)
!EL...to convert assim mol sec-1 to mol
    
      assim_accum = assim_accum + (sib%diag%tb_assim(i0)*time%dtsib) 

   enddo

   sib%diag%assim_d = assim_accum * 12.0 !multiplied by 12 to convert mol to g


!-----------------------------------------------------------
! allocation sheme for fractions for assimilate partitioning
!-----------------------------------------------------------   
!EL...1-roots, 2-leaves,3-stems,4-products (flowers and grains)

	if(sib%diag%gdd>=100.0 .and.sib%diag%gdd <500.0)then

		sib%diag%alloc(1)=0.5
		sib%diag%alloc(2)=0.25
		sib%diag%alloc(3)=0.25	
		sib%diag%alloc(4)=0.0	

    elseif(sib%diag%gdd>=500.0 .and. sib%diag%gdd <1000.0)then

        sib%diag%alloc(1)=0.5-(0.5-0.3)*(sib%diag%gdd-500)/(1000-500)
   		sib%diag%alloc(2)=0.25-(0.25-0.35)*(sib%diag%gdd-500)/(1000-500)
		sib%diag%alloc(3)=0.25-(0.25-0.35)*(sib%diag%gdd-500)/(1000-500)
		sib%diag%alloc(4)=0.0

		
	elseif(sib%diag%gdd>=1000.0 .and. sib%diag%gdd<1180.0)then

        sib%diag%alloc(1)=0.3-(0.3-0.3)*(sib%diag%gdd-1000)/(1180-1000)
		sib%diag%alloc(2)=0.35-(0.35-0.35)*(sib%diag%gdd-1000)/(1180-1000)
		sib%diag%alloc(3)=0.35-(0.35-0.35)*(sib%diag%gdd-1000)/(1180-1000)	
		sib%diag%alloc(4)=0.0

    elseif(sib%diag%gdd>=1180.0 .and. sib%diag%gdd<1360.0)then

        sib%diag%alloc(1)=0.3-(0.3-0.2)*(sib%diag%gdd-1180)/(1360-1180)
		sib%diag%alloc(2)=0.35-(0.35-0.45)*(sib%diag%gdd-1180)/(1360-1180)	
		sib%diag%alloc(3)=0.35-(0.3-0.2)*(sib%diag%gdd-1180)/(1360-1180)	
		sib%diag%alloc(4)=0.0-(0.0-0.15)*(sib%diag%gdd-1180)/(1360-1180)


	elseif(sib%diag%gdd>=1360.0 .and. sib%diag%gdd<1660.0)then

        sib%diag%alloc(1)=0.2-(0.2-0.05)*(sib%diag%gdd-1360)/(1660-1360)
		sib%diag%alloc(2)=0.45-(0.45-0.01)*(sib%diag%gdd-1360)/(1660-1360)	
		sib%diag%alloc(3)=0.2-(0.2-0.05)*(sib%diag%gdd-1360)/(1660-1360)		
		sib%diag%alloc(4)=0.15-(0.15-0.89)*(sib%diag%gdd-1360)/(1660-1360)	


	elseif(sib%diag%gdd>=1660.0 .and. sib%diag%gdd<2730.0)then

        sib%diag%alloc(1)=0.05-(0.05-0.05)*(sib%diag%gdd-1660)/(2730-1660)
		sib%diag%alloc(2)=0.01-(0.01-0.0)*(sib%diag%gdd-1660)/(2730-1660)
		sib%diag%alloc(3)=0.05-(0.05-0.0)*(sib%diag%gdd-1660)/(2730-1660)
		sib%diag%alloc(4)=0.89-(0.89-0.95)*(sib%diag%gdd-1660)/(2730-1660)	


	elseif(sib%diag%gdd < 100.0 .or. sib%diag%gdd >= 2730.0)then

		sib%diag%alloc(:)=0.0

        
        endif

!----------------------------------
!Calculate total weight allocation	
!----------------------------------
!EL...w_main is the daily amount of dry matter added, which is the basis for 
!EL...calculation of growth and maintenance respn

!EL...Multiplication factors below were derived using the info taken 
!EL...from past studies; mainly from de Vries et al., 1989.
    
    if((sib%diag%gdd >= 1300.0) .AND. (sib%diag%gdd < 2730.0)) then

       	sib%diag%w_main = sib%diag%assim_d /       &
               ((sib%diag%alloc(1) * 1.2214) +     &
                (sib%diag%alloc(2) * 1.2515) +     & 
                (sib%diag%alloc(3) * 1.2225) +     &
                (sib%diag%alloc(4) * 1.2095))

    endif

!EL...considering the fact that early corn growth is linearly related 
!EL...to temperature(Muchow and Carberry, 1989) &
!EL...considering total daily growth, based on the growth rate 
!EL...information given in de Vries et al.1989...

	if (sib%diag%gdd==100.0) then

		sib%diag%w_main=4.9 

!EL...(initial w_main g m-2 at emergence=mean of w_main from several 
!EL...iterations of the offline model using sib assimilation, 
!EL...and also based on the observed values from past studies, ),
!EL...and considering a planting density most commonly used (i.e. row 
!EL...spacing- 30", and plant spacing 6")

	endif

!EL...basic daily growth rate for the vegetative growth stages, depending 
!EL...on the average no. of days between planting and the end of V18 stage

 if (sib%diag%gdd>0.0001 .AND. sib%diag%gdd<=1300.0) then 
		drate=0.017	

!EL...g m-2; based on the range of results from original SiB 
!EL...simulations from the past
		max_wmain=4.9/0.3 

!EL...Temperatures and relevant fractions from the basic 
!EL...drates based on temperature 
!EL...were set based on the info in de Vries et al. 1989, 
!EL...and slightly modified  by looking 
!EL...at the observed LAI ranges for corn in Bondville.
  
	if (sib%diag%tempc<=8) then

		dgrowth = max_wmain*drate*0.01

	elseif (sib%diag%tempc >= 8 .and. sib%diag%tempc < 14) then

     	dgrowth = max_wmain * drate * (0.01-((0.01-0.2)      &
                  * (sib%diag%tempc - 8) / (14 - 8)))

	elseif (sib%diag%tempc >= 14 .and. sib%diag%tempc < 19) then

     	dgrowth = max_wmain * drate * (0.2 - ((0.2 - 0.6)      &
                 * (sib%diag%tempc - 14) / (19 - 14)))

	elseif (sib%diag%tempc >= 19 .and. sib%diag%tempc < 28) then

     	dgrowth = max_wmain * drate * (0.6-((0.6-1.0)      &
                  * (sib%diag%tempc - 19) / (28 - 19)))

	elseif (sib%diag%tempc >= 28 .and. sib%diag%tempc < 35) then

     	dgrowth = max_wmain * drate * (1.0 - ((1.0 - 0.9)      &
                   * (sib%diag%tempc - 28) / (35 - 28)))

	elseif (sib%diag%tempc>=35 .and. sib%diag%tempc<45) then

     	dgrowth = max_wmain * drate * (0.9-((0.9-0.01)      &
                  * (sib%diag%tempc - 35) / (45 - 35)))

	endif

   sib%diag%w_main = sib%diag%w_main+dgrowth

!EL..back calculation of the new assim_d:

	assimd_new = sib%diag%w_main *       &
                ((sib%diag%alloc(1) * 1.2214) +     &
                ( sib%diag%alloc(2) * 1.2515) +     & 
                ( sib%diag%alloc(3) * 1.2225) +     &
                ( sib%diag%alloc(4) * 1.2095))

    sib%diag%assim_d = assimd_new
 
endif

!itb_crop...harvest: reset 
if (time%doy > sib%diag%pd + 175)then

        sib%diag%w_main      = 0.0001
	    sib%diag%assim_d     = 0.0001
	    sib%diag%pd          = 0
        sib%diag%pd7         = 0        
        sib%diag%pd7_est     = 0
        sib%diag%pdindx7     = 0 
        sib%diag%phen_switch = 0

endif


!EL...Calculate w_main allocation (i.e. allocwt) to different plant parts
!EL.. calculates absolute allocation for roots(1),leaves(2),
!EL...stems(3), and products(4)
 
   do i = 1,4

        sib%diag%allocwt(i) = sib%diag%w_main * sib%diag%alloc(i) 

   enddo



!--------------------------
!Calculate maintanence respiration

!--------------------------

!EL...0.18 is the nonstructural C fraction of root C needing 
!EL...maintenance(calculations based on Brouquisse et al., 1998)
!EL...maint. coeff. info from Penning de Vries,1989, Amthor, 
!EL...1984, and Norman and Arkebauer, 1991)


       temp1 = (0.03 * 2.0 * 12.0 / 44.0) *                   &
             (2.0**((sib%diag%tempc - 20.0) / 10.0))
   
       sib%diag%phen_maintr(1) = sib%diag%cum_wt(1)       &
             * 0.18 * temp1  


!EL...(de Vries et al., 1989)
!EL.. 0.27 is the nonstructural fraction of leaf C needing 
!EL...maintenance (calculations based on Brouquisse et al., 1998)

       sib%diag%phen_maintr(2) = sib%diag%cum_wt(2) *     &
              0.27 * temp1


!EL... 0.24 is the nonstructural fraction of stem C  needing 
!EL...maintenance(calculations based on Brouquisse et al., 1998).

       temp1 = 0.01 * 2.0 * 12.0 / 44.0 * (2.0**((sib%diag%tempc-20.0) / 10.0))

       sib%diag%phen_maintr(3) = sib%diag%cum_wt(3)          &
              * 0.24 * temp1


!EL.. 0.7 is the nonstructural fraction of seed C  needing maintenance 
!EL...(calculations based on Thornton et al., 1969,Beauchemin et al., 1997)

       temp1 = 0.015 * 2.0 * 12.0 / 44.0 *      &
               (2.0 ** ((sib%diag%tempc - 20.0) / 10.0))

       sib%diag%phen_maintr(4) = sib%diag%cum_wt(4) * 0.7 * temp1




!---------------------------------------------    
!Calculate cumulative drywt.in each plant part 
!---------------------------------------------
!EL..cumulative dry weight will be used in calculation of maintenance respn

	 do j=1,4	 

             sib%diag%cum_wt(j) = sib%diag%cum_wt(j) + sib%diag%allocwt(j)
       

!itb_crop...this may have to be removed if crops are growing in
!itb_crop...Dec/Jan
	  	if ((time%doy-1)==0) then
          sib%diag%cum_wt(j)=0
	  	endif


 
     enddo
       

		
	  if (time%doy > sib%diag%pd + 175) then
      
        do i = 1,4
		
          sib%diag%cum_wt(i) = 0.0001

        enddo


	  endif

!----------------------------
!Calculate growth respiration
!----------------------------
 
        sib%diag%phen_growthr(1) = sib%diag%allocwt(1)*2*0.406 * 12/44
        sib%diag%phen_growthr(2) = sib%diag%allocwt(2)*2*0.461 * 12/44
        sib%diag%phen_growthr(3) = sib%diag%allocwt(3)*2*0.408 * 12/44
        sib%diag%phen_growthr(4) = sib%diag%allocwt(4)*2*0.384 * 12/44



!itb_crop...harvest: reset
		
	if (time%doy > sib%diag%pd+175) then
	
       do i = 1,4
         sib%diag%phen_maintr(i) = 0.0001
       enddo

	endif
		
!------------------------------
!Calculate dry weight change
!-----------------------------

      do i = 1,4
        sib%diag%wch(i) = sib%diag%allocwt(i) - sib%diag%phen_maintr(i)
      enddo
	

!itb_crop...harvest: reset

	  if (time%doy>sib%diag%pd+175) then
   
        do i = 1,4
		  sib%diag%wch(i)=0.0001
        enddo

	  endif

!--------------------------------------------------------------
!Calculate final cumulative dry weight (g C m-2) of each plant part
!--------------------------------------------------------------
	  do j=1,4	 
       
     	sib%diag%cum_drywt(j) = sib%diag%cum_drywt(j) + sib%diag%wch(j)

       
!itb_crop...may have to be removed for Dec/Jan crops
	  	if ((time%doy-1)==0) then
          sib%diag%cum_drywt(j)=0
	  	endif


 
      enddo
	  


!itb_crop...harvest: reset
	  if (time%doy > sib%diag%pd + 175) then

        do i = 1,4
          sib%diag%cum_drywt(i) = 0.0001
        enddo

	  endif


!------------------------------------------------
!final leaf weight (C g m-2) 
!-------------------------------------------------
!EL...i.e. after adjustment for senescence and harvest event
 
      if( sib%diag%gdd == 0.0 ) then

           sib%diag%leafwt_c =  0.01 
      
      elseif (sib%diag%gdd > 0.0001 .and. sib%diag%gdd < 2500) then

	       sib%diag%leafwt_c = sib%diag%cum_drywt(2)
	
      elseif (sib%diag%gdd >= 2500.0 .and. sib%diag%gdd < 2950.0) then

          sib%diag%leafwt_c = 0.95 * sib%diag%cum_drywt(2) -   &
          (0.95-0.1) * sib%diag%cum_drywt(2) *      &
                              ((sib%diag%gdd - 2500.0) / 450.0)
       
      elseif( sib%diag%gdd >= 2950.0 .or. time%doy>=sib%diag%pd+175) then

          sib%diag%leafwt_c=sib%diag%cum_drywt(2)*0.01

      endif

!--------------
!Calculate LAI
!-------------
!EL...convert to dry weight g m-2 and then multiply by SLA; 	
!EL...SLA determined by the averages based on several studies



      sib%diag%phen_LAI=sib%diag%leafwt_c*2*0.02 




!itb_crop...at the moment that growing degree days (gdd) passes
!itb_crop...100, we will initialize the LAI

    if(sib%diag%gdd >= 100.0_dbl_kind .AND. gdd_flag == 0) then

       sib%diag%phen_switch = 1

    endif


    if ((time%year==time%year + 1) .AND. sib%diag%doy==1) then
            sib%diag%assim_d=0.0001
    endif


print'(a,i8,4f12.4)','PD:',sib%diag%pd,(sib%diag%cum_wt(j),j=1,4)

end subroutine corn_phen





















!----------------------------------------------------------------
subroutine soy_phen
!----------------------------------------------------------------
!itb...local vars

   real(kind=dbl_kind) :: temp1

!--------------------------
!Calculate the planting date
!---------------------------

!EL...sib%diag%ndf_opt = no. of days with avg. temperature above 67F

	if (sib%diag%tempf < 67.0) then

    	sib%diag%ndf_opt = 0			

!EL...sib%diag%ndf_opt = no. of days with avg. temperature above 67F

	elseif (sib%diag%tempf >= 67.0) then
   	
           sib%diag%ndf_opt = sib%diag%ndf_opt + 1

	endif

	if (sib%diag%ndf_opt == 10 .AND. sib%diag%gdd < 200.0) then
 
!EL...gdd factor added to avoid taking any other pds later...

    	   sib%diag%pd = time%doy

	endif

!----------------------------
!Calculate growing degree days
!-----------------------------

!itb_crop...gdd flag to determine initial LAI on day that
!itb_crop...seeds emerge from ground

!   if(sib%diag%gdd > 100.0_dbl_kind) gdd_flag = 1


!EL...added to avoid gdd calculation before real planting date, 
!EL...since pd is printed out as 0 before the real planting 
!EL...date based on the above ndf_opt criterion

	if (sib%diag%pd    >  0                .AND.         & 
		time%doy       >= sib%diag%pd      .AND.         &
        sib%diag%tempf >  50.0             .AND.         &
        sib%diag%tempf <  86.0  )      then

      sib%diag%gdd=sib%diag%gdd + sib%diag%tempf - 50.0_dbl_kind
	
	endif

!EL...The following was added to avoid any occurrence of gdd at the 
!EL...beginning of the year, before planting

        if (time%doy < sib%diag%pd .AND. sib%diag%tempf < 50) then

            sib%diag%gdd = 0.0
            sib%diag%pd  = 0

        endif


!EL...to avoid gdd calculation after harvesting is done, allowing 
!EL...ample time between planting and harvesting

	if (time%doy > sib%diag%pd + 175) then

             sib%diag%gdd = 0.0001

	endif


!------------------------------
!	Reading and summing assimn
!-------------------------------

    assim_accum=0.0_dbl_kind
 

!EL...multiplied by the no. secs per each timestep (i.e. dtsib) 
!EL...to convert assim mol sec-1 to mol

	do i0 = 1, sib%diag%tb_indx
    
      assim_accum = assim_accum + (sib%diag%tb_assim(i0) * time%dtsib) 

    enddo


    sib%diag%assim_d = assim_accum * 12 !multiplied by 12 to convert mol to g


!-----------------------------------------------------------
! allocation sheme for fractions for assimilate partitioning
!-----------------------------------------------------------   
!EL...1-roots, 2-leaves,3-stems,4-products (flowers and grains)

	 if (sib%diag%pd>0 .AND. time%doy==(sib%diag%pd+10)) then

            sib%diag%alloc(1) = 0.5
			sib%diag%alloc(2) = 0.25
			sib%diag%alloc(3) = 0.25	
			sib%diag%alloc(4) = 0.0	
        
        else if (sib%diag%pd >  0                    .AND.    &
                 time%doy    >= (sib%diag%pd + 10)   .AND.    &
                 time%doy    <  (sib%diag%pd + 30) ) then

            sib%diag%alloc(1) = 0.5  - (0.5  - 0.5) *      &
                                (time%doy - (sib%diag%pd + 10)) / 20

            sib%diag%alloc(2) = 0.25 - (0.25 - 0.25) *     &
                                (time%doy - (sib%diag%pd + 10)) / 20 

            sib%diag%alloc(3) = 0.25 - (0.25 - 0.25) *     &
                                (time%doy - (sib%diag%pd + 10)) / 20

            sib%diag%alloc(4) = 0.0
         
         
         else if (sib%diag%pd > 0                       .AND.     &
                  time%doy    >= (sib%diag%pd + 30)     .AND.     &
                  time%doy    <  (sib%diag%pd + 40) )   then

            sib%diag%alloc(1) = 0.5  - (0.5  - 0.49) *      &
                                (time%doy - (sib%diag%pd + 30)) / 10

			sib%diag%alloc(2) = 0.25 - (0.25 - 0.26) *      &
                                (time%doy - (sib%diag%pd + 30)) / 10

			sib%diag%alloc(3) = 0.25 - (0.25 - 0.25) *      &
                                (time%doy - (sib%diag%pd + 30)) / 10
	
			sib%diag%alloc(4) = 0.0   

        else if (sib%diag%pd >  0                     .AND.     & 
                 time%doy    >= (sib%diag%pd + 40)    .AND.     &
                 time%doy    <  (sib%diag%pd + 50) )  then

            sib%diag%alloc(1) = 0.49 - (0.49 - 0.35) *     &
                                (time%doy - (sib%diag%pd + 40)) / 10

			sib%diag%alloc(2) = 0.26 - (0.26 - 0.40) *     &
                                (time%doy - (sib%diag%pd + 40)) / 10

			sib%diag%alloc(3) = 0.25 - (0.25 - 0.25) *     &
                                (time%doy - (sib%diag%pd + 40)) / 10
	
			sib%diag%alloc(4) = 0.0   

      
        else if (sib%diag%pd >  0                   .AND.     &
                 time%doy    >= (sib%diag%pd + 50)  .AND.     &
                 time%doy    <  (sib%diag%pd + 60)) then

            sib%diag%alloc(1) = 0.35 - (0.35 - 0.22) *      &
                                (time%doy - (sib%diag%pd + 50)) / 10

			sib%diag%alloc(2) = 0.40 - (0.40 - 0.53) *      &
                                (time%doy - (sib%diag%pd + 50)) / 10

			sib%diag%alloc(3) = 0.25 - (0.25 - 0.25) *      &
                                (time%doy - (sib%diag%pd + 50)) / 10
	
			sib%diag%alloc(4) = 0.0    
  

        else if (sib%diag%pd >  0                    .AND.    &
                 time%doy    >= (sib%diag%pd + 60)   .AND.    &
                 time%doy    <  (sib%diag%pd + 75))  then

            sib%diag%alloc(1) = 0.22 - (0.22 - 0.21) *      &
                                (time%doy - (sib%diag%pd + 60)) / 15

			sib%diag%alloc(2) = 0.53 - (0.53 - 0.50) *      &
                                (time%doy - (sib%diag%pd + 60)) / 15

			sib%diag%alloc(3) = 0.25 - (0.25 - 0.29) *      &
                                (time%doy - (sib%diag%pd + 60)) / 15	

			sib%diag%alloc(4) = 0.0


        else if (sib%diag%pd >  0                    .AND.     &
                 time%doy    >= (sib%diag%pd + 75)   .AND.     &
                 time%doy    <  (sib%diag%pd + 80))  then

            sib%diag%alloc(1) = 0.21 - (0.21 - 0.20) *      &
                                (time%doy - (sib%diag%pd + 75)) / 5

			sib%diag%alloc(2) = 0.50 - (0.50 - 0.35) *      &
                                (time%doy - (sib%diag%pd + 75)) / 5

			sib%diag%alloc(3) = 0.29 - (0.29 - 0.23) *      &
                                (time%doy - (sib%diag%pd + 75)) / 5	

			sib%diag%alloc(4) = 0.0 - (0.0-0.22) *      &
                                (time%doy - (sib%diag%pd + 75)) / 5 

        else if (sib%diag%pd >  0                      .AND.      &
                 time%doy    >= (sib%diag%pd + 80)     .AND.     &
                 time%doy    <  (sib%diag%pd + 89))    then

            sib%diag%alloc(1) = 0.20 

			sib%diag%alloc(2) = 0.35 - (0.35 - 0.20) *      &
                                (time%doy - (sib%diag%pd + 80)) / 9 

			sib%diag%alloc(3) = 0.23 - (0.23 - 0.20) *      &
                                (time%doy - (sib%diag%pd + 80)) / 9	

			sib%diag%alloc(4) = 0.22 - (0.22 - 0.40) *      &
                                (time%doy - (sib%diag%pd + 80)) / 9 

        else if (sib%diag%pd >  0                .AND.     &
                    time%doy >= (sib%diag%pd+89) .AND.     &
                    time%doy <  (sib%diag%pd+98)) then

            sib%diag%alloc(1) = 0.20 

			sib%diag%alloc(2) = 0.2 - 0.2 *      &
                                (time%doy - (sib%diag%pd + 89)) / 9

			sib%diag%alloc(3) = 0.2 - (0.2 - 0.015) *      &
                                (time%doy - (sib%diag%pd + 89)) / 9	

			sib%diag%alloc(4) = 0.4 - (0.4 - 0.785) *      &
                                (time%doy - (sib%diag%pd + 89)) / 9

        else if (sib%diag%pd >  0                .AND.     &
                 time%doy    >= (sib%diag%pd+98) .AND.     &
                 time%doy    <  (sib%diag%pd+108)) then

            sib%diag%alloc(1) = 0.2 - (0.2 - 0.1) *      &
                                (time%doy - (sib%diag%pd + 98)) / 10

			sib%diag%alloc(2) = 0.0 - (0.0 - 0.0) *      &
                                (time%doy - (sib%diag%pd + 98)) / 10

			sib%diag%alloc(3) = 0.015 - (0.015 - 0.0) *      &
                                (time%doy - (sib%diag%pd + 98)) / 10	

			sib%diag%alloc(4) = 0.785 - (0.755 - 0.93) *      &
                                (time%doy - (sib%diag%pd + 98)) / 10 

        else if (sib%diag%pd >  0                   .AND.     &
                 time%doy    >= (sib%diag%pd+108)   .AND.     &
                 time%doy    < (sib%diag%pd+130)) then

            sib%diag%alloc(1) = 0.1 
			sib%diag%alloc(2) = 0.0 
			sib%diag%alloc(3) = 0.0 
			sib%diag%alloc(4) = 0.9 


		elseif ( sib%diag%pd  > 0                   .AND.    &
                ((time%doy    < sib%diag%pd+10)     .OR.     &
                 (time%doy    > sib%diag%pd+130)))  then

			sib%diag%alloc(:)=0.0

        endif
	

!----------------------------------
!Calculate total weight allocation	
!----------------------------------
!EL...w_main is the daily amount of dry matter added, which is the basis 
!EL...for calculation of growth and maintenance respn

!EL...Multiplication factors below were derived using the info taken 
!EL...from past studies; mainly from de Vries et al., 1989.

    if (sib%diag%pd  >  0                       .AND.     &
        (time%doy    <  (sib%diag%pd + 10)      .OR.     &
         time%doy    >= (sib%diag%pd + 175)))   then 

        sib%diag%w_main =0.0001

    endif 


    if(sib%diag%pd >  0                       .AND.     &
       time%doy    >  (sib%diag%pd + 60)      .AND.     &
       time%doy    <= (sib%diag%pd + 130))    then

       	sib%diag%w_main = sib%diag%assim_d /      &
                             ((sib%diag%alloc(1 )* 1.2929) +    &
                              (sib%diag%alloc(2) * 1.431)  +    &
                              (sib%diag%alloc(3) * 1.2946) +    &
                              (sib%diag%alloc(4) *  1.6752))

    endif

  

!EL...considering the fact that early soybean growth is linearly 
!EL...related to temperature(Muchow and Carberry, 1989) &
!EL...considering total daily growth, based on the growth rate 
!EL...information given in de Vries et al.1989...

	if (sib%diag%pd > 0 .AND.time%doy == (sib%diag%pd + 10)) then

		sib%diag%w_main=2.8 

!EL...(initial w_main g m-2 at emergence=mean of w_main from several 
!EL...runs from the offline model using sib assimilation, 
!EL...and also based on the observed values from past studies, )

	endif

!EL...basic daily growth rate for the vegetative growth stages, 
!EL...depending on the average no. of days
!EL.. between planting and observed max LAI

 if (sib%diag%pd >  0                        .AND.      &
     time%doy    >  (sib%diag%pd + 10)       .AND.      &
     time%doy    <= (sib%diag%pd + 60))      then 

		drate=0.02	

!EL...g m-2; based on the range of results from original 
!EL...sib simulations from the past
		max_wmain=2.8/0.3 


!EL...Temperatures and relevant fractions from the basic 
!EL...drates based on temperature 
!EL...were set based on the info in de Vries et al. 1989, 
!EL...and slightly modified  by looking 
!EL...at the observed LAI ranges for soybean in Bondville.

    dgrowth = max_wmain * drate

	if (sib%diag%tempc<=8) then

		dgrowth = dgrowth * 0.01
	
	elseif (sib%diag%tempc >= 8.0 .AND. sib%diag%tempc<10.0) then

    	dgrowth = dgrowth * 0.01 - ((0.01 - 0.25)      &
                  * (sib%diag%tempc - 8.0) / (10.0 - 8.0))

	elseif (sib%diag%tempc >= 10.0 .AND. sib%diag%tempc < 20.0) then

    	dgrowth = dgrowth * (0.25 - ((0.25 - 0.7)      &
                  * (sib%diag%tempc - 10.0) / (20.0 - 10.0)))

	elseif (sib%diag%tempc >= 20.0 .AND. sib%diag%tempc < 27.0) then

    	dgrowth = dgrowth * (0.9 - ((0.9 - 1.0)        &
                  * (sib%diag%tempc - 20.0) / (27.0 - 20.0)))

	elseif (sib%diag%tempc >= 27.0 .and. sib%diag%tempc < 30.0) then

    	dgrowth = dgrowth * (1.0 - ((1.0 - 1.1)        &
                  * (sib%diag%tempc - 27.0) / (30.0 - 27.0)))

	elseif (sib%diag%tempc >= 30.0 .and. sib%diag%tempc < 40.0) then

    	dgrowth = dgrowth * (1.1 - ((1.1 -1.2 )        &
                  * (sib%diag%tempc - 30.0) / (40.0 - 30.0)))

	endif



   sib%diag%w_main=sib%diag%w_main+dgrowth

!EL..back calculation of the new assim_d:

	assimd_new= sib%diag%w_main *       &
               ((sib%diag%alloc(1) * 1.2929) +     &
                (sib%diag%alloc(2) * 1.431) +     & 
                (sib%diag%alloc(3) * 1.2946) +     &
                (sib%diag%alloc(4) * 1.6752))
    sib%diag%assim_d=assimd_new
 
 endif

!Calculate w_main allocation to different plant parts

!EL...calculates absolute allocation of biomass for roots(1),
!EL...leaves(2),stems(3),and products(4) using allocation fractions

        do i = 1,4
          sib%diag%allocwt(i) = sib%diag%w_main * sib%diag%alloc(i) 
        enddo




!--------------------------
!Calculate maintenance resp
!--------------------------
!EL...Q10 coefficient is 2.0 for soybean (Norman and Arkebauer, 1991);
!EL...0.32 is the nonstructural C fraction of root C needing maintenance
!EL...(calculations based on Allen et al., 1998 and Rogers et al., 2006)&
!EL...maint. coeff. info from Penning de Vries,1989, Amthor, 1984, 
!EL...and Norman and Arkebauer, 1991)
 coeff = 2.0_dbl_kind * 12.0_dbl_kind/ 44.0_dbl_kind
   
       temp1 = (0.03 * coeff) *                   &
             (2.0**((sib%diag%tempc - 20.0) / 10.0))
   
       sib%diag%phen_maintr(1) = sib%diag%cum_wt(1)       &
             * 0.32 * temp1  


!EL...(de Vries et al., 1989)
!EL.. 0.38 is the nonstructural fraction of leaf C needing maintenance
!EL...(calculations based on Allen et al., 1998 and Rogers et al., 2006).

       sib%diag%phen_maintr(2) = sib%diag%cum_wt(2) *     &
              0.38 * temp1


!EL... 0.32 is the nonstructural fraction of stem C  needing maintenance
!EL...(calculations based on Allen et al., 1998 and Rogers et al., 2006).

       temp1 = 0.01 * coeff * (2.0**((sib%diag%tempc-20.0) / 10.0))

       sib%diag%phen_maintr(3) = sib%diag%cum_wt(3)          &
              * 0.32 * temp1


!EL.. 0.46 is the nonstructural fraction of seed C  needing maintenance 
!EL...(calculations based on Allen et al., 1998 and Rogers et al., 2006).

       temp1 = 0.015 * coeff * (2.0**((sib%diag%tempc-20.0)/10.0))

       sib%diag%phen_maintr(4)= sib%diag%cum_wt(4) * 0.46*temp1


	
	  if (time%doy>sib%diag%pd+175) then
		sib%diag%phen_maintr(:)=0.0001
	  endif




!---------------------------------------------    
!Calculate cumulative drywt.in each plant part 
!---------------------------------------------
!EL...calculates cumulative dry weight before maint respn

  do j=1,4	 

     	sib%diag%cum_wt(j) = sib%diag%cum_wt(j) + sib%diag%allocwt(j)
       

!itb_crop...this may have to be removed for crops that are growing
!itb_crop...during December/January
	  if ((time%doy-1)==0) then
        sib%diag%cum_wt(j)=0
	  endif
 


  enddo
       

	 
	  if (time%doy > sib%diag%pd + 175) then
		sib%diag%cum_wt(:)=0.0001
	  endif


!----------------------------
!Calculate growth respiration
!----------------------------
 
        coeff = 2.0_dbl_kind * 12.0_dbl_kind/ 44.0_dbl_kind

        sib%diag%phen_growthr(1) = sib%diag%allocwt(1) * 0.537 * coeff
        sib%diag%phen_growthr(2) = sib%diag%allocwt(2) * 0.790 * coeff
        sib%diag%phen_growthr(3) = sib%diag%allocwt(3) * 0.540 * coeff
        sib%diag%phen_growthr(4) = sib%diag%allocwt(4) * 1.238 * coeff



!------------------------------
!Calculate dry weight change
!-----------------------------

      do i = 1,4
        sib%diag%wch(i) = sib%diag%allocwt(i) - sib%diag%phen_maintr(i)
      enddo

	 if (time%doy > sib%diag%pd + 175) then
		sib%diag%wch(:)=0.0001
	 endif

!--------------------------------------------------------------
!Recalculate final cumulative dry weight (g C m-2) of each plant part
!--------------------------------------------------------------

	  do j=1,4
	 
        sib%diag%cum_drywt(j) = sib%diag%cum_drywt(j) + sib%diag%wch(j)
       

!itb_crop...may have to be removed for Dec/Jan crops...
	  	if ((time%doy-1)==0) then
         sib%diag%cum_drywt(j)=0
	  	endif


 
          enddo


!itb_crop...harvest: reset
 	if (time%doy > sib%diag%pd + 175) then

		sib%diag%cum_drywt(:)=0.0001
	 
	endif

!------------------------------------------------
!final leaf weight (C g m-2) 
!-------------------------------------------------
!EL..i.e. after adjustment for senescence and harvest event
      
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
             time%doy    <  (sib%diag%pd+121))   then

               x=(time%doy - (sib%diag%pd+100)) / 21.0

	           sib%diag%leafwt_c = 0.85 * sib%diag%cum_drywt(2) -   &
                 (0.85-0.001) * sib%diag%cum_drywt(2) * x 

	 else

	       sib%diag%leafwt_c = sib%diag%cum_drywt(2) * 0.001

     endif

!--------------
!Calculate LAI
!-------------

!EL...convert to dry weight g m-2 and then multiply by SLA; 
!EL..SLA determined by the averages based on several studies

      sib%diag%phen_LAI = sib%diag%leafwt_c * 2.0 * 0.02 	


 	  sib%diag%tb_indx = 0	 !at the end of each day tb_index is set to zero


print*,sib%diag%allocwt(2),sib%diag%pd,sib%diag%w_main,sib%diag%phen_maintr(2),sib%diag%phen_growthr(2),sib%diag%cum_drywt(2),sib%diag%leafwt_c,sib%diag%phen_LAI

!itb_crop...at the moment that growing degree days (gdd) passes
!itb_crop...100, we will initialize the LAI

    if(sib%diag%pd >  0                     .AND.    &
       time%doy    >  (sib%diag%pd + 10)    .AND.    &
       gdd_flag    == 0)                    then

       sib%diag%phen_switch = 1

    endif

    if (time%doy < sib%diag%pd .and. time%doy>sib%diag%pd + 175)then

        sib%diag%w_main      = 0.0001
		sib%diag%assim_d     = 0.0001
	    sib%diag%pd          = 0
        sib%diag%phen_switch = 1

    endif


	
    if ((time%year   == time%year + 1) .AND.      &
         sib%diag%doy== 1) then

       sib%diag%assim_d=0.0001

    endif
	 


end subroutine soy_phen


end subroutine crop_accum
