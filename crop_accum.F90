
!==================SUBROUTINE CROP_ACCUM=======================================
subroutine crop_accum(sib,time,timevar)

use kinds
use sibtype
use timetype
use sib_bc_module
use physical_parameters, only: tice    
use sib_const_module, only:     latsib

implicit none

integer(kind=int_kind) :: i0,n,pd,i,j,ndf60,ndf65
real(kind=dbl_kind)    :: temp_accum,assim_accum,allocwt_accum,drate,max_assimd

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

integer(kind=int_kind) :: gdd_flag = 0


!----------------------------------------------------------------------
type(sib_t), intent(inout) :: sib
type(time_struct), intent(in) :: time
!----------------------------------------------------------------------  

   temp_accum = 0.0_dbl_kind
   	

   do i0 = 1, sib%diag%tb_indx
    
      temp_accum = temp_accum + sib%diag%tb_temp(i0) 

   enddo

   sib%diag%ta_bar = temp_accum / float(sib%diag%tb_indx)


!EL-conversion of avg. daily temperature (ta_bar) from Kelvin to Fahrenheit

sib%diag%tempf=((sib%diag%ta_bar-273.15)*1.8)+32.0	


!EL-conversion of avg. daily temperature (ta_bar) from Kelvin to Celcius

sib%diag%tempc=sib%diag%ta_bar - tice  !tice=273K


!EL-Calling for different phenology schemes based on the year and the crop

	if(mod(time%year,2)==0) then  
		call soy_phen
	else
		call corn_phen	
	endif


!itb_crop...need to calculate veg params daily...


   if(sib%diag%phen_switch > 0 ) then


      tempaerovar = aerovar(:,:,int(sib%param%biome))

      timevar%lai = sib%diag%phen_lai

      temptran(1,1) = sib%param%tran(1,1)
      temptran(1,2) = sib%param%tran(1,2)
      temptran(2,1) = sib%param%tran(2,1)
      temptran(2,2) = sib%param%tran(2,2)

      tempref(1,1) = sib%param%ref(1,1)
      tempref(1,2) = sib%param%ref(1,2)
      tempref(2,1) = sib%param%ref(2,1)
      tempref(2,2) = sib%param%ref(2,2)
        	
      call phen_mapper(                              &
         latsib(1),                             &
         time%mid_month(time%pmonth),           &
         sib%param%vcover,                  &
         sib%param%chil,                     &
         temptran,                              &
         tempref,                               & 
         morphtab(int(sib%param%biome)),          &
         tempaerovar,                           &
         laigrid,                               &
         fvcovergrid,                           &
         timevar)

!print*,timevar%lai		
			

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

contains

!------------------------------------------------------------------------------
subroutine corn_phen
!------------------------------------------------------------------------------


!itb...local vars

   real(kind=dbl_kind) :: temp1

!--------------------------
!Calculate the planting date
!---------------------------

   if (sib%diag%tempf<60.0) then

    ndf60=0			!ndf60= no. of days withe avg. temperature above 60F

   elseif (sib%diag%tempf>=60.0) then
  
   	ndf60=ndf60+1

   endif

   if (ndf60==5) then

     pd=time%doy

   endif

!----------------------------
!Calculate growing degree days
!-----------------------------

!itb_crop...gdd flag to determine initial LAI on day that
!itb_crop...seeds emerge from ground

   if(sib%diag%gdd > 100.0_dbl_kind) gdd_flag = 1


!EL...added to avoid gdd calculation before real planting date, 
!EL...since pd is printed out as 0 before the real planting 
!EL...date based on the above ndf60==5 criterion

	if (pd>0                  .AND.         & 
		time%doy >= pd        .AND.         &
        sib%diag%tempf>50.0   .AND.         &
        sib%diag%tempf<86.0)      then

    	sib%diag%gdd=sib%diag%gdd + sib%diag%tempf- 50.0_dbl_kind
	
	endif

	if (time%doy>300) then
        sib%diag%gdd=0.0001
	endif




!------------------------------
!	Reading and summing assimn
!-------------------------------
   assim_accum=0.0_dbl_kind

if (sib%diag%gdd<=1360.0) then 
	drate=0.025		!basic rate for the vegetative growth stages
!elseif (sib%diag%gdd>1360.0) then 
!drate=0.019		!basic rate for the reproductive growth stages
endif

if (sib%diag%gdd==100) then
	sib%diag%assim_d=4.13 !(initial assim_d at emergence=mean of assim_d from several runs from the offline model using sib assimilation)
endif
 
max_assimd=4.13/0.3 !g m-2; based on the results from offline run using sib asimilaiton and a modification of info from De Vries et al., 1989	

!if (sib%diag%gdd>=1560.0 .and.sib%diag%gdd<1660.0)then
!sib%diag%assim_d=max_assimd
!endif
  
if (sib%diag%tempc<=8) then
	assim_accum=max_assimd*drate*0.01

elseif (sib%diag%tempc>=8 .and. sib%diag%tempc<14) then
     assim_accum=max_assimd*drate*(0.01-((0.01-0.2)*(sib%diag%tempc-8)/(14-8)))

elseif (sib%diag%tempc>=14 .and. sib%diag%tempc<19) then
     assim_accum=max_assimd*drate*(0.2-((0.2-0.6)*(sib%diag%tempc-14)/(19-14)))

elseif (sib%diag%tempc>=19 .and. sib%diag%tempc<28) then
     assim_accum=max_assimd*drate*(0.6-((0.6-1.0)*(sib%diag%tempc-19)/(28-19)))

elseif (sib%diag%tempc>=28 .and. sib%diag%tempc<35) then
     assim_accum=max_assimd*drate*(1.0-((1.0-0.9)*(sib%diag%tempc-28)/(35-28)))

elseif (sib%diag%tempc>=35 .and. sib%diag%tempc<45) then
     assim_accum=max_assimd*drate*(0.9-((0.9-0.01)*(sib%diag%tempc-35)/(45-35)))

endif
! The above temperatures and relevant fractions from the basic derates based on temperature were set based on the info in de Vries et al. 1989, and slightly modified  by looking at the observed LAI ranges for corn in Bondville.

   sib%diag%assim_d=sib%diag%assim_d+assim_accum

 if (sib%diag%gdd>1360) then

   do i0 = 1, sib%diag%tb_indx

!EL...multiplied by the no. secs per each timestep (i.e. tbsib)
!EL...to convert assim mol sec-1 to mol
    
      assim_accum = assim_accum + (sib%diag%tb_assim(i0)*time%dtsib) 

   enddo

   sib%diag%assim_d = assim_accum * 12.0 !multiplied by 12 to convert mol to g
endif

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
		sib%diag%alloc(1)=0.0
		sib%diag%alloc(2)=0.0	
		sib%diag%alloc(3)=0.0	
		sib%diag%alloc(4)=0.0
        
    endif
!----------------------------------
!Calculate total weight allocation	
!----------------------------------
!EL...w_main is the (drywt+maintenancerespiration) together
!EL..the coefficients corresponds to growth resp & 
!EL...(drywt+maint R) ref Penning deVries (1989)
!EL.. here w_main is the (drywt+maint.respn)
    
    if((sib%diag%gdd >= 100.0) .AND. (sib%diag%gdd < 2730.0)) then

       	sib%diag%w_main = sib%diag%assim_d /       &
               ((sib%diag%alloc(1) * 1.2214) +     &
                (sib%diag%alloc(2) * 1.2515) +     & 
                (sib%diag%alloc(3) * 1.2225) +     &
                (sib%diag%alloc(4) * 1.2095))

    else
        sib%diag%w_main =0.0

    endif
 

!Calculate w_main allocation to different plant parts

!EL.. clculates absolute allocation for roots(1),leaves(2),stems(3), and products(4)
!EL..using allocation fraction for roots and assimilation  

         sib%diag%allocwt(1)=sib%diag%w_main*sib%diag%alloc(1) 

        sib%diag%allocwt(2)=sib%diag%w_main*sib%diag%alloc(2) 

        sib%diag%allocwt(3)=sib%diag%w_main*sib%diag%alloc(3)

        sib%diag%allocwt(4)=sib%diag%w_main*sib%diag%alloc(4)



!---------------------------------------------    
!Calculate cumulative drywt.in each plant part 
!---------------------------------------------
!EL..cumulative dry weight before respn

	  do j=1,4	 
       
     	sib%diag%cum_wt(time%doy,j)=sib%diag%cum_wt(time%doy-1,j)+sib%diag%allocwt(j)
       
	  	If ((time%doy-1)==0) then
     
        	sib%diag%cum_wt(time%doy,j)=0
     
	  	endif
 
      enddo
       
      sib%diag%cum_w(1)=sib%diag%cum_wt(time%doy,1)
	  sib%diag%cum_w(2)=sib%diag%cum_wt(time%doy,2)
	  sib%diag%cum_w(3)=sib%diag%cum_wt(time%doy,3)
	  sib%diag%cum_w(4)=sib%diag%cum_wt(time%doy,4)
		
	  if (time%doy>300) then

		sib%diag%cum_wt(time%doy,1)=0.0001
		sib%diag%cum_wt(time%doy,2)=0.0001
		sib%diag%cum_wt(time%doy,3)=0.0001
		sib%diag%cum_wt(time%doy,4)=0.0001

	  endif

!----------------------------
!Calculate growth respiration
!----------------------------
 
        sib%diag%phen_growthr(1)=sib%diag%allocwt(1)*2*0.406*12/44
        sib%diag%phen_growthr(2)=sib%diag%allocwt(2)*2*0.461*12/44
        sib%diag%phen_growthr(3)=sib%diag%allocwt(3)*2*0.408*12/44
        sib%diag%phen_growthr(4)=sib%diag%allocwt(4)*2*0.384*12/44

!--------------------------
!Calculate maintanence resp
!--------------------------

!EL...0.18 is the nonstructural C fraction of root C needing 
!EL...maintenance(calculations based on Brouquisse et al., 1998)
!EL...maint. coeff. info from Penning de Vries,1989, Amthor, 
!EL...1984, and Norman and Arkebauer, 1991)


       temp1 = (0.03 * 2.0 * 12.0 / 44.0) *                   &
             (2.0**((sib%diag%tempc - 20.0) / 10.0))
   
       sib%diag%phen_maintr(1) = sib%diag%cum_wt(time%doy-1,1)       &
             * 0.18 * temp1  


!EL...multiplied by 0.75 to incorporate that maintenance 
!EL...respiration during the daytime is half that of  nighttime.
!EL...(de Vries et al., 1989)
!EL.. 0.27 is the nonstructural fraction of leaf C needing 
!EL...maintenance (calculations based on Brouquisse et al., 1998)

       sib%diag%phen_maintr(2) = sib%diag%cum_wt(time%doy-1,2) *     &
              0.27 * temp1


!EL... 0.24 is the nonstructural fraction of stem C  needing 
!EL...maintenance(calculations based on Brouquisse et al., 1998).

       temp1 = 0.01 * 2.0 * 12.0 / 44.0 * (2.0**((sib%diag%tempc-20.0) / 10.0))

       sib%diag%phen_maintr(3) = sib%diag%cum_wt(time%doy-1,3)          &
              * 0.24 * temp1


!EL.. 0.4 is the nonstructural fraction of seed C  needing maintenance (calculations based on Beauchemin et al., 1997)

       temp1 = 0.015 * 2.0 * 12.0 / 44.0 * (2.0**((sib%diag%tempc-20.0)/10.0))

       sib%diag%phen_maintr(4)=sib%diag%cum_wt(time%doy-1,4) * temp1

 

		
	if (time%doy>300) then	
		sib%diag%phen_maintr(1)=0.0001
		sib%diag%phen_maintr(2)=0.0001
		sib%diag%phen_maintr(3)=0.0001
		sib%diag%phen_maintr(4)=0.0001
	endif
		
!------------------------------
!Calculate dry weight change
!-----------------------------

	  sib%diag%wch(1)=sib%diag%allocwt(1)- sib%diag%phen_maintr(1)
      sib%diag%wch(2)=sib%diag%allocwt(2)- sib%diag%phen_maintr(2)
      sib%diag%wch(3)=sib%diag%allocwt(3)- sib%diag%phen_maintr(3)
      sib%diag%wch(4)=sib%diag%allocwt(4)- sib%diag%phen_maintr(4)
	
	  if (time%doy>300) then
		sib%diag%wch(1)=0.0001
		sib%diag%wch(2)=0.0001
		sib%diag%wch(3)=0.0001
		sib%diag%wch(4)=0.0001
	  endif

!--------------------------------------------------------------
!Recalculate final cumulative dry weight (g C m-2) of each plant part
!--------------------------------------------------------------
	  do j=1,4	 
       
     	sib%diag%cum_drywt(time%doy,j) =            &
            sib%diag%cum_drywt(time%doy-1,j) + sib%diag%wch(j)
       
	  	If ((time%doy-1)==0) then
     
        	sib%diag%cum_drywt(time%doy,j)=0
     
	  	endif
 
      enddo
	  
 !Renamed to be output on a daily basis in a separate text file
      sib%diag%final_drywt(1)=sib%diag%cum_drywt(time%doy,1)
	  sib%diag%final_drywt(2)=sib%diag%cum_drywt(time%doy,2)
	  sib%diag%final_drywt(3)=sib%diag%cum_drywt(time%doy,3)
	  sib%diag%final_drywt(4)=sib%diag%cum_drywt(time%doy,4)

	  if (time%doy>300) then
		sib%diag%cum_drywt(time%doy,1)=0.0001
		sib%diag%cum_drywt(time%doy,2)=0.0001
		sib%diag%cum_drywt(time%doy,3)=0.0001
		sib%diag%cum_drywt(time%doy,4)=0.0001
	  endif
!------------------------------------------------
!final leaf weight (C g m-2) 
!-------------------------------------------------
!EL...after adjustment for senescence and harvest event
 
      if( sib%diag%gdd == 0.0 ) then

           sib%diag%leafwt_c =  0.01 
      
	  elseif (sib%diag%gdd > 0.0 .and. sib%diag%gdd < 2730) then

	       sib%diag%leafwt_c = sib%diag%cum_drywt(time%doy,2)
	
      elseif (sib%diag%gdd >= 2730.0 .and. sib%diag%gdd < 3300.0) then

          sib%diag%leafwt_c = 0.95 * sib%diag%cum_drywt(time%doy,2) -   &
          (0.95-0.2) * sib%diag%cum_drywt(time%doy,2) *      &
                              ((sib%diag%gdd - 2730.0) / 570.0)
       
	  elseif( sib%diag%gdd >= 3300.0) then

          sib%diag%leafwt_c=sib%diag%cum_drywt(time%doy,2)*0.01

       endif

!--------------
!Calculate LAI
!-------------
!EL...convert to dry weight g m-2 and then multiplied by SLA; 	
!EL...SLA determined by the averages based on several studies)



      sib%diag%phen_LAI=sib%diag%leafwt_c*2*0.02 

!print'(6g12.5)',sib%diag%phen_LAI,sib%diag%leafwt_c,sib%diag%cum_drywt(time%doy,2),   &
!                     sib%diag%cum_wt(time%doy,2),sib%diag%assim_d,   &
!                     sib%diag%alloc(2)
print*,sib%diag%assim_d,sib%diag%phen_LAI,timevar%lai


 	sib%diag%tb_indx = 0	 !at the end of each day tb_index is set to zero


!itb_crop...at the moment that growind degree days (gdd) passes
!itb_crop...100, we will initialize the LAI
    if(sib%diag%gdd >= 100.0_dbl_kind .AND. gdd_flag == 0) then

!       sib%param%zlt1    = sib%diag%zlt_crop_init
!       sib%param%zlt2    = sib%diag%zlt_crop_init
!       sib%param%zlt     = sib%diag%zlt_crop_init
!       sib%diag%phen_LAI = sib%diag%zlt_crop_init
!	   sib%diag%cum_drywt(2)=5.0

       sib%diag%phen_switch = 1

    endif

 
		write(20,'(i4.4,2x,i3.3,2x,46(1x,f11.2))')time%year,   &
            time%doy,sib%diag%tempf,sib%diag%tempc,            &
            sib%diag%gdd,sib%diag%assim_d,sib%diag%alloc(1:4) ,&
            sib%diag%w_main,sib%diag%allocwt(1:4),             &
            sib%diag%cum_w(1:4),sib%diag%phen_growthr(1:4),    &
            sib%diag%phen_maintr(1:4),sib%diag%wch(1:4),       &
            sib%diag%final_drywt(1:4),sib%diag%leafwt_c,       &
            sib%diag%phen_LAI,timevar%lai

end subroutine corn_phen




!-----------------------------------------------------------------------------------------------------------
subroutine soy_phen
!-----------------------------------------------------------------------------------------------------------

!--------------------------
!Calculate the planting date
!---------------------------

	if (sib%diag%tempf<65.0) then

    	ndf65=0			

!EL...ndf65= no. of days withe avg. temperature above 65F

	elseif (sib%diag%tempf>=65.0) then
   	
		 ndf65=ndf65+1

	endif

	if (ndf65==5) then

    	pd=time%doy

	endif

!----------------------------
!Calculate growing degree days
!-----------------------------

!itb_crop...gdd flag to determine initial LAI on day that
!itb_crop...seeds emerge from ground

!   if(sib%diag%gdd > 100.0_dbl_kind) gdd_flag = 1


!EL...added to avoid gdd calculation before real planting date, 
!EL...since pd is printed out as 0 before the real planting 
!EL...date based on the above ndf60==5 criterion

	if (pd>0                  .AND.         & 
		time%doy >= pd        .AND.         &
        sib%diag%tempf>50.0   .AND.         &
        sib%diag%tempf<86.0)      then

    	sib%diag%gdd=sib%diag%gdd + sib%diag%tempf- 50.0_dbl_kind
	
	endif

	if (time%doy>300) then
        sib%diag%gdd=0.0001
	endif



!------------------------------
!	Reading and summing assimn
!-------------------------------

		assim_accum=0.0_dbl_kind
 
	do i0 = 1, sib%diag%tb_indx
    
      	assim_accum = assim_accum + (sib%diag%tb_assim(i0)*time%dtsib) 

!EL...multiplied by the no. secs per each timestep (i.e. tbsib) to convert assim mol sec-1 to mol

   enddo

   sib%diag%assim_d = assim_accum*12 !multiplied by 12 to convert mol to g


!-----------------------------------------------------------
! allocation sheme for fractions for assimilate partitioning
!-----------------------------------------------------------   
!EL...1-roots, 2-leaves,3-stems,4-products (flowers and grains)

		if (time%doy==(pd+10)) then
	        sib%diag%alloc(1)=0.5
			sib%diag%alloc(2)=0.4
			sib%diag%alloc(3)=0.1	
			sib%diag%alloc(4)=0.0	
        
        else if (time%doy>=(pd+10).and.time%doy<(pd+30)) then
            sib%diag%alloc(1)=0.5-(0.5-0.48)*(time%doy-(pd+10))/20
            sib%diag%alloc(2)=0.4-(0.4-0.3)*(time%doy-(pd+10))/20
            sib%diag%alloc(3)=0.1-(0.1-0.22)*(time%doy-(pd+10))/20
            sib%diag%alloc(4)=0.0
         
         
         else if (time%doy>=(pd+30).and.time%doy<(pd+40)) then
            sib%diag%alloc(1)=0.48-(0.48-0.48)*(time%doy-(pd+30))/10
			sib%diag%alloc(2)=0.3-(0.3-0.3)*(time%doy-(pd+30))/10
			sib%diag%alloc(3)=0.22-(0.22-0.22)*(time%doy-(pd+30))/10	
			sib%diag%alloc(4)=0.0   

        else if (time%doy>=(pd+40).and.time%doy<(pd+50)) then
            sib%diag%alloc(1)=0.48-(0.48-0.37)*(time%doy-(pd+40))/10
			sib%diag%alloc(2)=0.3-(0.3-0.4)*(time%doy-(pd+40))/10
			sib%diag%alloc(3)=0.22-(0.22-0.23)*(time%doy-(pd+40))/10	
			sib%diag%alloc(4)=0.0   

      
        else if (time%doy>=(pd+50).and.time%doy<(pd+60)) then
            sib%diag%alloc(1)=0.37-(0.37-0.24)*(time%doy-(pd+50))/10
			sib%diag%alloc(2)=0.4-(0.4-0.54)*(time%doy-(pd+50))/10
			sib%diag%alloc(3)=0.23-(0.22-0.22)*(time%doy-(pd+50))/10	
			sib%diag%alloc(4)=0.0    
  

        else if (time%doy>=(pd+60).and.time%doy<(pd+75)) then
            sib%diag%alloc(1)=0.24-(0.24-0.2)*(time%doy-(pd+60))/15
			sib%diag%alloc(2)=0.54-(0.54-0.51)*(time%doy-(pd+60))/15
			sib%diag%alloc(3)=0.22-(0.22-0.29)*(time%doy-(pd+60))/15	
			sib%diag%alloc(4)=0.0


        else if (time%doy>=(pd+75).and.time%doy<(pd+80)) then
            sib%diag%alloc(1)=0.2-(0.2-0.2)*(time%doy-(pd+75))/5
			sib%diag%alloc(2)=0.51-(0.51-0.2)*(time%doy-(pd+75))/5
			sib%diag%alloc(3)=0.29-(0.29-0.38)*(time%doy-(pd+75))/5	
			sib%diag%alloc(4)=0.0-(0.0-0.22)*(time%doy-(pd+75))/5 

        else if (time%doy>=(pd+80).and.time%doy<(pd+89)) then
            sib%diag%alloc(1)=0.2-(0.2-0.2)*(time%doy-(pd+80))/9
			sib%diag%alloc(2)=0.2-(0.2-0.2)*(time%doy-(pd+80))/9
			sib%diag%alloc(3)=0.38-(0.38-0.2)*(time%doy-(pd+80))/9	
			sib%diag%alloc(4)=0.22-(0.22-0.4)*(time%doy-(pd+80))/9 

        else if (time%doy>=(pd+89).and.time%doy<(pd+98)) then
            sib%diag%alloc(1)=0.2-(0.2-0.2)*(time%doy-(pd+89))/9
			sib%diag%alloc(2)=0.2-(0.2-0.0)*(time%doy-(pd+89))/9
			sib%diag%alloc(3)=0.2-(0.2-0.015)*(time%doy-(pd+89))/9	
			sib%diag%alloc(4)=0.4-(0.4-0.785)*(time%doy-(pd+89))/9

        else if (time%doy>=(pd+98).and.time%doy<(pd+108)) then
            sib%diag%alloc(1)=0.2-(0.2-0.1)*(time%doy-(pd+98))/10
			sib%diag%alloc(2)=0.0-(0.0-0.0)*(time%doy-(pd+98))/10
			sib%diag%alloc(3)=0.015-(0.015-0.0)*(time%doy-(pd+98))/10	
			sib%diag%alloc(4)=0.785-(0.755-0.93)*(time%doy-(pd+98))/10 

        else if (time%doy>=(pd+108).and.time%doy<(pd+115)) then
            sib%diag%alloc(1)=0.1-(0.1-0.1)*(time%doy-(pd+108))/7
			sib%diag%alloc(2)=0.0-(0.0-0.0)*(time%doy-(pd+108))/7
			sib%diag%alloc(3)=0.0-(0.0-0.0)*(time%doy-(pd+108))/7	
			sib%diag%alloc(4)=0.9-(0.9-0.9)*(time%doy-(pd+108))/7

        else if (time%doy>=(pd+115).and.time%doy<=(pd+121)) then
            sib%diag%alloc(1)=0.1-(0.1-0.1)*(time%doy-(pd+115))/6
			sib%diag%alloc(2)=0.0-(0.0-0.0)*(time%doy-(pd+115))/6
			sib%diag%alloc(3)=0.0-(0.0-0.0)*(time%doy-(pd+115))/6	
			sib%diag%alloc(4)=0.9-(0.9-0.9)*(time%doy-(pd+115))/6

         else if (time%doy>=(pd+121).and.time%doy<=(pd+130)) then
            sib%diag%alloc(1)=0.1-(0.1-0.1)*(time%doy-(pd+115))/6
			sib%diag%alloc(2)=0.0-(0.0-0.0)*(time%doy-(pd+115))/6
			sib%diag%alloc(3)=0.0-(0.0-0.0)*(time%doy-(pd+115))/6	
			sib%diag%alloc(4)=0.9-(0.9-0.9)*(time%doy-(pd+115))/6

		elseif ((time%doy<pd).and.(time%doy>(pd+130))) then
			sib%diag%alloc(1)=0.0
			sib%diag%alloc(2)=0.0	
			sib%diag%alloc(3)=0.0	
			sib%diag%alloc(4)=0.0
        endif
	

!----------------------------------
!Calculate total weight allocation	
!----------------------------------
!EL..w_main is the (drywt+maintenancerespiration) together
!EL..the coefficients given corresponds to growth resp & (drywt+maint R) ref: deVries et al. (1989); 
!EL...here w_main is the (drywt+maint.respn)    

     if((time%doy>pd+10).and.(time%doy<(pd+130))) then

       	sib%diag%w_main=sib%diag%assim_d/((sib%diag%alloc(1)*1.2929)+(sib%diag%alloc(2)*1.431)& 

            +(sib%diag%alloc(3)*1.2946)+(sib%diag%alloc(4)*1.6752))

    else

        sib%diag%w_main =0.0

    endif


!Calculate w_main allocation to different plant parts

!EL... clculates absolute allocation for roots(1),leaves(2),stems(3),and products(4) using allocation fraction for roots and assimilation

        sib%diag%allocwt(1)=sib%diag%w_main*sib%diag%alloc(1) 
        sib%diag%allocwt(2)=sib%diag%w_main*sib%diag%alloc(2) 
        sib%diag%allocwt(3)=sib%diag%w_main*sib%diag%alloc(3) 
        sib%diag%allocwt(4)=sib%diag%w_main*sib%diag%alloc(4) 

!---------------------------------------------    
!Calculate cumulative drywt.in each plant part 
!---------------------------------------------
!EL...calculates cumulative dry weight before respn

  do j=1,4	 
       
     	sib%diag%cum_wt(time%doy,j)=sib%diag%cum_wt(time%doy-1,j)+sib%diag%allocwt(j)
       
	  If ((time%doy-1)==0) then
     
        sib%diag%cum_wt(time%doy,j)=0
     
	  endif
 
  enddo
       
      sib%diag%cum_w(1)=sib%diag%cum_wt(time%doy,1)
	  sib%diag%cum_w(2)=sib%diag%cum_wt(time%doy,2)
	  sib%diag%cum_w(3)=sib%diag%cum_wt(time%doy,3)
	  sib%diag%cum_w(4)=sib%diag%cum_wt(time%doy,4)
	 
	  if (time%doy>300) then
		sib%diag%cum_wt(time%doy,1)=0.0
		sib%diag%cum_wt(time%doy,2)=0.0
		sib%diag%cum_wt(time%doy,3)=0.0
		sib%diag%cum_wt(time%doy,4)=0.0
	  endif


!----------------------------
!Calculate growth respiration
!----------------------------
 
        sib%diag%phen_growthr(1)=sib%diag%allocwt(1)*2*0.537*12/44
        sib%diag%phen_growthr(2)=sib%diag%allocwt(2)*2*0.790*12/44
        sib%diag%phen_growthr(3)=sib%diag%allocwt(3)*2*0.540*12/44
        sib%diag%phen_growthr(4)=sib%diag%allocwt(4)*2*1.238*12/44


!--------------------------
!Calculate maintanence resp
!--------------------------
   
       sib%diag%phen_maintr(1)=sib%diag%cum_wt(time%doy-1,1)*0.32*(0.03*2*12/44)*(2.0**((sib%diag%tempc-20.)/10.))   

!EL...Q10 coefficient is 1.8 for soybean and 2.0 for corn (Norman and Arkebauer, 1991);
!EL...0.32 is the nonstructural C fraction of root C needing maintenance
!EL...(calculations based on Allen et al., 1998 and Rogers et al., 2006)&
!EL...maint. coeff. info from Penning de Vries,1989, Amthor, 1984, and Norman and Arkebauer, 1991)

       sib%diag%phen_maintr(2)=sib%diag%cum_wt(time%doy-1,2)*0.38*0.03*2*0.75*12/44*(2.0**((sib%diag%tempc-20.)/10.))

!EL..multiplied by 0.75 to incorporate that maintenance respiration during the daytime is half that of  nighttime (de Vries et al., 1989)
!EL.. 0.38 is the nonstructural fraction of leaf C needing maintenance (calculations based on Allen et al., 1998 and Rogers et al., 2006)

       sib%diag%phen_maintr(3)=sib%diag%cum_wt(time%doy-1,3)*0.32*0.01*2*12/44*(2.0**((sib%diag%tempc-20.)/10.))

!EL... 0.32 is the nonstructural fraction of stem C  needing maintenance
!EL..(calculations based on Allen et al., 1998 and Rogers et al., 2006).

       sib%diag%phen_maintr(4)=sib%diag%cum_wt(time%doy-1,4)*0.46*0.015*2*12/44*(2.0**((sib%diag%tempc-20.)/10.))

!EL.. 0.4 is the nonstructural fraction of seed C  needing maintenance
!EL..(calculations based on Allen et al., 1998 and Rogers et al., 2006) 

	
	  if (time%doy>300) then
		sib%diag%phen_maintr(1)=0.0001
		sib%diag%phen_maintr(2)=0.0001
		sib%diag%phen_maintr(3)=0.0001
		sib%diag%phen_maintr(4)=0.0001
	  endif

!------------------------------
!Calculate dry weight change
!-----------------------------

	  sib%diag%wch(1)=sib%diag%allocwt(1)- sib%diag%phen_maintr(1)
      sib%diag%wch(2)=sib%diag%allocwt(2)- sib%diag%phen_maintr(2)
      sib%diag%wch(3)=sib%diag%allocwt(3)- sib%diag%phen_maintr(3)
      sib%diag%wch(4)=sib%diag%allocwt(4)- sib%diag%phen_maintr(4)

	 if (time%doy>300) then
		sib%diag%wch(1)=0.0001
		sib%diag%wch(2)=0.0001
		sib%diag%wch(3)=0.0001
		sib%diag%wch(4)=0.0001
	 endif

!--------------------------------------------------------------
!Recalculate final cumulative dry weight (g C m-2) of each plant part
!--------------------------------------------------------------

	  do j=1,4	 
       
     		sib%diag%cum_drywt(time%doy,j)=sib%diag%cum_drywt(time%doy-1,j)+sib%diag%wch(j)
       
	  	If ((time%doy-1)==0) then
     
        	sib%diag%cum_drywt(time%doy,j)=0
     
	  	endif
 
      enddo

!EL..Renamed to be output on a daily basis in a separate text file

	  sib%diag%final_drywt(1)=sib%diag%cum_drywt(time%doy,1) 
	  sib%diag%final_drywt(2)=sib%diag%cum_drywt(time%doy,2)
	  sib%diag%final_drywt(3)=sib%diag%cum_drywt(time%doy,3)
	  sib%diag%final_drywt(4)=sib%diag%cum_drywt(time%doy,4)

 	if (time%doy>300) then

		sib%diag%cum_drywt(time%doy,1)=0.0
		sib%diag%cum_drywt(time%doy,2)=0.0
		sib%diag%cum_drywt(time%doy,3)=0.0
		sib%diag%cum_drywt(time%doy,4)=0.0
	 
	endif

!------------------------------------------------
!final leaf weight (C g m-2) 
!-------------------------------------------------
!EL..after adjustment for senescence and harvest event
      
	 if (time%doy>=(pd+10).and.time%doy<(pd+75)) then
	       sib%diag%leafwt_c=sib%diag%cum_drywt(time%doy,2)
	
     elseif (time%doy>=(pd+75).and.time%doy<(pd+105)) then
	       sib%diag%leafwt_c=sib%diag%cum_drywt(time%doy,2)*0.85 
 
     elseif (time%doy>=(pd+105).and.time%doy<(pd+121)) then
	       sib%diag%leafwt_c=sib%diag%cum_drywt(time%doy,2)*0.8 
 	 
	 elseif (time%doy>=(pd+121).and.time%doy<(pd+130)) then
	       sib%diag%leafwt_c=sib%diag%cum_drywt(time%doy,2)*0.6

	 elseif (time%doy>=(pd+131).and.time%doy<(pd+144)) then
	       sib%diag%leafwt_c=sib%diag%cum_drywt(time%doy,2)*0.1

	 else
	       sib%diag%leafwt_c=sib%diag%cum_drywt(time%doy,2)*0.01

     endif

!--------------
!Calculate LAI
!-------------

      sib%diag%phen_LAI=sib%diag%leafwt_c*2*0.02 	

!EL...convert to dry weight g m-2 and then multiplied by SLA; 
!EL..SLA determined by the averages based on several studies



 	sib%diag%tb_indx = 0	 !at the end of each day tb_index is set to zero


!print*,pd,sib%diag%assim_d,sib%diag%phen_LAI

 	write(20,'(i4.4,2x,i3.3,2x,43(1x,f11.2))')time%year,   &
            time%doy,sib%diag%tempf,sib%diag%tempc,            &
            sib%diag%gdd,sib%diag%assim_d,sib%diag%alloc(1:4) ,&
            sib%diag%w_main,sib%diag%allocwt(1:4),             &
            sib%diag%cum_w(1:4),sib%diag%phen_growthr(1:4),    &
            sib%diag%phen_maintr(1:4),sib%diag%wch(1:4),       &
            sib%diag%final_drywt(1:4),sib%diag%leafwt_c,       &
            sib%diag%phen_LAI 
	


sib%diag%tb_indx = 0	!at the end of each day tb_index is set to zero

	 
end subroutine soy_phen


end subroutine crop_accum
