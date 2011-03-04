!============SUBROUTINE CROP_ACCUM===================
!kdcorbin, 02/11 - removed pd7_est and pdindx7

subroutine crop_accum(sib,time,timevar)

use kinds
use sibtype
use timetype
use sib_bc_module
use physical_parameters, only: tice    
use sib_const_module

implicit none

integer(kind=int_kind) :: i0,n,i,j,sibpt
real(kind=dbl_kind)    :: drate,   &
                          dgrowth,dgrowth_opt,max_wmain,assimd_new,x

real(kind=dbl_kind)  :: coeff  !variable for carbon/CO2 conversion coefficient

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

integer(kind=int_kind) :: biomeclass
real(kind=real_kind) :: crop_doy

!-------------------------------------------------------------
type(sib_t), dimension(subcount),intent(inout) :: sib
type(time_struct), intent(in) :: time
!-------------------------------------------------------------  

coeff = 2.0_dbl_kind * 12.0_dbl_kind / 44.0_dbl_kind

do sibpt = 1,subcount  

   if (sib(sibpt)%param%biome >= 20) then
      if (sib(sibpt)%diag%tb_indx .gt. 1) then
         sib(sibpt)%diag%ta_bar = sib(sibpt)%diag%tb_temp / &
                                                 sib(sibpt)%diag%tb_indx
         sib(sibpt)%diag%rstfac_d = sib(sibpt)%diag%tb_rst / &
                                                 sib(sibpt)%diag%tb_indx

         !kdcorbin, 02/11 - moved accumulated assimn calculate here
         !...multiplied by the number of seconds per each timestep
         !...to convert assim mol/sec to mol; multiplied by 12 to convert mol to g
         sib(sibpt)%diag%assim_d = sib(sibpt)%diag%tb_assim*time%dtsib*12.0
      endif

        sib(sibpt)%diag%tb_temp = 0.
        sib(sibpt)%diag%tb_rst = 0.
        sib(sibpt)%diag%tb_assim = 0.
        sib(sibpt)%diag%tb_indx = 0.

!EL-conversion of avg. daily temperature (ta_bar) from Kelvin to Fahrenheit
sib(sibpt)%diag%tempf=((sib(sibpt)%diag%ta_bar-273.15)*1.8)+32.0

!EL-conversion of avg. daily temperature (ta_bar) from Kelvin to Celcius
sib(sibpt)%diag%tempc=sib(sibpt)%diag%ta_bar - tice  !tice=273K

!EL-Calling for different phenology schemes based on the biome type

!......BIOMES 20-29 reserved for crops
!.....20 - corn
!.....21 - soy
!.....22 - winter wheat
!.....23 - spring wheat

   if(sib(sibpt)%param%biome >= 20.0) then
        !Change temporary biome class depending on fallow or crops
	!kdcorbin, 02/11
        if (sib(sibpt)%diag%phen_switch == 0) then
	    biomeclass = temp_biome_bare
        else
            biomeclass = temp_biome_crop
	endif
   else
        biomeclass = int(sib(sibpt)%param%biome)
   endif

    if(sib(sibpt)%param%biome == (19 + corn_num)) then  
    	call corn_phen	
    elseif(sib(sibpt)%param%biome == (19 + soy_num)) then
	call soy_phen
    elseif (sib(sibpt)%param%biome == (19 + wwheat_num) .or. &
                  sib(sibpt)%param%biome == (19 + swheat_num)) then
	call wheat_phen
    endif

!itb_crop...need to calculate veg params daily...

    tempaerovar = aerovar(:,:,biomeclass)
    temptran(1,1) = sib(sibpt)%param%tran(1,1)
    temptran(1,2) = sib(sibpt)%param%tran(1,2)
    temptran(2,1) = sib(sibpt)%param%tran(2,1)
    temptran(2,2) = sib(sibpt)%param%tran(2,2)

    tempref(1,1) = sib(sibpt)%param%ref(1,1)
    tempref(1,2) = sib(sibpt)%param%ref(1,2)
    tempref(2,1) = sib(sibpt)%param%ref(2,1)
    tempref(2,2) = sib(sibpt)%param%ref(2,2)

    !kdcorbin,02/11 - changed time from pmonth
    !crop_doy = time%mid_month(time%pmonth)
    crop_doy = time%doy

    if(sib(sibpt)%diag%phen_switch == 1 ) then
         timevar%lai = sib(sibpt)%diag%phen_lai

         !kdcorbin, 02/11 - check to make sure LAI doesn't drop too low
          if (sib(sibpt)%diag%phen_lai .lt. min_lai_crop) then
               timevar%lai = min_lai_crop
          endif

          call phen_mapper(                              &
              latsib(sibpt),                             &
              crop_doy,                &
              sib(sibpt)%param%vcover,                           &
              sib(sibpt)%param%chil,                             &
              temptran,                                   &
              tempref,                                    & 
              morphtab(biomeclass),                       &
              tempaerovar,                                &
              laigrid,                                    &
              fvcovergrid,                                &
              timevar)	

         sib(sibpt)%param%aparc1 = timevar%fpar
         sib(sibpt)%param%zlt1 = timevar%lai
         sib(sibpt)%param%green1 = timevar%green
         sib(sibpt)%param%z0d1 = timevar%zo
         sib(sibpt)%param%zp_disp1 = timevar%zp_disp
         sib(sibpt)%param%rbc1 = timevar%rbc
         sib(sibpt)%param%rdc1 = timevar%rdc
         sib(sibpt)%param%gmudmu1 = timevar%gmudmu

   else  !phen_switch .ne. 1
        sib(sibpt)%param%aparc1       = sib(sibpt)%param%aparc2
        sib(sibpt)%param%zlt1         = sib(sibpt)%param%zlt2
        sib(sibpt)%param%green1       = sib(sibpt)%param%green2
        sib(sibpt)%param%z0d1         = sib(sibpt)%param%z0d2
        sib(sibpt)%param%zp_disp1     = sib(sibpt)%param%zp_disp2
        sib(sibpt)%param%rbc1         = sib(sibpt)%param%rbc2
        sib(sibpt)%param%rdc1         = sib(sibpt)%param%rdc2
        sib(sibpt)%param%gmudmu1      = sib(sibpt)%param%gmudmu2
        sib(sibpt)%param%d13cresp1    = sib(sibpt)%param%d13cresp2

         !kdcorbin, 02/11 - changed from pmonth
         call mapper(                              &
              latsib(sibpt),                             &
              crop_doy,           &
              min_ndvi_crop,             &
              min_ndvi_crop,             &
              min_fvcov_crop,            &
              sib(sibpt)%param%chil,                     &
              temptran,                              &
              tempref,                               & 
              morphtab(biomeclass),                  &
              tempaerovar,                           &
              laigrid,                               &
              fvcovergrid,                           &
              timevar)

      endif  !phen_switch test

       sib(sibpt)%param%aparc2 = timevar%fpar
       sib(sibpt)%param%zlt2 = timevar%lai
       sib(sibpt)%param%green2 = timevar%green
       sib(sibpt)%param%z0d2 = timevar%zo
       sib(sibpt)%param%zp_disp2 = timevar%zp_disp
       sib(sibpt)%param%rbc2 = timevar%rbc
       sib(sibpt)%param%rdc2 = timevar%rdc
       sib(sibpt)%param%gmudmu2 = timevar%gmudmu

 endif  !crop type
enddo

contains

!----------------------------------------------------------------
subroutine corn_phen
!----------------------------------------------------------------

implicit none

!Local Variables
real(kind=dbl_kind) :: temp1, litter
integer(kind=int_kind) :: dapd
real(kind=dbl_kind) :: vmax_factor

!--------------------------
!Calculate the planting date
!---------------------------

!EL...sib(sibpt)%diag%ndf_opt= no. of days with avg. temperature above 57F 
!EL..(i.e. avg warm enough  temp for considering planting)

    if (sib(sibpt)%diag%tempf<57.0) then
       sib(sibpt)%diag%ndf_opt=0	    
    else
       sib(sibpt)%diag%ndf_opt=sib(sibpt)%diag%ndf_opt+1
   endif

   if (sib(sibpt)%diag%ndf_opt==7 .AND.  &
       sib(sibpt)%diag%pd_annual == 0)  then
       sib(sibpt)%diag%pd = time%doy
       sib(sibpt)%diag%pd_annual = 1
   endif

   if (sib(sibpt)%diag%pd > 0 .AND.    &
       (time%doy == (sib(sibpt)%diag%pd+1)  .OR.  &
        time%doy == (sib(sibpt)%diag%pd+2)  .OR.  &
        time%doy == (sib(sibpt)%diag%pd+3)  .OR.  &
        time%doy == (sib(sibpt)%diag%pd+4)  .OR.  & 
        time%doy == (sib(sibpt)%diag%pd+5)  .OR.  &
        time%doy == (sib(sibpt)%diag%pd+6)  .OR.  &
        time%doy == (sib(sibpt)%diag%pd+7)) .AND.  &
        sib(sibpt)%diag%tempf < 53.0) then

           sib(sibpt)%diag%gdd = 0.0

   endif

!----------------------------
!Calculate growing degree days
!-----------------------------

!EL.. emergence at GDD=100.0
!EL.. sib(sibpt)%diag%nd_emerg= no. of days since emergence

        if (sib(sibpt)%diag%gdd<100.0_dbl_kind) then
            sib(sibpt)%diag%nd_emerg=0			
       elseif (sib(sibpt)%diag%gdd>=100.0_dbl_kind) then
            sib(sibpt)%diag%nd_emerg=sib(sibpt)%diag%nd_emerg+1
       endif

        if (sib(sibpt)%diag%nd_emerg ==1) then
           sib(sibpt)%diag%emerg_d = time%doy
        endif
  
!EL...added to avoid gdd calculation before real planting date, 
!EL...since pd is printed out as 0 before the real planting 
!EL...date based on the above ndf_opt criterion
!EL...Calculation of GDDs occurs between 50 and 86 F

	if (sib(sibpt)%diag%pd > 0                     .AND.         & 
	     time%doy  >=  sib(sibpt)%diag%pd   .AND.         &
             sib(sibpt)%diag%tempf > 50.0          .AND.         &
             sib(sibpt)%diag%tempf < 86.0)          then

  	     sib(sibpt)%diag%gdd=sib(sibpt)%diag%gdd + &
                          sib(sibpt)%diag%tempf- 50.0_dbl_kind
	endif


!itb_crop...harvest: reset
        !kdcorbin, 02/11 - added check for LAI and phen_switch
        !   If below the min value, then harvest   
	if (((time%doy>sib(sibpt)%diag%pd+175) .or. &
            ((time%doy>sib(sibpt)%diag%pd+125) .and. &
             (sib(sibpt)%param%zlt .lt. min_lai_crop))) .and. &
             (sib(sibpt)%diag%phen_switch == 1)) then

             sib(sibpt)%diag%gdd=0.0001
             sib(sibpt)%diag%pd=0
             sib(sibpt)%diag%w_main=0.0001
             sib(sibpt)%diag%w_main_pot=0.0001
             sib(sibpt)%diag%assim_d=0.0001
             sib(sibpt)%diag%phen_switch=0

             call set_ti(sib(sibpt))

             do i=1,4
                  sib(sibpt)%diag%phen_maintr(i) = 0.0001
                  sib(sibpt)%diag%cum_wt(i) = 0.0001
                  sib(sibpt)%diag%wch(i) = 0.0001
                  sib(sibpt)%diag%cum_drywt(i) = 0.0001
            enddo
	endif

!-----------------------------------------------------------
! allocation sheme for fractions for assimilate partitioning
!-----------------------------------------------------------   
!EL...1-roots, 2-leaves,3-stems,4-products (flowers and grains)
!EL...allocation to different growth stages after emergence 
!......(each stage given as a range of GDDs below) 
        if(sib(sibpt)%diag%gdd < 100.0) then
               sib(sibpt)%diag%alloc(:) = 0.0
	elseif(sib(sibpt)%diag%gdd <500.0) then
		sib(sibpt)%diag%alloc(1)=0.5
		sib(sibpt)%diag%alloc(2)=0.25
		sib(sibpt)%diag%alloc(3)=0.25	
		sib(sibpt)%diag%alloc(4)=0.0	
        elseif(sib(sibpt)%diag%gdd <1000.0) then
                sib(sibpt)%diag%alloc(1)=0.5-0.2*(sib(sibpt)%diag%gdd-500.0)/500.0
   		sib(sibpt)%diag%alloc(2)=0.25+0.1*(sib(sibpt)%diag%gdd-500.0)/500.0
		sib(sibpt)%diag%alloc(3)=0.25+0.1*(sib(sibpt)%diag%gdd-500.0)/500.0
		sib(sibpt)%diag%alloc(4)=0.0
        elseif(sib(sibpt)%diag%gdd<1180.0) then
                sib(sibpt)%diag%alloc(1)=0.3
		sib(sibpt)%diag%alloc(2)=0.35
		sib(sibpt)%diag%alloc(3)=0.35
		sib(sibpt)%diag%alloc(4)=0.0
        elseif(sib(sibpt)%diag%gdd<1360.0) then
                sib(sibpt)%diag%alloc(1)=0.3-0.1*(sib(sibpt)%diag%gdd-1180.0)/180.
		sib(sibpt)%diag%alloc(2)=0.35+0.1*(sib(sibpt)%diag%gdd-1180.0)/180.	
		sib(sibpt)%diag%alloc(3)=0.35-0.1*(sib(sibpt)%diag%gdd-1180.0)/180.	
		sib(sibpt)%diag%alloc(4)=0.15*(sib(sibpt)%diag%gdd-1180.0)/180.
	elseif(sib(sibpt)%diag%gdd<1660.0)then
                sib(sibpt)%diag%alloc(1)=0.2-0.15*(sib(sibpt)%diag%gdd-1360.0)/300.
		sib(sibpt)%diag%alloc(2)=0.45-0.44*(sib(sibpt)%diag%gdd-1360.0)/300.	
		sib(sibpt)%diag%alloc(3)=0.2-0.15*(sib(sibpt)%diag%gdd-1360.0)/300.	
		sib(sibpt)%diag%alloc(4)=0.15+0.74*(sib(sibpt)%diag%gdd-1360.0)/300.
	elseif(sib(sibpt)%diag%gdd<2730.0)then
                sib(sibpt)%diag%alloc(1)=0.05
		sib(sibpt)%diag%alloc(2)=0.01-0.01*(sib(sibpt)%diag%gdd-1660.0)/1070.
		sib(sibpt)%diag%alloc(3)=0.05-0.05*(sib(sibpt)%diag%gdd-1660.0)/1070.
		sib(sibpt)%diag%alloc(4)=0.89+0.06*(sib(sibpt)%diag%gdd-1660.0)/1070.	
	else
		sib(sibpt)%diag%alloc(:)=0.0
        endif

!----------------------------------
!Calculate total weight allocation	
!----------------------------------
!EL...w_main is the daily amount of dry matter added, which is the basis for 
!EL...calculation of growth and maintenance respn

!EL...Multiplication factors below were derived using the info taken 
!EL...from past studies; mainly from de Vries et al., 1989.

!EL...seedling weight at emergence  was determined based on 
!.......Richardson and Bacon (1993), Pinhero and Fletcher (1994), and Horri et al., 2007.
!EL...also considering a planting density most commonly used (i.e. row 
!......spacing- 30", and plant spacing 6")
!EL.. Carbon amount (0.37 g C m-2) was derived by multiplying the seedling 
!.......weight by 0.43

!EL..Calculating the w_main by using  assim, growth resp coefficients and 
!EL..allocation fractions from the original scheme

     if ((sib(sibpt)%diag%gdd >= 100.0) .AND. &
              (sib(sibpt)%diag%gdd < 2730.0)) then

            sib(sibpt)%diag%w_main = sib(sibpt)%diag%assim_d /   &
                  ((sib(sibpt)%diag%alloc(1) * 1.2214) +     &
                   (sib(sibpt)%diag%alloc(2) * 1.2515) +     & 
                   (sib(sibpt)%diag%alloc(3) * 1.2225) +     &
                   (sib(sibpt)%diag%alloc(4) * 1.2095))

    endif

!EL...considering the fact that early corn growth is linearly related 
!EL...to temperature(Muchow and Carberry, 1989) &
!EL...considering total daily growth, based on the growth rate 
!EL...information given in de Vries et al.1989...

!EL...basic daily growth rate for the initial growth phase, depending 
!EL...on the no. of days (and GDDs) that the seedling would depend 
!EL...on the seed carbon reserves

!EL...drate, in g m-2, is based on the range of results from original SiB 
!EL...simulations from the past and Goudriaan and Liar (1994)

!EL...Temperatures and relevant fractions of the basic 
!EL...daily growth rates based on temperature 
!EL...were set based on the growth rate info in de Vries et al. 1989, 
!EL...and slightly modified  by looking at the observed data.

 if (sib(sibpt)%diag%gdd>100.0 .AND. sib(sibpt)%diag%gdd<=190.0) then 
	drate=0.14
	max_wmain=8.0
        dgrowth_opt=(max_wmain-0.37)*drate

	if (sib(sibpt)%diag%tempc<8.0) then
               dgrowth = dgrowth_opt * sib(sibpt)%diag%rstfac_d * 0.01
	elseif (sib(sibpt)%diag%tempc < 14.0) then
                dgrowth = dgrowth_opt * sib(sibpt)%diag%rstfac_d * &
                              (0.01+.19 * (sib(sibpt)%diag%tempc-8)/6.)
	elseif (sib(sibpt)%diag%tempc < 19.0) then
                 dgrowth = dgrowth_opt * sib(sibpt)%diag%rstfac_d * &
                            (0.2 + 0.4 * (sib(sibpt)%diag%tempc-14.)/5.)
	elseif (sib(sibpt)%diag%tempc < 28.0) then
                  dgrowth = dgrowth_opt * sib(sibpt)%diag%rstfac_d * &
                                (0.6+0.4 * (sib(sibpt)%diag%tempc-19.)/9.)
	elseif (sib(sibpt)%diag%tempc < 35.0) then
                  dgrowth = dgrowth_opt * sib(sibpt)%diag%rstfac_d * &
                               (1.0 - 0.1 * (sib(sibpt)%diag%tempc-28.)/7.)
	elseif (sib(sibpt)%diag%tempc<45.0) then
                  dgrowth = dgrowth_opt * sib(sibpt)%diag%rstfac_d * &
                                (0.9-0.89 * (sib(sibpt)%diag%tempc-35.)/10.)
	endif          

        sib(sibpt)%diag%w_main_pot=sib(sibpt)%diag%w_main_pot+dgrowth

        if (time%doy == sib(sibpt)%diag%emerg_d) then
           sib(sibpt)%diag%w_main_pot=0.37+dgrowth
           sib(sibpt)%diag%w_main=sib(sibpt)%diag%w_main_pot
        endif

        if (sib(sibpt)%diag%gdd < 150.0) then
           sib(sibpt)%diag%w_main=max(sib(sibpt)%diag%w_main_pot, &
                                                           sib(sibpt)%diag%w_main) 
        endif

        !EL..back calculation of the new assim_d:
	assimd_new = sib(sibpt)%diag%w_main *       &
                ((sib(sibpt)%diag%alloc(1) * 1.2214) +     &
                ( sib(sibpt)%diag%alloc(2) * 1.2515) +     & 
                ( sib(sibpt)%diag%alloc(3) * 1.2225) +     &
                ( sib(sibpt)%diag%alloc(4) * 1.2095))
        sib(sibpt)%diag%assim_d = assimd_new
endif  !gdd > 1 and gdd < 190

!EL...Calculate w_main allocation (i.e. allocwt) to different plant parts
!EL.. calculates absolute allocation for roots(1),leaves(2),
!EL...stems(3), and products(4)
 
   do i = 1,4
        sib(sibpt)%diag%allocwt(i) = sib(sibpt)%diag%w_main * &
                                                     sib(sibpt)%diag%alloc(i) 
   enddo

!--------------------------
!Calculate maintanence respiration
!--------------------------

!EL...0.18 is the nonstructural C fraction of root C needing 
!EL...maintenance(calculations based on Brouquisse et al., 1998)
!EL...maint. coeff. info from Penning de Vries,1989, Amthor, 
!EL...1984, and Norman and Arkebauer, 1991).
!EL...final values can be also found in Lokupitiya et al., 2009)

       temp1 = (0.03 * coeff) *                   &
                  (2.0**((sib(sibpt)%diag%tempc - 20.0) / 10.0))
   
       sib(sibpt)%diag%phen_maintr(1) = sib(sibpt)%diag%cum_wt(1)       &
             * 0.18 * temp1  

!EL...(de Vries et al., 1989)
!EL.. 0.27 is the nonstructural fraction of leaf C needing 
!EL...maintenance (calculations based on Brouquisse et al., 1998)

       sib(sibpt)%diag%phen_maintr(2) = sib(sibpt)%diag%cum_wt(2) *     &
              0.27 * temp1

!EL... 0.24 is the nonstructural fraction of stem C  needing 
!EL...maintenance(calculations based on Brouquisse et al., 1998).

       temp1 = 0.01 *coeff * (2.0**((sib(sibpt)%diag%tempc-20.0) / 10.0))

       sib(sibpt)%diag%phen_maintr(3) = sib(sibpt)%diag%cum_wt(3)          &
              * 0.24 * temp1

!EL.. 0.7 is the nonstructural fraction of seed C  needing maintenance 
!EL...(calculations based on Thornton et al., 1969,Beauchemin et al., 1997)

       temp1 = 0.015 * coeff *      &
               (2.0 ** ((sib(sibpt)%diag%tempc - 20.0) / 10.0))

       sib(sibpt)%diag%phen_maintr(4) = sib(sibpt)%diag%cum_wt(4) * 0.7 * temp1

!---------------------------------------------    
!Calculate cumulative drywt.in each plant part 
!---------------------------------------------
!EL...calculates cumulative new biomass, which will be used in calculating
!..... the following day's maint respn

	 do j=1,4	 
             sib(sibpt)%diag%cum_wt(j) = sib(sibpt)%diag%cum_wt(j) + &
                                                           sib(sibpt)%diag%allocwt(j)
         enddo

!----------------------------
!Calculate growth respiration
!----------------------------
 
        sib(sibpt)%diag%phen_growthr(1) = &
                sib(sibpt)%diag%allocwt(1)*0.406 * coeff
        sib(sibpt)%diag%phen_growthr(2) = &
                sib(sibpt)%diag%allocwt(2)*0.461 * coeff
        sib(sibpt)%diag%phen_growthr(3) = &
               sib(sibpt)%diag%allocwt(3)*0.408 * coeff
        sib(sibpt)%diag%phen_growthr(4) = &
               sib(sibpt)%diag%allocwt(4)*0.384 * coeff

!------------------------------
!Calculate dry weight change
!-----------------------------

      do i = 1,4
          sib(sibpt)%diag%wch(i) = sib(sibpt)%diag%allocwt(i) - &
                                                sib(sibpt)%diag%phen_maintr(i)
      enddo

!--------------------------------------------------------------
!Calculate final cumulative dry weight (g C m-2) of each plant part
!--------------------------------------------------------------
	  do j=1,4	 
                sib(sibpt)%diag%cum_drywt(j) = sib(sibpt)%diag%cum_drywt(j) + &
                                                                sib(sibpt)%diag%wch(j)
           enddo
	  
          sib(sibpt)%diag%tot_biomass = sib(sibpt)%diag%cum_drywt(2) + &
                     sib(sibpt)%diag%cum_drywt(3) + sib(sibpt)%diag%cum_drywt(4)

!---------------------------------------------------
!
!Decreasing the vmax during the seed filling stage and senescence down
!  to 60% of original vmax value.  
!This is supported in Crafts-Brandner et al. (Plant Physiol., 1992), 
!   Jiang et al. (Plant Physiol., 1993), and Jiang et al. (Photosynthesis Research, 1997)
!
!kdcorbin, 02/11
!---------------------------------------------------

dapd = time%doy - sib(sibpt)%diag%pd
if (sib(sibpt)%diag%phen_switch == 1 .and.  &
    time%doy > sib(sibpt)%diag%pd + vmax_start(corn_num) .and. &
    time%doy < sib(sibpt)%diag%pd + vmax_stop(corn_num)) then
          vmax_factor = crop_vmax0a(corn_num) - crop_vmax0b(corn_num)
          sib(sibpt)%param%vmax0(2) = crop_vmax0a(corn_num) - &
                  vmax_factor * (dapd - vmax_start(corn_num)) / &
                    (vmax_stop(corn_num)-vmax_start(corn_num)-1)
endif

!------------------------------------------------
!final leaf weight (C g m-2) 
!-------------------------------------------------
!EL...i.e. after adjustment for senescence and harvest event
 
      if( sib(sibpt)%diag%gdd < 0.0001 ) then

           sib(sibpt)%diag%leafwt_c =  0.01 
      
      elseif (sib(sibpt)%diag%gdd < 2650) then

                sib(sibpt)%diag%leafwt_c = sib(sibpt)%diag%cum_drywt(2)
	
      elseif (sib(sibpt)%diag%gdd < 2730.0) then

                sib(sibpt)%diag%leafwt_c = 0.95 * sib(sibpt)%diag%cum_drywt(2) -   &
                         (0.95-0.5) * sib(sibpt)%diag%cum_drywt(2) *      &
                         ((sib(sibpt)%diag%gdd - 2650.0) / 80.0)

                !EL...allowing for carbon remobilization from senescing leaves
                !......to growing products
  
                sib(sibpt)%diag%cum_drywt(4)=sib(sibpt)%diag%cum_drywt(4)+ &
                       ((sib(sibpt)%diag%cum_drywt(2)- sib(sibpt)%diag%leafwt_c))*0.03
                litter=((sib(sibpt)%diag%cum_drywt(2)- sib(sibpt)%diag%leafwt_c))*0.97

       elseif (sib(sibpt)%diag%gdd < 2900.0) then 
         
          sib(sibpt)%diag%leafwt_c = 0.5 * sib(sibpt)%diag%cum_drywt(2) -   &
          (0.5-0.01) * sib(sibpt)%diag%cum_drywt(2) *      &
                              ((sib(sibpt)%diag%gdd - 2730.0) / 170.0)

          litter=((sib(sibpt)%diag%cum_drywt(2)- sib(sibpt)%diag%leafwt_c))

      !EL...introducing harvest towards  the end of field drying     
      elseif( sib(sibpt)%diag%gdd >= 2900.0) then
          sib(sibpt)%diag%leafwt_c=0.0001
          litter=sib(sibpt)%diag%cum_drywt(1) + sib(sibpt)%diag%cum_drywt(2) &
                    +sib(sibpt)%diag%cum_drywt(3)
      endif

!--------------
!Calculate LAI
!-------------
!EL...convert to dry weight g m-2 and then multiply by SLA; 	
!EL...SLA determined by the averages based on several studies

      sib(sibpt)%diag%phen_LAI=sib(sibpt)%diag%leafwt_c * 2.0 * 0.02

!EL...since the SLA for irrigated crops is higher, the following will be used
!EL...for irrigated corn (very slight increase, as no published data could
!EL...be found for corn).
     !sib(sibpt)%diag%phen_LAI=sib(sibpt)%diag%leafwt_c * 2.0 * 0.021

!itb_crop...at the moment that growing degree days (gdd) passes
!itb_crop...100, we will initialize the LAI

    if (sib(sibpt)%diag%gdd >= 100.0_dbl_kind .AND. &
       time%doy >= sib(sibpt)%diag%emerg_d .AND. &
       !kdcorbin, 02/11 - only set switch once
       sib(sibpt)%diag%phen_switch == 0) then
            sib(sibpt)%diag%phen_switch = 1
            call set_ti(sib(sibpt))
     endif

end subroutine corn_phen


!----------------------------------------------------------------
subroutine soy_phen
!----------------------------------------------------------------

implicit none

!Local Variables

real(kind=dbl_kind) :: temp1, litter
integer(kind=int_kind) :: dapd
real(kind=dbl_kind) :: vmax_factor

!--------------------------
!Calculate the planting date
!---------------------------

!EL...sib%diag%ndf_opt = no. of days with avg. temperature above 67F

	if (sib(sibpt)%diag%tempf < 66.5) then
              sib(sibpt)%diag%ndf_opt = 0			
	else   	
           sib(sibpt)%diag%ndf_opt = sib(sibpt)%diag%ndf_opt + 1
	endif

	if (sib(sibpt)%diag%ndf_opt == 7 .AND. &
            sib(sibpt)%diag%pd_annual == 0) then
              sib(sibpt)%diag%pd = time%doy
              sib(sibpt)%diag%pd_annual = 1
	endif

   if (sib(sibpt)%diag%pd > 0 .AND.    &
       (time%doy == (sib(sibpt)%diag%pd+1)  .OR.  &
        time%doy == (sib(sibpt)%diag%pd+2)  .OR.  &
        time%doy == (sib(sibpt)%diag%pd+3)  .OR.  &
        time%doy == (sib(sibpt)%diag%pd+4)  .OR.  &
        time%doy == (sib(sibpt)%diag%pd+5)) .AND.  &
        sib(sibpt)%diag%tempf < 50.0) then
          sib(sibpt)%diag%gdd = 0.0
   endif

!----------------------------
!Calculate growing degree days
!-----------------------------

!EL...added to avoid gdd calculation before real planting date, 
!EL...since pd is printed out as 0 before the real planting 
!EL...date based on the above ndf_opt criterion

	if (sib(sibpt)%diag%pd    >  0                .AND.         & 
            time%doy >= sib(sibpt)%diag%pd     .AND.         &
            sib(sibpt)%diag%tempf >  50.0         .AND.         &
            sib(sibpt)%diag%tempf <  86.0  )      then

             sib(sibpt)%diag%gdd=sib(sibpt)%diag%gdd + &
                     sib(sibpt)%diag%tempf - 50.0_dbl_kind
	endif

!EL...to avoid gdd calculation after harvesting is done, allowing 
!EL...ample time between planting and harvesting

        !kdcorbin, 02/11 - moved all harvest resets here
        !   and added phen_lai and phen_switch checks
	if (((time%doy > sib(sibpt)%diag%pd + 160) .or. &
             (  time%doy > sib(sibpt)%diag%pd + 120 .and. & 
                sib(sibpt)%diag%phen_lai < min_lai_crop )) .and. &
                sib(sibpt)%diag%phen_switch == 1) then
              sib(sibpt)%diag%gdd = 0.0001
              sib(sibpt)%diag%w_main = 0.0001
              sib(sibpt)%diag%assim_d = 0.0001
              sib(sibpt)%diag%pd = 0
              sib(sibpt)%diag%phen_switch = 0
              sib(sibpt)%diag%phen_lai = min_lai_crop

              call set_ti(sib(sibpt))

              do i=1,4
                   sib(sibpt)%diag%phen_maintr(i) = 0.0001
                   sib(sibpt)%diag%cum_wt(i) = 0.0001
                   sib(sibpt)%diag%wch(i) = 0.0001
                   sib(sibpt)%diag%cum_drywt(i) = 0.0001
              enddo
	endif

!-----------------------------------------------------------
! allocation sheme for fractions for assimilate partitioning
! modified to account for daily varying gmudmu - kdcorbin, 02/10
!-----------------------------------------------------------   
!EL...1-roots, 2-leaves,3-stems,4-products (flowers and grains)

         dapd = time%doy - sib(sibpt)%diag%pd

         if (sib(sibpt)%diag%pd<=0) then
             !Do nothing, not planted yet
         else
             if ((time%doy    < sib(sibpt)%diag%pd+10)     .OR.     &
                 (time%doy    > sib(sibpt)%diag%pd+121))  then
			sib(sibpt)%diag%alloc(:)=0.0
             elseif (time%doy < (sib(sibpt)%diag%pd + 30)) then
                        sib(sibpt)%diag%alloc(1) = 0.5
			sib(sibpt)%diag%alloc(2) = 0.25
			sib(sibpt)%diag%alloc(3) = 0.25	
			sib(sibpt)%diag%alloc(4) = 0.0	
             else if (time%doy <  (sib(sibpt)%diag%pd + 40) )   then
                        sib(sibpt)%diag%alloc(1) = 0.5  - 0.05 * (dapd - 30) / 10.0
			sib(sibpt)%diag%alloc(2) = 0.25 + 0.05 * (dapd - 30) / 10.0
			sib(sibpt)%diag%alloc(3) = 0.25
			sib(sibpt)%diag%alloc(4) = 0.0   
              else if (time%doy    <  (sib(sibpt)%diag%pd + 50) )  then
                        sib(sibpt)%diag%alloc(1) = 0.45 - 0.15 * (dapd - 40) / 10.0
			sib(sibpt)%diag%alloc(2) = 0.3 + .15 * (dapd - 40) / 10.0
			sib(sibpt)%diag%alloc(3) = 0.25
			sib(sibpt)%diag%alloc(4) = 0.0   
               else if (time%doy    <  (sib(sibpt)%diag%pd + 60)) then
                        sib(sibpt)%diag%alloc(1) = 0.3 - 0.08 * (dapd - 50) / 10.0
   			sib(sibpt)%diag%alloc(2) = 0.45 + 0.08 * (dapd - 50) / 10.0
			sib(sibpt)%diag%alloc(3) = 0.25
			sib(sibpt)%diag%alloc(4) = 0.0    
               else if (time%doy    <  (sib(sibpt)%diag%pd + 75))  then
                        sib(sibpt)%diag%alloc(1) = 0.22
			sib(sibpt)%diag%alloc(2) = 0.53 - 0.13 * (dapd - 60) / 15.0
			sib(sibpt)%diag%alloc(3) = 0.25 + .13 * (dapd - 60) / 15.0
			sib(sibpt)%diag%alloc(4) = 0.0
               else if (time%doy    <  (sib(sibpt)%diag%pd + 80))  then
                        sib(sibpt)%diag%alloc(1) = 0.21 - 0.07 * (dapd - 75) / 5.0
			sib(sibpt)%diag%alloc(2) = 0.40 - 0.1 * (dapd - 75) / 5.0
			sib(sibpt)%diag%alloc(3) = 0.38 - 0.2 * (dapd - 75) / 5.0
			sib(sibpt)%diag%alloc(4) = 0.37 * (dapd - 75) / 5.0
                else if (time%doy    <  (sib(sibpt)%diag%pd + 89))    then
                        sib(sibpt)%diag%alloc(1) = 0.17
			sib(sibpt)%diag%alloc(2) = 0.3 - 0.1 * (dapd - 80) / 9.0
			sib(sibpt)%diag%alloc(3) = 0.18
			sib(sibpt)%diag%alloc(4) = 0.37 + 0.08 * (dapd - 80) / 9.0
                else if (time%doy <  (sib(sibpt)%diag%pd+98)) then
                        sib(sibpt)%diag%alloc(1) = 0.17 
			sib(sibpt)%diag%alloc(2) = 0.2
			sib(sibpt)%diag%alloc(3) = 0.18 - 0.165 * (dapd - 89) / 9.0
			sib(sibpt)%diag%alloc(4) = 0.45 + 0.365 * (dapd - 89) / 9.0
                 else if (time%doy    <  (sib(sibpt)%diag%pd+108)) then
                        sib(sibpt)%diag%alloc(1) = 0.17 - 0.12 * (dapd - 98) / 10.0
			sib(sibpt)%diag%alloc(2) = 0.0
			sib(sibpt)%diag%alloc(3) = 0.015 - 0.015 * (dapd - 98) / 10.0
			sib(sibpt)%diag%alloc(4) = 0.815 + 0.135 * (dapd - 98) / 10.0
                 else if (time%doy    <= (sib(sibpt)%diag%pd+121)) then
                        sib(sibpt)%diag%alloc(1) = 0.05 
			sib(sibpt)%diag%alloc(2) = 0.0 
			sib(sibpt)%diag%alloc(3) = 0.0 
			sib(sibpt)%diag%alloc(4) = 0.95
                 else
                        sib(sibpt)%diag%alloc(:) = 0.0     
                 endif  !doy tests
          endif !pd >= test
	
!----------------------------------
!Calculate total weight allocation	
!----------------------------------
!EL...w_main is the daily amount of dry matter added, which is the basis 
!EL...for calculation of growth and maintenance respn

!EL...Multiplication factors below were derived using the info taken 
!EL...from past studies; mainly from de Vries et al., 1989.

!EL..Calculating the w_main by using  assim, growth resp coefficients and 
!EL..allocation fractions from the original scheme

    if (sib(sibpt)%diag%pd  >  0                                .AND.     &
        (time%doy    <  (sib(sibpt)%diag%pd + 10)      .OR.     &
         time%doy    >= (sib(sibpt)%diag%pd + 160)))   then 
          sib(sibpt)%diag%w_main =0.0001
    endif 

    if (sib(sibpt)%diag%pd >  0                               .AND.     &
       time%doy    >  (sib(sibpt)%diag%pd + 10)      .AND.     &
       time%doy    <= (sib(sibpt)%diag%pd + 121))    then
      	sib(sibpt)%diag%w_main = sib(sibpt)%diag%assim_d /      &
                             ((sib(sibpt)%diag%alloc(1)* 1.2929) +    &
                              (sib(sibpt)%diag%alloc(2) * 1.431)  +    &
                              (sib(sibpt)%diag%alloc(3) * 1.2946) +    &
                              (sib(sibpt)%diag%alloc(4) *  1.6752))
      endif

!EL...considering the fact that early soybean growth is linearly 
!EL...related to temperature(Muchow and Carberry, 1989) &
!EL...considering total daily growth, based on the growth rate 
!EL...information given in de Vries et al.1989...
!EL...seedling weight at emergence was determined based on 
!EL...Green and Sudia (1969) and Smiciklas et al. (1992).
!EL.. Carbon amount (0.37 g C m-2) was derived by multiplying 
!EL...the seedling weight by 0.43

	if (sib(sibpt)%diag%pd > 0 .AND.time%doy == (sib(sibpt)%diag%pd + 10)) then
		sib(sibpt)%diag%w_main=0.26 
	endif

!EL...basic daily growth rate for the initial growth phase, 
!EL...depending on the average no. of days
!EL.. between planting and observed max LAI
!EL...initial phase when the growth depends on the seed and 
!EL...cotyledon C stores was considered to extend from 
!EL...pd to pd+21(McWilliams et al., 1999)

!EL...Temperatures and relevant fractions from the basic 
!EL...drates based on temperature 
!EL...were set based on the info in de Vries et al. 1989, 
!EL...and slightly modified  by looking 
!EL...at the observed data.

        if (sib(sibpt)%diag%pd >  0 .AND.  &
            time%doy > (sib(sibpt)%diag%pd +10)  .AND.      &
            time%doy    <= (sib(sibpt)%diag%pd + 21)) then 
		drate=0.091	
		max_wmain=8.0
                dgrowth_opt = (max_wmain-0.26) * drate

          	if (sib(sibpt)%diag%tempc<=8) then
                    dgrowth = dgrowth_opt * sib(sibpt)%diag%rstfac_d * 0.01
         	elseif (sib(sibpt)%diag%tempc<10.0) then
                    dgrowth = dgrowth_opt * sib(sibpt)%diag%rstfac_d * 0.01 &
                                + (0.24 * (sib(sibpt)%diag%tempc-8.)/2.)
              	elseif (sib(sibpt)%diag%tempc < 20.0) then
                     dgrowth = dgrowth_opt * sib(sibpt)%diag%rstfac_d * &
                                (0.25 + (0.65 * (sib(sibpt)%diag%tempc-10.)/10.))
	        elseif (sib(sibpt)%diag%tempc < 27.0) then
                      dgrowth = dgrowth_opt * sib(sibpt)%diag%rstfac_d * &
                               (0.9 + (0.15 * (sib(sibpt)%diag%tempc-20.)/7.))
             	elseif (sib(sibpt)%diag%tempc < 30.0) then
                      dgrowth = dgrowth_opt * sib(sibpt)%diag%rstfac_d * &
                                (1.05 + (.06 * (sib(sibpt)%diag%tempc-27.)/3.))
	        elseif (sib(sibpt)%diag%tempc < 40.0) then
                      dgrowth = dgrowth_opt * sib(sibpt)%diag%rstfac_d * &
                                (1.1 + (0.1 * (sib(sibpt)%diag%tempc - 30.0)/10.))
         	endif

        sib(sibpt)%diag%w_main_pot=sib(sibpt)%diag%w_main_pot+dgrowth
        sib(sibpt)%diag%w_main=max(sib(sibpt)%diag%w_main_pot, &
                                                        sib(sibpt)%diag%w_main)
 
        !EL..back calculation of the new assim_d:
	assimd_new= sib(sibpt)%diag%w_main *       &
               ((sib(sibpt)%diag%alloc(1) * 1.2929) +     &
                (sib(sibpt)%diag%alloc(2) * 1.431) +     & 
                (sib(sibpt)%diag%alloc(3) * 1.2946) +     &
                (sib(sibpt)%diag%alloc(4) * 1.6752))
        sib(sibpt)%diag%assim_d=assimd_new
 
       endif

!Calculate w_main allocation to different plant parts

!EL...calculates absolute allocation of biomass for roots(1),
!EL...leaves(2),stems(3),and products(4) using allocation fractions

        do i = 1,4
          sib(sibpt)%diag%allocwt(i) = sib(sibpt)%diag%w_main * sib(sibpt)%diag%alloc(i) 
        enddo

!--------------------------
!Calculate maintenance resp
!--------------------------
!EL...Q10 coefficient is 2.0 for soybean (Norman and Arkebauer, 1991);
!EL...0.32 is the nonstructural C fraction of root C needing maintenance
!EL...(calculations based on Allen et al., 1998 and Rogers et al., 2006)&
!EL...maint. coeff. info from Penning de Vries,1989, Amthor, 1984, 
!EL...and Norman and Arkebauer, 1991)
!EL...the final values could also be found in Lokupiitya et al., 2009

       temp1 = (0.03 * coeff) *                   &
             (1.8**((sib(sibpt)%diag%tempc - 20.0) / 10.0))
       sib(sibpt)%diag%phen_maintr(1) =  sib(sibpt)%diag%cum_wt(1)       &
             * 0.32 * temp1  

!EL...(de Vries et al., 1989)
!EL.. 0.38 is the nonstructural fraction of leaf C needing maintenance
!EL...(calculations based on Allen et al., 1998 and Rogers et al., 2006).

       sib(sibpt)%diag%phen_maintr(2) =  sib(sibpt)%diag%cum_wt(2) *     &
              0.38 * temp1

!EL... 0.32 is the nonstructural fraction of stem C  needing maintenance
!EL...(calculations based on Allen et al., 1998 and Rogers et al., 2006).

       temp1 = 0.01 * coeff * (1.8**((sib(sibpt)%diag%tempc-20.0) / 10.0))
       sib(sibpt)%diag%phen_maintr(3) =  sib(sibpt)%diag%cum_wt(3)          &
              * 0.32 * temp1

!EL.. 0.46 is the nonstructural fraction of seed C  needing maintenance 
!EL...(calculations based on Allen et al., 1998 and Rogers et al., 2006).

       temp1 = 0.015 * coeff * (1.8**((sib(sibpt)%diag%tempc-20.0)/10.0))
       sib(sibpt)%diag%phen_maintr(4)=  sib(sibpt)%diag%cum_wt(4) * 0.46*temp1

!---------------------------------------------    
!Calculate cumulative drywt.in each plant part 
!---------------------------------------------
!EL...calculates cumulative new biomass, which will be used in calculating 
!EL...the following day's maint respn

     do j=1,4	 
       	sib(sibpt)%diag%cum_wt(j) = sib(sibpt)%diag%cum_wt(j) + &
                                                     sib(sibpt)%diag%allocwt(j)
    enddo
       
!----------------------------
!Calculate growth respiration
!----------------------------

        sib(sibpt)%diag%phen_growthr(1) = sib(sibpt)%diag%allocwt(1) * 0.537 * coeff
        sib(sibpt)%diag%phen_growthr(2) = sib(sibpt)%diag%allocwt(2) * 0.790 * coeff
        sib(sibpt)%diag%phen_growthr(3) = sib(sibpt)%diag%allocwt(3) * 0.540 * coeff
        sib(sibpt)%diag%phen_growthr(4) = sib(sibpt)%diag%allocwt(4) * 1.238 * coeff

!------------------------------
!Calculate dry weight change
!-----------------------------

      do i = 1,4
        sib(sibpt)%diag%wch(i) = sib(sibpt)%diag%allocwt(i) - &
                                               sib(sibpt)%diag%phen_maintr(i)
      enddo

!--------------------------------------------------------------
!Recalculate final cumulative dry weight (g C m-2) of each plant part
!--------------------------------------------------------------

	  do j=1,4
	       sib(sibpt)%diag%cum_drywt(j) = sib(sibpt)%diag%cum_drywt(j) + &
                                                                  sib(sibpt)%diag%wch(j)
          enddo

          sib(sibpt)%diag%tot_biomass = sib(sibpt)%diag%cum_drywt(2) + & 
                 sib(sibpt)%diag%cum_drywt(3)+ sib(sibpt)%diag%cum_drywt(4)

!---------------------------------------------------
!
!Decreasing the vmax during the seed filling stage and senescence down
!  to 60% of original vmax value.  
!This is supported in Crafts-Brandner et al. (Plant Physiol., 1992), 
!   Jiang et al. (Plant Physiol., 1993), and Jiang et al. (Photosynthesis Research, 1997)
!
!kdcorbin, 02/11
!---------------------------------------------------

if (sib(sibpt)%diag%phen_switch == 1 .and.  &
    time%doy > sib(sibpt)%diag%pd + vmax_start(soy_num) .and. &
    time%doy < sib(sibpt)%diag%pd + vmax_stop(soy_num)) then
          vmax_factor = crop_vmax0a(soy_num) - crop_vmax0b(soy_num)
          sib(sibpt)%param%vmax0(1) = crop_vmax0a(soy_num) - &
                  vmax_factor * (dapd - vmax_start(soy_num)) / &
                    (vmax_stop(soy_num)-vmax_start(soy_num)-1)
endif

!------------------------------------------------
!final leaf weight (C g m-2) 
!-------------------------------------------------

	 if (sib(sibpt)%diag%pd > 0) then
              if (time%doy < (sib(sibpt)%diag%pd+10)) then
                     sib(sibpt)%diag%leafwt_c = sib(sibpt)%diag%cum_drywt(2) * 0.01
              elseif (time%doy <  (sib(sibpt)%diag%pd+60)) then
                     sib(sibpt)%diag%leafwt_c = sib(sibpt)%diag%cum_drywt(2)
	      elseif (time%doy   <  (sib(sibpt)%diag%pd+90))   then
                    sib(sibpt)%diag%leafwt_c = sib(sibpt)%diag%cum_drywt(2) -   &
                          0.15 * sib(sibpt)%diag%cum_drywt(2) * (dapd - 60) / 30.0

                    !EL...allowing for carbon remobilization from senescing leaves
                    !EL...to growing products
                    sib(sibpt)%diag%cum_drywt(4)=sib(sibpt)%diag%cum_drywt(4)+ &
                        (sib(sibpt)%diag%cum_drywt(2)*0.15)*0.1
        
                    !EL.. the dead leaf pool
                    litter=(sib(sibpt)%diag%cum_drywt(2)*0.15)*0.9
         
             elseif (time%doy    <  (sib(sibpt)%diag%pd+121))   then
               
	           sib(sibpt)%diag%leafwt_c = sib(sibpt)%diag%cum_drywt(2)*0.85 -   &
                        0.84 * sib(sibpt)%diag%cum_drywt(2) * (dapd - 90) / 31.0 

                   !EL...allowing for carbon remobilization from senescing leaves  
                   !EL...to growing products
  
                   sib(sibpt)%diag%cum_drywt(4)=sib(sibpt)%diag%cum_drywt(4) + &
                         ((sib(sibpt)%diag%cum_drywt(2)- sib(sibpt)%diag%leafwt_c))*0.03

                   !EL.. the dead leaf pool
                   litter=(sib(sibpt)%diag%cum_drywt(2)-sib(sibpt)%diag%leafwt_c)*0.97

            else
                 sib(sibpt)%diag%leafwt_c = 0.0001
                 litter=sib(sibpt)%diag%cum_drywt(1)+sib(sibpt)%diag%cum_drywt(2) &
                          +sib(sibpt)%diag%cum_drywt(3)
            endif
	 else
             sib(sibpt)%diag%leafwt_c = sib(sibpt)%diag%cum_drywt(2) * 0.01 
         endif

!--------------
!Calculate LAI
!-------------

!EL...convert to dry weight g m-2 and then multiply by SLA; 
!EL..SLA determined by the averages based on several studies

      !kdcorbin, 02/11 - added test for minimum value
      if (sib(sibpt)%diag%leafwt_c .lt. min_lai_crop) then
           sib(sibpt)%diag%phen_LAI = min_lai_crop
      else
           sib(sibpt)%diag%phen_LAI = sib(sibpt)%diag%leafwt_c * 2.0 * 0.025
      endif

!EL..for irrigated soybean, the SLA is higher; so the LAI for irrigated soils, is as follows.
       !sib(sibpt)%diag%phen_LAI = sib(sibpt)%diag%leafwt_c * 2.0 * 0.032

!itb_crop...turn on phenology model 10 days after planting day
     if (sib(sibpt)%diag%pd >  0  .AND.    &
          time%doy    >=  (sib(sibpt)%diag%pd + 10) .AND. &
          !kdcorbin, 02/11 - only turning on switch once
          sib(sibpt)%diag%phen_switch == 0) then
          sib(sibpt)%diag%phen_switch = 1
          call set_ti(sib(sibpt))
    endif

end subroutine soy_phen


!--------------------------------------------------------------
subroutine wheat_phen
!--------------------------------------------------------------

implicit none

real(kind=dbl_kind) :: temp1
real(kind=dbl_kind) :: vmax_factor

!EL...calculation of the planting date of winterwheat was set to start 
!EL...after Aug 15 (by looking at USDA planting dates), and the
!EL..optimum temperature for germination and emergence is between 
!EL...15 and 25C (Lindstrom et al., 1976; Burleigh et al., 1965)
!EL...sib(sibpt)%diag%ndf_opt= no. of days withe avg. temperature between 20 and 25C

if (sib(sibpt)%param%biome == 22) then

    if (time%doy>=227 .and. sib(sibpt)%diag%tempc<20.0) then
         sib(sibpt)%diag%ndf_opt=0		
    elseif (time%doy>=227 .and. sib(sibpt)%diag%tempc>=20.0 .and. &
              sib(sibpt)%diag%tempc<25.0) then
         sib(sibpt)%diag%ndf_opt=sib(sibpt)%diag%ndf_opt+1
   endif

   if (sib(sibpt)%diag%ndf_opt==7 .AND. sib(sibpt)%diag%pd == 0)  then
        sib(sibpt)%diag%pd = time%doy
   endif

    if (sib(sibpt)%diag%tempc <= -20.0) then
         sib(sibpt)%diag%gdd=0.0
    endif

endif  !biome == 22

!EL...calculation of the planting date of spring wheat was set to start after 
!EL...April 01 (by looking at USDA planting dates), and the 
!EL..optimum temperature for planting is 12C (Canadian Organic Growers, Inc., 1992)
!EL...sib(sibpt)%diag%ndf_opt= no. of days withe avg. temperature between 12 and 25C

if (sib(sibpt)%param%biome == 23) then

    if (time%doy>=91 .and. sib(sibpt)%diag%tempc<12.0) then
        sib(sibpt)%diag%ndf_opt=0		
   elseif (time%doy>=91 .and. sib(sibpt)%diag%tempc>=12.0 .and.  &
              sib(sibpt)%diag%tempc<25.0) then
        sib(sibpt)%diag%ndf_opt=sib(sibpt)%diag%ndf_opt+1
   endif

   if (sib(sibpt)%diag%ndf_opt==7 .AND. sib(sibpt)%diag%pd == 0)  then
        sib(sibpt)%diag%pd = time%doy
   endif

    if (sib(sibpt)%diag%tempc <= -20.0) then
         sib(sibpt)%diag%gdd=0.0
    endif

endif  !biome == 23

!----------------------------
!Calculate growing degree days
!-----------------------------      

        if (sib(sibpt)%diag%gdd<105.0_dbl_kind) then
           sib(sibpt)%diag%nd_emerg=0   !nd_emerg= no. of days since emergence
        else
           sib(sibpt)%diag%nd_emerg=sib(sibpt)%diag%nd_emerg+1

           if (sib(sibpt)%diag%nd_emerg == 1) then
              sib(sibpt)%diag%emerg_d = time%doy
           endif
        endif

!EL...added to avoid gdd calculation before real planting date, 
!EL...since pd is printed out as 0 before the real planting 
!EL...date based on the above ndf_opt criterion

if (sib(sibpt)%param%biome == 22) then
        if (sib(sibpt)%diag%doy==1) then
	    	sib(sibpt)%diag%gdd=769.0
        elseif (sib(sibpt)%diag%pd == 0 .and. time%doy == 2) then
                sib(sibpt)%diag%gdd=769.0
	elseif (sib(sibpt)%diag%pd > 0  .AND.   &
            sib(sibpt)%diag%tempc > 0.0  .AND. sib(sibpt)%diag%tempc < 26.0)      then
                    sib(sibpt)%diag%gdd=sib(sibpt)%diag%gdd + sib(sibpt)%diag%tempc
         elseif (sib(sibpt)%diag%pd==0 .and. time%doy>2  .AND.  &
                   sib(sibpt)%diag%gdd>500.0 .AND. &
                   sib(sibpt)%diag%tempc > 0.0 .AND. sib(sibpt)%diag%tempc < 26.0 ) then
             sib(sibpt)%diag%gdd=sib(sibpt)%diag%gdd + sib(sibpt)%diag%tempc
        endif
endif  !biome == 22

if (sib(sibpt)%param%biome == 23) then
	if (sib(sibpt)%diag%pd > 0                     .AND.         & 
            time%doy  >=  sib(sibpt)%diag%pd   .AND.         &
            sib(sibpt)%diag%tempc > 0.0 .AND. sib(sibpt)%diag%tempc < 21.11) then
             	     sib(sibpt)%diag%gdd=sib(sibpt)%diag%gdd + &
                                                       sib(sibpt)%diag%tempc- 0.0_dbl_kind
        endif

        if (sib(sibpt)%diag%gdd >= 215.0) then
            if (sib(sibpt)%diag%pd > 0                     .AND.         & 
		time%doy  >=  sib(sibpt)%diag%pd    .AND.         &
                sib(sibpt)%diag%tempc > 0.0 .AND. sib(sibpt)%diag%tempc < 35.0 ) then
                     sib(sibpt)%diag%gdd=sib(sibpt)%diag%gdd + sib(sibpt)%diag%tempc
            endif
        endif
endif  !biome == 23

!kdcorbin, 02/11 - crop harvest
if (sib(sibpt)%diag%gdd .gt. 2300.) then
    sib(sibpt)%diag%gdd = 0.
    sib(sibpt)%diag%w_main = 0.0001
    sib(sibpt)%diag%w_main_pot = 0.0001
    sib(sibpt)%diag%assim_d = 0.0001
    sib(sibpt)%diag%pd = 0
    sib(sibpt)%diag%phen_switch = 0
    sib(sibpt)%diag%leafwt_c = 0.0001

    call set_ti(sib(sibpt))

    do i=1,4
         sib(sibpt)%diag%cum_drywt(i) = 0.0001
         sib(sibpt)%diag%phen_maintr(i) = 0.0001
         sib(sibpt)%diag%cum_wt(i) = 0.0001
         sib(sibpt)%diag%wch(i) = 0.0001
    enddo
endif

!-----------------------------------------------------------
! allocation sheme for fractions for assimilate partitioning
!-----------------------------------------------------------   
!EL...1-roots, 2-leaves,3-stems,4-products (flowers and grains)

	if(sib(sibpt)%diag%gdd < 105.0 .or. sib(sibpt)%diag%gdd > 2300.0) then
               sib(sibpt)%diag%alloc(:) = 0.0
        elseif (sib(sibpt)%diag%gdd <350.0) then
		sib(sibpt)%diag%alloc(1)=0.4
		sib(sibpt)%diag%alloc(2)=0.4
		sib(sibpt)%diag%alloc(3)=0.2	
		sib(sibpt)%diag%alloc(4)=0.0	
        elseif (sib(sibpt)%diag%gdd <680.0) then
                sib(sibpt)%diag%alloc(1)=0.4-0.1 * (sib(sibpt)%diag%gdd-350.0)/330.0
		sib(sibpt)%diag%alloc(2)=0.4-0.15 * (sib(sibpt)%diag%gdd-350.0)/330.0
		sib(sibpt)%diag%alloc(3)=0.2+0.25 * (sib(sibpt)%diag%gdd-350.0)/330.0
		sib(sibpt)%diag%alloc(4)=0.0
       elseif (sib(sibpt)%diag%gdd<769.0) then
                sib(sibpt)%diag%alloc(1)=0.3+0.1 * (sib(sibpt)%diag%gdd-680.0)/89.0
		sib(sibpt)%diag%alloc(2)=0.25-0.249 * (sib(sibpt)%diag%gdd-680.0)/89.0
		sib(sibpt)%diag%alloc(3)=0.45+0.149 * (sib(sibpt)%diag%gdd-680.0)/89.0
		sib(sibpt)%diag%alloc(4)=0.0
		sib(sibpt)%diag%alloc(4)=0.0
       elseif (sib(sibpt)%diag%gdd<910.0) then
                sib(sibpt)%diag%alloc(1)=0.4   
		sib(sibpt)%diag%alloc(2)=0.4   
		sib(sibpt)%diag%alloc(3)=0.2   
		sib(sibpt)%diag%alloc(4)=0.0
       elseif (sib(sibpt)%diag%gdd<1074.0) then 
    		sib(sibpt)%diag%alloc(1)=0.4 - 0.12 * (sib(sibpt)%diag%gdd - 910.0)/164.0
		sib(sibpt)%diag%alloc(2)=0.4 !- 0.05 * (sib(sibpt)%diag%gdd - 910.0)/164.0
		sib(sibpt)%diag%alloc(3)=0.2 + 0.12 * (sib(sibpt)%diag%gdd - 910.0)/164.0
		sib(sibpt)%diag%alloc(4)=0.0
      elseif (sib(sibpt)%diag%gdd<1569.0) then
        	sib(sibpt)%diag%alloc(1)= 0.28 - 0.05 * (sib(sibpt)%diag%gdd-1074.0)/495.0
                sib(sibpt)%diag%alloc(2) = 0.4 + 0.05*(sib(sibpt)%diag%gdd-1074.0)/495.0
                sib(sibpt)%diag%alloc(3) = 0.32 
		sib(sibpt)%diag%alloc(4)=0.0
	elseif (sib(sibpt)%diag%gdd<1629.0) then
                sib(sibpt)%diag%alloc(1)=0.23
		sib(sibpt)%diag%alloc(2)=0.45-0.05 * (sib(sibpt)%diag%gdd-1569.0)/60.0	
		sib(sibpt)%diag%alloc(3)=0.32-0.12 * (sib(sibpt)%diag%gdd-1569.0)/60.0	
		sib(sibpt)%diag%alloc(4)=0.17 * (sib(sibpt)%diag%gdd-1569.0)/60.0	
	elseif (sib(sibpt)%diag%gdd<1773.0) then
        	sib(sibpt)%diag%alloc(1)=0.23-0.13 * (sib(sibpt)%diag%gdd-1629.0)/144.0
		sib(sibpt)%diag%alloc(2)=0.4-0.25 * (sib(sibpt)%diag%gdd-1629.0)/144.0
		sib(sibpt)%diag%alloc(3)=0.2-0.1 * (sib(sibpt)%diag%gdd-1629.0)/144.0
		sib(sibpt)%diag%alloc(4)=0.17+0.48 * (sib(sibpt)%diag%gdd-1629.0)/144.0
	elseif (sib(sibpt)%diag%gdd < 2184.0) then
		sib(sibpt)%diag%alloc(1)=0.1-0.05 * (sib(sibpt)%diag%gdd-1773.0)/411.0	
		sib(sibpt)%diag%alloc(2)=0.15-0.1 * (sib(sibpt)%diag%gdd-1773.0)/411.0
		sib(sibpt)%diag%alloc(3)=0.1-0.05 * (sib(sibpt)%diag%gdd-1773.0)/411.0	
		sib(sibpt)%diag%alloc(4)=0.65 + 0.2 * (sib(sibpt)%diag%gdd-1773.0)/411.0
        !EL...winterwheat root fraction reaches 0.1 towards maturity (ref: Baret et al., 1992)
	elseif (sib(sibpt)%diag%gdd < 2300.0) then
		sib(sibpt)%diag%alloc(1)=0.05-0.02 * (sib(sibpt)%diag%gdd-2184.0)/116.0
		sib(sibpt)%diag%alloc(2)=0.05-0.02 * (sib(sibpt)%diag%gdd-2184.0)/116.0
		sib(sibpt)%diag%alloc(3)=0.05-0.02 * (sib(sibpt)%diag%gdd-2184.0)/116.0
		sib(sibpt)%diag%alloc(4)=0.85+0.06 * (sib(sibpt)%diag%gdd-2184.0)/116.0
	else
               print*,'Error with gdd and alloc calculations in crop_accum.F90'
               stop
        endif

!----------------------------------
!Calculate total weight allocation	
!----------------------------------
!EL...w_main is the daily amount of dry matter added, which is the basis for 
!EL...calculation of growth and maintenance respn

!EL...Multiplication factors below were derived using the info taken 
!EL...from past studies; mainly from de Vries et al., 1989.
!EL...initial seeding density (215-222 seeds/m2 or 20-30/ft2) was taken 
!EL...from the info from ag extension services.
!EL...Average seedling weight (5 mg each) was taken using several past studies 
!EL...(e.g.Blum et al., 1980,Hameed et al., 2003)
!EL...initial carbon in seedling (0.48 g C m-2) was obtained by multiplying the 
!EL...seedling weight by 0.43.

    if(sib(sibpt)%diag%gdd < 105.0_dbl_kind) then
                sib(sibpt)%diag%w_main=0.0
    elseif (sib(sibpt)%diag%gdd == 105.0_dbl_kind) then
                sib(sibpt)%diag%w_main=0.48    
       
     !EL..Calculating the w_main by using  assim, growth resp coefficients and 
     !EL..allocation fractions from the original scheme
    elseif(sib(sibpt)%diag%gdd < 2269.0) then
       	sib(sibpt)%diag%w_main = sib(sibpt)%diag%assim_d /       &
               ((sib(sibpt)%diag%alloc(1) * 1.2214) +     &
                (sib(sibpt)%diag%alloc(2) * 1.2515) +     & 
                (sib(sibpt)%diag%alloc(3) * 1.2225) +     &
                (sib(sibpt)%diag%alloc(4) * 1.1893))
    endif

!EL...considering the fact that early wheat growth is linearly related 
!EL...to temperature(Muchow and Carberry, 1989) &
!EL...considering total daily growth, based on the growth rate 
!EL...information given in de Vries et al.1989...

!EL...basic daily growth rate for the initial growth phase, depending 
!EL...on the average no. of days between planting until the tillers come out, 
!EL...and then from Jan. 01st to spring green.

     if (sib(sibpt)%diag%gdd>105.0 .AND. sib(sibpt)%diag%gdd<=310.0) then 
	drate=0.07
        max_wmain= 4.0

    elseif (sib(sibpt)%param%biome == 22 .and. sib(sibpt)%diag%gdd>769.0 &
                  .and. sib(sibpt)%diag%gdd<1074.0) then
        drate=0.015
        max_wmain= 8.0
    endif

!EL...Temperatures and relevant fractions from the basic 
!EL...drates based on temperature 
!EL...were set based on the info in de Vries et al. 1989 

       dgrowth_opt=(max_wmain-0.48)*drate

!EL.. initial growth scheme for winter wheat

    if (sib(sibpt)%param%biome == 22) then
         if ((sib(sibpt)%diag%gdd>=105.0 .AND. sib(sibpt)%diag%gdd<=310.0) &
               .or. (sib(sibpt)%diag%gdd>769.0 .and. &
                      sib(sibpt)%diag%gdd<1074.0)) then             
            	if (sib(sibpt)%diag%tempc<=2.0) then
                	dgrowth = 0.0
            	elseif (sib(sibpt)%diag%tempc >= 2.0 .and. &
                           sib(sibpt)%diag%tempc < 10.0) then
                     	dgrowth = dgrowth_opt * sib(sibpt)%diag%rstfac_d * &
                              (0.68 * (sib(sibpt)%diag%tempc-2.)/8.)
            	elseif (sib(sibpt)%diag%tempc >= 10.0 .and. &
                           sib(sibpt)%diag%tempc < 15.0) then
                    	dgrowth = dgrowth_opt * sib(sibpt)%diag%rstfac_d * &
                              (0.68+(0.21 * (sib(sibpt)%diag%tempc-10.)/5.))
             	elseif (sib(sibpt)%diag%tempc >= 15.0 .and. &
                          sib(sibpt)%diag%tempc < 20.0) then
                     	dgrowth = dgrowth_opt * sib(sibpt)%diag%rstfac_d * &
                              (0.89+(0.11 * (sib(sibpt)%diag%tempc-15.)/5.))
             	elseif (sib(sibpt)%diag%tempc >= 20.0 .and. &
                           sib(sibpt)%diag%tempc < 25.0) then
                     	dgrowth = dgrowth_opt * sib(sibpt)%diag%rstfac_d * &
                               (1.0+(0.03 * (sib(sibpt)%diag%tempc-20.)/5.))
            	elseif (sib(sibpt)%diag%tempc>=25.0 .and. &
                           sib(sibpt)%diag%tempc<30.0) then
                     	dgrowth = dgrowth_opt * sib(sibpt)%diag%rstfac_d * &
                              (1.03+(0.02 * (sib(sibpt)%diag%tempc-25)/5.0))
            	elseif (sib(sibpt)%diag%tempc>=30.0 .and. &
                          sib(sibpt)%diag%tempc<35.0) then
                      	dgrowth = dgrowth_opt * sib(sibpt)%diag%rstfac_d * &
                             (1.05+(0.01 * (sib(sibpt)%diag%tempc-30)/5.0))
	       endif

               if(sib(sibpt)%diag%gdd == 105.0_dbl_kind .or. &
                  sib(sibpt)%diag%gdd==769.0_dbl_kind) then
                      sib(sibpt)%diag%w_main_pot=0.48
                      sib(sibpt)%diag%w_main=0.48
               endif

              sib(sibpt)%diag%w_main_pot=sib(sibpt)%diag%w_main_pot+dgrowth
              sib(sibpt)%diag%w_main=max(sib(sibpt)%diag%w_main_pot, &
                                                               sib(sibpt)%diag%w_main) 
              sib(sibpt)%diag%w_main=min(sib(sibpt)%diag%w_main,max_wmain)
          endif  !gdd >= 105...

    endif  !biome == 22

!EL.. initial growth scheme for spring wheat
    if (sib(sibpt)%param%biome == 23) then

        if (sib(sibpt)%diag%gdd>=105.0 .AND. sib(sibpt)%diag%gdd<=310.0) then             
  	    if (sib(sibpt)%diag%tempc<-10.0) then
        	   dgrowth = 0.0
	    elseif (sib(sibpt)%diag%tempc < 0.0) then
          	  dgrowth = dgrowth_opt * sib(sibpt)%diag%rstfac_d * 0.01
	    elseif (sib(sibpt)%diag%tempc < 20.0) then
          	  dgrowth = dgrowth_opt * sib(sibpt)%diag%rstfac_d * &
                       (0.01 + (0.89 * sib(sibpt)%diag%tempc / 20.0))
	    elseif (sib(sibpt)%diag%tempc < 25.0) then
                  dgrowth = dgrowth_opt * sib(sibpt)%diag%rstfac_d * &
                       (0.9 + (0.1 * (sib(sibpt)%diag%tempc - 20.0) / 5.0))
	    elseif (sib(sibpt)%diag%tempc < 35.0) then
                  dgrowth = dgrowth_opt * sib(sibpt)%diag%rstfac_d * &
                       (1.0 + (0.2 * (sib(sibpt)%diag%tempc - 25.0) / 10.0 ))
            endif

            if(time%doy == sib(sibpt)%diag%emerg_d) then
                   sib(sibpt)%diag%w_main_pot=0.48 + dgrowth
                   sib(sibpt)%diag%w_main=sib(sibpt)%diag%w_main_pot
           endif

           sib(sibpt)%diag%w_main_pot=sib(sibpt)%diag%w_main_pot+dgrowth
           sib(sibpt)%diag%w_main=max(sib(sibpt)%diag%w_main_pot,  &
                    sib(sibpt)%diag%w_main) 
           sib(sibpt)%diag%w_main=min(sib(sibpt)%diag%w_main,max_wmain)
      
        endif  !gdd >= 105

   endif  !biome == 23

         !EL..back calculation of the new assim_d:
          assimd_new = sib(sibpt)%diag%w_main *       &
             ((sib(sibpt)%diag%alloc(1) * 1.2214) +     &
             ( sib(sibpt)%diag%alloc(2) * 1.2515) +     & 
             ( sib(sibpt)%diag%alloc(3) * 1.2225) +     &
             ( sib(sibpt)%diag%alloc(4) * 1.1893))

          sib(sibpt)%diag%assim_d = assimd_new
 
do i = 1,4
      sib(sibpt)%diag%allocwt(i) = sib(sibpt)%diag%w_main * sib(sibpt)%diag%alloc(i) 
enddo

!--------------------------
!Calculate maintanence respiration
!--------------------------

!EL...fractions of NSC and proteins were decided based on 
!EL...Collar and Aksland, 2001, and Blum, 1998.
!EL...maint. coeff. info from Norman and Arkebauer, 1991)

       temp1 = (0.03 * coeff) *                   &
             (2.0**((sib(sibpt)%diag%tempc - 20.0) / 10.0))
       sib(sibpt)%diag%phen_maintr(1) = sib(sibpt)%diag%cum_wt(1)       &
             * 0.2 * temp1  

!EL...Maintenance respn coefficients were based on 
!EL...de Vries et al., 1989 and Wang et al., 1992)

	temp1 = (0.016 * coeff) *                   &
             (2.0**((sib(sibpt)%diag%tempc - 20.0) / 10.0))
       sib(sibpt)%diag%phen_maintr(2) = sib(sibpt)%diag%cum_wt(2) *     &
              0.25 * temp1

       temp1 = 0.01 * coeff * (2.0**((sib(sibpt)%diag%tempc-20.0) / 10.0))
       sib(sibpt)%diag%phen_maintr(3) = sib(sibpt)%diag%cum_wt(3)          &
              * 0.3 * temp1

       temp1 = 0.015 * coeff *      &
               (2.0 ** ((sib(sibpt)%diag%tempc - 20.0) / 10.0))
       sib(sibpt)%diag%phen_maintr(4) = sib(sibpt)%diag%cum_wt(4) * 0.5 * temp1

!---------------------------------------------    
!Calculate cumulative drywt.in each plant part 
!---------------------------------------------
!EL..cumulative dry weight will be used in calculation of maintenance respn

	 do j=1,4	 
             sib(sibpt)%diag%cum_wt(j) = sib(sibpt)%diag%cum_wt(j) + &
                    sib(sibpt)%diag%allocwt(j)
         enddo
       
	  if (sib(sibpt)%diag%gdd>2290.0) then
              do i = 1,4
                   sib(sibpt)%diag%cum_wt(i) = 0.0001
              enddo
 	  endif

!----------------------------
!Calculate growth respiration
!----------------------------
 
        sib(sibpt)%diag%phen_growthr(1) = sib(sibpt)%diag%allocwt(1)*0.406 * coeff
        sib(sibpt)%diag%phen_growthr(2) = sib(sibpt)%diag%allocwt(2)*0.461 * coeff
        sib(sibpt)%diag%phen_growthr(3) = sib(sibpt)%diag%allocwt(3)*0.408 * coeff
        sib(sibpt)%diag%phen_growthr(4) = sib(sibpt)%diag%allocwt(4)*0.347 * coeff

!------------------------------
!Calculate dry weight change
!-----------------------------

      do i = 1,4
        sib(sibpt)%diag%wch(i) = sib(sibpt)%diag%allocwt(i) - &
                sib(sibpt)%diag%phen_maintr(i)
      enddo
	
!--------------------------------------------------------------
!Calculate final cumulative dry weight (g C m-2) of each plant part
!--------------------------------------------------------------
	  do j=1,4	 
                 sib(sibpt)%diag%cum_drywt(j) = sib(sibpt)%diag%cum_drywt(j) + &
                      sib(sibpt)%diag%wch(j)
          enddo

         sib(sibpt)%diag%tot_biomass= sib(sibpt)%diag%cum_drywt(2) + &
                  sib(sibpt)%diag%cum_drywt(3)+ sib(sibpt)%diag%cum_drywt(4)

!---------------------------------------------------
!
!Decreasing the vmax during the seed filling stage and senescence down
!  to 60% of original vmax value.  
!This is supported in Crafts-Brandner et al. (Plant Physiol., 1992), 
!   Jiang et al. (Plant Physiol., 1993), and Jiang et al. (Photosynthesis Research, 1997)
!
!kdcorbin, 02/11
!---------------------------------------------------

if (sib(sibpt)%diag%gdd >= vmax_start(wwheat_num) .and. &
     sib(sibpt)%diag%gdd < vmax_stop(wwheat_num)) then
          vmax_factor = crop_vmax0a(wwheat_num) - crop_vmax0b(wwheat_num)
          sib(sibpt)%param%vmax0(1) = crop_vmax0a(wwheat_num) - &
                  vmax_factor * (sib(sibpt)%diag%gdd - vmax_start(wwheat_num)) / &
                    (vmax_stop(wwheat_num)-vmax_start(wwheat_num)-1)
endif

!------------------------------------------------
!final leaf weight (C g m-2) 
!-------------------------------------------------
!EL...i.e. adjustment for leaf senescence and harvest event
 
      if( sib(sibpt)%diag%gdd <= 0.0 ) then
           sib(sibpt)%diag%leafwt_c =  0.01 
      elseif (sib(sibpt)%diag%gdd <= 2269.0) then
            sib(sibpt)%diag%leafwt_c = sib(sibpt)%diag%cum_drywt(2)
      elseif (sib(sibpt)%diag%gdd < 2440.0) then
            sib(sibpt)%diag%leafwt_c = sib(sibpt)%diag%cum_drywt(2) -   &
                     0.99 * sib(sibpt)%diag%cum_drywt(2) *      &
                    ((sib(sibpt)%diag%gdd - 2269.0) / 171.0)
      endif

!--------------
!Calculate LAI
!-------------
!EL...convert to dry weight g m-2 and then multiply by SLA; 	
!EL...SLA determined by the averages based on several studies

      sib(sibpt)%diag%phen_LAI=sib(sibpt)%diag%leafwt_c*2*0.02 

!itb_crop...at the moment that growing degree days (gdd) passes
!itb_crop...105, we will initialize the LAI

    if (sib(sibpt)%diag%gdd >= 105.0_dbl_kind .AND. &
        time%doy >= sib(sibpt)%diag%emerg_d  .AND.  &
        !kdcorbin, 02/11 - only set switch once
        sib(sibpt)%diag%phen_switch == 0) then
        sib(sibpt)%diag%phen_switch = 1
        call set_ti(sib(sibpt))
    endif

end subroutine wheat_phen

end subroutine crop_accum

