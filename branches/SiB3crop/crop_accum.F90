!==================SUBROUTINE CROP_ACCUM=======================================
subroutine crop_accum(sib,time)

use kinds
use sibtype
use timetype
use physical_parameters, only: tice    

implicit none

integer(kind=int_kind) :: i0,n,pd,i,j,ndf60,ndf65
real(kind=dbl_kind)    :: temp_accum,assim_accum,allocwt_accum


!----------------------------------------------------------------------
type(sib_t), intent(inout) :: sib
type(time_struct), intent(in) :: time
!----------------------------------------------------------------------  

   temp_accum = 0.0_dbl_kind
   	

   do i0 = 1, sib%diag%tb_indx
    
      temp_accum = temp_accum + sib%diag%tb_temp(i0) 
!print*,i0,time%doy,sib%diag%tb_temp(i0)
!pause
!      print*,i0,sib%diag%tb_temp(i0),temp_accum
 

   enddo

   sib%diag%ta_bar = temp_accum / float(sib%diag%tb_indx)

!conversion of avg. daily temperature (ta_bar) from Kelvin to Fahrenheit

sib%diag%tempf=((sib%diag%ta_bar-273.15)*1.8)+32.0	

!conversion of avg. daily temperature (ta_bar) from Kelvin to Celcius

sib%diag%tempc=sib%diag%ta_bar - tice  !tice=273K

!Calling for different phenology schemes based on the year and the crop
	if(mod(time%year,2)==0) then  
		call soy_phen
	else
		call corn_phen	
	endif
!print*,time%year,sib%diag%tempf,time%doy

contains

!------------------------------------------------------------------------------
subroutine corn_phen
!------------------------------------------------------------------------------

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

!print*,pd


!----------------------------
!Calculate growing degree days
!-----------------------------


	if (pd>0 							.and.		& !this line was added to avoid gdd calculation before real planting date, since pd is printed out as 0 before the real planting date based on the above ndf60==5 criterion
		time%doy >= pd            .and.          &
        sib%diag%tempf>50.0      .and.          &
        sib%diag%tempf<86.0)      then

    	sib%diag%gdd=sib%diag%gdd + sib%diag%tempf- 50.0_dbl_kind
	
	endif
	if (time%doy>300) then
        sib%diag%gdd=0.0
	endif
!   if(sib%diag%ta_bar > 20.0_dbl_kind + tice) then

!     sib%diag%gdd = sib%diag%gdd + sib%diag%ta_bar      &
!                                 - 20.0_dbl_kind + tice
!   endif


!------------------------------
!	Reading and summing assimn
!-------------------------------
assim_accum=0.0_dbl_kind
 do i0 = 1, sib%diag%tb_indx
    
      assim_accum = assim_accum + (sib%diag%tb_assim(i0)*time%dtsib) !multiplied by the no. secs per each timestep (i.e. tbsib) to convert assim mol sec-1 to mol

   enddo

   sib%diag%assim_d = assim_accum*12 !multiplied by 12 to convert mol to g


!print*,sib%diag%gdd,assim_d

!-----------------------------------------------------------
! allocation sheme for fractions for assimilate partitioning
!-----------------------------------------------------------   
	!1-roots, 2-leaves,3-stems,4-products (flowers and grains)

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


	elseif(sib%diag%gdd<100 .or. sib%diag%gdd>=2730.0)then
		sib%diag%alloc(1)=0.0
		sib%diag%alloc(2)=0.0	
		sib%diag%alloc(3)=0.0	
		sib%diag%alloc(4)=0.0
        
    endif
!----------------------------------
!Calculate total weight allocation	
!----------------------------------
   !w_main is the (drywt+maintenancerespiration) together
    
    if((sib%diag%gdd>100).and.((sib%diag%gdd<2730))) then
       	sib%diag%w_main=sib%diag%assim_d/((sib%diag%alloc(1)*1.2214)+(sib%diag%alloc(2)*1.2515)& 
            +(sib%diag%alloc(3)*1.2225)+(sib%diag%alloc(4)*1.2095))! the coefficients given corresponds to growth resp & (drywt+maint R) ref Penning deVries (1989); here w_main is the (drywt+maint.respn)
    else
        sib%diag%w_main =0.0
    endif
!print *,sib%diag%gdd,sib%diag%assim_d,sib%diag%w_main
!pause
 

!Calculate w_main allocation to different plant parts
  

         sib%diag%allocwt(1)=sib%diag%w_main*sib%diag%alloc(1) ! clculates absolute allocation for roots, using allocation fraction for roots and assimilation
        sib%diag%allocwt(2)=sib%diag%w_main*sib%diag%alloc(2) ! clculates absolute allocation for leaves, using allocation fraction for leaves and assimilation
        sib%diag%allocwt(3)=sib%diag%w_main*sib%diag%alloc(3) ! clculates absolute allocation for stem, using allocation fraction for stems and assimilation
        sib%diag%allocwt(4)=sib%diag%w_main*sib%diag%alloc(4) ! clculates absolute allocation for flowers and grains using allocation fraction for those and assimilation

!---------------------------------------------    
!Calculate cumulative drywt.in each plant part (before respn)
!---------------------------------------------

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

!print*,sib%diag%assim_d,sib%diag%allocwt(2),sib%diag%cum_w(time%doy,2)

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
   
       sib%diag%phen_maintr(1)=sib%diag%cum_wt(time%doy-1,1)*0.18*(0.03*2*12/44)*(2.0**((sib%diag%tempc-20.)/10.))   !Q10 coefficient is 1.8 for soybean and 2.0 for corn&
                          !(Norman and Arkebauer, 1991);0.18 is the nonstructural C fraction of root C needing maintenance(calculations based on Brouquisse et al., 1998)&
	                      !maint. coeff. info from Penning de Vries,1989, Amthor, 1984, and Norman and Arkebauer, 1991)
!print *,time%doy,(2.0**((sib%diag%tempc(time%doy)-20.)/10.))
!pause
       sib%diag%phen_maintr(2)=sib%diag%cum_wt(time%doy-1,2)*0.27*0.03*2*0.75*12/44*(2.0**((sib%diag%tempc-20.)/10.))
 !multiplied by 0.75 to incorporate that maintenance respiration during the daytime is half that of  nighttime.. (Penning de Vries, 1989)
                           ! 0.27 is the nonstructural fraction of leaf C needing maintenance (calculations based on Brouquisse et al., 1998)
       sib%diag%phen_maintr(3)=sib%diag%cum_wt(time%doy-1,3)*0.24*0.01*2*12/44*(2.0**((sib%diag%tempc-20.)/10.))! 0.24 is the nonstructural fraction of stem C  needing maintenance(calculations based on Brouquisse et al., 1998).
       sib%diag%phen_maintr(4)=sib%diag%cum_wt(time%doy-1,4)*0.4*0.015*2*12/44*(2.0**((sib%diag%tempc-20.)/10.))
! 0.4 is the nonstructural fraction of seed C  needing maintenance (calculations based on Beauchemin et al., 1997) 


!print*,sib%diag%cum_wt(time%doy,2),sib%diag%phen_maintr(2)
		
	if (time%doy>300) then	
		sib%diag%phen_maintr(1)=0.0
		sib%diag%phen_maintr(2)=0.0
		sib%diag%phen_maintr(3)=0.0
		sib%diag%phen_maintr(4)=0.0
	endif
		
!------------------------------
!Calculate dry weight change
!-----------------------------

	  sib%diag%wch(1)=sib%diag%allocwt(1)- sib%diag%phen_maintr(1)
      sib%diag%wch(2)=sib%diag%allocwt(2)- sib%diag%phen_maintr(2)
      sib%diag%wch(3)=sib%diag%allocwt(3)- sib%diag%phen_maintr(3)
      sib%diag%wch(4)=sib%diag%allocwt(4)- sib%diag%phen_maintr(4)
	
	  if (time%doy>300) then
		sib%diag%wch(1)=0.0
		sib%diag%wch(2)=0.0
		sib%diag%wch(3)=0.0
		sib%diag%wch(4)=0.0
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
	  
      sib%diag%final_drywt(1)=sib%diag%cum_drywt(time%doy,1) !Renamed to be output on a daily basis in a separate text file
	  sib%diag%final_drywt(2)=sib%diag%cum_drywt(time%doy,2)
	  sib%diag%final_drywt(3)=sib%diag%cum_drywt(time%doy,3)
	  sib%diag%final_drywt(4)=sib%diag%cum_drywt(time%doy,4)

!------------------------------------------------
!final leaf weight (C g m-2) (after adjustment for senescence and harvest event)
!-------------------------------------------------
       
		if (sib%diag%gdd<2730) then
	       sib%diag%leafwt_c=sib%diag%cum_drywt(time%doy,2)
	
       elseif (sib%diag%gdd>=2730.0 .and. sib%diag%gdd<3300.0) then
          sib%diag%leafwt_c=0.95*sib%diag%cum_drywt(time%doy,2)-(0.95-0.2)*sib%diag%cum_drywt(time%doy,2)*((sib%diag%gdd-2730.0)/570)
        
        else
          sib%diag%leafwt_c=sib%diag%cum_drywt(time%doy,2)*0.01

        endif

!--------------
!Calculate LAI
!-------------

      sib%diag%phen_LAI=sib%diag%leafwt_c*2*0.02 	!(convert to dry weight g m-2 and then multiplied by SLA; 
!SLA determined by the averages based on several studies)



 	sib%diag%tb_indx = 0	 !at the end of each day tb_index is set to zero


!print*,sib%diag%assim_d,sib%diag%phen_LAI



!	open(unit=20,file='phen_corn_test.dat',form='formatted')
 
		write(20,'(i4.4,2x,i3.3,2x,43(1x,f11.2))')time%year,time%doy,sib%diag%tempf,sib%diag%tempc,sib%diag%gdd,sib%diag%assim_d,sib%diag%alloc(1:4),&
             sib%diag%w_main,sib%diag%allocwt(1:4),sib%diag%cum_w(1:4),sib%diag%phen_growthr(1:4),sib%diag%phen_maintr(1:4),sib%diag%wch(1:4),&
sib%diag%final_drywt(1:4),sib%diag%leafwt_c,sib%diag%phen_LAI 

end subroutine corn_phen


!-----------------------------------------------------------------------------------------------------------
subroutine soy_phen
!-----------------------------------------------------------------------------------------------------------

if (sib%diag%tempf<65.0) then
    ndf65=0			!ndf65= no. of days withe avg. temperature above 65F
elseif (sib%diag%tempf>=65.0) then
    ndf65=ndf65+1
endif

if (ndf65==5) then

    pd=time%doy

endif

!print*,pd


!----------------------------
!Calculate growing degree days
!-----------------------------


	if (pd>0 							.and.		& !this line was added to avoid gdd calculation before real planting date, since pd is printed out as 0 before the real planting date based on the above ndf65==5 criterion
		time%doy >= pd            .and.          &
        sib%diag%tempf>50.0      .and.          &
        sib%diag%tempf<86.0)      then

    	sib%diag%gdd=sib%diag%gdd + sib%diag%tempf- 50.0_dbl_kind
	
	endif

	if (time%doy>300) then
        sib%diag%gdd=0.0
	endif

!   if(sib%diag%ta_bar > 20.0_dbl_kind + tice) then

!     sib%diag%gdd = sib%diag%gdd + sib%diag%ta_bar      &
!                                 - 20.0_dbl_kind + tice
!   endif


!------------------------------
!	Reading and summing assimn
!-------------------------------
assim_accum=0.0_dbl_kind
 do i0 = 1, sib%diag%tb_indx
    
      assim_accum = assim_accum + (sib%diag%tb_assim(i0)*time%dtsib) !multiplied by the no. secs per each timestep (i.e. tbsib) to convert assim mol sec-1 to mol

   enddo

   sib%diag%assim_d = assim_accum*12 !multiplied by 12 to convert mol to g


!print*,sib%diag%gdd,assim_d

!-----------------------------------------------------------
! allocation sheme for fractions for assimilate partitioning
!-----------------------------------------------------------   
	!1-roots, 2-leaves,3-stems,4-products (flowers and grains)
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
   !w_main is the (drywt+maintenancerespiration) together
    
     if((time%doy>pd+10).and.(time%doy<(pd+130))) then
       	sib%diag%w_main=sib%diag%assim_d/((sib%diag%alloc(1)*1.2929)+(sib%diag%alloc(2)*1.431)& 
            +(sib%diag%alloc(3)*1.2946)+(sib%diag%alloc(4)*1.6752))! the coefficients given corresponds to growth resp & (drywt+maint R) ref Penning deVries (1989); here w_main is the (drywt+maint.respn)
    else
        sib%diag%w_main =0.0
    endif
!print *,sib%diag%gdd,sib%diag%assim_d,sib%diag%w_main
!pause
 

!Calculate w_main allocation to different plant parts
  

         sib%diag%allocwt(1)=sib%diag%w_main*sib%diag%alloc(1) ! clculates absolute allocation for roots, using allocation fraction for roots and assimilation
        sib%diag%allocwt(2)=sib%diag%w_main*sib%diag%alloc(2) ! clculates absolute allocation for leaves, using allocation fraction for leaves and assimilation
        sib%diag%allocwt(3)=sib%diag%w_main*sib%diag%alloc(3) ! clculates absolute allocation for stem, using allocation fraction for stems and assimilation
        sib%diag%allocwt(4)=sib%diag%w_main*sib%diag%alloc(4) ! clculates absolute allocation for flowers and grains using allocation fraction for those and assimilation

!---------------------------------------------    
!Calculate cumulative drywt.in each plant part (before respn)
!---------------------------------------------
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

!	  do j=1,4	 
       
 !    	sib%diag%cum_w(time%doy,j)=sib%diag%cum_w(time%doy-1,j)+sib%diag%allocwt(j)
       
!	  If ((time%doy-1)==0) then
     
 !       sib%diag%cum_w(time%doy,j)=0
     
!	  endif
 
 !     enddo

!print*,sib%diag%assim_d,sib%diag%allocwt(2),sib%diag%cum_w(time%doy,2)

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
   
       sib%diag%phen_maintr(1)=sib%diag%cum_wt(time%doy-1,1)*0.32*(0.03*2*12/44)*(2.0**((sib%diag%tempc-20.)/10.))   !Q10 coefficient is 1.8 for soybean and 2.0 for corn&
                          !(Norman and Arkebauer, 1991);0.32 is the nonstructural C fraction of root C needing maintenance(calculations based on Allen et al., 1998 and Rogers et al., 2006)&
	                      !maint. coeff. info from Penning de Vries,1989, Amthor, 1984, and Norman and Arkebauer, 1991)
!print *,time%doy,(2.0**((sib%diag%tempc(time%doy)-20.)/10.))
!pause
       sib%diag%phen_maintr(2)=sib%diag%cum_wt(time%doy-1,2)*0.38*0.03*2*0.75*12/44*(2.0**((sib%diag%tempc-20.)/10.))
 !multiplied by 0.75 to incorporate that maintenance respiration during the daytime is half that of  nighttime.. (Penning de Vries, 1989)
                           ! 0.38 is the nonstructural fraction of leaf C needing maintenance (calculations based onAllen et al., 1998 and Rogers et al., 2006)
       sib%diag%phen_maintr(3)=sib%diag%cum_wt(time%doy-1,3)*0.32*0.01*2*12/44*(2.0**((sib%diag%tempc-20.)/10.))! 0.32 is the nonstructural fraction of stem C  needing maintenance(calculations based on Allen et al., 1998 and Rogers et al., 2006).
       sib%diag%phen_maintr(4)=sib%diag%cum_wt(time%doy-1,4)*0.46*0.015*2*12/44*(2.0**((sib%diag%tempc-20.)/10.))
! 0.4 is the nonstructural fraction of seed C  needing maintenance (calculations based on Allen et al., 1998 and Rogers et al., 2006) 

	
	  if (time%doy>300) then
		sib%diag%phen_maintr(1)=0.0
		sib%diag%phen_maintr(2)=0.0
		sib%diag%phen_maintr(3)=0.0
		sib%diag%phen_maintr(4)=0.0
	  endif
!print*,sib%diag%cum_w(time%doy,2),sib%diag%phen_maintr(2)

!------------------------------
!Calculate dry weight change
!-----------------------------

	  sib%diag%wch(1)=sib%diag%allocwt(1)- sib%diag%phen_maintr(1)
      sib%diag%wch(2)=sib%diag%allocwt(2)- sib%diag%phen_maintr(2)
      sib%diag%wch(3)=sib%diag%allocwt(3)- sib%diag%phen_maintr(3)
      sib%diag%wch(4)=sib%diag%allocwt(4)- sib%diag%phen_maintr(4)

	 if (time%doy>300) then
		sib%diag%wch(1)=0.0
		sib%diag%wch(2)=0.0
		sib%diag%wch(3)=0.0
		sib%diag%wch(4)=0.0
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

	  sib%diag%final_drywt(1)=sib%diag%cum_drywt(time%doy,1) !Renamed to be output on a daily basis in a separate text file
	  sib%diag%final_drywt(2)=sib%diag%cum_drywt(time%doy,2)
	  sib%diag%final_drywt(3)=sib%diag%cum_drywt(time%doy,3)
	  sib%diag%final_drywt(4)=sib%diag%cum_drywt(time%doy,4)

!------------------------------------------------
!final leaf weight (C g m-2) (after adjustment for senescence and harvest event)
!-------------------------------------------------
       
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

      sib%diag%phen_LAI=sib%diag%leafwt_c*2*0.02 	!(convert to dry weight g m-2 and then multiplied by SLA; 
!SLA determined by the averages based on several studies)



 	sib%diag%tb_indx = 0	 !at the end of each day tb_index is set to zero


print*,pd,sib%diag%assim_d,sib%diag%phen_LAI



	open(unit=20,file='phen_soy_test.dat',form='formatted')
 
		write(20,'(i4.4,2x,i3.3,2x,43(1x,f11.2))')time%year,time%doy,sib%diag%tempf,sib%diag%tempc,sib%diag%gdd,sib%diag%assim_d,sib%diag%alloc(1:4),&
             sib%diag%w_main,sib%diag%allocwt(1:4),sib%diag%cum_w(1:4),sib%diag%phen_growthr(1:4),sib%diag%phen_maintr(1:4),sib%diag%wch(1:4),&
sib%diag%final_drywt(1:4),sib%diag%leafwt_c,sib%diag%phen_LAI 
	


sib%diag%tb_indx = 0	 !at the end of each day tb_index is set to zero
end subroutine soy_phen


end subroutine crop_accum
