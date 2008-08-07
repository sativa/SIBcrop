!---------------------------------------------------------------------
subroutine init_var( sib)
!---------------------------------------------------------------------

!sets all variable, both global and local, initially to 0;
!called from SiBDRV.f90 after sib has been allocated.

use sibtype
use timetype
use sib_const_module
use sib_io_module
implicit none

! parameters
type(sib_t), dimension(subcount), intent(inout) :: sib
integer i
integer j
integer k

!set all diagnostic variables to 0.0;

!itb_crop
    sib%diag%tb_indx = 0
    sib%diag%pd      = 0
    sib%diag%ndf_opt = 0
    sib%diag%pd7     = 0
    sib%diag%pd7_est = 0
    sib%diag%pdindx7 = 0

    do j=1,20000
      sib%diag%tb_temp(j) = 0.0
      sib%diag%tb_assim(j) = 0.0 !added for daily accumulation of assimilation-EL
!El...the following added for testing..
     sib%diag%tb_rst(j) =0.0
    enddo
	
		sib%diag%tempf = 0.0
		sib%diag%tempc=0.0
		sib%diag%assim_d = 0.0
		sib%diag%w_main=0.0
!	        sib%diag%w_main_pot=0.0
!                sib%diag%cropht=0.0
                sib%diag%tot_biomass=0.0
!                sib%diag%tot_prod_an=0.0
!                sib%diag%tot_BM_an= 0.0
!                sib%diag%prodwt=0.0
		sib%diag%leafwt_c = 0.0
		sib%diag%phen_LAI=0.0
!     do i=1,13
		
!                sib%diag%tot_BM_an= 0.0
 !               sib%diag%tot_prod_an=0.0
  !   enddo


!	 do j=1,4

!			sib%diag%allocwt(i,j) = 0.0
!			sib%diag%cum_wt(i,j) = 0.0
!			sib%diag%cum_drywt(i,j) = 0.0
!	 enddo
!    enddo
      
     do j=1,4		


			sib%diag%cum_wt(j) = 0.0
                        sib%diag%cum_wt_P(j)=0.0
			sib%diag%cum_drywt(j) = 0.0

			sib%diag%alloc(j) = 0.0
			sib%diag%phen_maintr(j) = 0.0
			sib%diag%phen_growthr(j) = 0.0
			sib%diag%wch(j) = 0.0
			sib%diag%allocwt(j) = 0.0
    enddo


    sib%diag%ta_bar = 0.0
!EL..added for testing..
    sib%diag%rstfac_d = 0.0
    sib%diag%gdd    = 0.0
    sib%diag%gdd_c  = 0.0 !added temporarily

    sib%diag%use_phen = 0
	sib%diag%phen_switch=0


!itb_crop_end


	sib%diag%eastar = 0.0
	sib%diag%rha = 0.0
	sib%diag%psy = 0.0
	sib%diag%cas_cap_heat = 0.0
	sib%diag%cas_cap_vap = 0.0
	sib%diag%cas_cap_co2 = 0.0
	sib%diag%cas_e_storage = 0.0
	sib%diag%canex  = 0.0
	sib%diag%wc = 0.0
	sib%diag%wg  = 0.0
	
	sib%diag%areas  = 0.0
	sib%diag%a_areas = 0.0
	sib%diag%tsnow = 0.0
	sib%diag%snowmelt = 0.0
	sib%diag%www_tot_soil = 0.0
	sib%diag%roff = 0.0
	sib%diag%roffo  = 0.0
	sib%diag%qqq  = 0.0
	sib%diag%hr  = 0.0
	sib%diag%hrr  = 0.0
	sib%diag%respg = 0.0
	sib%diag%www_inflow = 0.0
	sib%diag%cu = 0.0
	sib%diag%ct = 0.0
	sib%diag%ustar = 0.0
	sib%diag%ventmf  = 0.0
	sib%diag%thvgm = 0.0
	sib%diag%ecmass  = 0.0
	sib%diag%egmass = 0.0
	sib%diag%chf = 0.0
	sib%diag%shf = 0.0
	sib%diag%ra  = 0.0
	sib%diag%rb = 0.0
	sib%diag%rc     = 0.0
	sib%diag%rd = 0.0
	sib%diag%rsoil  = 0.0
	sib%diag%rds  = 0.0
	sib%diag%thermk  = 0.0
	sib%diag%tgeff = 0.0
	sib%diag%thgeff = 0.0
	sib%diag%shgeff  = 0.0
	sib%diag%p0  = 0.0
	sib%diag%pcpg_rain = 0.0
	sib%diag%pcpg_snow = 0.0
	sib%diag%cuprt = 0.0
	sib%diag%lsprt = 0.0
	sib%diag%hg  = 0.0
	sib%diag%hc   = 0.0
	sib%diag%hs  = 0.0
	sib%diag%fss = 0.0
	sib%diag%fws = 0.0
	sib%diag%ec  = 0.0
	sib%diag%eg = 0.0
	sib%diag%es  = 0.0
	sib%diag%egi   = 0.0
	sib%diag%eci  = 0.0
	sib%diag%egs  = 0.0
	sib%diag%ess = 0.0
	sib%diag%ect  = 0.0
	sib%diag%aparkk  = 0.0
	sib%diag%pfd     = 0.0
	sib%diag%cflux = 0.0
	sib%diag%flux13c = 0.0
	sib%diag%flux12c  = 0.0
	sib%diag%flux_turb  = 0.0


	do i=1,13
		sib%diag%tot_an(i) = 0.0
		do j=1,nsoil
			sib%diag%tot_ss(i,j) = 0.0
		enddo
	enddo

	do j=1,nsoil
		
		sib%diag%soilscale(j) = 0.0
		sib%diag%soilq10(j) = 0.0
	enddo

	do i=1,6
		sib%diag%assimnp(i) = 0.0
		sib%diag%antemp(i) = 0.0
		sib%diag%ansqr(i) = 0.0
		sib%diag%omepot(i) = 0.0
		sib%diag%assimpot(i) = 0.0
		sib%diag%assimci(i) = 0.0
		sib%diag%wsfws(i) = 0.0
		sib%diag%wsfht(i) = 0.0
		sib%diag%wsflt(i) = 0.0
		sib%diag%wci(i) = 0.0
		sib%diag%whs(i) = 0.0
		sib%diag%wags(i) = 0.0
		sib%diag%wegs(i) = 0.0
		sib%diag%kiecps(i) = 0.0
		sib%diag%d13cassimn(i) = 0.0
		sib%diag%c13assimn(i) = 0.0
		sib%diag%c12assimn(i) = 0.0
		sib%diag%rcassimn(i) = 0.0
		sib%diag%ggl(i) = 0.0
		sib%diag%pco2i(i) = 0.0
		sib%diag%pco2c(i) = 0.0
		sib%diag%pco2s(i) = 0.0
		sib%diag%respc(i) = 0.0
		sib%diag%assim(i) = 0.0
		sib%diag%assimn(i) = 0.0
	enddo

	do i=1,4
		sib%diag%rstfac(i) = 0.0
	enddo

	do i=1,3
		sib%diag%snow_end(i) = 0.0
		sib%diag%radt(i) = 0.0
	enddo

	do i=1,2
		sib%diag%drag(i) = 0.0
		sib%diag%radc3(i) = 0.0
		do j=1,2
			sib%diag%salb(i,j) = 0.0
			do k=1,2
				sib%diag%radfac(i,j,k) = 0.0
			enddo
		enddo
	enddo



	!!  sib%diag%eff_poros(-nsnow+1:nsoil)

!------------------------------------------------------------------
!                   BOUNDARY CONDITION VARIABLES
!------------------------------------------------------------------

	sib%param%biome = 0.0
	sib%param%chil = 0.0
	sib%param%phc  = 0.0
	sib%param%z1  = 0.0
!	sib%param%z2 = 0.0
	sib%param%poros = 0.0
	sib%param%satco  = 0.0
	sib%param%bee  = 0.0
	sib%param%phsat = 0.0
	sib%param%slope = 0.0
	sib%param%vcover  = 0.0
	sib%param%zm = 0.0
	sib%param%wopt  = 0.0
	sib%param%wsat  = 0.0
	sib%param%sandfrac = 0.0
	sib%param%clayfrac = 0.0
	sib%param%vwcmin = 0.0
	sib%param%czc = 0.0
	sib%param%fieldcap  = 0.0
	sib%param%aparc = 0.0
	sib%param%aparc1 = 0.0
	sib%param%aparc2 = 0.0
	sib%param%zlt = 0.0
	sib%param%zlt1 = 0.0
	sib%param%zlt2  = 0.0
	sib%param%green  = 0.0
	sib%param%green1 = 0.0
	sib%param%green2 = 0.0
	sib%param%z0d  = 0.0
	sib%param%z0d1 = 0.0
	sib%param%z0d2  = 0.0
	sib%param%z0 = 0.0
	sib%param%z01 = 0.0
	sib%param%z02  = 0.0
	sib%param%zp_disp = 0.0
	sib%param%zp_disp1  = 0.0
	sib%param%zp_disp2 = 0.0
	sib%param%zpd_adj1 = 0.0
	sib%param%zpd_adj  = 0.0
	sib%param%zpd_adj2 = 0.0
	sib%param%cc1 = 0.0
	sib%param%cc11 = 0.0
	sib%param%cc12  = 0.0
	sib%param%cc2 = 0.0
	sib%param%cc21 = 0.0
	sib%param%cc22 = 0.0
	sib%param%rbc = 0.0
	sib%param%rbc1 = 0.0
	sib%param%rbc2 = 0.0
	sib%param%rdc = 0.0
	sib%param%rdc1 = 0.0
	sib%param%rdc2 = 0.0
	sib%param%gmudmu = 0.0
	sib%param%gmudmu1 = 0.0
	sib%param%gmudmu2 = 0.0
	sib%param%d13cresp = 0.0
	sib%param%d13cresp1 = 0.0
	sib%param%d13cresp2 = 0.0
	sib%param%d13cresp2 = 0.0








	do i=1,2
		sib%param%satcap(i) = 0.0
		sib%param%soref(i) = 0.0
		do j=1,2
			sib%param%tran(i,j) = 0.0
			sib%param%ref(i,j)  = 0.0
		enddo
	enddo
	
	do i=1,5
		sib%param%vmax0(i)  = 0.0
		sib%param%trop(i) = 0.0
		sib%param%trda(i) = 0.0
		sib%param%trdm(i) = 0.0
		sib%param%respcp(i) = 0.0
		sib%param%slti(i) = 0.0
		sib%param%shti(i) = 0.0
		sib%param%hltii(i) = 0.0
		sib%param%hhti(i)  = 0.0
		sib%param%effcon(i) = 0.0
		sib%param%binter(i) = 0.0
		sib%param%gradm(i) = 0.0
		sib%param%atheta(i) = 0.0
		sib%param%btheta(i) = 0.0
		sib%param%physfrac(i) = 0.0
		sib%param%physfrac1(i) = 0.0
		sib%param%physfrac2(i) = 0.0  
		sib%param%phystype(i) = 0.0
	enddo

	do i=1,nsoil
		sib%param%rootf(i) = 0.0
		sib%param%rootr(i) = 0.0
		sib%param%respfactor(i) = 0.0
		sib%param%tkmg(i) = 0.0
		sib%param%tksatu(i) = 0.0
		sib%param%tkdry(i) = 0.0
		sib%param%csolid(i) = 0.0
	enddo
	
		!sib%param%tksoil(-nsnow+1:nsoil) = 0.0
		!sib%param%slamda(-nsnow+1:nsoil) = 0.0
		!sib%param%shcap(-nsnow+1:nsoil) = 0.0


!------------------------------------------------------------------
!                   PROGNOSTIC VARIABLES
!------------------------------------------------------------------
	sib%prog%tg = 0.0
	sib%prog%ta = 0.0
	sib%prog%tc = 0.0
	sib%prog%tha = 0.0
	sib%prog%sha = 0.0
	sib%prog%ea = 0.0
	sib%prog%snow_veg = 0.0
	sib%prog%tke = 0.0
	sib%prog%snow_mass  = 0.0
	sib%prog%snow_depth  = 0.0
	sib%prog%snow_age = 0.0
	sib%prog%nsl   = 0.0
	sib%prog%pco2ap = 0.0
	sib%prog%pco2ap_old = 0.0
	sib%prog%cas = 0.0
	sib%prog%cas_old  = 0.0
	sib%prog%expand	 = 0.0
	sib%prog%pco2m = 0.0
	sib%prog%sw_dwn = 0.0
	sib%prog%sw_dwn1  = 0.0
	sib%prog%sw_dwn2  = 0.0
	sib%prog%radvbc  = 0.0
	sib%prog%radvdc 	 = 0.0
	sib%prog%radnbc = 0.0
	sib%prog%radndc  = 0.0
	sib%prog%dlwbot = 0.0
	sib%prog%dlwbot1 = 0.0
	sib%prog%dlwbot2 = 0.0
	sib%prog%vdcsav = 0.0
	sib%prog%tm = 0.0
	sib%prog%tm1 = 0.0
	sib%prog%tm2  = 0.0
	sib%prog%thm = 0.0
	sib%prog%sh = 0.0
	sib%prog%sh1  = 0.0
	sib%prog%sh2 = 0.0
	sib%prog%em = 0.0
	sib%prog%ps  = 0.0
	sib%prog%ps1 = 0.0
	sib%prog%ps2 = 0.0
	sib%prog%psb = 0.0
	sib%prog%zb = 0.0
	sib%prog%ros = 0.0
	sib%prog%cupr = 0.0
	sib%prog%cupr1 = 0.0
	sib%prog%cupr2 = 0.0
	sib%prog%lspr = 0.0
	sib%prog%lspr1 = 0.0
	sib%prog%lspr2 = 0.0
	sib%prog%spdm  = 0.0
	sib%prog%spdm1 = 0.0
	sib%prog%spdm2 = 0.0
	sib%prog%tcc1  = 0.0
	sib%prog%tcc2 = 0.0
	sib%prog%d13cca = 0.0
	sib%prog%d13cm = 0.0
	
	do i=1,2
		sib%prog%bps(i) = 0.0
		sib%prog%capac(i) = 0.0
	enddo
	
	do i=1,6
		sib%prog%rst(i) = 0.0
	enddo
	
	
		!sib%prog%dz(-nsnow+1:nsoil) = 0.0
		!sib%prog%layer_z(-nsnow:nsoil) = 0.0
		!sib%prog%node_z(-nsnow+1:nsoil) = 0.0
		!sib%prog%vol_ice(-nsnow+1:nsoil) = 0.0
		!sib%prog%vol_liq(-nsnow+1:nsoil) = 0.0
		!sib%prog%www_ice(-nsnow+1:nsoil) = 0.0
		!sib%prog%www_liq(-nsnow+1:nsoil) = 0.0
		!sib%prog%td(-nsnow+1:nsoil = 0.0	
	
!------------------------------------------------------------------
!                   STATUS VARIABLES
!------------------------------------------------------------------

	sib%stat%coszbar = 0.0
	sib%stat%cosz = 0.0
	sib%stat%dayflag = 0.0
	sib%stat%julday = 0
	sib%stat%pt_num = 0
	





end subroutine init_var
