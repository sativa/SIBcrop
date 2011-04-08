module cfrax

    !  CFRAX calculates 13C and 12C fluxes and concentrations in the canopy,
    !  mixed layer, carbon assimilation (photosynthate), respired soil carbon,
    !  assuming that discrimination against 13C during photosynthesis is a 
    !  function of discrimination during diffusion and the pCO2i/pCO2c ratio.
    !  C4 discrimination against 13C only results from diffusion. 


    !itb...modified 01 Oct 99 to try to settle some mass-balance problems
    !itb...we've been having - new implicit eqn's for canopy C12 and C13.

    use kinds
    use sibtype

    use sib_const_module, only: &
        pdb,       & 
        kieclfbl,  &
        kiecstom,  &
        kieclphas, &
        kiecdis,   &
        kiecrbsco, &
        tref,      &
        pref,      &
        dtt

    implicit none

    !...LOCAL VARIABLES...
    real(kind=dbl_kind) :: co2a_conc   ! CAS CO2 concentration (mol m^-3)
    real(kind=dbl_kind) :: co2m_conc   ! ref lev CO2 concentration (mol m^-3)
    real(kind=dbl_kind) :: rca         ! C13/C12 ratio of canopy CO2 
    real(kind=dbl_kind) :: rcm         ! C13/C12 ratio of ref level CO2
    real(kind=dbl_kind) :: c13cm       ! concentration of C13 in mixed 
    !  layer CO2
    real(kind=dbl_kind) :: c12cm       ! concentration of C12 in mixed 
    !  layer CO2 
    real(kind=dbl_kind) :: rcresp      ! C13/C12 ratio of respiratory CO2
    real(kind=dbl_kind) :: c13resp     ! flux of C13 in respiratory CO2
    !    (mol m^-2 sec^-1)
    real(kind=dbl_kind) :: c12resp     ! flux of C12 in respiratory CO2
    !    (mol m^-2 sec^-1)
    real(kind=dbl_kind) :: c13ca       ! concentration of C13 in canopy CO2 
    real(kind=dbl_kind) :: c12ca       ! concentration of C12 in canopy CO2
    

    contains

    subroutine cfrax_physloop(sib,i,c3)

        use kinds
        use sibtype

        implicit none

        !--------------------------------------------------------

        type(sib_t), intent(inout) :: sib

        !--------------------------------------------------------  

        !...INPUT VARIABLES
        integer(kind=int_kind),intent(in) :: i        ! phystype loop index
        real(kind=dbl_kind),intent(in)    :: c3       ! C3 flag

        !...convert canopy and mixed layer CO2 pressures from Pa to moles/m^3 
        !...uses PV = nRT; at STP and 44.6 moles of gas per m3.  

        co2a_conc = sib%prog%pco2ap * (tref/sib%prog%ta)*(44.6/pref)
        co2m_conc = sib%prog%pco2m  * (tref/sib%prog%ta)*(44.6/pref)


        !...d13Cca and d13Cm are converted to concentrations (moles/m3) of 13C 
        !...and 12C by first calculating isotope ratios (13C/12C) of the canopy     
        !...(Ca) and mixed layer (m). 
        rca   = ((sib%prog%d13cca * pdb) / 1000.) + pdb

        c13ca = (rca * co2a_conc) / (1. + rca)
        c12ca = co2a_conc / (1. + rca)

        rcm   = ((sib%prog%d13cm * pdb) / 1000.) + pdb

        c13cm = (rcm * co2m_conc) / (1. + rcm)
        c12cm = co2m_conc / (1. + rcm)


        !...13c and 12c fluxes (moles/m2/sec) arising from respiration are        
        !...calculated using conversions between delta notation and epsilon 
        !...notation and 13C/12C ratios.

        rcresp  = ((sib%param%d13cresp * pdb) / 1000.) + pdb
        c13resp = (rcresp * sib%diag%respg) / (1. + rcresp)
        c12resp = sib%diag%respg / (1. + rcresp)

        ! C13/C12 discrimination for C3 plants.  The isotope effect during C3
        ! photosynthesis is a function of a combination of the isotope effects
        ! associated with molecular transport of CO2 across the leaf boundary
        ! layer (lfbl), into the stoma (stom), dissolution to in mesophyll H2O
        ! (dis), and transport in the liquid phase (lphas).  The isotope effect 
        ! during C4 photosynthesis is only a function (for now) of transport into
        ! the stoma.    note: IECpsC3 is the isotope effect for carbon isotopic  
        ! discrimination during photosynthesis of C3 plants.  Similarly for IECpsC4,
        ! but for C4 plants. 


        if(sib%diag%assimn(i) > 0.0) then

            if(c3 == 1.0) then
	    
                sib%diag%kiecps(i) = ( kieclfbl * sib%prog%pco2ap + &
                    (kiecstom-kieclfbl) * sib%diag%pco2s(i) + &
                    (kiecdis + kieclphas - kiecstom) * sib%diag%pco2i(i) + &
                    (kiecrbsco-kiecdis-kieclphas) * sib%diag%pco2c(i) ) &
                    / sib%prog%pco2ap

            else !C4 plants given constant KIE = 4.4per mil  
	    
                sib%diag%kiecps(i) = kiecstom
		 
            endif



        else 

            ! We need values for when Assimn < 0.0, i.e. nighttime.  We revert
	    ! to the del13C of respiration both for convenience and because the
	    ! nighttime flux to the atmosphere reflects the 13C/12C ratio of
	    ! respiration

                sib%diag%kiecps(i) = sib%param%d13cresp - sib%prog%d13cm


        endif


        ! calculates d13C of carbon and fluxes of 13C and 12C assimilated 
	! by phystype i (moles/m2/sec)

        sib%diag%rcassimn(i)   =   rca / ((-sib%diag%kiecps(i) / 1000.) + 1.)

        sib%diag%d13cassimn(i) =  ((sib%diag%rcassimn(i) - pdb) / pdb ) * 1000.

        sib%diag%c13assimn(i)  = ( sib%diag%rcassimn(i) * sib%diag%assimn(i))/(1. + sib%diag%rcassimn(i))

        sib%diag%c12assimn(i)  =  sib%diag%assimn(i) / (1.+ sib%diag%rcassimn(i))


    end subroutine cfrax_physloop




    !ITB...this is the 'final' cfrax stuff to be done after the physiology 
    !ITB...loops are finished.

    subroutine cfrax_final(sib)

        use kinds
        use sibtype

        implicit none
        !------------------------------------------------------------

        type(sib_t), intent(inout) :: sib

        !------------------------------------------------------------  


        !...LOCAL VARIABLES...
        integer(kind=int_kind) :: i

        !...CFRAX...CFRAX...CFRAX...CFRAX...CFRAX...CFRAX...CFRAX...CFRAX
        sib%diag%c13assimn(6) = 0.0
        sib%diag%c12assimn(6) = 0.0

        phys_loop: do i=1,5

            !print*,'physfrac:',i,sib%param%physfrac(i)
            if(sib%param%physfrac(i) == 0.0) cycle phys_loop
             
            sib%diag%c13assimn(6) = sib%diag%c13assimn(6) + &
                sib%diag%c13assimn(i) * sib%param%physfrac(i)

            sib%diag%c12assimn(6) = sib%diag%c12assimn(6) + &
                sib%diag%c12assimn(i) * sib%param%physfrac(i)

        enddo phys_loop

        sib%diag%rcassimn(6) = sib%diag%c13assimn(6)/   &
                                                sib%diag%c12assimn(6)
 

        sib%diag%d13cassimn(6) =  ((sib%diag%rcassimn(6)/  &
                                   pdb) - 1.0_dbl_kind)*1000.0_dbl_kind


        !  Canopy concentrations at time n+1 is calculated using an implicit 
        !  scheme.
        c13ca = (c13ca + (dtt / sib%param%z2) * (c13resp - &
            sib%diag%c13assimn(6) + (c13cm / sib%diag%ra ))) &
            / (1. + (dtt / sib%diag%ra ) / sib%param%z2)

        c12ca = (c12ca + (dtt / sib%param%z2) * (c12resp - &
            sib%diag%c12assimn(6) + (c12cm / sib%diag%ra ))) &
            / (1. + (dtt / sib%diag%ra ) / sib%param%z2)



        !  del13C of canopy is recalculated using concentrations of 13Cca and 
        !  12Cca.  The fluxes (moles/m2/sec) of 13C and 12C out of the canopy 
        !  (the turbulent flux), and the del13C value (per mil vs PDB) of this 
        !  flux are calculated. 

        sib%prog%d13cca = ((c13ca/c12ca - pdb) / pdb) *1000.


        !  Use the following if you want the net flux from the canopy to 
        !  be based on differences in 12C and 13C net fluxes from respiration
        !  and photosynthesis

        !           sib%diag%flux13c = sib%diag%respg * rcresp / (1. + rcresp) -
        !     &  (sib%diag%assimn(6) * sib%diag%rcassimn / (1. + sib%diag%rcassimn))

        !           sib%diag%flux12c = sib%diag%respg / (1. + rcresp) -
        !     &  (sib%diag%assimn(6) / (1. + sib%diag%rcassimn))

        !  Use the following if you want the net flux from the canopy to 
        !  be based on differences  in concentration gradients between the 
        !  canopy and overlying atmosphere.

        sib%diag%flux13c = ( c13ca - c13cm) / sib%diag%ra 
        sib%diag%flux12c = ( c12ca - c12cm) / sib%diag%ra

        sib%diag%flux_turb = sib%diag%flux13c + sib%diag%flux12c

    end subroutine cfrax_final


end module cfrax
