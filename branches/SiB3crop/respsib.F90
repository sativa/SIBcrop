!----------------------------------------------------------------------
subroutine respsib(sib)
!----------------------------------------------------------------------
! respsib calculates ground respiration based on annual balance between 
! canopy net assimilation and respiration [Denning et al., 1996].
! Instantaneous ground respiration using the "R-star" approach of
! Denning et al (1996) scaled by soil temperature and moisture [Raich et al., 1991]
! and adapted for use with the Bonan/CLM 10-layer soil thermodynamics module.
!
! References:
!   Raich et al. [1991) Ecological Applications 1(4) pp 399-429 
!                     appendix I - Terrestrial Ecosystem Model
!   Denning, A.S., G.J. Collatz, C. Zhang,D.A. Randall, J.A. Berry, 
!                     P.J. Sellers, G.D. Colello D.A. Dazlich, 1996: 
!                     Simulations of Terrestrial Carbon Metabolism and 
!                     Atmospheric CO2 in  a General Circulation Model. 
!                     Part 1: Surace Carbon Fluxes. Tellus, 48B, 521-542.
! Modifications:
!   Scott Denning (9/14/95) changed soil Q10 value for soil respiration from 2.0 to 2.4
!      following Raich and Schelsinger (1992, Tellus 44B, 81-89)
!   Kevin Schaefer (8/4/2004) changed wetness exponent to use soil moisture 
!      fraction rather than percent to minimize cpu time
!   Kevin Schaefer (8/4/2004) changed wetness exp to use wfrac rather than www_liq
!   Lara Prihodko changed moist to use wet_exp for each layer (10/28/04)

    use kinds
    use sibtype
    use sib_const_module, only : denh2o,denice

    implicit none    

    type(sib_t), intent(inout) :: sib

!...LOCAL VARIABLES
    integer:: j

    real(kind=dbl_kind) :: &
        woptzm    ! optimal soil moisture fraction for resp (wopt) to the skewness exponent (zm)
    real(kind=dbl_kind),dimension(nsoil) :: &
        wet_exp, & ! wetness exponent (b in eqn 8 from Denning et al [1996])
        wfrac,   & ! soil moisture fraction of saturation for soil layer
        moist      ! soil moisture scaling factor for resp (eqn 1.14b from Raich et al [1991])
!
! calculate optimal soil moisture fraction for resp (wopt) to the skewness exponent (zm)
! only needed once, calculate once to save cpu time
    woptzm = (sib%param%wopt/100.) ** sib%param%zm
!
!...initialize ground respiration (respg) to zero
    sib%diag%respg = 0.0
!
! loop through soil layers
    do j=1,nsoil
!
! calculate soil moisture fraction of saturation
        wfrac(j) = sib%prog%www_liq(j) / (sib%prog%dz(j) * sib%param%poros * denh2o) 

! calculate wetness exponent
        wet_exp(j) = ((wfrac(j)**sib%param%zm-woptzm)/(1.-woptzm))**2
        wet_exp(j) = min(wet_exp(j),10.0_dbl_kind)
!
! calculate soil moisture respiration scaling factor
        moist(j) = 0.8*sib%param%wsat**wet_exp(j) + 0.2
!
! calculate soil temperature respiration scaling factor
        sib%diag%soilQ10(j) = exp(0.087547 * (sib%prog%td(j) - 298.15))
!
! calculate total soil respiration scaling factor
        sib%diag%soilscale(j) = sib%diag%soilQ10(j) * moist(j)
!
! calculate soil layer respiration; add to total ground respiration
        sib%diag%respg = sib%diag%respg + sib%param%respfactor(j) *   &
            sib%diag%soilscale(j)

    enddo  ! soil layers

end subroutine respsib
