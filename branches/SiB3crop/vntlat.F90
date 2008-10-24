subroutine vntlat(sib,sib_loc)


    use kinds
    use sibtype

    use sib_const_module, only:  &
        vkrmn,  &
        snofac, &
        dtt,    &
        zwind,  &
        ztemp

    implicit none

    !----------------------------------------------------------------------

    type(sib_t), intent(inout) :: sib
    type(sib_local_vars)     ,intent(inout) :: sib_loc
                                       ! variables local to SiB

    !----------------------------------------------------------------------  

    !...LOCAL VARIABLES

    real(kind=dbl_kind) :: ps
    real(kind=dbl_kind) :: u2 
    real(kind=dbl_kind) :: cuni 
    real(kind=dbl_kind) :: temv 
    real(kind=dbl_kind) :: zzwind
    real(kind=dbl_kind) :: zztemp
    real(kind=dbl_kind) :: eps
    real(kind=dbl_kind) :: epsc
    real(kind=dbl_kind) :: epsg 

    eps    = 1. / snofac                                   

    !
    !   calculate ventilation mass flux
    !

    zzwind = sib%param%z2 - sib%param%zpd_adj + zwind
    zztemp = sib%param%z2 - sib%param%zpd_adj + ztemp


    call vmfcalz(sib,zzwind,zztemp,cuni)

    !                                                                       
    !     AERODYNAMIC RESISTANCE                                            
    !                             

    sib%diag%ra    = sib%prog%ros / sib%diag%ventmf 
    temv = (sib%param%z2 - sib%param%zpd_adj) / sib%param%z0
!    print*,'vntlat,temv:',sib%param%z2,sib%param%zpd_adj,sib%param%z0
    !itb...PATCH...make sure that temv is not negative
    temv = max(temv,1.00_dbl_kind)
    temv = log(temv) 
!print*,sib%prog%spdm, cuni,vkrmn,temv, log(temv) 
    u2     = sib%prog%spdm / (cuni * vkrmn)
    u2 = u2 * temv

    !itb...HARDWIRE PATCH...keeping u2 from being zero. That's bad...
 !   u2 = MAX(u2,1.0_dbl_kind)

    sib%diag%drag(1) = sib%prog%ros * sib%diag%cu * sib%diag%ustar

    sib_loc%fc = 1.0
    sib_loc%fg = 1.0                                  

    !...calculate leaf surface-CAS and ground-CAS resistance
    call rbrd(sib,u2)

    epsc = 1.
    epsg = 1. 
    !   this only makes sense for canopy leaves, since
    !   there can only be water OR snow, not both. switching epsc
    !   epsc to eps makes the hltm adapt to freezing/fusion.
    if(sib%prog%snow_veg > 0.0) epsc = eps
    if(sib%prog%nsl      < 0)   epsg = eps
    !
    !   compute resistance terms
    !
    sib%diag%rc = sib%prog%rst(6) + sib%diag%rb + sib%diag%rb

    sib%diag%rds = sib%diag%rsoil * sib_loc%fg + sib%diag%rd

    sib_loc%gect =  (1. - sib%diag%wc) / sib%diag%rc
    sib_loc%geci = epsc * sib%diag%wc / (sib%diag%rb + sib%diag%rb)

    sib_loc%gegs =  (1. - sib%diag%wg) / sib%diag%rds
    sib_loc%gegi = epsg * sib%diag%wg / sib%diag%rd

    sib_loc%coc = sib_loc%gect + sib_loc%geci

    !...calculate ecmass -- canopy evapotranspiration
    !...temporary value to be passed into phosib
    sib%diag%ecmass = (sib_loc%etc - sib%prog%ea) * sib_loc%coc *  &
        sib%prog%ros  * 0.622e0_dbl_kind /sib%prog%ps * dtt

    !
    !   calculate soil respiration
    !
    call respsib(sib)
    !
    !   calculation of canopy conductance and photosynthesis
    !

    call phosib(sib,sib_loc)

    if(sib%prog%ea > sib_loc%etc) sib_loc%fc = 0.0
    if(sib%prog%ea > sib_loc%etg) sib_loc%fg = 0.0
    sib%diag%hrr = sib%diag%hr                                                
    if (sib_loc%fg < 0.5) sib%diag%hrr = 1.0                           

end subroutine vntlat
