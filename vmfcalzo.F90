
subroutine vmfcalzo(sib,zzwind,zztemp) 


!*****************************************************************************
!     VENTILATION MASS FLUX,Ustar, and transfer coefficients for momentum 
!     and heat fluxes, based on by Louis (1979, 1982), and revised by Holtslag
!     and Boville(1993), and by Beljaars and Holtslag (1991).              
!  
!     Rerences:
!       Beljars and Holtslag (1991): Flux parameterization over land surfaces
!              for atmospheric models. J. Appl. Meteo., 30, 327-341.
!       Holtslag and Boville (1993): Local versus nonlocal boundary-layer 
!              diffusion in a global climate model. J. of Climate, 6, 1825-
!              1842.
!       Louis, J. F., (1979):A parametric model of vertical eddy fluxes in
!              atmosphere. Boundary-Layer Meteo., 17, 187-202.
!       Louis, Tiedke, and Geleyn, (1982): A short history of the PBL
!              parameterization at ECMWF. Proc. ECMWF Workshop on Boundary-
!              Layer parameterization, ECMWF, 59-79.
!
!     General formulation:
!        surface_flux = transfer_coef.*U1*(mean_in_regerence - mean_at_sfc.) 
!     Transfer coefficients for mommentum and heat fluxes are:
!        CU = CUN*Fm, and
!        CT = CTN*Fh
!        where  CUN and CTN are nutral values of momentum and heat transfers,
!           and Fm and Fh are stability functions derived from surface
!           similarity relationships.     
!*****************************************************************************

    use kinds
    use sibtype

    use physical_parameters, only:  &
        grav

    use sib_const_module, only:  &
        vkrmn

    implicit none


    !----------------------------------------------------------------------

    type(sib_t), intent(inout) :: sib

    !----------------------------------------------------------------------  



    real(kind=dbl_kind),intent(in) :: zzwind
    real(kind=dbl_kind),intent(in) :: zztemp


    !...LOCAL VARIABLES...
    real(kind=dbl_kind) :: zrib
    real(kind=dbl_kind) :: z1z0u
    real(kind=dbl_kind) :: z1z0urt
    real(kind=dbl_kind) :: cun
    real(kind=dbl_kind) :: temv
    real(kind=dbl_kind) :: rib
    real(kind=dbl_kind) :: fmomn
    real(kind=dbl_kind) :: ribtemp
    real(kind=dbl_kind) :: dm

    ! constants for surface flux functions, according to Holtslag and
    ! Boville (1993, J. Climate)
    ! constants for unstable function
    real(kind=dbl_kind),parameter:: bunstablM = 10.  
    real(kind=dbl_kind),parameter:: cunstablM = 75.

                ! constants for stable function
    real(kind=dbl_kind),parameter:: bstabl = 8.   
    real(kind=dbl_kind),parameter:: cstabl = 10.


    zrib = zzwind**2.0 / zztemp

    !   Ratio of reference height (zwind/ztemp) and roughness length:
    z1z0u = zzwind/ 0.0002   ! oceanic roughness length
    z1z0urt = sqrt( z1z0U )
    z1z0u = log( z1z0U )

    !   Neutral surface transfers for momentum CUN and for heat/moisture CTN:

    cun = vkrmn*vkrmn / (z1z0u*z1z0u )   !neutral Cm & Ct

    !   SURFACE TO AIR DIFFERENCE OF POTENTIAL TEMPERATURE.            
    !   RIB IS THE BULK RICHARDSON NUMBER, between reference height and surface.

    temv = sib%prog%tha * sib%prog%spdm * sib%prog%spdm                          
    temv = max(1.0e-6_dbl_kind,temv)
    rib = -sib%diag%thvgm * grav * zrib / temv 

    !   The stability functions for momentum and heat/moisture fluxes as
    !   derived from the surface-similarity theory by Luis (1079, 1982), and
    !   revised by Holtslag and Boville(1993), and by Beljaars and Holtslag 
    !   (1991).

    if(rib >= 0.0_dbl_kind) then                                           

        !  THE STABLE CASE. RIB IS USED WITH AN UPPER LIMIT              

        rib = min( rib, 0.5_dbl_kind)                   
        fmomn = (1. + cstabl * rib * (1.+ bstabl * rib))
        fmomn = 1. / fmomn
        fmomn = max(1.0e-4_dbl_kind,fmomn)

    else                                  

        !  THE UNSTABLE CASE.    

        ribtemp = abs(rib)
        ribtemp = sqrt( ribtemp )
        dm = 1. + cunstablm * cun * z1z0urt * ribtemp
        fmomn = 1. - (bunstablm * rib ) / dm

    end if    

    !   surface-air transfer coefficients for momentum CU, for heat and 
    !   moisture CT.

    sib%diag%cu = cun * fmomn 

    !   Ustar and ventlation mass flux: note that the ustar and ventlation 
    !   are calculated differently from the Deardoff's methods due to their
    !   differences in define the CU and CT.

    sib%diag%ustar   = sib%prog%spdm * sib%prog%spdm * sib%diag%cu 
    sib%diag%ustar   = sqrt( sib%diag%ustar ) 
    sib%diag%drag(2) = sib%prog%ros * sib%diag%cu * sib%diag%ustar

    ! Note there is no CHECK FOR VENTMF EXCEEDS TOWNSENDS(1962) FREE CONVECTION  
    ! VALUE, like DEARDORFF EQ(40B), because the above CU and CT included
    ! free convection conditions. 



end subroutine vmfcalzo
