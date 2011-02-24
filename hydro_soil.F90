
subroutine hydro_soil(sib)

!-----------------------------------------------------------------
!
!     soil water - based on CLM_HYDRO_SOIL.F90 code...
!
!-----------------------------------------------------------------
!
!     following taken (almost) verbatim from CLM_HYDRO_SOIL comments...
!
!     main subroutine used to execute the calculation of water 
!     processes over soil.
!
!     1) water within soil (subroutine soilwater)
! 
!     2) runoff
!        The original code was provided by Robert E. Dickinson based on
!        following clues: exponential decrease of Ksat, a water table
!        level determination level including highland and lowland levels
!        and fractional area of wetland (water table above the surface).
!        Runoff is parameterized from the lowlands in terms of precip
!        incident on wet areas and a base flow, where these are estimated
!        using ideas from TOPMODEL.
!
!     The original scheme was modified by Z.-L. Yang and G.-Y. Niu,
!     *  using a new method to determine water table depth and the 
!        fractional wet area (fwsoil).
!     *  computing runoff (surface and subsurface) from this fraction and
!        the remaining fraction (i.e. 1-fwsoil)
!     *  for the 1-fwsoil part, using BATS1e method to compute surface and
!        subsurface runoff.
!
!     The original code on soil moisture and runoff was provided by R.E.
!     Dickinson in July 1996.
!
!     Revision History:
!     15 September 1999: Yongjiu Dai; Initial code
!     12 November  1999: Z.-L. Yang and G.-Y. Niu
!     15 December  1999: Paul Houser and Jon Radakovich; F90 revision
!     25 January   2002: Ian Baker; Integration into SiB
!
!----------------------------------------------------------------------

use kinds
use sibtype

use physical_parameters, only: &
    hltm

use sib_const_module, only: &
    nsoil,  &
    denice, &
    denh2o, &
    wtfact, &
    dtt,    &
    dti

implicit none

!----------------------------------------------------------------------

type(sib_t), intent(inout) :: sib

!----------------------------------------------------------------------  


!...LOCAL VARIABLES
integer(kind=int_kind) :: j,i      ! loop variable
real(kind=dbl_kind)    :: zwice    ! total ice mass in soil column
!   (kg m^-2)
real(kind=dbl_kind)    :: www_tot(1:nsoil)  
! total pore space occupied
!  by ice+liquid (-)
real(kind=dbl_kind)    :: dw_liq(1:nsoil)   
! change in layer liquid
!  -volumetric (m^3/m^3)
real(kind=dbl_kind)    :: dzmm(1:nsoil)
! soil layer thickness in mm (mm)
real(kind=dbl_kind)    :: zmm(1:nsoil)
! node depth in mm (mm)
real(kind=dbl_kind)    :: wmean    ! averaged soil wetness in top 
!  layers
real(kind=dbl_kind)    :: zmean    ! top soil layers contributing to 
!  runoff
real(kind=dbl_kind)    :: fz       ! factor (-)
real(kind=dbl_kind)    :: zwt      ! water table depth (m)
real(kind=dbl_kind)    :: infil    ! infiltration into soil 
!  (kg m^-2 sec^-1)
real(kind=dbl_kind)    :: hk(1:nsoil) 
! hydraulic conductivity 
!  (mm h20 sec^-1)
real(kind=dbl_kind)    :: dhkdw (1:nsoil)  
! d(hk)/d(water content)
real(kind=dbl_kind)    :: q_wet    ! subsurface runoff from 'wet' part
!  (mm H2O sec^-1)
real(kind=dbl_kind)    :: q_dry    ! subsurface runoff from 'dry' part
!  (mm H2O sec^-1)
real(kind=dbl_kind)    :: hksum    ! total hk for layers 6-9
real(kind=dbl_kind)    :: zsat     ! hk, weighted for soil thickness
real(kind=dbl_kind)    :: wsat     ! hk, weighted for soil wetness
real(kind=dbl_kind)    :: dzksum   ! hk, weighted for soil thickness
real(kind=dbl_kind)    :: fwsoil   ! saturation fraction
real(kind=dbl_kind)    :: watmin   ! minimum soil moisture 
real(kind=dbl_kind)    :: xs       ! excess soil moisture 
!  (kg m^-2) 


real(kind=dbl_kind)    :: roffo_t   ! placeholder for overland
! runoff (kg/m^2)

!---------------------------------------------------------------------


    zwice = 0.0  ! sum of ice mass of soil

    do j=1,nsoil
        zwice = zwice + sib%prog%www_ice(j)

        sib%prog%vol_ice(j) = min(sib%param%poros, sib%prog%www_ice(j)/(sib%prog%dz(j)*denice))
        sib%diag%eff_poros(j) = sib%param%poros - sib%prog%vol_ice(j)
        sib%prog%vol_liq(j) = min(sib%diag%eff_poros(j), &
            sib%prog%www_liq(j)/(sib%prog%dz(j)*denh2o))
        if (sib%prog%vol_liq(j) == 0.0 .and. sib%prog%www_liq(j) > 0.0 ) then
            sib%diag%roff = sib%diag%roff + sib%prog%www_liq(j) 
            sib%prog%www_liq(j) = 0.0
        endif

        www_tot(j) = min(1.0_dbl_kind,(sib%prog%vol_ice(j)+ &
            sib%prog%vol_liq(j))/sib%param%poros)
    enddo


    !...determine water table
    wmean = 0.0
    fz    = 1.0
    do j=1,nsoil
        wmean = wmean + www_tot(j)*sib%prog%dz(j)
    enddo
    zwt = fz * (sib%prog%layer_z(nsoil) - wmean)

    !...saturation fraction

    fwsoil = wtfact * min(1.0_dbl_kind,exp(-zwt))
    !      print'(2(a,g16.6))','water table depth, m=',zwt,' sat fraction=',fwsoil


    !...these soil calculations are hardwired for a 10-layer soil model
    wmean = 0.0
    zmean = 0.0

    do j=1,3
        zmean = zmean + sib%prog%dz(j)
        wmean = wmean + www_tot(j)*sib%prog%dz(j)
    enddo

    wmean = wmean/zmean

    !itb...modifying infiltration from CLM. We put precipitation/snowmelt 
    !itb...into surface interception store (sib%prog%capac(2)). Any infiltration
    !itb...will come from that reservoir. The hope is that we can hold
    !itb...onto puddles from one timestep to the next. Capac(2) is constrained
    !itb...by satcap(2), which is set in subroutine RADA2.


    !itb...THIS IS A PROBLEM-MESH BETWEEN SiB/CLM IS NOT CLEAN...

    !itb...old code...
    !      roffo_t  = max(0.0_dbl_kind,fwsoil*sib%diag%www_inflow) + &
    !                 max(0.0_dbl_kind,(1.0 - fwsoil)  &
    !               * min(1.0_dbl_kind,wmean**4*sib%diag%www_inflow))

    !itb...I think that roffo_t represents everything that can't 
    !itb...infiltrate during this timestep.
    !itb...new code
    roffo_t  = max(0.0_dbl_kind,fwsoil*sib%prog%capac(2)*dti) +  &
        max(0.0_dbl_kind,(1.0 - fwsoil) &
        * min(1.0_dbl_kind,wmean**4*sib%prog%capac(2)*dti))

    roffo_t  = roffo_t * dtt   ! kg/m^2

    !      sib%diag%roffo = sib%diag%roffo + roffo_t

    !...infiltration into surface soil layer 
    infil = sib%prog%capac(2) - roffo_t  ! units: kg/m^2

    sib%prog%capac(2) = sib%prog%capac(2) - infil 

    sib%diag%roffo = sib%diag%roffo + max(0.0_dbl_kind, &
        (sib%prog%capac(2) - sib%param%satcap(2)*denh2o))  ! kg/m^2

!    sib%prog%capac(2) = max(0.0_dbl_kind,sib%prog%capac(2) - &
!        sib%param%satcap(2)*denh2o)

    sib%prog%capac(2) = min(sib%prog%capac(2),sib%param%satcap(2)*denh2o)
    sib%prog%capac(2) = MAX(sib%prog%capac(2),0.0_dbl_kind)

    !...set up r, a, b, and c vectors (following Bonan's (1996) soil)
    !...for tridiagonal matrix.
    !...(length units will be millimeters)

    do j = 1,nsoil
        zmm(j)  = sib%prog%node_z(j) *1000.0
        dzmm(j) = sib%prog%dz(j) * 1000.0
    enddo

    !...need to convert infil from kg/m^2 to kg/m^2/sec
    if(infil > 0.0) then
        infil = infil * dti
    endif

    call soilwater( sib, infil, hk, dhkdw, dw_liq, dzmm, zmm )


    !...update the liquid water mass (kg/m^2)
    do j=1,nsoil
        sib%prog%www_liq(j) = sib%prog%www_liq(j) + &
            dw_liq(j)*sib%prog%dz(j)*denh2o
    enddo



    !...streamflow and total runoff
    !...The amount of streamflow is assumed to be maintained by flow 
    !...from the lowland water table with different levels contributing 
    !...according to their thickness and saturated hydraulic conductivity,
    !...i.e. a given layer below the water table interface loses water
    !...at a rate per unit depth given by drainage*hk/(sum over all layers 
    !...below this water table of hk*dz).  Because this is a slow smooth 
    !...process, and not strongly coupled to water in any one layer, it 
    !...should remain stable for explicit time differencing.  Hence, for 
    !...simplicity it is removed explicitly prior to the main soil water 
    !...calculation.  
    !...Another assumption: no subsurface runoff for ice mixed soil.
    !...Zong-Liang Yang and G.-Y. Niu

    sib%diag%qqq = 0.0
    q_wet  = 0.0
    q_dry  = 0.0
    hksum  = 0.0

    !...HARDWIRE...
    !...taking out streamflow...

    !      do j=6,nsoil-1
    !       hksum = hksum + hk(j)
    !      enddo

    !      if (zwice <= 0.0 .and. hksum > 0.0 )then
    !        zsat = 0.0
    !        wsat = 0.0
    !        dzksum = 0.0
    !        do j=6,nsoil-1
    !          zsat = zsat + sib%prog%dz(j) *hk(j)
    !          wsat = wsat + www_tot(j)*sib%prog%dz(j)*hk(j)
    !          dzksum = dzksum + hk(j)*sib%prog%dz(j)
    !        enddo
    !        wsat = wsat / zsat

    !        q_wet = fwsoil * 1.0e-5 * exp(-zwt)
    !        q_dry = (1.-fwsoil) * 4.0e-2 * wsat **(2.0*sib%param%bee+3)
    !        sib%diag%qqq = q_wet + q_dry

    !...streamflow...
    !        do j=6,nsoil-1
    !          sib%prog%www_liq(j) = sib%prog%www_liq(j) -dtt*sib%diag%qqq*sib%prog%dz(j)  &
    !                                            *hk(j)/dzksum
    !        enddo
    !      endif



    !...limit www_liq to be greater than or equal to watmin
    !...get water needed to bring www_liq equal to watmin from lower level

    watmin = 0.0
    !      watmin = sib%param%vwcmin
    do j=1,nsoil-1
        !        if(sib%prog%www_liq(j) < 0.0 ) then
        if(sib%prog%www_liq(j) < watmin*sib%prog%dz(j)*denh2o ) then
            xs = watmin * sib%prog%dz(j) * denh2o - sib%prog%www_liq(j)
        else
            xs = 0.0
        endif
        sib%prog%www_liq(j)   = sib%prog%www_liq(j)   + xs
        sib%prog%www_liq(j+1) = sib%prog%www_liq(j+1) - xs
    enddo

    j = nsoil
    if(sib%prog%www_liq(j) < watmin) then
        xs = watmin * sib%prog%dz(j) * denh2o - sib%prog%www_liq(j)
    else
        xs = 0.0
    endif
    sib%prog%www_liq(j) = sib%prog%www_liq(j) + xs
    sib%diag%qqq = sib%diag%qqq - xs

    !...determine water in excess of saturation

    !itb...I don't like how CLM included SATCAP here...
    !      xs = max(0.0_dbl_kind,sib%prog%www_liq(1)-                           &
    !                (sib%param%satcap(2)+sib%diag%eff_poros(1) * sib%prog%dz(1)*denh2o))
    !      if(xs > 0.0) sib%prog%www_liq(1) = &
    !               sib%param%satcap(2)+sib%diag%eff_poros(1) *sib%prog%dz(1)*denh2o

    xs = max(0.0_dbl_kind,sib%prog%www_liq(1)-sib%diag%eff_poros(1) &
        * sib%prog%dz(1)*denh2o)

    if(xs > 0.0) sib%prog%www_liq(1) = sib%diag%eff_poros(1) * &
        sib%prog%dz(1) * denh2o

    do j=2,nsoil
        xs = xs + max(sib%prog%www_liq(j)-sib%diag%eff_poros(j) &
            *sib%prog%dz(j)*denh2o,0.0_dbl_kind)

        sib%prog%www_liq(j) = min(sib%diag%eff_poros(j)*sib%prog%dz(j)*denh2o &
            ,sib%prog%www_liq(j))
    enddo

    !...sub-surface runoff and drainage
    sib%diag%qqq = sib%diag%qqq + xs                    
    sib%diag%roff = sib%diag%roff + sib%diag%qqq

    !...renew ice and liquid mass due to condensation
    !      if(sib%prog%nsl == 0) then
    !        sib%prog%www_liq(1) = sib%prog%www_liq(1) + dew*dtt
    !        sib%prog%www_ice(1) = sib%prog%www_ice(1) + (deposition -sublimation)*dtt
    !      endif


end subroutine hydro_soil
