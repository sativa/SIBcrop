!=================================================
subroutine mapper( lat,                 &
                   DOY,                 &
                   prevNDVI,            &
                   curNDVI,             &
                   fVCover,             &
                   ChiL,                &
                   LTran,               &
                   LRef,                &
                   MorphTab,            &
                   AeroVar,             &
                   LAIgrid,             &
                   fVCovergrid,         &
                   TimeVar )
!=================================================
! calculates time dependant boundary condition variables for SiB.

use kinds
IMPLICIT NONE

! begin input variables
real lat         ! center latitude of grid cell
real(kind=real_kind) :: curNDVI        ! FASIR NDVI values for a grid cell
real(kind=real_kind) :: prevNDVI    ! previous month's NDVI value
real(kind=dbl_kind) :: fVCover
real(kind=dbl_kind) :: ChiL
real LTran(2,2)
real LRef(2,2)

! begin input biome dependant, physical morphology variables
type biome_morph_var
   real (kind=real_kind) :: zc        ! Canopy inflection height (m)
   real (kind=real_kind) :: LWidth    ! Leaf width
   real (kind=real_kind) :: LLength   ! Leaf length	
   real (kind=real_kind) :: LAImax    ! Maximum LAI0.07

   real (kind=real_kind) :: stems     ! Stem area index
   real (kind=real_kind) :: NDVImax   ! Maximum NDVI
   real (kind=real_kind) :: NDVImin   ! Minimum NDVI
   real (kind=real_kind) :: SRmax     ! Maximum simple ratio
   real (kind=real_kind) :: SRmin     ! Minimum simple ratio
end type biome_morph_var

type(biome_morph_var) MorphTab

! begin input aerodynamic parameters
type aero_var
   real (kind=real_kind) :: zo       ! Canopy roughness coeff 
   real (kind=real_kind) :: zp_disp  ! Zero plane displacement
   real (kind=real_kind) :: RbC      ! RB Coefficient
   real (kind=real_kind) :: RdC      ! RC Coefficient
end type aero_var

type(aero_var) AeroVar(50,50) ! aerodynamic interpolation tables

real (kind=real_kind) :: LAIgrid(50)       ! grid of LAI values for lookup table
real (kind=real_kind) :: fVCovergrid(50)   ! grid of fVCover values for 
                                         !  interpolation table

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

! begin internal variables
real (kind=real_kind) ::DOY         ! Day of Year (DOY) of ndvi input map
real (kind=real_kind) ::prevfPAR    ! previous month's fPAR value
real(kind=real_kind), parameter :: fPARmax=0.95
                 ! Maximum possible FPAR corresponding to 98th percentile
real(kind=real_kind), parameter :: fPARmin=0.01
                 ! Minimum possible FPAR corresponding to 2nd percentile
!     For more information on fPARmin and fPARmax, see
!     Sellers et al. (1994a, pg. 3532); Los (1998, pg. 29, 37-39)

   !--------------------------------------------------------------
   ! Calculate time dependant variables
   !--------------------------------------------------------------
   ! Calculate first guess fPAR 
   ! use average of Simple Ratio (SR) and NDVI methods.

   !kdcorbin, 02/11 - initialize prevfPAR
   prevfPAR = 0.
   call AverageAPAR (prevNDVI, MorphTab%NDVImin, MorphTab%NDVImax,   &
                     MorphTab%SRmin, MorphTab%SRmax, fPARmax,        &
                     fParmin, prevfPAR)

   call AverageAPAR (curNDVI, MorphTab%NDVImin, MorphTab%NDVImax,    &
                     MorphTab%SRmin, MorphTab%SRmax, fPARmax,        &
                     fParmin, TimeVar%fPAR)


   ! Calculate leaf area index (LAI) and greeness fraction (Green)
   !   See S. Los et al 1998 section 4.2.

   !   Select previous month

   call laigrn (TimeVar%fPAR, prevfPAR, fPARmax, fVCover,         &
                MorphTab%stems, MorphTab%LAImax, TimeVar%Green,   &
                TimeVar%LAI)


   ! Interpolate to calculate aerodynamic, time varying variables
   call AeroInterpolate (TimeVar%LAI, fVCover, LAIgrid,fVCovergrid,   &
                         AeroVar, TimeVar%zo, TimeVar%zp_disp,        &
                         TimeVar%RbC, TimeVar%RdC)

   ! Calculate mean leaf orientation to par flux (gmudmu)
   call gmuder (lat, DOY, ChiL, TimeVar%gmudmu)

   ! recalculate fPAR adjusting for Sun angle, vegetation cover fraction,
   ! and greeness fraction, and LAI
   call aparnew (TimeVar%LAI, TimeVar%Green, LTran, LRef,   &
                 TimeVar%gmudmu, fVCover, TimeVar%fPAR,     &
                 fPARmax, fPARmin)


   return
end subroutine mapper


!-SUBROUTINE: averageapar--------------------------------------------

subroutine averageapar ( ndvi, ndvimin, ndvimax, srmin, srmax, &
        fparmax, fparmin, fpar )

    !----------------------------------------------------------------
    !
    ! Calucluates Canopy absorbed fraction of Photosynthetically
    ! Active Radiation (fPAR) using an average of the Simple Ratio (sr)
    ! and NDVI methods (Los et al. (1999), eqn. 5-6). The empirical
    ! SR method assumes a linear relationship between fPAR and SR.
    ! The NDVI assumes a linear relationship between fPAR and NDVI.
    !
    !----------------------------------------------------------------

    use kinds
    
    implicit none

    !-Parameters-----------------------------------------------------
    real(kind=real_kind) :: ndvi        ! normalized NDVI for vegetation type

    real(kind=real_kind) :: ndvimin     ! minimum NDVI for vegetation type
    real(kind=real_kind) :: ndvimax     ! maximum NDVI for vegetation type
    real(kind=real_kind) :: srmin       ! minimum SR for vegetation type
    real(kind=real_kind) :: srmax       ! maximum SR for vegetation type
    real(kind=real_kind) :: fparmax     ! maximum possible fPAR, corresponding
                                        !       to 98th percentile
    real(kind=real_kind) :: fparmin     ! minimum possible fPAR, corresponding
                                        !       to 2nd percentile

    real(kind=real_kind) :: fpar        ! canopy absorbed fraction of PAR

    !-Local Variables------------------------------------------------
    real(kind=real_kind) :: locndvi     ! local copy of ndvi
    real(kind=real_kind) :: sr          ! simple reatio of near IR and
                                        !       visible radiances
    real(kind=real_kind) :: ndvifpar    ! fPAR from NDVI method
    real(kind=real_kind) :: srfpar      ! fPAR from SR method

    ! Switch to local value of ndvi to prevent and changes
    locndvi = ndvi

    ! Insure calculated NDVI value falls within physical limits
    ! for vegetation type
    locndvi = max(locndvi, ndvimin)
    locndvi = min(locndvi, ndvimax)

    ! Calculate simple ratio (SR)
    sr=(1.+LocNDVI)/(1.-LocNDVI)

    ! Calculate fPAR using SR method (Los et al. (1999), eqn. 5)
    srfpar = (sr - srmin) * (fparmax - fparmin) / (srmax - srmin) + fparmin


    ! Calculate fPAR using NDVI method (Los et al. (1999), eqn. 6)
    ndvifpar = (locndvi - ndvimin) * (fparmax - fparmin) /  &
        (ndvimax - ndvimin) + fparmin

    ! Take the mean of the 2 methods
    fpar = 0.5 * (srfpar + ndvifpar)

    return

end subroutine averageapar


!==================================================
subroutine laigrn (fPAR,fPARm,fPARmax,fVCover,stems, LAImax,Green,LAI)
!==================================================
! calculate leaf area index (LAI) and greenness fraction (Green) from fPAR. 
! LAI is linear with vegetation fraction and exponential with fPAR.
! See Sellers et al (1994), Equations 7 through 13.

use kinds                                                                      
implicit none

! begin input variables
real fPAR     ! fraction of PAR absorbed by plants at current time
real fPARm    ! fraction of PAR absorbed by plants at previous time
real fPARmax  ! maximum possible FPAR corresponding to 98th percentile
real(kind=dbl_kind) :: fVCover  ! vegetation cover fraction
real stems    ! stem area index for the specific biome type
real LAImax   ! maximum total leaf area index for specific biome type

! begin output variables
real Green    ! greeness fraction of the total leaf area index
real LAI      ! area average total leaf area index

! begin internal variables
real LAIg     ! green leaf area index at current time
real LAIgm    ! green leaf area index at previous time
real LAId     ! dead leaf area index at current time

! Calculate current and previous green leaf area index (LAIg and LAIgm):
! LAIg is log-linear with fPAR.  Since measured fPAR is an area average, 
! divide by fVCover to get portion due to vegetation.  Since fVCover can
! be specified, check to assure that calculated fPAR does not exceed fPARMax.

   if(fPAR/fVCover.ge.fPARmax) then
      LAIg=LAImax
   else
      LAIg=alog(1.-fPAR/real(fVCover))*LAImax/alog(1-fPARmax)
   endif

   if(fPARm/fVCover.ge.fPARmax) then
      LAIgm=LAImax
   else
      LAIgm=alog(1.-fPARm/real(fVCover))*LAImax/alog(1-fPARmax)
   endif


   ! Calculate dead leaf area index (LAId):
   ! If LAIg is increasing or unchanged, the vegetation is in growth mode.
   ! LAId is then very small (very little dead matter).
   ! If LAIg is decreasing, the peak in vegetation growth has passed and
   ! leaves have begun to die off.  LAId is then half the change in LAIg,
   ! assuming half the dead leaves fall off.

   !     Growth mode dead leaf area index:
   if (LAIg.ge.LAIgm) LAId=0.0001

   !     die-off (post peak growth) dead leaf area index:
   if (LAIg.lt.LAIgm) LAId=0.5*(LAIgm-LAIg)

   ! Calculate area average, total leaf area index (LAI):
   LAI=(LAIg+LAId+stems)*fVCover

   ! Calculate greeness fraction (Green):
   ! Greeness fraction=(green leaf area index)/(total leaf area index)
   Green=LAIg/(LAIg+LAId+stems)

   return                                                                    
end subroutine laigrn



!-SUBROUTINE: aerointerpolate----------------------------------------

subroutine aerointerpolate ( lai, fvcover, laigrid, fvcovergrid, &
        aerovar, zo, zp_disp, rbc, rdc )

    !----------------------------------------------------------------
    !
    ! This subroutine calculates the aerodynamic parameters by
    ! bi-linear interpolation from a lookup tabke of previously
    ! calculated values.
    ! The interporation table is a numpts x numpts LAI/fVCover grid
    ! with LAI ranging from 0.02 to 10 and fVCover ranging from
    ! 0.01 to 1.
    !
    !----------------------------------------------------------------
    
    use kinds

    implicit none

    !-Parameters-----------------------------------------------------
    real(kind=real_kind) :: lai         ! actual area averages LAI
    real(kind=dbl_kind)  :: fvcover     ! vegetation cover fraction
    real(kind=real_kind) :: laigrid(50) ! grid of LAI values
    real(kind=real_kind) :: fvcovergrid(50) ! grid of fVCover values

    type aero_var
    real(kind=real_kind) :: zo          ! canopy roughness coefficient
    real(kind=real_kind) :: zp_disp     ! zero plane displacement
    real(kind=real_kind) :: rbc         ! RB coefficient
    real(kind=real_kind) :: rdc         ! RD coefficient
    end type aero_var

    type(aero_var) aerovar(50,50)       ! interpolation table

    real(kind=real_kind) :: rbc         ! interpolated RB coefficient
    real(kind=real_kind) :: rdc         ! interpolated RD coefficient
    real(kind=real_kind) :: zo          ! interpolated roughness length
    real(kind=real_kind) :: zp_disp     ! interpolated zero plane displacement

    !-Local Variables------------------------------------------------
    integer(kind=int_kind) :: i         ! index for LAI grid
    integer(kind=int_kind) :: j         ! index for fVCover grid
    real(kind=real_kind) :: loclai      ! local LAI
    real(kind=dbl_kind) :: locfvcover  ! local fVCover
    real(kind=dbl_kind) :: temp
    real(kind=real_kind) :: dlai        ! grid spacing between LAI values
    real(kind=real_kind) :: dfvcover    ! grid spacing between fVCover values

    ! Calculate grid spacing (assumed fixed)
    dlai = laigrid(2) - laigrid(1)
    dfvcover = fvcovergrid(2) - fvcovergrid(1)

    ! Assign input LAI and fVCover to local variables and make sure 
    ! they lie within the limits of the interpolation tables, assuring
    ! the LAI and fVConver values returned from the subroutine are not
    ! modified.
    loclai = max(lai, 0.02)
    temp = 0.01
    locfvcover = max(fvcover, temp)

    ! Determine the nearest array location for the desired LAI and fVCover
    i = int(loclai / dlai) + 1
    j = int(locfvcover / dfvcover) + 1
    j = min(j, 49)

    ! Interpolate RbC variable
    call interpolate( laigrid(i), loclai, dlai, &
        fvcovergrid(j), locfvcover, dfvcover,   &
        aerovar(i,j)%rbc, aerovar(i+1,j)%rbc,   &
        aerovar(i,j+1)%rbc, aerovar(i+1,j+1)%rbc, rbc )

    ! Interpolate RdC variable
    call interpolate( laigrid(i), loclai, dlai, &
        fvcovergrid(j), locfvcover, dfvcover,   &
        aerovar(i,j)%rdc, aerovar(i+1,j)%rdc,   &
        aerovar(i,j+1)%rdc, aerovar(i+1,j+1)%rdc, rdc )

    ! Interpolate roughness length
    call interpolate( laigrid(i), loclai, dlai, &
        fvcovergrid(j), locfvcover, dfvcover,   &
        aerovar(i,j)%zo, aerovar(i+1,j)%zo,     &
        aerovar(i,j+1)%zo, aerovar(i+1,j+1)%zo, zo )

    ! Interpolate zero plane displacement
    call interpolate( laigrid(i), loclai, dlai,       &
        fvcovergrid(j), locfvcover, dfvcover,         &
        aerovar(i,j)%zp_disp, aerovar(i+1,j)%zp_disp, &
        aerovar(i,j+1)%zp_disp, aerovar(i+1,j+1)%zp_disp, zp_disp )

    return

end subroutine aerointerpolate


!=================================================      
subroutine gmuder (Lat, DOY, ChiL, gmudmu)
!=================================================      
! calculates time mean leaf projection relative to the Sun.

use kinds
implicit none

! begin input variables
real Lat      ! latitude in degrees
real DOY      ! day-of-year (typically middle day of the month)
real(kind=dbl_kind) :: ChiL  ! leaf angle distribution factor

! begin output variables
real gmudmu   ! time mean projected leaf area normal to Sun

! begin internal variables
integer itime ! time counter
real gtime     ! time from 0:00 Greenwhich Mean Time (GMT)
real coshr    ! cosine of the Greenwhich Meridian (GM) Hour angle
real mu       ! cosine of the Solar zenith angle
real chiv     ! dummy variable for leaf angle distribution factor
real dec      ! declination of the Sun (Solar Declination)
real sin_dec  ! sine of the solar declination
real cos_dec  ! cosine of the solar declination
real pi       ! the constant pi
real pi180    ! conversion factor from degrees to radians
real aa       ! minimum possible LAI projection vs. cosine Sun angle
real bb       ! slope leaf area projection vs. cosine Sun angle
real cloud    ! (?) Cloud cover fraction
real fb       ! (?) mean cosine of Sun Angle
real swdown   ! (?) downward shortwave flux
real pardif   ! (?) PAR flux difracted into canopy by clouds
real pardir   ! (?) PAR flux directly onto canopy
real difrat   ! (?) fraction of shortwave flux defracted by clouds
real vnrat    ! (?) shortwave flux ratio of some sort
real tor      ! (?) TBD
real topint   ! (?) Total flux onto canopy adjusted for sun angle
real botint   ! (?) total PAR flux onto canopy during 24 hour period

! Assign values to constants
data pi /3.141592/

   pi180=pi/180. 
   cloud=0.5

   ! Calculate solar declination in radians
   dec=pi180*23.5*sin(1.72e-2*(DOY-80))

   ! Calculate sine and cosine of solar declination
   sin_dec=sin(dec)                                                         
   cos_dec=cos(dec)

   ! Begin time loop to average leaf projection over 24 hour period
   topint=0.
   botint=0.

   do itime=1, 48, 1                                                     

      ! Calculate time from zero Greenwhich Mean Time (GMT)
      gtime=0.5*real(itime) 

      ! Calculate cosine of hour angle of Grenwhich Meridion (GM)
      coshr=cos(-pi+gtime/24.*2.*pi)

      ! Calculate cosine of the Sun angle (mu)
      !     longitude=GM=0 degrees, latitude=Lat
      mu=sin(Lat*pi180)*sin_dec+cos(Lat*pi180)*cos_dec*coshr

      ! Ensure the cosine of Sun angle is positive, but not zero
      !     e.g., daylight, sun angle<=89.4 degrees (about start disc set/rise)
      mu=amax1(0.01, mu)                                            

      ! It looks like this code calculates the direct and difracted PAR based
      ! on the solar constant and a cloud fraction of 0.5 at the top and
      ! bottom of the atmosphere.  During night, mu=0.01, a constant.  These
      ! calculations do not match the definition of G(mu)/mu described in 
      ! Bonan (1996) and Sellers (1985).
      tor = 0.7**(1./mu)
      swdown = 1375.*mu*(tor+0.271-0.294*tor)
      difrat = 0.0604/(mu-0.0223)+0.0683
      difrat = max(difrat,0.)
      difrat = min(difrat,1.)
      difrat = difrat+(1.-difrat)*cloud
      vnrat  = (580.-cloud*464.)/((580.-cloud*499.) + (580.-cloud*464.))
      pardir = (1.-difrat)*vnrat*swdown
      pardif = difrat*vnrat*swdown
      topint = topint+pardir*mu+pardif*0.5
      botint = botint+pardir+pardif

   enddo                                                                 

   ! Calculate what looks like a mean value of something
   fb=topint/botint

   ! Calculate min and slope of LAI projection 
   chiv=ChiL                                                               
   if (abs(chiv) .le. 0.01) chiv=0.01
   !   calculate minimum value of projected leaf area
   aa=0.5-0.633*chiv-0.33*chiv*chiv
   !   calculate slope of projected leaf area wrt cosine sun angle
   bb=0.877*(1.-2.*aa)                                             

   ! Calculate mean projected LAI in Sun direction assuming fb approximates
   ! the mean cosine of the sun angle
   gmudmu=(aa+bb*fb)/fb                                            

   return   
                                                                    
end subroutine gmuder



!-SUBROIUTINE: aparnew-----------------------------------------------

subroutine aparnew ( lai, green, ltran, lref, &
       gmudmu, fvcover, fpar, fparmax, fparmin )

    !----------------------------------------------------------------
    !
    ! Recomputes the Canopy absorbed fraction of Photosynthetically
    ! Active Radiation (fPAR), adjusting for solar zenith angle and
    ! the vegetation cover fraction (fVCover) using a modified form
    ! of Beer's law.
    ! See Selleres et al. Part II (1996), eqns. 9-13.
    !
    !----------------------------------------------------------------

    use kinds

    implicit none

    !-Parameters-----------------------------------------------------
    real(kind=real_kind) :: lai            ! Leaf Area Index
    real(kind=real_kind) :: green       ! Greeness fraction of LAI
    real(kind=real_kind) :: fpar   ! Area average canopy absorbed fraction of PAR

    ! For LTran and LRef:
    !   (1,1) : shortwave, green plants
    !   (2,1) : longwave, green plants
    !   (1,2) : shortwave, brown plants
    !   (2,2) : longwave, brown plants
    real(kind=real_kind) :: ltran(2,2)  ! Leaf transmittance of 
                                        !       green/brown plants
    real(kind=real_kind) :: lref(2,2)   ! Leaf reflectance for
                                        !       green/brown plants

    real(kind=real_kind) :: scatp  ! Canopy transmittance + reflectance  wrt PAR
    real(kind=real_kind) :: scatg  ! Ground transmittance + reflectance wrt PAR
    real(kind=real_kind) :: park   ! Mean canopy absorption optical depth wrt PAR
    real(kind=real_kind) :: gmudmu      ! Daily time-mean canopy optical depth
    real(kind=dbl_kind)  :: fvcover     ! Canopy cover fraction
    real(kind=real_kind) :: fparmax     ! Maximum possible FPAR corresponding
                                        !       to 98th percentile
    real(kind=real_kind) :: fparmin     ! Minimum possible FPAR corresponding
                                        !       to 2nd percentile

    ! Calculate canopy transmittance + reflectance coefficient wrt PAR
    ! transmittance + reflectance coefficients = green plants + brown plants
    scatp = green * (ltran(1,1) + lref(1,1)) + (1. - green) * &
        (ltran(1,2) + lref(1,2))
    scatg = ltran(1,1) + lref(1,1)

    ! Caclulate PAR absorption optical depth in canopy adjusting for
    ! variance in projected leaf area wrt solar zenith angle
    ! (Sellers et al. Part II (1996), eqn. 13b)
    ! PAR absorption coefficient = (1 - scatp)
    park = sqrt(1. - scatp) * gmudmu

    ! Calculate the new fPAR (Sellers et al. Part II (1996), eqn 9)
    fpar = fvcover * (1. - exp(-park * lai / fvcover))

    ! Ensure calculated fPAR falls within physical limits
    fpar = amax1(fparmin, fpar)
    fpar = amin1(fparmax, fpar)

    return

end subroutine aparnew



!====================================================
subroutine interpolate(x1, x, Dx, y1, y, Dy, z11, z21, z12, z22, z)
!====================================================
! calculates the value of z=f(x,y) by linearly interpolating
! between the 4 closest data points on a uniform grid.  The subroutine
! requires a grid point (x1, y1), the grid spacing (Dx and Dy), and the 
! 4 closest data points (z11, z21, z12, and z22).

use kinds

! begin input variables
real x1  ! the x grid location of z11
real x   ! x-value at which you will interpolate z=f(x,y)
real Dx  ! grid spacing in the x direction
real y1  ! the y grid location of z11
!real y   ! y-value at which you will interpolate z=f(x,y)
real(kind=dbl_kind) :: y   ! y-value at which you will interpolate z=f(x,y)
real Dy  ! grid spacing in the y direction
real z11 ! f(x1, y1)
real z21 ! f(x1+Dx, y1)
real z12 ! f(x1, y1+Dy)
real z22 ! f(x1+Dx, y1+Dy)

! begin output variables
real z   ! f(x,y), the desired interpolated value

! begin internal variables
real zp  ! z'=first interpolated value at (x, y1)
real zpp ! z''=second interpolated value at (x, Y1+Dy)


   ! interpolate between z11 and z21 to calculate z' (zp) at (x, y1)
   zp=z11+(x-x1)*(z21-z11)/Dx

   ! interpolate between z12 and z22 to calculate z'' (zpp) at (x, Y1+Dy)
   zpp=z12+(x-x1)*(z22-z12)/Dx

   ! interpolate between zp and zpp to calculate z at (x,y)
   z=zp+(y-y1)*(zpp-zp)/Dy

   return
   
end subroutine interpolate
