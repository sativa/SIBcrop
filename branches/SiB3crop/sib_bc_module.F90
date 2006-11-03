module sib_bc_module

    !*************************************************************************
    ! This module contains boundary condition variable data that is used to  *
    ! calculate time variant boundary condition data. Values here must be    * 
    ! stored from month to month.                                            *
    !*************************************************************************

    use kinds
    implicit none

    ! Morph Table variables
    type biome_morph_var
        real(kind=real_kind) :: zc               ! Canopy inflection height (m)
        real(kind=real_kind) :: LWidth           ! Leaf width
        real(kind=real_kind) :: LLength          ! Leaf length 
        real(kind=real_kind) :: LAImax           ! Maximum LAI
        real(kind=real_kind) :: stems            ! Stem area index
        real(kind=real_kind) :: NDVImax          ! Maximum NDVI
        real(kind=real_kind) :: NDVImin          ! Minimum NDVI
        real(kind=real_kind) :: SRmax            ! Maximum simple ratio
        real(kind=real_kind) :: SRmin            ! Minimum simple ratio
    end type biome_morph_var

    type(biome_morph_var), allocatable :: MorphTab(:)

    ! Aerodynamic Interpolation table variables
    type aero_var
        real(kind=real_kind) :: zo          ! Canopy roughness coeff 
        real(kind=real_kind) :: zp_disp          ! Zero plane displacement
        real(kind=real_kind) :: RbC              ! RB Coefficient
        real(kind=real_kind) :: RdC             ! RC Coefficient
    end type aero_var

    type(aero_var),allocatable :: AeroVar(:,:,:)

    ! LAIgrid and FVCovergrid are used in mapper
    real(kind=real_kind) :: LAIgrid (50)        !
    real(kind=real_kind) :: fVCovergrid (50)    !

    ! NDVI data
    real(kind=real_kind), dimension(:), allocatable :: prevNDVI

end module sib_bc_module
