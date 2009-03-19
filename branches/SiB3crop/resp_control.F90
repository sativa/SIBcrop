subroutine respfactor_control( sib, time, rank )

use kinds
use timetype
use sib_const_module
use sib_io_module
use sibtype
use physical_parameters, only: tice  
implicit none

! parameters
type(sib_t), dimension(subcount), intent(inout) :: sib
type(time_struct), intent(in) :: time
integer(kind=int_kind), intent(in) :: rank

! local variables
integer(kind=int_kind) :: i, s, m
character(len=250) :: filename
real(kind=dbl_kind), dimension(:,:), allocatable :: resp

    ! add new values to sums of assimn and soilscale
    do i = 1, subcount
        sib(i)%diag%tot_an(time%month) = sib(i)%diag%tot_an(time%month) +  &
            sib(i)%diag%assimn(6) * time%dtsib
        do s = 1, nsoil
            sib(i)%diag%tot_ss(time%month,s) =  &
                sib(i)%diag%tot_ss(time%month,s) +   &
                sib(i)%diag%soilscale(s) * time%dtsib
        enddo
    enddo

!EL....added new info from crop_accum.F90
!    do i = 1, subcount
!       sib(i)%diag%tot_biomass(time%month) = sib(i)%diag%tot_biomass(time%month) +  &
!            sib(i)%diag%tot_biomass(time%doy)
 !      sib(i)%diag%prodwt(time%month) = sib(i)%diag%prodwt(time%month) +  &
 !           sib(i)%diag%prodwt(time%doy)
!    enddo


    ! if new year, calculate new respfactor
    if ( time%calc_respf ) then
        do i = 1, subcount
            sib(i)%diag%tot_an(13) = 0.0_dbl_kind
            sib(i)%diag%tot_ss(13,:) = 0.0_dbl_kind
!            sib(i)%diag%tot_BM_an(13)=0.0_dbl_kind
!            sib(i)%diag%tot_prod_an(13)=0.0_dbl_kind
        enddo
        do m = 1, 12
            do i = 1, subcount
                sib(i)%diag%tot_an(13) = sib(i)%diag%tot_an(13) +  &
                    sib(i)%diag%tot_an(m)
                do s = 1, nsoil
                    sib(i)%diag%tot_ss(13,s) = sib(i)%diag%tot_ss(13,s) +  &
                        sib(i)%diag%tot_ss(m,s)
                enddo
            enddo
!EL....added new info from crop_accum.F90
!            do i = 1, subcount 
!              sib(i)%diag%tot_BM_an(13)= sib(i)%diag%tot_BM_an +  &
!                    sib(i)%diag%tot_biomass(m)  
 !             sib(i)%diag%tot_prod_an(13)= sib(i)%diag%tot_prod_an +  &
!                    sib(i)%diag%prodwt(m)
!            enddo
        enddo
        


        call respire( sib )
        
        ! write out respfactor file
        if ( time%write_respf ) then
            write( filename, '(a,a,i4.4,a,i3.3,a)' ) trim(out_path),  &
                'CO2_respf_', time%year, 'p', rank, '_new'
            if ( drvr_type == 'single' ) then
                open( unit=37, file=trim(filename), status='unknown',  &
                    form='formatted')
                do s = 1, nsoil
                    write( 37,'(e10.5,i4,a)' ) sib(1)%param%respfactor(s), s,  &
                        '    respfactor, level'
                enddo
                close( 37 )
            else
                allocate( resp(nsib,nsoil) )
                resp(:,:) = 0.0_dbl_kind
                do i = 1, subcount
                    resp(subset(i),:) = sib(i)%param%respfactor(:)
                enddo
                open( unit=37, file=trim(filename), status='unknown',  &
                    form='unformatted')
                write( 37 ) nsib
                write( 37 ) nsoil
                write( 37 ) resp
                close( 37 )
                deallocate( resp )
            endif
        endif
        
    endif
    
end subroutine respfactor_control

!===============================================================================
!===============================================================================

subroutine respire( sib,time )

use kinds
use sibtype
use sib_const_module
use timetype
implicit none

!  calculate the annual respiration rate "respfactor" for each of 7
!   soil layer at each grid cell in the model, given monthly mean
!   maps of net carbon assimilation and "soilscale" at each level.

!  references:
!  denning et al., 1996a, tellus 48b, 521-542
!  denning at al., 1996b, tellus 48b, 543-567 

!  soilscale is the product of a temperature response and a moisture
!   response function, evaluated in each soil layer. it is also called
!   r* in denning et al (1996a), equations 6 through 9, and in
!   denning et al (1996b), equations 8 and 9

!  note: the soilscale qp3 from bugs has only 6 layers, so 
!  *** before calling respire, you must "fill in" the bottom (first)
!      layer of soilscale with zeros ****


! parameters
type(sib_t), dimension(subcount), intent(inout) :: sib
type(time_struct), intent(in) :: time                                
! local variables

integer :: n,l,n1                  ! looping indices
double precision :: xag         ! fraction of annual respiration from
                                !   above-ground sources
double precision :: xagmin      ! minimum above-ground resp fraction 
                                !   (marginal ecosystems, low gpp)
double precision :: xagmax      ! maximum above-ground resp fraction 
                                !   (robust ecosystems, high gpp)  
double precision :: anainflec   ! tot_an at inflection point for xag function
double precision :: kxag        ! exponential constant for xag function 
double precision :: xbg         ! fraction of annual respiration from
                                !   below-ground sources    
double precision :: roota       ! rooting distribution parameters
double precision :: rootb
real :: temp1,tempc_sib
parameter(xagmin = 0.10, xagmax = 0.75, anainflec = 1000.,  &
    kxag=5.e-3, roota = 5.0, rootb = 1.5)

    do n = 1, subcount

        ! apportion annual respiration flux due to above-ground and
        ! below-ground sources of organic matter at each grid point 

        ! above-ground fraction

!print*,'tot_an=',sib(n)%diag%tot_an(13),'anainflec=',anainflec, sib(n)%param%rootf(1),sib(n)%param%rootf(2)
     


        xag = xagmin + (xagmax - xagmin) /  &
           ( 1.0 + exp( -kxag * ( sib(n)%diag%tot_an(13) * 12. - anainflec ) ) )

        ! below-ground fraction

        xbg = 1.0 - xag

!for crops, multiply tot_an (below) by 0.78 to represent that 40% was removed in harvest..
        ! vertical distribution of respiration flux in the root zone   


        ! top two layers
        sib(n)%param%respfactor(1) = sib(n)%diag%tot_an(13) * 0.78 *( 0.5 * xag  +  &
            xbg * sib(n)%param%rootf(1))

        sib(n)%param%respfactor(2) = sib(n)%diag%tot_an(13) * 0.78 *( 0.5 * xag  +  &
            xbg * sib(n)%param%rootf(2))

        ! rooting layers (3..10)
        do l = 3, nsoil
            sib(n)%param%respfactor(l) = sib(n)%diag%tot_an(13)* 0.78 * xbg  &
                * sib(n)%param%rootf(l)
!print*,'rootf_l=3',sib(n)%param%rootf(l)
        end do

        !  divide by annual soilscale sum at each grid pt and layer


        do l = 1, nsoil
            sib(n)%param%respfactor(l) = sib(n)%param%respfactor(l) /  &
                sib(n)%diag%tot_ss(13,l)
!print*,'tot_ss',sib(n)%diag%tot_ss(13,l)
        end do ! (next layer)
 

    end do ! (next grid cell)
       

end subroutine respire
