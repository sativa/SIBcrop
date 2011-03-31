subroutine read_ndvi(filename, sib, time)
!--------------------------------------------------------------------
! opens current sib_bc file and checks if nsib matches, then reads data
! kdcorbin, 02/11 - modified input variables to include sib and time variables

use kinds
use sibtype
use timetype
use sib_const_module, only: nsib,     & ! total number of SiB points
                            subset,   & ! array mapping subdomain to nsib vector
                            subcount, & ! actual number of SiB points being simulated
                            physmax, & ! maximum # of physiology types
                            nper    ! number of ndvi composite periods - kdcorbin, 02/11
use sib_io_module, only: drvr_type,  & !jk flag to indicate driver data type
                         c4_source,  &
                         ndvi_source,  &
                         d13cresp_source, &
                         param_type !flag to indicate vegetation type
use netcdf
use typeSizes

#include "nc_util.h"

! define input variables
character *100, intent(in) :: filename ! Filename for sib_bc
type(sib_t), dimension(subcount), intent(inout) :: sib
type(time_struct) :: time

! define local variables
integer(kind=int_kind) :: i,j
integer(kind=int_kind) :: ntest1     ! value of nsib read from sib_bc file
integer(kind=int_kind) :: dimid
integer(kind=int_kind) :: status
integer(kind=int_kind) :: begin (2)  ! indices where to start and 
integer(kind=int_kind) :: finish (2) !   finish reading from file
integer(kind=int_kind) :: yyid       ! ndvi file id#
integer(kind=int_kind) :: ndvi_id    ! ndvi variable id#
integer(kind=int_kind) :: lai_id       ! LAI variable id# - kdcorbin, 02/11
integer(kind=int_kind) :: fpar_id     ! fPAR variable id# - kdcorbin, 02/11
integer(kind=int_kind) :: modis_time_id  ! modis time variable id# - kdcorbin, 02/11
integer(kind=int_kind) :: var_id       ! misc variable id #
integer(kind=int_kind) :: phys_id    ! physfrac variable id#
integer(kind=int_kind) :: d13_id     ! d13cresp variable id#
integer(kind=int_kind) :: x,y
real(kind=real_kind), dimension(nsib) :: ndvi ! temp variables for parameters
real(kind=real_kind), dimension(nsib,physmax)  :: frac
real(kind=real_kind), dimension(nsib)  :: d13
character(len=10) :: name

!kdcorbin, 02/11
real(kind=real_kind), dimension(:), allocatable :: mstart,mstop
real(kind=real_kind), dimension(nsib) :: lai     !temp lai variable
real(kind=real_kind), dimension(nsib) :: fpar   !temp fpar variable
real(kind=real_kind), dimension(nsib) :: modis_time   !temp modis time variable

    ! Open the sib_bc file
    CHECK( nf90_open( trim(filename), nf90_nowrite, yyid ) )
    CHECK( nf90_inq_dimid ( yyid, 'nsib', dimid ) )
    CHECK( nf90_inquire_dimension ( yyid, dimid, name, ntest1 ) )
    ! test if boundary condition files are correct
    if(ntest1 /= nsib) stop ' open: file sib_bc no match with model for nsib'

    !d13
    ENSURE_VAR( yyid, 'd13cresp', d13_id)
    CHECK( nf90_get_var ( yyid, d13_id, d13) )    

    !physfrac
    ENSURE_VAR( yyid, 'physfrac', phys_id)
    CHECK( nf90_get_var ( yyid, phys_id, frac) )

    ! read in global attributes that are passed on to output files
    CHECK( nf90_get_att( yyid, nf90_global, 'ndvi_source', ndvi_source ) )
    CHECK( nf90_get_att( yyid, nf90_global, 'c4_source', c4_source ) )
    CHECK( nf90_get_att( yyid, nf90_global, 'd13cresp_source', d13cresp_source ) )

    ! set the variables to sib values
    do i=1,subcount
         sib(i)%param%d13cresp = d13(subset(i))
         do j=1,physmax
              sib(i)%param%physfrac(j) = frac(subset(i),j)
         enddo
     enddo

    !kdcorbin, 02/11 - read in either NDVI or LAI/fPAR data:
    print*,'      param_type= ',param_type
    if ( param_type == 'ndvi' ) then
       ENSURE_VAR( yyid, 'ndvi', ndvi_id )
       ENSURE_VAR( yyid, 'd13cresp', d13_id )
       ENSURE_VAR( yyid, 'physfrac', phys_id )

       print*,'Not Done Yet!!!  Stopping.'
       stop

       ! read in one month's data for all nsib points
       !begin (1)  = time%ppmonth
       !begin (2)  = 1 
       !finish (1) = 1 
       !finish (2) = nsib 
       !CHECK( nf90_get_var ( yyid, ndvi_id, ndvi, begin, finish ) )
       !CHECK( nf90_get_var ( yyid, d13_id, d13 ) )
       !CHECK( nf90_get_var ( yyid, phys_id, frac ) )

       ! copy only points in subdomain
       !do i = 1, subcount
       !    ndvidata(i) = ndvi(subset(i))
       !    d13cresp(i) = d13(subset(i))
       !    do j = 1, physmax
       !        physfrac(i,j) = frac(subset(i),j)
       !    enddo
       !enddo

   else
       ENSURE_VAR( yyid, 'lai', lai_id )
       ENSURE_VAR( yyid, 'fpar', fpar_id)
       ENSURE_VAR( yyid, 'modis_time', modis_time_id)

       !get the number of composite periods, start and stop times
       ENSURE_VAR( yyid, 'mapsyear', var_id)
       CHECK( nf90_get_var ( yyid, var_id, nper ) )

       allocate(mstart(nper))
       ENSURE_VAR( yyid, 'modis_start', var_id)
       CHECK( nf90_get_var ( yyid, var_id, mstart ) )
       time%modis_start(1:nper) = mstart

       allocate(mstop(nper))
       ENSURE_VAR( yyid, 'modis_stop', var_id)
       CHECK( nf90_get_var ( yyid, var_id, mstop ) )
       time%modis_stop(1:nper) = mstop

       do i=1,nper
            if (time%doy >= time%modis_start(i)) then
                if (time%doy <= time%modis_stop(i)) time%bc_recnum=i
                if (i==nper .and. time%doy < time%modis_stop(i)+365.) time%bc_recnum=i
             endif
       enddo

       begin(1) = time%bc_recnum
       begin(2) = 1
       finish(1) = 1
       finish(2) = nsib

       CHECK( nf90_get_var ( yyid, lai_id, lai, begin, finish ) )
       CHECK( nf90_get_var ( yyid, fpar_id, fpar, begin, finish) )
       CHECK( nf90_get_var ( yyid, modis_time_id, modis_time, begin, finish) )

       !copy points, only using non-crop points
       do i=1,nsib
            if (sib(i)%param%biome < 20) then
                sib(i)%param%lai = lai(i)
                sib(i)%param%fpar = fpar(i)
                sib(i)%param%modis_time = modis_time(i)
                sib(i)%param%modis_period = time%bc_recnum
            endif   
       enddo
    endif  !ndvi or lai/fpar 

   ! close file
   CHECK( nf90_close( yyid ) )

end subroutine read_ndvi
