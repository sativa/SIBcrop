subroutine output_control( sib, time, rank )

use kinds
use sibtype
use timetype
use sib_const_module
use sib_io_module

use netcdf

#include "nc_util.h"

implicit none

! parameters
type(sib_t), dimension(subcount), intent(inout) :: sib
type(time_struct), intent(in) :: time
integer(kind=int_kind), intent(in) :: rank

    !itb...output routine here
    if(time%sec_year /= starttime) then
    	call diagnostic_output(sib,                            &
        	qpsib, qp3sib, pbpsib, pbp2sib, iiqpsib, iiqp3sib, &
        	iipbpsib, iipbp2sib, npbpsib, npbp2sib, ijtlensib, &
        	doqpsib, doqp3sib, nqpsib, nqp3sib, indxqp3sib,    &
        	indxqpsib, indxpbpsib, indxpbp2sib,time )
    endif

    ! switch qp file
    if ( time%switch_qp ) then
        call create_qp2( out_path, nqpsib, subcount, ihr, jhr, time%year,    &
                         time%month, longitude, latitude, sublon, sublat,    &
                         doqpsib, nameqpsib, listqpsib, qp2id, qp2varid,     &
                         qp2timeid, qp2charid, qp2startid, qp2endid,         &
                         qp2periodid, drvr_type, biome_source, soil_source,  &
                         soref_source, ndvi_source, c4_source, d13cresp_source,  &
                         rank )
        call create_qp3( out_path, nqp3sib, subcount, ihr, jhr,               &
                         time%year, time%month, nsoil, longitude, latitude,   &
                         sublon, sublat, doqp3sib, nameqp3sib,                &
                         listqp3sib, drvr_type, biome_source, soil_source,    &
                         soref_source, ndvi_source, c4_source, d13cresp_source,  &
                         qp3id, qp3varid, qp3timeid, qp3charid, qp3startid,   &
                         qp3endid, qp3periodid, rank )
    endif


    
    ! switch pbp file
    if ( time%switch_pbp .and. histpp ) then
        call create_pbp( ijtlensib, time%year, time%month, iipbpsib, npbpsib,  &
                         latpbp, lonpbp, dopbpsib, namepbpsib, listpbpsib,  &
                         indxpbpsib, drvr_type, biome_source, soil_source,  &
                         soref_source, ndvi_source, c4_source, d13cresp_source, &
                         out_path, pbptimeid, pbpcharid, pbpvarid,rank )
        call create_pbp2( ijtlensib, nsoil, time%year, time%month, iipbp2sib,  &
                          npbp2sib, latpbp, lonpbp, dopbp2sib,  &
                          namepbp2sib, listpbp2sib, indxpbp2sib, drvr_type,  &
                          biome_source, soil_source, soref_source, ndvi_source,  &
                          c4_source, d13cresp_source, out_path,  &
                          pbp2timeid, pbp2charid, pbp2varid,rank )
    endif

    ! output to pbp
    if ( time%write_pbp .and. histpp ) then
        pbpsib = pbpsib * time%pbp_incnt
        pbp2sib = pbp2sib * time%pbp_incnt
        call write_pbp( ijtlensib, time%year, time%month, time%day, time%sec_year,  &
                        iipbpsib, pbpsib, pbptimeid, pbpcharid,  &
                        pbpvarid, out_path, rank )
        call write_pbp2( ijtlensib, nsoil, time%year, time%month, time%day,  &
                         time%sec_year, iipbp2sib, pbp2sib, pbp2timeid,  &
                         pbp2charid, pbp2varid, out_path, rank )
        pbpsib(:,:) = 0.0
        pbp2sib(:,:,:) = 0.0
    endif

    ! output to qp
    if ( time%write_qp ) then
        qpsib = qpsib * time%qp_incnt
        qp3sib = qp3sib * time%qp_incnt
        call write_qp2( qp2id, qp2timeid, qp2startid, qp2endid, qp2periodid,  &
                        qp2charid, nqpsib, subcount, qp2varid, qpsib,         &
                        doqpsib, indxqpsib, time%year, time%month, time%day,  &
                        time%sec_year, time%end_period, time%period_length )
        call write_qp3( qp3id, qp3timeid, qp3startid, qp3endid, qp3periodid,  &
                        qp3charid, nqp3sib, subcount, nsoil, qp3varid,        &
                        qp3sib, doqp3sib, indxqp3sib, time%year, time%month,  &
                        time%day, time%sec_year, time%end_period,             &
                        time%period_length )
        qpsib(:,:) = 0.0
        qp3sib(:,:,:) = 0.0
    endif


end subroutine output_control


!-------------------------------------------------------------------------------
subroutine file_closer
!-------------------------------------------------------------------------------

  use kinds
  use sib_io_module
  use netcdf
  use typeSizes

  ! local variables
  integer(kind=int_kind) :: status

  ! close all output files
  status = nf90_close( qp2id )
  status = nf90_close( qp3id )
  status = nf90_close( pbpid )
  status = nf90_close( pbp2id )

end subroutine file_closer
