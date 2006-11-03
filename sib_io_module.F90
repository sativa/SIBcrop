module sib_io_module

    use kinds
    implicit none

    ! SiB diagnostics
    ! allocatable arrays

    integer(kind=int_kind) :: driver_id     ! driver data file id#

    real(kind=dbl_kind), dimension(:,:),   allocatable ::    &
        qpsib,   & ! horizontal only, time mean global diagnostics
        pbpsib     ! time series diagnostics at selected points, no vertical

    real(kind=dbl_kind), dimension(:,:,:), allocatable ::    &
        qp3sib,     & ! 3d time mean global diagnostics
        pbp2sib       ! time-height series diagnostics at selected points, 

    ! netcdf variable id #s
    integer(kind=int_kind) :: &
        qp2id,                &
        qp2timeid,            &
        qp2charid,            &
        qp2startid,           &
        qp2endid,             &
        qp2periodid,          &
        qp3id,                &
        qp3timeid,            &
        qp3charid,            &
        qp3startid,           &
        qp3endid,             &
        qp3periodid,          &
        pbpid,                &
        pbptimeid,            &
        pbpcharid,            &
        pbp2id,               &
        pbp2timeid,           &
        pbp2charid
        
    integer(kind=int_kind), dimension(:), allocatable ::  &
        qp2varid,       &
        qp3varid,       &
        pbpvarid,       &
        pbp2varid

    integer, dimension(:), allocatable ::    &
        indxqpsib,   & ! indices of saved qpsib fields
        indxqp3sib,  & ! indices of saved qp3sib fields
        indxpbpsib,  & ! indices of saved pbpsib fields
        indxpbp2sib, & ! indices of saved pbp2sib fields
        numqpsib,    & ! index number of qpsib field (from sib_qpopts)
        numqp3sib,   & ! index number of qp3sib field (from sib_qpopts)
        numpbpsib,   & ! index number of pbpsib field (from sib_pbpopts)
        numpbp2sib,  & ! index number of pbp2sib field (from sib_pbpopts)
        imultpbpsib    ! gridpoint indices where pbp fields saved

    character(len=80), dimension(:), allocatable ::    &
        listqpsib,  & ! descriptions of qpsib fields
        listpbpsib, & ! descriptions of qp3sib fields
        listqp3sib, & ! descriptions of pbpsib fields
        listpbp2sib   ! descriptions of pbp2sib fields

    character(len=16), dimension(:), allocatable ::    &
        nameqpsib,  & ! field names of qp fields
        nameqp3sib, & ! field names of qp3 fields
        namepbpsib, & ! field names of pbp fields
        namepbp2sib  ! field names of pbp2 fields

    logical, dimension(:), allocatable ::    &
        doqpsib,    & ! true for qpsib fields to be saved
        dopbpsib,   & ! true for qp3sib fields to be saved
        doqp3sib,   & ! true for pbpsib fields to be saved
        dopbp2sib     ! true for pbp2sib fields to be saved


    integer  ijtlensib,  & ! number of gridpoints where pbp fields are saved
        nqpsib,     & ! number of possible qpsib fields (from sib_qpopts)
        nqp3sib,    & ! number of possible qp3sib fields (from sib_qpopts)
        npbpsib,    & ! number of possible pbpsib fields (from sib_pbpopts)
        npbp2sib,   & ! number of possible pbp2sib fields (from sib_pbpopts)
        iiqpsib,    & ! number of saved qp fields
        iipbpsib,   & ! number of saved pbp fields
        iiqp3sib,   & ! number of saved qp3 fields
        iipbp2sib     ! number of saved pbp2 fields

    ! file path names read in from namel file
    character (len=256) ::     &
        param_path,   & !jk path to sib parameters file
        ic_path,      & !jk path to initial conditions file
        dr_format,    & !jk path to driver data in FORMAT form
                        !jk with provision to write *yyyymm*
                        !jk filenames.  Binary files need additional
                        !jk alpha format: *aaa_yyyymm* to specify data type
        out_path,     & !jk path for output files
        qp_path,      & !jk path to qpopts file
        pbp_path,     & !jk path to pbpopts file
        co2_path,     & !jk path to CO2_resp file
        grid_path,    & !jk path to gridmap file
        c4_path         !itb path to c4 fraction file-monthly bc

    ! jk flag to indicate driver data type
    character (len=8) :: drvr_type 
                        !jk 'ecmwf'   - ECMWF global 
                        !jk 'single'  - single point ASCII table
                        !KS 'ncep1'   - ncep1 global
                        !   'ncep2'   - ncep2 global
                        !   'geos4'   - geos4 global

    ! parameter attributes
    character (len=100) :: biome_source
    character (len=100) :: soil_source
    character (len=100) :: soref_source
    character (len=100) :: ndvi_source
    character (len=100) :: c4_source
    character (len=100) :: d13cresp_source

    !itb...LOGICAL OPTIONS
    logical (kind=log_kind)  ::   &
        hswtch = .false.      ,& !  flag to call switch routines
        qpintp = .true.       ,& !  flag to save qp info
        histpp = .false.      ,& !  flag to save pbp info
        debug = .false.       ,& !  run in debug mode?
        roll_respf               !  calculate a rolling respfactor?


end module sib_io_module
