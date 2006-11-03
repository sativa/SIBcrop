module kinds

    !-----------------------------------------------------------------------
    !
    !     This module defines variable precision for all common data
    !     types.
    !
    !-----------------------------------------------------------------------

    implicit none
    save

    !-----------------------------------------------------------------------

    integer, parameter :: &
        char_len  = 80,                    &
        int_kind  = kind(1),               &
        long_kind = selected_int_kind(18), &
        log_kind  = kind(.true.),          &
        real_kind = selected_real_kind(6), &
        dbl_kind  = selected_real_kind(13)

end module kinds
