!===============================================
subroutine diagnostic_output ( sib, qp2, qp3, pbp1, pbp2, nnqp2,      & 
        nnqp3, npbp1, npbp2, npbp1mx, npbp2mx, ijtlen, doqp2, doqp3,  &
        nnqp2mx, nnqp3mx, indxqp3, indxqp2, indxpbp1, indxpbp2, time)
!===============================================
! Calculates time averages for diagnostic output for single point (pbp) and
!  entire domain (qp)
!
! Modifications:
!  Kevin Schaefer corrected pbp2 from indxpbp1 to indxpbp2 (11/17/04)
!===============================================

    use kinds
    use sib_const_module
    use sibtype
    use timetype
    use physical_parameters, only:    &
           hltm


    use sib_io_module, only: qpintp, histpp, imultpbpsib

    implicit none


    integer(kind=int_kind),intent(in) :: nnqp2
    integer(kind=int_kind),intent(in) :: nnqp3
    integer(kind=int_kind),intent(in) :: npbp1
    integer(kind=int_kind),intent(in) :: npbp2
    integer(kind=int_kind),intent(in) :: ijtlen
    integer(kind=int_kind),intent(in) :: nnqp3mx
    integer(kind=int_kind),intent(in) :: nnqp2mx
    integer(kind=int_kind),intent(in) :: npbp2mx 
    integer(kind=int_kind),intent(in) :: npbp1mx

    logical(kind=log_kind),intent(in) :: doqp2(nnqp2mx) 
    logical(kind=log_kind),intent(in) :: doqp3(nnqp3mx)


    ! time average diagnostic fields (all points)
    real(kind=dbl_kind), intent(inout) :: qp2(subcount,nnqp2) 
    real(kind=dbl_kind), intent(inout) :: qp3(subcount, nsoil, nnqp3)

    ! index of saved diagnostics
    integer(kind=int_kind),intent(in) :: indxqp2(nnqp2mx) 
    integer(kind=int_kind),intent(in) :: indxqp3(nnqp3mx)        
    integer(kind=int_kind),intent(in) :: indxpbp1(npbp1mx)  
    integer(kind=int_kind),intent(in) :: indxpbp2(npbp2mx)


    ! time series diagnostic fields (select points)
    real(kind=dbl_kind), intent(inout) :: pbp1(npbp1+1,ijtlen) 
    real(kind=dbl_kind), intent(inout) :: pbp2(nsoil, npbp2, ijtlen)

    !-------------------------------------------------------------

    type(sib_t), dimension(subcount), intent(inout) :: sib
    type(time_struct), intent(in) :: time

    !-------------------------------------------------------------  

    !...LOCAL VARIABLES...
    integer(kind=int_kind) :: out_index,l,n
    real(kind=dbl_kind) :: discrim3 
    real(kind=dbl_kind) :: auxeadem   

    if(qpintp) then
        do out_index = 1,subcount
        !     save diagnostic output
        !------------------------------------------
        !  depth dependent diagnostics
            if(doqp3(1)) then 
                do l = 1,nsoil
                    qp3(out_index,l,indxqp3(1)) = qp3(out_index,l,indxqp3(1)) + &
                        sib(out_index)%diag%soilscale(l)
                enddo
            endif
            if(doqp3(2)) then 
                do l = 1,nsoil
                    qp3(out_index,l,indxqp3(2)) = qp3(out_index,l,indxqp3(2)) + &
                        sib(out_index)%prog%td(l)
                enddo
            endif
            if(doqp3(3)) then 
                do l = 1,nsoil
                    qp3(out_index,l,indxqp3(3)) = qp3(out_index,l,indxqp3(3)) + &
                        sib(out_index)%prog%www_liq(l) 
                enddo
            endif
            if(doqp3(4)) then 
                do l = 1,nsoil
                    qp3(out_index,l,indxqp3(4)) = qp3(out_index,l,indxqp3(4)) + &
                        sib(out_index)%prog%www_ice(l) 
                enddo
            endif
            if(doqp3(5)) then 
                do l = 1,nsoil
                    qp3(out_index,l,indxqp3(5)) = qp3(out_index,l,indxqp3(5)) + &
                        sib(out_index)%diag%soilq10(l)
                enddo
            endif

        !   depth independent diagnostics
            if(doqp2(1)) then 
                qp2(out_index,indxqp2(1)) = qp2(out_index,indxqp2(1)) +    &
                    sib(out_index)%diag%ventmf 
            endif
            if(doqp2(2)) then 
                qp2(out_index,indxqp2(2)) = qp2(out_index,indxqp2(2)) +    &
                    sib(out_index)%diag%ustar 
            endif
            if(doqp2(3)) then 
                qp2(out_index,indxqp2(3)) = qp2(out_index,indxqp2(3)) +    &
                    sib(out_index)%prog%td(sib(out_index)%prog%nsl+1) 
            endif
            if(doqp2(4)) then 
                qp2(out_index,indxqp2(4)) = qp2(out_index,indxqp2(4)) +    &
                    sib(out_index)%prog%tc 
            endif
            if(doqp2(5)) then 
                qp2(out_index,indxqp2(5)) = qp2(out_index,indxqp2(5)) +    &
                    sib(out_index)%diag%fws*day 
            endif
            if(doqp2(6)) then 
                qp2(out_index,indxqp2(6)) = qp2(out_index,indxqp2(6)) +    &
                    sib(out_index)%prog%snow_depth
            endif
            if(doqp2(7)) then 
                qp2(out_index,indxqp2(7)) = qp2(out_index,indxqp2(7)) +    &
                    sib(out_index)%prog%snow_veg 
            endif
            if(doqp2(8)) then 
                qp2(out_index,indxqp2(8)) = qp2(out_index,indxqp2(8)) +    &
                    sib(out_index)%diag%roffo 
            endif
            if(doqp2(9)) then 
                qp2(out_index,indxqp2(9)) = qp2(out_index,indxqp2(9)) +    &
                    sib(out_index)%diag%ggl(6) 
            endif
            if(doqp2(10)) then 
                if(sib(out_index)%prog%em.lt.1.e-6) then
                    auxeadem = 1.e6
                else
                    auxeadem = sib(out_index)%prog%ea / sib(out_index)%prog%em
                endif
                qp2(out_index,indxqp2(10)) = qp2(out_index,indxqp2(10)) +  &
                    sib(out_index)%diag%ggl(6) * AUXeadem
            endif
            if(doqp2(11)) then 
                qp2(out_index,indxqp2(11)) = qp2(out_index,indxqp2(11)) + &
                    sib(out_index)%diag%ggl(6) * sib(out_index)%diag%rha
            endif
            if(doqp2(12)) then 
                qp2(out_index,indxqp2(12)) = qp2(out_index,indxqp2(12)) + &
                    sib(out_index)%diag%antemp(6) / sib(out_index)%diag%rb
            endif
            if(doqp2(13)) then 
                qp2(out_index,indxqp2(13)) = qp2(out_index,indxqp2(13)) +    &
                    sib(out_index)%diag%assimn(6)
            endif
            if(doqp2(14)) then 
                qp2(out_index,indxqp2(14)) = qp2(out_index,indxqp2(14)) +    &
                    sib(out_index)%diag%ra
            endif
            if(doqp2(15)) then 
                qp2(out_index,indxqp2(15)) = qp2(out_index,indxqp2(15)) +    &
                    sib(out_index)%diag%rb
            endif
            if(doqp2(16)) then 
                qp2(out_index,indxqp2(16)) = qp2(out_index,indxqp2(16)) +    &
                    sib(out_index)%diag%rc
            endif
            if(doqp2(17)) then 
                qp2(out_index,indxqp2(17)) = qp2(out_index,indxqp2(17)) +    &
                    sib(out_index)%diag%rd
            endif
            if(doqp2(18)) then 
                qp2(out_index,indxqp2(18)) = qp2(out_index,indxqp2(18)) +    &
                    sib(out_index)%diag%rsoil
            endif
            if(doqp2(19)) then 
                qp2(out_index,indxqp2(19)) = qp2(out_index,indxqp2(19)) +    &
                    sib(out_index)%diag%rstfac(1)
            endif
            if(doqp2(20)) then 
                qp2(out_index,indxqp2(20)) = qp2(out_index,indxqp2(20)) +    &
                    sib(out_index)%diag%rstfac(2)
            endif
            if(doqp2(21)) then 
                qp2(out_index,indxqp2(21)) = qp2(out_index,indxqp2(21)) +    &
                    sib(out_index)%diag%rstfac(3)
            endif
            if(doqp2(22)) then 
                qp2(out_index,indxqp2(22)) = qp2(out_index,indxqp2(22)) +    &
                    sib(out_index)%diag%rstfac(4)
            endif
            if(doqp2(23)) then 
                qp2(out_index,indxqp2(23)) = qp2(out_index,indxqp2(23)) +    &
                    sib(out_index)%diag%aparkk
            endif
            if(doqp2(24)) then 
                qp2(out_index,indxqp2(24)) = qp2(out_index,indxqp2(24)) +    &
                    sib(out_index)%param%zlt
            endif
            if(doqp2(25)) then 
                qp2(out_index,indxqp2(25)) = qp2(out_index,indxqp2(25)) +    &
                    sib(out_index)%param%green
            endif
            if(doqp2(26)) then 
                qp2(out_index,indxqp2(26)) = qp2(out_index,indxqp2(26)) +    &
                    sib(out_index)%diag%chf
            endif
            if(doqp2(27)) then 
                qp2(out_index,indxqp2(27)) = qp2(out_index,indxqp2(27)) +    &
                    sib(out_index)%diag%shf
            endif
            if(doqp2(28)) then 
                qp2(out_index,indxqp2(28)) = qp2(out_index,indxqp2(28)) +    &
                    sib(out_index)%diag%hc * dti
            endif
            if(doqp2(29)) then 
                qp2(out_index,indxqp2(29)) = qp2(out_index,indxqp2(29)) +    &
                    sib(out_index)%diag%hg * dti
            endif
            if(doqp2(30)) then 
                qp2(out_index,indxqp2(30)) = qp2(out_index,indxqp2(30)) +    &
                    sib(out_index)%diag%ect * dti
            endif
            if(doqp2(31)) then 
                qp2(out_index,indxqp2(31)) = qp2(out_index,indxqp2(31)) +    &
                    sib(out_index)%diag%eci * dti
            endif
            if(doqp2(32)) then 
                qp2(out_index,indxqp2(32)) = qp2(out_index,indxqp2(32)) +    &
                    sib(out_index)%diag%egs * dti
            endif
            if(doqp2(33)) then 
                qp2(out_index,indxqp2(33)) = qp2(out_index,indxqp2(33)) +    &
                    sib(out_index)%diag%egi * dti
            endif
            if(doqp2(34)) then 
                qp2(out_index,indxqp2(34)) = qp2(out_index,indxqp2(34)) +    &
                    sib(out_index)%prog%ea
            endif
            if(doqp2(35)) then 
                qp2(out_index,indxqp2(35)) = qp2(out_index,indxqp2(35)) +    &
                    sib(out_index)%prog%ta
            endif
            if(doqp2(36)) then 
                qp2(out_index,indxqp2(36)) = qp2(out_index,indxqp2(36)) +    &
                    sib(out_index)%prog%em
            endif
            if(doqp2(37)) then 
                qp2(out_index,indxqp2(37)) = qp2(out_index,indxqp2(37)) +    &
                    sib(out_index)%diag%rha
            endif
            if(doqp2(38)) then 
                qp2(out_index,indxqp2(38)) = qp2(out_index,indxqp2(38)) +    &
                    sib(out_index)%diag%omepot(6)
            endif
            if(doqp2(39)) then 
                qp2(out_index,indxqp2(39)) = qp2(out_index,indxqp2(39)) +    &
                    sib(out_index)%diag%assimpot(6)
            endif
            if(doqp2(40)) then 
                qp2(out_index,indxqp2(40)) = qp2(out_index,indxqp2(40)) +    &
                    sib(out_index)%diag%assimnp(6)
            endif
            if(doqp2(41)) then 
                qp2(out_index,indxqp2(41)) = qp2(out_index,indxqp2(41)) +    &
                    sib(out_index)%diag%antemp(6)
            endif
            if(doqp2(42)) then 
                qp2(out_index,indxqp2(42)) = qp2(out_index,indxqp2(42)) +    &
                    sib(out_index)%diag%wsfws(6)
            endif
            if(doqp2(43)) then 
                qp2(out_index,indxqp2(43)) = qp2(out_index,indxqp2(43)) +    &
                    sib(out_index)%diag%wsfht(6)
            endif
            if(doqp2(44)) then 
                qp2(out_index,indxqp2(44)) = qp2(out_index,indxqp2(44)) +    &
                    sib(out_index)%diag%wsflt(6)
            endif
            if(doqp2(45)) then 
                qp2(out_index,indxqp2(45)) = qp2(out_index,indxqp2(45)) +    &
                    sib(out_index)%diag%wags(6)
            endif
            if(doqp2(46)) then 
                qp2(out_index,indxqp2(46)) = qp2(out_index,indxqp2(46)) +    &
                    sib(out_index)%diag%wegs(6)
            endif
            if(doqp2(47)) then 
                qp2(out_index,indxqp2(47)) = qp2(out_index,indxqp2(47)) +    &
                    sib(out_index)%param%aparc
            endif
            if(doqp2(48)) then 
                qp2(out_index,indxqp2(48)) = qp2(out_index,indxqp2(48)) +    &
                    sib(out_index)%diag%assimci(6)
            endif
            if(doqp2(49)) then 
                qp2(out_index,indxqp2(49)) = qp2(out_index,indxqp2(49)) +     &
                    sib(out_index)%diag%wci(6)
            endif
            if(doqp2(50)) then 
                qp2(out_index,indxqp2(50)) = qp2(out_index,indxqp2(50)) +    &
                    sib(out_index)%diag%pfd
            endif
            if(doqp2(51)) then 
                qp2(out_index,indxqp2(51)) = qp2(out_index,indxqp2(51)) +    &
                    (sib(out_index)%diag%ecmass + sib(out_index)%diag%egmass) &
                    * dti * 55.56
            endif
            if(doqp2(52)) then 
                qp2(out_index,indxqp2(52)) = qp2(out_index,indxqp2(52)) +    &
                    sib(out_index)%diag%assim(6)
            endif
            if(doqp2(53)) then 
                qp2(out_index,indxqp2(53)) = qp2(out_index,indxqp2(53)) +    &
                    sib(out_index)%diag%whs(6)
            endif
            if(doqp2(54)) then 
                qp2(out_index,indxqp2(54)) = qp2(out_index,indxqp2(54)) +    &
                    sib(out_index)%prog%capac(1)
            endif
            if(doqp2(55)) then 
                qp2(out_index,indxqp2(55)) = qp2(out_index,indxqp2(55)) +    &
                    sib(out_index)%prog%capac(2)
            endif
            if(doqp2(56)) then 
                qp2(out_index,indxqp2(56)) = qp2(out_index,indxqp2(56)) +    &
                    1./sib(out_index)%prog%rst(6)
            endif
            if(doqp2(57)) then 
                qp2(out_index,indxqp2(57)) = qp2(out_index,indxqp2(57)) +    &
                    sib(out_index)%diag%antemp(6) * sib(out_index)%prog%tc
            endif
            if(doqp2(58)) then 
                qp2(out_index,indxqp2(58)) = qp2(out_index,indxqp2(58)) +    &
                    sib(out_index)%diag%snowmelt
            endif
            if(doqp2(59)) then 
                qp2(out_index,indxqp2(59)) = qp2(out_index,indxqp2(59)) +    &
                    sib(out_index)%diag%ansqr(6)
            endif
            if(doqp2(60)) then 
                do l = 1,nsoil
                    qp2(out_index,indxqp2(60)) = qp2(out_index,indxqp2(60)) +    &
                        sib(out_index)%prog%www_liq(l) * 0.001 
                enddo
            endif
            if(doqp2(61)) then 
                qp2(out_index,indxqp2(61)) = qp2(out_index,indxqp2(61)) +    &
                    sib(out_index)%diag%fss
            endif
            if(doqp2(62)) then
                qp2(out_index,indxqp2(62)) = qp2(out_index,indxqp2(62)) +    &
                    sib(out_index)%prog%tm
            endif
            if(doqp2(63)) then
                qp2(out_index,indxqp2(63)) = qp2(out_index,indxqp2(63)) +    &
                    sib(out_index)%prog%thm
            endif
            if(doqp2(64)) then
                qp2(out_index,indxqp2(64)) = qp2(out_index,indxqp2(64)) +    &
                    sib(out_index)%prog%sh
            endif
            if(doqp2(65)) then
                qp2(out_index,indxqp2(65)) = qp2(out_index,indxqp2(65)) +    &
                    sib(out_index)%prog%radvbc
            endif
            if(doqp2(66)) then
                qp2(out_index,indxqp2(66)) = qp2(out_index,indxqp2(66)) +    &
                    sib(out_index)%prog%radvdc
            endif
            if(doqp2(67)) then
                qp2(out_index,indxqp2(67)) = qp2(out_index,indxqp2(67)) +    &
                    sib(out_index)%prog%radnbc
            endif
            if(doqp2(68)) then
                qp2(out_index,indxqp2(68)) = qp2(out_index,indxqp2(68)) +    &
                    sib(out_index)%prog%radndc
            endif
            if(doqp2(69)) then
                qp2(out_index,indxqp2(69)) = qp2(out_index,indxqp2(69)) +    &
                    sib(out_index)%prog%dlwbot
            endif
            if(doqp2(70)) then
                qp2(out_index,indxqp2(70)) = qp2(out_index,indxqp2(70)) +    &
                    sib(out_index)%prog%spdm
            endif
            if(doqp2(71)) then
                qp2(out_index,indxqp2(71)) = qp2(out_index,indxqp2(71)) +    &
                    sib(out_index)%prog%ps
            endif
            if(doqp2(72)) then
                qp2(out_index,indxqp2(72)) = qp2(out_index,indxqp2(72)) +    &
                    sib(out_index)%prog%lspr*3600.0
            endif
            if(doqp2(73)) then
                qp2(out_index,indxqp2(73)) = qp2(out_index,indxqp2(73)) +    &
                    sib(out_index)%prog%cupr*3600.0
            endif
            if(doqp2(74)) then 
                qp2(out_index,indxqp2(74)) = qp2(out_index,indxqp2(74)) +    &
                    sib(out_index)%diag%radc3(1)
            endif
            if(doqp2(75)) then 
                qp2(out_index,indxqp2(75)) = qp2(out_index,indxqp2(75)) +    &
                    sib(out_index)%diag%radc3(2)
            endif
            if(doqp2(76)) then
                discrim3 = -5.0 + 30.0 * sib(out_index)%diag%pco2i(6)   &
                    /sib(out_index)%prog%pco2ap
                qp2(out_index,indxqp2(76)) = qp2(out_index,indxqp2(76)) +    &
                    discrim3
            endif
            if(doqp2(77))then
                qp2(out_index,indxqp2(77)) = qp2(out_index,indxqp2(77)) +    &
                    sib(out_index)%prog%pco2ap
            endif
            if(doqp2(78))then
                qp2(out_index,indxqp2(78)) = qp2(out_index,indxqp2(78)) +    &
                    sib(out_index)%diag%pco2c(6)
            endif
            if(doqp2(79))then
                qp2(out_index,indxqp2(79)) = qp2(out_index,indxqp2(79)) +    &
                    sib(out_index)%diag%pco2i(6)
            endif
            if(doqp2(80))then
                qp2(out_index,indxqp2(80)) = qp2(out_index,indxqp2(80)) +    &
                    sib(out_index)%diag%pco2s(6)
            endif
            if(doqp2(81))then
                qp2(out_index,indxqp2(81)) = qp2(out_index,indxqp2(81)) +    &
                    sib(out_index)%prog%d13cca
            endif
            if(doqp2(82))then
                qp2(out_index,indxqp2(82)) = qp2(out_index,indxqp2(82)) +    &
                    sib(out_index)%prog%d13cm
            endif
            if(doqp2(83))then
                qp2(out_index,indxqp2(83)) = qp2(out_index,indxqp2(83)) +    &
                    sib(out_index)%param%d13cresp
            endif
            if(doqp2(84))then
                qp2(out_index,indxqp2(84)) = qp2(out_index,indxqp2(84)) +    &
                    sib(out_index)%diag%kiecps(1) !C3 plants
            endif
            if(doqp2(85))then
                qp2(out_index,indxqp2(85)) = qp2(out_index,indxqp2(85)) +    &
                    sib(out_index)%diag%kiecps(2) !C4 plants
            endif
            if(doqp2(86))then
                qp2(out_index,indxqp2(86)) = qp2(out_index,indxqp2(86)) +    &
                    sib(out_index)%diag%d13cassimn(1) !C3 plants
            endif
            if(doqp2(87))then
                qp2(out_index,indxqp2(87)) = qp2(out_index,indxqp2(87)) +    &
                    sib(out_index)%diag%d13cassimn(2) !C4 plants
            endif
            if(doqp2(88))then
                qp2(out_index,indxqp2(88)) = qp2(out_index,indxqp2(88)) +    &
                    sib(out_index)%diag%d13cassimn(6) !All plants summed
            endif
            if(doqp2(89))then
                qp2(out_index,indxqp2(89)) = qp2(out_index,indxqp2(89)) +    &
                    sib(out_index)%diag%flux13c
            endif

            if(doqp2(90))then
                qp2(out_index,indxqp2(90)) = qp2(out_index,indxqp2(90)) +    &
                    sib(out_index)%diag%flux12c
            endif

            if(doqp2(91)) then
                qp2(out_index,indxqp2(91)) = qp2(out_index,indxqp2(91)) +    &
                    sib(out_index)%diag%flux_turb
            endif
            if(doqp2(92)) then
                qp2(out_index,indxqp2(92)) = qp2(out_index,indxqp2(92)) +    &
                    sib(out_index)%diag%respg
            endif
            if(doqp2(93)) then
                qp2(out_index,indxqp2(93)) = qp2(out_index,indxqp2(93)) +    &
                    sib(out_index)%prog%sw_dwn
            endif
            if(doqp2(94)) then
                qp2(out_index,indxqp2(94)) = qp2(out_index,indxqp2(94)) +    &
                    sib(out_index)%stat%coszbar
            endif
            if(doqp2(95)) then
                qp2(out_index,indxqp2(95)) = qp2(out_index,indxqp2(95)) +    &
                    sib(out_index)%stat%cosz
            endif
            if(doqp2(96)) then
                qp2(out_index,indxqp2(96)) = qp2(out_index,indxqp2(96)) +    &
                    sib(out_index)%param%physfrac(2)
            endif
            if(doqp2(97)) then
                qp2(out_index,indxqp2(97)) = qp2(out_index,indxqp2(97)) +    &
                    sib(out_index)%diag%rcassimn(1) * sib(out_index)%diag%assimn(1)
            endif
            if(doqp2(98)) then
                qp2(out_index,indxqp2(98)) = qp2(out_index,indxqp2(98)) +    &
                    sib(out_index)%diag%rcassimn(2) * sib(out_index)%diag%assimn(2)
            endif
            if(doqp2(99)) then
                qp2(out_index,indxqp2(99)) = qp2(out_index,indxqp2(99)) +    &
                    sib(out_index)%diag%rcassimn(6) * sib(out_index)%diag%assimn(6)
            endif
            if(doqp2(100)) then
                qp2(out_index,indxqp2(100)) = qp2(out_index,indxqp2(100)) +    &
                    sib(out_index)%diag%assimn(1)
            endif
            if(doqp2(101)) then
                qp2(out_index,indxqp2(101)) = qp2(out_index,indxqp2(101)) +    &
                    sib(out_index)%diag%assimn(2)
            endif
            if(doqp2(102)) then
                qp2(out_index,indxqp2(102)) = qp2(out_index,indxqp2(102)) +    &
                    sib(out_index)%diag%assimn(6)*(sib(out_index)%diag%eastar-sib(out_index)%prog%ea)
            endif
            if(doqp2(103)) then
                qp2(out_index,indxqp2(103)) = qp2(out_index,indxqp2(103)) +    &
                    sib(out_index)%diag%antemp(1)*sib(out_index)%diag%kiecps(1) 
            endif
            if(doqp2(104)) then
                qp2(out_index,indxqp2(104)) = qp2(out_index,indxqp2(104)) +    &
                    sib(out_index)%diag%antemp(2)*sib(out_index)%diag%kiecps(2)
            endif
            if(doqp2(105)) then 
                qp2(out_index,indxqp2(105)) = qp2(out_index,indxqp2(105)) +    &
                    sib(out_index)%diag%antemp(1) 
            endif
            if(doqp2(106)) then 
                qp2(out_index,indxqp2(106)) = qp2(out_index,indxqp2(106)) +    &
                    sib(out_index)%diag%antemp(2)
            endif
            if(doqp2(107)) then
                qp2(out_index,indxqp2(107)) = qp2(out_index,indxqp2(107)) +    &
                    (sib(out_index)%diag%respg - sib(out_index)%diag%assimn(6))*1.0E6
            endif
            if(doqp2(108)) then
                qp2(out_index,indxqp2(108)) = qp2(out_index,indxqp2(108)) +    &
                    sib(out_index)%diag%cflux*1.0E6
            endif
            if(doqp2(109)) then
                qp2(out_index,indxqp2(109)) = qp2(out_index,indxqp2(109)) +    &
                    sib(out_index)%diag%fws
            endif
            if(doqp2(110)) then
                qp2(out_index,indxqp2(110)) = qp2(out_index,indxqp2(110)) +    &
                    sib(out_index)%diag%www_tot_soil
            endif
 
           
        enddo ! index
    endif

    if(histpp) then
        do n = 1,ijtlen
        !   depth dependent diagnostics
            do l = 1,nsoil
                pbp2(l,indxpbp2(1),n) = pbp2(l,indxpbp2(1),n) +    &
                    sib(imultpbpsib(n))%prog%www_liq(l)/(sib(imultpbpsib(n))%prog%dz(l) * &
                                sib(imultpbpsib(n))%param%poros * denh2o)
            enddo
            do l = 1,nsoil
                pbp2(l,indxpbp2(2),n) = pbp2(l,indxpbp2(2),n) +    &
                    sib(imultpbpsib(n))%prog%td(l)
            enddo
            do l = 1,nsoil
                pbp2(l,indxpbp2(3),n) = pbp2(l,indxpbp2(3),n) +    &
                    sib(imultpbpsib(n))%param%csolid(l)
            enddo
            do l = 1,nsoil
                pbp2(l,indxpbp2(4),n) = pbp2(l,indxpbp2(4),n) +    &
                    sib(imultpbpsib(n))%prog%node_z(l)
            enddo

            do l = 1,nsoil
                    pbp2(l,indxpbp2(5),n) = pbp2(l,indxpbp2(5),n) +    &
                    sib(imultpbpsib(n))%diag%soilscale(l)
            enddo
	    
        !   depth independent diagnostics
            out_index=1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%fss

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%fws

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%rst(6)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%assimn(6)*1.0E6

            out_index=out_index+1
            pbp1(indxpbp1(5),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%respg*1.E6

            out_index=out_index+1
            pbp1(indxpbp1(6),n) = pbp1(indxpbp1(out_index),n) +    &
                (sib(imultpbpsib(n))%diag%respg - &
                sib(imultpbpsib(n))%diag%assimn(6))*1.E6

            out_index=out_index+1
            pbp1(indxpbp1(7),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%cflux * 1.E6

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%td(sib(imultpbpsib(n))%prog%nsl+1)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%tc

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%ea

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%ta

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%em

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%rha

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%capac(1)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%capac(2)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%snow_veg

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) -    &
                sib(imultpbpsib(n))%prog%nsl

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%areas

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%snow_mass

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%snow_depth

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%pco2ap

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%pco2c(6)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%pco2i(6)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%pco2s(6)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%roff

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%qqq 

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%roffo

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%ra

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%rb

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%rc

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%rd

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%rsoil

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%rstfac(1)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%rstfac(2)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%rstfac(3)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%rstfac(4)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +   &
                sib(imultpbpsib(n))%param%scatp

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +   &
                sib(imultpbpsib(n))%param%scatg

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +   &
                sib(imultpbpsib(n))%param%park

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +   &
                sib(imultpbpsib(n))%param%gmudmu

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%aparkk

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%param%zlt

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%param%green

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%chf

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%shf

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%hc * dti

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%hg * dti

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%ect * dti

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%eci * dti

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%egs * dti

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%egi * dti

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%ess   

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%omepot(6)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%assimpot(6)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%assimnp(6) 

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%antemp(6) 

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%wsfws(6) 

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%wsfht(6)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    & 
                sib(imultpbpsib(n))%diag%wsflt(6)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%wags(6)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%wegs(6)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%param%aparc

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%assimci(6)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%wci(6)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%pfd

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) + &
                (sib(imultpbpsib(n))%diag%ecmass +    &
                sib(imultpbpsib(n))%diag%egmass) * dti * 55.56

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%assim(6)*1.e6

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%whs(6)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%cu

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%ct

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%ventmf

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%tm

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%thm

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%sh

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%radvbc

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%radvdc

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%radnbc

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%radndc

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%dlwbot

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%spdm

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%ps

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%lspr*3600.0

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%cupr*3600.0

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%radc3(1)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%radc3(2)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%stat%cosz

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%d13cca

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%d13cm

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%param%d13cresp

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%kiecps(1)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%kiecps(2)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%d13cassimn(1)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%d13cassimn(2)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%d13cassimn(6)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%flux13c

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    & 
                sib(imultpbpsib(n))%diag%flux12c

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%flux_turb

!itb...soil/snow temperatures (loop)

            out_index = out_index + 1 !index at: 98
            do l=-nsnow+1,nsoil

                pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                    sib(imultpbpsib(n))%prog%td(l)
                out_index = out_index + 1
            enddo


!itb...soil/snow liquid water (loop)

            do l=-nsnow+1,nsoil
                pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                    sib(imultpbpsib(n))%prog%www_liq(l)

                out_index = out_index + 1
            enddo


!itb...soil/snow ice water (loop)

            do l=-nsnow+1,nsoil
                pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                    sib(imultpbpsib(n))%prog%www_ice(l)

                out_index = out_index + 1
            enddo

          pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%www_tot_soil

!itb...volumetric soil water content (loop)

            out_index = out_index + 1 !index at: 144
            do l=1,nsoil

                pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) + &
                    sib(imultpbpsib(n))%prog%www_liq(l)/   &
                    (sib(imultpbpsib(n))%prog%dz(l)*denh2o)

                out_index = out_index + 1
            enddo


!itb...phystype-specific stomatal resistance (loop)

            !index at: 154
            do l=1,5

                pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                    sib(imultpbpsib(n))%prog%rst(l)

                out_index = out_index + 1
            enddo


!itb...phystype-specific net assimilation (loop)

             do l=1,5 
                pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                    sib(imultpbpsib(n))%diag%assimn(l)*1.0E6

                out_index = out_index + 1
            enddo


!itb...phystype-specific chloroplast CO2 partial pressure (loop)

            do l =1,5
                pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                    sib(imultpbpsib(n))%diag%pco2c(l)

                out_index = out_index + 1
            enddo


!itb...phystype-specific leaf internal CO2 partial pressure (loop)

            do l=1,5

                pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                    sib(imultpbpsib(n))%diag%pco2i(l)

                out_index = out_index + 1
            enddo

!itb...phystype-specific leaf surface CO2 partial pressure (loop)

            do l=1,5
                pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                    sib(imultpbpsib(n))%diag%pco2s(l)

                out_index = out_index + 1
            enddo

          pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%pco2m

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%stat%coszbar

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%sw_dwn2

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%prog%sw_dwn

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%radt(1)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%radt(2)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%radt(3)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%cas_e_storage

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%radfac(1,1,1)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%radfac(1,2,1)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%radfac(1,1,2)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%radfac(1,2,2)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%radfac(2,1,1)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%radfac(2,2,1)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%radfac(2,1,2)

            out_index=out_index+1
            pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%radfac(2,2,2)

            out_index=out_index+1 
            do l =1,6
                pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                    sib(imultpbpsib(n))%diag%respc(l)*1.E6
            
                out_index = out_index + 1
            enddo

          pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
               (sib(imultpbpsib(n))%diag%respg + &
               sib(imultpbpsib(n))%diag%aparkk* sib(imultpbpsib(n))%diag%respc(6))*1.E6

           out_index=out_index+1
          pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) + &
               sib(imultpbpsib(n))%param%clayfrac
           
            out_index=out_index+1
          pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) + &
               sib(imultpbpsib(n))%param%sandfrac

            out_index=out_index+1
         pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) + &
               sib(imultpbpsib(n))%param%biome

            out_index=out_index+1
          pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) + &
                sib(imultpbpsib(n))%param%vcover

          !Crop Variables (modified by kdcorbin, 02/11)
            out_index=out_index+1
          pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%gdd

            out_index=out_index+1
          pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%ta_bar

            out_index=out_index+1
          pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%rstfac_d

            out_index=out_index+1
         pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%assim_d

            out_index=out_index+1
          pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%w_main

            out_index=out_index+1
         pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%leafwt_c

            out_index=out_index+1
          pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%phen_lai

            out_index=out_index+1
            do l =1,4
                pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                    sib(imultpbpsib(n))%diag%alloc(l)
            
                out_index = out_index + 1
            enddo

         !index at: 218
            do l =1,4
                pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                    sib(imultpbpsib(n))%diag%cum_wt(l)
            
                out_index = out_index + 1
            enddo

            do l =1,4
                pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                    sib(imultpbpsib(n))%diag%cum_drywt(l)
            
                out_index = out_index + 1
            enddo

         !index at: 226
            do l =1,4
                pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                    sib(imultpbpsib(n))%diag%phen_maintr(l)
            
                out_index = out_index + 1
            enddo

            do l =1,4
                pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                    sib(imultpbpsib(n))%diag%phen_growthr(l)
            
                out_index = out_index + 1
            enddo
          pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%pd

            out_index=out_index+1
          pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%emerg_d

            out_index=out_index+1
          pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%ndf_opt

            out_index=out_index+1
          pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%phen_switch

            out_index=out_index+1
         pbp1(indxpbp1(out_index),n) = pbp1(indxpbp1(out_index),n) +    &
                sib(imultpbpsib(n))%diag%tot_biomass

        enddo
    endif

end subroutine diagnostic_output
