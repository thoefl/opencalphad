    module mod_th

    use liboceqplus
    implicit none
    contains

    pure function int_to_char(i) result(c)
        integer, intent(in) :: i
        character(:), allocatable :: c

        allocate(character(128) :: c)
        write(c, '(i0)') i
        c = trim(c)
    end function int_to_char

    pure subroutine errcheck_allocate(caller_name, err_msg, stat)
        character(*), intent(in) :: caller_name, err_msg
        integer, intent(in) :: stat
        if(stat /= 0) then
            error stop 'ERROR in procedure ' // caller_name // ': allocation returned error code '&
                // int_to_char(stat) // ': ' // err_msg
        end if
    end subroutine errcheck_allocate

    pure function linspace(startval, stopval, nval)
        real(dp), intent(in) :: startval, stopval
        integer, intent(in) :: nval
        real(dp) :: linspace(nval)
        integer :: i
    
        linspace = [(startval + real(i,dp)/(nval - 1)*(stopval - startval), i=0, nval-1)]
    end function linspace
    
    !/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

    !\addtotable subroutine equilph1d_th
    !\begin{verbatim}
    subroutine equilph1d_th(phase_tuple, temp, pressure, x_known, mu, n_endmemb, dmu, mobility, current_equi, silent_opt)
    ! equilibrates the constituent fractions of a phase for mole fractions x_known
    ! and calculates the Darken matrix and unreduced diffusivities
    ! phase_tuple is phase tuple
    ! tpval is T and P
    ! current_equi is a datastructure with all relevant thermodynamic data
    ! mu are the (calculated) chemical potentials
    ! silent is TRUE means no output
    ! n_endmemb is the number of values returned in dmu
    ! dmu are the derivatives of the chemical potentials wrt mole fractions??
    ! mobility are the mobilities
    type(gtp_phasetuple), pointer, intent(inout) :: phase_tuple
    real(dp), intent(in) :: temp, pressure
    real(dp), dimension(:), intent(in) :: x_known
    real(dp), dimension(:), intent(out) :: dmu, mobility, mu
    integer, intent(out) :: n_endmemb
    type(gtp_equilibrium_data), pointer, intent(inout) :: current_equi
    logical, intent(in), optional :: silent_opt

    logical silent
    !\end{verbatim} %+
    type(meq_setup) :: meqrec
    integer i

    if(present(silent_opt)) then
        silent = silent_opt
    else
        silent = .false.
    end if
    ! extract the current chemical potentials as start values
    do i = 1, noel()
        mu(i) = current_equi%cmuval(i)
    enddo
    if(gx%bmperr.ne.0) goto 1000
    ! create the meqrec structure
    call equilph1_meqrec(phase_tuple, meqrec, silent, current_equi)
    if(gx%bmperr.ne.0) goto 1000
    ! maybe we need RT ?
    current_equi%rtn = globaldata%rgas*temp
    ! iterate until equilibrium found for this phase
    call equilph1e_th(meqrec, meqrec%phr, temp, pressure, x_known, mu, silent, n_endmemb, dmu, mobility, current_equi)
1000 continue
    return
    end subroutine equilph1d_th

    !/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

    !\addtotable subroutine equilph1e_th
    !\begin{verbatim} %-
    subroutine equilph1e_th(meqrec, phr, temp, pressure, x_known, mu, silent, n_endmemb, dmu, mobility, current_equi)
    ! iterate constituent fractions of a phase for mole fractions x_known
    ! and calculate derivatives of MU and diffusion coefficients
    ! tpval is T and P
    ! x_known are mole fractions
    ! nrel is the number of components (elements)
    ! mu are the chemical potentials
    ! silent is TRUE if no output
    ! dmu is the derivatives of the chemical potentials wrt mole fractions??
    ! mobility are the mobilities
    ! current_equi is a datastructure with all relevant thermodynamic data

    type(meq_setup), intent(inout) :: meqrec
    type(meq_phase), dimension(*), target, intent(inout) :: phr
    real(dp), intent(in) :: temp, pressure
    real(dp), dimension(:), intent(in) :: x_known
    real(dp), dimension(:), intent(out) :: mu, mobility, dmu
    logical, intent(in) :: silent
    integer, intent(out) :: n_endmemb
    type(gtp_equilibrium_data), pointer, intent(inout) :: current_equi
    !\end{verbatim}
    real(dp) :: tpval(2), endmemb_stride
    integer :: i, nz1, nz2, converged, ierr, jj, nj, constit_idx, nl, is, jt, sublat_idx, component_idx, elem_idx
    integer :: lokph, n_constit(maxsubl), first(maxsubl+1), current(maxsubl), n_sublat, endmemb_idx
    integer :: deriv(maxsubl), ql, eq_idx, loksp
    type(meq_phase), pointer :: pmi
    real(dp), allocatable :: smat(:,:), svar(:), yarr(:), delta(:), coeff_mat(:,:)
    integer, allocatable, dimension(:) :: ipiv, selected_endmemb_idx
    ! dmu_endmemb_dy is derivatives of mu for endmembers wrt all constituents
    real(dp), dimension(:,:), allocatable :: dmu_endmemb_dy, dy_endmemb_dn, dmu_endmemb_dn,&
        dmu_constit_dy, dy_constit_dn, dmu_constit_dn
    real(dp), dimension(:), allocatable :: py, mu_endmemb, mu_constit, n_sites
    real(dp) :: chargefact, chargerr, pv, qq(5), ys, ycormax2, muall, sumsum, constit_mass, species_mass, species_charge
    integer :: mqindex, allocate_stat, lapack_stat, n_vacancy, n_eq, n_species_elem, n_extra, n_elem
    integer, dimension(maxel) :: species_elem_idx
    real(dp), dimension(maxel) :: stoich, extra
    character(256) :: allocate_errmsg
    character(24) :: component_name, constit_name
    character(24), allocatable, dimension(:,:) :: endmemb_constit_name
    character(*), parameter :: proc_name = 'equilph1e_th'

    tpval = [temp, pressure]
    ! mqindex is a constant set in gtpini in models/gtp3A.F90
    ! number of variables is number of components + one stable phase
    nz1 = meqrec%nrel + 1
    nz2 = nz1 + 1
    mobility = 0

    allocate(smat(nz1,nz2), stat=allocate_stat, errmsg=allocate_errmsg)
    call errcheck_allocate(proc_name, allocate_errmsg, allocate_stat)
    allocate(svar(nz1), stat=allocate_stat, errmsg=allocate_errmsg)
    call errcheck_allocate(proc_name, allocate_errmsg, allocate_stat)
    allocate(delta(phr(1)%ncc), stat=allocate_stat)
    call errcheck_allocate(proc_name, allocate_errmsg, allocate_stat)
    allocate(yarr(phr(1)%ncc), stat=allocate_stat)
    call errcheck_allocate(proc_name, allocate_errmsg, allocate_stat)

    chargefact = one
    chargerr = one
    ! we have just one phase in phr, phr must be TARGET
    pmi => phr(1)
    do
        converged=0
        smat = zero
        ! invert the phase matrix for pmi
        call meq_onephase(meqrec, pmi, current_equi)
        if(gx%bmperr.ne.0) goto 1000
        ! setup mass balance equations, note some components may be missing
        ! This is a simplified setup_equilmatrix using x_known as composition
        call setup_comp2cons(meqrec, phr, nz1, smat, tpval, x_known, converged, current_equi)
        if(gx%bmperr.ne.0) goto 1000
        ! solve the equilibrium matrix, some chemical potentials may be missing
        call lingld(nz1, nz2, smat, svar, nz1, ierr)
        if(ierr.ne.0) then
            write(*,*)'Error solving equilibrium matrix 3',ierr
            gx%bmperr=4203; goto 1000
        endif
        ! check that svar(1..meqrec%nrel) has converged
        do jj = 1, meqrec%nrel
            if(abs(svar(jj)-mu(jj)) > 1.0E1_dp*current_equi%xconv) then
                converged = 7
                cerr%mconverged = converged
                if(cerr%nvs.lt.10) then
                    cerr%nvs=cerr%nvs+1
                    cerr%typ(cerr%nvs)=7
                    cerr%val(cerr%nvs)=svar(jj)
                    cerr%dif(cerr%nvs)=mu(jj)
                endif
            endif
            mu(jj)=svar(jj)
        enddo
        !    write(*,111)'svar3: ',0,(svar(jj),jj=1,nz1)
111     format(a,i2,6(1pe12.4))
        ! check dxmol ... seems OK
        !    do constit_idx=1,phr(1)%ncc
        !       write(*,111)'dxmol: ',constit_idx,(phr(1)%dxmol(nl,constit_idx),nl=1,meqrec%nrel)
        !    enddo
        ! update constituent fractions in just one phase
        !    lap: do jj=1
        jj=1
        ! The current chemical potentials are in current_equi%cmuval(i) svar(1..n)
        ! jj is stable, increment kk but do not make it larger than meqrec%nstph
        ! save index in meqrec%stphl in jph !!!!!!!!!!! kk never used !!!!!!!!!
        !    jph=kk
        !    kk=min(kk+1,meqrec%nstph)
        ! if phr(jj)%xdone=1 then phase has no composition variation
        if(phr(jj)%xdone.eq.1) goto 1000
        !----------------------------------------------------
        ycormax2=zero
        !    write(*,*)'cc: ',jj
        ! loop for all constituents
        !    write(*,112)'Y0: ',meqrec%noofits,converged,(yarr(nj),nj=1,phr(jj)%ncc)
        moody: do nj=1,phr(jj)%ncc
            ys=zero
            do constit_idx=1,phr(jj)%ncc
                pv=zero
                do nl=1,meqrec%nrel
                    ! current_equi%cmuval(nl) is the chemical potential of element nl (divided by RT)
                    ! When a chemical potential is fixed use meqrec%mufixval
                    ! phr(jj)%dxmol(nl,constit_idx) is the derivative of component nl
                    ! wrt constituent constit_idx
                    !?             pv=pv+current_equi%complist(nl)%chempot(1)/current_equi%rtn*phr(jj)%dxmol(nl,constit_idx)
                    !             pv=pv+current_equi%cmuval(nl)*phr(jj)%dxmol(nl,constit_idx)
                    !             pv=pv+svar(nl)*phr(jj)%dxmol(nl,constit_idx)
                    pv=pv+mu(nl)*phr(jj)%dxmol(nl,constit_idx)
                enddo
                pv=pv-phr(jj)%curd%dgval(1,constit_idx,1)
                ys=ys+phr(jj)%invmat(nj,constit_idx)*pv
                !          write(*,111)'pv: ',nj,ys,pv,phr(1)%curd%dgval(1,constit_idx,1),&
                !               phr(1)%invmat(nj,constit_idx)
            enddo
            if(phr(jj)%chargebal.eq.1) then
                ! For charged phases add a term
                ! phr(jj)%invmat(phr(jj)%idim,phr(jj)%idim)*Q
                ys=ys-chargefact*phr(jj)%invmat(nj,phr(jj)%idim)*&
                    phr(jj)%curd%netcharge
            endif
            delta(nj)=ys
            if(abs(delta(nj)).gt.ycormax2) then
                ycormax2=delta(nj)
            endif
            if(abs(ys).gt.current_equi%xconv) then
                ! if the change in any constituent fraction larger than xconv continue iterate
                if(converged.lt.4) then
                    ! large correction in fraction of constituent fraction of stable phase
                    !             write(*,*)'mm converged 4C: ',jj,nj,ys
                    converged=4
                    !             yss=ys
                    !             yst=phr(jj)%curd%yfr(nj)
                endif
                !       elseif(phr(jj)%stable.eq.1) then
                ! check to find good convergence criteria in Re-V test case
                !          if(abs(delta(nj)).gt.ysmm) then
                !             jmaxy=jj
                !             ysmm=abs(delta(nj))
                !             ysmt=phr(jj)%curd%yfr(nj)
                !           endif
            endif
            yarr(nj)=phr(jj)%curd%yfr(nj)+delta(nj)
        enddo moody
        ! >>>>>>>>>>>>>>>>>> HERE the new constitution is set <<<<<<<<<<<<<<<<<<<<<
        !    write(*,112)'YC: ',jj,(delta(nj),nj=1,phr(jj)%ncc)
        !    write(*,112)'YY: ',meqrec%noofits,converged,(yarr(nj),nj=1,phr(jj)%ncc)
112     format(a,2i3,8F8.5)
        !    write(*,*)'MM calling set_constitution 5:',phr(jj)%iph,phr(jj)%ics
        call set_constitution(phr(jj)%iph,phr(jj)%ics,yarr,qq,current_equi)
        if(gx%bmperr.ne.0) goto 1000
        !-------------------------- end of iteration
        ! check convergence
        meqrec%noofits=meqrec%noofits+1
        if(converged.gt.3) then
            if(meqrec%noofits.le.current_equi%maxiter) cycle
            gx%bmperr=4204
            !       write(*,*)'MM Too many iterations',current_equi%maxiter
            goto 1000
        elseif(meqrec%noofits.lt.6) then
            cycle
        else
            if(.not.btest(meqrec%status, MMQUIET)) write(*,202) meqrec%noofits
202         format('Calculation required ', i4, ' its')
            exit
        endif
    end do

    do is=1,meqrec%nrel
        mu(is)=svar(is)
    enddo

    !----------------------------------------------------------
    ! When the calculation converged we calculate dmu and interdiffusivites
    ! A nontrival expression:
    !
    ! dmu_i/dx_j = 1/N (d2G/dx_i/dx_j - \sum_k x_k (d2G/dx_k/dx_i + d2G/dx_k/dx_j)+
    !                              \sum_k\sum_m x_k x_m d2G/dx_k/dx_m )
    !
    ! NOTE THIS IS SYMMETRICAL, dmu_i/dx_j = dmu_j/dx_i.
    ! If the phase is ideal then d2G/dx_i/dx_j = RT/x_i if i=j, otherwise zero
    ! This gives for
    ! dmu_i/dx_i = RT/N * (1-x_i)/x_i
    ! dmu_i/dx_j = - RT/N                  (i not equal to j)
    !
    ! We calc             sum_k (x_k*d2G/dx_k/dx_i)   in delta(i)
    !         sum_m x_m ( sum_k (x_k*d2G/dx_k/dx_m))  in sumsum
    !
    ! new use of delta !!!
    delta=zero
    muall=pmi%curd%gval(1,1)
    sumsum=zero
    ! Here we calculate delta(is) =           \sum_jt y(jt)*d2G/dy_jt/dy_is and
    !                   sumsum = \sum_m y(is) \sum_jt y(jt)*d2G/dy_jt/dy_is
    ! The loop of is is for all constituents
    do is=1,phr(1)%ncc
        ! The loop for jt are for all constituents in all sublattices
        do jt=1,phr(1)%ncc
            ! STRANGE that d2G/dy_Va/dy_Va is zero ... should be 1 (*RT) ...does not matter
            !          if(is.gt.jt) stop "wrong order 1"
            ! keep ixsym here as I do not know if jt>is or not
            delta(is)=delta(is)+pmi%curd%yfr(jt)*pmi%curd%d2gval(ixsym(is,jt),1)
            !          write(*,*)'d2G/dy/dy: ',is,jt,pmi%curd%d2gval(ixsym(is,jt),1)
        enddo
        sumsum=sumsum+pmi%curd%yfr(is)*delta(is)
        muall=muall-pmi%curd%dgval(1,is,1)*pmi%curd%yfr(is)
    enddo
    ! muall    = G_m - \sum_i y_i dG/dy_i
    ! delta(i) = \sum_j y_j d2G/dy_i/dy_j             sum for all y_j for one y_i
    ! sumsum   = \sum_i \sum_j y_i y_j d2G/dy_i/dy_j  sum for all y_i and y_j
    !-------------------- summations over all constituents in all sublattices
    ! now we must generate the endmembers, loop over all sublattices
    ! but sublattics and number of constituents in each are in the phase record
    ! and protected ... use a subroutine ...
    lokph=pmi%curd%phlink
    call get_phase_structure(lokph,n_sublat,n_constit)
    if(gx%bmperr.ne.0) goto 1000
    ! ---------------------------------------------------------------------
    substitutional: if(n_sublat.eq.1) then
        ! specially simple if n_sublat=1 (substitutional)
        n_endmemb=n_constit(1)
        allocate(mu_endmemb(n_endmemb), stat=allocate_stat, errmsg=allocate_errmsg)
        call errcheck_allocate(proc_name, allocate_errmsg, allocate_stat)
        ! calculate just mu(endmember)
        !       loop1: do endmemb_idx=1,n_endmemb
        !          mu_endmemb(endmemb_idx)=muall+pmi%curd%dgval(1,endmemb_idx,1)
        !          loop2: do jt=1,n_endmemb
        ! the chemical potential has the derivative of the constituent
        !             mu_endmemb(endmemb_idx)=mu_endmemb(endmemb_idx)+pmi%curd%dgval(1,jt,1)
        !          enddo loop2
        !       enddo loop1
        ! now we calculate dmu(end)/dy_is (just for substitutional)
        allocate(dmu_endmemb_dy(n_endmemb,pmi%ncc), stat=allocate_stat, errmsg=allocate_errmsg)
        call errcheck_allocate(proc_name, allocate_errmsg, allocate_stat)
        dmu_endmemb_dy=zero
        ! For a substitutional solution:
        ! dmu_i/dx_j = 1/N ( d2G/dx_i/dx_j -
        !                    \sum_k x_k d2G/dx_k/dx_i - \sum_k x_k d2G/dx_k/dx_j+
        !                    \sum_k\sum_m x_k x_m d2G/dx_k/dx_m )
        ! NOTE THIS SHOULD BE SYMMETRICAL, dmu_i/dx_j = dmu_j/dx_i.
        ! use delta(i) and sumsum calculated above
        !       write(*,*)'Derivatives of chemical potentials',n_endmemb
        nl=0
        loop3: do is=1,n_endmemb
            mu_endmemb(is)=muall+pmi%curd%dgval(1,is,1)
            loop4: do jt=1,n_endmemb
                !             if(is.gt.jt) stop "wrong order 2"
                ! keep using ixsym here as I do not know if jt>is
                dmu_endmemb_dy(is,jt)=pmi%curd%d2gval(ixsym(is,jt),1)-&
                    delta(is)-delta(jt)+sumsum
                !             write(*,775)'dd1:',1,is,jt,&
                !                  dmu_endmemb_dy(is,jt),pmi%curd%d2gval(ixsym(is,jt),1),&
                !                  delta(is),delta(jt),sumsum
                nl=nl+1
                dmu(nl)=dmu_endmemb_dy(is,jt)*current_equi%rtn
            enddo loop4
            !          write(*,777)'dd: ',(current_equi%rtn*dmu_endmemb_dy(is,jt),jt=1,n_endmemb)
            !777       format(a,6(1pe12.4))
        enddo loop3
        ! UNFINISHED ?? I do not divide by N
        !       write(*,777)'mu: ',(mu_endmemb(is),is=1,n_endmemb)
        !-------------------
    else ! not substitutional below (2 or more sublattices)
        ! now we have to handle sublattices and endmembers
        ! n_sublat is number of sublattices and n_constit(1..n_sublat) the number of const in each
        n_endmemb=1
        is=1
        first=0
        do nl=1,n_sublat
            ! endmemb_idx is number of endmembers
            ! here first and current are set to first constituent index in each sublattice
            n_endmemb=n_endmemb*n_constit(nl)
            first(nl)=is
            current(nl)=is
            deriv(nl)=is
            is=is+n_constit(nl)
        enddo
        ! we need this to indicate when we reached the end
        first(n_sublat+1)=is
        allocate(mu_endmemb(n_endmemb), stat=allocate_stat, errmsg=allocate_errmsg)
        call errcheck_allocate(proc_name, allocate_errmsg, allocate_stat)
        allocate(py(n_endmemb), stat=allocate_stat, errmsg=allocate_errmsg)
        call errcheck_allocate(proc_name, allocate_errmsg, allocate_stat)
        py=one
        !       write(*,611)'first: ',n_sublat,(first(nj),nj=1,n_sublat)
        !611    format(a,i2,2x,10i3)
        ! all partials have this term
        mu_endmemb=muall
        !       write(*,*)'MM muall: ',muall,pmi%curd%gval(1,1)
        ! The partial Gibbs energy, for each sublattice add one dG/dy_is
        endmemb_idx=0
        nj=0
        allpg: do while(nj.le.n_sublat)
            endmemb_idx=endmemb_idx+1
            ! the partials constituents, G_I, are in current(1..n_sublat)
            nlloop: do nl=1,n_sublat
                is=current(nl)
                ! endmembers like 1:1:1, 1:1:2, 1:2:1, 1:2:2, 2:1:1, 2:1:2, 2:2:1, 2:2:2 =8
                ! constituents are in current(1..n_sublat)
                mu_endmemb(endmemb_idx)=mu_endmemb(endmemb_idx)+pmi%curd%dgval(1,is,1)
            enddo nlloop
            ! generate a new set of constituents in current
            nj=1
888         continue
            current(nj)=current(nj)+1
            if(current(nj).eq.first(nj+1)) then
                ! note first(n_sublat+1) is the end of all constituents
                current(nj)=first(nj)
                nj=nj+1
                if(nj.le.n_sublat) goto 888
            endif
        enddo allpg
        if(.not.silent) then
            write(*,881)(mu_endmemb(jt),jt=1,n_endmemb)
881         format('Calculated potentials for all endmembers/RT: '/6(1x,1pe12.4))
        endif
        !-----------------------------------------------------------------------
        ! FIXME_TH: DOCUMENTATION GOES HERE
        ! the part below is messy and unfinished
        !---------------- now the derivative of the partial Gibbs energy
        ! The partial Gibbs energy, for each sublattice add one dG/dy_is
        ! the derivative of the partial Gibbs energy wrt all other endmembers ....
        ! dG_i/dn_J = 1/N_J( \sum_s (d2G/dy_is/dy_js - delta(is) - delta(js)) + sumsum )
        ! delta(is) = \sum_s \sum_k y_k d2G/dy_is/dy_k
        ! sumsum    = \sum_k \sum_m y_k y_m d2G/dy_k/dy_m (already added above)
        !---------------------------------------------------
        n_elem = size(pmi%xmol)
        endmemb_stride = (endmemb_idx - 1.0_dp)/(n_elem - 1)
        selected_endmemb_idx = nint(linspace(1.0_dp, real(n_endmemb, dp), n_elem))
        !allocate(endmemb_constit_name(n_sublat, n_endmemb), stat=allocate_stat, errmsg=allocate_errmsg)
        !call errcheck_allocate(proc_name, allocate_errmsg, allocate_stat)
        allocate(dmu_endmemb_dy(n_elem, pmi%ncc), stat=allocate_stat, errmsg=allocate_errmsg)
        call errcheck_allocate(proc_name, allocate_errmsg, allocate_stat)
        allocate(dy_endmemb_dn(n_elem, pmi%ncc), stat=allocate_stat, errmsg=allocate_errmsg)
        call errcheck_allocate(proc_name, allocate_errmsg, allocate_stat)
        allocate(n_sites(n_sublat), stat=allocate_stat, errmsg=allocate_errmsg)
        call errcheck_allocate(proc_name, allocate_errmsg, allocate_stat)

        dmu_endmemb_dy = 0
        dy_endmemb_dn = 0
        delta = 0
        do is=1,phr(1)%ncc
            do jt=1,phr(1)%ncc
                delta(is) = delta(is) + pmi%curd%yfr(jt)*pmi%curd%d2gval(ixsym(is,jt),1)
            enddo
        enddo
        
        call get_sublattice_structure(pmi%iph, pmi%ics, n_sublat, n_constit, n_sites, current_equi)
        
        allocate(coeff_mat(n_elem, n_elem), stat=allocate_stat, errmsg=allocate_errmsg)
        call errcheck_allocate(proc_name, allocate_errmsg, allocate_stat)
        allocate(ipiv(n_elem), stat=allocate_stat, errmsg=allocate_errmsg)
        call errcheck_allocate(proc_name, allocate_errmsg, allocate_stat)
        
        coeff_mat = 0
        
        !subroutine get_constituent_location(lokph,cno,loksp)
        ! returns the location of the species record of a constituent
        ! requred for ionic liquids as phlista is private
        
         !subroutine get_species_component_data(loksp,nspel,compnos,stoi,smass,qsp,ceq)
        ! return species data, loksp is from a call to find_species_record
        ! Here we return stoichiometry using components 
        ! nspel: integer, number of components in species
        ! compno: integer array, component (species) indices
        ! stoi: double array, stoichiometric factors
        ! smass: double, mass of species
        ! qsp: double, charge of the species
         
        !subroutine get_species_data(loksp,nspel,ielno,stoi,smass,qsp,nextra,extra)
        ! return species data, loksp is from a call to find_species_record
        ! nspel: integer, number of elements in species
        ! ielno: integer array, element indices
        ! stoi: double array, stoichiometric factors
        ! smass: double, mass of species
        ! qsp: double, charge of the species
        ! nextra, integer, number of additional values
        ! extra: double, some additional values like UNIQUAC volume and area
        
        ! loop for all partial Gibbs energies G_I
        do endmemb_idx = 1, n_endmemb
            eq_idx = findloc(selected_endmemb_idx, endmemb_idx, 1)
            ! loop for all sublattices
            if(eq_idx > 0) then
                mu(eq_idx) = mu_endmemb(endmemb_idx)
                dmu_endmemb_dy(eq_idx,:) = -delta
                dy_endmemb_dn(eq_idx,:) = -pmi%curd%yfr
                do sublat_idx = 1, n_sublat
                    is = current(sublat_idx)
                        call get_phase_record(pmi%iph, lokph)
                        call get_constituent_location(lokph, is, loksp)
                        ! Skip if the species is a vacancy, as defined in gtp3A.F90
                        if(loksp > 1) then
                            call get_species_data(&
                                loksp, n_species_elem, species_elem_idx, stoich, species_mass, species_charge, n_extra, extra)
                            do concurrent(i = 1:n_species_elem)
                                elem_idx = species_elem_idx(i)
                                coeff_mat(eq_idx,elem_idx) = n_sites(sublat_idx)*stoich(i)
                            end do
                        end if
                    !endmemb_constit_name(sublat_idx, endmemb_idx) = constit_name
                    ! loop for all constituents in all sublattices
                    do constit_idx = 1, phr(1)%ncc
                        dmu_endmemb_dy(eq_idx,constit_idx) =&
                            dmu_endmemb_dy(eq_idx,constit_idx) + pmi%curd%d2gval(ixsym(is,constit_idx),1)
                    end do
                    dy_endmemb_dn(eq_idx,is) = dy_endmemb_dn(eq_idx,is) + 1
                end do
            end if
            
            ! set constituent indices in sublattices for next endmember
            do sublat_idx = 1, n_sublat
                current(sublat_idx) = current(sublat_idx) + 1
                if(current(sublat_idx) == first(sublat_idx + 1)) then
                    ! note first(n_sublat) is the end of all constituents
                    current(sublat_idx) = first(sublat_idx)
                else
                    exit
                end if
            end do
        end do
        
        call dgetrf(n_elem, n_elem, coeff_mat, n_elem, ipiv, lapack_stat)
        call dgetrs('N', n_elem, pmi%ncc, coeff_mat, n_elem, ipiv, dmu_endmemb_dy, n_elem, lapack_stat) 
        call dgetrs('N', n_elem, pmi%ncc, coeff_mat, n_elem, ipiv, dy_endmemb_dn, n_elem, lapack_stat)
        call dgetrs('N', n_elem, 1, coeff_mat, n_elem, ipiv, mu, n_elem, lapack_stat)
        dmu = reshape(matmul(dmu_endmemb_dy, transpose(dy_endmemb_dn)), [n_elem**2])
        dmu = dmu*globaldata%rgas*temp
        mu = mu*globaldata%rgas*temp
        
        ! Calculate component chemical potentials and partial Gibbs energy derivatives from endmembers, if possible
        ! This currently only works for phases containing all components selected from the database
        
    endif substitutional
1000 continue
    return
    end subroutine equilph1e_th

    !/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\!/!\

    !\begin{verbatim}
      subroutine tqgdmat_th(phase_tuple_idx, temp, pressure, x_known, mu, dmu, current_equi, silent_opt)
    ! equilibrates the constituent fractions of a phase for mole fractions xknown
    ! and calculates the Darken matrix and unreduced diffusivities
    ! phtup is phase tuple
    ! tpval is T and P
    ! ceq is a datastructure with all relevant thermodynamic data
    ! cpot are the (calculated) chemical potentials
    ! tyst is TRUE means no outut
    ! nend is the number of values returned in mugrad
    ! mugrad are the derivatives of the chemical potentials wrt mole fractions??
    ! mobval are the mobilities
                    !  type(gtp_phasetuple), pointer, intent(inout) :: phase_tuple
                    !real(dp), intent(in) :: temp, pressure
                    !real(dp), dimension(:), intent(in) :: x_known
                    !real(dp), dimension(:), intent(out) :: dmu, mobility, mu
                    !integer, intent(out) :: n_endmemb
                    !type(gtp_equilibrium_data), pointer, intent(inout) :: current_equi
                    !logical, intent(in), optional :: silent_opt
        integer, intent(in) :: phase_tuple_idx
        real(dp), intent(in) :: temp, pressure
        real(dp), dimension(:), intent(in) :: x_known
        real(dp), dimension(:), intent(out) :: dmu, mu
        type(gtp_equilibrium_data), pointer, intent(inout) :: current_equi
        logical, intent(in), optional :: silent_opt
        logical :: silent
        type(gtp_phasetuple), pointer :: phase_tuple
        integer :: n_endmemb
        real(dp), dimension(:), allocatable :: mobility
    !\end{verbatim}

        if(present(silent_opt)) then
            silent = silent_opt
        else
            silent = .false.
        end if
             
        phase_tuple => phasetuple(phase_tuple_idx)
        allocate(mobility(0))
        call equilph1d_th(phase_tuple, temp, pressure, x_known, mu, n_endmemb, dmu, mobility, current_equi, silent)
      end subroutine tqgdmat_th

    !\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\!/!!\
    
    end module mod_th