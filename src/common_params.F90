module common_params
    
! Single precision real numbers
integer, parameter :: single_prec = selected_real_kind(6)

! Double precision real numbers
integer, parameter :: double_prec = selected_real_kind(15)
    
end module