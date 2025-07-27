submodule(fmm) cheb
    contains


    module subroutine setup_ch_pols()
        integer ::k
        do k = 0, p-1
            ch_pols(k)%k = k
            call compute_chebyshev_poly_denom(ch_pols(k))
        end do
    end subroutine

    module subroutine compute_chebyshev_nodes()  !!Computes the chebyshev nodes
        integer :: k
        do k = 0,p-1
            nodes(k) = cos((k+0.5)*pi/p) 
        end do
    end subroutine

    subroutine compute_chebyshev_poly_denom(this)
        type(ch_pol) :: this
        integer :: m
        real(kind) ::prod
        prod = 1
        do m = 0, p-1
            if(this%k /= m) then
                prod = prod * (nodes(this%k)- nodes(m))
            endif
        end do
        this%invdenom=1.0_kind/prod
    end subroutine

    pure elemental real(kind) module function eval_chebyshev(this, x) result(res)
        type(ch_pol), intent(in) :: this
        real(kind), intent(in) :: x
        integer :: m
        real(kind) :: prod
        prod = 1
        do m = 0, p-1
            if(m/=this%k) then
                prod = prod * (x - nodes(m))
            endif
        end do
        prod = prod * this%invdenom
        res = prod
    end function

end submodule 