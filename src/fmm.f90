module fmm
    use constants
    implicit none
    private
    integer, parameter :: p = 5 !!order of interpolation
    type ch_pol
        real(kind) :: invdenom
        integer :: k
    end type

    real(kind) :: nodes(0:p-1)
    type(ch_pol) :: ch_pols(0:p-1)
    interface
        module subroutine setup_ch_pols()
        end subroutine

        module subroutine compute_chebyshev_nodes()
        end subroutine

        pure elemental module function eval_chebyshev(this, x) result(res)
            type(ch_pol), intent(in) :: this
            real(kind), intent(in) :: x
            real(kind) :: res
        end function
    end interface



    contains

    subroutine setup_ch()
        call compute_chebyshev_nodes()
        call setup_ch_pols()
    end subroutine

end module