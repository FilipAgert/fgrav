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


    type cluster
        integer :: N
        real(kind), allocatable :: poles(:,:), weights(:)
    end type

    type sph_harm_coeff !!Data structure for storing spherical harmonic coefficients c_n^m
    !since each n stores (2n+1) coefficients, a 2d array is inefficient. 1d array is better.
        real(kind), private, allocatable :: data(:)
        integer, private, allocatable :: m_ptr(:)
        contains
            procedure :: alloc => alloc_sph_coeff
            procedure :: get => get_sph_coeff
            procedure :: set => set_sph_coeff
            procedure, private :: get_m_ptr
    end type
    contains
    real(kind) function  get_sph_coeff(this, n, m) result(val)
        class(sph_harm_coeff), intent(inout) :: this
        integer, intent(in) :: n, m
        integer :: start, shift
        start = this%get_m_ptr(n)
        !order -m, -(m-1), ... m-1, m
        shift = n + m !!= 0 if n = -m, 
        val = this%data(shift+start)
    end function
    subroutine set_sph_coeff(this, n, m, val)
        class(sph_harm_coeff), intent(inout) :: this
        integer, intent(in) :: n, m
        real(kind), intent(in) ::val
        integer :: start, shift
        start = this%get_m_ptr(n)
        !order -m, -(m-1), ... m-1, m
        shift = n + m !!= 0 if n = -m, 
        this%data(shift+start) = val
    end subroutine
    pure integer elemental function get_m_ptr(this, n)
        class(sph_harm_coeff), intent(in) :: this
        integer, intent(in) :: n
        get_m_ptr = n*n + 1
        ! n = 0 => 0. n = 1 => 1.    n = 2 => 4 (+ 5 = 9) (+7 = 16) (+9 = 25)
    end function
    subroutine alloc_sph_coeff(this,n)
        class(sph_harm_coeff) :: this
        integer :: n, numdata
        numdata = (n+1)**2
        allocate(this%data(numdata))
        allocate(this%m_ptr(n+2))
        ! Each level (2n+1). Sum k=0 N (2k+1)
    end subroutine




    function compute_clustg(clust) result(g_lr)
        type(cluster), intent(inout) :: clust
        integer :: xx,yy,zz, aa
        real(kind) :: x, y,z, xn,yn,zn
        real(kind) :: g_lr(0:p-1, 0:p-1, 0:p-1)
        g_lr = 0

        do zz = 0, p-1
            z = nodes(zz)
            do yy = 0, p-1
                y = nodes(yy)
                do xx = 0, p-1
                    x = nodes(xx)
                    do aa = 1, clust%N
                        xn = clust%poles(1,aa)
                        yn = clust%poles(2,aa)
                        zn = clust%poles(3,aa)
                        g_lr(xx,yy,zz) = g_lr(xx,yy,zz) + clust%weights(aa)/sqrt((x-xn)**2 + (y-yn)**2 + (z-zn)**2)
                    end do
                end do
            end do
        end do
    end function

    real(kind) function eval_lr(r, g_lr) result(V)
        real(kind), intent(in) :: g_lr(0:p-1, 0:p-1, 0:p-1)
        real(kind), intent(in) :: r(3)

    end function

    subroutine setup_ch()
        call compute_chebyshev_nodes()
        call setup_ch_pols()
    end subroutine

end module