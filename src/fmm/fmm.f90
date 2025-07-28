module fmm
    use constants
    implicit none
    private
    public :: cluster, calc_mulp_exp_c, eval_mulp_exp_c, Ynm, sph_harm_coeff_d, sph_harm_coeff_c
    integer, parameter :: p = 5 !!order of interpolation
    type sph_harm_coeff_d !!Data structure for storing spherical harmonic coefficients c_n^m
    !since each n stores (2n+1) coefficients, a 2d array is inefficient. 1d array is better.
        real(kind), private :: data((p+1)**2)
        contains
            procedure :: get => get_sph_coeff
            procedure :: set => set_sph_coeff
            procedure :: add => add_sph_coeff
            procedure :: mulval => mul_sph_coeff
    end type

    type sph_harm_coeff_c !!Data structure for storing spherical harmonic coefficients c_n^m
    !since each n stores (2n+1) coefficients, a 2d array is inefficient. 1d array is better.
        complex(kind), private :: data((p+1)**2)
        contains
            procedure :: get => get_sph_coeff_c
            procedure :: set => set_sph_coeff_c
            procedure :: add => add_sph_coeff_c
           ! procedure :: mulval => mul_sph_coeff_c
    end type

    type cluster
        real(kind), allocatable :: pos(:,:) ! r, theta, phi
        real(kind), allocatable :: weights(:)
        type(sph_harm_coeff_c) :: mp_exp!!multipole expansion coefficients
    end type

     
    interface
        module real(kind) function get_sph_coeff(this, n, m)
            class(sph_harm_coeff_d), intent(in) :: this
            integer, intent(in) :: n, m
        end function
        module complex(kind) function get_sph_coeff_c(this, n, m)
            class(sph_harm_coeff_c), intent(in) :: this
            integer, intent(in) :: n, m
        end function


        module subroutine set_sph_coeff(this, n,m,val)
            class(sph_harm_coeff_d), intent(inout) :: this
            integer, intent(in) :: n, m
            real(kind), intent(in) ::val
        end subroutine

        module subroutine set_sph_coeff_c(this, n,m,val)
            class(sph_harm_coeff_c), intent(inout) :: this
            integer, intent(in) :: n, m
            complex(kind), intent(in) ::val
        end subroutine

        ! M(n,m) += val
        module subroutine add_sph_coeff(this, n,m,val)
            class(sph_harm_coeff_d), intent(inout) :: this
            integer, intent(in) :: n, m
            real(kind), intent(in) ::val
        end subroutine
        module subroutine add_sph_coeff_c(this, n,m,val)
            class(sph_harm_coeff_c), intent(inout) :: this
            integer, intent(in) :: n, m
            complex(kind), intent(in) ::val
        end subroutine


        ! M(n,m) *= val
        module subroutine mul_sph_coeff(this, n,m,val)
            class(sph_harm_coeff_d), intent(inout) :: this
            integer, intent(in) :: n, m
            real(kind), intent(in) ::val
        end subroutine

        module subroutine Ynm(Y,theta, phi)
            type(sph_harm_coeff_c), intent(out) ::Y !!Contains factor part of function fnm = Nnm * Pnm(cos theta)
            real(kind), intent(in) :: theta,phi
        end subroutine

    end interface


    contains

    subroutine calc_mulp_exp_c(clust)!!Calculates multipole expansion coefficients from a cluster.
        type(cluster), intent(inout) :: clust
        integer :: n, m, i
        real(kind) :: rho, alpha, beta, rhoraised, w
        complex(kind) :: z
        type(sph_harm_coeff_c) :: Y
        clust%mp_exp%data=0

        do i = 1, size(clust%weights,1)
            rho = clust%pos(1,i)
            alpha = clust%pos(2,i)
            beta = clust%pos(3,i)
            w = clust%weights(i)
            call Ynm(Y, alpha, beta)
            do n = 0,p
                rhoraised = rho**n
                do m = -n,n
                    z = Y%get(n,-m)
                    call clust%mp_exp%add(n, m, z*rhoraised*w)
                end do
            end do
        end do
    end subroutine

    function eval_mulp_exp_c(clust, r) result(Psi)
        type(cluster), intent(in) :: clust
        real(kind), intent(in) :: r(3)
        type(sph_harm_coeff_c) :: Y
        real(kind) ::Psi, invrpow
        complex(kind) :: z
        integer :: n,m
        call Ynm(Y, r(2), r(3))
        Psi = 0
        do n = 0,p
            z = 0

            do m = -n, n
                z = z + Y%get(n,m)*clust%mp_exp%get(n,m)
            end do
            invrpow = 1.0_kind/r(1)**(n+1)
            Psi = Psi + REAL(z,kind)*invrpow
        end do
    end function


end module