module fmm
    use constants
    implicit none
    private
    integer, parameter :: p = 5 !!order of interpolation
    type sph_harm_coeff !!Data structure for storing spherical harmonic coefficients c_n^m
    !since each n stores (2n+1) coefficients, a 2d array is inefficient. 1d array is better.
        real(kind), private :: data((p+1)**2)
        contains
            procedure :: get => get_sph_coeff
            procedure :: set => set_sph_coeff
            procedure :: add => add_sph_coeff
            procedure :: mulval => mul_sph_coeff
    end type

    type cluster
        real(kind), allocatable :: pos(:,:) ! r, theta, phi
        real(kind), allocatable :: weights(:)
        type(sph_harm_coeff) :: mp_exp!!multipole expansion coefficients
    end type

     
    interface
        module real(kind) function get_sph_coeff(this, n, m)
            class(sph_harm_coeff), intent(in) :: this
            integer, intent(in) :: n, m
        end function

        module subroutine set_sph_coeff(this, n,m,val)
            class(sph_harm_coeff), intent(inout) :: this
            integer, intent(in) :: n, m
            real(kind), intent(in) ::val
        end subroutine

        ! M(n,m) += val
        module subroutine add_sph_coeff(this, n,m,val)
            class(sph_harm_coeff), intent(inout) :: this
            integer, intent(in) :: n, m
            real(kind), intent(in) ::val
        end subroutine

        ! M(n,m) *= val
        module subroutine mul_sph_coeff(this, n,m,val)
            class(sph_harm_coeff), intent(inout) :: this
            integer, intent(in) :: n, m
            real(kind), intent(in) ::val
        end subroutine

        module subroutine alloc_sph_coeff(this, n)
            class(sph_harm_coeff), intent(inout) :: this
            integer, intent(in) :: n
        end subroutine

        module subroutine Ynm(Y,theta, phi)
            type(sph_harm_coeff), intent(out) ::Y !!Contains factor part of function fnm = Nnm * Pnm(cos theta)
            real(kind), intent(in) :: theta,phi
        end subroutine

    end interface


    contains

        subroutine calc_mulp_exp_c(clust)!!Calculates multipole expansion coefficients from a cluster.
            type(cluster), intent(inout) :: clust
            integer :: n, m, i
            real(kind) :: s, rho, alpha, beta, rhoraised, w
            type(sph_harm_coeff) :: Y
            clust%mp_exp%data=0

            do i = 1, size(clust%weights,1)
                rho = clust%pos(1,i)
                alpha = clust%pos(2,i)
                beta = clust%pos(3,i)
                w = clust%weights(i)
                call Ynm(Y, alpha, beta)
                do n = 0,p
                    rhoraised = rho**n
                    do m = -p,p
                        s = Y%get(n,m)
                        call clust%mp_exp%add(n, m, s*rhoraised*w)
                    end do
                end do
            end do
        end subroutine


end module