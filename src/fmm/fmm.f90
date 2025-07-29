module fmm
    use constants, only: kind, ckind
    implicit none
    private
    public :: cluster, calc_mulp_exp, eval_mulp_exp_c, Ynm, sph_harm_coeff_d, sph_harm_coeff_c, eval_grad_mulpexp, calc_mulp_gl_exp
    integer, parameter :: p = 6 !!order of interpolation
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
        complex(ckind), private :: data((p+1)**2)
        contains
            procedure :: get => get_sph_coeff_c
            procedure :: set => set_sph_coeff_c
            procedure :: add => add_sph_coeff_c
           ! procedure :: mulval => mul_sph_coeff_c
    end type

    type cluster
        real(kind), pointer :: pos(:,:) ! x, y, z global coordinates.
        real(kind) :: cluster_pos(3) ! cluster center
        real(kind), pointer :: weights(:)
        integer :: startidx = 0, stopidx=0 !!start(inclusive) and stop(exclusive) index of cluster's own values in the pos and weights arrays
        type(sph_harm_coeff_c) :: mp_exp!!multipole expansion coefficients
        type(sph_harm_coeff_c) :: mp_gl_exp !!Global coordinate system multipole expansion coefficient
    end type

     
    interface
        module real(kind) function get_sph_coeff(this, n, m)
            class(sph_harm_coeff_d), intent(in) :: this
            integer, intent(in) :: n, m
        end function
        module complex(ckind) function get_sph_coeff_c(this, n, m)
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
            complex(ckind), intent(in) ::val
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
            complex(ckind), intent(in) ::val
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

         real(kind) module function fac(n)
            integer, intent(in) :: n
        end function
    end interface


    contains
    pure function toSpherical(cart) result(sph)
        real(kind), intent(in) :: cart(3)
        real(kind) :: sph(3)
        sph(1) = sqrt(sum(cart*cart,1))
        sph(2) = acos(cart(3)/sph(1))
        sph(3) = atan2(cart(2),cart(1))
    end function
    subroutine calc_mulp_exp(clust)!!Calculates multipole expansion coefficients from a cluster.
        type(cluster), intent(inout) :: clust
        integer :: n, m, i
        real(kind) :: rhoraised, w, sph(3)
        complex(ckind) :: z
        type(sph_harm_coeff_c) :: Y
        clust%mp_exp%data=0
        write(*,*) clust%startidx, clust%stopidx
        do i = clust%startidx, clust%stopidx - 1
            sph = toSpherical(clust%pos(:,i)-clust%cluster_pos)!psherical coordinates wrp center of cluster
            w = clust%weights(i)
            call Ynm(Y, sph(2), sph(3))
            rhoraised = 1.0_kind/sph(1)
            do n = 0,p
                rhoraised = rhoraised*sph(1)!r^n
                do m = -n,n
                    z = Y%get(n,-m)
                    call clust%mp_exp%add(n, m, z*rhoraised*w)
                end do
            end do
        end do
    end subroutine

    subroutine calc_mulp_gl_exp(clust)!!Calculates global multipole expansion coefficients from a cluster given the local one is already calculated.
        type(cluster), intent(inout) :: clust
        integer :: n, m, k, j
        real(kind) :: rho, alpha, beta, rhopow, sph(3)
        complex(ckind) :: z, zm
        type(sph_harm_coeff_c) :: Y
        clust%mp_gl_exp%data=0
        sph = toSpherical(clust%cluster_pos)
        rho = sph(1)
        alpha= sph(2)
        beta = sph(3)
        call Ynm(Y, alpha, beta)
        do j = 0,p
            do k = -j, j
                z = 0
                do n = 0,j
                    if(n ==0) then
                        rhopow = 1
                    else 
                        rhopow = rhopow * rho
                    endif
                    zm = 0
                    do m = -n, n
                        zm = zm + clust%mp_exp%get(j-n,k-m) * complex(0,1.0_kind)**(abs(k)-abs(m)-abs(k-m))*Anm(n,m)*Anm(j-n,k-m)*Y%get(n,-m)
                    end do
                    z = z + zm*rhopow
                end do
                z = z/Anm(k,j)
                call clust%mp_gl_exp%set(j,k,z)
            end do
        end do
    end subroutine

    function eval_mulp_exp_c(clust, r) result(Psi)
        type(cluster), intent(in) :: clust
        real(kind), intent(in) :: r(3) !in spherical coordinates
        type(sph_harm_coeff_c) :: Y
        real(kind) ::Psi, rpow
        complex(ckind) :: z
        integer :: n,m
        call Ynm(Y, r(2), r(3))
        Psi = 0
        rpow = 1
        do n = 0,p
            z = 0

            do m = -n, n
                z = z + Y%get(n,m)*clust%mp_exp%get(n,m)
            end do
            rpow = rpow*r(1)!r^n+1
            Psi = Psi + REAL(z,kind)/rpow
        end do
    end function

    !Evaluates the gradient of the potential given by coefficients g at point r.
   function eval_grad_mulpexp(cf, r) result(g)
        type(sph_harm_coeff_c),intent(in) :: cf
        type(sph_harm_coeff_c) :: Y
        real(kind), intent(in) :: r(3) !!point in spherical coordinates
        real(kind) :: g(3) !!gradient in spherical coordinates
        real(kind) :: cott, isint, rpow
        complex(ckind) :: zr, zt, zp, phase
        integer :: m,n
        call Ynm(Y, r(2), r(3))
        cott = cotan(r(2))
        isint = 1.0_kind/sin(r(2))
        rpow = r(1)
        phase = exp(complex(0.0_kind,-1.0_kind) * r(3))
        write(*,*) cott, isint
        g=0
        do n = 0,p
            zr = 0
            zt = 0
            zp = 0

            do m = -n, n
                zr = zr + cf%get(n,m)* Y%get(n,m)
                zp = zp + cf%get(n,m)*complex(0.0_kind, 1.0_kind) * m * isint *  Y%get(n,m)
                !zt = zt + cf%get(n,m)*dYdt%get(n,m)
                ! dYnm/dtheta = n cotan(theta) Ynm - 1/(sin theta) * Y_(n-1)m * sqrt(n+|m|)*sqrt(n-|m|)
                zt = zt + cf%get(n,m) * (n*cott * Y%get(n,m) - isint*Y%get(n-1,m)*sqrt(1.0_kind*(n*n-m*m)))
            end do
            rpow = rpow*r(1)!r^n+2
            g(1) = g(1)+ (-n-1)*real(zr,kind)/rpow
            g(2) = g(2)+ real(zt,kind)/rpow
            !write(*,'(6e15.3)') zr,zt,zp
            g(3) = g(3)+ real(zp,kind)/rpow
        end do

    end function

    real(kind) function Anm(n,m)
        logical, save :: first_time = .true.
        type(sph_harm_coeff_d), save :: A
        integer, intent(in) :: n,m
        integer :: nn,mm

        if(first_time) then
            do nn = 0,p
                do mm = -nn,nn
                    call A%set(nn,mm, (-1)**n *1.0_kind / sqrt(fac(n-m)*fac(n+m)))
                end do
            end do
            first_time = .false.
        endif
        Anm = A%get(n,m)
    end function


end module