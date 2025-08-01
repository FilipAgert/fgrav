module mulpol
    use constants, only: kind, ckind
    implicit none
    private
    public :: cluster, calc_mulp_exp, eval_mulp_exp_c, Ynm, sph_harm_coeff_d, sph_harm_coeff_c, eval_grad_mulpexp, calc_mulp_gl_exp, Gmulpol
    integer, parameter :: p = 8 !!order of interpolation
    real(kind) :: Gmulpol
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
        module pure elemental real(kind) function get_sph_coeff(this, n, m)
            class(sph_harm_coeff_d), intent(in) :: this
            integer, intent(in) :: n, m
        end function
        module pure elemental complex(ckind) function get_sph_coeff_c(this, n, m)
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

    pure function cart2sph(cart) result(sph)
        !!Converts cartesian to spherical coordinates
        real(kind), intent(in) :: cart(3)!!Cartesian coordinates
        real(kind) :: sph(3) !!spherical coordinates
        sph(1) = sqrt(sum(cart*cart,1))
        sph(2) = acos(cart(3)/sph(1))
        sph(3) = atan2(cart(2),cart(1))
    end function


    pure function vec_sph2cart(vecSph,sphpos) result(vecCart)
        !!Converts a vector from spherical to cartesian. E.g. acceleration vector
        !!https://en.wikipedia.org/wiki/Del_in_cylindrical_and_spherical_coordinates
        real(kind), intent(in) :: vecSph(3)!! vector in spherical coordinates
        real(kind), intent(in) :: sphpos(3)!! position in spherical coordinates
        real(kind) :: vecCart(3) !!Vector in cartesian coordinates
        real(kind):: st, ct, sp, cp
        st = sin(sphpos(2))
        ct = cos(sphpos(2))
        sp = sin(sphpos(3))
        cp = cos(sphpos(3))
        vecCart(1) = st*cp*vecSph(1) + ct*cp*vecSph(2) - sp*vecSph(3)
        vecCart(2) = st*sp*vecSph(1) + ct*sp*vecSph(2) + cp*vecSph(3)
        vecCart(3) = ct*vecSph(1) - st*vecSph(2)
    end function

    
    subroutine calc_mulp_exp(clust)!!Calculates multipole expansion coefficients from a cluster.
        type(cluster), intent(inout) :: clust
        integer :: n, m, i
        real(kind) :: rhoraised, w, sph(3)
        complex(ckind) :: z
        type(sph_harm_coeff_c) :: Y
        clust%mp_exp%data=0
        do i = clust%startidx, clust%stopidx - 1
            sph = cart2sph(clust%pos(:,i)-clust%cluster_pos)!psherical coordinates wrp center of cluster
            w = clust%weights(i)*Gmulpol
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
        integer :: m_start,m_end, jnm
        clust%mp_gl_exp%data=0
        sph = cart2sph(clust%cluster_pos)
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
                    m_start = max(-n, k - (j - n))
                    m_end   = min( n, k + (j - n))
                    do m = m_start, m_end
                        
                        zm = zm +clust%mp_exp%get(j-n,k-m) * complex(0,1.0_kind)**(abs(k)-abs(m)-abs(k-m))*Anm(n,m)*Anm(j-n,k-m)*Y%get(n,-m)

                        ! if(m < m_start) then
                        !     write(*,'(a,i3,1x,i3)') "m<m_start: j-n, k-m", j-n, k-m
                        !     write(*,'(a,i3,1x,i3)') "m_start,m_stop: j-n, k-m", m_start, m_end
                        ! else if(m > m_end) then 
                        !     write(*,'(a,i3,1x,i3)') "m>m_end: j-n, k-m", j-n, k-m
                        ! endif
                    end do
                    z = z + zm*rhopow

                end do
                z = z/Anm(j,k)
                !write(*,*) j,k,z
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
        real(kind) :: sph(3) 
        real(kind) :: g(3) !!gradient in cartesian coordinates
        real(kind) :: gsph(3) !!gradient in spherical coordinates
        real(kind), intent(in) :: r(3)!!point in cartesian coordinates
        real(kind) :: cott, isint, rpow
        complex(ckind) :: zr, zt, zp, phase, Ynm, cfnm
        complex(ckind), parameter :: imag = complex(0.0_kind, 1.0_kind)
        integer :: m,n
        sph = cart2sph(r)
        call Ynm(Y, sph(2), sph(3))
        cott = cotan(sph(2))
        isint = 1.0_kind/sin(sph(2))
        rpow = sph(1)
        phase = exp(complex(0.0_kind,-1.0_kind) * sph(3))
        !write(*,*) cott, isint
        gsph=0
        do n = 0,p
            zr = 0
            zt = 0
            zp = 0

            do m = -n, n-1
                Ynm = Y%get(n,m)
                cfnm = cf%get(n,m)
                zr = zr + cfnm* Ynm
                zp = zp + cfnm*imag * m * isint *  Ynm
                !zt = zt + cf%get(n,m)*dYdt%get(n,m)
                ! dYnm/dtheta = n cotan(theta) Ynm - 1/(sin theta) * Y_(n-1)m * sqrt(n+|m|)*sqrt(n-|m|)

                zt = zt +cfnm * (n*cott * Ynm - isint*Y%get(n-1,m)*sqrt(1.0_kind*(n*n-m*m)))
            end do
            m = n
            Ynm = Y%get(n,m)
            cfnm = cf%get(n,m)
            zr = zr + cfnm* Ynm
            zp = zp + cfnm*imag * m * isint * Ynm
            zt = zt + cfnm* (n*cott * Ynm)

            rpow = rpow*sph(1)!r^n+2
            gsph(1) = gsph(1)+ (-n-1)*real(zr,kind)/rpow
            gsph(2) = gsph(2)+ real(zt,kind)/rpow
            !write(*,'(6e15.3)') zr,zt,zp
            gsph(3) = gsph(3)+ real(zp,kind)/rpow
        end do
        g = vec_sph2cart(gsph, sph)
    end function

    real(kind) function Anm(n,m)
        logical, save :: first_time = .true.
        type(sph_harm_coeff_d), save :: A
        integer, intent(in) :: n,m
        integer :: nn,mm

        if(first_time) then
            do nn = 0,p
                do mm = -nn,nn
                    call A%set(nn,mm, ((-1)**nn) *1.0_kind / sqrt(fac(nn-mm)*fac(nn+mm)))
                end do
            end do
            first_time = .false.
        endif
        !error stop
        Anm = A%get(n,m)
    end function


end module