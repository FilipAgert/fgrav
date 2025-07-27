submodule(fmm) math
    contains

    real(kind) module function  get_sph_coeff(this, n, m) result(val)
        class(sph_harm_coeff), intent(in) :: this
        integer, intent(in) :: n, m
        integer :: start, shift
        start = get_m_ptr(n)
        !order -m, -(m-1), ... m-1, ms
        shift = n + m !!= 0 if n = -m, 
        val = this%data(shift+start)
    end function
    module subroutine set_sph_coeff(this, n, m, val)
        class(sph_harm_coeff), intent(inout) :: this
        integer, intent(in) :: n, m
        real(kind), intent(in) ::val
        integer :: start, shift
        start = get_m_ptr(n)
        !order -m, -(m-1), ... m-1, m
        shift = n + m !!= 0 if n = -m, 
        this%data(shift+start) = val
    end subroutine

    module subroutine add_sph_coeff(this, n, m, val)
        class(sph_harm_coeff), intent(inout) :: this
        integer, intent(in) :: n, m
        real(kind), intent(in) ::val
        integer :: start, shift
        start = get_m_ptr(n)
        !order -m, -(m-1), ... m-1, m
        shift = n + m !!= 0 if n = -m, 
        this%data(shift+start) = this%data(shift+start) + val
    end subroutine

    module subroutine mul_sph_coeff(this, n, m, val)
        class(sph_harm_coeff), intent(inout) :: this
        integer, intent(in) :: n, m
        real(kind), intent(in) ::val
        integer :: start, shift
        start = get_m_ptr(n)
        !order -m, -(m-1), ... m-1, m
        shift = n + m !!= 0 if n = -m, 
        this%data(shift+start) = this%data(shift+start) * val
    end subroutine

    integer function get_m_ptr(n)
        integer, intent(in) :: n
        integer, dimension(0:p), save :: ptr
        integer :: i
        logical, save :: first_time = .true.
        if(first_time) then
            do i = 0,p
                ptr(i) = i*i+1
            end do
        endif
        get_m_ptr = ptr(n)
        ! n = 0 => 0. n = 1 => 1.    n = 2 => 4 (+ 5 = 9) (+7 = 16) (+9 = 25)
    end function

    module subroutine Ynm(Y,theta, phi)  !!Gets spherical harmonics from n=0..p evaluated at theta,phi
        real(kind), intent(in) :: theta,phi
        type(sph_harm_coeff), save :: Nnm !!normalization constant
        type(sph_harm_coeff), intent(out) :: Y
        type(sph_harm_coeff) :: Pnm
        logical, save :: first_time = .true.
        integer :: nn, mm
        real(kind) :: norm, val
        if(first_time) then
            do nn = 0,p
                do mm = -nn, nn
                    norm = sqrt((2.0_kind*nn+1.0_kind) * fac(nn-mm)) / sqrt(4.0_kind*pi * fac(nn+mm))
                    call Nnm%set(nn,mm,norm)
                end do
            end do
            first_time = .false.
        endif
        call compute_legendre_cos_gamma(theta, Pnm)
        do nn = 0,p
            do mm = -nn,nn
                val = 1
                if(mm < 0) then
                    val =sin(-mm*phi)
                else if (mm > 0) then
                    val = cos(mm*phi)
                endif
                call Y%set(nn,mm, Nnm%get(nn,mm) * Pnm%get(nn,mm) * val)
            end do
        end do
    end subroutine

    real(kind) function fac(n)
        integer, intent(in) :: n
        logical, save :: first_time = .true.
        integer, parameter :: maxf = 170
        real(kind), save :: f(0:maxf)
        integer ::i
        if(first_time) then
            f(0)= 1
            do i =0,maxf-1
                f(i+1) = f(i) * (i+1)
            end do
            first_time = .false.
        endif
        fac = f(n)
    end function

    real(kind) function ffac(n)!!double factorial
        integer, intent(in) :: n
        logical, save :: first_time = .true.
        integer, parameter :: fmax = 301
        real(kind), save :: f(0:fmax)
        integer ::i
        if(first_time) then
            f(0)= 1
            f(1)=1
            do i =1,fmax-1
                f(i+1) = f(i-1) * (i+1)
            end do
            first_time = .false.
        endif
        ffac = f(n)
    end function

    subroutine compute_legendre_cos_gamma(gamma, Pnm)
        real(kind), intent(in) :: gamma
        real(kind) :: x, val, y
        type(sph_harm_coeff), intent(inout) :: Pnm
        integer :: n,m
        call Pnm%set(0,0,1.0_kind)

        x = cos(gamma)
        y = sin(gamma) !! == sqrt(1-x^2)
        !https://en.wikipedia.org/wiki/Associated_Legendre_polynomials#Recurrence_formula 
        do n = 1, p
            val = -(2*n-1) * y * Pnm%get(n-1,n-1) !recurrence to get top of chain.
            call Pnm%set(n,n, val)
            val = x*(2*n-1)*Pnm%get(n-1,n-1) !recurrence to get one lower.
            call Pnm%set(n-1,n, val)
            !top two m values are set. use recurrence to find the others
            do m = n-1, -(n-1), -1
                val = (2.0_kind*m*x*Pnm%get(n,m)/y - Pnm%get(n,m+1))/((n+m)*(n-m+1.0_kind))
                call Pnm%set(n,m-1,val)
            end do
        end do
    end subroutine

end submodule 