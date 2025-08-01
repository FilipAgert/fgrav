submodule(mulpol) math
    integer ::i
    integer, parameter :: m_ptr(0:p) = [(i*i+1,i=0,p)]
    contains

    real(kind) pure elemental module function  get_sph_coeff(this, n, m) result(val)
        class(sph_harm_coeff_d), intent(in) :: this
        integer, intent(in) :: n, m
        integer :: start, shift
        if(abs(m) > n) then
            val = 0
            return
        else
            start = m_ptr(n)
            !order -m, -(m-1), ... m-1, ms
            shift = n + m !!= 0 if n = -m, 
            val = this%data(shift+start)
        endif
    end function
    complex(ckind) pure elemental module function  get_sph_coeff_c(this, n, m) result(val)
        class(sph_harm_coeff_c), intent(in) :: this
        integer, intent(in) :: n, m
        integer :: start, shift
        if(abs(m) > n) then
            val = 0
            return
        endif
        start = m_ptr(n)
        !order -m, -(m-1), ... m-1, ms
        shift = n + m !!= 0 if n = -m, 
        val = this%data(shift+start)
    
    end function
    module subroutine set_sph_coeff(this, n, m, val)
        class(sph_harm_coeff_d), intent(inout) :: this
        integer, intent(in) :: n, m
        real(kind), intent(in) ::val
        integer :: start, shift
        start = m_ptr(n)
        !order -m, -(m-1), ... m-1, m
        shift = n + m !!= 0 if n = -m, 
        this%data(shift+start) = val
    end subroutine
    module subroutine set_sph_coeff_c(this, n, m, val)
        class(sph_harm_coeff_c), intent(inout) :: this
        integer, intent(in) :: n, m
        complex(ckind), intent(in) ::val
        integer :: start, shift
        start = m_ptr(n)
        !order -m, -(m-1), ... m-1, m
        shift = n + m !!= 0 if n = -m, 
        this%data(shift+start) = val
    end subroutine

    module subroutine add_sph_coeff(this, n, m, val)
        class(sph_harm_coeff_d), intent(inout) :: this
        integer, intent(in) :: n, m
        real(kind), intent(in) ::val
        integer :: start, shift
        start = m_ptr(n)
        !order -m, -(m-1), ... m-1, m
        shift = n + m !!= 0 if n = -m, 
        this%data(shift+start) = this%data(shift+start) + val
    end subroutine
    module subroutine add_sph_coeff_c(this, n, m, val)
        class(sph_harm_coeff_c), intent(inout) :: this
        integer, intent(in) :: n, m
        complex(ckind), intent(in) ::val
        integer :: start, shift
        start = m_ptr(n)
        !order -m, -(m-1), ... m-1, m
        shift = n + m !!= 0 if n = -m, 
        this%data(shift+start) = this%data(shift+start) + val
    end subroutine

    module subroutine mul_sph_coeff(this, n, m, val)
        class(sph_harm_coeff_d), intent(inout) :: this
        integer, intent(in) :: n, m
        real(kind), intent(in) ::val
        integer :: start, shift
        start = m_ptr(n)
        !order -m, -(m-1), ... m-1, m
        shift = n + m !!= 0 if n = -m, 
        this%data(shift+start) = this%data(shift+start) * val
    end subroutine

    module subroutine Ynm(Y,theta, phi)  !!Gets spherical harmonics from n=0..p evaluated at theta,phi
        real(kind), intent(in) :: theta,phi
        type(sph_harm_coeff_d), save :: Nnm !!normalization constant
        type(sph_harm_coeff_c), intent(out) :: Y
        type(sph_harm_coeff_d) :: Pnm
        logical, save :: first_time = .true.
        integer :: nn, mm
        real(kind) :: norm, cm(-p:p), sm(-p:p), arg
        complex(ckind):: val
        if(first_time) then
            do nn = 0,p
                do mm = -nn, nn
                    !norm = sqrt((2.0_kind*nn+1.0_kind) * fac(nn-abs(mm))) / sqrt(4.0_kind*pi * fac(nn+abs(mm)))
                    norm =  sqrt(fac(nn-abs(mm))) / sqrt(fac(nn+abs(mm)))
                    !write(*,'(a,2i3,a,f10.3)') "n,m:",nn,mm, ", norm:", norm
                    call Nnm%set(nn,mm,norm)
                end do
            end do
            first_time = .false.
        endif
        call compute_legendre_cos_gamma(theta, Pnm)
        cm(0) = 1
        sm(0) = 0
        do mm = 1, p
            arg = mm*phi
            cm(mm) = cos(arg)
            sm(mm) = sin(arg)
        end do
        do mm = -p, -1
            cm(mm) =cm(-mm)
            sm(mm) = -sm(-mm)
        end do
        do nn = 0,p
            do mm = -nn,nn
                val = cmplx(cm(mm),sm(mm), kind=ckind)
                call Y%set(nn,mm, Nnm%get(nn,mm) * Pnm%get(nn,abs(mm)) * val)
            end do
        end do
    end subroutine

    real(kind) module function fac(n)
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
        type(sph_harm_coeff_d), intent(inout) :: Pnm
        integer :: n,m
        call Pnm%set(0,0,1.0_kind)

        x = cos(gamma)
        y = sin(gamma) !! == sqrt(1-x^2)
        !https://en.wikipedia.org/wiki/Associated_Legendre_polynomials#Recurrence_formula 
        if(y==0.0_kind) then!! x is either -1 or +1
            do n = 0,p
                do m = 0, n
                    if(m/=0) then
                        call Pnm%set(n,m,0.0_kind)
                    else
                        if(mod(n,2)==0) then !If even, Pnm is 1 at endpoints for m =0
                            call Pnm%set(n,m,1.0_kind)
                        else
                            call Pnm%set(n,m,x)!If odd, Pnm is -1 at x = -1 and +1 at x = +1
                        endif
                    endif
                end do
            end do
            return
        endif
        do n = 1, p!! y is nonzero.
            val = -(2*n-1) * y * Pnm%get(n-1,n-1) !recurrence to get top of chain.
            call Pnm%set(n,n, val) !(1,1)
            val = x*(2*n-1)*Pnm%get(n-1,n-1) !recurrence to get one lower.
            call Pnm%set(n,n-1, val)!(1,0)
            !top two m values are set. use recurrence to find the others
   
            !write(*,'(2I3,2f10.3)')n,n, Pnm%get(n,n), x
            !write(*,'(2I3,2f10.3)')n,n-1,Pnm%get(n,n-1), x
            do m = n-1, 0, -1
                if(y/=0.0_kind) then
                    val = -(2.0_kind*m*x*Pnm%get(n,m)/y + Pnm%get(n,m+1))/((n+m)*(n-m+1.0_kind))
                else
                    error stop "Should not be here. Legendre error. "
                endif

                call Pnm%set(n,m-1,val)
                !write(*,'(2I3,2f10.3)')n,m-1,Pnm%get(n,m-1), x
            end do
        end do
    end subroutine

end submodule 