module sim
    use constants
    use body_module
    implicit none


    contains

    subroutine adaptive_step(bodies, dt, tol)
        type(body), dimension(:), intent(inout) :: bodies
        real(kind), intent(inout) :: dt
        real(kind), intent(in) :: tol
        type(body), dimension(size(bodies)) :: bodies_half, bodies_full
        real(kind) :: err
        integer :: i
        logical :: accept
        real(kind) :: ep(3), p0(3,size(bodies)), dp(3)

        accept = .false.
        do i = 1, size(bodies)
            p0(:,i) = bodies(i)%position
        end do
        do while (.not. accept)
            bodies_half = bodies
            bodies_full = bodies

            ! Full step
            call step(bodies_full, dt)

            ! Two half steps
            call step(bodies_half, dt/2.0)
            call step(bodies_half, dt/2.0)

            err = 0.0
            do i = 1, size(bodies)
                ep = bodies_full(i)%position - bodies_half(i)%position
                dp = bodies_half(i)%position

                err = max(err, maxval(abs(ep/dp)))
            end do

            if (err < tol) then
                accept = .true.
                bodies = bodies_half  ! use better estimate
                if (err < tol / 4.0) dt = dt * 2.0  ! increase step size
            else
                dt = dt / 2.0  ! reduce step size
            end if
        end do
    end subroutine


    subroutine step(bodies, dt)
        type(body), dimension(:), intent(inout) :: bodies
        real(kind), intent(in) :: dt
        real(kind), dimension(3,size(bodies)*2)::y, dy  
        real(kind), dimension(size(bodies)) :: M
        integer ::bi
        real(kind), dimension(3,2*size(bodies)) ::k1,k2,k3,k4
        integer :: N
        logical, dimension(size(bodies)) :: negligble
        N = size(bodies)
        do bi = 1, N
            y(:,bi) = bodies(bi)%position
            y(:,bi+N) = bodies(bi)%velocity
            M(bi) = bodies(bi)%mass
            negligble(bi) = bodies(bi)%negligble
        end do
        k1 = f(y,M,N, negligble)
        k2 = f(y+dt*k1/2.0,M,N, negligble)
        k3 = f(y+dt*k2/2.0, M,N, negligble)
        k4 = f(y+dt*k3,M,N, negligble)
        dy = dt * (k1+2.0*k2+2.0*k3+k4) / 6.0_kind
        y = y + dy

        do bi = 1, N
            bodies(bi)%position = y(:,bi)
            bodies(bi)%velocity = y(:,bi+N)
            !write(*,'(A,A,2x,6E10.1)') "position change for ",  bodies(bi)%name, dy(1:2,bi)
        end do

    end subroutine

    function f(y,M,N, negligble)
        real(kind), dimension(3,2*N) ,intent(in) :: y
        logical, dimension(N), intent(in) :: negligble
        integer, intent(in) ::N
        real(kind), intent(in) :: M(N)
        real(kind) :: f(3,N*2)
        integer :: i, j

        do i = 1, N
            f(:,i) = y(:,i + N)
        end do

        do i = 1, N
            f(:,i+N) = 0
            do j = 1,N
                if(i /= j .and. .not. negligble(j)) then
                    f(:,i+N) = f(:,i+N) + a_grav(y(:,i), y(:,j), M(j))
                end if
            end do
        end do

    end function


    

    !!Gravitational acceleration between two bodies
    pure function a_grav(x1,x2,m2)
        real(kind), dimension(3), intent(in) :: x1,x2
        real(kind), intent(in) :: m2
        real(kind) :: r
        real(kind) :: a_grav(3)
        r = sqrt(sum((x1-x2)**2,1))
        a_grav = (x2-x1)*G * m2 / r**3
    end function a_grav


end module