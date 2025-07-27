program main
    use constants
    use body_module
    use sim
    use out
    implicit none

    real(kind) :: dt = 1.0 ! days
    integer :: s, numsteps, numdays
    real(kind) :: dayctr
    real(kind), dimension(3) :: pos, vel
    integer :: i
    real :: r, angle, speed
    real, parameter :: max_r = 4.0, max_v = 0.02
 
    type(body), dimension(4) :: bodies
    real(kind) :: COM(3)
    real(kind) :: E0, T, V,alpha, scale, positions(size(bodies),3), lasttime, tp
    integer :: ii

    
    call bodies(2)%initialize("Moon", 0.012*3e-6_kind, 0.0_kind)
    call bodies(1)%initialize("Earth", 3.003e-6_kind, 0.0_kind)
    call bodies(3)%initialize("Asteroid", 0.0_kind, 0.0_kind, .true.)
    call bodies(4)%initialize("x", 3e-7_kind, 0.0_kind)

    bodies(2)%position = [0.002578,0.0,0.0]
    bodies(2)%velocity = [0.0, 5.88e-4, 0.0] ! circular orbit around Earth
    bodies(3)%position = [0.0,0.004,0.0]; bodies(3)%velocity = [-2.5e-4, 0.0, 0.0] 
    bodies(4)%position = [0.0,0.012,0.0]; bodies(4)%velocity = [5e-3, -5e-4, 0.0]
    
    ! Simulation setup
    numdays = 365 * 15 ! simulate 15 years
    numsteps = int(numdays / dt)
    dayctr = 0.0
    E0 = Vpot(bodies)+Tkin(bodies)
    call write_solar_system(bodies, dayctr)
    s = 0
    do while(dayctr < numdays)
        s = s +1


        if (mod(s, 1) == 0) then
            call write_solar_system(bodies, dayctr)
            COM = 0
            do ii = 1, size(bodies)
                COM=COM+bodies(ii)%position*bodies(ii)%mass
            end do
            COM = com/size(bodies)
            do ii = 1,size(bodies)
                positions(ii,:) = bodies(ii)%position-COM
            end do
            
            scale = 0.005!int(maxval(sqrt(sum(positions**2,1)))/2.0 + 1)

            call print_solar_system(bodies, scale, 'abs')

            write(*,'(A,f8.1,A,A,f8.2,a,I8)') "T = ", dayctr, " d", " dt = ", dt, " d, steps: ", s
            T = Tkin(bodies)
            V = Vpot(bodies)
            alpha = T/abs(V)
            write(*,'(A,f10.3,a,f10.3)') "Tkin/V", alpha, " E/E0:", (T+V)/E0
            tp=(dayctr-lasttime)/3
            call pause(tp)
            lasttime = dayctr
        end if

        call adaptive_step(bodies, dt,1e-12_kind)
        dayctr = dayctr + dt
    end do

contains

    subroutine pause(t)
        use iso_c_binding
        real(kind) :: t
        interface
            subroutine usleep(microseconds) bind(C, name="usleep")
                import :: c_int
                integer(c_int), value :: microseconds
            end subroutine usleep
        end interface
        call usleep(int(t*1e6))  ! microseconds
    end subroutine pause

    real(kind) function Vpot(bodies)
        type(body), intent(in) :: bodies(:)
        integer :: ii,jj
        real(kind) :: d(3)
        real(kind), parameter ::AU = 1.496e11
        Vpot = 0
        do ii = 1,size(bodies)
            do jj = ii+1,size(bodies)
                d = bodies(ii)%position-bodies(jj)%position
                Vpot = Vpot - 6.67430e-11*bodies(ii)%mass*bodies(jj)%mass*M_sun*M_sun/(AU*sqrt(sum(d**2)))
            end do
        end do


    end function
    real(kind) function Tkin(bodies)
        type(body), intent(in) :: bodies(:)
        integer :: ii
        Tkin = 0
        do ii = 1,size(bodies)
            Tkin = Tkin + bodies(ii)%mass * M_sun * sum((bodies(ii)%velocity*1.731e6)**2)/2.0
        end do

    end function

end program main
