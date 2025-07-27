module out
    use constants
    use body_module
    implicit none
    contains

    subroutine write_solar_system(bodies, dayctr)
        implicit none
        type(body), dimension(:), intent(in) :: bodies
        real(kind), intent(in) :: dayctr
        integer ::ii, unit
        logical :: first_time = .true.
        character(len=200) :: file = "out/run.dat"
        unit = 3
        if(first_time) then
            open(unit = unit, file=file, status='replace', action='write')
            write(unit,'(A)')"BODY      MASS (M_SUN)    RAD(km)       x(AU)      y(AU)       z(AU)      vx (AU/d)      vy (AU/d)     vz(AU/d)"
            first_time = .false.
        else
            open(unit = unit, file=file, status='old', action='readwrite', position='append')
        endif
        write(unit,'(a,f10.1,a)')'T= ',dayctr," d"
        do ii = 1, size(bodies)
            write(unit, '(A20, 1x,2e10.3,1x, 6e20.8)') bodies(ii)%name, bodies(ii)%mass, bodies(ii)%radius, bodies(ii)%position(:), bodies(ii)%velocity
        end do
        close(unit)
    end subroutine

subroutine print_solar_system(bodies, scale, frame)
    implicit none
    type(body), dimension(:), intent(in) :: bodies
    character(len=3) :: frame
    real(kind), intent(in) :: scale
    integer, parameter :: width = 71, height = 35
    character(len=width) :: line
    character(len=1), dimension(width, height) :: grid
    real(kind) :: x, y
    integer :: i, gx, gy
    integer :: cx, cy, row
    real(kind) :: total_mass
    real(kind), dimension(3) :: com  ! center of mass

    ! Compute center of mass
    total_mass = 0.0
    com = 0.0
    do i = 1, size(bodies)
        total_mass = total_mass + bodies(i)%mass
        com = com + bodies(i)%mass * bodies(i)%position
    end do
    com = com / total_mass
    if(frame/='com') com =0

    ! Clear screen and move cursor home
    write(*,'(a)',advance='no') char(27)//'[2J'//char(27)//'[H'

    cx = width / 2
    cy = height / 2
    grid = ' '

    ! Draw axes
    grid(:, cy) = '-'
    do row = 1, height
        grid(cx, row) = '|'
    end do
    grid(cx, cy) = '+'

    ! Place bodies relative to COM
    do i = 1, size(bodies)
        x = bodies(i)%position(1) - com(1)
        y = bodies(i)%position(2) - com(2)
        gx = nint((x / scale) * (width / 2)) + cx
        gy = nint((y / scale) * (height / 2)) + cy

        if (gx >= 1 .and. gx <= width .and. gy >= 1 .and. gy <= height) then
            grid(gx, gy) = bodies(i)%name(1:1)
        end if
    end do

    ! Print all rows
    do row = height, 1, -1
        line = ''
        do i = 1, width
            line(i:i) = grid(i, row)
        end do
        write(*,'(a)') line
    end do
    write(*,'(a,f6.2,a)') "Scale: ±", scale, " AU (centered on COM)"
end subroutine



end module