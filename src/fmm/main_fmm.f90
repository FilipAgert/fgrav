program main
    use constants
    use mulpol
    use tree_module
    use fmm

    implicit none
    integer ::numparticles =1000
    real(kind),allocatable :: cords(:,:), weights(:), acc(:,:)
    integer :: ii, status
    character(len=100) :: str
    call get_command_argument(1, str,status=status)
    if (status .eq. 0) read(str,'(i100)') numparticles
    allocate(cords(3,numparticles), weights(numparticles), acc(3,numparticles))
    call rand(cords)
    weights = 1
    acc = eval_acc(cords, weights, G_Au)
    ! do ii = 1, numparticles
    !     write(*,'(a,3(f6.3,1x),a,3(f6.3,1x),a)') 'x: [', cords(:,ii),'], a=[',acc(:,ii),']'
    ! end do

    contains

    subroutine rand(arr)
        real(kind), intent(inout) :: arr(:,:)


        ! Initialize random seed
        call random_seed()
        ! Fill array with random numbers in [0,1)
        call random_number(arr)
    end subroutine
end program