program main
    use constants
    use fmm
    use tree_module

    implicit none
    integer, parameter ::numparticles =4000
    real(kind) :: cords(3,numparticles), weights(numparticles)
    type(tree), pointer :: root
    allocate(root)
    call rand(cords)
    root%bounds(:,1) = 0
    root%bounds(:,2) = 1
    root%width = 1
    allocate(root%clust)
    allocate(root%clust%pos(3,numparticles))
    allocate(root%clust%weights(numparticles))
    root%clust%pos = 0.5
    root%clust%pos = cords
    root%clust%weights = 1
    call root%split()
    write(*,*)"number of children:", root%numChild
    !call print(root, 0) 


    contains

    subroutine rand(arr)
        real(kind), intent(inout) :: arr(:,:)


        ! Initialize random seed
        call random_seed()
        ! Fill array with random numbers in [0,1)
        call random_number(arr)
    end subroutine
end program