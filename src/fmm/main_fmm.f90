program main
    use constants
    use fmm
    use tree_module

    implicit none
    integer, parameter ::numparticles =50
    real(kind) :: cords(3,numparticles), weights(numparticles)
    type(tree), pointer :: root
    call rand(cords)
    call root_constructor(root, cords, weights)
    call root%split()
    write(*,*)"number of children:", root%numChild
    call print(root, 0) 


    contains

    subroutine rand(arr)
        real(kind), intent(inout) :: arr(:,:)


        ! Initialize random seed
        call random_seed()
        ! Fill array with random numbers in [0,1)
        call random_number(arr)
    end subroutine
end program