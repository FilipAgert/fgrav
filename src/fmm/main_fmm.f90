program main
    use constants
    use fmm


    implicit none
    type(cluster) :: cl
    real(kind) :: pos(3), weight
    pos = [1.0_kind, pi/2,0.0_kind]
    weight = 1.0_kind
    allocate(cl%pos(3,1), cl%weights(1))
    cl%pos(:,1) = pos
    cl%weights = weight
    call calc_mulp_exp_c(cl)
    
    write(*,*) "hello world!"

end program