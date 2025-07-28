program main
    use constants
    use fmm


    implicit none
    type(cluster) :: cl
    real(kind) :: pos(3), weight, r(3), P1,P2
    complex(kind):: SUM
    type(sph_harm_coeff_c) :: Yt, Ys
    integer :: n,m
    pos = [1.0_kind, 2.0_kind,1.0_kind]
    weight = 1.0_kind
    allocate(cl%pos(3,1), cl%weights(1))
    cl%pos(:,1) = pos
    cl%weights = weight
    call calc_mulp_exp_c(cl)

    r = [50.0_kind, 2.0_kind, 1.0_kind]
    P1 = eval_mulp_exp_c(cl, r)
    P2 = 1.0_kind/abs(r(1)-pos(1))
    call Ynm(Ys, pos(2),pos(3))
    call Ynm(Yt, r(2), r(3))
    do n = 1,1
        SUM = 0
        do m = -n,n
            Sum = Sum+Yt%get(n,m)*Ys%get(n,-m)
        end do
        write(*,*) "Sum:", sum
    end do

    

    write(*,'(a,f10.6,a,f10.4)') "Exp method:", P1, ", actual:", P2

end program