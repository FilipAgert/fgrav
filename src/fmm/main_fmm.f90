program main
    use constants
    use fmm


    implicit none
    type(cluster) :: cl
    real(kind) :: pos(3), weight, r(3), P1,P2, a(3),x,y,z
    complex(kind):: SUM
    type(sph_harm_coeff_c) :: Yt, Ys
    integer :: n,m
    pos = [1.0_kind, 2.0_kind, pi/4]
    weight = 1.0_kind
    allocate(cl%pos(3,1), cl%weights(1))
    cl%pos(:,1) = pos
    cl%weights = weight
    call calc_mulp_exp(cl)

    

    r = [5.0_kind, 3.0_kind, pi/5]
    P1 = eval_mulp_exp_c(cl, r)
    a= eval_grad_mulpexp(cl%mp_exp,r)

    P2 = 1.0_kind/&
    sqrt(r(1)**2 + pos(1)**2 - 2*r(1)*pos(1)*(sin(r(2))*sin(pos(2))*cos(r(3)-pos(3))  + cos(r(2))*cos(pos(2))))


    
    write(*,'(a,f10.6,a,f10.6)') "Potential: Exp method:", P1, ", actual:", P2
    write(*,'(a,e16.8,a,e16.8)') "Acceleration: Exp method:", sqrt(sum(a*a)), ", actual:", P2*P2

end program