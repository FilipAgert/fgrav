module fmm
    use constants
    use tree_module
    use mulpol, only: eval_grad_mulpexp, Gmulpol
    implicit none
    private
    public :: eval_acc


    contains

    function eval_acc(pos, weights, G) result(acc)
        real(kind), dimension(:,:), intent(inout) ::pos!!in: position in world coorindates. Array will change order
        real(kind), dimension(:), intent(inout) :: weights !!in: weights in world coordinates. Array will change order
        real(kind), dimension(3,size(weights,1)) :: acc!!acceleration in cartesian world coordinates corresponding to positions.
        real(kind) :: scaling, translation(3)
        real(kind), intent(in) :: G !!gravitational constant
        type(tree), pointer :: root
        Gmulpol = G
        call setup(scaling, translation, pos, weights)!rescale from po
        call root_constructor(root, pos, weights)
        call root%split()
        call eval_tree_acc(root,acc)



        !!take coordinates back to world coordinates
        pos = (pos - spread(translation,2,ncopies=size(weights,1)))/scaling
        acc = acc/scaling
        Gmulpol = Gmulpol/scaling**3
    end function

    subroutine setup(scaling, translation, pos, weights)
        real(kind), dimension(:,:), intent(inout) ::pos!!in: position in world coorindates. Out: in [0,1] coordinates
        real(kind), dimension(:), intent(inout) :: weights !!in: weights in world coordinates. Out: in [0,1] coordinates
        real(kind), intent(out) :: translation(3) !!Shift to bring world coords to [0,1]
        real(kind), intent(out) :: scaling !!scaling used to bring world coords to [0,1]
        type(tree), pointer :: root
        call world_to_tree(translation, scaling, pos)
        pos= pos* scaling + spread(translation,2,ncopies=size(weights,1))
        Gmulpol =Gmulpol*scaling**3 !Gravitational constant has units [L^3] so it needs to be scaled appropriately.
        call root_constructor(root, pos, weights)
        call root%split() !Subdivide root.
    end subroutine



end module