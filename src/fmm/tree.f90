module tree_module
    use constants, only: kind
    use fmm, only:cluster
    implicit none
    private
    integer, parameter :: clust_ceil = 10 !! 10 particles per cluster
    public world_to_tree, subdivide

    type tree
        real(kind) :: bounds(3,2)!!xyz minmax values
        contains


    end type
    contains

    subroutine world_to_tree(translation, scaling, coords) !Convert world coordinates to [0,1]^3
        real(kind), intent(in) ::coords(:,:)
        real(kind) :: bounds(3,2)
        real(kind),intent(out) :: translation(3), scaling !!translation and scaling to put coordinates from world to [0,1]
        bounds(1,1) = minval(coords(1,:),1)
        bounds(1,2) = maxval(coords(1,:),1)
        bounds(2,1) = minval(coords(2,:),1)
        bounds(2,2) = maxval(coords(2,:),1)
        bounds(3,1) = minval(coords(3,:),1)
        bounds(3,2) = maxval(coords(3,:),1)

        scaling = 1.0_kind/maxval(bounds(:,2)-bounds(:,1),1)![bounds(3,2)-bounds(3,1), bounds(2,2)-bounds(2,1), bounds(1,2)-bounds(1,1)])
        !multiply world coordinates by scaling to convert world into [a,b] where b-a = 1
        translation = -bounds(:,1)*scaling
    end subroutine

    !!Subdivide a cube into 8 subcubes
    function subdivide(bounds) result(subs)
        real(kind), intent(in) :: bounds(3,2)
        real(kind) :: subs(3,2,8)
        real(kind) :: mids(3)!division midpoints.
        integer :: i
        integer :: shft(3)
        real(kind) ::ps(3,3) !x,y,z for L,M,U
        mids = (bounds(:,2)+bounds(:,1))/2.0_kind
        ps(:,1) = bounds(:,1)!lb
        ps(:,2) = mids       !mid
        ps(:,3) = bounds(:,2)!ub

        do i =1,8
            !The ISHFTS and IAND are truth tables for 3 variables. x rightmost, y, then z
            !   0 0 1
            !   0 1 0
            !   0 1 1 etc
            shft = [IAND(i-1,1), ISHFT(IAND(i-1,2),-1), ISHFT(IAND(i-1,4),-2)] !figure out if x,y,z are lb:mid or mid:ub.
            !write(*,'(4I4)') i, u
            subs(1,:,i) = ps(1,shft(1)+1:shft(1)+2) 
            subs(2,:,i) = ps(2,shft(2)+1:shft(2)+2) 
            subs(3,:,i) = ps(3,shft(3)+1:shft(3)+2) 
        end do
    end function


end module