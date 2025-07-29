module tree_module
    use constants, only: kind
    use fmm, only:cluster, calc_mulp_exp, calc_mulp_gl_exp
    implicit none
    private
    integer, parameter :: clust_ceil = 10 !! 10 particles per cluster
    
    public world_to_tree, subdivide, tree

    !When we create a subcube => the subcubes must recursively walk the tree
    !and get pointers to all nodes in the tree with near field
    !and all nodes with far fields
    type :: tree
        real(kind) :: bounds(3,2)!!xyz minmax values
        real(kind) :: width !!Length of box sides.
        integer :: numChild = 0 !!Number of chilren and their children etc
        type(treeptr), allocatable :: child_arr(:)
        type(tree), pointer :: parent => null()
        type(cluster) :: clust
        contains

        procedure, pass(this) :: isRoot
        procedure, pass(this) :: isLeaf
        procedure, pass(this) :: isEmpty
        procedure, pass(this) :: split
    end type

    type :: treeptr
        type(tree), pointer :: child
    end type



    contains

    
    logical function isEmpty(this)
        class(tree), intent(in) :: this
        isEmpty = .not. allocated(this%clust%pos)
        if(.not. isEmpty) then
            isEmpty = size(this%clust%weights,1) == 0
        endif
    end function

    logical function isRoot(this)
        class(tree), intent(in) ::this
        isRoot = .not. associated(this%parent)
    end function

    logical function isLeaf(this)
        class(tree), intent(in) ::this
        isLeaf = .not. allocated(this%child_arr)
    end function


    recursive subroutine split(this)!!
        class(tree), intent(inout), target :: this
        integer :: ii
        real(kind) :: subs(3,2,8)
        if(size(this%clust%weights,1) > clust_ceil ) then
            subs = subdivide(this%bounds)
            this%numChild = 8
            if(.not. this%isRoot()) this%parent%numChild = this%parent%numChild + 8
            allocate(this%child_arr(8))
            do ii = 1,8
                allocate(tree :: this%child_arr(ii)%child)
                this%child_arr(ii)%child%parent => this
                this%child_arr(ii)%child%bounds = subs(:,:,ii)
                this%child_arr(ii)%child%width = this%width*0.5_kind
                call assign_clust(this%child_arr(ii)%child, this%clust)
                call calc_clust(this%child_arr(ii)%child)!!Precompute the multipole expansion for this cluster.
                call this%child_arr(ii)%child%split()
            end do
        endif
    end subroutine

    recursive function getRoot(this) result(root)
        type(tree), intent(in), pointer :: this
        type(tree), pointer :: root

        if(this%isRoot()) then
            root => this
        else
            root => getRoot(this%parent)
        endif
    end function

    subroutine startDescentNearFar(origin, nears, nearIdx, fars, farsIdx)
        type(tree), intent(in), pointer :: origin
        type(treeptr), intent(out), allocatable :: nears(:), fars(:)
        integer, intent(out) :: nearIdx, farsIdx
        integer :: numnodes
        type(tree), pointer :: root
        root => getRoot(origin)
        numnodes = root%numChild
        allocate(nears(numNodes), fars(numnodes))
        nearIdx = 0
        farsIdx = 0
        call descentAddNearFar(origin, root, nears, nearIdx, fars, farsIdx)
    end subroutine

    recursive subroutine descentAddNearFar(origin, node, nears, nearIdx, fars, farsIdx)
        type(tree), intent(in), pointer :: node
        type(tree), intent(in), pointer :: origin
        type(treeptr), intent(inout) :: nears(:), fars(:)
        integer, intent(inout) ::nearIdx, farsIdx
        integer :: ii
        if(node%isEmpty()) return
        if(associated(origin, node)) return!origin == node. skip.
        if(.not. isNearField(origin, node)) then!far field
            farsIdx = farsIdx + 1
            fars(farsIdx)%child => node
            !Stop descending if farfield as all children will be farfield
        else!near field
            if(node%isLeaf()) then !if its a leaf, add as near field
                nearIdx = nearIdx + 1
                nears(nearIdx)%child => node
            else!if not a leaf, one of the children might be far field, so recursively call this function.
                do ii = 1,8
                    call descentAddNearFar(origin, node%child_arr(ii)%child, nears, nearIdx, fars, farsIdx)
                end do
            endif
        endif

    end subroutine

    !!Near field if this box touches other box (including diagonal)
    pure logical function isNearField(this, other)
        type(tree), intent(in) :: this, other
        logical t(3)
        t = this%bounds(:,1) <= other%bounds(:,2) .and. other%bounds(:,1) <= this%bounds(:,2)
        isNearField = ALL(t)
        !Touching if my min is smaller than their max AND their min is smaller than my max
        !e.g. for A = [0.25, 0.5] B=[0.5, 0.75]. A.min <= B.max .and. B.min <= A.max
    end function

    !!From a parents cluster, take relevant points into this cluster's 
    subroutine assign_clust(this, parent_clust)
        type(tree), intent(inout) :: this
        type(cluster), intent(in) ::parent_clust
        integer :: ctr, ii
        ctr = 0
        do ii = 1, size(parent_clust%weights,1)
            if (contains_point(this%bounds, parent_clust%pos(:,ii))) ctr = ctr + 1
        end do
        if(ctr > 0) then
            allocate(this%clust%pos(3,ctr), this%clust%weights(ctr))
            ctr = 0
            this%clust%cluster_pos = (this%bounds(:,2) + this%bounds(:,1) )/2.0_kind
            do ii = 1, size(parent_clust%weights,1)
                if (contains_point(this%bounds, parent_clust%pos(:,ii)))then
                    this%clust%pos(:,ctr) = parent_clust%pos(:,ii)
                    this%clust%weights(ctr) = parent_clust%weights(ii)
                endif
            end do
            !Positions are all relative center of cluster.
            this%clust%pos = this%clust%pos - (this%clust%cluster_pos - parent_clust%cluster_pos)
        endif
    end subroutine

    subroutine calc_clust(this)
        type(tree), intent(inout) :: this
        call calc_mulp_exp(this%clust)
        call calc_mulp_gl_exp(this%clust)
    end subroutine

    pure logical function contains_point(bounds, point)
        real(kind), intent(in) :: bounds(3,2)
        real(kind), intent(in) :: point(3)
        logical :: inside(3)
        inside = point >= bounds(:,1) .and. point < bounds(:,2)!Lower edge is included. upper edge is exlcluded.
        inside = inside .or. (point==1.0_kind .and. bounds(:,2) == 1.0_kind) !Upper limit of box edge is included.
        contains_point = ALL(inside,1)
    end function

    subroutine world_to_tree(translation, scaling, coords) !Convert world coordinates to [0,1]^3
        real(kind), intent(in) ::coords(:,:)
        real(kind) :: bounds(3,2)
        real(kind),intent(out) :: translation(3), scaling !!translation and scaling to put coordinates from world to [0,1]
        bounds(:,1) = minval(coords,2) !minvals for x,y,z
        bounds(:,2) = maxval(coords,2) !maxvals for x,y,z
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
            !   0 0 0-ä
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