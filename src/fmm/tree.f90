module tree_module
    use constants, only: kind
    use mulpol, only:cluster, calc_mulp_exp, calc_mulp_gl_exp, eval_grad_mulpexp, Gmulpol
    implicit none
    private
    integer, parameter :: clust_ceil = 20 !! particles per cluster
    integer, parameter :: maxnumfar = 200, maxnumclose = 27
    
    public world_to_tree, subdivide, tree, print, root_constructor, eval_tree_acc

    !When we create a subcube => the subcubes must recursively walk the tree
    !and get pointers to all nodes in the tree with near field
    !and all nodes with far fields
    type :: tree
        real(kind) :: bounds(3,2)!!xyz minmax values
        real(kind) :: width !!Length of box sides.
        integer :: numChild = 0 !!Number of chilren and their children etc
        type(treeptr), allocatable :: child_arr(:)
        type(tree), pointer :: parent => null()
        type(cluster), allocatable :: clust
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

    subroutine root_constructor(this, pos, weights)
        type(tree), pointer :: this
        real(kind), dimension(:,:) ::pos
        real(kind), dimension(:) :: weights
        allocate(this)
        this%bounds(:,1) = 0
        this%bounds(:,2) = 1
        this%width = 1
        allocate(this%clust)
        allocate(this%clust%pos(3,size(weights)))
        allocate(this%clust%weights(size(weights)))
        this%clust%cluster_pos =0.5
        this%clust%pos = pos
        this%clust%weights = weights
        this%clust%startidx=1
        this%clust%stopidx = size(weights) + 1
    end subroutine

    logical function isEmpty(this)
    !!Checks if node has particles
        class(tree), intent(in) :: this
        isEmpty = .not. allocated(this%clust)
        if(.not. isEmpty) then
            isEmpty = this%clust%startidx == this%clust%stopidx
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
        integer :: ii, index
        real(kind) :: subs(3,2,8)
        real(kind) :: mids(3)
        if(this%clust%stopidx-this%clust%startidx > clust_ceil ) then
            subs = subdivide(this%bounds)
            this%numChild = 8
            if(.not. this%isRoot()) this%parent%numChild = this%parent%numChild + 8
            allocate(this%child_arr(8))
            mids = (this%bounds(:,2) + this%bounds(:,1))/2.0_kind
            call quicksort_octant(this%clust%pos, this%clust%weights, mids, this%clust%startidx,this%clust%stopidx-1)
            index = this%clust%startidx
            do ii = 1,8
                allocate(tree :: this%child_arr(ii)%child)
                this%child_arr(ii)%child%parent => this
                this%child_arr(ii)%child%bounds = subs(:,:,ii)
                this%child_arr(ii)%child%width = this%width*0.5_kind
                call assign_clust(this%child_arr(ii)%child, this%clust, index)
                if(this%child_arr(ii)%child%isEmpty()) cycle !If child has no points, dont compute anything.
                !call shrinkwrap(this%child_arr(ii)%child) !!Reduce bounds of child to smallest possible
                call calc_clust(this%child_arr(ii)%child)!!Precompute the multipole expansion for this cluster.
                call split(this%child_arr(ii)%child)
            end do
        else
            !write(*,*) "Didnt split, size of weights:", size(this%clust%weights,1)
        endif
        !else, do nothing as this node only has a few particles inside. Dont need to subdivide further.
    end subroutine

    recursive subroutine quicksort_octant(pos, weights, mids, idxstart, idxstop)
        real(kind), pointer :: pos(:,:)   ! shape (3, N)
        real(kind), pointer :: weights(:) ! length N
        real(kind), intent(in) :: mids(3) ! midpoints for octant calc
        integer, intent(in) :: idxstart, idxstop
        integer :: i, j, pivot, tmpOct
        real(kind) :: tmpW
        real(kind) :: tmpPos(3)

        if (idxstart >= idxstop) return

        ! Choose pivot = octant of middle element
        pivot = octant(pos, (idxstart + idxstop) / 2, mids)

        i = idxstart
        j = idxstop

        do
            do while (octant(pos, i, mids) < pivot)
                i = i + 1
            end do
            do while (octant(pos, j, mids) > pivot)
                j = j - 1
            end do

            if (i <= j) then
                ! Swap positions
                tmpPos = pos(:, i)
                pos(:, i) = pos(:, j)
                pos(:, j) = tmpPos

                ! Swap weights
                tmpW = weights(i)
                weights(i) = weights(j)
                weights(j) = tmpW

                i = i + 1
                j = j - 1
            end if

            if (i > j) exit
        end do

        if (idxstart < j) call quicksort_octant(pos, weights, mids, idxstart, j)
        if (i < idxstop)  call quicksort_octant(pos, weights, mids, i, idxstop)

    end subroutine

    !!gets which octant coordinate is in based on midpoints of bounds
    integer function octant(pos, idx, mids)
        real(kind), pointer :: pos(:,:)
        integer, intent(in) :: idx
        real(kind), intent(in) ::mids(3)
        integer :: bits(3)
        bits = merge(1,0,pos(:,idx) >= mids)
        octant = 1 + bits(1) + 2*bits(2) + 4*bits(3)
        !

    end function
    recursive function getRoot(this) result(root)
        type(tree), intent(in), pointer :: this
        type(tree), pointer :: root

        if(this%isRoot()) then
            root => this
        else
            root => getRoot(this%parent)
        endif
    end function

    !!In the tree, adds all nodes which are nearfield to nears and all nodes which are far field to fars.
    !!Excludes self in both of these since we dont want selfinteraction for a particle.
    subroutine startDescentNearFar(origin, nears, nearIdx, fars, farsIdx)
        type(tree), intent(in), pointer :: origin
        type(treeptr), intent(out) :: nears(maxnumclose), fars(maxnumfar)
        integer, intent(out) :: nearIdx, farsIdx
        integer :: numnodes
        type(tree), pointer :: root
        root => getRoot(origin)
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
            !if(farsIdx > maxnumfar) error stop "too many fars"
            fars(farsIdx)%child => node
            !Stop descending if farfield as all children will be farfield
        else!near field
            if(node%isLeaf()) then !if its a leaf, add as near field
                nearIdx = nearIdx + 1
                !if(nearIdx > maxnumclose) error stop "too many close"
                nears(nearIdx)%child => node
            else!if not a leaf, one of the children might be far field, so recursively call this function.
                do ii = 1,8
                    call descentAddNearFar(origin, node%child_arr(ii)%child, nears, nearIdx, fars, farsIdx)
                end do
            endif
        endif

    end subroutine

    recursive subroutine eval_tree_acc(node, acc)
    !!evaluates acceleration for all subnodes of input node.
        type(tree), pointer :: node !!in: root of tree. 
        real(kind), dimension(:,:), intent(inout) :: acc
        integer :: ii
        if(node%isLeaf()) then
            if(node%isEmpty()) return
            !write(*,*) "Not empty leaf."
            call eval_node_acc(node, acc)
        else
            if(node%isEmpty()) then
                !write(*,*) "Empty node."
                return
            else !has children
                !write(*,*) "do for children."
                do ii = 1,8
                    call eval_tree_acc(node%child_arr(ii)%child,acc)
                end do
            endif
        endif
    end subroutine

    subroutine eval_node_acc(node, acc)
        type(tree), pointer :: node
        real(kind), dimension(:,:), intent(inout) :: acc
        type(treeptr) :: nears(maxnumclose), fars(maxnumfar)
        integer :: nearIdx, farIdx, ii
        call startDescentNearFar(node, nears, nearIdx, fars, farIdx) !!gets interaction list for node.
        do ii = 1, nearIdx
            call near_interaction(node, nears(ii)%child,acc)
        end do
        do ii = 1, farIdx
            call far_interaction(node, fars(ii)%child, acc)
        end do  
        call self_interaction(node, acc)

    end subroutine

    subroutine near_interaction(node1, node2, acc)
        type(tree), pointer:: node1, node2
        real(kind), dimension(:,:), intent(inout) :: acc
        real(kind), dimension(3) :: locAcc
        real(kind) :: r(3)
        integer :: ii,jj
        locAcc = 0
        do ii = node1%clust%startidx, node1%clust%stopidx -1!iterate over all particles in node1
            do jj = node2%clust%startidx, node2%clust%stopidx - 1
                r = node2%clust%pos(:,jj) - node1%clust%pos(:,ii) !vector from node1 to node2. Acceleration on node1 particle.
                locAcc = r*Gmulpol/(sum(r*r,1)**(3.0/2.0)) 
                acc(:,ii) = acc(:,ii) + locAcc*node2%clust%weights(jj)!!acceleration on node1 particle.
                acc(:,jj) = acc(:,jj) - locAcc*node1%clust%weights(ii)!!acceleration on node2 particle.
            end do
        end do
    end subroutine

    subroutine self_interaction(node, acc)
        type(tree), pointer:: node
        real(kind), dimension(:,:), intent(inout) :: acc
        real(kind), dimension(3) :: locAcc
        real(kind) :: r(3)
        integer :: ii,jj
        locAcc = 0
        do ii = node%clust%startidx, node%clust%stopidx -1!iterate over all particles in node1
            do jj = ii + 1, node%clust%stopidx -1
                r = node%clust%pos(:,jj) - node%clust%pos(:,ii) !vector from p1 to p2. Acceleration on p1 particle.
                locAcc = r*Gmulpol/(sum(r*r,1)**(3.0/2.0)) 
                acc(:,ii) = acc(:,ii) + locAcc*node%clust%weights(jj)!!acceleration on p1 particle.
                acc(:,jj) = acc(:,jj) - locAcc*node%clust%weights(ii)!!acceleration on p2 particle.
            end do
        end do
    end subroutine

    subroutine far_interaction(node1, node2, acc) 
        !!acceleration to particles in node1 from node2
        type(tree), pointer ::node1, node2
        real(kind), dimension(:,:), intent(inout) :: acc
        real(kind) :: g(3)
        integer :: ii
        do ii = node1%clust%startidx, node1%clust%stopidx -1!iterate over all particles in node1
            g = eval_grad_mulpexp(node2%clust%mp_gl_exp, node1%clust%pos(:,ii))
            acc(:,ii) = acc(:,ii) + g/node1%clust%weights(ii)
        end do


    end subroutine

    pure logical function isNearField(this, other)
        type(tree), intent(in) :: this, other
        real(kind) :: centerdist
        real(kind), parameter :: distToCorner = sqrt(3.0_kind)/2.0_kind !distance to corner from center
        !of a unit cube
        centerdist = sqrt(sum((this%clust%cluster_pos-other%clust%cluster_pos)**2,1))
        isNearField = centerdist < (this%width + other%width)*distToCorner
    end function

    !!From a parents cluster, take relevant points into this cluster's 
    subroutine assign_clust(this, parent_clust, idx)
        type(tree), intent(inout) :: this
        type(cluster), intent(in) ::parent_clust
        integer, intent(inout) :: idx !on input: index pointing to first place in positions array. on output:
        !where next child needs to start searching
        integer :: ctr, ii
        real(kind) :: shift(3)
        logical :: allocated
        do
            if(contains_point(this%bounds, parent_clust%pos(:,idx))) then
                allocate(this%clust)
                this%clust%startidx = idx
                this%clust%pos => parent_clust%pos
                this%clust%weights => parent_clust%weights
                this%clust%cluster_pos = (this%bounds(:,2) + this%bounds(:,1) )/2.0_kind
                idx = idx + 1
                exit
            else
                return
            endif
        end do

        do
            if(contains_point(this%bounds, parent_clust%pos(:,idx))) then
                idx = idx + 1
            else
                exit
            endif
        end do
        idx = idx !exclusive
        this%clust%stopidx = idx
    end subroutine

    !!Adjust bounds of node to smallest possible
    subroutine shrinkwrap(this)
        type(tree), intent(inout) ::this
        real(kind) :: width, bounds(3,2), slide(3,2), llup, uld
        integer :: a,b
        a = this%clust%startidx
        b = this%clust%stopidx
        this%bounds(:,1) = minval(this%clust%pos(:,a:b-1),2) !minvals for x,y,z
        this%bounds(:,2) = maxval(this%clust%pos(:,a:b-1),2) !maxvals for x,y,z
        slide(:,1) = bounds(:,1)-this%bounds(:,1)
        slide(:,2) = this%bounds(:,2)-bounds(:,2)
        llup = minval(slide(:,1),1)
        uld = minval(slide(:,2),1)
        this%bounds(:,1) = this%bounds(:,1) + llup
        this%bounds(:,2) = this%bounds(:,2) - uld
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

    subroutine world_to_tree(translation, scaling, coords) 
        !!Convert world coordinates to [0,1]^3
        !!To use: converted coordinates are cords*scaling + translation
        real(kind), intent(in) ::coords(:,:)
        real(kind) :: bounds(3,2)
        real(kind),intent(out) :: translation(3)!!translation to shift scaled coordinates to [0,1]
        real(kind), intent(out) :: scaling !!scaling to scale coordinates to [a,b] where b-a = 1
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

    recursive subroutine dealloc(this)
        type(tree), intent(inout), pointer :: this
        integer :: ii
        if(allocated(this%clust)) deallocate(this%clust)
        if(isLeaf(this)) then

        else
            do ii = 1,8
                call dealloc(this%child_arr(ii)%child)
            end do
        endif
    end subroutine

    recursive subroutine print(this, level)
        type(tree), pointer ::this
        integer, intent(in)::level
        integer :: ii
        write(*,*) 
        write(*,'(a,i3)') "Level:", level
        write(*,'(a,i4)') "number of particles: ", this%clust%stopidx-this%clust%startidx
        write(*,'(a,f5.3,1x,f5.3,a,f5.3,1x,f5.3,a,f5.3,1x,f5.3)') "x:", this%bounds(1,:), ", y:", this%bounds(2,:), ", z:", this%bounds(3,:)
        do ii = this%clust%startidx, this%clust%stopidx-1
            write(*,'(a,f5.3,1x,f5.3,1x,f5.3,a)') "[",this%clust%pos(:,ii),"]"
        end do

        if(.not. this%isLeaf()) then
            do ii = 1, 8
                if(.not. this%child_arr(ii)%child%isEmpty()) call print(this%child_arr(ii)%child,level + 1)
            end do

        endif

        
    end subroutine
end module