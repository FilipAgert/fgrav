module body_module
    use constants
    implicit none

    type :: body
        character(len=20) :: name
        real(kind) :: mass ! in solar masses
        real(kind) :: radius ! in pc
        real(kind) :: position(3) ! in pc
        real(kind) :: velocity(3) ! in km/s
        logical :: negligble = .false.
    contains
        procedure :: initialize => init_body
    end type body

    contains
        subroutine init_body(this, name, mass, radius, neg)
            class(body), intent(inout) :: this
            character(len=*), intent(in) :: name
            real(kind), intent(in) :: mass, radius
            logical, optional :: neg
            logical :: acneg
            acneg = .false.
            if(present(neg)) acneg =neg

            this%name = name
            this%mass = mass
            this%radius = radius
            this%position = 0
            this %velocity = 0
            this%negligble = acneg
        end subroutine init_body

        subroutine update(this, dx, dv) !update position and velocity of body
            class(body), intent(inout) ::this
            real(kind), intent(in), dimension(3) :: dx, dv
            this%velocity = this%velocity + dv
            this%position = this%position + dx
        end subroutine

end module