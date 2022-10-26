module gsh_mod

    !*****************************************************************************************************************
    !
    !  OVERVIEW
    !
    ! The gsh module defines all properties and procedures
    ! associated to gas structuration history structure
    ! This structure allows to "follow" gas structuration process
    ! to the different binned scales
    !
    !*****************************************************************************************************************

    use parameters  ! Acces to global defintions and properties
    use model_mod   ! Acces to model parameters and integration scheme configuration
    use scale_mod   ! Acces to scale properties and procedures

    implicit none

    public

    type gsh
        real(kind=rkd)             :: mass          ! Total mass of the structure
        real(kind=rkd)             :: l             ! Current largest scale (at least one Cloud)
        type(scale), allocatable :: cascade(:)    ! Gas structuration cascade, set of scale
    contains
        procedure  :: create => gsh_create        ! Create the gas structuration (cascade)
        procedure  :: delete => gsh_delete        ! Delete the gas structuration (cascade)
        procedure  :: copy => gsh_copy            ! Copy a gsh object
        procedure  :: evolve => gsh_evolve        ! Evolve the complete gsh struture by dt
        procedure  :: status => gsh_status        ! Return the current status (outRate) of the gsh
        procedure  :: solve => gsh_solve          ! Solve the input output system during dt
        procedure  :: Vesc => gsh_escape_velocity ! Return the escape velocity of the gs
    end type gsh

    ! INTERFACE OPERATOR DECLARATIONS

    interface assignment (=)  ! allows to copy a gsh component by using the symbol '='
        module procedure gsh_copy_
    end interface assignment (=)

contains

    !
    ! SUBROUTINES
    !

    ! **********************************
    subroutine gsh_init()

        ! Initialize the gsh module and dependancies

        implicit none

        call scale_init()

    end subroutine gsh_init

    ! **********************************
    subroutine gsh_create(this)

        ! Create (allocate) a gas hstory structure

        implicit none

        integer(kind=ikd)  :: i

        real(kind=rkd)     :: l
        class(gsh)       :: this

        this%mass = 0.d0                ! Total mass of the structure
        allocate(this%cascade(nScales)) ! Create the cascade with nScales cells

        ! Create gas element in each cells
        do i = 1, nScales
            l = gsh_i2l(i)
            call this%cascade(i)%create(l)
        end do

        ! Set the current largest scale
        this%l = this%cascade(nScales)%l
        i = gsh_l2i(this%l)

        return
    end subroutine gsh_create

    ! **********************************
    subroutine gsh_delete(this)

        ! Delete a gas structuration history structure

        implicit none

        integer(kind=ikd)  :: il
        class(gsh)       :: this

        this%mass = 0.d0     ! Total mass of the structure
        this%l    = 0.d0     ! Highest occupied scale

        ! Delete gas element in each cells
        do il = 1, nScales
            call this%cascade(il)%delete()
        end do

        deallocate(this%cascade) ! Delete the cascade

    end subroutine gsh_delete

    ! **********************************
    subroutine gsh_copy(gs1, gs2)

        ! Copy the gsh object gs2 into the gsh object gs1

        implicit none

        integer(kind=ikd)                  :: s

        type(gsh), intent(in)            :: gs2

        class(gsh), intent(inout)        :: gs1

        ! Copy
        gs1%mass = gs2%mass
        gs1%l    = gs2%l
        ! 
        ! Allocate cascade
        if (.not. allocated(gs1%cascade)) then
            ! Create the cascade with nScales cells
            allocate(gs1%cascade(nScales))
        end if
        ! Copy each scale
        do s = 1, nScales
            gs1%cascade(s) = gs2%cascade(s)
        end do
    
    end subroutine gsh_copy

    ! **********************************
    subroutine gsh_copy_(gs1, gs2)

        ! Interface procedure to copy

        implicit none

        type(gsh), intent(in)      :: gs2
        class(gsh), intent(inout)  :: gs1

        call gsh_copy(gs1, gs2)

    end subroutine gsh_copy_

    ! **********************************
    subroutine gsh_solve(this, dt, inRate, l)

        ! Solve the evolution of the current gas structuration
        !  according to a constant external input rate "inRate"

        implicit none

        integer(kind=ikd)          :: s
        integer(kind=ikd)          :: il

        real(kind=rkd), intent(in) :: dt         ! The gsh is evolved during dt
        real(kind=rkd), intent(in) :: l          ! Scale injection

        type(gas), intent(in)    :: inRate       ! The (dt-)constant input rate
        type(gas), allocatable   :: inRates(:)   ! The set of input rates
        type(gas), allocatable   :: outRates(:)  ! The corrected set of output rates
        type(gas), allocatable   :: outRates1(:) ! Intermediates set of output rates
        type(gas), allocatable   :: outRates2(:)
        type(gas), allocatable   :: outRates3(:)
        type(gas), allocatable   :: outRates4(:)

        type(gsh)                :: gsh_tmp  ! Intermediate states of the scale

        class(gsh)               :: this     ! The current scale

        ! Init outRates
        allocate(outRates(nScales))
        allocate(inRates(nScales))
        !
        ! Init inRates and add external input rates
        il = gsh_l2i(l)
        do s = 1, nScales
            call inRates(s)%create()
            if (s == il) inRates(s) = inRate
        end do
        !
        ! Get curent status
        outRates1 = this%status()
        ! Performed evolution for dt/2
        gsh_tmp = this%evolve(dt/2.d0, inRates, outRates1)
        !
        select case (trim(intScheme))
            !
            case ('RK4')
                ! Range-Kutta 4th order
                ! Get curent status
                outRates2 = gsh_tmp%status()
                ! Get second intermediate state
                gsh_tmp = this%evolve(dt/2.d0, inRates, outRates2)
                outRates3 = gsh_tmp%status()
                ! Get last intermediate state
                gsh_tmp = this%evolve(dt, inRates, outRates3)
                outRates4 = gsh_tmp%status()
                ! Perform complete (corrected) evolution
                do s = 1, nScales
                    outRates(s) = real(1.d0/6.d0,kind=rkd)*(outRates1(s) + &
                                                            real(2.d0, kind=rkd)*outRates2(s) + &
                                                            real(2.d0, kind=rkd)*outRates3(s) + &
                                                            outRates4(s))
                end do
                ! Final evolution with corrected output rates
                this = this%evolve(dt, inRates, outRates)
            case default
                ! Range-Kutta 2d order
                ! Get intermedite status
                outRates = gsh_tmp%status()
                ! Perform complete evolution
                this = this%evolve(dt, inRates, outRates)
        end select
        !
        ! Delete temporary gas object
        do s = 1, nScales
            call inRates(s)%delete()
            call outRates(s)%delete()
            call outRates1(s)%delete()
            if (allocated(outRates2)) call outRates2(s)%delete()
            if (allocated(outRates3)) call outRates3(s)%delete()
            if (allocated(outRates4)) call outRates4(s)%delete()
        end do
        ! Deallocate arrays
        deallocate(inRates, outRates, outRates1)
        if (allocated(outRates2)) deallocate(outRates2)
        if (allocated(outRates3)) deallocate(outRates3)
        if (allocated(outRates4)) deallocate(outRates4)

    end subroutine gsh_solve

    !
    ! FUNCTIONS
    !

    ! **********************************
    function gsh_i2l(i) result(l)

        ! Return the ith scale of the cascade

        implicit none

        integer(kind=ikd), intent(in)  :: i

        real(kind=rkd)                 :: l

        l = stepFactor**(i-1)*lStar

    end function gsh_i2l

    ! **********************************
    function gsh_l2i(l) result(i)

        ! Return the index (in the cascade array) associated to the scale l

        implicit none

        integer(kind=ikd)  :: i

        real(kind=rkd)     :: l    ! The scale

        i =  int(log(l/lStar)/log(stepFactor), kind=ikd) + int(1, kind=ikd)
        i = min(i, nScales)

    end function gsh_l2i

    ! **********************************
    function gsh_escape_velocity(this) result(Vesc)

        ! Return the escape velovity of the current gas structuration
        ! The total mass (gsh%mass) of the gas is distributed into Nclouds formed
        ! at the largest occupied scale (gsh%l)

        implicit none

        integer(kind=ikd)    :: il

        real(kind=rkd)       :: Vesc
        real(kind=rkd)       :: NClouds
        real(kind=rkd)       :: mass

        class(gsh)         :: this

        il = gsh_l2i(this%l)  ! Get the index of the largest occupied scale
        
        ! The mass is distributed in the NClouds formed 
        ! at the largest occupied scale
        NClouds = this%cascade(il)%NClouds()
        Vesc = 0.d0 ! Init
        if (NClouds > 0.d0) then
            mass = this%mass / NClouds
            !
            Vesc = sqrt(2.d0*GCst_CU*mass/this%l)  ! [CU]
        end if

    end function gsh_escape_velocity

    ! **********************************
    function gsh_status(this) result(rates)

        implicit none

        integer(kind=ikd)        :: s

        type(gas), allocatable :: rates(:)  ! The set of output rates (all scales)
        
        class(gsh)             :: this
        !
        ! Create the set of output rates
        allocate(rates(nScales))
        ! Set values
        do s = 1, nScales
            rates(s) = this%cascade(s)%status()
        end do

    end function gsh_status

    ! **********************************
    function gsh_evolve(this, dt, inRates, outRates) result(gs)

        ! Evolve a gas structuration history by dt
        ! The gas structuration cascade evolves with:
        !  - A set of input rates "inRates", 
        !    (external input at the scale "l" is already included)
        ! - A set of output rates "outRates"

        implicit none
    
        integer(kind=ikd)                      :: s

        real(kind=rkd), intent(in)          :: dt          ! The time-step
        real(kind=rkd)                      :: mass        ! The updated total mass of the gsh
        real(kind=rkd)                      :: ll          ! The largest occupied scale

        type(gas), intent(in), allocatable  :: inRates(:)  ! The input rate
        type(gas), intent(in), allocatable  :: outRates(:) ! The set of output rates
        type(gas)                           :: sInRate     ! local input rate

        type(gsh)                           :: gs          ! The evolved gsh

        class(gsh)                          :: this        ! The current gsh

        ! Init gs structure from current
        gs = this
        ! Init sInRate
        call sInRate%create()
        !
        ! Loop over scales from the largest to the lowest
        mass = 0. ! Init the total mass
        do s = nScales, 1, -1
            !
            ! Take into account internal transfer rates
            ! The largest scale if only fed by external input
            if (s < nScales) then
                sInRate = InRates(s) + outRates(s+1)
            end if
            ! Apply evolution at scale s
            gs%cascade(s) = this%cascade(s)%evolve(dt, sInRate, outRates(s))
            !
            ! Update total mass
            mass = mass + gs%cascade(s)%gas%mass
            !
            ! Update the largest occupied scale
            if ((mass > 0.) .and. (ll < gs%cascade(s)%l)) then
                ll = gs%cascade(s)%l
            end if
        end do
        !
        ! Set new total gas mass
        gs%mass = mass
        ! Set new largest occupied scale
        gs%l = ll

    end function gsh_evolve

end module gsh_mod