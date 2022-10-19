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
        real(kind=8)             :: mass       ! Total mass of the structure
        type(scale), allocatable :: cascade(:) ! Gas structuration cascade, set of scale
    contains
        procedure  :: create => gsh_create     ! Create the gas structuration (cascade)
        procedure  :: delete => gsh_delete     ! Delete the gas structuration (cascade)
        procedure  :: copy => gsh_copy         ! Copy a gsh object
        procedure  :: evolve => gsh_evolve     ! Evolve the complete gsh struture by dt
        procedure  :: status => gsh_status     ! Return the current status (outRate) of the gsh
        procedure  :: solve => gsh_solve       ! Solve the input output system during dt
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

        integer(kind=4)  :: i

        real(kind=8)     :: l
        class(gsh)       :: this

        this%mass = 0.d0                ! Total mass of the structure
        allocate(this%cascade(nScales)) ! Create the cascade with nScales cells

        ! Create gas element in each cells
        do i = 1, nScales
            l = gsh_i2l(i)
            call this%cascade(i)%create(l)
        end do

    end subroutine gsh_create

    ! **********************************
    subroutine gsh_delete(this)

        ! Delete a gas structuration history structure

        implicit none

        integer(kind=4)  :: il
        class(gsh)       :: this

        this%mass = 0.d0     ! Total mass of the structure

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

        integer(kind=4)                  :: s

        type(gsh), intent(in)            :: gs2

        class(gsh), intent(inout)        :: gs1

        ! Copy mass
        gs1%mass = gs2%mass
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

        integer(kind=4)          :: s

        real(kind=8), intent(in) :: dt           ! The gsh is evolved during dt
        real(kind=8), intent(in) :: l            ! Scale injection

        type(gas), intent(in)    :: inRate       ! The (dt-)constant input rate
        type(gas), allocatable   :: outRates(:)  ! The corrected set of output rates
        type(gas), allocatable   :: outRates1(:) ! Intermediates set of output rates
        type(gas), allocatable   :: outRates2(:)
        type(gas), allocatable   :: outRates3(:)
        type(gas), allocatable   :: outRates4(:)

        type(gsh)                :: gsh_tmp  ! Intermediate states of the scale

        class(gsh)               :: this     ! The current scale

        ! Init OutRates
        allocate(outRates(nScales))
        !
        ! Get curent status
        outRates1 = this%status()
        ! Performed evolution for dt/2
        gsh_tmp = this%evolve(dt/2.d0, l, inRate, outRates1)
        !
        select case (trim(intScheme))
            !
            case ('RK4')
                ! Range-Kutta 4th order
                ! Get curent status
                outRates2 = gsh_tmp%status()
                ! Get second intermediate state
                gsh_tmp = this%evolve(dt/2.d0, l, inRate, outRates2)
                outRates3 = gsh_tmp%status()
                ! Get last intermediate state
                gsh_tmp = this%evolve(dt, l, inRate, outRates3)
                outRates4 = gsh_tmp%status()
                ! Perform complete (corrected) evolution
                do s = 1, nScales
                    outRates(s) = 1.d0/6.d0*(outRates1(s) + &
                                            2.d0*outRates2(s) + &
                                            2.d0*outRates3(s) + &
                                            outRates4(s))
                end do
                ! Final evolution with corrected output rates
                this = this%evolve(dt, l, inRate, outRates)
            case default
                ! Range-Kutta 2d order
                ! Get intermedite status
                outRates = gsh_tmp%status()
                ! Perform complete evolution
                this = this%evolve(dt, l, inRate, outRates)
        end select
        !
        ! Delete temporary gas object
        do s = 1, nScales
            call outRates(s)%delete()
            call outRates1(s)%delete()
            if (allocated(outRates2)) call outRates2(s)%delete()
            if (allocated(outRates3)) call outRates3(s)%delete()
            if (allocated(outRates4)) call outRates4(s)%delete()
        end do
        ! Deallocate arrays
        deallocate(outRates, outRates1)
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

        integer(kind=4), intent(in)  :: i

        real(kind=8)                 :: l

        l = stepFactor**(i-1)*lStar

    end function gsh_i2l

    ! **********************************
    function gsh_l2i(l) result(i)

        ! Return the index (in the cascade array) associated to the scale l

        implicit none

        integer(kind=4)  :: i

        real(kind=8)     :: l    ! The scale

        i =  int(log(l/lStar)/log(stepFactor), kind=4) + 2

    end function gsh_l2i

    ! **********************************
    function gsh_status(this) result(rates)

        implicit none

        integer(kind=4)        :: s

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
    function gsh_evolve(this, dt, l, inRate, outRates) result(gs)

        ! Evolve a gas structuration history by dt
        ! The gas structuration cascade is fed
        ! at the scale "l" by an input rate "inRate"
        ! and evolve trought a set of output rates "outRates"

        implicit none
    
        integer(kind=4)                     :: s
        integer(kind=4)                     :: is

        character(MAXPATHSIZE)              :: calledBy

        real(kind=8), intent(in)            :: l           ! The injection scale
        real(kind=8), intent(in)            :: dt      

        type(gas), intent(in)               :: inRate      ! The input rate
        type(gas), intent(in), allocatable  :: outRates(:) ! The set of output rates
        type(gas)                           :: sInRate     ! Total input rate at scale s

        class(gsh)                          :: this        ! The current gsh
        type(gsh)                           :: gs          ! The evolved gsh

        write(calledBy, '(a)') 'gsh_evolve'

        ! Init gs structure from current
        gs = this
        ! Init sInRate
        call sInRate%create()
        ! 
        ! Gas injection will be performed at scale indexed is
        is = gsh_l2i(l)
        !
        ! Loop over scales from the largest to the lowest
        do s = nScales, 1, -1
            !
            ! Largest scale case, no input from larger scale, only external input
            if (s < nScales) then
                call sInRate%copy(outRates(s+1))  ! Input from larger scale
            end if
            !
            ! External injection at specific scale
            if (s == is) then
                sInRate = sInRate + inRate
            end if
            !
            ! Apply evolution at scale s
            gs%cascade(s) = this%cascade(s)%evolve(dt, sInRate, outRates(s))
            !
        end do
        !
        ! sInRate is delete
        call sInRate%delete()

    end function gsh_evolve

end module gsh_mod