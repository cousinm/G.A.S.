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
    use status_mod  ! Acces to transfer rates status
    use scale_mod   ! Acces to scale properties and procedures

    implicit none

    public

    type gsh
        real(kind=rkd)             :: mass          ! Total mass of the structure
        real(kind=rkd)             :: l             ! Current largest scale (at least one Cloud)
        type(gas)                  :: sfr           ! Instantaneous star formation rate (the tranfer rate of the lowest scale)
        type(scale), allocatable   :: cascade(:)    ! Gas structuration cascade, set of scale
    contains
        procedure  :: create => gsh_create        ! Create the gas structuration (cascade)
        procedure  :: delete => gsh_delete        ! Delete the gas structuration (cascade)
        procedure  :: copy => gsh_copy            ! Copy a gsh object
        procedure  :: solve => gsh_solve          ! Solve, for dt, the next step of the integration scheme
        procedure  :: stevolve => gsh_stevolve    ! Evolve, for dt, the next step of the integration scheme
        procedure  :: evolve => gsh_evolve        !         according to internal/external input/output
        procedure  :: Vesc => gsh_escape_velocity ! Return the escape velocity of the gs
    end type gsh

    ! TEMPORARY INTERMEDIATE STATUS
    ! Targets
    type(status), allocatable, target     :: gshStatus(:)

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

        ! Init intermediate status
        allocate(gshStatus(nScales))

    end subroutine gsh_init

    ! **********************************
    subroutine gsh_create(this)

        ! Create (allocate) a gas history structure

        implicit none

        integer(kind=ikd)  :: i

        real(kind=rkd)     :: l
        class(gsh)       :: this

        this%mass = 0.d0                ! Total mass of the structure
        call this%sfr%create()          ! Create SFR, set to null

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

        this%mass = 0.d0       ! Total mass of the structure
        this%l    = 0.d0       ! Highest occupied scale

        call this%sfr%delete() ! SFR

        ! Delete gas element in each cells
        do il = 1, nScales
            call this%cascade(il)%delete()
        end do

        deallocate(this%cascade) ! Delete the cascade

    end subroutine gsh_delete

    ! **********************************
    subroutine gsh_copy(gsh1, gsh2)

        ! Copy the gsh object gsh2 into the gsh object gsh1

        implicit none

        integer(kind=ikd)           :: s

        type(gsh), intent(in)       :: gsh2

        class(gsh), intent(inout)   :: gsh1

        ! Copy
        gsh1%mass = gsh2%mass
        gsh1%l    = gsh2%l
        gsh1%sfr  = gsh2%sfr
        ! 
        ! Allocate cascade
        if (.not. allocated(gsh1%cascade)) then
            ! Create the cascade with nScales cells
            allocate(gsh1%cascade(nScales))
        end if
        ! Copy each scale
        do s = 1, nScales
            gsh1%cascade(s) = gsh2%cascade(s)
        end do
    
    end subroutine gsh_copy

    ! **********************************
    subroutine gsh_copy_(gsh1, gsh2)

        ! Interface procedure to copy

        implicit none

        type(gsh), intent(in)      :: gsh2
        class(gsh), intent(inout)  :: gsh1

        call gsh_copy(gsh1, gsh2)

    end subroutine gsh_copy_

    ! **********************************
    subroutine gsh_stevolve(this, dt, st, l, inRate, outRates, pStatus, pGsh)

        ! Compute current status of the gas structuration
        ! according to :
        ! - Internal evolution processes
        ! - A constant input rate (due to external process) "inRate" and
        ! - A set of constant output rates (due to external process) "outRates"
        ! Then
        ! Apply the next step "st" of the complete integration scheme

        implicit none
    
        integer(kind=ikd), intent(in)         :: st
        integer(kind=ikd)                     :: il
        integer(kind=ikd)                     :: s

        real(kind=rkd), intent(in)            :: dt          ! The time-step
        real(kind=rkd), intent(in)            :: l           ! The injection scale
        real(kind=rkd)                        :: mass        ! mass contains in the complete cascade

        type(gas), intent(in)                 :: inRate      ! The external input rate
        type(gas), intent(inout), allocatable :: outRates(:) ! The set of output rates
        type(gas)                             :: sclInRate

        type(scale), pointer                  :: pScl        ! Pointer to the current scale
        type(gsh), pointer                    :: pGsh        ! pointer to the current evolved state

        type(status), pointer                 :: pStatus(:)  ! Pointer to the final status
        type(status), pointer                 :: pSclStatus

        class(gsh)                            :: this        ! The current gsh

        ! Get the injection scale bin index
        il = gsh_l2i(l)

        ! Init sclInRate
        call sclInRate%create()

        ! Init total mass
        mass = real(0.d0, kind=rkd)
        !
        ! Reset instantaneous SFR
        call pGsh%sfr%create()

        ! Run through the different scales,
        ! from the largest one to the lowest one
        do s = nScales, 1, -1 
            ! Compute current global status due to
            ! Internal process + external input/output
            ! The largest scale can only be fed by external input
            ! Here sInRate is null for s = nScale
            if (s == il) sclInRate = sclInRate + inRate
            !
            ! Get pointers
            if (st == 1) then
                ! Reset "final" status to null
                call pStatus(s)%reset()
            end if
            pSclStatus => pStatus(s)
            pScl => pGsh%cascade(s)
            !
            if (pScl%gas%mass > 0.d0 .or. sclInRate%mass > 0.d0) then
                ! Apply the new solver step
                ! In output of stevolve, sclInrate contains the transfer rate
                ! from the current scale to the next lower one
                call this%cascade(s)%stevolve(dt, st, sclInRate, outRates(s), pSclStatus, pScl)
                !
                ! Update total cascade mass
                mass = mass + pGsh%cascade(s)%gas%mass
            end if
        end do
        !
        ! Set total mass
        pGsh%mass = mass
        !
        ! Set instantaneous star formation rate
        ! Use the transfer rate of the lowest scale
        if (pStatus(1)%tr%isCreated()) then
            pGsh%sfr = pStatus(1)%tr
        end if
        !
        ! Deallocate local structure
        call sclInRate%delete()

    end subroutine gsh_stevolve

    ! **********************************
    subroutine gsh_evolve(this, dt, l, inRate, outRates)

        ! Evolve, during dt, the current gas structuration
        ! according to :
        ! - Internal evolution processes
        ! - A constant input rate (due to external process) "inRate" and
        ! - A set of constant output rates (due to external process) "outRates"

        implicit none

        integer(kind=ikd)                     :: st          ! Step index of evolution scheme

        real(kind=rkd), intent(in)            :: dt          ! The scale is evolve during dt
        real(kind=rkd), intent(in)            :: l           ! The injection scale

        type(gas), intent(in)                 :: inRate      ! The (dt-)constant input rate
        type(gas), intent(inout), allocatable :: outRates(:) ! The set of output rates

        type(gsh), target                     :: aGsh        ! The current intermediate evolved gas structuration
        type(gsh), pointer                    :: pGsh        ! Pointer to the current evolved gas structuration

        type(status), pointer                 :: pStatus(:)  ! Pointer to the final status

        class(gsh)                            :: this        ! The current gas structuration

        ! Init intermediate status with the current status
        aGsh = this
        pGsh => aGsh
        ! Define complete corrected status
        pStatus => gshStatus
        do st = 1, nSolverStep
            call this%stevolve(dt, st, l, inRate, outRates, pStatus, pGsh)
        end do
        ! Save complete evolved state
        this = aGsh

    end subroutine gsh_evolve

    !
    ! FUNCTIONS
    !

    ! **********************************
    function gsh_i2l(i) result(l)

        ! Return the ith scale of the cascade
        ! in the cell i, gas is structured at scale (l) equal or lower than l(i)
        ! e.g cascade(i=1) stores gas structured at scale lower or equal than lstar = 0.1pc

        implicit none

        integer(kind=ikd), intent(in)  :: i

        real(kind=rkd)                 :: l

        l = stepFactor**(i-1)*lStar
        return

    end function gsh_i2l

    ! **********************************
    function gsh_l2i(l) result(i)

        ! Return the index (in the cascade array) associated to the scale l

        implicit none

        integer(kind=ikd)  :: i, s

        real(kind=rkd)     :: l    ! The scale

        do s = 1, nScales
            if (gsh_i2l(s) >= l) then
                exit
            end if
        end do
        i = min(s, nScales)

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
    function gsh_solve(this, it, dt, aStatusTab, pStatus) result(aGsh)

        ! Apply one solver step

        implicit none

        integer(kind=ikd), intent(in)               :: it                 ! Current solver iteration
        integer(kind=ikd)                           :: s

        real(kind=rkd), intent(in)                  :: dt                 ! Full time step

        type(gsh)                                   :: aGsh               ! The new intermediate scale evolution
        type(status), intent(in), allocatable       :: aStatusTab(:)      ! Set of current global evolution rate
        type(status), pointer                       :: pStatus(:)         ! Pointer to the current corrected scale status

        class(gsh)                                  :: this               ! The current scale

        aGsh = this
        do s = 1, nScales
            aGsh%cascade(s) = this%cascade(s)%solve(it, dt, aStatusTab(s), pStatus(s))
        end do

    end function gsh_solve

end module gsh_mod