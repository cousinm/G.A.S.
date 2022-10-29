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
        procedure  :: solve => gsh_solve          ! Solve, for dt, the next step of the integration scheme
        procedure  :: stevolve => gsh_stevolve    ! Evolve, for dt, the next step of the integration scheme
        procedure  :: evolve => gsh_evolve        !         according to internal/external input/output
        procedure  :: Vesc => gsh_escape_velocity ! Return the escape velocity of the gs
    end type gsh

    ! TEMPORARY INTERMEDIATE STATUS
    ! Targets
    type(gas), allocatable, target     :: gshStatus(:)

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

        integer(kind=ikd)    :: s

        call scale_init()

        ! Init intermediate status
        allocate(gshStatus(nScales))
        do s = 1, nScales
            call gshStatus(s)%create()
        end do

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
    subroutine gsh_stevolve(this, dt, st, l, inRate, outRates, pStatus, pGs)

        ! Compute current status of the gas structuration
        ! according to :
        ! - Internal evolution processes
        ! - A constant input rate (due to external process) "inRate" and
        ! - A set of constant output rates (due to external process) "outRates"
        ! Then
        ! Apply the next step "st" of the complete integration scheme

        implicit none
    
        integer(kind=ikd), intent(in)       :: st
        integer(kind=ikd)                   :: il
        integer(kind=ikd)                   :: s

        real(kind=rkd), intent(in)          :: dt          ! The time-step
        real(kind=rkd), intent(in)          :: l           ! The injection scale
        real(kind=rkd)                      :: mass        ! mass contains in the complete cascade

        type(gas), intent(in)               :: inRate      ! The external input rate
        type(gas), intent(in), allocatable  :: outRates(:) ! The set of output rates
        type(gas), pointer                  :: pStatus(:)  ! Pointer to the complete corrected status
        type(gas), pointer                  :: pSclStatus
        type(gas)                           :: sInRate

        type(scale), pointer                :: pScl        ! Pointer to the current scale
        type(gsh), pointer                  :: pGs         ! pointer to the current evolved state

        class(gsh)                          :: this        ! The current gsh

        ! Get the injection scale bin index
        il = gsh_l2i(l)

        ! Init sInRate
        call sInRate%create()

        ! Init total mass
        mass = real(0.d0, kind=rkd)
        do s = nScales, 1, -1
            ! Compute current global status due to
            ! Internal process + external input/output
            ! The largest scale can only be fed by external input
            ! Here sInRate is null for s = nScale
            if (s == il) sInRate = sInRate + inRate
            !
            ! Apply the new solver step
            pSclStatus => pStatus(s)
            pScl => pGs%cascade(s)
            call this%cascade(s)%stevolve(dt, st, sInRate, outRates(s), pSclStatus, pScl)
            !
            ! Save the inRate of the next lower scale
            sInRate = this%cascade(s)%status()
            !
            ! Update total cascade mass
            mass = mass + pGs%cascade(s)%gas%mass
        end do
        ! Set total mass
        pGs%mass = mass

    end subroutine gsh_stevolve

    ! **********************************
    subroutine gsh_evolve(this, dt, l, inRate, outRates)

        ! Evolve, during dt, the current gas structuration
        ! according to :
        ! - Internal evolution processes
        ! - A constant input rate (due to external process) "inRate" and
        ! - A set of constant output rates (due to external process) "outRates"

        implicit none

        integer(kind=ikd)                   :: st          ! Step index of evolution scheme

        real(kind=rkd), intent(in)          :: dt          ! The scale is evolve during dt
        real(kind=rkd), intent(in)          :: l           ! The injection scale

        type(gas), intent(in)               :: inRate      ! The (dt-)constant input rate
        type(gas), intent(in), allocatable  :: outRates(:) ! The set of output rates
        type(gas), pointer                  :: pStatus(:)  ! Pointer to the corrected status

        type(gsh), target                   :: gs          ! The current intermediate evolved scale
        type(gsh), pointer                  :: pGsh        ! Pointer to the curretn evolved scale

        class(gsh)                          :: this        ! The current scale

        ! Init intermediate status with the current status
        gs = this
        pGsh => gs
        ! Define complete corrected status
        pStatus => gshStatus
        do st = 1, nSolverStep
            call this%stevolve(dt, st, l, inRate, outRates, pStatus, pGsh)
        end do
        ! Save complete evolved state
        this = gs

    end subroutine gsh_evolve

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

        ! Return the current status of the gas structuration
        ! according to its internal evolution only

        implicit none
  
        integer(kind=ikd)       :: s

        type(gas), allocatable  :: rates(:)   ! Set of output rates (can be negative)
        type(gas)               :: sInRate    ! in at s = out at s+1

        class(gsh)              :: this
        !
        ! Create the set of global evolution rates
        allocate(rates(nScales))
        ! Init
        call sInRate%create()
        ! Set values, run from the largest to the lowest scale
        do s = nScales, 1, -1
            ! The largest scale can only be fed by external input
            ! Here sInRate is null for s = nScale
            rates(s) = sInRate - this%cascade(s)%status()
            ! Bkp current scale status (current output) for the next scale
            sInRate = this%cascade(s)%status()
        end do

    end function gsh_status

    ! **********************************
    function gsh_solve(this, it, dt, status, pStatus) result(gs)

        ! Apply one solver step

        implicit none

        integer(kind=ikd), intent(in)       :: it         ! Current solver iteration
        integer(kind=ikd)                   :: s

        real(kind=rkd), intent(in)          :: dt         ! Full time step

        type(gas), intent(in), allocatable  :: status(:)  ! Set of current global evolution rate
        type(gas), pointer                  :: pStatus(:) ! Pointer to the current corrected scale status

        type(gsh)                           :: gs         ! The new intermediate scale evolution

        class(gsh)                          :: this       ! The current scale

        gs = this
        do s = 1, nScales
            gs%cascade(s) = this%cascade(s)%solve(it, dt, status(s), pStatus(s))
        end do

    end function gsh_solve

end module gsh_mod