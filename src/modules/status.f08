module status_mod

    !*****************************************************************************************************************
    !
    !  OVERVIEW
    !
    ! The status module defines all properties and procedures
    ! associated to intermediate and final transfer rates status
    ! A status is defined for scale and ssp object according to 
    ! a set of gas objects, dedicated to the current input and ouptut rate
    !
    !*****************************************************************************************************************
    !
    use gas_mod     ! Acces to gas properties and procedures
    use model_mod   ! Acces to model parameters and configurations
    use solver_mod  ! Acces to solver parameters and methods

    implicit none

    public
    !
    ! Define evolution status structure
    type status
        real(kind=rkd)  :: mass   ! Reference mass (set at st = 1)
        type(gas)       :: in     ! Input rate
        type(gas)       :: tr     ! Internal transfer rate
        type(gas)       :: out    ! Output rate
        
    contains
        procedure  :: create => status_create
        procedure  :: update => status_update
        procedure  :: set => status_set
        procedure  :: delete => status_delete
        procedure  :: dtMax => status_dtmax
    end type status

contains

    !
    ! SUBROUTINES
    !
    ! **********************************

    subroutine status_create(this)

        implicit none

        class(status)      :: this           ! The current status

        call this%in%create()                ! Create input rate
        call this%tr%create()                ! Create internal transfer rate
        call this%out%create()               ! Create output rate

    end subroutine status_create

    ! **********************************
    subroutine status_set(this, in, tr, out)

        ! Update the current status from the new one

        implicit none

        type(gas), intent(in)   :: in       ! Input rate
        type(gas), intent(in)   :: tr       ! Internal transfer rate
        type(gas), intent(in)   :: out      ! Output rate
        

        class(status)           :: this     ! The current status

        this%in = in        ! Set input rate
        this%tr = tr        ! Set internal tranfer rate
        this%out = out      ! Set output rate

    end subroutine status_set

    ! **********************************
    subroutine status_delete(this)

        ! Delete and deallocate all associated object

        class(status)       :: this     ! The current status

        call this%in%delete()
        call this%tr%delete()
        call this%out%delete()

    end subroutine status_delete

        ! **********************************
    subroutine status_update(this, st, aStatus)

        ! Update the current status from the new one

        implicit none

        integer(kind=ikd), intent(in) :: st        ! Current solver step

        real(kind=rkd)                :: w         ! Weight

        type(status), target          :: aStatus   ! The new status

        class(status), target         :: this      ! The current status

        ! Get status ponderation
        w = solver_w(st)
        !
        ! With 'RK2' scheme, pStatus is alaways reset
        ! With 'RK4' scheme, pStatus is only reset at the first step
        if (st == 1 .or. solver == 'RK2') then
            call this%create()   ! Reset status
        end if
        !
        ! Update the "final" status
        this%in = this%in + w * aStatus%in    ! Input rate
        this%tr = this%tr + w * aStatus%tr    ! Internal transfer rate
        this%out = this%out + w * aStatus%out ! Output rate

    end subroutine status_update

    !
    ! FUNCTIONS
    !

    ! **********************************
    function status_dtmax(this, dt, mass) result(dt_max)

        ! Return the maximal time step value possible to keep mass variations positive

        real(kind=rkd), intent(in)    :: dt
        real(kind=rkd), intent(in)    :: mass
        real(kind=rkd)                :: dmOut, dmIn, dmTr ! Output, input, and transfered masses
        real(kind=rkd)                :: dm
        real(kind=rkd)                :: dt_max

        class(status)                 :: this     ! The current status
        !
        ! Init to current dt
        dt_max = dt
        !
        ! Adaptive time-step
        ! Compute mass variation
        dmIn = dt * this%in%mass
        dmTr = dt * this%tr%mass
        dmOut = dt * this%out%mass
        !
        ! Compute complete variation
        dm = dmIn - dmTr - dmOut
        !
        if (dm < 0.) then
            ! Mass variation is strickly negative
            ! dt should be adapted
            dt_max = min(dt_max, mass / abs(this%in%mass - this%tr%mass - this%out%mass))
        end if
    end function status_dtmax

end module status_mod