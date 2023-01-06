module status_mod

    !*****************************************************************************************************************
    !
    !  OVERVIEW
    !
    ! The status module defines all properties and procedures
    ! associated to intermediate and final transfer rates status
    ! A status is definedfor scale and ssp object according to 
    ! a set of two gas obejcts, dedicated to the current input and ouptut rate
    !
    !*****************************************************************************************************************
    !
    use gas_mod     ! Acces to gas properties and procedures

    implicit none

    public
    !
    ! Define evolution status structure
    type status
        type(gas)   :: in     ! Input rate
        type(gas)   :: tr     ! Internal transfer rate
        type(gas)   :: out    ! Output rate
    contains
        procedure  :: reset => status_reset
        procedure  :: resetFrom => status_resetFrom
        procedure  :: updateFrom => status_updateFrom
        procedure  :: set => status_set
        procedure  :: delete => status_delete
    end type status

contains

    !
    ! SUBROUTINES
    !
    ! **********************************

    subroutine status_reset(this)

        implicit none

        class(status)      :: this      ! The current status

        call this%in%create()    ! Reset input rate
        call this%tr%create()    ! Reset internal transfer rate
        call this%out%create()   ! Reset output rate

    end subroutine status_reset

    subroutine status_resetFrom(this, aStatus, w)

        ! Reset the current status fomr the new one

        implicit none

        real(kind=rkd), intent(in)  :: w         ! weight

        type(status), intent(in)    :: aStatus   ! The new status

        class(status)               :: this      ! The current status

        this%in = w * aStatus%in     ! Reset input rate
        this%tr = w * aStatus%tr     ! Reset internal transfer rate
        this%out = w * aStatus%out   ! Reset output rate

    end subroutine status_resetFrom

    ! **********************************

    subroutine status_updateFrom(this, aStatus, w)

        ! Update the current status from the new one

        implicit none

        real(kind=rkd), intent(in)  :: w        ! weight

        type(status), intent(in)    :: aStatus  ! The new status

        class(status)               :: this     ! The current status

        if (.not. this%in%isCreated()) then
            call status_resetFrom(this, aStatus, w)
        else
            this%in = this%in + w * aStatus%in      ! Reset input rate
            this%tr = this%tr + w * aStatus%tr      ! Reset internal transfer rate
            this%out = this%out + w * aStatus%out   ! Reset output rate
        end if


    end subroutine status_updateFrom

    ! **********************************

    subroutine status_set(this, in, tr, out)

        ! Update the current status from the new one

        implicit none

        type(gas), intent(in)  :: in       ! Input rate
        type(gas), intent(in)  :: tr       ! Internal transfer rate
        type(gas), intent(in)  :: out      ! Output rate

        class(status)          :: this     ! The current status

        this%in = in      ! Set input rate
        this%tr = tr      ! Set internal tranfer rate
        this%out = out    ! Set output rate

    end subroutine status_set

    ! **********************************
    subroutine status_delete(this)

        ! Delete and deallocate all associated object

        class(status)       :: this     ! The current status

        call this%in%delete()
        call this%tr%delete()
        call this%out%delete()

    end subroutine status_delete

end module status_mod