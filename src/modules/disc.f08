module disc_mod

    use gsh_mod
    use sp_mod

    !*****************************************************************************************************************
    !
    !  OVERVIEW
    !
    !  The disc module defines all properties and procedures assossiated to a disc object.
    !  In G.A.S. a disc contains mainly a gas structuration history structure (gsh)
    !  strongly linked to the stellar population (sp) evolution
    !  Incomming gas feeds the gsh, in cascade output new stars are formed
    !  and evolved in the sp structure
    !  This module starts with the definition of this disc type
    !
    !*****************************************************************************************************************

    implicit none

    public

    type disc
        type(gsh)  :: myGsh       ! The gas structuration history
        type(sp)   :: mySp        ! The stellar population
    contains
        procedure  :: create => disc_create     ! Create a disc
        procedure  :: delete => disc_delete     ! Delete a disc
        procedure  :: stevolve => disc_stevolve ! Apply the next integration scheme step
        procedure  :: evolve => disc_evolve     ! Evolve a disc structure during dt
    end type disc

    ! INTERFACE OPERATOR DECLARATIONS

    interface assignment (=)  ! allows to copy a disc component by using the symbol '='
        module procedure disc_copy_
    end interface assignment (=)

contains

    !
    ! SUBROUTINES
    !

    ! **********************************
    subroutine disc_init()

        ! Initialize the disc module and dependancies

        implicit none

        ! Init gas structuration history
        call gsh_init()
        ! Init stellar population
        call sp_init()

    end subroutine disc_init

    ! **********************************
    subroutine disc_create(this)

        ! Create a disc structure

        implicit none

        class(disc)     :: this

        ! Create the gas structuration history
        call this%myGsh%create()
        !
        ! Create the stellar popualtion
        call this%mySp%create()

    end subroutine disc_create

    ! **********************************
    subroutine disc_delete(this)

        ! Delete a disc structure

        implicit none

        class(disc)     :: this

        ! Create the gas structuration history
        call this%myGsh%delete()
        !
        ! Create the stellar popualtion
        call this%mySp%delete()

    end subroutine disc_delete

    ! **********************************
    subroutine disc_copy(disc1, disc2)

        ! Copy the disc object disc2 into the disc object disc1

        implicit none

        type(disc), intent(in)       :: disc2

        class(disc), intent(inout)   :: disc1

        ! Copy
        ! gsh
        disc1%myGsh = disc2%myGsh
        ! sp
        disc1%mySp = disc2%mySp
    
    end subroutine disc_copy

    ! **********************************
    subroutine disc_copy_(disc1, disc2)

        ! Interface procedure to copy

        implicit none

        type(disc), intent(in)       :: disc2

        class(disc), intent(inout)   :: disc1

        call disc_copy(disc1, disc2)

    end subroutine disc_copy_

    ! **********************************
    subroutine disc_stevolve(this, dt, st, l, inRate, outRate, pDisc)

        ! Compute current status of the disc
        ! according to :
        ! - Internal evolution processes
        ! - A constant input rate (fresh gas accretion) "inRate"
        ! Then
        ! Apply the next step "st" of the complete integration scheme
        !
        ! In output, outRate in the gas ejection rate produced

        implicit none
    
        integer(kind=ikd), intent(in)   :: st

        real(kind=rkd), intent(in)      :: dt               ! The time-step
        real(kind=rkd), intent(in)      :: l                ! The gas injection scale
        real(kind=rkd)                  :: Qsn              ! SN gas disruption power

        type(gas), intent(in)           :: inRate           ! The external input rate
        type(gas), intent(out)          :: outRate          ! The output rate
        type(gas)                       :: spOutRate        ! The stellar population outRate (ejecta)
        type(gas)                       :: iSFR             ! The instantaneous SFR
        type(gas)                       :: gshInRate        ! Total gsh input rate, fresh input + ejected gas from sp
        type(gas)                       :: ejR

        type(status), pointer           :: pGshStatus(:)    ! Pointer to the gsh status
        type(status), pointer           :: pSpStatus(:, :)  ! Pointer to the sp status

        type(gsh), pointer              :: pGsh             ! pointer to the current gsh stage
        type(sp), pointer               :: pSp              ! pointer to the current stage of the sp

        type(disc), pointer             :: pDisc            ! Pointer to the current stage

        class(disc)                     :: this             ! The current gsh

        ! Get current status for the gas structuration history
        ! and the stellar population
        pGshStatus => gshStatus
        pSpStatus => spStatus
        pGsh => pDisc%myGsh
        pSp => pDisc%mySp
        !
        ! Init gas objects
        call gshInRate%create()
        !
        ! 1- Get instantaneous SN disruption power, the energy per time unit generes by SNs
        Qsn = pSp%iSNP()
        !
        ! 2- Get instantaneous gas ejection rate produced by the stellar population
        !    Add to the current fresh gas input rate
        ejR = pSp%iEJR()
        gshInRate = inRate + ejR
        !
        ! 3- Get the instantaneous SFR
        iSFR = pGsh%iSFR()
        !
        ! 4- Evolve gas structuration history
        call this%myGsh%stevolve(dt, st, l, Qsn, gshInRate, outRate, pGshStatus, pGsh)
        !
        ! 5- Evolve the stellar population
        call this%mySp%stevolve(dt, st, iSFR, spOutRate, pSpStatus, pSp)

    end subroutine disc_stevolve

    ! **********************************
    subroutine disc_evolve(this, dt, l, inRate, outRate)

        ! Evolve, during dt, the current disc
        ! according to :
        ! - Internal evolution processes
        ! - A constant input rate (due to external process) "inRate" and
        ! Gas is injected at scale l (deduced from dark-matter properties)
        !
        ! In output, outRate save the gas ejection rate

        implicit none

        integer(kind=ikd)              :: st          ! Step index of evolution scheme

        real(kind=rkd), intent(in)     :: dt          ! The disc evolves during dt
        real(kind=rkd), intent(in)     :: l           ! The disc gas injection scale

        type(gas), intent(in)          :: inRate      ! The (dt-)constant input rate
        type(gas), intent(out)         :: outRate     ! The gas ejection rate

        type(disc), target             :: aDisc       ! The current stage
        type(disc), pointer            :: pDisc       ! Pointer to the current stage

        class(disc)                    :: this        ! The current disc

        ! Init intermediate status with the current status
        aDisc = this
        pDisc => aDisc
        ! Apply solver steps
        do st = 1, nSolverStep
            call this%stevolve(dt, st, l, inRate, outRate, pDisc)
        end do
        ! Save complete evolved state
        this = aDisc
        ! Delete the tmp disc structure
        call aDisc%delete()

    end subroutine disc_evolve

end module disc_mod