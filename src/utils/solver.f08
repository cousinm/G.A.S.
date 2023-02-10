module solver_mod

    !*****************************************************************************************************************
    !
    !  OVERVIEW
    !
    !  This module contains all method associated to integration scheme
    !
    !*****************************************************************************************************************

    use parameters  ! Acces to global defintions and properties

    implicit none

    public

    ! Define integration scheme
    character(MAXPATHSIZE)   :: solver

    integer(kind=ikd)        :: nSolverStep

    real(kind=rkd)           :: solver_dt

contains

    !
    ! SUBROUTINES
    !

    subroutine solver_init(solver_name)

        ! Initialize the solver

        character(*), intent(in)   :: solver_name
        !
        ! Set
        solver = solver_name
        !
        ! Set parameter depending to the solver
        select case (trim(solver))
        case ('RK4')
            ! Range-Kutta 4th order
            nSolverStep = 4
        case default
            ! Range Kutte 2d order
            nSolverStep = 2
        end select

    end subroutine solver_init

    !
    ! FUNCTIONS
    !

    function solver_isFinalStep(st) result(isFinalStep)

        ! Return true if the current step is the final one

        implicit none

        integer(kind=ikd), intent(in)   :: st

        logical                         :: isFinalStep

        if (st == nSolverStep) then
            isFinalStep = .TRUE.
        else
            isFinalStep = .FALSE.
        end if
    end function solver_isFinalStep

    ! **********************************
    function solver_w(st) result(w)

        ! Return the "status" ponderation factor associated to the current integrator step

        implicit none

        integer(kind=ikd), intent(in)  :: st

        real(kind=rkd)                 :: w

        select case (trim(solver))
        case ('RK4')
            ! Range-Kutta 4th order
            select case (st)
            case (2, 3)
                ! Intermediate step
                w = real(2.d0, kind=rkd)/real(6.d0, kind=rkd)
            case default
                ! First or final step
                w = real(1.d0, kind=rkd)/real(6.d0, kind=rkd)
            end select
        case default
            ! Range-Kutta 2d order
            w = real(1.d0, kind=rkd)
        end select
    end function solver_w

    ! **********************************
    function solver_wdt(st) result(wdt)

        ! Return the "time-step" ponderation factor associated to the current integrator step

        implicit none

        integer(kind=ikd), intent(in)  :: st

        real(kind=rkd)                 :: wdt

        select case (trim(solver))
        case ('RK4')
            ! Range-Kutta 4th order
            select case (st)
            case (1, 2, 3)
                ! First and intermediate step
                wdt = real(1.d0, kind=rkd)/real(2.d0, kind=rkd)
            case default
                ! Final step
                wdt = real(1.d0, kind=rkd)
            end select
        case default
            ! Range-Kutta 2d order
            select case (st)
            case (1)
                ! First step
                wdt = real(1.d0, kind=rkd)/real(2.d0, kind=rkd)
            case default
                ! Final step
                wdt = real(1.d0, kind=rkd)
            end select
        end select
    end function solver_wdt

end module solver_mod
