module scale_tests_mod

    !*****************************************************************************************************************
    !
    !  OVERVIEW
    !
    !  This module is dedicated to scale module tests.
    !
    !*****************************************************************************************************************

    use parameters ! Acces to global defintions and properties
    use config_mod ! Acces to configurations parameters (path)
    use scale_mod  ! Acces to scale properties and procedures
    use log_mod    ! Acces to logging procedures

    implicit none

    public

contains

    !
    ! SUBROUTINES
    !

    ! **********************************
    subroutine scale_all_tests()

        ! Run all tests for scale module
    
        implicit none

        integer(kind=ikd), parameter :: u = 333       ! file unit

        logical                    :: isValid

        character(MAXPATHSIZE)     :: filename

        real(kind=rkd)             :: tstart, tend

        ! Open log file for these tests
        write(filename,'(a, a)') trim(validPath), '/scale_tests.log'
        open(unit=u, file=filename, status='new')

        call cpu_time(tstart)
        isValid = test_scale_constant_injection()
        call cpu_time(tend)
        write(u, '(a,f7.3,a,l)') 'test_scale_constant_injection (', tend-tstart, ' sec): ', isValid

        ! Close log file
        close(u)

        return
    end subroutine scale_all_tests

    !
    ! FUNCTIONS
    !

    ! **********************************
    function test_scale_constant_injection() result(isValid)

        implicit none

        integer(kind=ikd), parameter   :: u = 334          ! file unit

        logical                      :: isValid

        character(MAXPATHSIZE)       :: filename

        real(kind=rkd), parameter    :: l = 1.d2*pc2kpc  ! [pc] Injection scale
        real(kind=rkd), parameter    :: dt = 1.d-4       ! CU [Gyr]
        real(kind=rkd), parameter    :: evolTime = 1.d0  ! CU [Gyr]
        real(kind=rkd)               :: t
        real(kind=rkd)               :: solution, diff

        type(gas)                    :: inRate  ! The constant injection rate
        type(gas)                    :: outRate ! The constant ejection rate
        type(gas)                    :: ejGas   ! Gas ejected from the scale

        type(scale)                  :: scl     ! A gas scale

        isValid = .TRUE.

        inRate = 1.d1 * MassRate_CU * initAbund(nMetBins)  ! 10Msun/yr in CU
        call outRate%create()  ! null

        ! Create a scale
        call scl%create(l)
        ! Create gas reservoir
        call ejGas%create()

        ! Open data file for this test
        write(filename,'(a, a)') trim(validPath), '/scale_test_constant_injection.dat'
        open(unit=u, file=filename, status='new')

        ! Evolution
        t = 0.d0         ! init
        solution = 0.d0  ! init
        diff = 0.d0      ! init
        write(u, '(a)') 'nClouds | scl mass [CU] | gas mass [CU] | solution | diff'
        do while (t < evolTime)
            !
            write(u, *) scl%nClouds(), scl%gas%mass, ejGas%mass, solution, diff
            !
            ! Compute evolution
            call scl%evolve(dt, inRate, outRate)
            !
            ! Compute real solution
            solution = solution + inRate%mass * dt
            !
            ! Update ejected gas reservoir
            ejGas = ejGas + dt * outRate
            !
            ! Test, mass conservation
            diff = abs(ejGas%mass + scl%gas%mass - solution)
            if (diff > num_accuracy) then
                isValid = .FALSE.
            end if
            t = t + dt
        end do

        ! Close data file
        close(u)

        ! Delete structure
        call scl%delete()
        call ejGas%delete()

    end function test_scale_constant_injection

end module scale_tests_mod
