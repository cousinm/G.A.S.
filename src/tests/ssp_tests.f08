module ssp_tests_mod

    !*****************************************************************************************************************
    !
    !  OVERVIEW
    !
    !  This module is dedicated to ssp module tests.
    !
    !*****************************************************************************************************************

    use parameters ! Acces to global defintions and properties
    use config_mod ! Acces to configurations parameters (path)
    use ssp_mod    ! Acces to ssp properties and procedures
    use log_mod    ! Acces to logging procedures

    implicit none

    public

contains

    !
    ! SUBROUTINES
    !

    ! **********************************
    subroutine ssp_all_tests()

        ! Run all tests for ssp module
    
        implicit none

        integer(kind=ikd), parameter :: u = 333       ! file unit

        logical                      :: isValid

        character(MAXPATHSIZE)       :: filename

        real(kind=rkd)               :: tstart, tend

        ! Open log file for these tests
        write(filename,'(a, a)') trim(validPath), '/ssp_tests.log'
        open(unit=u, file=filename, status='new')

        call cpu_time(tstart)
        isValid = test_ssp_constant_injection()
        call cpu_time(tend)
        write(u, '(a,f7.3,a,l)') 'test_ssp_constant_injection (', tend-tstart, ' sec): ', isValid

        ! Close log file
        close(u)

        return
    end subroutine ssp_all_tests

    !
    ! FUNCTIONS
    !

    ! **********************************
    function test_ssp_constant_injection() result(isValid)

        implicit none

        integer(kind=ikd), parameter :: u = 334          ! file unit
        integer(kind=ikd), parameter :: iAge = 10

        logical                      :: isValid

        character(MAXPATHSIZE)       :: filename

        real(kind=rkd), parameter    :: l = 1.d2*pc2kpc  ! [pc] Injection ssp
        real(kind=rkd), parameter    :: dt = 1.d-4       ! CU [Gyr]
        real(kind=rkd)               :: adt
        real(kind=rkd), parameter    :: evolTime = 1.d0  ! CU [Gyr]
        real(kind=rkd)               :: t
        real(kind=rkd)               :: solution, diff

        type(gas)                    :: inRate  ! The constant injection rate
        type(gas)                    :: ejGas   ! Gas ejected from the ssp

        type(ssp)                    :: aSsp    ! A ssp

        isValid = .TRUE.

        inRate = 1.d1 * MassRate_CU * initAbund(nMetBins)  ! 10Msun/yr in CU

        ! Create a ssp
        call aSsp%create(iAge, nMetBins)
        ! Create gas reservoir
        call ejGas%create()

        ! Open data file for this test
        write(filename,'(a, a)') trim(validPath), '/ssp_test_constant_injection.dat'
        open(unit=u, file=filename, status='new')

        ! Evolution
        t = 0.d0         ! init
        solution = 0.d0  ! init
        diff = 0.d0      ! init
        write(u, '(a)') 't | ssp age [CU] | ssp mass [CU] | gas mass [CU] | solution | diff'
        do while (t < evolTime)
            !
            write(u, *) t, aSsp%avgAge, aSsp%mass, ejGas%mass, solution, diff
            !
            ! Compute evolution
            adt = dt
            call aSsp%evolve(adt, inRate)
            !
            ! Compute real solution
            solution = solution + inRate%mass * adt
            !
            ! Update ejected gas reservoir
            ejGas = ejGas + adt * (mySspStatus%out + mySspStatus%tr)
            !
            ! Test, mass conservation
            diff = abs(ejGas%mass + aSsp%mass - solution)
            if (diff > num_accuracy) then
                isValid = .FALSE.
            end if
            t = t + adt
        end do

        ! Close data file
        close(u)

        ! Delete structure
        call ejGas%delete()

    end function test_ssp_constant_injection

end module ssp_tests_mod
