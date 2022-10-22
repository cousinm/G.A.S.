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

        integer(kind=4), parameter :: u = 333       ! file unit

        logical                    :: isValid

        character(MAXPATHSIZE)     :: filename

        ! Open log file for these tests
        write(filename,'(a, a)') trim(validPath), '/scale_tests.log'
        open(unit=u, file=filename, status='new')

        isValid = test_scale_constant_injection()
        write(u, '(a,l)') 'test_scale_constant_injection: ', isValid

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

        integer(kind=4), parameter   :: u = 334          ! file unit

        logical                      :: isValid

        character(MAXPATHSIZE)       :: filename

        real(kind=8), parameter      :: l = 1.d2*pc2kpc  ! [pc] Injection scale
        real(kind=8), parameter      :: dt = 1.d-4       ! CU [Gyr]
        real(kind=8), parameter      :: evolTime = 1.d0  ! CU [Gyr]
        real(kind=8)                 :: t

        type(gas)                    :: inRate  ! The constant injection rate
        type(scale)                  :: scl     ! A gas scale

        isValid = .FALSE.

        inRate = 1.d1 * MassRate_CU * initAbund(nMetBins)  ! 10Msun/yr in CU

        ! Create a scale
        call scl%create(l)

        ! Open data file for this test
        write(filename,'(a, a)') trim(validPath), '/scale_test_constant_injection.dat'
        open(unit=u, file=filename, status='new')

        ! Evolution
        t = 0.d0 ! init
        write(u, '(a)') 'nClouds | mass | mZ | Eint | mH | mHe | mC | mN | mO | mFe'
        do while (t < evolTime)
            write(u, *) scl%nClouds(), scl%gas%mass, scl%gas%mZ, scl%gas%Eint, &
                        scl%gas%elts(1), scl%gas%elts(2), scl%gas%elts(3), &
                        scl%gas%elts(4), scl%gas%elts(5), scl%gas%elts(6)
            call scl%solve(dt, inRate)
            t = t + dt
        end do

        ! Close data file
        close(u)

        ! Delete structure
        call scl%delete()

        isValid = .TRUE.

    end function test_scale_constant_injection

end module scale_tests_mod
