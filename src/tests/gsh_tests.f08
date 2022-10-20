module gsh_tests_mod

    !*****************************************************************************************************************
    !
    !  OVERVIEW
    !
    !  This module is dedicated to gas structuration history module tests.
    !
    !*****************************************************************************************************************

    use gsh_mod    ! Acces to gas structuration history module
    use config_mod ! Acces to configurations parameters (path)
    use log_mod    ! Acces to logging procedures

    implicit none

    public

contains

    !
    ! SUBROUTINES
    !

    ! **********************************
    subroutine gsh_all_tests()

        ! Run all tests for gsh module
    
        implicit none

        integer(kind=4), parameter :: u = 333       ! file unit

        logical                    :: isValid

        character(MAXPATHSIZE)     :: filename

        ! Open log file for these tests
        write(filename,'(a, a)') trim(validPath), '/gsh_tests.log'
        open(unit=u, file=filename, status='new')

        isValid = test_gsh_create()
        write(u, '(a,l)') 'test_gsh_create: ', isValid

        isValid = test_gsh_delete()
        write(u, '(a,l)') 'test_gsh_delete: ', isValid

        isValid = test_gsh_constant_injection()
        write(u, '(a,l)') 'test_gsh_constant_injection: ', isValid

        ! Close log file
        close(u)

        return
    end subroutine gsh_all_tests

    !
    ! FUNCTIONS
    !

    ! **********************************
    function test_gsh_create() result(isValid)

        logical           :: isValid

        type(gsh)         :: g

        ! Create a gsh object
        call g%create()

        ! Tests
        isValid = .TRUE.
        if ((g%mass < 0.) .or. (g%mass > 0.)) then
            isValid = .FALSE.
        else if (.not. allocated(g%cascade)) then
            isValid = .FALSE.
        end if

    end function test_gsh_create

    ! **********************************
    function test_gsh_delete() result(isValid)

        logical           :: isValid

        type(gsh)         :: g

        ! Create a gsh object
        call g%create()
        ! Delete it
        call g%delete()

        ! Tests
        isValid = .TRUE.
        if ((g%mass < 0.) .or. (g%mass > 0.)) then
            isValid = .FALSE.
        else if (allocated(g%cascade)) then
            isValid = .FALSE.
        end if

    end function test_gsh_delete

        ! **********************************
    function test_gsh_constant_injection() result(isValid)

        implicit none

        integer(kind=4)              :: is
        integer(kind=4), parameter   :: u = 334          ! file unit

        logical                      :: isValid

        character(MAXPATHSIZE)       :: filename

        real(kind=8), parameter      :: l = 1.d0*pc2kpc  ! [pc] Injection scale
        real(kind=8), parameter      :: dt = 1.d-4       ! CU [Gyr]
        real(kind=8), parameter      :: evolTime = 1.d0  ! CU [Gyr]
        real(kind=8)                 :: t

        type(gas)                    :: inRate  ! The constant injection rate
        type(gsh)                    :: gs      ! gas structuration history

        isValid = .FALSE.

        inRate = 1.d1 * MassRate_CU * initAbund(nMetBins)  ! 10Msun/yr in CU

        ! Create a scale
        call gs%create()

        ! Open data file for this test
        write(filename,'(a, a)') trim(validPath), '/gsh_test_constant_injection.dat'
        open(unit=u, file=filename, status='new')

        ! Evolution
        is = gsh_l2i(l)
        t = 0.d0 ! init
        write(u, '(a)') 'mass | nClouds(nScale) | nClouds(1) | mass(nScales) | mass(1)'
        do while (t < evolTime)
            write(u, *) gs%mass, gs%cascade(is)%nClouds(), gs%cascade(1)%nClouds(), &
                        gs%cascade(is)%gas%mass, gs%cascade(1)%gas%mass
            call gs%solve(dt, inRate, l)
            t = t + dt
        end do

        ! Close data file
        close(u)

        ! Delete structure
        call gs%delete()

        isValid = .TRUE.

    end function test_gsh_constant_injection

end module gsh_tests_mod
