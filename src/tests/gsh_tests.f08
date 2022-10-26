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

        integer(kind=ikd), parameter :: u = 333       ! file unit

        logical                    :: isValid

        real(kind=rkd)             :: tstart, tend

        character(MAXPATHSIZE)     :: filename

        ! Open log file for these tests
        write(filename,'(a, a)') trim(validPath), '/gsh_tests.log'
        open(unit=u, file=filename, status='new')

        isValid = test_gsh_create()
        write(u, '(a,l)') 'test_gsh_create: ', isValid

        isValid = test_gsh_delete()
        write(u, '(a,l)') 'test_gsh_delete: ', isValid

        call cpu_time(tstart)
        isValid = test_gsh_constant_injection()
        call cpu_time(tend)
        write(u, '(a,f7.3,a,l)') 'test_gsh_constant_injection (', tend-tstart, ' sec): ', isValid

        call cpu_time(tstart)
        isValid = test_gsh_constant_and_stop_injection()
        call cpu_time(tend)
        write(u, '(a,f7.3,a,l)') 'test_gsh_constant_and_stop_injection (', tend-tstart, ' sec): ', isValid

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

        integer(kind=ikd)              :: is
        integer(kind=ikd)              :: i
        integer(kind=ikd), parameter   :: u = 10000        ! file unit

        logical                      :: isValid

        character(MAXPATHSIZE)       :: filename

        real(kind=rkd), parameter      :: l = 1.d2*pc2kpc  ! [pc] Injection scale
        real(kind=rkd), parameter      :: dt = 1.d-4       ! CU [Gyr]
        real(kind=rkd), parameter      :: evolTime = 1.d0  ! CU [Gyr]
        real(kind=rkd)                 :: t

        type(gas)                    :: inRate  ! The constant injection rate
        type(gsh)                    :: gs      ! gas structuration history

        isValid = .FALSE.

        inRate = 1.d1 * MassRate_CU * initAbund(nMetBins)  ! 10Msun/yr in CU

        ! Create a scale
        call gs%create()

        ! Open data files for this test
        ! Save each scale in a dedicated file
        do i = 1, nScales
            write(filename, '(a,a,i2.2,a)') trim(validPath), '/gsh_test_constant_injection_s', i, '.dat'
            open(unit=u+i, file=filename, status='new')
            write(u+i, '(a)') '# t | mass | ll | Vesc [km/s] | nClouds(s) | mass(s)'
        end do

        ! Evolution
        is = gsh_l2i(l)
        t = 0.d0 ! init
        do while (t < evolTime)
            ! Save each scale in a dedicated file
            do i = 1, nScales
                write(u+i, *) t, gs%mass, gs%l, gs%Vesc()*Velocity_km_s, gs%cascade(i)%nClouds(), gs%cascade(i)%gas%mass
            end do
            call gs%solve(dt, inRate, l)
            t = t + dt
        end do

        ! Close data files
        do i = 1, nScales
            close(u+i)
        end do

        ! Delete structure
        call gs%delete()

        isValid = .TRUE.

    end function test_gsh_constant_injection

    ! **********************************
    function test_gsh_constant_and_stop_injection() result(isValid)

        implicit none

        integer(kind=ikd)              :: is
        integer(kind=ikd)              :: i
        integer(kind=ikd), parameter   :: u = 10000        ! file unit

        logical                      :: isValid

        character(MAXPATHSIZE)       :: filename

        real(kind=rkd), parameter      :: l = 1.d2*pc2kpc  ! [pc] Injection scale
        real(kind=rkd), parameter      :: dt = 1.d-4       ! CU [Gyr]
        real(kind=rkd), parameter      :: evolTime = 2.d0  ! CU [Gyr]
        real(kind=rkd)                 :: t

        type(gas)                    :: inRate  ! The constant injection rate
        type(gsh)                    :: gs      ! gas structuration history

        isValid = .FALSE.

        inRate = 1.d1 * MassRate_CU * initAbund(nMetBins)  ! 10 Msun/yr in CU

        ! Create a scale
        call gs%create()

        ! Open data files for this test
        ! Save each scale in a dedicated file
        do i = 1, nScales
            write(filename, '(a,a,i2.2,a)') trim(validPath), '/gsh_test_constant_and_stop_injection_s', i, '.dat'
            open(unit=u+i, file=filename, status='new')
            write(u+i, '(a)') '# t | mass | ll | Vesc [km/s] | nClouds(s) | mass(s)'
        end do

        ! Evolution
        is = gsh_l2i(l)
        t = 0.d0 ! init
        do while (t < evolTime)
            ! Save each scale in a dedicated file
            do i = 1, nScales
                write(u+i, *) t, gs%mass, gs%l, gs%Vesc()*Velocity_km_s, gs%cascade(i)%nClouds(), gs%cascade(i)%gas%mass
            end do
            if (t > 1.0) then
                inRate = real(0.d0, kind=rkd) * initAbund(nMetBins)
            end if
            call gs%solve(dt, inRate, l)
            t = t + dt
        end do

        ! Close data files
        do i = 1, nScales
            close(u+i)
        end do

        ! Delete structure
        call gs%delete()

        isValid = .TRUE.

    end function test_gsh_constant_and_stop_injection

end module gsh_tests_mod
