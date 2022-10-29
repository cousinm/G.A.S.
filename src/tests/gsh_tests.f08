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

        integer(kind=ikd)              :: s
        integer(kind=ikd), parameter   :: u = 10000        ! file unit

        logical                      :: isValid

        character(MAXPATHSIZE)       :: filename

        real(kind=rkd), parameter    :: l = real(1.d2, kind=rkd)*pc2kpc  ! [pc] Injection scale
        real(kind=rkd), parameter    :: dt = real(1.d-4, kind=rkd)       ! CU [Gyr]
        real(kind=rkd), parameter    :: evolTime = real(1.d0, kind=rkd)  ! CU [Gyr]
        real(kind=rkd)               :: t

        type(gas)                    :: inRate  ! The constant injection rate at scale l
        type(gas), allocatable       :: outRates(:) 
        type(gsh)                    :: gs      ! gas structuration history

        isValid = .FALSE.

        ! Create input rate
        inRate = real(1.d1, kind=rkd) * MassRate_CU * initAbund(nMetBins)  ! 10Msun/yr in CU

        ! Create set of output rates (null)
        allocate(outRates(nScales))
        do s = 1, nScales
            call outRates(s)%create()
        end do

        ! Create a scale
        call gs%create()

        ! Open data files for this test
        ! Save each scale in a dedicated file
        do s = 1, nScales
            write(filename, '(a,a,i2.2,a)') trim(validPath), '/gsh_test_constant_injection_s', s, '.dat'
            open(unit=u+s, file=filename, status='new')
            write(u+s, '(a)') '# t | mass | ll | Vesc [km/s] | nClouds(s) | mass(s)'
        end do

        ! Evolution
        t = 0.d0 ! init
        do while (t < evolTime)
            ! Save each scale in a dedicated file
            do s = 1, nScales
                write(u+s, *) t, gs%mass, gs%l, gs%Vesc()*Velocity_km_s, gs%cascade(s)%nClouds(), gs%cascade(s)%gas%mass
            end do
            call gs%evolve(dt, l, inRate, outRates)
            t = t + dt
        end do

        ! Close data files
        do s = 1, nScales
            close(u+s)
        end do

        ! Delete structures
        call gs%delete()
        call inRate%delete()
        do s = 1, nScales
            call outRates(s)%delete()
        end do
        deallocate(outRates)

        isValid = .TRUE.

    end function test_gsh_constant_injection

    ! **********************************
    function test_gsh_constant_and_stop_injection() result(isValid)

    implicit none

        integer(kind=ikd)              :: s
        integer(kind=ikd), parameter   :: u = 10000        ! file unit

        logical                      :: isValid

        character(MAXPATHSIZE)       :: filename

        real(kind=rkd), parameter    :: l = real(1.d2, kind=rkd)*pc2kpc  ! [pc] Injection scale
        real(kind=rkd), parameter    :: dt = real(1.d-4, kind=rkd)       ! CU [Gyr]
        real(kind=rkd), parameter    :: evolTime = real(1.d0, kind=rkd)  ! CU [Gyr]
        real(kind=rkd)               :: t

        type(gas)                    :: inRate  ! The constant injection rate at scale l
        type(gas), allocatable       :: outRates(:) 
        type(gsh)                    :: gs      ! gas structuration history

        isValid = .FALSE.

        ! Create input rate
        inRate = real(1.d1, kind=rkd) * MassRate_CU * initAbund(nMetBins)  ! 10Msun/yr in CU

        ! Create set of output rates (null)
        allocate(outRates(nScales))
        do s = 1, nScales
            call outRates(s)%create()
        end do

        ! Create a scale
        call gs%create()

        ! Open data files for this test
        ! Save each scale in a dedicated file
        do s = 1, nScales
            write(filename, '(a,a,i2.2,a)') trim(validPath), '/gsh_test_constant_and_stop_injection_s', s, '.dat'
            open(unit=u+s, file=filename, status='new')
            write(u+s, '(a)') '# t | mass | ll | Vesc [km/s] | nClouds(s) | mass(s)'
        end do

        ! Evolution
        t = 0.d0 ! init
        do while (t < evolTime)
            ! Save each scale in a dedicated file
            do s = 1, nScales
                write(u+s, *) t, gs%mass, gs%l, gs%Vesc()*Velocity_km_s, gs%cascade(s)%nClouds(), gs%cascade(s)%gas%mass
            end do
            if (t > 8.d-1) inRate = real(0.d1, kind=rkd) * initAbund(nMetBins)
            call gs%evolve(dt, l, inRate, outRates)
            t = t + dt
        end do

        ! Close data files
        do s = 1, nScales
            close(u+s)
        end do

        ! Delete structures
        call gs%delete()
        call inRate%delete()
        do s = 1, nScales
            call outRates(s)%delete()
        end do
        deallocate(outRates)

        isValid = .TRUE.

    end function test_gsh_constant_and_stop_injection

end module gsh_tests_mod
