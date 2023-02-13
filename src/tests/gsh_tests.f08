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

        logical                      :: isValid

        real(kind=rkd)               :: tstart, tend

        character(MAXPATHSIZE)       :: filename

        ! Open log file for these tests
        write(filename,'(a, a)') trim(validPath), '/gsh_tests.log'
        open(unit=u, file=filename, status='new')

        call cpu_time(tstart)
        isValid = test_gsh_i2l_and_l2i()
        call cpu_time(tend)
        write(u, '(a,f7.3,a,l)') 'test_gsh_i2l_and_l2i (', tend-tstart, ' sec): ', isValid

        isValid = test_gsh_create()
        write(u, '(a,l)') 'test_gsh_create: ', isValid

        isValid = test_gsh_delete()
        write(u, '(a,l)') 'test_gsh_delete: ', isValid

        call cpu_time(tstart)
        isValid = test_gsh_constant_injection_and_stop()
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

    function test_gsh_i2l_and_l2i() result(isValid)

        implicit none

        integer(kind=ikd) :: il, il_test

        logical           :: isValid

        real(kind=rkd)    :: li, li_test

        isValid = .TRUE.
        il = 1
        li = gsh_i2l(il)
        il_test = gsh_l2i(li)
        if (il_test /= il) then
            isValid = .FALSE.
            return
        end if

        li = real(1.d2, kind=rkd)*pc2kpc
        il = gsh_l2i(li)
        li_test = gsh_i2l(il)
        if (il < nScales .and. li_test < li)  then
            isValid = .FALSE.
            return
        end if

        return
    end function test_gsh_i2l_and_l2i

    ! **********************************
    function test_gsh_constant_injection_and_stop() result(isValid)

        implicit none

        integer(kind=ikd)              :: s
        integer(kind=ikd), parameter   :: u = 10000        ! file unit

        logical                      :: isValid

        character(MAXPATHSIZE)       :: filename

        real(kind=rkd), parameter    :: l = real(1.d2, kind=rkd)*pc2kpc  ! [pc] Injection scale
        real(kind=rkd), parameter    :: dt = real(1.d-4, kind=rkd)       ! CU [Gyr]
        real(kind=rkd), parameter    :: evolTime = real(2.d0, kind=rkd)  ! CU [Gyr]
        real(kind=rkd), parameter    :: Q = real(0.d0, kind=rkd)         ! Gas injection power (here = 0.)
        real(kind=rkd)               :: t
        real(kind=rkd)               :: adt
        real(kind=rkd)               :: solution, diff

        type(gas)                    :: inRate  ! The constant injection rate at scale l
        type(gas)                    :: ejGas   ! Gas ejected from the gsh

        type(gsh)                    :: aGsh    ! gas structuration history

        isValid = .TRUE.

        ! Create input rate
        inRate = real(1.d1, kind=rkd) * MassRate_CU * initAbund(nMetBins)  ! 10Msun/yr in CU

        ! Create a gas structuration history
        call aGsh%create()
        ! Create gas reservoir
        call ejGas%create()

        ! Open data files for this test
        ! Save each scale in a dedicated file
        do s = 1, nScales
            write(filename, '(a,a,i2.2,a)') trim(validPath), '/gsh_results_constant_injection_and_stop_s', s, '.dat'
            open(unit=u+s, file=filename, status='new')
            write(u+s, '(a)') '# time | mass | il | Vesc [km/s] | nClouds(s) | mass(s) | Ms '
        end do
        ! Create a other main test file for mass conservation test
        write(filename, '(a,a)') trim(validPath), '/gsh_test_constant_injection_and_stop.dat'
        open(unit=u, file=filename, status='new')
        write(u, '(a)') '# time | Gsh mass | Ej gas mass | Solution | Error'

        ! Evolution
        t = 0.d0         ! init
        solution = 0.d0  ! init
        diff = 0.d0      ! init
        do while (t < evolTime)
            !
            if (t > 2.d-1) then
                call inRate%create()
            end if
            !
            ! Save each scale in a dedicated file
            do s = 1, nScales
                write(u+s, *) t, aGsh%mass, aGsh%ih, aGsh%Vesc()*Velocity_km_s, aGsh%cascade(s)%nClouds(), aGsh%cascade(s)%gas%mass, ejGas%mass
            end do
            write(u, *) t, aGsh%mass, ejGas%mass, solution, diff
            !
            ! Compute evolution
            adt = dt
            call aGsh%evolve(adt, l, Q, inRate)
            !
            ! Compute real solution
            solution = solution + inRate%mass * adt
            !
            ! Update ejected gas reservoir
            ejGas = ejGas + adt * (aGsh%myStatus%tr + aGsh%myStatus%out)
            !
            ! Test, mass conservation
            diff = ejGas%mass + aGsh%mass - solution
            if (diff > num_accuracy) then
                isValid = .FALSE.
            end if
            t = t + adt
        end do

        ! Close data files
        do s = 1, nScales
            close(u+s)
        end do
        close(u)

        ! Delete structures
        call aGsh%delete()
        call inRate%delete()
        call ejGas%delete()

    end function test_gsh_constant_injection_and_stop

end module gsh_tests_mod
