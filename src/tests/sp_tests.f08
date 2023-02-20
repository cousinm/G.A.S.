module sp_tests_mod

    !*****************************************************************************************************************
    !
    !  OVERVIEW
    !
    !  This module is dedicated to stellar population module tests.
    !
    !*****************************************************************************************************************

    use sp_mod     ! Acces to stellar population procedure and parameters
    use config_mod ! Acces to configurations parameters (path)
    use log_mod    ! Acces to logging procedures

    implicit none

    public

contains

    !
    ! SUBROUTINES
    !

    ! **********************************
    subroutine sp_all_tests()

        ! Run all tests for sp module
    
        implicit none

        integer(kind=ikd), parameter :: u = 333       ! file unit

        real(kind=rkd)             :: tstart, tend

        logical                    :: isValid

        character(MAXPATHSIZE)     :: filename

        ! Open log file for these tests
        write(filename,'(a, a)') trim(validPath), '/sp_tests.log'
        open(unit=u, file=filename, status='new')

        isValid = test_sp_new()
        write(u, '(a,l)') 'test_sp_new: ', isValid

        call cpu_time(tstart)
        isValid = test_sp_instantaneous_burst()
        call cpu_time(tend)
        write(u, '(a,f7.3,a,l)') 'test_sp_instantaneous_burst (', tend-tstart, ' sec): ', isValid

        call cpu_time(tstart)
        isValid = test_sp_constant_injection_and_stop()
        call cpu_time(tend)
        write(u, '(a,f7.3,a,l)') 'test_sp_constant_injection_and_stop (', tend-tstart, ' sec): ', isValid

        ! Close log file
        close(u)

        return
    end subroutine sp_all_tests

    !
    ! FUNCTIONS
    !

    ! **********************************
    function test_sp_new() result(isValid)

        implicit none

        logical                     :: isValid

        integer(kind=ikd)           :: i

        real(kind=rkd)              :: im(nMetBins)

        type(gas)                   :: g
        type(gas)                   :: gOut
        type(gas), allocatable      :: m(:)

        isValid = .TRUE.

        ! Test with mass only in bin 3
        im = [0., 0., 1., 0., 1., 1., 0.]
        call g%create()
        do i = 1, nMetBins
            g = g + real(im(i)*1.d-3, kind=rkd)*initAbund(i)
        end do
        ! Solve
        m = g2s(g)
        ! Test distribution
        call gOut%create()
        do i = 1, nMetBins
            if ((m(i)%mass > 0.) .and. (im(i) < 1.)) then
                isValid = .FALSE.
                return
            end if
            gOut = gOut + m(i)
        end do
        ! Test total mass
        if (abs(gOut%mass - g%mass) .gt. num_accuracy) then
            isValid = .FALSE.
            return
        end if

        deallocate(m)

    end function test_sp_new

    ! **********************************
    function test_sp_instantaneous_burst() result(isValid)

        ! Test stellar population evolution following an instantaneous burst at t = 0

        implicit none

        integer(kind=ikd), parameter   :: u = 10000        ! file unit

        logical                        :: isValid

        character(MAXPATHSIZE)         :: filename

        real(kind=rkd), parameter      :: dt = real(1.d-4, kind=rkd)  ! CU [Gyr]
        real(kind=rkd), parameter      :: M = real(1.d-5, kind=rkd)   ! CU [10^11Msun]
        real(kind=rkd), parameter      :: evolTime = 3.d0             ! CU [Gyr]
        real(kind=rkd)                 :: t
        real(kind=rkd)                 :: adt
        real(kind=rkd)                 :: solution, diff

        type(gas)                      :: inRate   ! The constant SFR
        type(gas)                      :: ejGas    ! Gas ejected by stellar population
        type(sp)                       :: aSp      ! A stellar population

        isValid = .TRUE.

        ! Init inRate
        call inRate%create()

        ! Create the stellar population
        call aSp%create()
        ! Create a 1e6 Msun ssp in the first age bin
        aSp%mass = M
        aSp%sfh(1, 3)%mass = M
        !
        ! Create gas reservoir
        call ejGas%create()

        ! Open data files for this test
        write(filename, '(a,a,i2.2,a)') trim(validPath), '/sp_test_instantaneous_burst.dat'
        open(unit=u, file=filename, status='new')
        write(u, '(a)') '# time | Stellar mass | Gas mass | Solution | Error '

        ! Evolution
        t = 0.d0 ! init
        solution = M
        diff = 0.
        do while (t < evolTime)
            !
            write(u, *) t, aSp%mass, ejGas%mass, solution, diff
            !
            ! Compute evolution
            adt = dt
            call aSp%evolve(adt, inRate)
            !
            ! Update ejected gas reservoir
            ejGas = ejGas + adt * aSp%myStatus%out
            !
            ! Test, mass conservation
            diff = ejGas%mass + aSp%mass - solution
            if (diff > num_accuracy) then
                isValid = .FALSE.
            end if
            t = t + adt
        end do

        ! Close data files
        close(u)

        ! Delete structure
        call aSp%delete()
        call inRate%delete()

    end function test_sp_instantaneous_burst

    ! **********************************
    function test_sp_constant_injection_and_stop() result(isValid)

        ! Test stellar population evolution according to a constant SFR

        implicit none

        integer(kind=ikd), parameter   :: u = 10000        ! file unit

        logical                        :: isValid

        character(MAXPATHSIZE)         :: filename

        real(kind=rkd), parameter      :: dt = real(1.d-4, kind=rkd)        ! CU [Gyr]
        real(kind=rkd), parameter      :: evolTime = real(2.d0, kind=rkd)   ! CU [Gyr]
        real(kind=rkd)                 :: t
        real(kind=rkd)                 :: adt
        real(kind=rkd)                 :: solution, diff

        type(gas)                      :: inRate   ! The constant SFR
        type(gas)                      :: ejGas    ! Gas ejected by stellar population
        type(sp)                       :: aSp      ! A stellar population

        isValid = .TRUE.

        ! Init inRate
        inRate = real(10., kind=rkd) * Msun_Yr2MassRateCU * initAbund(3)   ! 10Msun/yr in CU

        ! Create the stellar population
        call aSp%create()
        ! Create gas reservoir
        call ejGas%create()

        ! Open data files for this test
        write(filename, '(a,a,i2.2,a)') trim(validPath), '/sp_test_constant_injection_and_stop.dat'
        open(unit=u, file=filename, status='new')
        write(u, '(a)') '# time | Stellar mass | Gas mass | Solution | Error'

        ! Evolution
        solution = 0.d0
        diff = 0.d0
        t = 0.d0 ! init
        do while (t < evolTime)
            !
            if (t > 0.5) call inRate%create()
            !
            write(u, *) t, aSp%mass, ejGas%mass, solution, diff
            !
            ! Compute evolution
            adt = dt
            call aSp%evolve(adt, inRate)
            !
            ! Compute real solution
            solution = solution + inRate%mass * adt
            !
            ! Update ejected gas reservoir
            ejGas = ejGas + adt * aSp%myStatus%out
            !
            ! Test, mass conservation
            diff = ejGas%mass + aSp%mass - solution
            if (diff > num_accuracy) then
                isValid = .FALSE.
            end if
            t = t + adt
        end do

        ! Close data files
        close(u)

        ! Delete structure
        call inRate%delete()
        call ejGas%delete()
        call aSp%delete()

    end function test_sp_constant_injection_and_stop

end module sp_tests_mod
