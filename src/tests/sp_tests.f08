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

        !call cpu_time(tstart)
        !isValid = test_sp_instantaneous_burst()
        !call cpu_time(tend)
        !write(u, '(a,f7.3,a,l)') 'test_sp_instantaneous_burst (', tend-tstart, ' sec): ', isValid

        call cpu_time(tstart)
        isValid = test_sp_constant_SFR()
        call cpu_time(tend)
        write(u, '(a,f7.3,a,l)') 'test_sp_constant_SFR (', tend-tstart, ' sec): ', isValid

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
        real(kind=rkd), parameter      :: evolTime = 3.d0             ! CU [Gyr]
        real(kind=rkd)                 :: t
        real(kind=rkd)                 :: adt
        real(kind=rkd)                 :: solution, diff

        type(gas)                      :: inRate   ! The constant SFR
        type(gas)                      :: outRate  ! wind/sn ejection rate
        type(gas)                      :: ejGas    ! Gas ejected by stellar population
        type(sp)                       :: aSp      ! A stellar population

        isValid = .TRUE.

        ! Init inRate
        inRate = real(10.d0 * MassRate_CU, kind=rkd) * initAbund(3)  ! 10Msun/yr in CU
        ! Init outRate
        call outRate%create()

        ! Create the stellar population
        call aSp%create()
        ! Create gas reservoir
        call ejGas%create()

        ! Open data files for this test
        write(filename, '(a,a,i2.2,a)') trim(validPath), '/sp_test_instantaneous_burst.dat'
        open(unit=u, file=filename, status='new')
        write(u, '(a)') '# t | stellar mass [CU] | gas mass [CU] | solution | diff | mAge'

        ! Evolution
        t = 0.d0 ! init
        solution = 0.
        diff = 0.
        do while (t < evolTime)
            !
            write(u, *) t, aSp%mass, ejGas%mass, solution, diff, aSp%mAge
            !
            if (t > 0.) then
                inRate = real(0.d0 * MassRate_CU, kind=rkd) * initAbund(3)  ! 10Msun/yr in CU
            end if
            !
            ! Compute real solution
            solution = solution + inRate%mass * dt
            !
            ! Compute evolution
            adt = dt
            call aSp%evolve(adt, inRate, outRate)
            !
            ! Update ejected gas reservoir
            ejGas = ejGas + adt * outRate
            !
            ! Test, mass conservation
            diff = abs(ejGas%mass + aSp%mass - solution)
            if (diff > num_accuracy) then
                isValid = .FALSE.
            end if
            t = t + adt
        end do

        ! Close data files
        close(u)

        ! Delete structure
        call aSp%delete()

    end function test_sp_instantaneous_burst

    ! **********************************
    function test_sp_constant_SFR() result(isValid)

        ! Test stellar population evolution according to a constant SFR

        implicit none

        integer(kind=ikd), parameter   :: u = 10000        ! file unit

        logical                        :: isValid

        character(MAXPATHSIZE)         :: filename

        real(kind=rkd), parameter      :: dt = real(1.d-4, kind=rkd)  ! CU [Gyr]
        real(kind=rkd), parameter      :: evolTime = 1.d0             ! CU [Gyr]
        real(kind=rkd)                 :: t
        real(kind=rkd)                 :: adt
        real(kind=rkd)                 :: solution, diff

        type(gas)                      :: inRate   ! The constant SFR
        type(gas)                      :: outRate  ! SN ejection rate
        type(gas)                      :: ejGas    ! Gas ejected by stellar population
        type(sp)                       :: aSp      ! A stellar population

        isValid = .TRUE.

        ! Init inRate
        inRate = real(1.d1 * MassRate_CU, kind=rkd) * initAbund(3)  ! 10Msun/yr in CU
        ! Init outRate
        call outRate%create()

        ! Create the stellar population
        call aSp%create()
        ! Create gas reservoir
        call ejGas%create()

        ! Open data files for this test
        write(filename, '(a,a,i2.2,a)') trim(validPath), '/sp_test_constant_SFR.dat'
        open(unit=u, file=filename, status='new')
        write(u, '(a)') '# t | stellar mass [CU] | gas mass [CU] | solution | diff | mAge'

        ! Evolution
        solution = 0.d0
        diff = 0.d0
        t = 0.d0 ! init
        do while (t < evolTime)
            !
            write(u, *) t, aSp%mass, ejGas%mass, solution, diff, aSp%mAge
            !
            ! Compute evolution
            adt = dt
            call aSp%evolve(adt, inRate, outRate)
            !
            ! Compute real solution
            solution = solution + inRate%mass * adt
            !
            ! Update ejected gas reservoir
            ejGas = ejGas + adt * outRate
            !
            ! Test, mass conservation
            diff = abs(ejGas%mass + aSp%mass - solution)
            if (diff > num_accuracy) then
                isValid = .FALSE.
            end if
            t = t + adt
        end do

        ! Close data files
        close(u)

        ! Delete structure
        call aSp%delete()

    end function test_sp_constant_SFR

end module sp_tests_mod
