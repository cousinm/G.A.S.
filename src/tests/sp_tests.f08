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
    use solver_mod ! Acces to solver parameters

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

        ! call cpu_time(tstart)
        ! isValid = test_sp_constant_injection_and_stop()
        ! call cpu_time(tend)
        ! write(u, '(a,f7.3,a,l)') 'test_sp_constant_injection_and_stop (', tend-tstart, ' sec): ', isValid

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

        integer(kind=ikd)           :: i, n
        integer(kind=ikd)           :: nTests

        real(kind=rkd)              :: im(nMetBins)
        real(kind=rkd), parameter   :: totalMass = real(1.d2, kind=rkd)

        type(gas)                   :: g
        type(gas)                   :: gOut
        type(gas), allocatable      :: m(:)

        isValid = .TRUE.

        nTests = 4
        !
        ! Run ton test cases
        do n = 1, nTests
            select case (n)
            case(1)
                ! 1. Test with mass in only one bin
                im = [0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0]
                
            case(2)
                ! 2. Test with mass in tree bins
                im = [1.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0]
            case(3)
                ! 3. Test with mass in all bins
                im = [1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0]
            case(4)
                ! 4. Test with inhomogeneous distribution, realistic case
                im = [0.d0, 2.d-1, 3.d0, 1.d-1, 8.d-4, 0.d0, 0.d0]
            end select
            !
            ! Normalisation
            im = im / sum(im) * totalMass
            call g%create()
            do i = 1, nMetBins
                g = g + real(im(i), kind=rkd) * initAbund(i)
            end do
            !
            ! Solve
            m = g2s(g)
            !
            ! Test distribution
            call gOut%create()
            do i = 1, nMetBins
                if (abs(m(i)%mass - im(i))/im(i) > 1d-4) then
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
        end do
        
    end function test_sp_new

    ! **********************************
    function test_sp_instantaneous_burst() result(isValid)

        ! Test stellar population evolution following an instantaneous burst at t = 0

        implicit none

        integer(kind=ikd), parameter   :: u = 334          ! file unit
        integer(kind=ikd), parameter   :: v = 335          ! file unit
        integer(kind=ikd), parameter   :: w = 336          ! file unit

        logical                        :: isValid

        character(MAXPATHSIZE)         :: filename

        real(kind=rkd), parameter      :: M = real(1.d-5, kind=rkd)        ! CU [10^11Msun]
        real(kind=rkd), parameter      :: evolTime = real(1.d1, kind=rkd)  ! CU [Gyr]
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
        !
        ! Create gas reservoir
        call ejGas%create()

        ! Open data files for this test
        write(filename, '(a,a,i2.2,a)') trim(validPath), '/test_sp_instantaneous_burst_mass_conservation.dat'
        open(unit=u, file=filename, status='new')
        write(u, '(a)') '# time | Stellar mass | Gas mass | Solution | Error '
                !
        write(filename,'(a, a)') trim(validPath), '/test_sp_instantaneous_burst_rates.dat'
        open(unit=v, file=filename, status='new')
        write(v, '(a)') '# time | In Rate | Out Rate '
        !
        write(filename,'(a, a)') trim(validPath), '/test_sp_instantaneous_burst_timescales.dat'
        open(unit=w, file=filename, status='new')
        write(w, '(a)') '# time | Average Age '

        ! Evolution
        t = 0.d0 ! init
        solution = M
        diff = 0.
        do while (t < evolTime)
            !
            write(u, *) t, aSp%mass, ejGas%mass, solution, diff
            write(v, *) t, aSp%myStatus%in%mass, aSp%myStatus%out%mass
            write(w, *) t, aSp%mAge
            !
            adt = solver_dt
            !
            if (t < num_accuracy) then
                ! Add mass in a very short event
                adt = real(solver_dt / 10., kind=rkd)
                inRate = M / adt * initAbund(3)
            else
                ! No new mass
                call inRate%create()
            end if
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
        close(v)
        close(w)

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

        real(kind=rkd), parameter      :: dt = real(5.d-4, kind=rkd)        ! CU [Gyr]
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
