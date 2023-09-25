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
    use solver_mod ! Acces to solver parameters

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
        isValid = test_ssp_instantaneous_burst()
        call cpu_time(tend)
        write(u, '(a,f7.3,a,l)') 'test_ssp_instantaneous_burst (', tend-tstart, ' sec): ', isValid

        ! call cpu_time(tstart)
        ! isValid = test_ssp_constant_injection_and_stop()
        ! call cpu_time(tend)
        ! write(u, '(a,f7.3,a,l)') 'test_ssp_constant_injection_and_stop (', tend-tstart, ' sec): ', isValid

        ! Close log file
        close(u)

        return
    end subroutine ssp_all_tests

    !
    ! FUNCTIONS
    !

    ! **********************************
    function test_ssp_instantaneous_burst() result(isValid)

        implicit none

        integer(kind=ikd), parameter :: u = 334          ! file unit
        integer(kind=ikd), parameter :: v = 335          ! file unit
        integer(kind=ikd), parameter :: w = 336          ! file unit
        integer(kind=ikd), parameter :: iAge = 65
        integer(kind=ikd), parameter :: jMet = 3

        logical                      :: isValid

        character(MAXPATHSIZE)       :: filename

        real(kind=rkd)               :: adt
        real(kind=rkd), parameter    :: evolTime = 1.d1  ! CU [Gyr]
        real(kind=rkd), parameter    :: M = real(1.d-1, kind=rkd)   ! CU [10^11Msun]
        real(kind=rkd)               :: t
        real(kind=rkd)               :: solution, diff

        type(gas)                    :: inRate  ! The constant injection rate
        type(gas)                    :: ejGas   ! Gas ejected from the ssp
        type(gas)                    :: trGas   ! Gas transfered to the next age bin

        type(ssp)                    :: aSsp    ! A ssp

        isValid = .TRUE.

        ! Create a ssp
        call aSsp%create(iAge, jMet)
        ! Set inRate to null
        call inRate%create()
        ! Create gas reservoir
        call ejGas%create()
        call trGas%create()

        ! Open data files for this test
        write(filename,'(a, a)') trim(validPath), '/test_ssp_instantaneous_burst_mass_conservation.dat'
        open(unit=u, file=filename, status='new')
        write(u, '(a)') '# time | Ssp mass | Ej gas mass | Tr gas mass | Solution | Error'
        !
        write(filename,'(a, a)') trim(validPath), '/test_ssp_instantaneous_burst_rates.dat'
        open(unit=v, file=filename, status='new')
        write(v, '(a)') '# time | In Rate | Tr Rate | Out Rate '
        !
        write(filename,'(a, a)') trim(validPath), '/test_ssp_instantaneous_burst_timescales.dat'
        open(unit=w, file=filename, status='new')
        write(w, '(a)') '# time | Age Bin | Next Age Bin  | Average Age | Formation Age'


        ! Evolution
        t = 0.d0      ! init
        solution = M  ! init
        diff = 0.d0   ! init
        
        do while (t < evolTime)
            !
            write(u, *) t, aSsp%mass, ejGas%mass, trGas%mass, solution, diff
            write(v, *) t, sspFinalStatus%in%mass, sspFinalStatus%tr%mass, sspFinalStatus%out%mass
            write(w, *) t, ageBins(iAge), ageBins(iAge + 1), aSsp%avgAge, aSsp%tform
            !
            adt = solver_dt
            !
            if (t < num_accuracy) then
                ! Add mass in a very short event
                adt = real(solver_dt / 10., kind=rkd)
                inRate = M / adt * initAbund(jMet)
            else
                ! No new mass
                call inRate%create()
            end if
            !
            ! Compute evolution
            call aSsp%evolve(adt, inRate)
            !
            ! Update ejected gas reservoir
            ejGas = ejGas + adt * sspFinalStatus%out
            ! Update transfered mass reservoir
            trGas = trGas + adt * sspFinalStatus%tr
            !
            ! Test, mass conservation
            diff = abs(ejGas%mass + aSsp%mass + trGas%mass - solution)
            if (diff > num_accuracy) then
                isValid = .FALSE.
            end if
            t = t + adt
        end do

        ! Close data file
        close(u)
        close(v)
        close(w)

        ! Delete structure
        call inRate%delete()
        call ejGas%delete()
        call trGas%delete()

    end function test_ssp_instantaneous_burst

    ! **********************************
    function test_ssp_constant_injection_and_stop() result(isValid)

        implicit none

        integer(kind=ikd), parameter :: u = 334          ! file unit

        integer(kind=ikd), parameter :: iAge = 30

        logical                      :: isValid

        character(MAXPATHSIZE)       :: filename

        real(kind=rkd)               :: adt
        real(kind=rkd), parameter    :: evolTime = 1.d0  ! CU [Gyr]
        real(kind=rkd)               :: t
        real(kind=rkd)               :: solution, diff

        type(gas)                    :: inRate  ! The constant injection rate
        type(gas)                    :: ejGas   ! Gas ejected from the ssp
        type(gas)                    :: trGas   ! Gas transfered to the next age bin

        type(ssp)                    :: aSsp    ! A ssp

        isValid = .TRUE.

        inRate = 1.d1 * initAbund(nMetBins) * Msun_Yr2MassRateCU  ! 10 Msun/yr in CU

        ! Create a ssp
        call aSsp%create(iAge, nMetBins)
        ! Create gas reservoir
        call ejGas%create()
        call trGas%create()

        ! Open data file for this test
        write(filename,'(a, a)') trim(validPath), '/ssp_test_constant_injection_and_stop_mass_conservation.dat'
        open(unit=u, file=filename, status='new')

        ! Evolution
        t = 0.d0         ! init
        solution = 0.d0  ! init
        diff = 0.d0      ! init
        write(u, '(a)') '# time | Ssp mass | Ej gas mass | Tr gas mass | Solution | Error'
        do while (t < evolTime)
            !
            if (t > 0.2) then
                call inRate%create()
            end if
            !
            write(u, *) t, aSsp%mass, ejGas%mass, trGas%mass, solution, diff
            !
            ! Compute evolution
            adt = solver_dt
            call aSsp%evolve(adt, inRate)
            !
            ! Compute real solution
            solution = solution + inRate%mass * adt
            !
            ! Update ejected gas reservoir
            ejGas = ejGas + adt * sspFinalStatus%out
            ! Update transfered mass reservoir
            trGas = trGas + adt * sspFinalStatus%tr
            !
            ! Test, mass conservation
            diff = abs(ejGas%mass + aSsp%mass + trGas%mass - solution)
            if (diff > num_accuracy) then
                isValid = .FALSE.
            end if
            t = t + adt
        end do

        ! Close data file
        close(u)

        ! Delete structure
        call ejGas%delete()
        call inRate%delete()
        call trGas%delete()

    end function test_ssp_constant_injection_and_stop

end module ssp_tests_mod
