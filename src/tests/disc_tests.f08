module disc_test_mod

    use config_mod ! Acces to configurations parameters (path)
    use log_mod    ! Acces to logging procedures
    use disc_mod   ! Acces to disc methods and procedures

    !*****************************************************************************************************************
    !
    !  OVERVIEW
    !
    !  This module is dedicated to disc module tests.
    !
    !*****************************************************************************************************************

    implicit none

    public

    contains

    !
    ! SUBROUTINES
    !

    ! **********************************
    subroutine disc_all_tests()

        ! Run all tests for disc module
    
        implicit none

        integer(kind=ikd), parameter :: u = 333       ! file unit

        logical                      :: isValid

        real(kind=rkd)               :: tstart, tend

        character(MAXPATHSIZE)     :: filename

        ! Open log file for these tests
        write(filename,'(a, a)') trim(validPath), '/disc_tests.log'
        open(unit=u, file=filename, status='new')

        call cpu_time(tstart)
        isValid = test_disc_constant_injection_and_stop()
        call cpu_time(tend)
        write(u, '(a,f7.3,a,l)') 'test_disc_constant_injection_and_stop (', tend-tstart, ' sec): ', isValid

        ! Close log file
        close(u)

        return
    end subroutine disc_all_tests
    !
    ! FUNCTIONS
    !

    ! **********************************
    function test_disc_create() result(isValid)

        logical           :: isValid

        type(disc)         :: d

        ! Create a disc object
        call d%create()

        ! Tests
        isValid = .TRUE.

    end function test_disc_create

    ! **********************************
    function test_disc_delete() result(isValid)

        logical           :: isValid

        type(disc)         :: d

        ! Create a disc object
        call d%create()
        ! Delete it
        call d%delete()

        ! Tests
        isValid = .TRUE.

    end function test_disc_delete

    ! **********************************
    function test_disc_constant_injection_and_stop() result(isValid)

        implicit none

        integer(kind=ikd), parameter   :: u = 10000        ! file unit
        integer(kind=ikd), parameter   :: v = 10001        ! file unit

        logical                      :: isValid

        character(MAXPATHSIZE)       :: filename

        real(kind=rkd), parameter    :: l = real(1.d2, kind=rkd)*pc2LenCU  ! [pc] Injection scale
        real(kind=rkd), parameter    :: dt = real(1.d-4, kind=rkd)       ! CU [Gyr]
        real(kind=rkd)               :: adt
        real(kind=rkd), parameter    :: evolTime = real(1.d0, kind=rkd)  ! CU [Gyr]
        real(kind=rkd)               :: t
        real(kind=rkd)               :: solution, diff

        type(gas)                    :: inRate       ! The constant injection rate at scale l
        type(gas)                    :: ejGas        ! Ejected gas

        type(disc)                   :: aDisc        ! A disc

        isValid = .TRUE.

        ! Create input rate
        inRate = real(10., kind=rkd) * initAbund(4) * Msun_Yr2MassRateCU ! 10Msun/yr in CU

        ! Create a disc
        call aDisc%create()
        ! Create gas reservoir
        call ejGas%create()

        ! Open data files for this test
        ! Mass conservation data file
        write(filename, '(a,a,i2.2,a)') trim(validPath), '/disc_test_constant_injection_and_stop.dat'
        open(unit=u, file=filename, status='new')
        write(u, '(a)') '# time | Gas mass | Star mass | Ejected gas mass | Solution | Error'
        !
        ! Results file
        write(filename, '(a,a,i2.2,a)') trim(validPath), '/disc_results_constant_injection_and_stop.dat'
        open(unit=v, file=filename, status='new')
        write(v, '(a)') '# time | inRate | SFR | outRate | Solution | Error'
        !
        ! Evolution
        t = 0.d0         ! init
        solution = 0.d0  ! init
        diff = 0.d0      ! init
        do while (t < evolTime)
            !
            if (t > 0.05) then
                call inRate%create()
            end if
            !
            ! Save
            write(u, *) t, aDisc%myGsh%mass, aDisc%mySp%mass, ejGas%mass, solution, diff
            write(v, *) t, aDisc%myGsh%myStatus%in%mass, &
                        aDisc%mySp%myStatus%in%mass, aDisc%myGsh%myStatus%out%mass, solution, diff
            !
            ! Compute evolution
            adt = dt
            call aDisc%evolve(adt, l, inRate)
            if (adt < dt) write(v, *) 'adatative time-step:', adt
            !
            ! Compute real solution
            solution = solution + inRate%mass * adt
            !
            ! Update ejected gas reservoir
            ejGas = ejGas + adt * aDisc%myGsh%myStatus%out
            !
            ! Test, mass conservation
            diff = aDisc%myGsh%mass + aDisc%mySp%mass + ejGas%mass - solution
            if (diff > num_accuracy) then
                isValid = .FALSE.
            end if
            t = t + adt
        end do

        ! Close data file
        close(u)
        close(v)

        ! Delete structures
        call aDisc%delete()
        call inRate%delete()
        call ejGas%delete()

    end function test_disc_constant_injection_and_stop

    ! **********************************

end module disc_test_mod