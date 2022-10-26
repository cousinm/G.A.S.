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

        ! Init module
        call sp_init()

        ! Open log file for these tests
        write(filename,'(a, a)') trim(validPath), '/sp_tests.log'
        open(unit=u, file=filename, status='new')

        isValid = test_sp_new()
        write(u, '(a,l)') 'test_sp_new: ', isValid

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
        real(kind=rkd), allocatable :: m(:)

        type(gas)                   :: g

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
        do i = 1, nMetBins
            if ((m(i) > 0.) .and. (im(i) < 1.)) then
                isValid = .FALSE.
                return
            end if
        end do
        ! Test total mass
        if (abs(sum(m) - g%mass) .gt. num_accuracy) then
            isValid = .FALSE.
            return
        end if

        deallocate(m)

        isValid = .TRUE.

    end function test_sp_new

    ! **********************************
    function test_sp_constant_SFR() result(isValid)

        ! Test stellar population evolution according to a constant SFR

        implicit none

        integer(kind=ikd), parameter   :: u = 10000        ! file unit

        logical                        :: isValid

        character(MAXPATHSIZE)         :: filename

        real(kind=rkd), parameter      :: dt = 1.d-4       ! CU [Gyr]
        real(kind=rkd), parameter      :: evolTime = 1.d0  ! CU [Gyr]
        real(kind=rkd)                 :: t

        type(gas)                      :: inRate  ! The constant SFR
        type(sp)                       :: esp     ! The evolvinf stellar population

        isValid = .FALSE.

        inRate = real(1.d1 * MassRate_CU, kind=rkd) * initAbund(3)  ! 10Msun/yr in CU

        ! Create a scale
        call esp%create()

        ! Open data files for this test
        write(filename, '(a,a,i2.2,a)') trim(validPath), '/sp_test_constant_SFR.dat'
        open(unit=u, file=filename, status='new')
        write(u, '(a)') '# t | mass [CU] |'

        ! Evolution
        t = 0.d0 ! init
        do while (t < evolTime)
            write(u, *) t, esp%mass
            call esp%solve(dt, inRate)
            t = t + dt
        end do

        ! Close data files
        close(u)

        ! Delete structure
        call esp%delete()

        isValid = .TRUE.

    end function test_sp_constant_SFR

end module sp_tests_mod
