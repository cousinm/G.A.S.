module gas_tests_mod

    !*****************************************************************************************************************
    !
    !  OVERVIEW
    !
    !  This module is dedicated to gas module tests.
    !
    !*****************************************************************************************************************

    use gas_mod    ! Acces to gas module
    use log_mod    ! Acces to logging procedures

    implicit none

    public

contains

    !
    ! SUBROUTINES
    !

    ! **********************************
    subroutine gas_all_tests()

        ! Run all tests for gas module
    
        implicit none

        integer(kind=4), parameter :: u = 333       ! file unit

        logical                    :: isValid

        character(MAXPATHSIZE)     :: filename

        ! Open log file for these tests
        write(filename,'(a, a)') trim(logPath), '/gas_tests.log'
        open(unit=u, file=filename, status='new')

        isValid = test_gas_create()
        write(u, '(a,l)') 'gas_create: ', isValid

        isValid = test_gas_temperature()
        write(u, '(a,l)') 'test_gas_temperature: ', isValid

        isValid = test_gas_copy()
        write(u, '(a,l)') 'gas_copy: ', isValid

        isValid = test_gas_add()
        write(u, '(a,l)') 'gas_add: ', isValid

        isValid = test_gas_scalar_multiply_left()
        write(u, '(a,l)') 'test_gas_scalar_multiply_left: ', isValid

        isValid = test_gas_scalar_multiply_right()
        write(u, '(a,l)') 'test_gas_scalar_multiply_right: ', isValid

        isValid = test_gas_sub()
        write(u, '(a,l)') 'gas_sub: ', isValid

        ! Close log file
        close(u)

        return
    end subroutine gas_all_tests

    !
    ! FUNCTIONS
    !

    ! **********************************
    function test_gas_create() result(isValid)

        ! Test the gas creation procedure

        implicit none

        integer(kind=4)   :: e

        logical           :: isValid

        type(gas)         :: g

        ! Create a gas object
        call g%create()

        ! Tests
        isValid = .TRUE.
        if ((g%mass < 0.) .or. (g%mass > 0.)) then
            isValid = .FALSE.
        else if ((g%mZ < 0.) .or. (g%mZ > 0.)) then
            isValid = .FALSE.
        else if ((g%Eint < 0.) .or. (g%Eint > 0.)) then
            isValid = .FALSE.
        else if (.not. allocated(g%elts)) then
            isValid = .FALSE.
        else
            do e = 1, nElts
                if ((g%elts(e) < 0.) .or. (g%elts(e) > 0.)) then
                    isValid = .FALSE.
                end if
            end do
        end if
    end function test_gas_create

    ! **********************************
    function test_gas_temperature() result(isValid)

        ! Test temperature setting and the associated total energy conversion

        implicit none

        logical                 :: isValid

        real(kind=8), parameter :: T_ref = 1.d5
        real(kind=8)            :: T

        type(gas)               :: g

        ! Define a gas object
        g = initAbund(nMetBins)
        ! Set the gas temperature to T_ref
        call g%setTemperature(T_ref)

        ! Tests 
        isValid = .TRUE.
        if (g%Eint <= 0.d0) then
            ! The total internal energy is not valid
            isValid = .FALSE.
        else
            ! Get back the temperature
            T = g%temperature()
            if ((T > T_ref) .or. (T < T_ref)) then
                isValid = .FALSE.
            end if
        end if

    end function test_gas_temperature

    ! **********************************
    function test_gas_copy() result(isValid)
    
        ! Test the gas copy procedure

        implicit none
        
        integer(kind=4)         :: e

        logical                 :: isValid

        real(kind=8), parameter :: T_ref = 1.d5  ! [K]

        type(gas)               :: g1
        type(gas)               :: g2

        ! Create a gas object from initial abundance gas object
        g1 = initAbund(nMetBins)
        ! Set the temperature to T_ref
        call g1%setTemperature(T_ref)

        ! Copy
        g2 = g1

        ! Tests
        isValid = .TRUE.
        if ((g2%mass < g1%mass) .or. (g2%mass > g1%mass)) then
            isValid = .FALSE.
        else if ((g2%mZ < g1%mZ) .or. (g2%mZ > g1%mZ)) then
            isValid = .FALSE.
        else if ((g2%Eint < g1%Eint) .or. (g2%Eint > g1%Eint)) then
            isValid = .FALSE.
        else if (.not. allocated(g2%elts)) then
            isValid = .FALSE.
        else
            do e = 1, nElts
                if ((g2%elts(e) < g1%elts(e)) .or. (g2%elts(e) > g1%elts(e))) then
                    isValid = .FALSE.
                end if
            end do
        end if

    end function test_gas_copy

    ! **********************************
    function test_gas_add() result(isValid)

        ! Test the gas add procedure
        ! Start from a empty gas and add the initial composition
        ! Of the highest metalicity bin (all elements mass > 0.)

        implicit none

        integer(kind=4)  :: e

        logical          :: isValid

        type(gas)        :: g1
        type(gas)        :: g2

        ! Create a gas object
        call g1%create()

        g2 = g1 + initAbund(nMetBins)

        ! Tests
        isValid = .TRUE.
        if ((g2%mass < initAbund(nMetBins)%mass) .or. (g2%mass > initAbund(nMetBins)%mass)) then
            isValid = .FALSE.
        else if ((g2%mZ < initAbund(nMetBins)%mZ) .or. (g2%mZ > initAbund(nMetBins)%mZ)) then
            isValid = .FALSE.
        else if ((g2%Eint < initAbund(nMetBins)%Eint) .or. (g2%Eint > initAbund(nMetBins)%Eint)) then
            isValid = .FALSE.
        else if (.not. allocated(g2%elts)) then
            isValid = .FALSE.
        else
            do e = 1, nElts
                if ((g2%elts(e) < initAbund(nMetBins)%elts(e)) .or. (g2%elts(e) > initAbund(nMetBins)%elts(e))) then
                    isValid = .FALSE.
                end if
            end do
        end if

    end function test_gas_add

    ! **********************************
    function test_gas_scalar_multiply_left() result(isValid)

        ! Test the scalar multiplication of a gas object: a * g

        implicit none

        integer(kind=4)         :: e

        logical                 :: isValid

        real(kind=8), parameter :: a = 3.
        real(kind=8), parameter :: T_ref = 1.d5  ! [K]

        type(gas)               :: g1
        type(gas)               :: g2

        ! Create a gas object from initial abundance gas object
        g1 = initAbund(nMetBins)
        ! Set the temperature to T_ref
        call g1%setTemperature(T_ref)

        g2 = a * g1

        ! Tests
        isValid = .TRUE.
        if ((g2%mass < a*g1%mass) .or. (g2%mass > a*g1%mass)) then
            isValid = .FALSE.
        else if ((g2%mZ < a*g1%mZ) .or. (g2%mZ > a*g1%mZ)) then
            isValid = .FALSE.
        else if ((g2%Eint < a*g1%Eint) .or. (g2%Eint > a*g1%Eint)) then
            isValid = .FALSE.
        else if (.not. allocated(g2%elts)) then
            isValid = .FALSE.
        else
            do e = 1, nElts
                if ((g2%elts(e) < a*g1%elts(e)) .or. (g2%elts(e) > a*g1%elts(e))) then
                    isValid = .FALSE.
                end if
            end do
        end if

    end function test_gas_scalar_multiply_left

    ! **********************************
    function test_gas_scalar_multiply_right() result(isValid)

        ! Test the scalar multiplication of a gas object: g * a

        implicit none

        integer(kind=4)         :: e

        logical                 :: isValid

        real(kind=8), parameter :: a = 3.

        type(gas)               :: g1
        type(gas)               :: g2

        g1 = initAbund(nMetBins)
        g2 = g1 * a

        ! Tests
        isValid = .TRUE.
        if ((g2%mass < a*g1%mass) .or. (g2%mass > a*g1%mass)) then
            isValid = .FALSE.
        else if ((g2%mZ < a*g1%mZ) .or. (g2%mZ > a*g1%mZ)) then
            isValid = .FALSE.
        else if ((g2%Eint < a*g1%Eint) .or. (g2%Eint > a*g1%Eint)) then
            isValid = .FALSE.
        else if (.not. allocated(g2%elts)) then
            isValid = .FALSE.
        else
            do e = 1, nElts
                if ((g2%elts(e) < a*g1%elts(e)) .or. (g2%elts(e) > a*g1%elts(e))) then
                    isValid = .FALSE.
                end if
            end do
        end if

    end function test_gas_scalar_multiply_right

    ! **********************************
    function test_gas_sub() result(isValid)

        implicit none

        integer(kind=4)         :: e

        logical                 :: isValid

        real(kind=8), parameter :: a = 2.

        type(gas)               :: g1
        type(gas)               :: g2

        g1 = initAbund(nMetBins)
        g2 = a*g1 - g1

        ! Tests g2 = g1
        isValid = .TRUE.
        if ((g2%mass < g1%mass) .or. (g2%mass > g1%mass)) then
            isValid = .FALSE.
        else if ((g2%mZ < g1%mZ) .or. (g2%mZ > g1%mZ)) then
            isValid = .FALSE.
        else if ((g2%Eint < g1%Eint) .or. (g2%Eint > g1%Eint)) then
            isValid = .FALSE.
        else if (.not. allocated(g2%elts)) then
            isValid = .FALSE.
        else
            do e = 1, nElts
                if ((g2%elts(e) < g1%elts(e)) .or. (g2%elts(e) > g1%elts(e))) then
                    isValid = .FALSE.
                end if
            end do
        end if

    end function test_gas_sub

end module gas_tests_mod
