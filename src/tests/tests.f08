program GAS_tests

    ! Run G.A.S tests

    use config_mod     ! Acces to configurations parameters (path)
    use gas_tests_mod  ! gas tests

    implicit none

    ! Initialize the configuration
    call config_init()
    ! Run all "gas" object tests
    call gas_all_tests()
    
end program GAS_tests