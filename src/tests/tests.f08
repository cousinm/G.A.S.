program GAS_tests

    ! Run G.A.S tests

    use config_mod      ! Acces to configurations parameters (path)
    use model_mod       ! Acces to model parameters and integration scheme configuration
    use gas_tests_mod   ! gas tests
    use scale_tests_mod ! scale tests
    use gsh_tests_mod   ! gas tests
    use sp_tests_mod    ! stellar population

    implicit none

    ! Initialize the configuration
    call config_init()
    ! Initialize model configuration and integration scheme
    call model_init()
    ! Initialize modules
    call gsh_init()
    ! Run all "gas" object tests
    !call gas_all_tests()
    ! Run all "scale" object tests
    !call scale_all_tests()
    ! Run all "gsh" object tests
    !call gsh_all_tests()
    ! Run all "sp" object tests
    call sp_all_tests()
    
end program GAS_tests