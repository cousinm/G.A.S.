program G_A_S

    ! main program

    use config_mod ! Acces to configurations parameters (path)
    use model_mod  ! Acces to model parameters and configurations
    use gsh_mod    ! Acces to gas structuration history module
    use sp_mod     ! Acces to stellar population module

    implicit none 

    ! Initialize the configuration
    call config_init()
    ! Initialize model configuration and integration scheme
    call model_init()
    ! Initialize modules
    call gsh_init()
    call sp_init()

end program G_A_S
