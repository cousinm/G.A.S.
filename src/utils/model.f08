module model_mod

    !*****************************************************************************************************************
    !
    !  OVERVIEW
    !
    !  This module allows to define and set all model parameters
    !  and integration scheme configuration
    !
    !*****************************************************************************************************************

    use parameters  ! Acces to global defintions and properties
    use solver_mod  ! Acces to solver parameters and methods
    use config_mod  ! Acces to configurations parameters (path)
    use log_mod     ! Acces to logging procedures

    implicit none

    public
    !
    ! Define scale specific parameters
    real(kind=rkd)    :: lStar           ! Star formation scale (~ 0.1pc)
    real(kind=rkd)    :: muStar          ! Specific gas surface density at the star forming scale [CU]
    real(kind=rkd)    :: sigmaStar       ! Specific gas velocity dispersion at the star forming scale [CU]
    real(kind=rkd)    :: mu_slope        ! Slope index for the gas surface density scalling relation
    real(kind=rkd)    :: sigma_slope     ! Slope index for the velocity dispersion scalling relation
    real(kind=rkd)    :: ETRV            ! Energy transfert rate per volume unit
    !
    ! Define gas structuration history parameters
    integer(kind=ikd) :: nScales         ! Number of structuration scales
    real(kind=rkd)    :: stepFactor      ! Stepping factor of the gas structuration cascade
    !
    ! Define stellar population parameters
    character(MAXPATHSIZE)  :: IMF       ! Initial Mass Function Reference
    !
    ! Define feedback parameters
    real(kind=rkd)     :: gasDisruptionEfficiency  ! Gas disruption efficiency [0, 1]
    real(kind=rkd)     :: Vwind                    ! Velocity of gas ejected from the gsh

contains

    !
    ! SUBROUTINES
    !

    ! **********************************
    subroutine model_init()

        ! Initialize the GAS model configuration

        implicit none

        ! Read the main parameter file
        call model_read_configuration_file()

    end subroutine model_init

    ! **********************************
    subroutine model_read_configuration_file()
    
        ! Read the model configuration file containing
        ! model paramters and integration scheme configuration
    
        implicit none
    
        integer(kind=ikd)         :: loop
    
        character(MAXPATHSIZE)  :: filename
        character(MAXPATHSIZE)  :: message
        character(MAXPATHSIZE)  :: line
        character(MAXPATHSIZE)  :: name
        character(MAXPATHSIZE)  :: val
    
        ! Get user parameter file path
        call getarg(2, filename)
        !
        open(inputModelFileUnit, file=filename, status='old')
        do
            read(inputModelFileUnit, '(a)', end=2) line
            loop = scan(line, '=')
            ! '#' is the symbol for header comments in all input data files
            if ((loop == 0) .or. (line(1:1) == '#')) cycle
            !
            name = trim(adjustl(line(:loop-1)))
            val  = trim(adjustl(line(loop+1:)))
            ! '!' is the symbol for comments in all input data files
            loop = scan(val, '!')
            if (loop /= 0) val = trim(adjustl(val(:loop-1))) ! cut the comment part of the readed line
            !
            ! Run through the list of all keyword defined in the file
            ! and treat each one
            select case (trim(name))
                !
                ! Integration Scheme
                case ('solver')
                    read(val, '(a)') solver
                    call solver_init(solver)
                case ('dt_max', 'dt')
                    read(val, *) solver_dt
                !
                ! Scale
                case ('lStar')
                    read(val, *) lStar
                case ('sigmaStar')
                    read(val, *) sigmaStar 
                case ('muStar')
                    read(val, *) muStar
                case ('mu_slope')
                    read(val, *) mu_slope
                !
                ! Gsh
                case ('nScales')
                    read(val, *) nScales
                case ('stepFactor')
                    read(val, *) stepFactor
                !
                ! Stellar populations
                case ('IMF')
                    read(val, '(a)') IMF
                !
                ! Feedback
                case ('disrupt_eff', 'disrupt_efficiency')
                    read(val, *) gasDisruptionEfficiency
                case ('Vwind', 'vWind', 'VWind')
                    read(val, *) Vwind
                    ! conversion in code unit
                    Vwind = Vwind * km_s2kpc_Gyr
                case default
                    write(message, '(a,a,a)') 'Model parameter ', trim(name), ' is unknown'
                    call log_message(message, logLevel=LOG_ERROR)
            end select
        end do
2       close(11)
        !
        write(message, *) 'Read model configuration file : ',trim(filename)//new_line('a')// &
                        '| INTEGRATOR'//new_line('a')// &
                        '|> intScheme : '//trim(solver)//new_line('a')// &
                        '| SCALE CONFIGURATION'//new_line('a')// &
                        '|> lStar     : ', lStar, new_line('a')// &
                        '|> sigmaStar : ', sigmaStar, new_line('a')// &
                        '|> muStar    : ', muStar, new_line('a')// &
                        '|> mu_slope  : ', mu_slope, new_line('a')// &
                        '| GSH CONFIGURATION'//new_line('a')// &
                        '|> lStar           : ', nScales, new_line('a')// &
                        '|> stepFactor      : ', stepFactor, new_line('a')// &
                        '| STELLAR POPULATION'//new_line('a')// &
                        '|> IMF             : ', trim(IMF), new_line('a')
        call log_message(message)
    end subroutine model_read_configuration_file

end module model_mod
