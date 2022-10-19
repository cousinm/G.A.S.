module config_mod

    !*****************************************************************************************************************
    !
    !  OVERVIEW
    !
    !  This module allows to define and set all environment parameters
    !  like path directories by loading the main GAS parameters files
    !
    !*****************************************************************************************************************

    use parameters  ! Acces to global defintions and properties
    use log_mod     ! Acces to logging procedures

    implicit none

    public

    ! File units for input files
    integer(kind=4), parameter   :: inputParameterFileUnit = 110  ! Used to read input parameter file

    ! Main directories path
    character(MAXPATHSIZE)       :: librariesPath
    character(MAXPATHSIZE)       :: stellarPopulationPath
    character(MAXPATHSIZE)       :: treesPath
    character(MAXPATHSIZE)       :: outputPath

contains

    !
    ! SUBROUTINES
    !

    ! **********************************
    subroutine config_init()

        ! Initialize the GAS environment configuration

        implicit none

        ! Read the main parameter file
        call config_read_main_configuration_file()

    end subroutine config_init

    ! **********************************
    subroutine config_read_main_configuration_file()
    
        ! Read the main user simulation configuration file containing
        ! input, output pathes all model and tree cleaning parameters
    
        implicit none
    
        integer(kind=4)         :: loop
    
        character(MAXPATHSIZE)  :: filename
        character(MAXPATHSIZE)  :: message
        character(MAXPATHSIZE)  :: line
        character(MAXPATHSIZE)  :: name
        character(MAXPATHSIZE)  :: val
        character(MAXPATHSIZE)  :: lPath
    
        ! Get user parameter file path
        call getarg(1, filename)
        !
        open(inputParameterFileUnit, file=filename, status='old')
        do
            read(inputParameterFileUnit, '(a)', end=2) line
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
                case ('logPath')
                    read(val, '(a)') lPath          ! tree files path
                    ! Init log file using lPath
                    call log_init(lPath)
                !
                case ('librariesPath')
                    read(val, '(a)') librariesPath  ! libraries input path
                !
                case ('treesPath')
                    read(val, '(a)') treesPath      ! tree files path
                !
                case ('outputPath')
                    read(val, '(a)') outputPath     ! Main output directory path
            end select
        end do
        !
2       close(11)
        !
        write(message, '(a,a)') 'Read configuration file : ', trim(filename)//new_line('a')// &
                                '| PATH'//new_line('a')// &
                                '|> logPath         : '//trim(lPath)//new_line('a')// &
                                '|> librariesPath   : '//trim(librariesPath)//new_line('a')// &
                                '|> treesPath       : '//trim(treesPath)//new_line('a')// &
                                '|> outputPath      : '//trim(outputPath)//new_line('a')
        call log_message(message)
    end subroutine config_read_main_configuration_file

    ! **********************************
    subroutine config_read_cosmology
    
      ! Read cosmological parameters associated to the DM simulation

      implicit none

  end subroutine config_read_cosmology
    
end module config_mod
