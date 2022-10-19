module log_mod

    !*****************************************************************************************************************! 
    ! OVERVIEW
    !
    ! This module defined log printing procedures
    !
    !*****************************************************************************************************************
    
    use parameters  ! Acces to global definitions and properties
    use PrDi_mod    ! Acces to multi processing distribution methods

    implicit none
   
    public

    ! Define some constants and parameters
    !
    ! Define the main log directories
    character(MAXPATHSIZE)     :: logPath
    !
    ! Log level
    integer(kind=4), parameter :: LOG_INFO = 1
    integer(kind=4), parameter :: LOG_WARNING = 2
    integer(kind=4), parameter :: LOG_ERROR = 3
    !
    ! Common file unit for log
    !   the rank of the current process is added to return
    !   the specific file unit of the current process
    integer(kind=4), parameter :: logProcessUnit = 2000

contains

    !
    ! SUBROUTINES
    !

    ! **********************************
    subroutine log_init(lPath)

        ! Initialize and open the log file of the current process

        implicit none

        integer(kind=4)                     :: u        ! log file unit

        character(MAXPATHSIZE), intent(in)  :: lPath    ! log directory
        character(MAXPATHSIZE)              :: filename ! log filename

        ! Set the log directory
        logPath = lPath

        ! Get file unit associated to the current process myRank
        u = log_get_process_file_unit()

        write(filename,'(a,a,i3.3,a)') trim(logPath), '/Process_', myRank, '.log'
        open(unit=u, file=filename, status='new')
        call log_message('Log process initialized')

    end subroutine log_init

    ! **********************************
    subroutine log_finalize()

        ! Finalize and close the log file of the current process

        implicit none

        integer(kind=4) :: u  ! log file unit

        ! Get file unit associated to the current process myRank
        u = log_get_process_file_unit()

        close(u)

    end subroutine log_finalize

    ! **********************************
    subroutine log_message(message, logLevel, componentRank, &
                           paramNames, realParams, intParams, calledBy)

        ! Print a standard log message

        implicit none

        integer(kind=4)                          :: i                ! loop over physical parameters that should be printed
        integer(kind=4)                          :: u                ! file unit
        integer(kind=4), intent(inout), optional :: componentRank    ! component indentation for log printing
        integer(kind=4)                          :: compRank
        integer(kind=4), intent(in), optional    :: logLevel         ! log Level (info, warning, error)
        integer(kind=4), intent(in), optional    :: intParams(:)     ! values of physical parameters (integer format)

        character(*), intent(in)                 :: message          ! the message
        character(*), intent(in), optional       :: calledBy         ! trace back of the error message
        character(30)                            :: space            ! indentation
        character(18)                            :: currentTime      ! current time
        character(12)                            :: lLevel           ! log level
        character(25), intent(in), optional      :: paramNames(:)    ! name of physical parameters

        real(kind=8), intent(in), optional       :: realParams(:)    ! values of physical parameters (real format)
        real(kind=8)                             :: time             ! current time
        real(kind=8)                             :: hour, min, sec   !

        call cpu_time(time)                   ! current time
        hour = floor(time/60./60.)            ! number of hours
        min  = floor((time-hour*60.*60.)/60.) ! number of min
        sec  = time-hour*60.*60-min*60.       ! number of sec
        write(currentTime,'(i2.2,a,i2.2,a,f5.2,a)') int(hour), 'h ', int(min), 'min ', sec, 'sec'
        
        if (present(componentRank)) then
            compRank = componentRank
        else
            compRank = 1
        end if
        space(compRank:) = ' |-> '
        
        ! Get file unit associated to the current process myRank
        u = log_get_process_file_unit()

        lLevel = 'INFO | '
        if (present(logLevel)) then
            if (logLevel .eq. LOG_ERROR) then
                lLevel = 'ERROR | '
            else if (logLevel .eq. LOG_WARNING) then
                lLevel = 'WARNING | '
            end if
        end if

        write(u, *) trim(currentTime), trim(space), trim(lLevel), trim(message)

        if (present(paramNames)) then
            do i = 1, size(paramNames)
                if (present(realParams)) then
                    write(u, *) '> ', paramNames(i), ' = ', realParams(i)
                else if (present(intParams)) then
                    write(u, *) '> ', paramNames(i), ' = ', intParams(i)
                end if
            end do
        end if

        if (present(calledBy)) then
            write(u,'(a,a)') 'Called by: ', calledBy
        end if

        if (present(logLevel)) then
            if (logLevel .eq. LOG_ERROR) then
                stop
            end if
        end if

        return
    end subroutine log_message

    !
    ! FUNCTIONS
    !

    ! **********************************
    function log_get_process_file_unit()
  
        ! Return the log file unit associated to the current process

        implicit none 
        
        integer(kind=4)  :: log_get_process_file_unit
        
        log_get_process_file_unit = logProcessUnit + myRank
        
        return
    end function log_get_process_file_unit
    
    ! **********************************
    function log_calledBy(calledBy, from)

        ! Format the calledBy message storing calling sequence

        implicit none

        character(MAXPATHSIZE), intent(in), optional  :: calledBy
        character(MAXPATHSIZE), intent(in), optional  :: from
        character(MAXPATHSIZE)                        :: log_calledBy

        write(log_calledBy, *) ''
        if (present(calledBy)) then
            write(log_calledBy, *) trim(calledBy)
        end if
        if (present(from)) then
            write(log_calledBy, *) trim(log_calledBy), new_line('a')//'From ', trim(from)
        end if

    end function log_calledBy

end module log_mod
