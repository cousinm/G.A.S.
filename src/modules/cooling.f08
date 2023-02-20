module cooling_mod

    !*****************************************************************************************************************
    !
    !  OVERVIEW
    !
    ! The cooling module defines all properties and procedures
    ! associated to gas cooling (efficiency and thermal instabilities aspect)
    ! Those parameters are defined according to metalicity and temperature
    !
    !*****************************************************************************************************************
    
    use parameters  ! Acces to global defintions and properties
    use tools_mod   ! Acces to tool methods (locate2D)
    use config_mod  ! Acces to configurations parameters (path)
    use gas_mod     ! Acces to gas object and methods

    implicit none

    public

    ! Cooling class definition

    type cooling
        real(kind=rkd)    :: clock  ! Gas clock
        type(gas)         :: gas    ! The gas stored
    contains
        procedure :: t_cool => cooling_t_cool
    end type cooling

    !
    ! Define specific constants and parameters
    integer(kind=ikd)           :: cooling_nTempBins                   ! Number of temperature bins
    integer(kind=ikd)           :: cooling_nMetBins                    ! Number of metallicity bins
  
    real(kind=rkd),allocatable  :: cooling_tempBins(:)                 ! List of temperature bins
    real(kind=rkd),allocatable  :: cooling_metBins(:)                  ! List of metallicity bins 
    real(kind=rkd),allocatable  :: log_cooling_efficiency_table(:, :)  ! log of the cooling_efficiency (erg/s.cm^3) (first dim = Temp, second = Met)
    real(kind=rkd),allocatable  :: thermal_instabillity_table(:, :)    ! termal instabillity parameter: 2. - dln Lambda / dln T log (first dim = Temp, second = Met)
    real(kind=rkd)              :: dlTemp, dlMet                       ! Temp and Met elementary log step
    real(kind=rkd)              :: lTempMin, lMetMin                   ! minimum lTemp and lMet bin (use for localisation process in tables)
    
contains

    !
    ! SUBROUTINES
    !
    ! **********************************

    subroutine cooling_init()

        ! Initialize the cooling module
        !   Load cooling efficiency and thermal instability parameter tables

        implicit none

        call cooling_read_cooling_efficiency()
        call cooling_read_thermal_instability()

    end subroutine cooling_init

    ! **********************************
    subroutine cooling_finalize
    
        ! Deallocate all cooling tables
        
        implicit none
        
        if (allocated(cooling_tempBins))             deallocate(cooling_TempBins)              ! Temperature bins
        if (allocated(cooling_metBins))              deallocate(cooling_MetBins)               ! Metallicity bins 
        if (allocated(log_cooling_efficiency_table)) deallocate(log_cooling_efficiency_table)  ! first dim = Temp, second = Met
        if (allocated(thermal_instabillity_table))   deallocate(thermal_instabillity_table)    ! first dim = Temp, second = Met

    end subroutine cooling_finalize

    ! **********************************
    subroutine cooling_read_cooling_efficiency

        ! Load cooling efficiency data
        ! These data are pre-tabultaed accorind to metalicities and temperatures

        implicit none

        integer(kind=ikd)        :: i,j             ! Loop indexes under temperature and metallicity  

        character(MAXPATHSIZE)   :: filename
        character(MAXPATHSIZE)   :: line            ! A red line
        character(MAXPATHSIZE)   :: message         ! A message

        ! Build the input filename
        write(filename,'(a, a)') trim(librariesPath), '/cooling_efficiency.in'
        write(message,'(a, a)') 'Load data from: ', trim(filename)
        call log_message(message)

        ! Open the file
        open(unit=cooling_unit, file=trim(filename), status='old')

        ! Read and load data
        do
            read(cooling_unit, '(a)', end=2) line
            if (trim(line) .eq. '----') then
                !
                ! Read
                read(cooling_unit, *) cooling_nTempBins, cooling_nMetBins ! Read numbers of bins for temperatures and metallicities
                read(cooling_unit, *) dlTemp, dlMet                       ! elementary temperature and metalicity log step
                !
                ! allocate arrays
                allocate(cooling_TempBins(cooling_nTempBins))                               ! Temperatures scale (in log)
                allocate(cooling_MetBins(cooling_nMetBins))                                 ! Metalicity scale (in log)
                allocate(log_cooling_efficiency_table(cooling_nTempBins, cooling_nMetBins)) ! Cooling efficiency table (T, Z) (in log)
                !
                ! Load scales
                read(cooling_unit, *) cooling_TempBins ! Read Temperatures
                read(cooling_unit, *) cooling_MetBins  ! Read metallicities
                !
                ! Define minimum bin values
                lTempMin = minval(cooling_tempBins)
                lMetMin  = minval(cooling_metBins)
                !
                ! Read the lastest comment line #
                read(cooling_unit, *)
                !
                ! Load data
                do i = 1, cooling_nTempBins
                    read(cooling_unit, *) (log_cooling_efficiency_table(i, j), j=1, cooling_nMetBins) ! [erg/s.cm3]
                end do
                !
                exit  ! quit do loop
            end if
            if (line(1:1) .eq. '#') then
                ! Header or something like this (skip)
                cycle 
            else
                write(message, '(a, a)') 'Impossible to read a line in ', trim(filename)
                call log_message(message, logLevel=LOG_ERROR, &
                                calledBy='cooling_read_cooling_efficiency')
            end if  
        end do
2       close(cooling_unit)

    end subroutine cooling_read_cooling_efficiency

    ! **********************************
    subroutine cooling_read_thermal_instability
    
        ! Load thermal instability parameter
        ! These data are pre-tabultaed accorind to metalicities and temperatures

        implicit none

        integer(kind=ikd)        :: i,j             ! Loop indexes under temperature and metallicity  

        character(MAXPATHSIZE)   :: filename,line   ! A red line
        character(MAXPATHSIZE)   :: message         ! A message

        ! Build the input filename
        write(filename,'(a, a)') trim(librariesPath), '/thermal_instability.in'
        write(message,'(a, a)') 'Load data from: ', trim(filename)
        call log_message(message)

        ! Open the file
        open(unit=cooling_unit, file=trim(filename), status='old')

        ! Read and load data
        do
            read(cooling_unit, '(a)', end=2) line
            if (trim(line) .eq. '----') then
                !
                ! Read
                read(cooling_unit, *) ! Cooling_nTempBins, cooling_nMetBins ! already loaded
                read(cooling_unit, *) ! dlTemp, dlMet                       ! already loaded
                !
                ! allocate array
                allocate(thermal_instabillity_table(cooling_nTempBins, cooling_nMetBins))
                !
                read(cooling_unit, *) ! cooling_TempBins ! already loaded
                read(cooling_unit, *) ! cooling_MetBins  ! already loaded
                !
                ! Read the lastest comment line #
                read(cooling_unit, *)
                !
                ! Load data
                do i = 1, cooling_nTempBins
                    read(cooling_unit, *) (thermal_instabillity_table(i,j), j=1, cooling_nMetBins)
                end do  
                exit  ! Quit do loop
            end if
            if (line(1:1) .eq. '#') then
                !
                cycle ! Header or something like this (skip)
            else
                write(message, '(a, a)') 'Impossible to read a line in ', trim(filename)
                call log_message(message, logLevel=LOG_ERROR, &
                                calledBy='cooling_read_thermal_instability')
            end if  
        end do
2       close(cooling_unit)

    end subroutine cooling_read_thermal_instability

    !
    ! FUNCTIONS
    !
    ! **********************************

    function Lambda(T, Z) result(val)

        ! Return the cooling efficiency for a given temperature and metalicity
        ! The efficiency is given in J/s.m3 (SI)

        implicit none

        character(MAXPATHSIZE)      :: message         ! A message

        real(kind=rkd),intent(in)   :: T, Z        ! temperature and metallicity given in input
        real(kind=rkd)              :: val         ! The cooling efficiency
        real(kind=rkd)              :: logT, logZ  ! log values for the temperature and for the metallicity

        ! Temperature, init to minimal value
        logT = lTempMin
        if (T .gt. 0.d0) then
            logT = log10(T)
        else
            write(message, *) 'Temperature is not valid: ', T
            call log_message(message, logLevel=LOG_ERROR, &
                             calledBy='Lambda')
        end if
        !
        ! Metalicity, init to minimal value
        logZ = lMetMin
        if (Z .gt. 0.d0) then
            logZ = log10(Z)
        else
            write(message, *) 'Metalicity is not valid: ', Z
            call log_message(message, logLevel=LOG_ERROR, &
                             calledBy='Lambda')
        end if
        
        ! Get current cooling efficiency
        val = locate2D(logT, logZ, log_cooling_efficiency_table, &
                       cooling_nTempBins, cooling_nMetBins, &
                       lTempmin, lMetmin, &
                       dlTemp, dlMet)
        
        val = (1.d1)**val  ! Here the value is given in erg.cm^3/s
        
        ! conversion in Internationnal unit J/s.m^3
        ! 1 erg = 10-7 Joules
        val = val*1.d-7
        ! 1 m = 1.d-2 cm
        val = val*(1.d-2)**3.   ! Here val is given in J/s m^3 = kg m^5 / s^3

    end function Lambda

    ! **********************************
  
    function dlnLambda_dlnT(T,Z) result(val)

        ! Return the thermal instability parameter for a given tempertaure and metalicity

        implicit none

        character(MAXPATHSIZE)      :: message     ! A message

        real(kind=rkd),intent(in)   :: T, Z        ! temperature and metallicity given in input
        real(kind=rkd)              :: val         ! The thermal instability parameter
        real(kind=rkd)              :: logT, logZ  ! log values for the temperature and for the metallicity

        ! Temperature, init to minimal value
        logT = lTempMin
        if (T .gt. 0.d0) then
            logT = log10(T)
        else
            write(message, *) 'Temperature is not valid: ', T
            call log_message(message, logLevel=LOG_ERROR, &
                             calledBy='dlnLambda_dlnT')
        end if
        !
        ! Metalicity, init to minimal value
        logZ = lMetMin
        if (Z .gt. 0.d0) then
            logZ = log10(Z)
        else
            write(message, *) 'Metalicity is not valid: ', Z
            call log_message(message, logLevel=LOG_ERROR, &
                             calledBy='dlnLambda_dlnT')
        end if
        
        ! Load the thermal instability parameter
        val = locate2D(logT, logZ, thermal_instabillity_table, &  
                       cooling_nTempBins, cooling_nMetBins, &
                       lTempmin, lMetmin, &
                       dlTemp,dlMet)

    end function dlnLambda_dlnT

    ! **********************************

    function cooling_t_cool(this, T_, rho) result(t_cool)

        ! Return the effective cooling time to reach T_target for a gas with
        !   - tempertaure T
        !   - metalicity: Z [metal mass fraction]
        !   - molecular mass: mu
        !   - average density: rho

        implicit none
    
        real(kind=rkd), intent(in) :: T_        ! Target temperature
        real(kind=rkd), intent(in) :: rho       ! Gas density
        real(kind=rkd)             :: T, Z, mu
        real(kind=rkd)             :: ntot, ne
        real(kind=rkd)             :: t_cool    ! The effective cooling timescale

        class(cooling), intent(in) :: this      ! The cooling gas component

        ! Get gas properties
        T = this%gas%temperature()
        Z = this%gas%metalicity()
        mu = this%gas%molecular_mass()
        !
        ntot = rho * MassCU2kg / LenCU2m**3.                    ! Conversion to SI (kg/m^3)
        ntot = ntot / (real(1.4, kind=rkd) * mu * mp)           ! ntot = rho/(1.4*mu*mp)
        !
        ne = real(1.2, kind=rkd) * ntot
        !
        ! Compute effective cooling time
        t_cool = real(3./2., kind=rkd) * ntot * kb * (T - T_)   ! in J*kg = kg^2*m^2/s^2                                                          
        t_cool = t_cool / (ne**2. * Lambda(T, Z))               ! in second (Lambda is given in J/s.m^3 = kg*m^5/s^3)
        t_cool = t_cool * s2TimeCU                              ! in code unit (Gyr)

    end function cooling_t_cool

end module cooling_mod