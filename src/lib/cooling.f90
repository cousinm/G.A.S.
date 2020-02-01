module cooling

  use gas     ! Contains gas structure and gas exploitation functions

  public

  !*****************************************************************************************************************
  ! 
  ! OVERVIEW
  !
  ! SUBROUTINES IN THIS MODULE
  !
  ! cooling_read_cooling_efficiency  : Read the cooling efficiency table
  !
  ! cooling_read_thermal_instability : Read the thermal instability efficiency table
  ! 
  !
  ! cooling_function_finalize        : Deallocate all data array used by the cooling function module
  !
  ! FUNCTION IN THIS MODULE 
  !
  ! Lambda(T,Z)                      : Return the cooling efficiency at a given temperature "T" and a given metallicity "Z" 
  !
  ! dlnLambda_dlnT(T,Z)              : Return the thermal instability parameter  at a given temperature "T" and a given metallicity "Z"                  
  !
  ! t_cool                           : Return the cooling timescale associated to a given gas component
  !
  ! MFF                              : Return the structured/fragmented gas fraction associated to a given t_cool_eff / cooling time-scale
  !
  ! dMFF                             : Return the fraction (in mass) of diffuse/unstructured gas that can effectivelly condensate 
  !                                        as a function of the effective cooling time
  !
  !*****************************************************************************************************************

  integer(kind=4)           :: cooling_nTempBins                   ! Number of temperature bins
  integer(kind=4)           :: cooling_nMetBins                    ! Number of metallicity bins
  
  real(kind=8),allocatable  :: cooling_TempBins(:)                 ! Temperature bins
  real(kind=8),allocatable  :: cooling_MetBins(:)                  ! Metallicity bins 
  real(kind=8),allocatable  :: log_cooling_efficiency_table(:,:)   ! log of the cooling_efficiency (erg.cm^3/s) (first dim = Temp, second = Met)
  real(kind=8),allocatable  :: thermal_instabillity_table(:,:)     ! termal instabillity parameter: 2. - dln Lambda / dln T log (first dim = Temp, second = Met)
  real(kind=8)              :: dlTemp, dlMet                       ! Temp and Met log bin size
  real(kind=8)              :: lTempmin, lMetmin                   ! minimum lTemp and lMet bin (use for localisation process in tables)

  contains

  !*****************************************************************************************************************
  
  subroutine cooling_read_cooling_efficiency
  
    ! LOAD COOLING EFFICIENCIES
    ! Pre-tabulated temperature and metallicity-dependent values are load from 'cooling_efficiency.in'

    implicit none

    integer(kind=4)          :: i,j             ! loop indexes under temperature and metallicity  

    character(MAXPATHSIZE)   :: filename,line   ! line is a tmp character variable used to read the file
    character(MAXPATHSIZE)   :: message         ! a message  
    
    call IO_print_message('cooling_read_cooling_efficiency')
    
    if (main_process .or. physical_process) then
      ! 
      ! Build the input filename
      write(filename,'(a,a)') trim(input_path), '/cooling_efficiency.in'
      write(message,'(a,a)') 'Load data from : ', 'cooling_efficiency.in'
      call IO_print_message(message)
      !
      ! open the file, using cooling_unit
      open(unit = cooling_unit,file = trim(filename),status = 'old')
      do
        read(cooling_unit, '(a)', end = 2) line
        if (trim(line) .eq. 'START') then
          !
          ! Read
          read(cooling_unit,*) cooling_nTempBins, cooling_nMetBins ! read numbers of bins in temperature and in metallicity
          read(cooling_unit,*) dlTemp, dlMet                       ! log bin for temperature and metalicity
          !
          if (physical_process) then
            !
            ! allocate arrays
            allocate(cooling_TempBins(cooling_nTempBins))                              ! Temperature in log scale
            allocate(cooling_MetBins(cooling_nMetBins))                                ! Metalicity in log scale
            allocate(log_cooling_efficiency_table(cooling_nTempBins,cooling_nMetBins)) ! cooling efficiency table (T,Z) in log scale
            !
            read(cooling_unit,*) cooling_TempBins ! read Temperature log bins
            read(cooling_unit,*) cooling_MetBins  ! read metallicity log bins
            !
            ! define minimum bin values
            lTempmin = minval(cooling_TempBins)
            lMetmin  = minval(cooling_MetBins)
            !
            read(cooling_unit,*) ! read a comment line #
            !
            do i = 1, cooling_nTempBins
              read(cooling_unit,*) (log_cooling_efficiency_table(i,j),j=1,cooling_nMetBins)
            end do  
          end if
          !
          exit  ! quit do loop
        end if
        if (line(1:1) .eq. '#') then
          !
          cycle ! header or something like this (skip)
        else
          !
          call IO_print_error_message('Impossible to read the line',only_rank=rank,called_by='cooling_read_cooling_efficiency') 
          write(*,*) line
          stop ! stop the program
        end if  
      end do
2     close(cooling_unit)
    end if

    return
  end subroutine cooling_read_cooling_efficiency

  !*****************************************************************************************************************
  
    subroutine cooling_read_thermal_instability
    
    ! LOAD THERMAL INSTABILITY EFFICIENCIES
    ! Pre-tabulated temperature and metallicity-dependent values are load from 'termal_instability.in'

    implicit none

    integer(kind=4)          :: i,j             ! loop indexes under temperature and metallicity  

    character(MAXPATHSIZE)   :: filename,line   ! line is a tmp character variable used to read the file
    character(MAXPATHSIZE)   :: message         ! a message  
    
    call IO_print_message('cooling_read_thermal_instability')
    
    if (main_process .or. physical_process) then
      !
      ! Build the input filename
      write(filename,'(a,a)') trim(input_path), '/thermal_instability.in'
      write(message,'(a,a)') 'Load data from : ', 'thermal_instability.in'
      call IO_print_message(message)
      !
      ! open the file, using cooling_unit
      open(unit = cooling_unit,file = trim(filename),status = 'old')
      !
      do
        read(cooling_unit, '(a)', end = 2) line
        if (trim(line) .eq. 'START') then
          !
          ! Read
          read(cooling_unit,*) ! cooling_nTempBins, cooling_nMetBins ! already loaded 
          read(cooling_unit,*) ! dlTemp, dlMet                       ! already loaded 
          !
          if (physical_process) then
              !
              ! allocate array
              allocate(thermal_instabillity_table(cooling_nTempBins,cooling_nMetBins))
              !
              read(cooling_unit,*) ! cooling_TempBins ! already loaded 
              read(cooling_unit,*) ! cooling_MetBins  ! already loaded 
              !
              read(cooling_unit,*) ! read a comment line #
              !
              do i = 1, cooling_nTempBins
                read(cooling_unit,*) (thermal_instabillity_table(i,j),j=1,cooling_nMetBins)
              end do  
          end if
          exit  ! quit do loop
        end if
        if (line(1:1) .eq. '#') then
          !
          cycle ! header or something like this (skip)
        else
          !
          call IO_print_error_message('Impossible to read the line',only_rank=rank,called_by='cooling_read_thermal_instability') 
          stop ! stop the program
        end if  
      end do
2     close(cooling_unit)
    end if

    return
  end subroutine cooling_read_thermal_instability

  !*****************************************************************************************************************

  subroutine cooling_function_finalize
    
    ! DEALLOCATE ALL DATA ARRAY USED BY THE COOLING FUNCTION LIBRARY
    
    implicit none
    
    if (allocated(cooling_TempBins))             deallocate(cooling_TempBins)              ! Temperature bins
    if (allocated(cooling_MetBins))              deallocate(cooling_MetBins)               ! Metallicity bins 
    if (allocated(log_cooling_efficiency_table)) deallocate(log_cooling_efficiency_table)  ! first dim = Temp, second = Met
    if (allocated(thermal_instabillity_table))   deallocate(thermal_instabillity_table)    ! first dim = Temp, second = Met
    
    return
  end subroutine cooling_function_finalize

  !*****************************************************************************************************************
  !
  ! FUNCTIONS
  !
  !*****************************************************************************************************************

  function Lambda(T,Z)  

    ! RETURN THE COOLING EFFICIENCY FOR A GIVEN TEMPERATURE AND METALLICITY
    ! The efficiency is given in Joule/s * m3 (SI)

    implicit none

    real(kind=8),intent(in)   :: T, Z        ! temperature and metallicity given in input
    real(kind=8)              :: Lambda      ! The cooling efficiency
    real(kind=8)              :: logT, logZ  ! log values for the temperature and for the metallicity

    ! Temperature 
    logT = 4.0
    if (T .gt. 0.d0) logT = log10(T)
    !
    ! Metalicity
    logZ = -3.d0
    if (Z .gt. 0.d0) logZ = log10(Z)
    
    ! Load the log(Lambda) 
    Lambda = locate2D(logT,logZ,log_cooling_efficiency_table, &  
                            cooling_nTempBins,cooling_nMetBins, &
                            lTempmin,lMetmin, &
                            dlTemp,dlMet)
    
    Lambda = (1.d1)**Lambda  ! Here the value is given in erg.cm^3/s
    
    ! conversion in Internationnal unit Joules/s.m^3
    ! 1 erg = 10-7 Joules
    Lambda = Lambda*1.d-7
    ! 1 m = 1.d-2 cm
    Lambda = Lambda*(1.d-2)**3.   ! Here lambda is given in J/s m^3 = kg m^5 / s^3
   
    return
  end function Lambda
 
  !*****************************************************************************************************************
  
  function dlnLambda_dlnT(T,Z)  

    ! RETURN THE THERMAL INSTABILITY EFFICIENCY FOR A GIVEN TEMPERATURE AND METALLICITY 

    implicit none

    real(kind=8),intent(in)   :: T, Z           ! temperature and metallicity given in input
    real(kind=8)              :: dlnLambda_dlnT ! the thermal instability parameter
    real(kind=8)              :: logT, logZ     ! log values for the temperature and for the metallicity

    ! Temperature 
    logT = 4.d0
    if (T .gt. 0.d0) logT = log10(T)
    !
    ! Metalicity
    logZ = -3.0
    if (Z .gt. 0.d0) logZ = log10(Z)
    
    ! Load the thermal instability parameter
    dlnLambda_dlnT  = locate2D(logT,logZ,thermal_instabillity_table, &  
                            cooling_nTempBins,cooling_nMetBins, &
                            lTempmin,lMetmin, &
                            dlTemp,dlMet)
    return
  end function dlnLambda_dlnT
  
  !*****************************************************************************************************************  
  
  function t_cool(T_g,rho_g,Z_g)

    ! RETURN THE COOLING TIME OF A GAS WITH A TEMPERATURE T_g, A DENSITY rho_g AND A METALLICITY Z_g 
        
     implicit none

     real(kind=8),intent(in) :: T_g,Z_g,rho_g  ! cooling process properties (temperature, metalicity and density)
     real(kind=8)            :: t_cool         ! the cooling timescale

     ! Cornuault their Eq 4. 
     t_cool = 4.3d0*(5.d0/2.d0)*mu*kb*(T_g)                      ! in Joules*kg  = kg^2*m^2/s^2 
     t_cool = t_cool/(rho_g*mass_code_unit_in_kg/(kpc_in_m**3.)) ! in Joules*m^3 = kg*m^5/s^2 / conversion of the gas density in internationnal unit (kg/m^3)   
     t_cool = t_cool/Lambda(T_g,Z_g)                             ! in sec  (because Lambda is in Joules/s.m^3 = kg*m^5/s^3)
     t_cool = t_cool/Gyr_in_s                                    ! in Gyr

     return
  end function t_cool
 
  !*****************************************************************************************************************
  
  function MFF(x)
  
    ! RETURN THE STRUCTURED/FRAGMENTED GAS FRACTION ASSOCIATED TO A GIVEN t_cool_eff / cooling time-scale
    
    real(kind=8),intent(in)   :: x
    real(kind=8)              :: x_
    real(kind=8),parameter    :: x_t = 5.51885003749d-1
    real(kind=8),parameter    :: s   = 1.33172792468d-1
    real(kind=8)              :: MFF
    
    if (x .gt. 1.d1) then
        MFF = 1.d0
        return
    end if
    
    if (x .lt. 1.d-2) then
        MFF = 0.d0
        return
    end if
    
    x_ = log10(x/x_t)/sqrt(s)
    
    MFF = 5.d-1*(1.0d0+erf(x_))  
    
    return
  end function MFF
   
  !*****************************************************************************************************************
  
  function dMFF_dt(x)
  
    ! RETURN THE FRACTION (in mass) OF UNSTRUCTURED (diffuse) GAS THAT CAN CONDENSATE 
    ! AS FUNCTION OF THE EFFECTIVE COOLING TIME
    
    implicit none
    
    real(kind=8),intent(in)   :: x
    real(kind=8)              :: x_
    real(kind=8),parameter    :: amp = 1.44887662202d0
    real(kind=8),parameter    :: x_t = 3.87234590598d-1
    real(kind=8),parameter    :: s   = 2.58508609094d-1
    real(kind=8)              :: dMFF_dt
    
    ! t_ = effective cooling time (clock) / cooling timescale (10⁴, Z)
  
    x_ = log10(x/x_t)**2./(2.d0*s**2.)
    
    dMFF_dt = amp*exp(-x_)
    
    return
  end function dMFF_dt
  
  !*****************************************************************************************************************

end module cooling
