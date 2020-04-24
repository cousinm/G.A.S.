module stellar_population_library

  use gas    ! Contains gas structure and gas exploitation functions 

  public
  
  !*****************************************************************************************************************
  ! 
  ! OVERVIEW
  !
  ! Here are define all functions or subroutines dedicated to stellar population evolution
  ! In the header of this module are defined some global variables as the mass_loss_rate tables containing 
  ! ejecta_rates linked to a stellar population with a given age and metallicity.
  !
  ! SUBROUTINE IN THIS MODULE
  !
  ! stellar_population_read_mass_loss_rates      : Allows to read stellar ejecta rates for 
  !   called by : main program                         a simple stellar population (SSP) of a given age and metallicity
  !  
  ! stellar_population_read_SN_rates             : Allows to read global SN(SNIa+SNII) evenement rates for 
  !                                                    a simple stellar population (SSP) of a given age and metallicity
  !
  ! stellar_population_read_stellar_spectra      : Allows to read stellar spectra associated to SSPs
  !
  ! stellar_population_finalize                  : Deallocate all data array used by a stellar population
  !
  !*****************************************************************************************************************

  ! DEFINITION OF GLOBAL VARIABLES LINKED TO THE STELLAR POPULATION EVOLUTION ************

  integer(kind=4)            :: nAgeBins                         ! nb of age steps in the stellar evolution model
  integer(kind=4)            :: nWaves                           ! number wavelenghts used in stellar spectra 
  
  real(kind=8)               :: StellarTimeStep                  ! minimum age step separation
  real(kind=8),allocatable   :: AgeBins(:)                       ! stellar ages corresponding to age step
  real(kind=8),allocatable   :: SN_rates(:,:)                    ! the global SN (SNIa + SNII) event rates [nb/Gyr/Msun]
  real(kind=8),allocatable   :: LumBins(:,:)                     ! elementary luminosity in a given age/metallicity bin in Lsun/Msun
  real(kind=4),allocatable   :: Waves(:)                         ! wavelengh table [in micron 1e-6 m]
  real(kind=4),allocatable   :: sb_sed(:,:,:)                    ! starburst spectrum (nWaves, nAgeBins, nMetBins)
  real(kind=8)               :: ISRF_ref                         ! Reference factor for ISRF erg/s/cm^2 G0 = 1
  real(kind=8),allocatable   :: ISRFBins(:)                      ! scalling factor of ISRF
  
  type(gas_type),allocatable :: mass_loss_rates(:,:)             ! mass_loss_rates[a][Z] for given age 'a' and a metallicity 'Z'
                                                                 ! each element of mass_loss_rate is a gas object
                                                                 ! the total ejecta rate at the age 'a' and for the metallicity bin 'Z' is therefore given by
                                                                 ! - mass_loss_rates[a][Z]%mass
                                                                 ! the metals ejected rate of at the age 'a' and for the metallicity bin 'Z' is given by
                                                                 ! - mass_loss_rates[a][Z]%mZ
                                                                 ! and the ejecta rate associated to a element 'e' is given by 
                                                                 ! mass_loss_rates[a][Z]%elts[e]
  
  ! wavelengh reference indexes
  integer(kind=4)           :: i_FUV              ! index (in Waves) of the FUV filter central wavelenght (0.15mu m) 
                                                  ! setted in dust_read_abs_sca_properties
  integer(kind=4)           :: i_V, i_B           ! index (in Waves) of the Visible filter central wavelenght (0.55 mu m) 
                                                  !                     the Blue filter central wavelenght    (0.43 mu m) 
                                                  ! setted in dust_read_abs_sca_properties
  integer(kind=4)           :: i_8mic, i_1000mic  ! index (in Waves) of the IR range [8:1000] microns
                                                  ! setted in dust_read_dust_SEDs
  integer(kind=4)           :: i_14eV, i_6eV      ! index (in Waves) of the FUV (6ev < h*nu < 13.6eV) band
                                                  ! settled in dust_read_abs_sca_properties
                                                  
  ! define limit of the young stellar population
  real(kind=8),parameter     :: young_stars_MaxAge     =  5.d-2 ! [Gyr] = 5.e7 yr
                                                 
  ! model parameters associated to SN feedback processes
  real(kind=8),parameter     :: SN_kinetic_fraction    = 2.d0/3.d0 ! fraction of the instantaneous SN explosion power converted in kinetic power  
  real(kind=8)               :: SN_thermal_fraction    = 9.9d-1     ! fraction of the instantaneous non-kinetic SN explosion power converted in thermal power  
  !real(kind=8),parameter     :: p_sn                   = 3.d4                                             ! SN momentum (Msun km/s) produced by one SN 
                                                                                                          ! assuming that each SN produce 10M_sun of material moving a v = 3000 km.s-1 
  !real(kind=8),parameter     :: p_sn_code_unit         = p_sn/mass_code_unit_in_M_Sun/kpc_in_km*Gyr_in_s  ! Supernovae momentum in code unit : 10^11 Msun * kpc / Gyr
  real(kind=8),parameter     :: Esn_erg                = 1.d51                                            ! SN energy (erg)
  real(kind=8),parameter     :: Esn_J                  = Esn_erg*1.e-7                                    ! SN energy (Joule)
  real(kind=8),parameter     :: Esn_code_unit          = Esn_J*J_in_code_unit                             ! SN energy (code unit : 10^11 Msun * kpc^2 / Gyr^2) 
  
contains

  !*****************************************************************************************************************
  !
  ! SUBROUTINES
  !
  !*****************************************************************************************************************

  subroutine stellar_population_read_mass_loss_rates

    ! Read the mass loss rates for a stellar population of given age and metallicity
    ! the model used is linked to the initial mass function used

    implicit none
    
    integer(kind=4)          :: nMetBins_loc ! local nMetBins, 
                                             !   allow to compared the value read in stellar_population_mass_loss_rate file with 
                                             !   the global nMetBins parameter saved in the gas module
    integer(kind=4)          :: nElts_loc    ! local nElts
                                             !   allow to compared the value read in stellar_population_mass_loss_rate file with 
                                             !   the global nElts parameter saved in the gas module
    integer(kind=4)          :: i            ! loop index under stellar age
    integer(kind=4)          :: j            ! loop index under metalicity
    integer(kind=4)          :: e            ! loop index under elements
    
    character(MAXPATHSIZE)   :: filename
    character(MAXPATHSIZE)   :: message      ! a message  
    character(MAXPATHSIZE)   :: line

    call IO_print_message('stellar_population_read_mass_loss_rates')
     
    if (main_process .or. physical_process) then
       !
       ! Build the input filename, The initial mass function (IMF) is given in input
       write(filename,'(a,a,a,a)') trim(stellar_population_input_path), '/stellar_population_library_mass_loss_rates[BC03]', trim(IMF), '.in'
       write(message,'(a,a,a,a)') 'Load data from : ', 'stellar_population_library_mass_loss_rates[BC03]', trim(IMF), '.in'
       call IO_print_message(message)
       ! Open the library file
       open(unit = massloss_unit, file = filename, status = 'old')
       !
       do
          !
          read(massloss_unit, '(a)', end = 2) line
          if (trim(line) .eq. 'START') then
             !
             ! all lines before the START keywork are header lines: they contain informations about data 
             ! read nAgeBins  : the number of stellar age bin used in the stellar population model
             ! read nMetBins  : the number of metallicity bins used in the stellar population model
             ! read nElts     : the number of elements followed by the stellar population model (C,N,O,Fe)
             read(massloss_unit,*) nAgeBins, nMetBins_loc, nElts_loc
             ! test nMetBins and nElts already saved in the gas module
             if (nMetBins_loc .ne. nMetBins) then
                !
                call IO_print_error_message('nMetBins(stellar) != nMetBins(gas) ',called_by = 'stellar_population_read_mass_loss_rates')
                call IO_print_message('use',param_name=(/'nMetBins(stellar)        ','nMetBins(gas)            '/), &
                     int_param_val=(/nMetBins_loc,nMetBins/))
                stop ! stop the program
             endif
             !
             if (nElts_loc .ne. nElts) then
                !
                call IO_print_error_message('nElts(stellar) != nElts(gas) ',called_by = 'stellar_population_read_mass_loss_rates')
                call IO_print_message('use',param_name=(/'nElts(stellar)           ','nElts(gas)               '/), &
                     int_param_val=(/nElts_loc,nElts/))
                stop ! stop the program
             end if
             !
             call IO_print_message('use',param_name=(/'nAgeBins                 ','nMetBins                 ','nElts                    '/), &
                                         int_param_val=(/nAgeBins, nMetBins, nElts/))
             !
             read(massloss_unit,*) ! skip the line with Metallicity bins, these data are already saved in the gas module
             !
             ! allocate arrays  
             ! MetBins and Elt_names are already saved in the gas module
             ! stellar age (Gyr) corresponding to each step of the stellar model used
             allocate(AgeBins(nAgeBins))                   
             ! the mass_loss_rate array is composed of gas objects, 
             ! the ejecta rates associated to elements (H, He, C, N, O and Fe) are stored following the gas object structure
             allocate(mass_loss_rates(nAgeBins,nMetBins)) ! ejected mass fraction/Gyr for each age and metallicity bin
             !
             do i = 1, nAgeBins
                !
                do j = 1, nMetBins
                  !                      
                  call gas_void(mass_loss_rates(i,j)) ! initialize the mass_loss_rates array
                end do                 
                ! For each stellar age bin read : the age of the population, the mass loss rate and the yield of each elements take into acount
                read(massloss_unit,*) AgeBins(i), (mass_loss_rates(i,j)%mass,j=1,nMetBins), (mass_loss_rates(i,j)%mZ,j=1,nMetBins), ((mass_loss_rates(i,j)%elts(e),j=1,nMetBins),e=1,nElts)            
             end do
             !       
             ! Define StellarTimeStep
             ! StellarTimeStep is the minimal step in the stellar population evolution model
             ! It is therfore the maximum time step that the evolution model can used to evolve stellar population
             StellarTimeStep = minval(AgeBins(2:nAgeBins)-AgeBins(:nAgeBins-1)) 
             !
             if ((main_process) .and. (nbproc .gt. 1)) then
                !
                ! the main process doesn't not compute physical evolution
                if (allocated(mass_loss_rates)) deallocate(mass_loss_rates)
                call IO_print_message('use',param_name=(/'StellarTimeStep          '/),real_param_val=(/StellarTimeStep/))
             end if
             exit  ! quit do loop
          end if
          if (line(1:1) .eq. '#') then
            !                   
            cycle ! header or something like this (skip)
          else
             !
             call IO_print_error_message('Impossible to read the line',only_rank=rank,called_by='stellar_population_read_massloss') 
             stop ! stop the program
          end if
       end do
2      close(massloss_unit)
    end if

    return
  end subroutine stellar_population_read_mass_loss_rates

  !*****************************************************************************************************************

  subroutine stellar_population_read_SN_rates

    implicit none

    integer(kind=4)          :: nAgeBins_loc ! local nAgeBins,
                                             !   allow to compared the value read in stellar_population_SN_rates file with 
                                             !   the global nAgeBins parameter read previously in stellar_population_read_mass_loss_rates
    integer(kind=4)          :: nMetBins_loc ! local nMetBins, 
                                             !   allow to compared the value read in stellar_population_mass_loss_rate file with 
                                             !   the global nMetBins parameter saved in the gas module
    integer(kind=4)          :: i            ! loop index under stellar age
    integer(kind=4)          :: j            ! loop index under metalicity

    character(MAXPATHSIZE)   :: filename  
    character(MAXPATHSIZE)   :: message      ! a message  
    character(MAXPATHSIZE)   :: line

    real(kind=8)             :: agebin       ! a local variable used to read the stellar age

    call IO_print_message('stellar_population_read_SN_rates')

    if (main_process .or. physical_process) then
       !
       ! Build the input filename, The initial mass function (IMF) is given in input
       write(filename,'(a,a,a,a)') trim(stellar_population_input_path), '/stellar_population_library_SN_rates[BC03]', trim(IMF), '.in'
       write(message,'(a,a,a,a)') 'Load data from : ', 'stellar_population_library_SN_rates[BC03]', trim(IMF), '.in'
       call IO_print_message(message)
       ! Open the library file
       open(unit = snrate_unit, file = filename, status = 'old')
       
       do
          !
          read(snrate_unit, '(a)', end = 2) line
          if (trim(line) .eq. 'START') then
             !
             ! all lines before the START keywork are header lines: they contain informations about data 
             ! read nAgeBins  : the number of stellar age bin used in the stellar population model
             ! read nMetBins  : the number of metallicity bins used in the stellar population model
             read(snrate_unit,*) nAgeBins_loc, nMetBins_loc
             ! test nAgeBins and nMetBins already saved previously
             if (nAgeBins_loc .ne. nAgeBins) then
                !
                call IO_print_error_message('nAgeBins(SN) != nAgeBins(stellar) ',called_by = 'stellar_population_read_SN_rates')
                call IO_print_message('use',param_name=(/'nAgeBins(SN)             ','nAgeBins(stellar)        '/), &
                     int_param_val=(/nAgeBins_loc,nAgeBins/))
                stop ! stop the program
             end if
             !
             if (nMetBins_loc .ne. nMetBins) then
                !
                call IO_print_error_message('nMetBins(SN) != nMetBins(stellar) ',called_by = 'stellar_population_read_SN_rates')
                call IO_print_message('use',param_name=(/'nMetBins(SN)             ','nMetBins(stellar)        '/), &
                     int_param_val=(/nMetBins_loc,nMetBins/))
                stop ! stop the program
             endif
             read(snrate_unit,*) ! skip the line with Metallicity bins, these data are already saved in the gas module
             !
             if (physical_process) then
                !
                ! allocate arrays  
                ! AgeBins are already save in stellar_population_read_mass_loss_rates                
                ! SN_rate number of SN events per Gyr and per Msun [nb/Gyr/Msun]
                allocate(SN_rates(nAgeBins,nMetBins)) ! 
                SN_rates(:,:) = 0.d0                  ! initialize 
                do i = 1, nAgeBins
                    !
                    ! For each stellar age bin read : the age of the stellar population and the SN event rate
                    read(snrate_unit,*) agebin, (SN_rates(i,j),j=1,nMetBins)
                end do
                !
                ! convert in code unit
                SN_rates(:,:) = SN_rates(:,:)*mass_code_unit_in_M_Sun   ! nb/Gyr/(10^11Msun)
             end if
             exit  ! quit do loop
          end if
          if (line(1:1) .eq. '#') then
             !
             cycle ! header or something like this (skip)
          else
             !
             call IO_print_error_message('Impossible to read the line',only_rank=rank,called_by='stellar_population_read_SN_rates') 
             stop ! stop the program
          end if
       end do
2      close(massloss_unit)
    end if

    return
  end subroutine stellar_population_read_SN_rates

  !*****************************************************************************************************************
  
  subroutine stellar_population_read_stellar_spectra
    
    ! READ STELLAR SPECTRA AND ASSOCIATED INSTANTANEOUS LUMINOSITY
    !
    ! Physical processes have to read the library !!
    !  They need to have an acces to instantaneous luminosities of stellar populations
    !     to compute lum-weighted age or metallicities
    
    implicit none
    
    integer(kind=4)          :: i               ! index loop 
    integer(kind=4)          :: j               ! index loop 
    integer(kind=4)          :: nAgeBins_spec   ! nb of age steps in the stellar spectra table
    integer(kind=4)          :: nMetBins_spec  
     
    character(MAXPATHSIZE)   :: filename
    character(MAXPATHSIZE)   :: message         ! a message  
    character(MAXPATHSIZE)   :: line
    
    real(kind=8),allocatable :: AgeBins_spec(:)
    real(kind=8),allocatable :: MetBins_spec(:)
    
    call IO_print_message('stellar_population_read_stellar_spectra')

    ! Build the input filename, The initial mass function (IMF) is given in input    
    write(filename,'(a,a,a,a)') trim(stellar_population_input_path), '/stellar_spectra[BC03]', trim(IMF), '.in'
    write(message,'(a,a,a,a)') 'Load data from : ', 'stellar_spectra[BC03]', trim(IMF), '.in'
    call IO_print_message(message)
    ! Open the library file
    open(unit = stellarspectra_unit, file = filename, status = 'old')
       
    do
      !
      read(stellarspectra_unit, '(a)', end = 2) line
      if (trim(line) .eq. 'START') then
        !
        ! all lines before the START keywork are header lines: they contain informations about data 
        !
        if (main_process .or. luminous_process) then
           !
           read(stellarspectra_unit,*) nAgeBins_spec, nMetBins_spec, nWaves
        else
           !
           read(stellarspectra_unit,*) ! skip the line
        end if
        !
        if (main_process) then
          !
          ! a luminous process doesn't read nMetBins and nAgeBins
          ! but the main process 'yes' therfore the main process checks if nMetBins_spec = nMetBins and if nAgeBins_spec = nAgeBins.
          ! If it is not the case the main process kills the run
          if (nMetBins_spec .ne. nMetBins) then
            !
            call IO_print_error_message('nMetBins_spec /= nMetBins',only_rank=rank,called_by='stellar_population_read_stellar_spectra') 
            stop ! stop the program 
          end if
          !
          if (nAgeBins_spec .ne. nAgeBins) then
            !
            call IO_print_error_message('nAgeBins_spec /= nAgeBins',only_rank=rank,called_by='stellar_population_read_stellar_spectra') 
            call IO_print_message('use',param_name=(/'nAgeBins(spec)           ','nAgeBins(stellar)        '/), &
                    int_param_val=(/nAgeBins_spec,nAgeBins/))
            stop ! stop the program
          end if
        end if
        !
        if (luminous_process) then
          !
          ! If nMetBins_spec /= nMetBins and if nAgeBins_spec /= nAgeBins the main process kills the run
          ! Therefore if G.A.S. still runs we can set: 
          nAgeBins = nAgeBins_spec
          nMetBins = nMetBins_spec
        end if
        !
        ! allocate arrays
        if (main_process .or. luminous_process) then
           !
           allocate(AgeBins_spec(nAgeBins))
           if (.not. allocated(AgeBins)) allocate(AgeBins(nAgeBins))
           allocate(MetBins_spec(nMetBins))
           if (.not. allocated(MetBins)) allocate(MetBins(nMetBins))
           allocate(Waves(nWaves))   ! [in micron 1.e-6 m]
           allocate(sb_sed(nWaves,nAgeBins,nMetBins))
        end if
        !
        allocate(LumBins(nAgeBins,nMetBins))
        !
        if (main_process .or. luminous_process) then
          !
          ! read stellar ages
          read(stellarspectra_unit,*) AgeBins_spec
          ! Check age bins
          if (main_process) then
            !
            ! a luminous process doesn't read MetBins and AgeBins table
            ! but the main process 'yes' therfore the main process checks if AgeBins_spec(i) = AgeBins(i) for all i in nAgeBins
            ! If it is not the case the main process kills the run
            do i = 1, nAgeBins
               !
               if (AgeBins_spec(i) .ne. AgeBins(i)) then
                  !  
                  call IO_print_error_message('AgeBins_spec(i) /= AgeBins(i)',only_rank=rank,called_by='stellar_population_read_stellar_spectra') 
                  call IO_print_message('with',only_rank=rank, &
                            param_name = (/'iAge                     '/), int_param_val = (/i/))
                  call IO_print_message('with',only_rank=rank, &
                            param_name = (/'AgeBins_spec(i)          ','AgeBins(i)               '/), real_param_val = (/AgeBins_spec(i),AgeBins(i)/))
                  stop ! stop the program
               end if
            end do
          end if
          !
          ! Read metallicities
          read(stellarspectra_unit,*) MetBins_spec
          ! Read wavelenght
          read(stellarspectra_unit,*) Waves(1:nWaves)
        else
          ! 
          ! I am a physical process
          read(stellarspectra_unit,*)  ! skip the line
          read(stellarspectra_unit,*)  ! skip the line
          read(stellarspectra_unit,*)  ! skip the line
        end if
        !
        ! read elementary luminosities
        do i = 1, nAgeBins
          !
          read(stellarspectra_unit,*) LumBins(i,1:nMetBins)  ! [Lsun / Msun]
        end do
        ! convert in code unit
        LumBins = LumBins*mass_code_unit_in_M_Sun            ! [Lsun / (1.e11*Msun)]
        !
        if (luminous_process) then
          !
          ! If AgeBins_spec /= AgeBins the main process kills the run
          ! Therefore if G.A.S. still runs we can set: 
          AgeBins = AgeBins_spec
          MetBins = MetBins_spec
          !
          ! skip a line
          read(stellarspectra_unit,*)
          !
          ! read stellar spectra  
          do i = 1, nAgeBins
             !
             do j = 1, nWaves
                !
                read(stellarspectra_unit,*) sb_sed(j,i,1:nMetBins) ! [Lsun / Msun]
            end do
          end do
        end if
        !
        if (main_process) then
          !
          ! Print some informations 
          call IO_print_message('use', &
                param_name=(/'nWaves                   ','nAgeBins_spec            ','nMetBins_spec            '/),&
                int_param_val=(/nWaves,nAgeBins_spec,nMetBins_spec/))
       
        end if
        exit  ! quit do loop    
      end if
      if (line(1:1) .eq. '#') then
         !
         cycle ! header or something like this (skip)
      else
         !
         call IO_print_error_message('Impossible to read the line',only_rank=rank,called_by='stellar_population_read_stellar_spectra') 
         stop ! stop the program
      end if
    end do
    !
    ! erase local tables
    if (allocated(AgeBins_spec)) deallocate(AgeBins_spec)     
    if (allocated(MetBins_spec)) deallocate(MetBins_spec) 
    !    
2   close(stellarspectra_unit)
    
    return  
  end subroutine stellar_population_read_stellar_spectra
  
  !*****************************************************************************************************************

  subroutine stellar_population_finalize

    implicit none

    if (allocated(AgeBins))         deallocate(AgeBins)         ! ages corresponding to age step
    if (allocated(mass_loss_rates)) deallocate(mass_loss_rates) ! mass loss rate for given age and metallicity
    if (allocated(SN_rates))        deallocate(SN_rates)        ! SN event rate for a given age and metallicity

    if (allocated(Waves))   deallocate(Waves)   ! uv, opt, near ir wavelengh, used in stellar spectra (Waves) 
    if (allocated(sb_sed))  deallocate(sb_sed)  ! starburst spectrum (Waves, nAgeBins, nMetBins)
    
    return
  end subroutine stellar_population_finalize

  !*****************************************************************************************************************

end module stellar_population_library


