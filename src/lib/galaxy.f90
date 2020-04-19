module galaxy

  use disc    ! Contains disc structure and disc exploitation functions
  use filters ! acces to filters and subroutines to compute magnitudes

  public

  !*****************************************************************************************************************
  !
  !  OVERVIEW
  !
  !  galaxy_module defines the galaxy data structure univ(ts)%halo(ih)%galaxy
  !  This module contains the definitions and all functions and subroutines associated to a galaxy component
  !  A galaxy is composed of a disc and a bulge (if a merger event occurs in the past evolution of the galaxy)
  !  The galaxy component allows to link the disc(bulge) components to the baryon halo phase
  !  Global feedback is computed at the galaxy scale by taking into account all processes acting at the disc scale
  !  In the header of the module are defined output properties of a galaxy component (labels, units and formats)
  !
  !  MPI PROCESSES IN THIS MODULE
  !
  !   galaxy_send_data                               : send specific informations about a galaxy in a given halo
  !
  !   galaxy_receive_data                            : receive specific informations about a galaxy in a given halo
  !
  !  SUBROUTINES IN THIS MODULE
  !
  !   galaxy_set_reference_mass                      : initialize reference mass
  !      called by : halo_set_reference_parameters
  !
  !   galaxy_void                                    : init all galaxy properties
  !      called by : halo_void
  !
  !   galaxy_copy                                    : copy the galaxy properties into an other
  !      called by : halo_copy & halo_save
  !
  !   galaxy_deallocate                              : deallocate all arrays allocated in galaxy components (disc and bulge)
  !      called by : halo_deallocate
  !
  !   galaxy_compute_galaxy_feedback_activities      : returns all feedback properties associated to a given galaxy (SN, AGN)
  !      called by : galaxy_evolve
  !                : halo_evolve
  !
  !   galaxy_evolve                                  : evolve a galaxy component during dt
  !     called by : halo_evolve
  !
  !   galaxy_merge                                   : merge two galaxies, redistribute stars and gas, merge gsh_tab of the two progenitors
  !      called by : halo_merge
  !
  !   galaxy_spectrum                                : Build the complete spectral energy distribution of a given galaxy, compute AB Mags
  !
  !  FUNCTIONS IN THIS MODULE
  !
  !   galaxy_mass                                    : return the galaxy mass encolse in the radius r
  !
  !   galaxy_frac_mass_radius                        : returns the radius which enclose a given fraction of the galaxy mass
  !
  !   galaxy_gas_signature                           : Return the gas signature (elements mass fraction) of a given galaxy component (or a sub-component)
  !

  !   galaxy_compute_escape_velocity                 : returns the escape velocity of the galaxy + dm dynamical system
  !
  !  PRINTING PROCEDURES
  !
  !   galaxy_load_gal_data                           : create the galaxy output properties list
  !
  !   galaxy_print                                   : print galaxy properties in a binary file
  !      called by : halo_print
  !
  !*****************************************************************************************************************

  type galaxy_type
    ! Evolution
    real(kind=8)             :: life_time           ! time since galaxy formation (oldest progenitors)
    real(kind=8)             :: age_form            ! age of the universe when the galaxy (the first oldest prog) has been formed
    !
    ! main components
    type(disc_type)          :: disc                ! a disc
    type(bulge_type)         :: bulge               ! and a bulge
    ! galaxy global properties
    real(kind=8)             :: R50                 ! radius which enclose 50% of the galaxy mass
    real(kind=8)             :: R50_stars           ! radius which enclose 50% of the galaxy stellar mass
    real(kind=8)             :: Vesc                ! escape velocity
    real(kind=8)             :: Vwind               ! wind velocity
    real(kind=8)             :: mu_tot              ! mass ratio of galaxy mass and dark matter at the last merger event
    real(kind=8)             :: mu_gas              ! gas mass ratio of the two last merger progenitors
    ! mass assembly properties
    integer(kind=4)          :: nb_merger           ! number of merger
    integer(kind=4)          :: nb_major_merger     ! number of major merger
    real(kind=8)             :: accreted_gas_mass   ! all gas mass accreted by galaxy over the galaxy's history
    real(kind=8)             :: merger_gas_mass     ! all gas mass gained by merger event over the galaxy's history
    ! stellar properties
    real(kind=8)             :: Stel_Age(2)         ! Stellar Age (mass and lum weighted)
    real(kind=8)             :: Stel_dAge(2)        ! Error on the stellar Age (mass and lum weighted)
    real(kind=8)             :: Stel_Z(2)           ! Stellar metallicity (mass and lum weighted)
    real(kind=8)             :: Stel_dZ(2)          ! Error on the stellar metallicity (mass and lum weighted)
    ! spectra
    real(kind=4),allocatable :: spectra(:,:)        ! spectra of the global stellar population,
                                                    ! spectra(:,1) young stars, spectra(:,2) old stars, spectra(:,3) full galaxy spectrum
    ! Mags
    real(kind=4),allocatable :: Mags(:,:)           ! Magnitudes in the observer-frame and rest-frame
    real(kind=4),allocatable :: Mags_disc(:)        ! Magnitudes in the observer-frame for the disc only
    real(kind=4),allocatable :: dMags(:)            ! First order derivative magnitudes in the observer-frame
    real(kind=4),allocatable :: Non_Ext_Mags(:,:)   ! Non-extincted Magnitudes in the observer-frame and rest-frame
    !
    ! additionnal luminous properties
    real(kind=8)             :: Lbol_stars          ! bolometric luminosity of the whole stellar population
    real(kind=8)             :: Mdust               ! Mass of dust deduced from the total IR luminosity
    real(kind=8)             :: Lir(3)              ! Infrared luminosity [8 : 1000] micron 1: young stars contribution, 2: full galaxy, 3: AGN
    real(kind=8)             :: ISRF(2)             ! InterStellar Radiation Field [G0 Habing unit]
                                                    ! computed for two stellar population ages 1: young stars only / 2: old stars only
    real(kind=8)             :: AV(2)               ! A(V) computed in 1: the disc ISM, 2: the disc BC
    real(kind=8)             :: E_BV(2)             ! E(B-V) computed in 1: the disc ISM, 2: the disc BC
  end type galaxy_type

  ! hdu reference for galaxy structure
  integer(kind=4)           :: hdu_gal
  ! printable physical properties for galaxy structure
  integer(kind=4),parameter :: nb_gal_field = 17    ! Number of gal properties saved
  ! Name of each output colomn data
  character(len=ttype_len),dimension(nb_gal_field) :: ttype_gal = (/'nb_merger             ','nb_major_merger       ','tilt                  ',&
                                                                    'accreted_mass         ','merger_mass           ','m_stars               ',&
                                                                    'mu_tot                ','mu_gas                ','R50_stars             ',&
                                                                    'disc_bulge_ratio      ','gal_stars_age_MW      ','gal_stars_dage_MW     ',&
                                                                    'gal_stars_age_LW      ','gal_stars_dage_LW     ','gal_Z_stars_MW        ',&
                                                                    'gal_Z_stars_LW        ','gal_Vesc              '/)
  ! Physical unit of each column data
  character(len=tunit_len),dimension(nb_gal_field) :: tunit_gal = (/'#           ','#           ','radian      ','M_sun       ',&
                                                                    'M_sun       ','M_sun       ','w_o_unit    ','w_o_unit    ',&
                                                                    'kpc         ','w_o_unit    ','Gyr         ','Gyr         ',&
                                                                    'Gyr         ','Gyr         ','Z_sun       ','Z_sun       ',&
                                                                    'km/s        '/)
  ! Data type of each column data
  character(len=tform_len),dimension(nb_gal_field) :: tform_gal = (/'1J','1J','1E','1E','1E','1E','1E','1E','1E', &
                                                                    '1E','1E','1E','1E','1E','1E','1E','1E'/)
  !
  ! printable properties for spectra structures
  integer(kind=4),parameter             :: nb_spec_field = 1
  ! hdu reference for stellar population spectra
  integer(kind=4)                       :: hdu_ystars   ! for young stars
  integer(kind=4)                       :: hdu_ostars   ! for old stars
  integer(kind=4)                       :: hdu_fullspec ! for the full galaxy spectra
  ! header informations
  character(len=ttype_len)              :: ttype_spec
  character(len=tunit_len)              :: tunit_spec
  character(len=5)                      :: tform_spec
  !
  ! printable properties for Magnitudes
  integer(kind=4)                       :: nb_Mag_field
  ! hdu reference Magnitudes
  integer(kind=4)                       :: hdu_NE_Mag_obs, hdu_NE_Mag_gal, hdu_Mag_obs, hdu_Mag_disc_obs, hdu_dMag_obs, hdu_Mag_gal, hdu_Props
  ! header informations (depends of the number of filters in input)
  character(len=ttype_len),allocatable  :: ttype_mag(:)
  character(len=tunit_len),allocatable  :: tunit_mag(:)
  character(len=tform_len),allocatable  :: tform_mag(:)

  ! other printable luminous properties (Lir, ... )
  integer(kind=4),parameter             :: nlumprops = 11
  ! header informations
  character(len=ttype_len),dimension(nlumprops)  :: ttype_lumprops = (/'Lbol_stars            ','ISRF_disc_ISM         ','ISRF_BC               ', &
                                                                       'Lir                   ','Lir_young_stars       ','Lir_AGN               ',&
                                                                       'Mdust                 ','AV_disc_ISM           ','AV_BC                 ', &
                                                                       'E_BV_disc_ISM         ','E_BV_BC               '/)
  character(len=tunit_len),dimension(nlumprops)  :: tunit_lumprops = (/'Lsun        ','G0 unit     ','G0 unit     ','Lsun        ','Lsun        ',&
                                                                       'Lsun        ','Msun        ','Mag         ','Mag         ','Mag         ','Mag         '/)
  character(len=tform_len),dimension(nlumprops)  :: tform_lumprops = (/'1E','1E','1E','1E','1E','1E','1E','1E','1E','1E','1E'/)

contains

  !*****************************************************************************************************************
  !
  ! MPI PROCESSES
  !
  !*****************************************************************************************************************

#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------
  subroutine galaxy_send_data(gal)

    ! SEND SPECIFIC INFORMATIONS ABOUT A GALAXY IN A GIVEN HALO (spectrum)

    implicit none

    logical                      :: go_down

    type(galaxy_type),intent(in) :: gal      ! a galaxy

    ! data are sent by physical process (odd) and receive by luminous process (even)

    go_down = .false. ! init
    if (galaxy_mass(gal,component='stars') .gt. 0.d0) then
        ! the galaxy exist and it have a stellar population
        !
        go_down = .true. ! init
        ! send data
        call MPI_SEND(go_down,1,MPI_LOGICAL,rank+1,gal_tag+1,MPI_COMM_WORLD,ierror)
        !
        ! send stellar population of the disc
        call disc_send_data(gal%disc)
        !
        ! send stellar population of the bulge
        call bulge_send_data(gal%bulge)
        !
    else
        ! send exit loop order
        call MPI_SEND(go_down,1,MPI_LOGICAL,rank+1,gal_tag+1,MPI_COMM_WORLD,ierror)
    end if

    return
  end subroutine galaxy_send_data

  !*****************************************************************************************************************

  subroutine galaxy_receive_data(gal)

    ! RECEIVE SPECIFIC INFORMATIONS ABOUT A GALAXY IN A GIVEN HALO

    implicit none

    logical                         :: go_down ! exit loop order (if = .false. returnto tree level)

    type(galaxy_type),intent(inout) :: gal     ! a galaxy

    call galaxy_void(gal)                      ! init the galaxy

    ! receive exit loop order
    call MPI_RECV(go_down,1,MPI_LOGICAL,rank-1,gal_tag+1,MPI_COMM_WORLD,statut,ierror)
    !
    if (go_down) then
        ! the galaxy exists and has a global stellar population
        !
        ! receive stellar population of the disc
        call disc_receive_data(gal%disc)
        !
        ! receive stellar population of the bulge
        call bulge_receive_data(gal%bulge)
        !
    end if

    return
  end subroutine galaxy_receive_data
! -------------------------------------------------
#endif
! LUMINOUS_PROCESSES

  !*****************************************************************************************************************
  !
  ! SUBROUTINES
  !
  !*****************************************************************************************************************

  subroutine galaxy_set_reference_mass

    ! INITIALIZE REFERENCE PARAMETERS LINKED TO THE DARK-MATTER N-BODY SIMULATION

    implicit none

    ! initialize M_stars_min and young_stars_iMaxAge
    call stars_set_reference_mass_and_age

    return
  end subroutine galaxy_set_reference_mass

  !*****************************************************************************************************************

  subroutine galaxy_void(gal)

    ! INIT OR VOID A GALAXY COMPONENT

    implicit none

    type(galaxy_type),intent(inout) :: gal

    ! Evolution
    gal%life_time           = 0.d0    ! will be incremented in galaxy evolve
    gal%age_form            = 0.d0    ! will be set in galaxy_evolve during the fisrt evolution step
    !
    ! main components
    call disc_void(gal%disc)          ! void the disc component
    call bulge_void(gal%bulge)        ! void the bulge component
    ! galaxy global properties
    gal%R50                 = 0.d0    ! all other parameters are set to null values
    gal%Vesc                = 0.d0
    gal%Vwind               = 0.d0
    gal%mu_tot              = 0.d0
    gal%mu_gas              = 0.d0
    ! mass assembly properties
    gal%nb_merger           = 0       ! will be updated in galaxy_merge
    gal%nb_major_merger     = 0       ! will be updated in galaxy_merge
    gal%accreted_gas_mass   = 0.d0    ! will be incremented in galaxy_evolve
    gal%merger_gas_mass     = 0.d0    ! will be incremented in galaxy_merge
    ! stellar properties
    gal%Stel_Age            = 0.d0
    gal%Stel_dAge           = 0.d0
    gal%Stel_Z              = 0.d0
    gal%Stel_dZ             = 0.d0
    ! spectra
    if (allocated(gal%spectra))      gal%spectra = 0.d0
    ! Mags
    if (allocated(gal%Mags))         gal%Mags = 0.d0
    if (allocated(gal%Mags_disc))    gal%Mags_disc = 0.d0
    if (allocated(gal%dMags))        gal%dMags = 0.d0
    if (allocated(gal%Non_Ext_Mags)) gal%Non_Ext_Mags = 0.d0
    ! additionnal luminous properties
    gal%Lbol_stars          = 0.d0
    gal%Mdust               = 0.d0
    gal%Lir                 = 0.d0
    gal%ISRF                = 0.d0
    gal%AV                  = -1.d0
    gal%E_BV                = -1.d0

    return
  end subroutine galaxy_void

  !*****************************************************************************************************************

  subroutine galaxy_copy(gal1,gal2)

    ! COPY gal2 IN gal1

    implicit none

    type(galaxy_type),intent(inout) :: gal1
    type(galaxy_type),intent(in)    :: gal2

    ! Evolution
    gal1%life_time           = gal2%life_time
    gal1%age_form            = gal2%age_form
    !
    ! main properties
    call disc_copy(gal1%disc,gal2%disc)       ! copy the disc component
    call bulge_copy(gal1%bulge,gal2%bulge)    ! copy the bulge component
    gal1%R50                 = gal2%R50
    gal1%Vesc                = gal2%Vesc
    gal1%Vwind               = gal2%Vwind
    gal1%mu_tot              = gal2%mu_tot
    gal1%mu_gas              = gal2%mu_gas
    ! mass assembly properties
    gal1%nb_merger           = gal2%nb_merger
    gal1%nb_major_merger     = gal2%nb_major_merger
    gal1%accreted_gas_mass   = gal2%accreted_gas_mass
    gal1%merger_gas_mass     = gal2%merger_gas_mass
    ! stellar properties
    gal1%Stel_Age            = gal2%Stel_Age
    gal1%Stel_dAge           = gal2%Stel_dAge
    gal1%Stel_Z              = gal2%Stel_Z
    gal1%Stel_dZ             = gal2%Stel_dZ
    ! spectra
    if (allocated(gal2%spectra)) then
        if (.not. allocated(gal1%spectra)) then
            allocate(gal1%spectra(nWaves,2))
            gal1%spectra = gal2%spectra
        else
            call IO_print_error_message('Overwrite an existing gal%spectra',only_rank=rank,called_by='galaxy_copy')
            stop
        end if
    end if
    ! Mags
    if (allocated(gal2%Mags)) then
        if (.not. allocated(gal1%Mags)) then
            allocate(gal1%Mags(nfilters,2))
            gal1%Mags = gal2%Mags
        else
            call IO_print_error_message('Overwrite an existing gal%Mags',only_rank=rank,called_by='galaxy_copy')
            stop
        end if
    end if
    ! Mags_disc
    if (allocated(gal2%Mags_disc)) then
        if (.not. allocated(gal1%Mags_disc)) then
            allocate(gal1%Mags_disc(nfilters))
            gal1%Mags_disc = gal2%Mags_disc
        else
            call IO_print_error_message('Overwrite an existing gal%Mags_disc',only_rank=rank,called_by='galaxy_copy')
            stop
        end if
    end if
    if (allocated(gal2%dMags)) then
        if (.not. allocated(gal1%dMags)) then
            allocate(gal1%dMags(nfilters))
            gal1%dMags = gal2%dMags
        else
            call IO_print_error_message('Overwrite an existing gal%dMags',only_rank=rank,called_by='galaxy_copy')
            stop
        end if
    end if
    if (allocated(gal2%Non_Ext_Mags)) then
        if (.not. allocated(gal1%Non_Ext_Mags)) then
            allocate(gal1%Non_Ext_Mags(nfilters,2))
            gal1%Non_Ext_Mags = gal2%Non_Ext_Mags
        else
            call IO_print_error_message('Overwrite an existing gal%Non_Ext_ags',only_rank=rank,called_by='galaxy_copy')
            stop
        end if
    end if
    ! additionnal luminous properties
    gal1%Lbol_stars          = gal2%Lbol_stars
    gal1%Mdust               = gal2%Mdust
    gal1%Lir                 = gal2%Lir
    gal1%ISRF                = gal2%ISRF
    gal1%AV                  = gal2%AV
    gal1%E_BV                = gal2%E_BV

    return
  end subroutine galaxy_copy

  !*****************************************************************************************************************

  subroutine galaxy_deallocate(gal)

    ! DEALLOCATE ALL ARRAYS ALLOCATED IN THE DIFFERENTS GALAXY COMPONENT

    implicit none

    type(galaxy_type),intent(inout) :: gal  ! a galaxy

    ! spectra
    if (allocated(gal%spectra))      deallocate(gal%spectra)
    ! Mags
    if (allocated(gal%Mags))         deallocate(gal%Mags)
    if (allocated(gal%Mags_disc))    deallocate(gal%Mags_disc)
    if (allocated(gal%dMags))        deallocate(gal%dMags)
    if (allocated(gal%Non_Ext_Mags)) deallocate(gal%Non_Ext_Mags)
    ! galaxy components
    call disc_deallocate(gal%disc)
    call bulge_deallocate(gal%bulge)

    return
  end subroutine galaxy_deallocate

  !*****************************************************************************************************************

  subroutine galaxy_compute_galaxy_feedback_activities(gal,dm,ejecta_rate,agn_acc_rate,Vwind,Qtherm,Qturb,Qrad,f_in)

    ! RETURN THE TOTAL GALAXY EJECTA RATE COME FROM SN FEEDBACK AND/OR AGN FEEDBACK PROCESSES
    ! compute global ejecta rate, thermal power, radiative power and mean wind velocity
    ! In this new version, all ejecta are due to processes into the disc structure, the bulge is considered as a passive component
    ! Even if the creation and the evolution of a SMBH is linked to bulge properties,
    ! all effects produced by the SMBH are computed in the context of the disc structure (gas structuration and gas reservoir)

    implicit none

    real(kind=8),intent(out),optional :: ejecta_rate      ! the global ejecta rate generated by SN explosion and/or SMBH activity
    real(kind=8)                      :: ej_rate
    real(kind=8),intent(out),optional :: agn_acc_rate     ! the accretion rate onto the SMBH
    real(kind=8),intent(out),optional :: Vwind            ! the rate-weighted velocity of the wind produced by feedback processes
    real(kind=8)                      :: Vw               ! the rate-weighted velocity of the wind produced by feedback processes
    real(kind=8),intent(out),optional :: Qtherm           ! thermal power (allow to compute mean wind temperature)
    real(kind=8)                      :: Qt
    real(kind=8),intent(out),optional :: Qturb            ! turbulent heating power (allow to reduced the effective accretion rate)
    real(kind=8),intent(out),optional :: Qrad             ! non kinetic and non-thermal power (radiations: allow to reduced the effective cooling rate)
    real(kind=8),intent(out),optional :: f_in             ! predicted fraction of galaxy ejected-mass that is catched by the dark matter potentiel well
    real(kind=8)                      :: T_ej             ! predicted temperature of galaxy ejecta
    real(kind=8)                      :: Vesc             ! escape velocity associated to the dark matter halo

    type(gas_type)                    :: ej               ! ejecta
    type(galaxy_type),intent(in)      :: gal              ! the galaxy component
    type(dm_type),intent(in)          :: dm               ! the dm component

    ! disc_compute_disc_feedback_activities(disc,Vesc,ejecta_rate,agn_acc_rate,Vwind,Qtherm,Qrad)
    call disc_compute_disc_feedback_activities(gal%disc,gal%Vesc, &
          ejecta_rate=ej_rate,agn_acc_rate=agn_acc_rate,Vwind=Vw,Qtherm=Qt,Qturb=Qturb,Qrad=Qrad)
    !
    ! set ejecta_rate
    if (present(ejecta_rate)) ejecta_rate = ej_rate
    !
    ! set Vwind
    if (present(Vwind)) Vwind = Vw
    !
    ! set Qterm
    if (present(Qtherm)) Qtherm = Qt
    !
    if (present(f_in)) then
        f_in = 0.d0 ! init
        if (ej_rate .gt. 0.d0) then
            ! we assume a constant ejection process during StellarTimeStep
            ! compute ejected_mass
            ej = ej_rate*StellarTimeStep*disc_gas_signature(gal%disc,apply_as='rate_builder')
            ! set temperature of ejecta
            call gas_inject_termal_energy(ej,Qt*StellarTimeStep)
            T_ej = gas_temp(ej)
            ! by using Vesc, Vwind and T_ej we can compute a prediction for f_in
            Vesc = dm_escape_velocity(dm%r_core,dm)
            f_in = min(1.d0,max(0.d0,Ronbint(Maxwell_Boltzman_Vdist_shifted,0.d0,Vesc,(/T_ej,Vw/))))
            if (f_in .le. 1.d-2)  f_in = 0.d0
            if (f_in .gt. 9.9d-1) f_in = 1.d0
        end if
    end if

    return
  end subroutine galaxy_compute_galaxy_feedback_activities

  !*****************************************************************************************************************

  subroutine galaxy_evolve(gal,dm,dt,fresh_gas_acc_rate,gal_stripping_rate,gal_ejecta_rate_Wd,gal_ejecta_rate_Wu,f_in_Wd,f_in_Wu, &
        gal_ejecta_rate,Vwind,dt_optim,stop_before_end)

    ! EVOLVE A GALAXY COMPONENT
    !     - compute the evolution of the disc and the bulge component
    !         in respecting dynamical times of these components
    !     - compute galaxy global properties

    implicit none

    logical,intent(out)              :: stop_before_end            ! = .true. if the evolution has been stopped before the end of initial time-step (dt_optim < dt)
    !
    character(MAXPATHSIZE)           :: message                    ! a message to display
    !
    real(kind=8),intent(in)          :: dt                         ! time-step dt
    real(kind=8),intent(in)          :: f_in_Wd                    ! lower value of f_in
    real(kind=8),intent(in)          :: f_in_Wu                    ! upper value of f_in
    real(kind=8),intent(out)         :: Vwind                      ! (time/rate)-mean galaxy wind velocity during dt_optim
    real(kind=8),intent(out)         :: dt_optim                   ! max value of dt which respect the quasi-static equilibrium criterion for all galaxy component
    real(kind=8)                     :: time_before_end            ! time before the end of the time-step dt
    ! gal
    real(kind=8)                     :: gal_time                   ! galaxy evolution time (defined when disc and bulge are synchronized)
    real(kind=8)                     :: gal_dt_optim               ! optimal time-step for the galaxy
    real(kind=8)                     :: inst_gal_ejecta_rate       ! instantaneous ejecta rate of the galaxy
    real(kind=8)                     :: gal_ejecta_rate_test       ! (time)-mean galaxy ejecta rate (test version)
    real(kind=8)                     :: M_ej
    real(kind=8)                     :: Wd, Wu
    real(kind=8)                     :: inst_Vwind                 ! instantaneous wind velocity
    real(kind=8)                     :: inst_Qtherm                ! instantaneous thermal power
    real(kind=8)                     :: Qtherm                     ! (time)-mean thermal power
    real(kind=8)                     :: f_in                       ! fraction of ejecta that stay in the dark-matter potential well
    ! disc
    real(kind=8)                     :: disc_dt_optim              ! optimal adatative time-step for the disc component of the galaxy
    real(kind=8)                     :: disc_time                  ! evolve time of the disc
    real(kind=8)                     :: disc_next_stop             ! step mark for disc evolution
    ! bulge
    real(kind=8)                     :: bulge_dt_optim             ! optimal adaptative time-step for the bulge component
    real(kind=8)                     :: bulge_time                 ! evolve time of the bulge
    real(kind=8)                     :: bulge_next_stop            ! step mark for bulge evolution

    type(gas_type),intent(in)        :: fresh_gas_acc_rate         ! in presence of an accretion rate "cold_inflow_rate" due to the cold phase
                                                                   ! and of an accretion rate "cooling_inflow_rate" due to the colling of the hot phase
    type(gas_type),intent(in)        :: gal_ejecta_rate_Wd         ! ejecta rate produced by the galaxy can not be smaller than this first barrier
    type(gas_type),intent(in)        :: gal_ejecta_rate_Wu         ! and can not be larger than this second barrier
    type(gas_type),intent(inout)     :: gal_stripping_rate         ! the average gas stripping rate acting during dt
    type(gas_type),intent(out)       :: gal_ejecta_rate            ! (time)-mean ejecta rate produced by the galaxy (SN and AGN)
    type(gas_type)                   :: gal_ejecta                 ! cumulative mass of ejecta produced by the galaxy during dt
    type(gas_type)                   :: gal_stripped_gas           ! cumulative mass of stripped gas

    type(stars_type)                 :: stars_tmp                  ! A local copy of the global disc and bulge stellar population

    type(dm_type),intent(in)         :: dm                         ! in a DM halo dm

    type(galaxy_type),intent(inout)  :: gal                        ! Evolve the galaxy

#ifdef PRINTALL
! -------------------------------------------------
  call IO_print_message('galaxy_evolve',only_rank=rank,component='galaxy')
! -------------------------------------------------
#endif

    call gas_void(gal_ejecta_rate)       ! init
    call gas_void(gal_ejecta)            ! init
    call gas_void(gal_stripped_gas)      ! init
    !
    inst_gal_ejecta_rate  = 0.d0         ! init
    M_ej                  = 0.d0         ! init
    !
    time_before_end       = dt           ! init, the galaxy must evolve until dt
    disc_time             = 0.d0         ! init
    bulge_time            = 0.d0         ! init
    disc_dt_optim         = 0.d0         ! init
    bulge_dt_optim        = 0.d0         ! init
    disc_next_stop        = 0.d0         ! init
    bulge_next_stop       = 0.d0         ! init
    gal_time              = 0.d0         ! init
    Vwind                 = 0.d0         ! init the (time)-mean velocity wind
    Qtherm                = 0.d0         ! init the (time)-mean themal power associated to the wind phase
    stop_before_end       = .false.      ! init
    !
    ! AGE FORM
    if ((galaxy_mass(gal) .eq. 0.d0) .and. (gas_mass(fresh_gas_acc_rate) .gt. 0.d0)) then
       ! set gal%age_form
       gal%age_form = dm%age_form + dm%life_time
       ! print the starting state of the galaxy
       if (FOLLOW_UP .and. PR_FOLLOW_UP) call galaxy_print(follow_up_unit(current_index),'fits','phy',gal)
    end if

    ! time loop, the evolution process have to be run until dt
    do while (abs(time_before_end) .gt. num_accuracy)
      !
      ! ****************************************************************************************************
      ! PREDICTOR PART : compute input and output rate for the disc and the bulge component
      !                  compute optimal evolution time-step for these two components
      !                  these computation has perfomed by two dedicated subroutines: disc_evolve_I and bulge_evolve_I
      ! ****************************************************************************************************
      !
      ! ******************************************
      ! DISC
      ! the disc receives the fresh accreted gas, forms new stars etc ...
      ! ******************************************
      !
      if (disc_time .lt. dt) then
        if (abs(disc_time - gal_time) .lt. num_accuracy) then
          !
          if (disc_mass(gal%disc) .ge. 0.d0) then
            !
            ! compute optimal integration time for the disc component
            call disc_evolve_I(gal%disc,dm,gal%Vesc,fresh_gas_acc_rate,gal_stripping_rate,disc_dt_optim)
            !
            if (disc_dt_optim .gt. 0.d0) then
              disc_dt_optim = min((dt-disc_time),disc_dt_optim)
            else
              disc_dt_optim = dt-disc_time
            end if
          end if
          !
          ! FEEDBACK ACTIVITIES
          ! in the current model, only the disc host feedback activities
          ! the bulge is passive, the instantaneous ejecta rate is therefore given by the disc structure
          inst_gal_ejecta_rate = gas_mass(gal%disc%ejecta_rate)
          M_ej = gas_mass(gal_ejecta)
          Wu   = gas_mass(gal_ejecta_rate_Wu)
          Wd   = gas_mass(gal_ejecta_rate_Wd)
          !
          ! compute mean ejecta rate
          gal_ejecta_rate_test = (M_ej+inst_gal_ejecta_rate*disc_dt_optim)/(disc_time+disc_dt_optim)
          ! in the previous formulae, M_ej is the total gas mass already ejected (during disc_time)
          !
          if ((gal_ejecta_rate_test .gt. Wu) .or. (gal_ejecta_rate_test .lt. Wd)) then
            ! the halo equilibrium is modified
            stop_before_end = .true.
            exit
          end if
          ! @ this time the (time)-average galaxy ejecta rate is compatible
          ! with the two ejecta barriers given in input
          ! and this, until the next disc_next_stop flag
          !
          if (inst_gal_ejecta_rate .gt. 0.d0) then
             !
             ! compute (time-)mean values for feedback effects
             ! there are some ejected mass
             ! galaxy_compute_galaxy_feedback_activities(gal,dm,[ejecta_rate,agn_acc_rate,Vwind,Qtherm,Qturb,Qrad])
             call galaxy_compute_galaxy_feedback_activities(gal,dm,Vwind=inst_Vwind,Qtherm=inst_Qtherm,f_in=f_in)
             !
             ! test f_in
             if ((f_in .lt. f_in_Wd) .or. (f_in .gt. f_in_Wu)) then
                ! the halo equilibrium is modified
                stop_before_end = .true.
                exit
             end if
             ! @ this time f_in is compatible
             ! with the predicted value
             ! and this, until the next disc_next_stop flag
             !
             ! update integrated value
             gal_ejecta = gal_ejecta + gal%disc%ejecta_rate*disc_dt_optim               ! total ejected mass
             Vwind      = Vwind      + inst_Vwind*inst_gal_ejecta_rate*disc_dt_optim    ! wind velocity
             Qtherm     = Qtherm     + inst_Qtherm*disc_dt_optim                        ! thermal power
             !
          end if
          !
          ! STRIPPING OF THE DIFFUSE GAS
          ! compute the total stripped gas mass
          gal_stripped_gas = gal_stripped_gas + gal%disc%stripping_rate*disc_dt_optim
        end if
      end if
      !
      ! ******************************************
      ! BULGE
      ! The bulge is considered as a passive component
      ! The evolution rate is only based on stellar passive cycle (ISM enrichment but no star formation)
      ! ******************************************
      !
      if (bulge_time .lt. dt) then
        if (abs(bulge_time - gal_time) .lt. num_accuracy) then
          !
          if (bulge_mass(gal%bulge) .gt. 0.d0) then
            !
            ! compute optimal integration time for the bulge component
            call bulge_evolve_I(gal%bulge,bulge_dt_optim)
            !
            if (bulge_dt_optim .gt. 0.d0) then
              bulge_dt_optim = min((dt-bulge_time),bulge_dt_optim)
            else
              bulge_dt_optim = dt-bulge_time
            end if
          else ! there are no bulge
            bulge_dt_optim  = dt-bulge_time ! init
          end if
        end if
      end if
      !
      ! ********************************************
      ! COMPUTE OPTIMAL TIME-STEP
      ! ********************************************
      !
      ! disc
      disc_next_stop  = min(dt,disc_time + disc_dt_optim)     ! stop evolve mark
      disc_dt_optim   = disc_next_stop - disc_time            ! computed corresponding time-step
      !
      if (disc_dt_optim .le. 0.d0) then
        call IO_print_error_message('disc_dt_optim <= 0.',only_rank=rank,called_by='galaxy_evolve')
        call IO_print_message('used',only_rank=rank,component='disc', &
                param_name = (/'disc_dt_optim            ','dt                       ','time_before_end          ', &
                               'disc_next_stop           ','disc_time                ','gal_time                 ', &
                               'bulge_time               '/), &
                real_param_val  = (/disc_dt_optim,dt,time_before_end,disc_next_stop,disc_time,gal_time,bulge_time/))
        write(*,*) 'sbe: ', stop_before_end
        stop
      end if
      !
      ! bulge
      bulge_next_stop = min(dt,bulge_time + bulge_dt_optim)   ! stop evolve mark
      bulge_dt_optim  = bulge_next_stop - bulge_time          ! computed corresponding time-step
      !
      if (bulge_dt_optim .le. 0.d0) then
        call IO_print_error_message('bulge_dt_optim <= 0.',only_rank=rank,called_by='galaxy_evolve')
#ifdef PRINT_WARNING
! -------------------------------------------------
        call IO_print_message('used',only_rank=rank,component='disc', &
                param_name = (/'bulge_dt_optim           ','dt                       ', &
                               'bulge_next_stop          ','bulge_time               '/), &
                real_param_val  = (/bulge_dt_optim,dt,bulge_next_stop,bulge_time/))
! -------------------------------------------------
#endif
! PRINT_WARNING
        stop
      end if
      !
      ! @ this point the next time-step ddt can be add to the global evolution
      ! In the current model, feedback processes are only present in the disc structure
      ! If, during dt, the mean galaxy ejecta rate violates one of the two barriers, the evolution process must be stopped
      ! In this context the bulge_time cannot be larger than the disc_time
      gal_dt_optim = min(disc_next_stop,bulge_next_stop) - gal_time
      gal_time     = gal_time + gal_dt_optim
      !
      ! ****************************************************************************************************
      ! CORRECTOR PART : apply input and output rate for the disc and the bulge component
      !                  onto the optimal evolution time-step pre-computed for these two components
      !                  these computation has perfomed by two dedicated subroutines: disc_evolve_II and bulge_evolve_II
      ! ****************************************************************************************************
      !
      ! ********************************************
      ! DISC
      ! ********************************************
      !
      if (abs(disc_next_stop - gal_time) .lt. num_accuracy) then
        call disc_evolve_II(gal%disc,dm,gal%bulge,disc_dt_optim)
        disc_time = disc_next_stop
        gal_time  = disc_next_stop
      end if
      !
      ! ********************************************
      ! BULGE
      ! ********************************************
      !
      if (abs(bulge_next_stop - gal_time) .lt. num_accuracy) then
        call bulge_evolve_II(gal%bulge,bulge_dt_optim)
        bulge_time = bulge_next_stop
        gal_time   = bulge_next_stop
      end if
      !
      ! update time_before_end
      time_before_end = dt - gal_time
      !
    end do
    ! END TIME EVOLUTION LOOP
    !
    ! ********************************************
    ! REAL TIME STEP EVOLUTION AND CHECKS
    ! ********************************************
    !
    dt_optim = 0.d0  ! init
    if (time_before_end .le. num_accuracy) then
      ! THE EVOLUTION HAS BEEN PERFORM OVER ALL dt
      dt_optim = dt
    else
      ! THE EVOLUTION HAS BEEN STOPPED (ejecta effects)
      dt_optim = dt - time_before_end
    end if
    ! Check bulge evolution
    ! The time evolution of the passive bulge can be smaller than the disc evolution time
    if (bulge_time .lt. disc_time) then
      ! evolve the bulge component until disc_time
      call bulge_evolve_II(gal%bulge,(disc_time-bulge_time))
      bulge_time = disc_time
    end if
    if (dt_optim .gt. dt) then
      call IO_print_error_message('dt_optim > dt',only_rank=rank,called_by='galaxy_evolve')
      call IO_print_message('used',only_rank=rank,component='galaxy', &
                param_name=(/'dt_optim                 ','dt                       '/), &
                real_param_val=(/dt_optim,dt/))
      stop ! stop the program
    end if
    if (dt_optim .le. 0.d0) then
      call IO_print_error_message('dt_optim <= 0.d0',only_rank=rank,called_by='galaxy_evolve')
      call IO_print_message('used',only_rank=rank,component='galaxy', &
                param_name=(/'dt                       ','gal_dt_optim             ','gal_time                 ', &
                             'disc_next_stop           ','disc_dt_optim            ','disc_time                ', &
                             'bulge_next_stop          ','bulge_dt_optim           ','bulge_time               ', &
                             'stripping_rate           ','inst_gal_ejecta_rate     '/), &
                real_param_val=(/dt,gal_dt_optim,gal_time,disc_next_stop, &
                                 disc_dt_optim,disc_time,bulge_next_stop,bulge_dt_optim,bulge_time, &
                                 gal%disc%stripping_rate%mass, inst_gal_ejecta_rate/))
      write(*,*) 'sbe: ', stop_before_end
      stop ! stop the program
    end if
    !
    ! compare disc and bulge evolution time
    if (abs(disc_time - bulge_time) .gt. num_accuracy) then
      call IO_print_error_message('disc_time != bulge_time',only_rank=rank,called_by='galaxy_evolve')
      if (stop_before_end) call IO_print_message('galaxy evolution has been stop !! ',only_rank=rank,component='galaxy')
      call IO_print_message('used',only_rank=rank,component='galaxy', &
                param_name=(/'dt                       ','gal_dt_optim             ','gal_time                 ', &
                             'disc_next_stop           ','disc_dt_optim            ','disc_time                ', &
                             'bulge_next_stop          ','bulge_dt_optim           ','bulge_time               '/), &
                real_param_val=(/dt,gal_dt_optim,gal_time,disc_next_stop, &
                                 disc_dt_optim,disc_time,bulge_next_stop,bulge_dt_optim,bulge_time/))
      stop ! stop the program
    end if
    !
    ! ********************************************
    ! UPDATE GLOBAL FEEDBACK PROPERTIES
    ! ********************************************
    !
    ! EJECTA
    if (gas_mass(gal_ejecta) .gt. 0.d0) then
        ! mean wind velocity
        Vwind  = Vwind / gas_mass(gal_ejecta)
        !
        ! Set temperature of the ejecta
        ! We use gal_ejecta, the integrated gas mass ejected during dt_optim
        call gas_inject_termal_energy(gal_ejecta,Qtherm)
        if (gas_temp(gal_ejecta) .lt. diffuse_gas_temp) then
            write(message,'(a)') 'T_ej < T_cool_gas; set T_ej gas to T_cool'
            call IO_print_warning_message(message,only_rank=rank,called_by='galaxy_evolve')
            call gas_set_component(gal_ejecta,diffuse_gas_temp,component='Temp')
        end if
        !
        ! set gal_ejecta_rate, the mean effective ejecta rate acting during dt_optim
        gal_ejecta_rate = (1.d0/dt_optim)*gal_ejecta
    end if
    !
    if ((gas_mass(gal_ejecta_rate) .gt. 0.d0) .and. ((Qtherm .le. 0.d0) .or. (Vwind .le. 0.d0))) then
      call IO_print_error_message('No consistency m_ejected < 0. and Vwind(or Qtherm) > 0.',only_rank=rank,called_by='galaxy_evolve')
      call IO_print_message('used',only_rank=rank,component='galaxy', &
                param_name=(/'dt                       ','gal_ejecta_rate          ', &
                             'Vwind                    ','Qtherm                   '/), &
                real_param_val=(/dt,gas_mass(gal_ejecta_rate),Vwind,Qtherm/))
      stop ! stop the program
    end if
    !
    ! STRIPPING
    if (gas_mass(gal_stripped_gas) .gt. 0.d0) then
        ! set gal_stripping_rate, the mean effective stripping rate acting during dt_optim
        gal_stripping_rate = (1.d0/dt_optim)*gal_stripped_gas
    end if
    !
    ! ********************************************
    ! UPDATE GALAXY PROPERTIES
    ! ********************************************
    !
    if (galaxy_mass(gal) .gt. 0.d0) then
        ! galaxy properties
        ! Evolution
        ! update life_time
        gal%life_time = gal%life_time + dt_optim
        ! galaxy global properties
        gal%R50_stars = galaxy_frac_mass_radius(gal,5.d-1,component='stars',called_by='galaxy_evolve') ! radius which enclose 50% of the galaxy stellar mass
        gal%R50       = galaxy_frac_mass_radius(gal,5.d-1,called_by='galaxy_evolve')                   ! radius which enclose 50% of the galaxy mass (stars + gas)
        gal%Vesc      = galaxy_compute_escape_velocity(gal,dm)                                         ! escape velocity of the galaxy
        gal%Vwind     = Vwind                                                                          ! (time-)mean ejecta wind velocity
        ! mass assembly properties
        ! update the gas mass accreted by galaxy over the galaxy's history
        gal%accreted_gas_mass = gal%accreted_gas_mass + gas_mass(fresh_gas_acc_rate*dt_optim)
        ! global stellar properties (Age and Z)
        call stars_void(stars_tmp)                   ! Create the global stellar population
        stars_tmp = gal%disc%stars + gal%bulge%stars !
        gal%Stel_Age  = stars_tmp%Age                ! Copy the stellar age
        gal%Stel_dAge = stars_tmp%dAge               ! Error on the stellar age
        gal%Stel_Z    = stars_tmp%Z                  ! Stellar metallicity
        gal%Stel_dZ   = stars_tmp%dZ                 ! Error on the stellar metallicity
        call stars_void(stars_tmp)                   ! erase the local copy
        !
        if (FOLLOW_UP .and. PR_FOLLOW_UP) call galaxy_print(follow_up_unit(current_index),'fits','phy',gal)
    end if

#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('galaxy_evolve ... done',only_rank=rank,component='galaxy')
! -------------------------------------------------
#endif

    return
  end subroutine galaxy_evolve

  !*****************************************************************************************************************

  subroutine galaxy_merge(gal1,gal2,dm1,dm2)

    ! MERGE TWO GALAXIES
    ! the remnant galaxy, result of the merger, is gal1
    ! this galaxy evolves in the dark-matter halo dm1

    implicit none

    real(kind=8)                    :: r1,r2             ! galaxy half mass radius
    real(kind=8)                    :: M1,M2             ! galaxy + dark matter halo mass
    real(kind=8)                    :: Mg,Ms             ! gas mass and stellar mass
    real(kind=8)                    :: mu, mu_gas        ! progenitor mass ratio (gal +dm) and (gas)
    real(kind=8)                    :: f                 ! mass fraction
    real(kind=8)                    :: Mgal
    real(kind=8)                    :: Mprog1, Mprog2
    real(kind=8)                    :: error

    type(gas_type)                  :: unstr_from_bulges ! gas transfered from bulges (1 &2) to the new gaseous disc
    type(gas_type)                  :: unstr_in_torus    ! gas added to the AGN torus
    type(gas_type)                  :: gas               ! a gas component
    type(stars_type)                :: stars_tmp         ! A local copy of the global disc and bulge stellar population

    type(dm_type),intent(inout)     :: dm1               ! the DM haloes of galaxy1
    type(dm_type),intent(in)        :: dm2               ! the DM haloes of galaxy2

    type(galaxy_type),intent(inout) :: gal1              ! a galaxy
    type(galaxy_type),intent(in)    :: gal2              ! an other galaxy
    type(galaxy_type)               :: gal               ! a local copy of the remnant galaxy

#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('galaxy_merge',only_rank=rank,component='gal')
! -------------------------------------------------
#endif

    ! save mass of progenitors
    Mprog1 = galaxy_mass(gal1)
    Mprog2 = galaxy_mass(gal2)
    ! init
    call gas_void(unstr_from_bulges)
    call gas_void(unstr_in_torus)
    !
    if (Mprog1 .gt. 0.d0) then
      ! the first galaxy exist
      if (Mprog2 .gt. 0.d0) then
        ! the second galaxy exist, it is a real merger event
        ! create a local remnant galaxy
        call galaxy_void(gal)
        ! compute the merger factor
        ! if mu > epsilon_merge the merger is considered as a major merger
        r1 = gal1%R50
        r2 = gal2%R50
        ! check stellar and star forming gas
        M1 = galaxy_mass(gal1,r=r1) + 2.d0*dm_mass(dm1,r=r1)   ! dynamical mass of the first system
        M2 = galaxy_mass(gal2,r=r2) + 2.d0*dm_mass(dm2,r=r2)   ! dynamical mass of the second system
        mu = min(M1,M2)/max(M1,M2)                             ! merger factor
        ! computed gas merger factor
        Mg = galaxy_mass(gal1,component='gas') + galaxy_mass(gal2,component='gas')
        Ms = galaxy_mass(gal1,component='stars') + galaxy_mass(gal2,component='stars')
        mu_gas = Mg/(Mg+Ms)
        !
        ! MERGER PROPERTIES
        ! associated to the dark-matter halo and to the galaxy
        ! update the merger_gas_mass, integration of the gas mass gained by merger
        ! in the current model, it exist three gas component, two in the disc and one in the bulge
        gal%merger_gas_mass = gal1%merger_gas_mass + galaxy_mass(gal2,component='gas')
        gal%nb_merger       = max(gal1%nb_merger,gal2%nb_merger) +1
        gal%nb_major_merger = max(gal1%nb_major_merger,gal2%nb_major_merger) ! init
        ! update the number of major merger event if mu > epsilon_merge
        if (mu .ge. epsilon_merge) gal%nb_major_merger = gal%nb_major_merger +1
        ! updated the accreted gas mass
        ! for the remnant structure, the accreted gas mass is the mass accreted by the firts progenitor
        gal%accreted_gas_mass = gal1%accreted_gas_mass
        if (gal1%disc%t_since_last_merger .eq. 0.d0) then
          !
          ! gal1%disc%t_since_last_merger = 0., it is a multi-merger event
          ! we save the most important merger events of the sequence
          gal%mu_tot = max(mu,gal1%mu_tot)
          gal%mu_gas = max(mu_gas,gal1%mu_gas)
          !
        else
          !
          ! gal1%t_since_last_merger > 0.
          ! it the first merger (may be in a multi-merger event but at this time it is impossible to kwown)
          gal%mu_tot = mu
          gal%mu_gas = mu_gas
        end if
        !
        if (mu .gt. epsilon_merge) then
          !
          ! ****************************************************
          ! MAJOR MERGER
          ! the morphology of the stellar component is strongly affected
          ! stars are distributed in a spheroid
          ! gas coming from the two progenitors are added and form a new disc
          ! the geometrical properties of this new disc are computed used the current dark-matter halo
          ! the gas structuration is strongly modified ...
          ! ****************************************************
          !
          ! REMNANT GALAXY PROPERTIES
          ! compute typical Hernquist size of the spheroid
          gal%bulge%rb = ((M1+M2)**2./(M1**2./r1+M2**2./r2+2.5d0*M1*M2/(r1+r2))) ! half mass radius
          gal%bulge%rb = gal%bulge%rb/(1.d0+sqrt(2.d0))                          ! bulge scale radius
          ! compute exponential radius of the remnant disc
          gal%disc%rd = (M1*gal1%disc%rd + M2*gal2%disc%rd)/(M1+M2)    ! mass weighted exponential radius
          gal%disc%rd = max(gal%disc%rd,dm1%spin*dm1%R_vir/sqrt(2.d0)) ! max value with dm spin properties
          ! compute the new disc inclination
          call disc_compute_inclination(gal%disc,dm1,disc1=gal1%disc,disc2=gal2%disc,mu=mu)
          !
          ! STARS
          ! 1th step: all stellar components are added to the remnant bulge component
          gal%bulge%stars = gal1%bulge%stars + gal1%disc%stars + gal2%disc%stars + gal2%bulge%stars
          ! 2th step: young stars are transfered to the disc structure
          call stars_transfer_young_stars(gal%bulge%stars,gal%disc%stars)
          !
        else
          !
          ! ****************************************************
          ! MINOR MERGER
          ! during a minor merger the global morphology of the galaxy is not strongly affected
          ! stellar discs and speroids are kept but all the gas (evolving in the different components (bulges and disc)
          ! are added and associated to the larger progenitor disc
          ! ****************************************************
          !
          ! REMNANT GALAXY PROPERTIES
          ! compute exponential radius of the remnant disc
          gal%disc%rd = (M1*gal1%disc%rd + M2*gal2%disc%rd)/(M1+M2)  ! mass weighted exponential radius
          gal%bulge%rb = max(gal1%bulge%rb,gal2%bulge%rb)            ! the larger bulge survives
          ! compute the new disc inclination
          call disc_compute_inclination(gal%disc,dm1,disc1=gal1%disc,disc2=gal2%disc,mu=mu)
          !
          ! STARS
          ! i) stellar population of the two discs are added
          gal%disc%stars = gal1%disc%stars + gal2%disc%stars
          ! ii) stellar population of the two bulges are added
          gal%bulge%stars = gal1%bulge%stars + gal2%bulge%stars
        end if
        !
        ! GAS
        ! @ this point the new caracteristics of the remnant disc are computed and save in gal%disc
        ! in the two cases, minor or major mergers,
        ! all gas components coming from the two progenitors are added. They form a new gas disc
        unstr_from_bulges = gal1%bulge%gas + gal2%bulge%gas
        !
        ! AGN
        if ((agn_mass(gal1%disc%agn) .gt. 0.d0) .or. (agn_mass(gal2%disc%agn) .gt. 0.d0)) then
            ! One or the two progenitors have already a SMBH
            ! Merge
            gal%disc%agn = agn_merge(gal1%disc%agn,gal2%disc%agn)
        else
            ! No pre-existing SMBH
            ! try to create it
            if (bulge_mass(gal%bulge) .gt. 1.d4*M_BH_min) then
                ! compute the fraction of the stellar mass which will be turned into a black-hole
                f = M_BH_min / bulge_mass(gal%bulge,component='stars')
                ! the the black hole mass is composed of stellar stars
                gal%bulge%stars = gal%bulge%stars*(1.d0-f)
                !
                ! create agn component
                call agn_create_agn(gal%disc%agn,M_BH_min)
            end if
        end if
        if (agn_mass(gal%disc%agn) .gt. 0.d0) then
            ! a SMBH exist or has been created
            ! reset the time elapsed since the last merger event
            call agn_reset_t_since_last_merger(gal%disc%agn)
            !
            ! add mass to the torus
            gas =  disc_mass(gal1%disc,r=3.d0*r_torus,component='unstr')*disc_gas_signature(gal1%disc,component='unstr') + &
                    disc_mass(gal2%disc,r=3.d0*r_torus,component='unstr')*disc_gas_signature(gal2%disc,component='unstr')
            if (gas_mass(gas) .gt. M_BH_min) then
                unstr_in_torus = mu*mu_gas*gas
                ! this gas will be substracted in disc_update_gas_struct_history
                ! feed the gas torus
                call agn_add_torus_mass(gal%disc%agn,unstr_in_torus)
            end if
            !
            ! set the new dynamical time of the accretion process
            call agn_set_dynamical_time(gal%disc%agn,r_torus/disc_velocity(r_torus,dm1,gal%disc,gal%bulge))
            !
        end if
        !
        ! Gas structuration history
        call disc_update_gas_struct_history(gal%disc,disc1=gal1%disc,disc2=gal2%disc, &
                    unstr_from_bulges=unstr_from_bulges, &
                    unstr_in_torus=unstr_in_torus)
        !
        ! OTHER PROPERTIES
        ! for disc
        ! the time life of the remnant disc is set to the maximum value of the two progenitors
        gal%disc%life_time = max(gal1%disc%life_time,gal2%disc%life_time)
        ! age_form is set to the minimum value of the two progenitors
        gal%disc%age_form  = min(gal1%disc%age_form,gal2%disc%age_form)
        ! we reset t_since_last_merger
        gal%disc%t_since_last_merger = 0.d0
        ! update cooling clock
        call disc_update_cooling_clock(gal%disc,disc1=gal1%disc,disc2=gal2%disc)
        ! computed the orbital velocity at the half mass radius
        gal%disc%V         = disc_velocity(1.68d0*gal%disc%rd,dm1,gal%disc,gal%bulge)   
        ! compute the epicyclic frequency at the disc half mass radius
        gal%disc%kappa     = disc_kappa(1.68d0*gal%disc%rd,dm1,gal%disc,gal%bulge)      
        ! update the structuration fraction
        gal%disc%f_str     = disc_gas_fraction(gal%disc,component='structured')
        ! computed dynamical time
        gal%disc%t_dyn     = disc_dynamical_time(gal%disc,dm1,gal%bulge)
        ! compute velocity dispersion
        gal%disc%dV        = disc_update_velocity_dispersion(gal1%disc,disc2=gal2%disc)
        ! compute the scale height of the disc
        gal%disc%h         = disc_scale_height(gal%disc)
        ! update inertial cascade properties t_emp, t_form, t_cascade, ngc
        call disc_update_inertial_cascade(gal%disc,disc1=gal1%disc,disc2=gal2%disc)
        !
        ! for bulge
        ! set the bulge formation time and life time
        if (bulge_mass(gal%bulge) .gt. 0.d0) then
          ! the remant galaxy has a bulge
          if (bulge_mass(gal1%bulge) .gt. 0.d0) then
            ! gal1 has a bulge
            if (bulge_mass(gal2%bulge) .gt. 0.d0) then
              ! gal2 has a bulge
              ! the two progenitors have a bulge
              gal%bulge%age_form  = min(gal1%bulge%age_form,gal2%bulge%age_form)
              gal%bulge%life_time = max(gal1%bulge%life_time,gal2%bulge%life_time)
            else
              ! only gal1 has a bulge
              gal%bulge%age_form  = gal1%bulge%age_form
              gal%bulge%life_time = gal1%bulge%life_time
            end if
          else
            ! gal1 has no bulge
            if (bulge_mass(gal2%bulge) .gt. 0.d0) then
              ! gal2 has a bulge
              gal%bulge%age_form  = gal2%bulge%age_form
              gal%bulge%life_time = gal2%bulge%life_time
            else
              ! no progenitor with a bulge
              gal%bulge%age_form = gal1%age_form + gal1%life_time ! = gal2%age_form + gal2%life_time = dm1%age_form + dm1%life_time
              gal%bulge%life_time = 0.d0
            end if
          end if
        end if
        !
        ! GALAXY PROPERTIES
        ! the time life of the galaxy is set to the maximum value of the two progenitors
        gal%life_time = max(gal1%life_time,gal2%life_time)
        ! age_form is set to the minimum value of the two progenitors
        gal%age_form  = min(gal1%age_form,gal2%age_form)
        ! compute other global properties
        ! the half-mass radius
        gal%R50_stars = galaxy_frac_mass_radius(gal,5.d-1,called_by='galaxy_merge',component='stars')
        gal%R50       = galaxy_frac_mass_radius(gal,5.d-1,called_by='galaxy_merge')
        ! the escape velocity
        gal%Vesc      = galaxy_compute_escape_velocity(gal,dm1)
        ! global stellar properties (Age and Z)
        call stars_void(stars_tmp)                   ! Create the global stellar population
        stars_tmp = gal%disc%stars + gal%bulge%stars !
        gal%Stel_Age  = stars_tmp%Age                ! Copy the stellar age
        gal%Stel_dAge = stars_tmp%dAge               ! Error on the stellar age
        gal%Stel_Z    = stars_tmp%Z                  ! Stellar metallicity
        gal%Stel_dZ   = stars_tmp%dZ                 ! Error on the stellar metallicity
        call stars_void(stars_tmp)                   ! erase the local copy
        !
        if (FOLLOW_UP) then
          if (PR_FOLLOW_UP) then
            ! in this case, the gal2 is the followed galaxy,
            ! we have to erase properties associated to gal1 and replace them by gal2 properties
            gal%nb_merger         = gal2%nb_merger +1
            gal%nb_major_merger   = gal2%nb_major_merger                             ! init
            if (mu .ge. epsilon_merge)  gal%nb_major_merger = gal%nb_major_merger +1 ! update
            gal%life_time         = gal2%life_time
            gal%age_form          = gal2%age_form
            gal%disc%life_time    = gal2%disc%life_time
            gal%disc%age_form     = gal2%disc%age_form
            gal%bulge%life_time   = gal2%bulge%life_time
            gal%bulge%age_form    = gal2%bulge%age_form  ! init
            if ((bulge_mass(gal%bulge) .gt. 0.d0) .and. (bulge_mass(gal2%bulge) .eq. 0.d0)) then
                ! formation of a bulge component in the followed galaxy
                ! set the formation age of this new bulge
                gal%bulge%age_form = gal2%age_form + gal2%life_time
            end if
            gal%accreted_gas_mass = gal2%accreted_gas_mass
            gal%merger_gas_mass   = gal2%merger_gas_mass + galaxy_mass(gal1,component='gas')
          end if
        end if
        !
        ! erase gal1
        call galaxy_void(gal1)
        ! copy the remnant galaxy (gal) in gal1
        call galaxy_copy(gal1,gal)
        ! erase the local copy gal
        call galaxy_void(gal)
        !
      end if ! gal2 > 0.
    else
       ! it the first copy in the descendent merger tree halo
       call galaxy_copy(gal1,gal2)
    end if
    !
    ! check merger events stats
    if (gal1%nb_merger .lt. gal1%nb_major_merger) then
      call IO_print_error_message('nb mergers < nb major mergers',only_rank=rank,called_by='galaxy_merge')
      call IO_print_message('used',only_rank=rank,component='galaxy', &
                param_name=(/'nb_merger                ','nb_major_merger          '/), &
                int_param_val=(/gal1%nb_merger,gal1%nb_major_merger/))
      stop ! stop the program
    end if
    !
    ! check mass conservation
    Mgal = galaxy_mass(gal1)
    error = abs(Mgal-(Mprog1+Mprog2))/(Mprog1+Mprog2)
    if (error .gt. num_accuracy) then
       call IO_print_error_message('Mgal != Mprogs', &
                 only_rank = rank, called_by = 'galaxy_merge')
       call IO_print_message('with',only_rank=rank,component='gal',&
                 param_name=(/'dM                       ','error                    ', &
                              'Mgal                     ','Mprogs                   ', &
                              'Mprog1                   ','Mprog2                   '/), &
                 real_param_val =(/Mgal-(Mprog1+Mprog2),error,Mgal,Mprog1+Mprog2,Mprog1,Mprog2/))
       stop ! stop the program
    end if

#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('galaxy_merge ... done ',only_rank=rank,component='gal')
! -------------------------------------------------
#endif

    return
  end subroutine galaxy_merge

#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------

  subroutine galaxy_spectrum(gal,z)

    ! BUILD THE GALAXY SPECTRUM (STARS + DUST)

    implicit none

    real(kind=8),intent(in)               :: z                       ! redshift of the galaxy
    real(kind=8)                          :: Lbol_before, Lbol_after
    real(kind=8)                          :: min_test
    real(kind=8)                          :: ISRF
    real(kind=8)                          :: Lir, Lir_AGN
    real(kind=8)                          :: Mdust, Chi
    real(kind=8)                          :: spec_min, spec_max
    real(kind=4),allocatable              :: bulge_spectra(:,:)      ! young stars and old stars spectra --> bulge component
    real(kind=4),allocatable              :: disc_spectra(:,:)       ! young stars and old stars spectra --> disc component
    real(kind=4),allocatable              :: tmp_IR_spectrum(:)      !
    real(kind=4),allocatable              :: tmp_agn_eff_ext(:)      ! agn effective extinction (disc ISM + bulge component)

    type(galaxy_type),intent(inout)       :: gal                     ! a galaxy

#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('galaxy_spectrum',only_rank=rank,component='gal')
! -------------------------------------------------
#endif

    !******************************
    ! 1: GALAXY SPECTRA
    !******************************

    ! create galaxy%spectra
    allocate(gal%spectra(nWaves,3))  ! young stars, old stars, full sppectrum
    gal%spectra(:,:) = 0.d0          ! init
    !
    allocate(tmp_IR_spectrum(nWaves))
    tmp_IR_spectrum = 0.d0           ! init
    !
    allocate(tmp_agn_eff_ext(nWaves))
    tmp_agn_eff_ext = 0.d0           ! init
    !
    Mdust = 0.d0
    !
    ! ******************
    ! 1-1: DISC
    ! ******************
    !
    if (stars_mass(gal%disc%stars) .gt. 0.d0) then
        !
        ! ************************
        ! 1-1-1: DISC - OLD STARS
        ! ************************
        !
        ! Build spectrum associated to the disc stellar population
        call stars_build_stellar_spectra(gal%disc%stars,disc_spectra)
        ! build galaxy stellar spectra
        gal%spectra(:,1:2) = gal%spectra(:,1:2) + disc_spectra ! Lsun --> lamb*Flamb
        !
        ! compute average ISRF associated to the component
        ! For old stars in the diffuse ISM we assume :
        !   - an homogeneous distribution of the power
        !   - an average distance between stars and gas/dust D = h/4 (h the disc scale height)
        !   - each side of the disc is illuminated by half of the power
        ! These hypothesis lead to very simple prescription
        ISRF = gal%disc%stars%ISRF(2)/(2.d0*pi*1.1d1*gal%disc%rd)**2./ISRF_ref    ! [G0 unit]
        ! reset
        gal%disc%stars%ISRF(2) = ISRF  ! [G0 unit]
        !
        ! If extinction is turn ON
        if (gal%disc%dust(2)%tau .gt. 0.d0) then
            !
            ! Compute bolometric luminosity of old stars in the disc
            ! use Flamb
            Lbol_before = trap(disc_spectra(:,2)/Waves,Waves)
            !
            ! compute dust extinction
            call dust_compute_effective_extinction(gal%disc%dust(2))
            ! save ISM extinction curve, will be added to bulge component extinction and apply onto un-extinguish AGN SED
            tmp_agn_eff_ext = gal%disc%dust(2)%eff_ext
            !
            ! apply extinction on the old stellar population
            disc_spectra(:,2) = disc_spectra(:,2)*gal%disc%dust(2)%eff_ext
            !
            ! test minimum value
            min_test = minval(disc_spectra(:,2))
            if (min_test .lt. 0.d0) then
                call IO_print_warning_message('Extinction process (disc old): negative value in spectra ',only_rank = rank, called_by = 'galaxy_spectrum')
                call IO_print_message('with',only_rank=rank,component='gal', &
                          param_name = (/'min_test                 '/), real_param_val = (/min_test/))
                where (disc_spectra(:,2) .lt. 0.d0)
                    disc_spectra(:,2) = 0.d0
                endwhere
            end if
            !
            ! build full galaxy spectrum
            gal%spectra(:,3) = gal%spectra(:,3) + disc_spectra(:,2) ! only old stars of the disc, Lsun --> lamb*Flamb
            !
            ! Compute bolometric luminosity of old stars after extinction
            ! use Flamb
            Lbol_after = trap(disc_spectra(:,2)/Waves,Waves)
            ! test energy conservation
            if (Lbol_after .gt. Lbol_before) then
                call IO_print_warning_message('Extinction process (disc old): Lbol_after > Lbol_before ',only_rank = rank, called_by = 'galaxy_spectrum')
                call IO_print_message('with',only_rank=rank,component='gal', &
                          param_name = (/'Lbol_after               ','Lbol_before              '/), real_param_val = (/Lbol_after, Lbol_before/))
            end if
            !
            ! build dust spectrum
            call dust_build_dust_spectrum(gal%disc%dust(2),ISRF,dspt='old stars') ! ISRF old stars
            !
            ! compute A(V) and E(B-V) in the disc ISM
            gal%AV(2)   = dust_AV(gal%disc%dust(2))
            gal%E_BV(2) = dust_E_BV(gal%disc%dust(2))
            !
            ! build full galaxy spectrum
            ! add dust component associated to the old stellar population of the disc, Lsun --> lamb*Flamb
            ! normalised dust spectrum to 1 Lsun and apply luminosity transfer
            ! compute Lir --> energy conservation
            Lir = Lbol_before - Lbol_after
            !
            ! Normalization factor
            Chi = trap(gal%disc%dust(2)%spectrum/Waves,Waves) ! [Lsun/Msun]
            ! deduce Mdust associated to the component
            gal%disc%dust(2)%mass = Lir / Chi * (gal%disc%dust(2)%MZ_MH / MZ_MH_solar)
            Mdust = Mdust + gal%disc%dust(2)%mass
            gal%spectra(:,3) = gal%spectra(:,3) + real((Lir/Chi),4)*gal%disc%dust(2)%spectrum ! Lsun --> lamb*Flamb
        else
            ! No extinction, copy only old star SED
            gal%spectra(:,3) = gal%spectra(:,3) + disc_spectra(:,2)
        end if
        !
        ! **************************
        ! 1-1-2: DISC - YOUNG STARS
        ! **************************
        !
        ! compute average ISRF associated to the component
        ! For young stars in the diffuse ISM we assume :
        !   - young stars distribute in disc%ngc giant molecular clouds
        !   - an average distance between stars and dust D = 3h/8 (h the disc scale height)
        ! These hypothesis lead to:
        ISRF = gal%disc%stars%ISRF(1)/real(gal%disc%ngc,8)/(4.d0*pi*(3.d0*gal%disc%h/8.d0)**2.)/ISRF_ref  ! [G0 unit]
        ! reset
        gal%disc%stars%ISRF(1) = ISRF ! [G0 unit]
        !
        if (gal%disc%dust(1)%tau .gt. 0.d0) then
            !
            ! Compute bolometric luminosity of young stars
            ! use Flamb
            Lbol_before = trap(disc_spectra(:,1)/Waves,Waves)
            !
            ! apply extinction from ISM
            if (gal%disc%dust(2)%tau .gt. 0.d0) then
                disc_spectra(:,1) = disc_spectra(:,1)*gal%disc%dust(2)%eff_ext
                ! test minimum value
                min_test = minval(disc_spectra(:,1))
                if (min_test .lt. 0.d0) then
                    call IO_print_warning_message('Extinction process (disc young ISM): negative value in spectra',only_rank = rank, called_by = 'galaxy_spectrum')
                    call IO_print_message('with',only_rank=rank,component='gal', &
                              param_name = (/'min_test                 '/), real_param_val = (/min_test/))
                    where (disc_spectra(:,1) .lt. 0.d0)
                        disc_spectra(:,1) = 0.d0
                    endwhere
                end if
            end if
            !
            ! apply additional extinction from burst cloud
            ! compute dust extinction
            call dust_compute_effective_extinction(gal%disc%dust(1))
            disc_spectra(:,1) = disc_spectra(:,1)*gal%disc%dust(1)%eff_ext
            ! test minimum value
            min_test = minval(disc_spectra(:,1))
            if (min_test .lt. 0.d0) then
                call IO_print_warning_message('Extinction process (disc young BC): negative value in spectra',only_rank = rank, called_by = 'galaxy_spectrum')
                call IO_print_message('with',only_rank=rank,component='gal', &
                          param_name = (/'min_test                 '/), real_param_val = (/min_test/))
                where (disc_spectra(:,1) .lt. 0.d0)
                    disc_spectra(:,1) = 0.d0
                endwhere
            end if
            !
            ! build full galaxy spectrum
            gal%spectra(:,3) = gal%spectra(:,3) + disc_spectra(:,1) ! only young stars, Lsun --> lamb*Flamb
            !
            ! Compute bolometric luminosity of young stars after complete extinction (ISM + BC)
            ! use Flamb
            Lbol_after = trap(disc_spectra(:,1)/Waves,Waves)
            ! test energy conservation
            if (Lbol_after .gt. Lbol_before) then
                call IO_print_warning_message('Extinction process (disc young): Lbol_after > Lbol_before ',only_rank = rank, called_by = 'galaxy_spectrum')
                call IO_print_message('with',only_rank=rank,component='gal', &
                          param_name = (/'Lbol_after               ','Lbol_before              '/), real_param_val = (/Lbol_after, Lbol_before/))
            end if
            !
            ! build dust spectrum
            call dust_build_dust_spectrum(gal%disc%dust(1),ISRF,dspt='young stars') ! ISRF young stars
            !
            ! compute A(V) and E(B-V) in disc BCs
            if (gal%AV(2) .gt. 0.) then
                gal%AV(1)   = dust_AV(gal%disc%dust(1)) + gal%AV(2)
            else
                gal%AV(1)   = dust_AV(gal%disc%dust(1))
            end if
            if (gal%E_BV(2) .gt. 0.d0) then
                gal%E_BV(1) = dust_E_BV(gal%disc%dust(1)) + gal%E_BV(2)
            else
                gal%E_BV(1) = dust_E_BV(gal%disc%dust(1))
            end if
            !
            ! build full galaxy spectrum
            ! add dust burst cloud component of the disc, Lsun --> lamb*Flamb
            ! normalised dust spectrum to 1 Lsun and apply luminosity transfer
            ! compute Lir --> energy conservation
            Lir = Lbol_before - Lbol_after
            !
            ! Normalization factor
            Chi = trap(gal%disc%dust(1)%spectrum/Waves,Waves) ! [Lsun/Msun]
            ! deduce Mdust associated to the component
            gal%disc%dust(1)%mass = Lir / Chi * (gal%disc%dust(1)%MZ_MH / MZ_MH_solar)
            Mdust = Mdust + gal%disc%dust(1)%mass
            gal%spectra(:,3) = gal%spectra(:,3) + real((Lir/Chi),4)*gal%disc%dust(1)%spectrum ! Lsun --> lamb*Flamb
            !
            ! build the full spectra associated to young stars
            ! young stars spectrum + dust spectrum of BC
            tmp_IR_spectrum = disc_spectra(:,1) + real((Lir/Chi),4)*gal%disc%dust(1)%spectrum
            ! compute Lir [8:1000 mu m] associated to young stars
            gal%Lir(1) = trap(tmp_IR_spectrum(i_8mic:i_1000mic)/Waves(i_8mic:i_1000mic),Waves(i_8mic:i_1000mic)) ! Lsun
            !
        else
            ! No extinction: copy young stars SED
            gal%spectra(:,3) = gal%spectra(:,3) + disc_spectra(:,1)
        end if
    end if
    !
    !******************************
    ! 1-2: DISC MAGNITUDES
    !******************************
    !
    ! @ this point galaxy spectrum contains only disc stellar popualtion and dust
    ! create gal%Mags_disc
    allocate(gal%Mags_disc(nfilters))  ! observer-frame
    ! Init
    gal%Mags_disc = 9.9999999999
    !
    ! compute magnitudes
    if (maxval(gal%spectra(:,3)) .gt. 0.d0) then
        call filters_spec2mag(gal%spectra(:,3),z,gal%Mags_disc,frame='observer')
    end if
    !
    ! ******************
    ! 1-3: BULGE
    ! ******************
    ! In the current model, bulges do not host young stars, they have only old stellar populations
    !
    if (stars_mass(gal%bulge%stars) .gt. 0.d0) then
        !
        ! Build spectrum associated to the bulge stellar population
        call stars_build_stellar_spectra(gal%bulge%stars,bulge_spectra)
        !
        ! build galaxy stellar spectrum
        gal%spectra(:,1:2) = gal%spectra(:,1:2) + bulge_spectra  ! Lsun --> lamb*Flamb
        !
        ! compute average ISRF associated to the component
        ! For old stars in the bulge component we assume :
        !   - half of the power in the half mass radius of the bulge
        !   - an average distance between stars and dust D = r50 (r50 the bulge half mass radius)
        ! These hypothesis lead to
        ISRF = gal%bulge%stars%ISRF(2)/2.d0/(4.d0*pi*((1.d0+sqrt(2.d0))*gal%bulge%rb)**2.)/ISRF_ref  ! [G0 unit]
        ! reset
        gal%bulge%stars%ISRF(2) = ISRF ! [G0 unit]
        !
        if (gal%bulge%dust%tau .gt. 0.d0) then
            !
            ! Compute bolometric luminosity of stars (old) in the bulge
            ! use Flamb
            Lbol_before = trap(bulge_spectra(:,2)/Waves,Waves)
            !
            ! compute dust extinction
            call dust_compute_effective_extinction(gal%bulge%dust)
            ! take into account dust extinction from bulge
            tmp_agn_eff_ext = tmp_agn_eff_ext*gal%bulge%dust%eff_ext
            !
            ! apply extinction
            bulge_spectra(:,2) = bulge_spectra(:,2)*gal%bulge%dust%eff_ext
            ! test minimum value
            min_test = minval(bulge_spectra(:,2))
            if (min_test .lt. 0.d0) then
                call IO_print_warning_message('Extinction process (bulge old): negative value in spectra',only_rank = rank, called_by = 'galaxy_spectrum')
                call IO_print_message('with',only_rank=rank,component='gal', &
                          param_name = (/'min_test                 '/), real_param_val = (/min_test/))
                where (bulge_spectra(:,2) .lt. 0.d0)
                    bulge_spectra(:,2) = 0.d0
                endwhere
            end if
            !
            ! build full galaxy spectrum
            gal%spectra(:,3) = gal%spectra(:,3) + bulge_spectra(:,2) ! only stars (old) of the bulge, Lsun --> lamb*Flamb
            !
            ! Compute bolometric luminosity of stars (old) after extinction
            ! use Flamb
            Lbol_after = trap(bulge_spectra(:,2)/Waves,Waves)
            ! test energy conservation
            if (Lbol_after .gt. Lbol_before) then
                call IO_print_warning_message('Extinction process (bulge old): Lbol_after > Lbol_before ',only_rank = rank, called_by = 'galaxy_spectrum')
                call IO_print_message('with',only_rank=rank,component='gal', &
                          param_name = (/'Lbol_after               ','Lbol_before              '/), real_param_val = (/Lbol_after, Lbol_before/))
            end if
            !
            ! build dust spectrum
            call dust_build_dust_spectrum(gal%bulge%dust,ISRF,dspt='old stars')
            !
            ! build full galaxy spectrum
            ! add dust component of the bulge, Lsun --> lamb*Flamb
            ! normalised dust spectrum to 1 Lsun and apply luminosity transfer
            ! compute Lir --> energy conservation
            Lir = Lbol_before - Lbol_after
            !
            ! Normalization factor
            Chi = trap(gal%bulge%dust%spectrum/Waves,Waves) ! [Lsun/Msun]
            ! deduce Mdust associated to the component
            gal%bulge%dust%mass = Lir / Chi
!~          Mdust = Mdust + gal%bulge%dust%mass
            gal%spectra(:,3) = gal%spectra(:,3) + real((Lir/Chi),4)*gal%bulge%dust%spectrum ! Lsun --> lamb*Flamb
        else
            ! No extinction: copy old stars SED
            gal%spectra(:,3) = gal%spectra(:,3) + bulge_spectra(:,2)
        end if
    end if
    !
    ! ******************
    ! 1-1-3: AGN
    ! ******************
    !
    Lir_AGN = 0.d0  ! init
    !
#ifdef AGN_SPECTRUM
! ------------------------------------------------
    if (agn_mass(gal%disc%agn) .gt. 0.d0) then
        ! compute the sepctrum associated to the AGN
        call agn_spectrum(gal%disc%agn)
        ! Apply extinction from disc ISM and bulge component
        gal%disc%agn%spectrum = gal%disc%agn%spectrum*tmp_agn_eff_ext
        ! compute IR component of the AGN emission
        Lir_AGN = trap(gal%disc%agn%spectrum(i_8mic:i_1000mic)/Waves(i_8mic:i_1000mic),Waves(i_8mic:i_1000mic)) ! Lsun
        ! add to the full gaalxy spectrum
        gal%spectra(:,3) = gal%spectra(:,3) + gal%disc%agn%spectrum
    end if
! ------------------------------------------------
#endif
! AGN_SPECTRUM
    !
    gal%Lir(3) = Lir_AGN
    !
    !******************************
    ! 2 : GALAXY MAGNITUDES
    !******************************
    !
    ! create gal%Mags and gal%Non_Ext_Mags
    allocate(gal%Mags(nfilters,2))          ! observer-frame and rest-frame
    allocate(gal%Non_Ext_Mags(nfilters,2))  ! observer-frame and rest-frame
    ! Init
    gal%Mags = 9.9999999999
    gal%Non_Ext_Mags = 9.9999999999
    ! create and init dMags
    allocate(gal%dMags(nfilters))           ! observer-frame
    gal%dMags = 0.d0
    !
    ! compute magnitudes
    if (maxval(gal%spectra(:,3)) .gt. 0.d0) then
        call filters_spec2mag(gal%spectra(:,3),z,gal%Mags(:,1),frame='observer')
        call filters_spec2mag(gal%spectra(:,3),z,gal%Mags(:,2),frame='galaxy')
        ! precompute first order derivative magnitudes
        call filters_spec2mag(gal%spectra(:,3),z-2.d-1,gal%dMags,frame='observer')
    end if
    ! compute non-extinced magnitudes, used only young stars and old stars
    if (maxval(gal%spectra(:,1)+gal%spectra(:,2)) .gt. 0.d0) then
        call filters_spec2mag(gal%spectra(:,1)+gal%spectra(:,2),z,gal%Non_Ext_Mags(:,1),frame='observer')
        call filters_spec2mag(gal%spectra(:,1)+gal%spectra(:,2),z,gal%Non_Ext_Mags(:,2),frame='galaxy')
    end if
    !
    ! compute first order derivative magnitudes
    where ((gal%Mags(:,1) < 9.d2) .and. (abs(gal%dMags) > 0.d0))
        gal%dMags = (gal%dMags-gal%Mags(:,1))/(-0.2)
    elsewhere
        gal%dMags = 0.
    end where

    !***********************************
    ! 3 : ADDITIONNAL PROPERTIES
    !***********************************

    if (maxval(gal%spectra(:,3)) .gt. 0.d0) then
        !
        ! Total mass of dust
        gal%Mdust = Mdust
        !
        ! bolometric luminosity of the whole stellar population
        gal%Lbol_stars = trap((gal%spectra(:,1)+gal%spectra(:,2))/Waves,Waves)   ! Lsun
        !
        ! ISRF young stars
        if (gal%disc%stars%ISRF(1) .gt. 0.d0) then
            gal%ISRF(1) = gal%disc%stars%ISRF(1)  ! G0 unit
        end if
        ! ISRF old stars
        if (gal%disc%stars%ISRF(2) .gt. 0.d0) then
            gal%ISRF(2) = gal%disc%stars%ISRF(2)  ! G0 unit
        end if
        !
        ! Infrared luminosity [8 : 1000] mu m "full"
        gal%Lir(2) = trap(gal%spectra(i_8mic:i_1000mic,3)/Waves(i_8mic:i_1000mic),Waves(i_8mic:i_1000mic)) ! Lsun
        !
        ! test Lir
        if (gal%Lir(2)-Lir_AGN .gt. gal%Lbol_stars) then
            spec_min = minval(gal%spectra(i_8mic:i_1000mic,3))
            spec_max = maxval(gal%spectra(i_8mic:i_1000mic,3))
            call IO_print_error_message('Too much IR emission !',only_rank = rank, called_by = 'galaxy_spectrum')
            call IO_print_message('with',only_rank=rank,component='gal', &
                      param_name = (/'Lbol_stars               ','Lir "full"               ', &
                                     'spec_min                 ','spec_max                 '/), &
                                     real_param_val = (/real(gal%Lbol_stars,8), real(gal%Lir(2),8),spec_min,spec_max/))
            stop ! stop the program
        end if
    end if

    ! deallocate tmp array
    if (allocated(tmp_IR_spectrum)) deallocate(tmp_IR_spectrum)
    if (allocated(tmp_agn_eff_ext)) deallocate(tmp_agn_eff_ext)

#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('galaxy_spectrum ... done',only_rank=rank,component='gal')
! -------------------------------------------------
#endif

    return
  end subroutine galaxy_spectrum

! -------------------------------------------------
#endif
! LUMINOUS_PROCESSES

  !*****************************************************************************************************************
  !
  ! FUNCTIONS
  !
  !*****************************************************************************************************************

  function galaxy_mass(gal,r,component)

    ! RETURN THE TOTAL MASS OF THE GALAXY
    ! component is constistent with option given in disc or bulge component
    ! if r is provided, return mass of gal within r

    implicit none

    character(*),intent(in),optional  :: component    ! allows to select a component (gas, stars ...)

    real(kind=8),intent(in),optional  :: r            ! the radius
    real(kind=8)                      :: galaxy_mass  ! the mass of tha galaxy

    type(galaxy_type),intent(in)      :: gal          ! a galaxy component


    galaxy_mass = disc_mass(gal%disc) + bulge_mass(gal%bulge)

    if (galaxy_mass .le. 0.d0) return ! no galaxy

    if (present(r)) then
      if (max(gal%disc%rd, gal%bulge%rb) .le. 0.d0) then
        call IO_print_error_message('r = 0 with massive galaxy', &
                only_rank = rank, called_by = 'galaxy_mass')
        stop ! stop the program
      end if
      galaxy_mass = disc_mass(gal%disc,r=r,component=component) + bulge_mass(gal%bulge,r=r,component=component)
    else
      galaxy_mass = disc_mass(gal%disc,component=component) + bulge_mass(gal%bulge,component=component)
    end if

    return
  end function galaxy_mass

  !*****************************************************************************************************************

  function galaxy_frac_mass_radius(gal,frac,component,called_by)

    ! RETURN THE RADIUS WHICH ENCLOSE 100*frac % OF THE GALAXY MASS
    ! component is constistent with option given in disc or bulge component

    implicit none

    character(*),intent(in),optional  :: component                 ! allow to select the gaalxy component (disc or bulge)
    character(*),intent(in),optional  :: called_by                 ! allow to known the calling function (or subroutine)
    character(MAXPATHSIZE)            :: message                   ! a message to display

    real(kind=8),intent(in)           :: frac                      ! fraction of the galaxy mass
    real(kind=8)                      :: galaxy_frac_mass_radius   ! the radius which enclose the gioven fraction of the galaxy mass
    real(kind=8)                      :: mass, r                   ! local mass and radius
    real(kind=8)                      :: target_mass               !
    real(kind=8)                      :: r_min, r_max              ! interval in which a solution for galaxy_half_mass_radius is searched

    type(galaxy_type),intent(in)      :: gal                       ! a galaxy component

    if (present(called_by)) write(message,'(a,a)') 'Function <<galaxy_frac_mass_radius>> is called by: ', trim(called_by)

    galaxy_frac_mass_radius = 0.d0
    !
    ! init the target mass
    target_mass = frac*galaxy_mass(gal,component=component)
    if (target_mass .le. 0.d0) return
    !
    r_min = 0.                                   ! init the minima radius
    r_max = 1.d2*max(gal%disc%rd,gal%bulge%rb)   ! init the maximal radius
    !
    if (r_max .le. 0.d0) then
      call IO_print_error_message('r_max <= 0', &
            only_rank = rank, called_by = 'galaxy_frac_mass_radius')
      if (present(called_by)) call IO_print_message(trim(message),only_rank=rank,component='gal')
      call IO_print_message('used',only_rank=rank,component='gal',&
        param_name=(/'disc%rd                  ','bulge%rb                 ','gal%mass                 '/), &
        real_param_val =(/gal%disc%rd,gal%bulge%rb,galaxy_mass(gal)/))
      stop ! stop the program
    end if
    !
    ! test with r_max and crash the code if the optimal radius is larger than this barrier.
    mass        = galaxy_mass(gal,r=r_max,component=component)
    !
    if (mass .lt. target_mass) then
      call IO_print_error_message('Optimal r > r_max', &
            only_rank = rank, called_by = 'galaxy_frac_mass_radius')
      if (present(called_by)) call IO_print_message(trim(message),only_rank=rank,component='gal')
      call IO_print_message('used',only_rank=rank,component='gal',&
        param_name=(/'disc%rd                  ','bulge%rb                 ','gal%mass                 '/), &
        real_param_val =(/gal%disc%rd,gal%bulge%rb,galaxy_mass(gal)/))
      stop ! stop the program
    end if

    r = -1.0d0
    do while ((abs(mass-target_mass)/target_mass .gt. num_precision) &
      .and. (r_min .ne. r_max))
       r = 5.d-1*(r_min+r_max)
       mass = galaxy_mass(gal,r=r,component=component)
       ! the galaxy mass is an increasing function of the galactic radius
       if (mass .le. target_mass) then
          r_min = r
       else
          r_max = r
       end if
    end do

    ! test the result and crash the code if the optimal radius is smaller than 0.
    if (r .gt. 0.d0) then
       galaxy_frac_mass_radius = r
    else
      call IO_print_error_message('Optimal r < 0', &
            only_rank = rank, called_by = 'galaxy_frac_mass_radius')
      if (present(called_by)) call IO_print_message(trim(message),only_rank=rank,component='gal')
        call IO_print_message('used',only_rank=rank,component='gal',&
            param_name=(/'disc%rd                  ','bulge%rb                 ','gal%mass                 '/), &
            real_param_val =(/gal%disc%rd,gal%bulge%rb,galaxy_mass(gal)/))
      stop ! stop the program
    end if

    return
  end function galaxy_frac_mass_radius

  !*****************************************************************************************************************

  function galaxy_gas_signature(gal,component,subcomponent,apply_as,called_by)

    ! RETURN A GAS SIGNATURE (A NORMALIZED GAS OBJECT)

    implicit none

    character(*),intent(in),optional  :: component             ! allow to select the galaxy component (disc, galaxy)
    character(*),intent(in),optional  :: subcomponent          ! allow to select the subcomponent (sfg, non-sfg)
    character(*),intent(in),optional  :: called_by             ! name of the function which has called this function
    character(*),intent(in),optional  :: apply_as              ! the gas signature function can be use in different cases
                                                               ! when this function is used as an output rate builder, we have to check critical mass
    character(MAXPATHSIZE)            :: message               ! a message to display

    type(gas_type)                    :: galaxy_gas_signature  ! the gas signature

    type(galaxy_type),intent(in)      :: gal                   ! a galaxy component

    if (present(component)) then
        select case (trim(component))
        case ('disc')
            galaxy_gas_signature = disc_gas_signature(gal%disc,component=subcomponent,apply_as=apply_as,called_by=called_by)
        case ('bulge')
            galaxy_gas_signature = bulge_gas_signature(gal%bulge,component=subcomponent)
        case ('galaxy')
            galaxy_gas_signature = disc_mass(gal%disc,component='gas')*disc_gas_signature(gal%disc,component='gas',apply_as=apply_as,called_by=called_by)
            galaxy_gas_signature = galaxy_gas_signature + bulge_mass(gal%bulge,component='gas')*bulge_gas_signature(gal%bulge,component='gas')
            galaxy_gas_signature = gas_signature(galaxy_gas_signature)
        case default
            write(message,'(a,a,a)') 'Keyword ', trim(component), ' not defined'
            call IO_print_error_message(message,only_rank=rank,called_by='galaxy_gas_signature')
            stop  ! stop the program
        end select
    else
        ! use all gas as defaut value
        galaxy_gas_signature = disc_mass(gal%disc,component='gas')*disc_gas_signature(gal%disc,component='gas',apply_as=apply_as,called_by=called_by)
        galaxy_gas_signature = galaxy_gas_signature + bulge_mass(gal%bulge,component='gas')*bulge_gas_signature(gal%bulge,component='gas')
        galaxy_gas_signature = gas_signature(galaxy_gas_signature)
    end if

    return
  end function galaxy_gas_signature

  !*****************************************************************************************************************

  function galaxy_compute_escape_velocity(gal,dm)

    ! RETURNS THE ESCAPE VELOCITY OF THE GALAXY (in code unit)
    ! This function have to be called after the computation of gal%R50

    implicit none

    real(kind=8)                    :: galaxy_compute_escape_velocity
    real(kind=8)                    :: N_clumps

    type(galaxy_type),intent(in)    :: gal   ! the galaxy component
    type(dm_type),intent(in)        :: dm    ! the dark matter component

    galaxy_compute_escape_velocity = -1.d0   ! init

    if (galaxy_mass(gal) .le. 0.d0) return   ! no galaxy

    if (gal%R50 .gt. 0.d0) then
        N_clumps = min(3.d1,max(2.d0,real(gal%disc%ngc,4)))
        galaxy_compute_escape_velocity = sqrt(2.d0*gravconst_code_unit*( &
                            disc_mass(gal%disc,r=gal%R50,component='stars') + &
                            disc_mass(gal%disc,r=gal%R50,component='gas')/N_clumps + &
                            bulge_mass(gal%bulge,r=gal%R50) &
                            + dm_mass(dm,r=gal%R50))/gal%R50)
    end if

    return
  end function

  !*****************************************************************************************************************
  !
  ! PRINTING PROCEDURES
  !
  !*****************************************************************************************************************

  subroutine galaxy_load_gal_data(gal,gal_data,data_type)

    ! CREATE THE GAL OUTPUT PROPERTIES LIST

    implicit none

    character(*),intent(in)                :: data_type              ! phy or lum data

    real(kind=8),allocatable,intent(inout) :: gal_data(:)
    real(kind=8)                           :: Msd,Msb,dbr

    type(galaxy_type),intent(in)           :: gal   ! gal component

    select case (trim(data_type))
    case('phy','physical','phy_data')
        !
        ! PHYSICAL PROPERTIES
        !
        ! For information
        ! 'nb_merger             ','nb_major_merger       ','tilt                  ',&
        ! 'accreted_mass         ','merger_mass           ','m_stars               ',&
        ! 'mu_tot                ','mu_gas                ','R50_stars             ',&
        ! 'disc_bulge_ratio      ','gal_stars_age_MW      ','gal_stars_dage_MW     ',&
        ! 'gal_stars_age_LW      ','gal_stars_dage_LW     ','gal_Z_stars_MW        ',&
        ! 'gal_Z_stars_LW        ','Vesc                  '
        !
        ! allocate gal_data array
        allocate(gal_data(nb_gal_field))
        !
        ! compute the disc to bulge stellar mass ratio
        Msd = disc_mass(gal%disc,component='stars')
        Msb = bulge_mass(gal%bulge,component='stars')
        dbr = -1.d0 ! init
        if (Msd + Msb .gt. 0.d0) dbr = Msd/(Msd+Msb)
        !
        gal_data = (/real(gal%nb_merger,8),real(gal%nb_major_merger,8),gal%disc%incl, &
                     mass_code_unit_in_M_Sun*gal%accreted_gas_mass, mass_code_unit_in_M_Sun*gal%merger_gas_mass , &
                     mass_code_unit_in_M_Sun*galaxy_mass(gal,component='stars'),gal%mu_tot, gal%mu_gas, gal%R50_stars, dbr, &
                     gal%Stel_Age(1),gal%Stel_dAge(1),gal%Stel_Age(2),gal%Stel_dAge(2),gal%Stel_Z(1)/Z_sun,gal%Stel_Z(2)/Z_sun, &
                     gal%Vesc/)
        !
    case('lum','luminous','lum_data')
        !
        ! LUMINOUS PROPERTIES
        !
        ! For information
        ! non-extincted magnitudes for each input filter (obs-frame and rest-frame)
        ! extincted magnitudes for each input filter (obs-frame)
        ! extincted magnitudes associated to the disc, for each input filter (obs-frame only)
        ! extincted magnitudes for each input filter (rest-frame)
        ! first order derivative magnitudes for each input filter (obs-frame only)
        !   + additionnal luminosity poperties Lbol_stars, Lir ...
        !
        ! allocate gal_data array
        ! set nb_Mag_field
        allocate(gal_data(6*nfilters + nlumprops))
        gal_data = 0.d0
        !
        ! non-extincted Magnitudes
        gal_data(1:nfilters)              = gal%Non_Ext_Mags(:,1)  ! observer-frame
        gal_data(nfilters+1:2*nfilters)   = gal%Non_Ext_Mags(:,2)  ! galaxy-frame
        ! Magnitudes
        gal_data(2*nfilters+1:3*nfilters) = gal%Mags(:,1)          ! observer-frame
        gal_data(3*nfilters+1:4*nfilters) = gal%Mags_disc(:)       ! disc only observer-frame
        gal_data(4*nfilters+1:5*nfilters) = gal%Mags(:,2)          ! galaxy-frame
        ! dMags
        gal_data(5*nfilters+1:6*nfilters) = gal%dMags              ! observer-frame
        !
        ! add additionnal luminous properties (Lir, Lbol ...)
        ! for information
        !'Lbol_stars            ','ISRF_disc_ISM         ','ISRF_BC               ', &
        !'Lir                   ','Lir_young_stars       ','Lir_AGN               ',&
        !'Mdusct                ','AV_disc_ISM           ','AV_BC                 ', &
        !'E_BV_disc_ISM         ','E_BV_BC
        gal_data(6*nfilters+1:) = (/gal%Lbol_stars,gal%ISRF(2),gal%ISRF(1),gal%Lir(2),gal%Lir(1),gal%Lir(3),gal%Mdust, &
                                    gal%AV(2),gal%AV(1),gal%E_BV(2),gal%E_BV(1)/)
        !
    case default
        call IO_print_error_message('Unknwon data type', &
                only_rank = rank, called_by = 'galaxy_load_gal_data')
        stop ! stop the program
    end select

    return
  end subroutine galaxy_load_gal_data

  !*****************************************************************************************************************

  subroutine galaxy_print(unit,form,data_type,gal,go_down)

    ! PRINT GALAXY-PROPERTIES IN eGALICS OUTPUT FILE

    implicit none

    integer(kind=4),intent(in)   :: unit              ! file unit
    integer(kind=4)              :: status,i,hdutype  ! fits write process information

    logical,optional              :: go_down          ! = .true. if higher level printing procedures have to be called
                                                      ! only in merger tree printing procedure

    character(*)                 :: form              ! fits or tmp_bin
    character(*),intent(in)      :: data_type         ! phy or lum data

    real(kind=8),allocatable     :: gal_data(:)       !

    type(galaxy_type),intent(in) :: gal               ! the galaxy component

    select case (trim(data_type))
    case('phy','physical','phy_data')
        !
        ! PHYSICAL PROPERTIES
        !
        call galaxy_load_gal_data(gal,gal_data,data_type)
        ! select the output format
        select case (trim(form))
        case ('tmp_bin')
            write(unit) gal_data  ! directly write data in the tmp binary output file
        case ('fits')
            ! move to gal extension
            call ftmahd(unit,hdu_gal,hdutype,status)
            if (status .gt. 0) then
                call IO_print_error_message('ftmahd status', &
                    only_rank = rank, called_by = 'gal_print')
                stop ! stop the program
            end if
            ! init
            call ftirow(unit,0,1,status)
            if (status .gt. 0) then
                call IO_print_error_message('ftirow status', &
                    only_rank = rank, called_by = 'gal_print')
                stop ! stop the program
            end if
            ! write data in the dm entension
            call ftpcld(unit,1,1,1,1,gal%age_form+gal%life_time,status)
            do i=2, nb_gal_field+1
                call ftpcld(unit,i,1,1,1,gal_data(i-1),status)
                if (status .gt. 0) then
                    call IO_print_error_message('ftpcld status', &
                        only_rank = rank, called_by = 'gal_print')
                    stop ! stop the program
                end if
            end do
        case default
            call IO_print_error_message('Unknwon output data format', &
                only_rank = rank, called_by = 'gal_print')
            stop ! stop the program
        end select
        !
        if (present(go_down)) then
            if (go_down) then
                ! write disc and bulge properties
                ! disc_print(unit,form,disc)
                call disc_print(unit,form,gal%disc)   ! print disc properties
                ! bulge_print(unit,form,bulge)
                call bulge_print(unit,form,gal%bulge) ! print bulge properties
            end if
        end if
        !
    case('lum','luminous','lum_data')
        !
        ! LUMINOUS PROPERTIES
        !
        call galaxy_load_gal_data(gal,gal_data,data_type)
        !
        ! For information
        ! non-extincted magnitudes for each input filter (obs-frame and rest-frame)
        ! extincted magnitudes for each input filter(obs-frame and rest-frame)
        ! first order derivative magnitudes for each input filter (observer frame only)
        !   + additionnal luminosity poperties Lbol_stars, Lir ...
        write(unit) real(gal_data,4)
#ifdef GAL_SPECTRA
! -------------------------------------------------
        ! write the galaxy spectrum associated to young stars
        write(unit) gal%spectra(:,1)
        ! write the galaxy spectrum associated to old stars
        write(unit) gal%spectra(:,2)
        ! write the full galaxy spectrum
        write(unit) gal%spectra(:,3)
! -------------------------------------------------
#endif
! GAL_SPECTRA

    case default
        call IO_print_error_message('Unknwon data type', &
                only_rank = rank, called_by = 'galaxy_print')
        stop ! stop the program
    end select

    return
  end subroutine galaxy_print

  !*****************************************************************************************************************

end module galaxy
