module disc

  use agn      ! Contains agn structure definition and agn exploitation functions
  use bulge    ! Contains bulge structure definition and bulge exploitation functions

  public

  !*****************************************************************************************************************
  !
  !  OVERVIEW
  !
  !  disc_module defines the disc data structure univ(ts)%halo(ih)%galaxy%disc
  !
  !  This module contains the definition of the disc structure and all functions and subroutines associated.
  !  In the header of the module are defined output properties of the disc component (labels, units and formats)
  !  The disc is the active part of the galaxy. It accreates diffuse gas. This diffuse gas is then progressivelly converted
  !  into dense/fragmented gas and then into star-forming gas
  !  New stars are only formed into the disc. Fresh gas accretion are only supported by the disc component
  !  SN that explose into the disc, allows to head the gas, eject it but allows feed turbulent motion in the gas
  !  Gas fragmentation process is mainly followed according to the diffuse gas condensation rate and the scale-height of the disc
  !  The disc evolution procedure is divided in two routines, disc_evolve_I and disc_evolve_II
  !
  !  MPI PROCESSES IN THIS MODULE
  !
  !  disc_send_data                                  : send specific informations about a disc
  !
  !  disc_receive_data                               : receive specific informations about a disc
  !
  !  SUBROUTINES IN THIS MODULE
  !
  !   disc_void                                      : init all properties of the disc structure
  !
  !   disc_copy                                      : copy a disc component to an other one
  !
  !   disc_evolve_I                                  : first part (PREDICTOR) of the disc evolution scheme
  !      called by : galaxy_evolve                     compute input, output and evolution rate associated to the gas structuration history
  !      contains  : disrupt_str_gas_rate              return the instantaneous disruption rate of the structured gas phase
  !                  gas_str_rate                      return the instantaneous strucuration rate of teh diffuse gas phase
  !
  !   disc_evolve_II                                 : second part (CORRECTOR) of the disc evolution scheme
  !       called by : galaxy_evolve                    apply to the gas structuration table (gsh_tab) evolution rates
  !                                                    computed in disc_evolve_I
  !
  !   disc_compute_disc_feedback_activities          : compute global effets of the feedback into the disc structure (SN + AGN)
  !      called by : disc_evolve_I
  !      contained : SN_ejecta_rate                  : Return the gas ejection rate associated to the SN kinetic energy injection
  !                  SN_wind_velocity                : Return the velocity of wind produced by SN kinetic energy injection
  !                  SN_non_kinetic_power            : return the non-kenitic power produced by SN
  !                  SN_thermal_power                : return the thermal power produced by SN
  !                  SN_rad_power                    : return the radiative power produced by SN (residual luminosity)
  !                  SN_turbulent_heating_power      : return the SN power injected in the gas turbulnet motion
  !
  !   disc_update_gas_struct_history                 : update the gas structuration history table associated to a given disc
  !      called by : disc_evolve_II
  !                : galaxy_merge
  !
  !   disc_compute_inclination                       : compute the inclination of the disc, and the angular momentum
  !      called by : disc_evolve_II
  !                : galaxy_merge
  !
  !   disc_deallocate                                : deallocate all allocated array in a disc structure
  !
  !   disc_update_cooling_clock                      : follow the effective cooling time of the diffuse gas contained into the disc
  !
  !  FUNCTIONS IN THIS MODULE
  !
  !   disc_SFR                                       : Return the star formation rate associated to a given disc
  !
  !   disc_mass_surf_density_                        : return the disc mass surface density
  !
  !   disc_mass                                      : return the disc mass
  !
  !   disc_gas_signature                             : return the signature of a given gas component
  !
  !   disc_exponential_radius                        : return the exponential radius of a disc component
  !
  !   disc_scale_height                              : return the disc scale height of a disc component
  !
  !   disc_escape_velocity                           : return the disc escape velocity
  !
  !   disc_mass_weighted_orbital_velocity            : return the mass-weighted orbital velocity of the disc
  !
  !   disc_dynamical_time                            : return the dynamical time of the disc
  !
  !   disc_star_formation_timescale                  : return the star formation timescale associated to a given disc
  !
  !   disc_emptying_timescale                        : return the time need to completly void the gas structured/fragmented
  !                                                         gas reservoir according to the instantaneous fragmentation timescale
  !
  !   disc_gas_fraction                              : return the gas fraction (in mass) of a given gas disc component
  !
  !   isolated_disc_velocity (INTERFACE VERSION)     : return the rotationnal velocity of a isolated disc (without bulge and dm structures)
  !
  !   isolated_disc_velocity_ (CORE VERSION)
  !
  !   disc_velocity (INTERFACE VERSION)              : return the total rotationnal velocity of a disc
  !
  !   disc_velocity_ (CORE VERSION)
  !
  !   disc_velocity_dispersion                       : return the velocity dispersion of a disc component
  !
  !   disc_kappa (INTERFACE VERSION)                 : return the epicyclic frequency (at a given radius) of a disc
  !
  !   disc_kappa_ (CORE VERSION)
  !
  !   disc_Toomre_parameter                          : return the overall disc Toomre parameter associated to the diffuse gas phase
  !
  !   kappa_sigma_r                                  : return the product of the epicyclic frequency, the mass surface density and the orbital radius
  !                                                    this function allows to compute the mass-weighted epiclyclic frequency
  !   disc_mass_weighted_kappa                       : return the mass-weighted epicyclic frequency
  !
  !   V_sigma_r                                      : return the product of the orbital velocity, the mass surface density and the orbital radius
  !                                                    the function allows to compute the mass-weighted rotational velocity of the disc
  !
  !  PRINTING PROCEDURES
  !
  !   disc_load_disc_data
  !      called by  : disc_print                     : create the output property list of the disc component
  !
  !   disc_print                                     : print disc-properties in output files
  !       called by : disc_evolve_II
  !
  !*****************************************************************************************************************

  type disc_type
    ! Evolution
    real(kind=8)        :: life_time           ! the time life of this disc component
    real(kind=8)        :: age_form            ! age of the universe when the disc has been formed
    real(kind=8)        :: t_since_last_merger ! time since the last merger event between two discs
    !
    ! morphology
    character(6)        :: morpho              ! relaxed disc or clumpy
    !
    ! main components
    type(gsh_type)      :: gsh_tab             ! the gas structuration history table [m_gas,in_rate,out_rate]
    type(stars_type)    :: stars               ! a stellar population
    type(dust_type)     :: dust(2)             ! dust, associated to 1) burst cloud, 2) diffuse ISM
    type(agn_type)      :: agn                 ! A SMBH and this gas torus
    !
    ! disc properties
    integer(kind=4)     :: nlevels             ! number of gas structuration levels (-1 if no gas structuration)
    integer(kind=4)     :: ngc                 ! number of giant clouds
    real(kind=8)        :: f_str               ! fraction of structured/fragmented gas in the disc
    real(kind=8)        :: rd                  ! the disc has an exponential profile and its size is described by the exp radius, in [kpc]
    real(kind=8)        :: h                   ! scale height of the disc, in [kpc]
    real(kind=8)        :: incl                ! inclination of the disc
    real(kind=8)        :: L(3)                ! angular momentum (Lx, Ly, Lz)
    real(kind=8)        :: kappa               ! epicyclic frequency
    real(kind=8)        :: V                   ! orbital velocity, computed ar 2.2rd (Pelliccia+17)
    real(kind=8)        :: t_dyn               ! dynamical time of the disc (2.2rd/V) in [Gyr]
    real(kind=8)        :: t_cool              ! cooling clock of the unstructured gas
    real(kind=8)        :: cooling_timescale   ! cooling timescale of the diffuse gas
    real(kind=8)        :: dV                  ! mean velocity dispersion at the disc scale height (in the diffuse gas)
    real(kind=8)        :: Q                   ! Toomre parameter
    real(kind=8)        :: sfr_burst           ! save value of the instantaneous sfr at the last merger event
    real(kind=8)        :: Qturb_unstr         ! the instantaneous turbulent power injected in the diffuse gas phase
    !
    ! interaction with the main halo hot gas (stripping)
    type(gas_type)      :: stripping_rate      ! stripping rate undergoes by the diffuse gas phase of the disc
    !
    ! transfer rates into the disc structure
    type(gas_type)      :: fresh_gas_acc_rate  ! the global accretion rate onto the disc (cold-streams + cooling)
    type(gas_type)      :: sfr                 ! star formation rate
    type(gas_type)      :: str_rate            ! gas structuration rate, transfer from the diffuse to the fragmented non star forming gas phases
    type(gas_type)      :: disrupt_rate        ! gas disruption rate, transfer from the fragmented non star formaing and star-forming gas to the diffuse gas
    type(gas_type)      :: ejecta_rate         ! ejecta rate of the disc (large scale wind generated from SN + AGN)
  end type disc_type

  real(kind=8)          :: disc_dt_min_use                          ! minimal time-step used at the disc scale

  ! hdu reference for disc structure
  integer(kind=4)           :: hdu_disc
  ! printable properties for disc structure
  integer(kind=4),parameter :: nb_disc_field = 38 ! Number of disc properties saved
  ! Name of each output colomn data
  character(len=ttype_len),dimension(nb_disc_field) :: ttype_disc = (/'t_since_last_merger   ','disc_unstr_gas        ','disc_unstr_mZ         ',&
                                                                      'disc_unstr_mH         ','disc_unstr_mC         ','disc_unstr_mN         ',&
                                                                      'disc_unstr_mO         ','disc_unstr_mFe        ','disc_unstr_OH_index   ',&
                                                                      'disc_str_gas          ','disc_str_mZ           ','disc_str_mH           ',&
                                                                      'disc_str_mC           ','disc_str_mN           ','disc_str_mO           ',&
                                                                      'disc_str_mFe          ','disc_str_OH_index     ','disc_f_str            ',&
                                                                      'disc_ngc              ','disc_rd               ','disc_h                ',&
                                                                      'disc_V                ','disc_dV               ','disc_Q                ',&
                                                                      'disc_cool_timescale   ','disc_t_dyn            ','disc_t_cool           ',&
                                                                      'disc_t_cascade        ','disc_t_form           ','disc_t_eq             ',&
                                                                      'disc_t_emp            ','disc_gas_acc_rate     ','disc_stripping_rate   ',&
                                                                      'disc_str_rate         ','disc_disrupt_rate     ','disc_ejecta_rate      ',&
                                                                      'disc_sfr              ','disc_sfr_burst        '/)
  ! Physical unit of each column data
  character(len=tunit_len),dimension(nb_disc_field) :: tunit_disc = (/'Gyr         ','M_sun       ','M_sun       ','M_sun       ','Msun        ','M_sun       ',&
                                                                      'M_sun       ','M_sun       ','w_o_unit    ','M_sun       ','Msun        ','M_sun       ',&
                                                                      'M_sun       ','M_sun       ','M_sun       ','M_sun       ','w_o_unit    ','w_o_unit    ',&
                                                                      'w_o_unit    ','kpc         ','kpc         ','km/s        ','km/s        ','wounit      ',&
                                                                      'Gyr         ','Gyr         ','Gyr         ','Gyr         ','Gyr         ','Gyr         ',&
                                                                      'Gyr         ','M_sun/yr    ','M_sun/yr    ','M_sun/yr    ','M_sun/yr    ','M_sun/yr    ',&
                                                                      'M_sun/yr    ','M_sun/yr    '/)
  ! Data type of each column data
  character(len=tform_len),dimension(nb_disc_field) :: tform_disc = (/'1E','1E','1E','1E','1E','1E','1E','1E','1E','1E','1E','1E',&
                                                                      '1E','1E','1E','1E','1E','1E','1J','1E','1E','1E','1E','1E','1E',&
                                                                      '1E','1E','1E','1E','1E','1E','1E','1E','1E','1E','1E','1E','1E'/)
  ! *** STARS ***
  ! Name of each output colomn data
  character(len=ttype_len),dimension(nb_stars_field) :: ttype_starsd = (/'disc_stars_mass       ','disc_stars_f_young    ', &
                                                                         'disc_stars_age_MW     ','disc_stars_dage_MW    ', &
                                                                         'disc_stars_age_LW     ','disc_stars_dage_LW    ', &
                                                                         'disc_stars_Z_MW       ','disc_stars_Z_LW       '/)
  ! Physical unit of each column data
  character(len=tunit_len),dimension(nb_stars_field) :: tunit_starsd = (/'M_sun       ','wounit      ','Gyr         ','Gyr         ', &
                                                                         'Gyr         ','Gyr         ','Z_sun       ','Z_sun       '/)
  ! Data type of each column data
  character(len=tform_len),dimension(nb_stars_field) :: tform_starsd = (/'1E','1E','1E','1E','1E','1E','1E','1E'/)

  ! *** DUSTS ***
  ! Name of each output colomn data
  character(len=ttype_len),dimension(nb_dust_field) :: ttype_dust_ISMd = (/'disc_dust_ISM_f_pah   ','disc_dust_ISM_f_bg    ', &
                                                                           'disc_dust_ISM_tau     ','disc_dust_ISM_mZmH    '/)
  ! Name of each output colomn data
  character(len=ttype_len),dimension(nb_dust_field) :: ttype_dust_BCd  = (/'disc_dust_BC_f_pah    ','disc_dust_BC_f_bg     ', &
                                                                           'disc_dust_BC_tau      ','disc_dust_BC_mZmH     '/)

  ! Physical unit of each column data
  character(len=tunit_len),dimension(nb_dust_field) :: tunit_dustd = (/'w_o_unit    ','w_o_unit    ','w_o_unit    ','w_o_unit    '/)

  ! Data type of each column data
  character(len=tform_len),dimension(nb_dust_field) :: tform_dustd = (/'1E','1E','1E','1E'/)

contains

  !*****************************************************************************************************************
  !
  ! MPI PROCESSES
  !
  !*****************************************************************************************************************

#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------
  subroutine disc_send_data(disc)

    ! SEND SPECIFIC INFORMATIONS ABOUT A DISC FROM ONE PROCESS TO AN OTHER

    implicit none

    logical                      :: go_down

    type(disc_type),intent(in)   :: disc

    ! data are sent by a physical process and receive by a luminous process

    go_down = .false.

    if (disc_mass(disc,component='stars') .gt. 0.d0) then
        !
        ! the disc have a stellar population
        go_down = .true.
        call MPI_SEND(go_down,1,MPI_LOGICAL,rank+1,disc_tag+1,MPI_COMM_WORLD,ierror)
        ! send disc exponential radius
        call MPI_SEND(disc%rd,1,MPI_REAL8,rank+1,disc_tag+2,MPI_COMM_WORLD,ierror)
        ! send disc scale height
        call MPI_SEND(disc%h,1,MPI_REAL8,rank+1,disc_tag+3,MPI_COMM_WORLD,ierror)
        ! send number of giant molecular clouds
        call MPI_SEND(disc%ngc,1,MPI_INTEGER4,rank+1,disc_tag+4,MPI_COMM_WORLD,ierror)
        ! send stellar population of the disc
        call stars_send_data(disc%stars)
        ! send disc's dust component for disc BC extinction
        call dust_send_data(disc%dust(1))
        ! send disc's dust component for disc ISM extinction
        call dust_send_data(disc%dust(2))
        ! send agn component
        call agn_send_data(disc%agn)
    else
        !
        ! send exit loop order
        call MPI_SEND(go_down,1,MPI_LOGICAL,rank+1,disc_tag+1,MPI_COMM_WORLD,ierror)
    end if

    return
  end subroutine disc_send_data

  !*****************************************************************************************************************

  subroutine disc_receive_data(disc)

    ! RECEIVE SPECIFIC INFORMATIONS ABOUT A DISC

    implicit none

    logical                      :: go_down

    type(disc_type),intent(out)  :: disc

    ! data are sent by a physical process and receive by a luminous process

    call disc_void(disc)  ! init the disc

    ! receive exit loop order
    call MPI_RECV(go_down,1,MPI_LOGICAL,rank-1,disc_tag+1,MPI_COMM_WORLD,statut,ierror)

    if (go_down) then
        !
        ! receive disc exponential radius
        call MPI_RECV(disc%rd,1,MPI_REAL8,rank-1,disc_tag+2,MPI_COMM_WORLD,statut,ierror)
        ! receive disc scale height
        call MPI_RECV(disc%h,1,MPI_REAL8,rank-1,disc_tag+3,MPI_COMM_WORLD,statut,ierror)
        ! receice number of giant molecular clouds
        call MPI_RECV(disc%ngc,1,MPI_INTEGER4,rank-1,disc_tag+4,MPI_COMM_WORLD,statut,ierror)
        ! receive stellar population of the disc
        call stars_receive_data(disc%stars)
        ! receive disc's dust component, BC
        call dust_receive_data(disc%dust(1))
        ! receive disc's dust component, ISM
        call dust_receive_data(disc%dust(2))
        ! receive agn's component
        call agn_receive_data(disc%agn)
    end if

    return
  end subroutine disc_receive_data
! -------------------------------------------------
#endif
! LUMINOUS_PROCESSES

  !*****************************************************************************************************************
  !
  ! SUBROUTINES
  !
  !*****************************************************************************************************************

  subroutine disc_void(disc)

    ! INIT OR VOID A DISC COMPONENT

    implicit none
    type(disc_type),intent(inout) :: disc   ! a disc component

    ! Evolution
    disc%life_time           = 0.d0                       ! will be incremented in disc_evolve_II
    disc%age_form            = 0.d0                       ! will be set in disc_evolve_II during the first evolution step
    disc%t_since_last_merger = 0.d0                       ! will be set in galaxy_merge
    !
    ! main components
    !
    call gas_void_gsh(disc%gsh_tab)                       ! void the gas structuration history table
    call stars_void(disc%stars)                           ! void the stellar population hosted by the disc
    call dust_void(disc%dust(1),init_geom='BC    ')       ! void the dust component associated to the young stellar population: BC
    call dust_void(disc%dust(2),init_geom='slab  ')       ! void the dust component associated to the old stelalr population: ISM
    call agn_void(disc%agn)                               ! void the agn component
    !
    ! disc properties
    disc%ngc               = 0
    disc%f_str             = -1.d0                        ! all other fields are set to null value
    disc%rd                = -1.d0
    disc%h                 = -1.d0
    disc%incl              = -1.d0
    disc%L                 = (/-1.d0,-1.d0,-1.d0/)
    disc%kappa             = -1.d0
    disc%V                 = -1.d0
    disc%dV                = -1.d0
    disc%Q                 =  1.d0
    disc%t_dyn             = -1.d0
    disc%t_cool            =  0.d0
    disc%cooling_timescale =  0.d0
    disc%Qturb_unstr       =  0.d0
    !
    ! interaction with the main halo hot gas (stripping)
    call gas_void(disc%stripping_rate)
    !
    ! transfer rates between gas phases
    call gas_void(disc%fresh_gas_acc_rate)
    call gas_void(disc%sfr)
    call gas_void(disc%str_rate)
    call gas_void(disc%disrupt_rate)
    call gas_void(disc%ejecta_rate)

    return
  end subroutine disc_void

  !*****************************************************************************************************************

  subroutine disc_copy(disc1,disc2)

    ! COPY THE disc2 COMPONENT IN THE disc1 COMPONENT

    implicit none

    type(disc_type),intent(inout) :: disc1
    type(disc_type),intent(in)    :: disc2

    ! Evolution
    disc1%life_time           = disc2%life_time
    disc1%age_form            = disc2%age_form
    disc1%t_since_last_merger = disc2%t_since_last_merger
    !
    ! main components
    disc1%gsh_tab = disc2%gsh_tab               ! copy the gas structuration history tab
    disc1%stars = disc2%stars                   ! copy the stellar population
    call dust_copy(disc1%dust(1),disc2%dust(1)) ! copy the dust component
    call dust_copy(disc1%dust(2),disc2%dust(2)) ! copy the dust component
    call agn_copy(disc1%agn,disc2%agn)          ! copy the SMBH component
    !
    ! disc properties
    disc1%nlevels            = disc2%nlevels    ! copy all other fields
    disc1%ngc                = disc2%ngc
    disc1%f_str              = disc2%f_str
    disc1%rd                 = disc2%rd
    disc1%h                  = disc2%h
    disc1%incl               = disc2%incl
    disc1%L                  = disc2%L
    disc1%kappa              = disc2%kappa
    disc1%V                  = disc2%V
    disc1%dV                 = disc2%dV
    disc1%Q                  = disc2%Q
    disc1%t_dyn              = disc2%t_dyn
    disc1%t_cool             = disc2%t_cool
    disc1%cooling_timescale  = disc2%cooling_timescale
    disc1%sfr_burst          = disc2%sfr_burst
    disc1%Qturb_unstr        = disc2%Qturb_unstr
    !
    ! interaction with the main halo hot gas (stripping)
    disc1%stripping_rate     = disc2%stripping_rate
    !
    ! transfer rates into the disc structure
    disc1%fresh_gas_acc_rate = disc2%fresh_gas_acc_rate
    disc1%sfr                = disc2%sfr
    disc1%str_rate           = disc2%str_rate
    disc1%disrupt_rate       = disc2%disrupt_rate
    disc1%ejecta_rate        = disc2%ejecta_rate

    return
  end subroutine disc_copy

  !*****************************************************************************************************************

  subroutine disc_evolve_I(disc,dm,bulge,r50,fresh_gas_acc_rate,gal_stripping_rate,dt_optim)

    ! COMPUTE FIRST PART OF DISC EVOLUTION SCHEME (PREDICTOR)
    ! COMPUTE INPUT, OUTPUT AND EVOLUTION RATES ASSOCIATED TO THE GAS STRUCTURATION HISTORY OF A GIVEN DISC
    ! These new evolution rates are saved into the disc structure 'disc'
    ! the impacts on masses are then computed and saved in 'disc_evolve_II'
    ! disc_evolve_I also compute optimal dt for disc evolution

    implicit none

    character(MAXPATHSIZE)         :: message                   ! a message to display

	real(kind=8),intent(in)        :: r50                       ! galaxy half mass radius
    real(kind=8)                   :: Vesc                      ! galaxy escape velocity
    real(kind=8),intent(out)       :: dt_optim                  ! optimal time-step for the disc evolution scheme
    real(kind=8)                   :: dt_agn                    ! optimal timestep for the agn component
    real(kind=8)                   :: dt_gas, dt_stars          ! optimal time-step for gas and stellar component
    real(kind=8)                   :: ej                        ! instantaneous ejecta rate
    real(kind=8)                   :: Qturb                     ! SN and AGN turbulent power injected into the fragmented gas phase
    real(kind=8)                   :: Qturb_unstr               ! residual Qturb injected in the unstructured/diffuse gas
    real(kind=8)                   :: f_unstr                   ! unstructured/diffuse mass fraction
    real(kind=8)                   :: f_str_in                  ! structured/fragmented mass fraction of the fresh accreted gas

    type(gas_type),intent(in)      :: fresh_gas_acc_rate        ! the overall accretion rate (= free-fall + cooling)
    type(gas_type),intent(in)      :: gal_stripping_rate        ! the stripping rate acting on the diffuse gas phase
    type(gas_type)                 :: rate                      ! gas transfer rate
    type(gas_type)                 :: tmp_rate                  ! a tmp rate
    type(gas_type)                 :: sfr                       ! instantaneous star formation rate
    type(gas_type)                 :: ejecta_rate               ! global ejecta rate (gas object form)
    type(gas_type)                 :: stripping_rate            ! the instantaneous stripping rate (gas object form)
    type(gas_type)                 :: str_rate                  ! structuring/fragmenting rate (from the diffuse to the fragmented non star-forming phases)
    type(gas_type)                 :: disrupt_rate              ! disrupting rate (from the fragmented to the diffuse gas phases)

    type(dm_type),intent(in)       :: dm                        ! a dark matter component
	type(bulge_type),intent(in)    :: bulge                     ! a bulge component
    type(disc_type),intent(inout)  :: disc                      ! a disc component

#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('disc_evolve_I',only_rank=rank,component='disc')
! -------------------------------------------------
#endif

    call gas_void(rate)           ! init
    call gas_void(ejecta_rate)    ! init
    call gas_void(str_rate)       ! init
    call gas_void(disrupt_rate)   ! init
    call gas_void(stripping_rate) ! init
    call gas_void(tmp_rate)       ! init
    call gas_void(sfr)            ! init
    !
    ! AGE FORM
    if ((disc%life_time .eq. 0.d0) .and. (gas_mass(fresh_gas_acc_rate) .gt. 0.d0)) then
        !
        ! set disc%age_form
        disc%age_form = dm%age_form + dm%life_time
    end if

    ! ***********************
    ! GAS
    ! ***********************
    call disc_update_gas_struct_history(disc) ! use for initialization, input and output rates are setled to 0.
    !
    !************************
    ! EVOLUTION RATES
    !************************
    !
    ! set the fraction of fragmented/structured gas in the fresh newly accreted gas
    f_str_in = fresh_gas_acc_rate%f_str
    !
    ! compute:
    ! - the star formation rate
    sfr = disc_SFR(disc)*disc_gas_signature(disc,component='structured',apply_as='rate_builder',called_by='sfr')
    ! - the ejecta rate produced by SN and AGN: ej
    ! the overall ejecta rate is stored in ej, the gas_object (ejecta_rate)
    ! will be progressively built as functions of the gas reservoirs
    ! - the turbulent power: Qturb
    Vesc = disc_escape_velocity(disc,bulge,dm,r50)
    call disc_compute_disc_feedback_activities(disc,Vesc,ejecta_rate=ej,Qturb=Qturb)
    !
    !************************
    ! MASS TRANSFER RATES
    !************************
    !
    ! RUN THROUGH GAS PHASES
    !
    ! OUTPUTS OF THE DIFFUSE GAS
    !
    call gas_void(rate) ! init, set to 0.
    !
    ! The diffuse gas is progressively converted in fragmented/structured gas
    ! the structuration rate is scaled according to the effective cooling time of the diffuse gas
    str_rate = gas_str_rate(disc,cooling_timescale=disc%cooling_timescale)* &
                    disc_gas_signature(disc,component='unstructured',apply_as='rate_builder',called_by='str_rate')
    rate = str_rate ! rate is initialize to str_rate
    !
    ! The diffuse gas can be ejected from the disc due to SN+AGN kinetic energy injection
    f_unstr = disc_gas_fraction(disc,component='unstructured') ! fraction of diffuse gas
    ! ejecta rate is progressivelly built accoring to the distribution of the gas between the different gas phases
    ! ejecta rate is initialized with the diffuse gas phase
    ejecta_rate = f_unstr*ej*disc_gas_signature(disc,component='unstructured',apply_as='rate_builder',called_by='ejecta_rate/unstr') ! init
    rate = rate + ejecta_rate ! rate is updated
    !
    ! The diffuse gas can be stripped by the interaction with an hot gas atmosphere
    ! the stripping_rate (flotting form) is converted into a gas object
    stripping_rate = gas_mass(gal_stripping_rate)*disc_gas_signature(disc,component='unstructured',apply_as='rate_builder',called_by='stripping_rate/unstr') ! init
    rate = rate + stripping_rate ! rate is updated
    !
    ! set the output rate of the diffuse gas phase
    disc%gsh_tab%out_rate(1) = rate
    !
    ! INPUTS OF THE STRUCTURED/FRAGMENTED GAS
    !
    ! The input gas is composed of :
    ! 1/ The gas transferred (condensed) from the diffuse to the structured/fragmented gas phases
    ! 2/ A fraction "f_str_in" of the fresh newly accreted gas, considered as already fragmented/structured
    str_rate = str_rate + f_str_in*fresh_gas_acc_rate
    ! 3/ The gas ejected by stellar winds and SN explosions is also strictly added to the structured/fragmented gas
    ! the stellar loss rate is pre-computed during stellar evolution in disc_evolve_II
    ! we use here the value computed and saved in star_evolve_II
    ! Metals contained into the stellar winds and produced by SN explosions are immediatly impacted by the disruption process
    ! and a fraction is therefore transferred to the diffuse gas phase
    ! set the input rate of the structured/fragmented gas phase
    disc%gsh_tab%in_rate(2) = str_rate + disc%stars%loss_rate
    !
    ! OUPUT OF THE STRUCTURED/FRAGMENTED GAS
    !
    call gas_void(rate) ! init, set to 0.
    !
    ! The structured/fragmented gas is progressively converted into stars 
    rate = sfr
    !
    ! The structured/fragmented gas can be disrupted by SN+AGN kinetic energy injection
    ! A fraction of the enery injected into the structured/fragmented gas allows to disrupt it
    ! The residual turbulent energy (non used to disrupt the dense gas) is saved in Qturb_unstr
    disrupt_rate = disrupt_str_gas_rate(disc,Qturb,Qturb_unstr=Qturb_unstr)* &
                    disc_gas_signature(disc,component='structured',apply_as='rate_builder',called_by='disrupt_rate')
    rate = rate + disrupt_rate
    !
    ! The structured/fragmented gas can be ejected from the disc due to SN+AGN kinetic energy injection
    ! update ejecta_rate
    call gas_void(tmp_rate)
    tmp_rate = (1.d0-f_unstr)*ej*disc_gas_signature(disc,component='structured',apply_as='rate_builder',called_by='ejecta_rate/str')
    ! tmp_rate is used to update
    ejecta_rate = ejecta_rate + tmp_rate
    ! and update
    rate = rate + tmp_rate
    ! set the output rate of the structured/fragmented gas phase
    disc%gsh_tab%out_rate(2) = rate
    !
    ! INPUTS OF THE UNSTRUCTURED/DIFFUSE GAS
    !
    ! the diffuse gas is composed by:
    ! 1/ the diffuse fresh gas accreted by the disc
    ! 2/ the dense/fragmented gas disrupted by SN+AGN
    disc%gsh_tab%in_rate(1) = (1.d0 - f_str_in)*fresh_gas_acc_rate + disrupt_rate
    !
    !************************************
    ! TIME-STEP ASSOCIATED TO COMPONENTS
    !************************************
    !
    ! *** gas ***
    ! compute optimal time-step for gas components,
    ! the lower one is the selected for the evolution
    dt_gas = 1.d15  ! init
    ! diffuse gas
    write(message,'(a)') 'disc_evolve_I / diffuse gas'
    dt_gas = min(dt_gas,gas_dt_optim(disc%gsh_tab%gas(1),disc%gsh_tab%in_rate(1),disc%gsh_tab%out_rate(1),called_by=trim(message)))
    ! structured/fragmented
    write(message,'(a)') 'disc_evolve_I / fragmented gas'
    dt_gas = min(dt_gas,gas_dt_optim(disc%gsh_tab%gas(2),disc%gsh_tab%in_rate(2),disc%gsh_tab%out_rate(2),called_by=trim(message)))
    !
    ! *** stars ***
    ! compute optimal time-step for the stellar component,
    dt_stars = -1.d0 ! init
    dt_stars = stars_dt_optim(disc%stars,sfr)
    !
    ! *** agn ***
    ! compute optimal time-step for the AGN component,
    ! set effective accretion rate
    call agn_evolve_I(disc%agn,dt_agn)
    !
    ! *****************************
    ! COMPUTE THE OPTIMAL TIME-STEP
    ! *****************************
    !
    dt_optim = -1.d0 ! init
    if (dt_gas .gt. 0.d0) dt_optim = dt_gas
    if (dt_stars .gt. 0.d0) then
        !
        if (dt_optim .gt. 0.d0) then
            !
            dt_optim = min(dt_stars,dt_optim)
        else
            !
            dt_optim = dt_stars
        end if
    end if
    !
    if (dt_agn .gt. 0.d0) then
        !
        if (dt_optim .gt. 0.d0) then
            !
            dt_optim = min(dt_agn,dt_optim)
        else
            !
            dt_optim = dt_agn
        end if
    end if
    !
    ! *****************************
    ! SET DISC PROPERTIES (rates, Qturb)
    ! *****************************
    !
    ! fresh gas accretion rate
    disc%fresh_gas_acc_rate = fresh_gas_acc_rate
    ! interaction with environment
    disc%stripping_rate = stripping_rate
    ! set Qturb_unstr
    disc%Qturb_unstr = Qturb_unstr
    ! the transfer rates are updated here
    ! mass reservoirs will be updated in disc_evolve_II
    ! set new value of the star formation rate
    disc%sfr = sfr
    if (disc%t_since_last_merger .eq. 0.d0) then
        !
        ! it is the first step after a merger event
        ! save the sfr at the merger time
        disc%sfr_burst = gas_mass(sfr)
    end if
    !
    ! set new value of the structuration rate
    disc%str_rate     = str_rate
    ! set new value of the disruption rate
    disc%disrupt_rate = disrupt_rate
    ! set new value of ejecta rate
    disc%ejecta_rate  = ejecta_rate

#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('disc_evolve_I ... done',only_rank=rank,component='disc')
! -------------------------------------------------
#endif

    return

    contains

    ! *************************************************

    function gas_str_rate(disc,cooling_timescale)

        ! COMPUTE AND RETURN THE GAS STRUCTURATION RATE
        ! We assume that - the structuration rate is governed by the cooling timescale of the diffuse gas
        !                - the diffuse gas is at 10^4 K
        ! At each time step we compute the fraction of unstructured gas that can cool and feed the GMCs

        implicit none

        real(kind=8)                      :: gas_str_rate      ! the structuration rate of the diffuse gas
        real(kind=8)                      :: M_unstr           ! the diffuse gas mass
        real(kind=8)                      :: f_unstr, f_Q      ! some efficiency coefficients
        real(kind=8)                      :: rho_unstr         ! unstructured/diffuse gas density
        real(kind=8)                      :: Z_unstr           ! the unstructured/diffuse gas metallicity
        real(kind=8)                      :: t_cool_ts
        real(kind=8),intent(out),optional :: cooling_timescale ! the cooling timescale at 10^4 K
        real(kind=8)                      :: x

        type(disc_type),intent(in)        :: disc              ! the disc component

        gas_str_rate      = 0.d0 ! init
        cooling_timescale = 0.d0 ! init

        if (disc%t_cool .le. 0.d0) return    ! cooling effective time is null

        M_unstr = disc_mass(disc,component='unstructured')
        !
        ! fraction of unstructured/diffuse gas
        f_unstr = disc_gas_fraction(disc,component='unstructured')
        !
        ! Take into account unstability factor
        f_Q = 1.d0 ! init
        if ((disc%Q .gt. 0.d0) .and. (disc%Q .le. 1.d0)) f_Q = Q_crit/disc%Q
        !
        ! compute cooling timescale associated to a 10^4 gas at metallicity Z
        ! compute the gas metallicity of the unstructured gas phase
        Z_unstr = gas_metallicity(disc_gas_signature(disc,component='unstructured'))
        ! compute the density of the unstructured/diffuse gas phase
        rho_unstr = 5.d-1*M_unstr / (pi*(1.68d0*disc%rd)**2.) / disc%h  ! [code unit]
        ! cooling timescale
        t_cool_ts = t_cool(diffuse_gas_temp,rho_unstr,Z_unstr)   ! [Gyr]
        if (present(cooling_timescale)) cooling_timescale = t_cool_ts
        !
        x = 0.d0
        if (t_cool_ts .gt. 0.d0) then
            !
            ! deduce the normalised cooling timescale
            x = disc%t_cool/t_cool_ts
            if (x .gt. 0.d0) then
                !
                gas_str_rate = f_unstr*f_Q*M_unstr*MFF(x)/t_cool_ts
            end if
        end if
        !
        if (is_NaN(gas_str_rate)) then
            call IO_print_error_message('gas_str_rate is NAN',only_rank=rank,called_by='gas_str_rate')
            stop ! stop the program
        end if

        return
    end function gas_str_rate

    ! *************************************************

    function disrupt_str_gas_rate(disc,Qturb,Qturb_unstr)

        ! COMPUTE THE DISRUPTION RATE OF THE STRUCTURED/FRAGMENTED GAS ASSOCIATED TO SN+AGN KINETIC INJECTION ENERGY

        implicit none

        real(kind=8),intent(in)            :: Qturb                ! turbulent power, injected by SN and AGN
        real(kind=8),intent(out),optional  :: Qturb_unstr          ! the residual turbulent power, then injected in the diffuse gas
        real(kind=8)                       :: disrupt_str_gas_rate ! the disruption rate of the structured/frgamented gas
        real(kind=8)                       :: f_str                ! the structured/frgamented mass fraction

        type(disc_type),intent(in)         :: disc                 ! the disc component

        disrupt_str_gas_rate = 0.d0  ! init
        if (present(Qturb_unstr)) Qturb_unstr = 0.d0  ! init

        f_str = disc_gas_fraction(disc,component='structured')

        if (Qturb .gt. 0.d0) then
            !
            ! Some turbulent energy is injected into the dense/fragmented gas
            disrupt_str_gas_rate = min(1.d5,2.d0*f_str*Qturb/disc%dV**2.)
            !
            if (present(Qturb_unstr)) then
                !
                ! In input Qturb is the full turbulent power
                ! If asked, in output Qturb_unstr saves the turbulent power
                ! that will be injected into the diffuse/unstructured gas phase
                Qturb_unstr = (1.d0-f_str)*Qturb
            end if
        end if

        return
    end function disrupt_str_gas_rate

    ! *************************************************
    ! end contains disc_evolve_I

  end subroutine disc_evolve_I

  !*****************************************************************************************************************

  subroutine disc_evolve_II(disc,dm,bulge,r50,dt)

    ! COMPUTE SECOND PART OF DISC EVOLUTION (CORRECTOR) :
    ! this subroutine evolve gas and stars component (mass) by using the optimal time-step (dt) and transfer rates computed previously in disc_evolve_I
    ! Concerning the gas structuration history, all informations about its evolution is stored in the gsh_tab
    ! For the gas, disc_evolve_II has to evolve the diffuse and the fragmented gas component
    ! For the stars, disc_evolve_II has to evolve stars component

    implicit none

    real(kind=8),intent(in)         :: r50                    ! galaxy half mass radius
    real(kind=8),intent(in)         :: dt                     ! effective evolving time of the disc
    real(kind=8)                    :: Uunstr,Ustr,dU         ! mass variation test
    real(kind=8)                    :: Vesc                   ! galaxy escape velocity

    type(gas_type)                  :: loss_rate_prior        ! mass loss rate used as prior value in disc_evolve_I
    type(gas_type)                  :: loss_rate              ! effective stellar mass loss rate
    type(gas_type)                  :: gas_return             ! return gas associated to star formation process
    type(gas_type)                  :: gas                    ! a gas component
    type(gas_type)                  :: dMacc                  ! accreted mass
    type(gas_type)                  :: dMdisrupt              ! disrupted mass
    type(gas_type)                  :: dMcond                 ! condensed mass

    type(dm_type),intent(in)        :: dm                     ! The dm component

    type(disc_type),intent(inout)   :: disc                   ! The disc component

    type(bulge_type),intent(in)     :: bulge                  ! The bulge component

#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('disc_evolve_II',only_rank=rank,component='disc')
! -------------------------------------------------
#endif

    if (dt .le. 0.d0) then
      !
      call IO_print_error_message('dt <= 0.',only_rank=rank,called_by = 'disc_evolve_II')
      call IO_print_message('use: ',only_rank=rank,param_name=(/'dt                       '/),real_param_val=(/dt/))
      stop ! stop the program
    end if

    if ((disc%life_time .eq. 0.d0) .and. (disc%age_form .gt. 0.d0)) then
        !
        ! print the starting state of the disc_component
        if (FOLLOW_UP .and. PR_FOLLOW_UP) then
            !
            call disc_print(follow_up_unit(current_index),'fits',disc)
        end if
    end if
    !
    ! *****************************
    ! STARS
    ! *****************************
    !
    ! save the previous mass loss rate, used as prior in disc_evolve_I
    loss_rate_prior = disc%stars%loss_rate
    call gas_void(gas_return)
    call stars_evolve(disc%stars,dt,loss_rate,sfr=disc%sfr,gas_return=gas_return)
    ! loss_rate is the real mass loss rate
    ! this value can be slightly different to loss_rate_prior, indeed the real stellar evolution takes AgeBin jumps into account
    ! the difference between the two values is taken into account after the gas reservoir update
    !
    ! *****************************
    ! GAS
    ! *****************************
    !
    ! save previous value for the diffuse gas phase
    Uunstr = disc_mass(disc,component='unstructured')
    ! save previous value for structured gas phase (without the star forming gas)
    Ustr = disc_mass(disc,component='structured')
    ! evolve the gas structuration history
    call disc_update_gas_struct_history(disc,dt=dt)
    !
    ! test mass variation
    if (Uunstr .gt. M_gas_min) then
        !
        dU = abs(disc_mass(disc,component='unstructured')-Uunstr)
        if ((abs(dU)/Uunstr) .gt. physical_precision) then
          !
          call IO_print_error_message('No quasi-static state of the diffuse gas component',only_rank=rank,called_by='disc_evolve_II')
          call IO_print_message('with',only_rank=rank,component='disc', &
                      param_name = (/'dU/U (%)                 ','U                        ','dU                       ',&
                                     'loss_rate_prior          ','loss_rate                '/), &
                      real_param_val = (/1.d2*dU/Uunstr,Uunstr,dU,gas_mass(loss_rate_prior),gas_mass(loss_rate)/))
                      write(*,*) 'dt = ', dt
          stop  ! stop the program
        end if
    end if
    !
    ! test mass variation
    if (Ustr .gt. M_gas_min) then
        !
        dU = abs(disc_mass(disc,component='structured')-Ustr)
        if ((abs(dU)/Ustr) .gt. physical_precision) then
          !
          call IO_print_error_message('No quasi-static state of the fragmented gas component',only_rank=rank,called_by='disc_evolve_II')
          call IO_print_message('with',only_rank=rank,component='disc', &
                      param_name = (/'dU/U (%)                 ','U                        ','dU                       ',&
                                     'sfr                      '/), &
                      real_param_val = (/1.d2*dU/Ustr,Ustr,dU,gas_mass(disc%sfr)/))
          stop  ! stop the program
        end if
    end if
    !
    ! Add gas_return to the fragmented component
    call gas_add(disc%gsh_tab%gas(2),gas_return)
    !
    ! *****************************
    ! SET disc PROPERTIES
    ! *****************************
    !
    if (disc_mass(disc) .gt. 0.d0) then
      !
      ! Some updating procedures use the latest accreted mass and/or the latest disrupted mass
      ! At this point new accreted, ejected of disrupted mass have been updated, therefore we have to
      ! evaluate them with associated rates
      ! dMacc
      dMacc = disc%fresh_gas_acc_rate*dt
      ! dMdisrupt
      dMdisrupt = disc%disrupt_rate*dt
      ! dMcond
      dMcond = disc%str_rate*dt
      !
      ! add time to the disc evolution time
      disc%life_time = disc%life_time + dt
      ! update t_since_last_merger
      disc%t_since_last_merger = disc%t_since_last_merger + dt
      ! update effective cooling clock
      call disc_update_cooling_clock(disc,dt=dt,gas=dMacc+dMdisrupt)
      ! compute the new exponential radius of the disc
      disc%rd        = disc_exponential_radius(disc,dm,accreted_mass=dMacc)
      ! compute the new disc inclination
      call disc_compute_inclination(disc,dm,accreted_mass=dMacc)
      ! computed the orbital velocity at the half mass radius
      disc%V         = disc_velocity(1.68d0*disc%rd,dm,disc,bulge)
      ! compute the epicyclic frequency at the disc half mass radius
      disc%kappa     = disc_kappa(1.68d0*disc%rd,dm,disc,bulge)
      ! computed dynamical time
      disc%t_dyn     = disc_dynamical_time(disc,dm,bulge)
      ! compute velocity dispersion
      Vesc = disc_escape_velocity(disc,bulge,dm,r50)
      disc%dV        = disc_update_velocity_dispersion(disc,Vesc,dm=dm,accreted_mass=dMacc,dt=dt)
      ! compute the scale height of the disc
      disc%h         = disc_scale_height(disc)
      ! compute Toomre parameter
      disc%Q         = disc_Toomre_parameter(disc)
      ! update inertial cascade properties t_emp, t_form, t_cascade, ngc
      call disc_update_inertial_cascade(disc,dt=dt,condensed_mass=dMcond)
      ! compute stracuration fraction
      disc%f_str     = disc_gas_fraction(disc,component='structured')
      !
      ! *****************************
      ! AGN
      ! *****************************
      !
      ! evolve the agn component
      call agn_evolve_II(disc%agn,dt,disc%incl)
      ! compute the dynamical acctetion time [r = 5pc / Vdisc( 5 pc )]
      call agn_set_dynamical_time(disc%agn,r_torus/disc_velocity(r_torus,dm,disc,bulge))
      !
      ! *****************************
      ! DUST
      ! *****************************!
      !
      ! evolve dust properties
      ! dust_evolve(d,gas,rc,incl=incl)
      ! BC
      ! Fragmented gas and dusts are assumed to be distributed homogeneously into the different GMCs
      ! build the complete gas object
      gas = disc_mass(disc,component='str')*disc_gas_signature(disc,component='str')
      ! Scale the gas component to a BC mass
      if (disc%ngc .gt. 0) then
        !
        gas = min(gas_mass(gas),Bonnor_Ebert_Mass(1.d0/disc%h))*gas_signature(gas)
        if (gas_mass(gas) > 0.d0) then
            !
            ! We compute here the effective extinction associated to 1 GMC
            call dust_evolve(disc%dust(1),gas,disc%h/2.d0)
        end if
      endif
      !
      ! ISM
      gas = disc_mass(disc,component='unstr')*disc_gas_signature(disc,component='unstr')
      if ((gas_mass(gas) .gt. 0.d0) .and. (disc%rd .gt. 0.d0)) then
        !
        call dust_evolve(disc%dust(2),gas,1.68d0*disc%rd,incl=disc%incl)
      end if
      !
      ! print instantaneous disc properties of the followed disc
      if (FOLLOW_UP .and. PR_FOLLOW_UP) then
        !
        call disc_print(follow_up_unit(current_index),'fits',disc)
      end if
    end if

#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('disc_evolve_II ... done',only_rank=rank,component='disc')
! -------------------------------------------------
#endif

    return
  end subroutine disc_evolve_II

  !*****************************************************************************************************************

  subroutine disc_compute_disc_feedback_activities(disc,Vesc,ejecta_rate,agn_acc_rate,Vwind,Qtherm,Qturb,Qrad)

    ! COMPUTE ALL PROPERTIES LINKED TO FEEDBACK PROCESSES INTO THE DISC STRUCTURE
    ! - the global gas ejection rate,
    ! - the global SMBH accretion rate,
    ! - the instantaneous velocity of the global wind
    ! - the instantaneous thermal, turbulent heating and radiative power generate by SN explosion and SMBH activity

    implicit none

    real(kind=8),intent(in)           :: Vesc             ! the escape velocity of the galaxy (in code unit)
    real(kind=8),intent(out),optional :: ejecta_rate      ! the global ejecta rate generated by SN explosion and/or SMBH activity
    real(kind=8),intent(out),optional :: agn_acc_rate     ! the real accretion rate onto the SMBH
    real(kind=8),intent(out),optional :: Vwind            ! the rate-weighted velocity of the wind produced by feedback processes
    real(kind=8),intent(out),optional :: Qtherm           ! thermal power (allow to compute mean wind temperature)
    real(kind=8),intent(out),optional :: Qturb            ! turbulent heating power
    real(kind=8),intent(out),optional :: Qrad             ! non kinetic and non-thermal power
    real(kind=8)                      :: SN_ej_rate       ! ejecta rate produced by SN explosions
    real(kind=8)                      :: AGN_ej_rate      ! ejecta rate produced by the SMBH activity
    real(kind=8)                      :: Qturb_SN         ! turbulent power produced by SN explosions
    real(kind=8)                      :: Qturb_AGN        ! turbulent power produced by the SMBH activity
    real(kind=8)                      :: Vwind_SN
    real(kind=8)                      :: Vwind_AGN
    real(kind=8)                      :: f_esc, Qt, T_ej, Vw

    type(gas_type)              :: ej                     ! a gas component
    type(disc_type),intent(in)  :: disc                   ! a disc component

    if (present(ejecta_rate)) then
        !
        ! init global ejecta rate
        ejecta_rate = 0.d0
    end if
    !
    if (present(Vwind)) then
      !
      ! init global (rate-weighted) velocity of the wind
      Vwind         = 0.d0
    end if
    !
    if (present(agn_acc_rate)) then
      !
      ! init the real accretion onto the SMBH
      agn_acc_rate  = 0.d0
    end if
    !
    if (present(Qtherm)) then
      !
      ! init the instantaneous thermal power
      Qtherm        = 0.d0
    end if
    !
    if (present(Qturb)) then
      !
      ! init the instantaneous turbulent power
      Qturb         = 0.d0
    end if
    !
    if (present(Qrad)) then
      !
      ! init the instantaneous radiative power
      Qrad          = 0.d0
    end if
    !
    ! Compute SN ejecta rate
    SN_ej_rate = SN_ejecta_rate(disc,Vesc)
    ! compute AGN ejecta rate
    AGN_ej_rate = agn_compute_agn_ejecta_rate(disc%agn,Vesc)
    if (present(ejecta_rate)) then
        !
        ! compute global ejecta rate
        ejecta_rate = SN_ej_rate + AGN_ej_rate
    end if
    !
    ! compute global (rate-weighted) velocity of the wind
    Vw = 0.d0  ! init
    if ((SN_ej_rate + AGN_ej_rate) .gt. 0.d0) then
      !
      Vwind_SN  = SN_wind_velocity(disc)
      Vwind_AGN = agn_compute_agn_wind_velocity(disc%agn,Vesc)
      Vw        = (SN_ej_rate*Vwind_SN + AGN_ej_rate*Vwind_AGN)/(SN_ej_rate + AGN_ej_rate)
    end if
    !
    if (present(Vwind)) then
      !
      Vwind = Vw
    end if
    !
    if (present(agn_acc_rate)) then
      !
      ! compute the real accretion onto the SMBH
      agn_acc_rate = agn_compute_AGN_accretion_rate(disc%agn)
    end if
    !
    Qt = SN_thermal_power(disc) + agn_compute_AGN_thermal_power(disc%agn)
    if (present(Qtherm)) then
      !
      ! compute the instantaneous thermal power
      Qtherm = Qt
    end if
    !
    if (present(Qrad)) then
      !
      ! compute the instantaneous radiative power
      Qrad  = SN_rad_power(disc) + agn_compute_AGN_rad_power(disc%agn)
    end if
    !
    ! compute the instantaneous turbulent heating power
    ! from SN
    Qturb_SN  = SN_turbulent_heating_power(disc)
    ! from AGN
    Qturb_AGN = agn_compute_AGN_turbulent_heating_power(disc%agn)
    if (present(Qturb)) then
       !
       Qturb = Qturb_SN + Qturb_AGN
    end if
    !
    ! apply correction factor to the ejecta rate due to escape galaxy velocity
    f_esc = 0.d0 ! init
    if (present(ejecta_rate)) then
       !
       if (ejecta_rate .gt. 0.d0) then
         !
         ! we assume a constant ejection process during a StellarTimeStep
         ! compute ejected_mass
         ej = ejecta_rate*StellarTimeStep*disc_gas_signature(disc,apply_as='rate_builder')
         if (gas_mass(ej) .gt. M_gas_null) then
            ! set temperature of ejecta
            call gas_inject_termal_energy(ej,Qt*StellarTimeStep)
            T_ej = gas_temp(ej)
            f_esc = 1.d0 - min(9.99d-1,max(1.d-3,Ronbint(Maxwell_Boltzman_Vdist_shifted,0.d0,Vesc,(/T_ej,Vw/),called_by='disc_compute_disc_feedback_activities')))
            ejecta_rate = f_esc*(SN_ej_rate + AGN_ej_rate)
            !
            if (present(Qturb)) then
                !
                ! Add power non used in real ejecta (1-fesc) into the turbulent motion
                Qturb = Qturb + (1.d0 - f_esc)*5.d-1*(SN_ej_rate + AGN_ej_rate)*Vw**2.
            end if
         else
            ! too few ejected mass 
			ejecta_rate = 0.d0
         end if
       end if
    end if

    return

    contains

    ! *************************************************

    function SN_ejecta_rate(disc, Vesc)

        ! COMPUTE EJECTA RATE OF THE DISC COMPONENT

        implicit none

		real(kind=8), intent(in)   :: Vesc    ! galaxy escpae velocity
        real(kind=8)               :: tot_gas ! the total gas mass in the disc
        real(kind=8)               :: Vwind   ! wind velocity
        real(kind=8)               :: eta_sn  ! SN event rate
        real(kind=8)               :: SN_ejecta_rate

        type(disc_type),intent(in) :: disc    ! the disc component

        SN_ejecta_rate = 0.d0 ! init

        tot_gas = disc_mass(disc,component='gas')

        if (disc_ejecta_efficiency .le. 0.d0) return  ! no ejecta takes into account
        !
        ! compute wind velocity 
        Vwind = SN_wind_velocity(disc)  ! in km/s
        !
        if (Vwind .le. 0.d0) return
        !
        eta_sn = disc%stars%SN_rate
        ! Energy driven
        ! in code unit 10^11 Msun / Gyr
        SN_ejecta_rate = 2.d0*disc_ejecta_efficiency*eta_sn*SN_kinetic_fraction*Esn_code_unit/Vwind**2
        !
        ! Check the result and crash the code if:
        ! disc_compute_disc_SN_ejecta_rate is NAN
        if (is_NaN(SN_ejecta_rate)) then
          !
          call IO_print_error_message('SN ejecta rate is NAN', &
                only_rank = rank, called_by = 'SN_ejecta_rate')
          call IO_print_message('used',only_rank=rank,component='disc',&
            param_name=(/'eta_sn                   ','tot_gas                  ','Vwind                    '/), &
            real_param_val =(/eta_sn,tot_gas,Vwind/))
          stop ! stop the program
        end if
        !
        ! SN_ejecta_rate is smaller than 0.d0
        if (SN_ejecta_rate .lt. 0.d0) then   ! check positivity of the ejecta
          !
          call IO_print_error_message('SN ejecta rate < 0.', &
                only_rank = rank, called_by = 'SN_ejecta_rate')
          call IO_print_message('used',only_rank=rank,component='disc',&
            param_name=(/'ejecta rate              ','eta_sn                   ', &
                         'tot_gas                  ','Vwind                    '/), &
            real_param_val =(/SN_ejecta_rate,eta_sn,tot_gas,Vwind/))
          stop ! stop the program
        end if

        return
      end function SN_ejecta_rate

      ! *************************************************

      function SN_wind_velocity(disc)

         ! RETURN WIND VELOCITY OF DISC EJECTA
         ! For a given power attributed to the ejecta process, two trends are possibles
         !  - a small amount of gas ejected with an high velocity
         !  - a large amount of gas ejected with a small velocity
         ! the balance between the two cases are evaluated according to the mass filling factor 
         ! more condensed/concentred the mass, larger the mass ejected (and smaller teh velocity)

         implicit none

         real(kind=8)                :: SN_wind_velocity  ! the velocity (in code unit) of the ejecta wind
         real(kind=8)                :: tot_gas
         real(kind=8),parameter      :: ref_mass_concentration = 1.d-2  ! 1Msun/pc3 in code unit
         real(kind=8),parameter      :: ref_SN_wind_velocity = 2.5d2/vel_code_unit_2_kmPerSec
         real(kind=8),parameter      :: min_SN_wind_velocity = 1.d2/vel_code_unit_2_kmPerSec
         real(kind=8)                :: mass_concentration
         real(kind=8)                :: M_unstr, M_str
         real(kind=8)                :: f_unstr
         real(kind=8)                :: scaling
         
         type(disc_type),intent(in)  :: disc              ! the disc component

         SN_wind_velocity = 0.d0 ! init

         if (disc_ejecta_efficiency .le. 0.d0)  return ! no ejecta 
         
         tot_gas = disc_mass(disc, component='gas')
         
         if ((disc%h .gt. 0.d0) .and. (tot_gas .gt. 0.d0)) then
			!
			! Evaluate the effective mass concentration
			M_unstr = disc_mass(disc, component='unstructured')
			M_str = disc_mass(disc, component='structured')
			f_unstr = disc_gas_fraction(disc, component='unstructured')
			mass_concentration = f_unstr*5.d-1*M_unstr / (pi*(1.68d0*disc%rd)**2.*disc%h) & ! diffuse gas
                                 + (1.d0-f_unstr)*M_str / (real(max(1,disc%ngc),kind=8)*(4.d0/3.d0*pi*(disc%h/2.d0)**3.)) ! fragmented gas
                                 
            !
            scaling = ref_mass_concentration/mass_concentration
			SN_wind_velocity = max(min_SN_wind_velocity,ref_SN_wind_velocity*scaling)                     
         else
			return
         end if 
         
         !SN_wind_velocity = 2.d2/vel_code_unit_2_kmPerSec

         return
      end function SN_wind_velocity

      ! *************************************************

      function SN_non_kinetic_power(disc)

        ! RETURN THE INSTANTANEOUS NON-KINETIC POWER PRODUCED BY SN EXPLOSIONS

        implicit none

        real(kind=8)               :: eta_sn                ! SN event rate
        real(kind=8)               :: SN_non_kinetic_power  ! in code unit (10^11.Msun.kpc^2/Gyr^3)

        type(disc_type),intent(in) :: disc                  ! the disc component

        SN_non_kinetic_power = 0.d0  ! init

        eta_sn = disc%stars%SN_rate
        SN_non_kinetic_power = &     ! in code unit (10^11.Msun.kpc^2/Gyr^3)
              eta_sn*(1.d0 - SN_kinetic_fraction)*Esn_code_unit

        return
      end function SN_non_kinetic_power

      ! *************************************************

      function SN_thermal_power(disc)

        ! RETURN THE INSTANTANEOUS THERMAL POWER PRODUCED BY SN EXPLOSIONS

        implicit none

        real(kind=8)               :: SN_thermal_power  ! in code unit (10^11.Msun.kpc^2/Gyr^3)

        type(disc_type),intent(in) :: disc              ! the disc component

        SN_thermal_power = SN_thermal_fraction*SN_non_kinetic_power(disc)

        return
      end function SN_thermal_power

      ! *************************************************

      function SN_rad_power(disc)

        ! RETURN THE INSTANTANEOUS RADIATIVE POWER (NON THERMAL AND NON KINETIC POWER) PRODUCED BY THE SN EXPLOSION

        implicit none

        real(kind=8)               :: SN_rad_power   ! in code unit (10^11.Msun.kpc^2/Gyr^3)

        type(disc_type),intent(in) :: disc                        ! the disc component

        SN_rad_power = (1.d0 - SN_thermal_fraction)*SN_non_kinetic_power(disc)

        return
      end function SN_rad_power

      ! *************************************************

      function SN_turbulent_heating_power(disc)

        ! COMPUTE THE TOTAL INSTANTANEOUS KINETIC POWER INJECTED IN THE TUBULENT GAS MOTION

        implicit none

        real(kind=8)               :: eta_sn                      ! SN event rate
        real(kind=8)               :: SN_turbulent_heating_power  ! in code unit (10^11.Msun.Mpc^2/Gyr^3)

        type(disc_type),intent(in) :: disc  ! the disc component

        SN_turbulent_heating_power = 0.d0  ! init

        if (disc_ejecta_efficiency .le. 0.d0) return  ! no ejecta takes into account

        eta_sn = disc%stars%SN_rate
        SN_turbulent_heating_power = &     ! in code unit (10^11.Msun.kpc^2/Gyr^3)
              (1.d0-disc_ejecta_efficiency)*eta_sn*SN_kinetic_fraction*Esn_code_unit

        return
      end function SN_turbulent_heating_power

      ! *************************************************
      ! end contains disc_compute_disc_feedback_activities

  end subroutine disc_compute_disc_feedback_activities

  !*****************************************************************************************************************

  subroutine disc_update_gas_struct_history(disc,dt,disc1,disc2,unstr_from_bulges,unstr_in_torus)

    ! UPDATE THE GAS STRUCTURATION HISTORY OF A GIVEN DISC
    ! Warning ! This subroutine must be called after the update of :
    !                             - the caracteristic exponential radius of the disc
    !                             - the disc scale height

    implicit none

    integer(kind=4)                      :: n                   ! loop index (structuration levels)

    character(MAXPATHSIZE)               :: message             ! a message to display

    real(kind=8),intent(in),optional     :: dt                  ! evolution time-step
    real(kind=8)                         :: M_comp              ! Mass of the component

    type(gas_type),intent(in),optional   :: unstr_from_bulges   ! nosfg coming from bulges (1 & 2)
    type(gas_type),intent(in),optional   :: unstr_in_torus      ! nosfg added to the AGN torus

    type(disc_type),intent(inout)        :: disc                ! The remnant disc (result of the merger)
    type(disc_type),intent(in),optional  :: disc1               ! The first progenitor disc component
    type(disc_type),intent(in),optional  :: disc2               ! The other progenitor disc component

#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('disc_update_gas_struct_history',only_rank=rank,component='disc')
! -------------------------------------------------
#endif
    !
    ! @ this point :
    ! - the new geometrical caracteristics of the (remnant) disc are computed and saved in disc
    ! - the remnant disc is "disc"
    !
    ! This subroutine is used in two different cases
    ! - A simple restructuration of a given disc
    ! - A merger of two differents gsh_tab, coming from two disc progenitors
    !
    if (present(disc1) .and. present(disc2) .and. present(unstr_from_bulges) .and. present(unstr_in_torus)) then
      !
      ! MERGER OF TWO DISC COMPONENTS
      !
      ! disc%gsh_tab is void
      ! disc1%gsh_tab contains gas structuration history of the first progenitor
      ! disc2%gsh_tab contains gas structuration history of the second progenitor
      !
      ! merge progenitors disc structuration histories
      disc%gsh_tab%gas = disc1%gsh_tab%gas + disc2%gsh_tab%gas
      !
      ! add diffuse/unstructured gas coming from bulge
      call gas_add(disc%gsh_tab%gas(1), unstr_from_bulges)
      !
      ! substract the nosfg added to the agn torus
      call gas_sub(disc%gsh_tab%gas(1), unstr_in_torus)
      !
    else
      !
      if (present(dt)) then
          !
          ! RESTRUCTURATION OF A GIVEN DISC
          !
          ! update gas structuration
          do n = 1, 2
              !
              ! apply input rate
              !
              disc%gsh_tab%gas(n) = disc%gsh_tab%gas(n) + disc%gsh_tab%in_rate(n)*dt
              !
              M_comp = gas_mass(disc%gsh_tab%gas(n))
              !
              if (M_comp .gt. 0.d0) then
                !
                ! apply output rate
                !
                if (gas_mass(disc%gsh_tab%out_rate(n)*dt) .lt. gas_mass(disc%gsh_tab%gas(n))) then
                   write(message,'(a,i2.2)') 'apply output rate on gsh_tab%gas(n), n = ', n
                   call gas_sub(disc%gsh_tab%gas(n),disc%gsh_tab%out_rate(n)*dt,called_by=message)
                else
                   call IO_print_error_message('gsh_tab%out_rate*dt > gsh_tab%gas',only_rank=rank,called_by = 'disc_update_gas_struct_history')
                   call IO_print_message('use: ',only_rank=rank, &
                               param_name=(/'dt                       ','M_gas_min                ','M_gas_crit               ','out_rate*dt              ','gas                      '/), &
                               real_param_val=(/dt,M_gas_min,M_gas_crit,gas_mass(disc%gsh_tab%out_rate(n)*dt),gas_mass(disc%gsh_tab%gas(n))/))
                   call IO_print_message('use: ',only_rank=rank, &
                               param_name=(/'n                        '/), int_param_val=(/n/))
                   stop ! stop the program
                end if
                !
              else
                if (gas_mass(disc%gsh_tab%out_rate(n)*dt) .gt. 0.d0) then
                  call IO_print_error_message('Try to substract mass from an empty component',called_by = 'disc_update_gas_struct_history')
                  stop ! stop the program
                end if
              end if
          end do
       else
         !
         ! INITIALISATION CASE
         !
         ! erase input and output rates of the gsh_tab
         call gas_void_array(disc%gsh_tab%in_rate)
         call gas_void_array(disc%gsh_tab%out_rate)
       end if
    end if
    !
#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('disc_update_gas_struct_history ... done',only_rank=rank,component='disc')
! -------------------------------------------------
#endif

    return
  end subroutine disc_update_gas_struct_history

  !*****************************************************************************************************************

  subroutine disc_update_cooling_clock(disc,dt,gas,disc1,disc2)

    ! UPDATE THE EFFECTIVE COOLING TIME OF THE UNSTRUCTURED GAS PHASE OF THE DISC
    ! The gas is assumed to be continously added to the reservoir,
    ! We take into account the half of the time-step (dt/2.d0),
    ! WARNING: When this routine is called the mass has already been added to the hot reservoir

    implicit none

    real(kind=8),intent(in),optional    :: dt     ! the lastest timestep
    real(kind=8)                        :: M1, M2 ! The unstructured gas mass: already present into the disc, just accreted

    type(gas_type),intent(in),optional  :: gas    ! the latest accreted and disrupted gas
    type(disc_type),intent(inout)       :: disc   ! a disc component
    type(disc_type),intent(in),optional :: disc1  ! first porgenitor
    type(disc_type),intent(in),optional :: disc2  ! second progenitor

    if ((present(disc1) .and. present(disc2))) then
        !
        ! merger case
        M1 = disc_mass(disc1,component='unstructured')  ! The mass of unstructured gas into the disc
        M2 = disc_mass(disc2,component='unstructured') ! The mass of unstructured gas into the disc2
        !
        if ((M1 + M2) .gt. 0.d0) then
            !
            disc%t_cool = (M1*disc1%t_cool + M2*disc2%t_cool)/(M1+M2)
        end if
    else
        !
        ! secular evolution
        M1 = disc_mass(disc,component='unstructured') ! The mass of unstructured gas into the disc
        M2 = gas_mass(gas)                            ! The mass of unstructured newly accreted gas
        !
        disc%t_cool = (max(0.d0,M1-M2)*(disc%t_cool + dt) + 5.d-1*dt*M2)/M1
    end if

    return
  end subroutine disc_update_cooling_clock

  !*****************************************************************************************************************

  subroutine disc_update_inertial_cascade(disc,dt,condensed_mass,disc1,disc2)

    ! UPDATE PROPERTIES OF THE GAS INERTIAL CASCADE
    ! timescales, global effective clock and the number of giant birth clouds

    implicit none

    type(gas_type),intent(in),optional  :: condensed_mass  ! the latest condensed gas
    type(disc_type),intent(inout)       :: disc            ! a disc component
    type(disc_type),intent(in),optional :: disc1           ! first porgenitor
    type(disc_type),intent(in),optional :: disc2           ! second progenitor

    real(kind=8),intent(in),optional    :: dt              ! the lastest timestep
    real(kind=8)                        :: M1, M2          ! 
    real(kind=8)                        :: MBE1, MBE2      ! 
    real(kind=8)                        :: t_emp           ! instantaneous emptying timescale
    real(kind=8)                        :: t_form          ! instantaneous formation timescale
    real(kind=8)                        :: t_eq            ! instantaneous equilibrium timescale
    real(kind=8)                        :: t_prod          ! instantaneous production timescale

    if ((present(disc1) .and. present(disc2))) then
        !
        ! merger case
        M1 = disc_mass(disc1,component='structured')  ! The mass of structured gas into the disc1
        M2 = disc_mass(disc2,component='structured')  ! The mass of structured gas into the disc2
        !
        ! timescales and clock are just updated according to a mass weighted computation
        if ((M1 + M2) .gt. 0.d0) then
          !
          disc%gsh_tab%t_cascade = (M1*disc1%gsh_tab%t_cascade + M2*disc2%gsh_tab%t_cascade)/(M1+M2)
          !
          disc%gsh_tab%t_eq = (M1*disc1%gsh_tab%t_eq + M2*disc2%gsh_tab%t_eq)/(M1+M2)
          !
          disc%gsh_tab%t_form = (M1*disc1%gsh_tab%t_form + M2*disc2%gsh_tab%t_form)/(M1+M2)
          disc%gsh_tab%t_form = min(disc%gsh_tab%t_form, 9.9d-1*disc%gsh_tab%t_eq)
          !
          disc%gsh_tab%t_emp = (M1*disc1%gsh_tab%t_emp + M2*disc2%gsh_tab%t_emp)/(M1+M2)
          !
          disc%gsh_tab%t_prod = (M1*disc1%gsh_tab%t_prod + M2*disc2%gsh_tab%t_prod)/(M1+M2)
          disc%gsh_tab%t_prod = min(disc%gsh_tab%t_prod, 9.9d-1*disc%gsh_tab%t_emp)
          !
          MBE1 = min(M1, Bonnor_Ebert_Mass(1./disc1%h))
          MBE2 = min(M2, Bonnor_Ebert_Mass(1./disc2%h))
          disc%ngc = 0
          if ((disc%h > 0.d0) .and. ((MBE1 + MBE2) > 0.)) then
             disc%ngc = max(1,int((MBE1*float(disc1%ngc) + MBE2*float(disc2%ngc)) / (MBE1 + MBE2)))
          end if
       else
			! Reset clock and timescales
			disc%gsh_tab%t_cascade = 0.d0
			disc%gsh_tab%t_eq      = t_eq_max
			disc%gsh_tab%t_form    = 9.9d-1*t_eq_max
			disc%gsh_tab%t_emp     = t_str_max
			disc%gsh_tab%t_prod    = 9.9d-1*t_str_max
			disc%ngc = 0
        end if
    else
        !
        ! secular evolution
        M1 = disc_mass(disc,component='structured')   ! The mass of structured gas into the disc
        M2 = gas_mass(condensed_mass)                 ! The mass of newly condensed gas
        !
        if (M1 .gt. M_gas_crit) then
          ! update the glogal effective clock of the inertial cascade
          ! it is assumed that the new incoming gas is continuously added
          ! therefore, in average the new gas has already spent dt/2 in the cascade
          disc%gsh_tab%t_cascade = (max(0.d0,M1-M2)*(disc%gsh_tab%t_cascade + dt) + dt/2.d0*M2)/M1
          !
          ! update the effective equilibrimum timescale of the inertial cascade
          t_eq = gas_timescales(disc%str_rate,disc%h,timescale='equilibrium')
          disc%gsh_tab%t_eq = (max(0.d0,M1-M2)*disc%gsh_tab%t_eq + M2*t_eq)/M1
          !
          ! update the effective formation timescale of the inertial cascade
          t_form = gas_timescales(disc%str_rate,disc%h,timescale='formation')
          disc%gsh_tab%t_form = (max(0.d0,M1-M2)*disc%gsh_tab%t_form + M2*t_form)/M1
          disc%gsh_tab%t_form = min(disc%gsh_tab%t_form, 9.9d-1*disc%gsh_tab%t_eq)
          !
          ! update the effective emptying timescale of the inertial cascade
          t_emp = gas_timescales(disc%str_rate,disc%h,timescale='structuration')
          disc%gsh_tab%t_emp = (max(0.d0,M1-M2)*disc%gsh_tab%t_emp + M2*t_emp)/M1
          !
          ! update the effective production timescale of the inertial cascade
          t_prod = Bonnor_Ebert_Mass(1./disc%h) / gas_mass(disc%str_rate)
          disc%gsh_tab%t_prod = (max(0.d0,M1-M2)*disc%gsh_tab%t_prod + M2*t_prod)/M1
          disc%gsh_tab%t_prod = min(disc%gsh_tab%t_prod, 9.9d-1*disc%gsh_tab%t_emp)
          !
          ! update the number of giant birth clouds
          if (disc%h > 0.d0) then
             MBE1 = min(M1, Saturation_Mass(1./disc%h, disc%str_rate))
             disc%ngc = int(max(1,ceiling(MBE1/Bonnor_Ebert_Mass(1./disc%h))))
          end if
        else
			! Reset clock and timescales
			disc%gsh_tab%t_cascade = 0.d0
			disc%gsh_tab%t_eq      = t_eq_max
			disc%gsh_tab%t_form    = 9.9d-1*t_eq_max
			disc%gsh_tab%t_emp     = t_str_max
			disc%gsh_tab%t_prod    = 9.9d-1*t_str_max
			disc%ngc = 0
        end if
    end if

	! Check NaN Values
	if (is_NaN(disc%gsh_tab%t_cascade) .or. is_NaN(disc%gsh_tab%t_eq) .or. \
	    is_NaN(disc%gsh_tab%t_emp) .or. is_NaN(disc%gsh_tab%t_prod) .or. \
	    is_NaN(disc%gsh_tab%t_form)) then
		!
        call IO_print_error_message('gas timescale(s) is(are) NAN',only_rank=rank,called_by='disc_update_inertial_cascade')
        call IO_print_message('used',only_rank=rank,component='disc',& 
        param_name=(/'t_cascade                ','t_eq                     ','t_emp                    ',&
                     't_prod                   ','t_form                   ','M1                       ',& 
                     'M2                       ','h                        '/), &
        real_param_val =(/disc%gsh_tab%t_cascade,disc%gsh_tab%t_eq,disc%gsh_tab%t_emp, \
                          disc%gsh_tab%t_prod,disc%gsh_tab%t_form,M1,M2,disc%h/))
        stop ! stop the program
	endif

    return
  end subroutine disc_update_inertial_cascade

  !*****************************************************************************************************************

  subroutine disc_compute_inclination(disc,dm,accreted_mass,disc1,disc2,mu)

     implicit none

     type(gas_type),intent(in),optional  :: accreted_mass ! the new accreted mass (used in secular evolution computation)

     type(disc_type),intent(inout)       :: disc          ! a disc component
     type(disc_type),intent(in),optional :: disc1         ! the first disc component (used in disc merger computation)
     type(disc_type),intent(in),optional :: disc2         ! the second disc component (used in disc merger computation)

     type(dm_type),intent(in)            :: dm            ! a dark-matter component

     real(kind=8),intent(in),optional    :: mu            ! merger mass factor
     real(kind=8)                        :: norm
     real(kind=8)                        :: new_incl
     real(kind=8)                        :: Md1, Md2, Mgas

     disc%incl = -1.0  ! init

     if (disc%rd .le. 0.d0) return     ! no disc
     !
     ! Update disc angular momentum
     if (present(disc1) .and. present(disc2) .and. present(mu)) then
        !
        ! MERGER CASE
        Md1 = disc_mass(disc1)
        Md2 = disc_mass(disc2)
        if ((Md1 .gt. 0.d0) .and. (Md2 .gt. 0.d0)) then
            !
            if (mu .gt. epsilon_merge) then
                !
                ! Major merger
                ! reset on the dark matter structure properties
                disc%L = dm%L
            else
                !
                ! Minor merger
                ! mass-weighted angular momentum of the two progenitor discs
                disc%L = (Md1*disc1%L + Md2*disc2%L)/(Md1 + Md2)
            end if
        end if
     else
        !
        ! SECULAR EVOLUTION CASE
        if (present(accreted_mass)) then
            Mgas   = gas_mass(accreted_mass)
            Md1    = max(0.d0,disc_mass(disc)-Mgas)          ! we have to remove the accreted mass
            disc%L = (Md1*disc%L + Mgas*dm%L)/(Md1 + Mgas)   ! mass-weighted angular momentum
        end if
     end if
     !
     ! Update inclination
     norm = sqrt(sum(disc%L**2.))
     if (norm .gt. 0.d0) then
        !
        new_incl = acos(disc%L(3)/norm) ! radian
     else
        !
        call IO_print_error_message('norm <= 0. !',only_rank=rank,called_by='inclination')
     end if

     ! TEST new inclination and crash the code if is NAN
     if (is_NaN(new_incl)) then
        !
        call IO_print_error_message('new disc inclination is NAN',only_rank=rank,called_by='inclination')
        stop ! stop the program
     end if

     disc%incl = new_incl  ! radian

     return
  end subroutine disc_compute_inclination

  !*****************************************************************************************************************

  subroutine disc_deallocate(disc)

    implicit none

    type(disc_type),intent(inout) :: disc   ! a disc component

    call gas_deallocate_gsh(disc%gsh_tab) ! deallocate the gas structuration history table
    call stars_deallocate(disc%stars)     ! deallocate the star formation history table

    return
  end subroutine disc_deallocate

  !*****************************************************************************************************************
  !
  ! FUNCTIONS
  !
  !*****************************************************************************************************************

  function disc_SFR(disc)

    ! COMPUTE STAR FORMATION RATE IN A DISC COMPONENT
    implicit none

    real(kind=8)               :: disc_SFR  ! the instantaneous star formation rate
    real(kind=8)               :: M_str     ! Mass of structured gas
    real(kind=8)               :: dt
    real(kind=8)               :: a, b 

    type(disc_type),intent(in) :: disc      ! the disc component

    disc_SFR = 0.d0 ! init
    M_str = disc_mass(disc, component='str')

    if ((M_str .gt. 0.d0) .and. (disc%gsh_tab%t_cascade .gt. 0.d0) .and. (disc%h .gt. 0.d0) .and. (disc%ngc .gt. 0)) then
        !
        ! The fragmented NON star-forming gas is converted in fragmented star-forming gas
        if (disc%gsh_tab%t_cascade .gt. disc%gsh_tab%t_eq) then
           ! The cascade reaches its equilibrium
           dt = disc%gsh_tab%t_prod
        else
           if (disc%gsh_tab%t_prod .lt. disc%gsh_tab%t_emp) then
              a = (disc%gsh_tab%t_prod - disc%gsh_tab%t_emp) / (disc%gsh_tab%t_eq - disc%gsh_tab%t_form)
              b = disc%gsh_tab%t_emp - a * disc%gsh_tab%t_form
              dt = a * disc%gsh_tab%t_cascade + b
           else
              dt = disc%gsh_tab%t_prod
           end if
        end if
        if (dt .gt. 0.d0) then
           disc_SFR = min(M_str, Bonnor_Ebert_Mass(1.d0/disc%h)) / dt
        else
           !
		   call IO_print_error_message('dt < 0',only_rank=rank,called_by='disc_SFR')
           call IO_print_message('used',only_rank=rank,component='disc',& 
           param_name=(/'t_cascade                ','t_eq                     ','t_emp                    ',&
                        't_prod                   ','t_form                   ','M_str                    ',& 
                        'h                        '/), &
           real_param_val =(/disc%gsh_tab%t_cascade,disc%gsh_tab%t_eq,disc%gsh_tab%t_emp, \
                             disc%gsh_tab%t_prod,disc%gsh_tab%t_form,M_str,disc%h/))
           stop ! stop the program
        end if
    end if

    return
  end function disc_SFR

  !*****************************************************************************************************************

  function disc_mass_surf_density_(r,param)

    ! RETURN THE DISC MASS SURFACE DENSITY

    implicit none

    real(kind=8),intent(in)   :: r            ! the disc radius
    real(kind=8),intent(in)   :: param(2)     ! parameter array : param(1) = disc%rd
                                              !                   param(2) = disc_mass(disc,component=component)

    real(kind=8)              :: disc_mass_surf_density_  ! the mass surface density

    disc_mass_surf_density_ = param(2)/(2.d0*pi*param(1)**2.)*exp(-r/param(1))

    return
  end function disc_mass_surf_density_

  !*****************************************************************************************************************

  function disc_mass(disc,r,component,element)

    ! RETURN THE DISC MASS
    ! if r is present, return the mass enclosed in the radius r

    implicit none

    character(*),intent(in),optional :: component   ! allows to select the disc component (diffuse, fragmented, stars, agn)
    character(*),intent(in),optional :: element     ! allows to select the gas element (H, He, C, ... )
    character(MAXPATHSIZE)           :: message     ! a message to display

    real(kind=8),intent(in),optional :: r           ! the disc radius
    real(kind=8)                     :: disc_mass   ! the disc mass
    real(kind=8)                     :: x           ! addimentionnal radius

    type(disc_type),intent(in)       :: disc        ! the disc component

    disc_mass = 0.d0

    ! compute total mass
    ! sum over structuration level
    disc_mass = gas_mass(disc%gsh_tab%gas(1),component=element)             ! diffuse gas
    disc_mass = disc_mass + gas_mass(disc%gsh_tab%gas(2),component=element) ! structured/fragemented gas
    ! add stellar and agn mass
    disc_mass = disc_mass + stars_mass(disc%stars) + agn_mass(disc%agn)

    if (disc_mass .le. 0.d0) return ! no disc
    ! @ this point the disc contains mass

    if (present(component)) then
        !
        ! A specific component of the gas is selected
        select case (trim(component))
        case ('gas')
            disc_mass = gas_mass(disc%gsh_tab%gas(1),component=element)             ! diffuse gas
            disc_mass = disc_mass + gas_mass(disc%gsh_tab%gas(2),component=element) ! structured/fragemented gas
        case ('agn','SMBH','black-hole')
            disc_mass = agn_mass(disc%agn,component=component)
        case ('without_SMBH','no_SMBH','no_AGN')
            disc_mass = disc_mass - agn_mass(disc%agn)
        case('structured','str','fragmented')
            disc_mass = gas_mass(disc%gsh_tab%gas(2),component=element)
        case('unstructured','unstr','diffuse')
            disc_mass = gas_mass(disc%gsh_tab%gas(1),component=element)
        case ('stars')
            disc_mass = stars_mass(disc%stars)
        case default
            write(message,'(a,a,a)') 'Keyword ', trim(component), ' not defined'
            call IO_print_error_message(message,only_rank=rank,called_by='disc_mass')
            stop  ! stop the program
        end select
    else
        !
        ! No specific component
        if (present(element)) then
            !
            ! A specific element of the gas is selected
            disc_mass = 0.d0
            ! sum over structuration level
            disc_mass = gas_mass(disc%gsh_tab%gas(1),component=element)              ! diffuse gas
            disc_mass = disc_mass + gas_mass(disc%gsh_tab%gas(2),component=element)  ! structured/fragemented gas
        end if
    end if

    if (disc_mass .le. 0.d0) then ! no disc
      !
      disc_mass = 0.d0
      return
    end if
    !
    ! ADAPT TO THE RADIUS R
    if (present(r)) then
      !
      if (r .le. 0.d0) then
        !
        call IO_print_error_message('r <= 0.',only_rank=rank,called_by='disc_mass')
        call IO_print_message('with',only_rank=rank,component='disc', &
               param_name = (/'r                        ','rd                       ','disc_mass                '/), &
               real_param_val  = (/r,disc%rd,disc_mass/))
        stop  ! stop the program
      end if
      !
      if (disc%rd .le. 0.d0) then
        !
        call IO_print_error_message('rd <= 0.',only_rank=rank,called_by='disc_mass')
        call IO_print_message('with',only_rank=rank,component='disc', &
               param_name = (/'r                        ','rd                       '/), &
               real_param_val  = (/r,disc%rd/))
        stop  ! stop the program
      endif
      !
      x     = r/disc%rd                                ! compute adimentionnal radius
      disc_mass = disc_mass*(1.0d0-exp(-x)-x*exp(-x))  ! apply geometrical function
    end if
    !
    ! CHECK the results and crash the code if disc_mass is NAN
    if (is_NaN(disc_mass)) then
      !
      call IO_print_error_message('disc_mass is NaN',only_rank=rank,called_by='disc_mass')
      stop  ! stop the program
    end if
    !
    return
  end function disc_mass

  !*****************************************************************************************************************

  function disc_gas_signature(disc,component,apply_as,called_by)

    ! RETURN A GAS SIGNATURE (A NORMALIZED GAS OBJECT)

    implicit none

    character(*),intent(in),optional  :: component          ! allow to select the disc component (sfg, non-sfg)
    character(*),intent(in),optional  :: called_by          ! name of the function which has called this function
    character(*),intent(in),optional  :: apply_as           ! the gas signature function can be use in different cases
                                                            ! when this function is used as an output rate builder, we have to check critical mass
    character(MAXPATHSIZE)            :: message            ! a message to display

    real(kind=8)                      :: f_unstr            ! mass fraction of unstructured/diffuse gas
    
    type(gas_type)                    :: disc_gas_signature ! the gas signature

    type(disc_type),intent(in)        :: disc               ! a disc component

    call gas_void(disc_gas_signature)

    f_unstr = disc_gas_fraction(disc, component='unstr')
    ! gas signature of the complete disc
    ! WARNING, as in disc_evolve_I, the gas signature has to be applied to each component 
    ! (because of the M_gas_crit threshold)
    disc_gas_signature = f_unstr*gas_signature(disc%gsh_tab%gas(1),apply_as=apply_as,called_by=called_by)  ! diffuse gas
    disc_gas_signature = disc_gas_signature + (1.d0-f_unstr)*gas_signature(disc%gsh_tab%gas(2),apply_as=apply_as,called_by=called_by)  ! diffuse gas ! structured/fragmented gas

    if (present(component)) then
        !
        select case (trim(component))
        case('unstructured','unstr','diffuse')
            disc_gas_signature = gas_signature(disc%gsh_tab%gas(1),apply_as=apply_as,called_by=called_by)
        case('structured','str','fragmented')
            disc_gas_signature = gas_signature(disc%gsh_tab%gas(2),apply_as=apply_as,called_by=called_by)
        case ('all','gas')
            disc_gas_signature = disc_gas_signature
        case default
            write(message,'(a,a,a)') 'Keyword ', trim(component), ' not defined'
            call IO_print_error_message(message,only_rank=rank,called_by='disc_gas_signature')
            stop  ! stop the program
        end select
    end if

    return
  end function disc_gas_signature

  !*****************************************************************************************************************

  function disc_exponential_radius(disc,dm,accreted_mass,disc2)

     ! RETURN THE EXPONENTIAL RADIUS OF THE DISC COMPONENT
     ! WARNING !! When this subroutine is called, the new accreted mass is already added to the previous disc mass

     implicit none

     real(kind=8)                        :: new_rd        ! the new disc exponential radius
     real(kind=8)                        :: old_rd        ! backup of the old exponential radius of the disc
     real(kind=8)                        :: dm_rd         ! disc exponential radisu deduced of the dark-matter structure (used for new accretion)
     real(kind=8)                        :: old_mass      ! old mass of the disc (without the diffuse gas = stars + fragmented gas)
     real(kind=8)                        :: acc_mass      ! the newly accreted mass
     real(kind=8)                        :: old_diff_mass ! the old mass of diffuse gas (before the newly accreted)
     real(kind=8)                        :: pi_Sigma, W
     real(kind=8)                        :: disc_exponential_radius

     type(gas_type),intent(in),optional  :: accreted_mass ! the new accreted mass (used in secular evolution computation)

     type(disc_type),intent(in)          :: disc          ! a disc component
     type(disc_type),intent(in),optional :: disc2         ! a second disc component (used in disc merger computation)

     type(dm_type),intent(in)            :: dm            ! a dark-matter component

     disc_exponential_radius = 0.d0

     if (present(disc2)) then
       !
       ! merger case
       !
       ! In this case, the exponential radius of the remanent disc is setted to
       ! the larger exponential radius of the two progenitors.
       new_rd = max(disc%rd,disc2%rd)
     else
       !
       ! Secular evolution
       !
       if (disc_mass(disc) .le. 0.d0) return     ! no disc
       !
       ! save the previous value of the exponential disc radius
       old_rd = disc%rd
       !
       ! compute the exponential radius linked to the new structure of the dark matter halo
       dm_rd  = dm%spin*dm%R_halo/sqrt(2.d0)
       !
       if (present(accreted_mass)) then
          !
          acc_mass = gas_mass((1.d0-accreted_mass%f_str)*accreted_mass) ! the mass of the new accreted unstructured gas
       else
          !
          acc_mass = 0.d0
       end if
       !
       ! computed the current fragmented/structured disc mass
       ! this mass takes into account all the disc mass without the diffuse gas
       ! the diffuse gas is expected to rich a disc size 2 times larger than the stellar disc
       old_mass = max(0.d0,disc_mass(disc) - disc_mass(disc,component='unstr'))
       !
       ! compute the diffuse gas mass contained into the disc before the accretion
       ! substarct the diffuse accreted mass (already added in the diffuse gas reservoir in disc_evolve_II)
       old_diff_mass = max(0.d0,disc_mass(disc,component='unstr') - acc_mass)
       !
       ! compute new rd
       ! The new exponential radius is computed by using mass-weighted computations
       ! we consider that the new accreted gas forms a disc
       ! with a exponential radius computed as function of the new properties of the dark-matter halo
       pi_Sigma = 0.
       if (dm_rd .gt. 0.d0)  pi_Sigma = pi_Sigma + (acc_mass/dm_rd)**2.
       !
       if (old_rd .gt. 0.d0) then
            !
            W = 2.d0
            if (2.d0*old_rd .ge. dm_rd) W = 1.0
            pi_Sigma = pi_Sigma + (old_mass/old_rd)**2. + (old_diff_mass/min(2.d0*old_rd,dm_rd))**2.
       end if
       !
       if (pi_Sigma .gt. 0.d0) then
           !
           new_rd = disc_mass(disc)**2.*sqrt(1.d0/pi_Sigma)/(W*old_diff_mass + acc_mass + old_mass)
       else
           !
           call IO_print_error_message('pi_Sigma <= 0.',only_rank=rank,called_by='disc_exponential_radius')
           stop ! stop the program
       end if
     end if
     !
     ! TEST new rd and crash the code if the value is smaller than 0.
     if ((new_rd .le. 0.d0) .and. (disc_mass(disc) .gt. 0.d0)) then
       !
       call IO_print_error_message('new disc rd <= 0.',only_rank=rank,called_by='disc_exponential_radius')
       stop ! stop the program
     end if
     !
     ! TEST new rd and crash the code if new_size is NAN
     if (is_NaN(new_rd)) then
       !
       call IO_print_error_message('new disc rd is NAN',only_rank=rank,called_by='disc_exponential_radius')
       call IO_print_message('used',only_rank=rank,component='disc', &
               param_name=(/'Md                       ','new_rd                   ', &
                            'old_rd                   ','accreted mass            '/), &
               real_param_val=(/disc_mass(disc),new_rd,old_rd,acc_mass/))
       stop ! stop the program
     end if
     !
     ! set
     disc_exponential_radius = max(5.d-2,new_rd)

     return
  end function disc_exponential_radius

  !*****************************************************************************************************************

  function disc_scale_height(disc)

    ! RETURN THE SCALE HEIGHT OF THE DISC COMPONENT

    implicit none

    real(kind=8)                  :: disc_scale_height

    type(disc_type),intent(inout) :: disc                           ! a disc component

    ! We simply inverse the larson velocity dispersion relation
    disc_scale_height = l_star*(sig_star/disc%dV)**(-1./larson_sig_slope) ! [code unit]

    return
  end function disc_scale_height

  !*****************************************************************************************************************
  
  function disc_escape_velocity(disc,bulge,dm,r50)
  
	! RETURN THE ESCAPE VELOCITY OF THE DISC
	
	real(kind=8),intent(in)       :: r50     ! galaxy half mass radius (including dark-matter)
	real(kind=8)                  :: r_esc   ! radius of escape velocity computation
	
	type(disc_type),intent(in)    :: disc    ! the disc component
    type(bulge_type),intent(in)   :: bulge   ! the bulge component
    type(dm_type),intent(in)      :: dm      ! the dark matter halo component
     
    r_esc = max(1.68d0*disc%rd,r50)
   
    disc_escape_velocity = -1.d0
    
    if (r_esc .gt. 0.d0) then 
		disc_escape_velocity = sqrt(2.d0*gravconst_code_unit*( &
                                    disc_mass(disc,r=r_esc) + &
                                    bulge_mass(bulge,r=r_esc) + &
                                    dm_mass(dm,r=r_esc))/r_esc)
    end if
    
    return
  end function disc_escape_velocity  
   	
  !*****************************************************************************************************************
  
  function disc_mass_weighted_orbital_velocity(disc,dm,bulge)

     implicit none

     real(kind=8)                  :: disc_mass_weighted_orbital_velocity
     real(kind=8)                  :: r_ext   ! external radius (set as 11*disc%rd, enclose 99.999% of the disc mass)

     type(disc_type),intent(in)    :: disc    ! the disc component
     type(bulge_type),intent(in)   :: bulge   ! the bulge component
     type(dm_type),intent(in)      :: dm      ! the dark matter halo component

     disc_mass_weighted_orbital_velocity = 0.d0

     if (disc%rd .gt. 0.d0) then
        !
        r_ext = 1.d1*disc%rd
        disc_mass_weighted_orbital_velocity = 1.d0/disc_mass(disc)* &
                        Ronbint(V_sigma_r,r_sf,r_ext, &
                                    (/dm%rho_core,dm%r_core,disc%rd,disc_mass(disc),bulge%rb,bulge_mass(bulge)/), &
                                    called_by = 'disc_mass_weighted_orbital_velocity')
     endif

     return
  end function disc_mass_weighted_orbital_velocity

  !*****************************************************************************************************************

  function disc_dynamical_time(disc,dm,bulge)

     ! RETURN THE DYNAMICAL TIME OF THE DISC (rd/V)
     ! WARNING !! This function must be called after disc_exponential_radius

     implicit none

     real(kind=8)                :: disc_dynamical_time
     real(kind=8)                :: r

     type(disc_type),intent(in)  :: disc    ! the disc component
     type(bulge_type),intent(in) :: bulge   ! the bulge component
     type(dm_type),intent(in)    :: dm      ! the dark matter halo component

     if (disc%rd .gt. 0.d0) then
        !
        r = 2.2d0*disc%rd
        disc_dynamical_time = r/disc_velocity(r,dm,disc,bulge)
     end if

     return
  end function disc_dynamical_time

 !*****************************************************************************************************************
 
  function disc_gas_fraction(disc,component)

     ! RETURN THE MASS FRACTION OF A GIVEN COMPONENT

     implicit none

     integer(kind=4)                   :: e                  ! loop indexes

     character(MAXPATHSIZE)            :: message            ! a message to display
     character(*),intent(in),optional  :: component          ! allow to select the bulge component (sf, str, stars)
     character(MAXPATHSIZE)            :: subcomponent

     real(kind=8)                      :: M_gas              ! total gas mass in the disc
     real(kind=8)                      :: disc_gas_fraction  ! fraction of mass structured in the disc

     type(disc_type),intent(in)        :: disc               ! the disc component

     disc_gas_fraction = 0.d0 ! init

     M_gas = disc_mass(disc,component='gas')
     if (M_gas .le. 0.d0) return

     if (present(component)) then
        !
        subcomponent = 'all' ! the default values
        if ((trim(component) .eq. 'Metals') .or. (trim(component) .eq. 'metals')) then
          !
          subcomponent = trim(component)
        else
          !
          do e = 1, nElts
            if (trim(component) .eq. trim(Elt_names(e))) then
                subcomponent = trim(Elt_names(e))
                exit
            end if
          end do
        end if
        !
        select case (trim(component))
        
        case ('gas')
            !
            disc_gas_fraction = gas_mass(disc%gsh_tab%gas(1),component=subcomponent) ! diffuse gas
            disc_gas_fraction = disc_gas_fraction + gas_mass(disc%gsh_tab%gas(2),component=subcomponent) ! structured gas
        case('structured','str','fragmented')
            !
            disc_gas_fraction = gas_mass(disc%gsh_tab%gas(2),component=subcomponent)
        case('unstructured','unstr','diffuse')
            !
            disc_gas_fraction = gas_mass(disc%gsh_tab%gas(1),component=subcomponent)
        case default
            if (trim(subcomponent) .eq. trim(component)) then
                !
                disc_gas_fraction = gas_mass(disc%gsh_tab%gas(1),component=subcomponent) ! diffuse gas
                disc_gas_fraction = disc_gas_fraction + gas_mass(disc%gsh_tab%gas(2),component=subcomponent) ! structured gas
            else
                !
                write(message,'(a,a,a)') 'Keyword ', trim(component), ' not defined'
                call IO_print_error_message(message,only_rank=rank,called_by='disc_gas_fraction')
                stop  ! stop the program
            end if
         end select
     end if

     disc_gas_fraction = disc_gas_fraction / M_gas
     return
  end function disc_gas_fraction

  !*****************************************************************************************************************

  function isolated_disc_velocity(r,disc) ! INTERFACE VERSION

    ! COMPUTE THE ROTATIONNAL VELOCITY OF A ISOLATED DISC
    ! use FREEMAN+1970 eqs. 10 & 11

    implicit none

    real(kind=8), intent(in)    :: r                      ! radius (in code unit: Mpc)
    real(kind=8)                :: isolated_disc_velocity ! in kpc/Gyr

    type(disc_type), intent(in) :: disc                   ! a disc component

    isolated_disc_velocity = isolated_disc_velocity_(r,(/disc%rd,disc_mass(disc,component='no_AGN')/))

    return
  end function isolated_disc_velocity

  !*****************************************************************************************************************

  function isolated_disc_velocity_(r,param)

    ! COMPUTE THE ROTATIONNAL VELOCITY OF A ISOLATED DISC
    ! use FREEMAN+1970 eqs. 10 & 11

    implicit none

    real(kind=8), intent(in)    :: r          ! radius (in code unit: kpc)
    real(kind=8),intent(in)     :: param(2)   ! parameter array : param(1) = disc%rd
                                              !                   param(2) = disc_mass(disc,component='no_AGN')

    real(kind=8)                :: v                       ! disc velocity
    real(kind=8)                :: isolated_disc_velocity_ ! in kpc/Gyr
    real(kind=8)                :: u                       ! a tmp variable

    u = r/(2.d0*param(1))
    v = gravconst_code_unit*param(2)/(2.d0*param(1)**3.)*r**2. ! in code unit
    v = v*(BESSI0(u)*BESSK0(u) - BESSI1(u)*BESSK1(u))          ! in code unit
    isolated_disc_velocity_ = sqrt(v)                          ! in kpc/Gyr

    return
  end function isolated_disc_velocity_

  !*****************************************************************************************************************

  function disc_velocity(r,dm,disc,bulge) ! INTERFACE VERSION

    ! COMPUTE THE GLOBAL ROTATIONNAL VELOCITY OF THE DISC
    ! Take into account the dark matter and the bulge structure.

    implicit none

    real(kind=8), intent(in)     :: r             ! radius (in kpc)
    real(kind=8)                 :: disc_velocity ! in kpc/Gyr

    type(dm_type),intent(in)     :: dm            ! a dark matter component
    type(disc_type), intent(in)  :: disc          ! a disc component
    type(bulge_type), intent(in) :: bulge         ! a bulge component

    disc_velocity = disc_velocity_(r,(/dm%rho_core,dm%r_core,disc%rd,disc_mass(disc,component='no_AGN'),agn_mass(disc%agn),bulge%rb,bulge_mass(bulge,r=r)/))

    return
  end function disc_velocity

  !*****************************************************************************************************************

  function disc_velocity_(r,param)

    ! COMPUTE THE GLOBAL ROTATIONNAL VELOCITY OF THE DISC
    ! Take into account the dark matter and the bulge structure.

    implicit none

    real(kind=8), intent(in)   :: r              ! radius (in kpc)
    real(kind=8),intent(in)    :: param(7)       ! parameter array : param(1) = dm%rho_core
                                                 !                   param(2) = dm%r_core
                                                 !                   param(3) = disc%rd
                                                 !                   param(4) = disc_mass(disc,component='no_AGN')
                                                 !                   param(5) = agn_mass(disc%agn)
                                                 !                   param(6) = bulge%rb
                                                 !                   param(7) = bulge_mass(bulge)
    real(kind=8)               :: disc_velocity_ ! in kpc/Gyr

    disc_velocity_ = sqrt(isolated_disc_velocity_(r,param(3:4))**2. + &    ! disc
                          isolated_agn_velocity_(r,param(5))**2. + &       ! agn
                          isolated_bulge_velocity_(r,param(6:7))**2. + &   ! bulge
                          dm_velocity_(r,param(1:2))**2.)                  ! dm

    return
  end function disc_velocity_

  !*****************************************************************************************************************

  function disc_update_velocity_dispersion(disc,Vesc,disc2,dm,accreted_mass,dt)

     ! RETURN THE DISC VELOCITY DISPERSION
     ! The velocity dispersion is produced by injection of
     !   -a fraction of the kinetic energy of the accreted mass
     !   -a fraction of the kinetic energy generated byt SN and AGN
     !   -a fraction of the gravitatial interaction energy
     ! A fraction of the turbulent energy (associated to the velocity dispersion) is lost during each galaxy rotation cycle
     ! this routine has to be called after the update of disc%V, disc%t_dyn and disc%gsh_tab

     implicit none

	 real(kind=8),intent(in)              :: Vesc ! The galaxy escape velocity 
     real(kind=8),intent(in),optional     :: dt
     real(kind=8)                         :: disc_update_velocity_dispersion
     real(kind=8)                         :: M1, M2
     real(kind=8)                         :: Edisp1, Edisp2, Eint, Eturb
     real(kind=8)                         :: Egrav, Ekin
     real(kind=8)                         :: racc, t_disp
     real(kind=8)                         :: Vacc, dV
     real(kind=8)                         :: fdisp

     type(gas_type),intent(in),optional   :: accreted_mass ! the new accreted mass (used in secular evolution computation)

     type(disc_type),intent(in)           :: disc          ! a disc component (the first progenitor)
     type(dm_type),intent(in),optional    :: dm            ! a dark matter component
     type(disc_type),intent(in),optional  :: disc2         ! a disc component (the second progenitor

     disc_update_velocity_dispersion = disc%dV ! init to the previous value

     if (present(dt)) then
        !
        ! secular evolution case
        !
        ! accreted gaseous disc
        M1        = gas_mass(accreted_mass)
        ! Gravitational energy of the accreted gas
        Egrav = gravconst_code_unit*dm%M_vir * M1 / dm%R_vir
        ! half mass radius of the accreted gaseous disc
        racc      = dm%spin*dm%R_halo/sqrt(2.d0)
        ! orbital velocity of the accreted gaseous disc
        Vacc = disc_velocity_(racc,(/dm%rho_core,dm%r_core,1.68d0*racc,M1,0.d0,0.d0,0.d0/))
        ! kinetic energy of the accreted gaseous disc
        Ekin = 5.d-1*M1*Vacc**2.
        ! deduce the Turbulent residual energy
        Eturb = max(0.d0,Egrav-Ekin)
        fdisp = 0.d0 ! initially no dissipation
        if (Eturb .gt. 0.d0) then
            ! we assume that a fraction of the turbulent energy is lost during a dissipation timescale
            t_disp = min(1.d0,Ekin/Eturb)*(racc / Vacc)
            fdisp  = min(1.d0,dt/t_disp)
        end if
        Edisp1 = (1.d0-fdisp)*Eturb
        !
        ! pre-formed disc
        ! this routine is called after the update of disc%gsh_tab, we have to remove the latest accreted mass
        ! focus only on unstructured gas
        M1 = (1.d0 - accreted_mass%f_str)*M1
        M2 = max(0.d0,disc_mass(disc,component='unstr') - M1)
        t_disp = disc%t_dyn
        fdisp  = min(1.d0,dt/t_disp)
        Edisp2 = (1.d0-fdisp)*5.d-1*M2*disc%dV**2.
        !
        ! No gravitational energy in that case
        Eint = 0.d0
        !
        ! Add turbulence energy
        t_disp = disc%h / disc%dV
        fdisp  = min(1.d0,dt/t_disp)
        Eturb = (1.d0-fdisp)*disc%Qturb_unstr*dt
     end if
     !
     if (present(disc2)) then
        !
        ! merger case
        !
        ! gas mass of the first progenitor
        M1     = disc_mass(disc,component='unstr')
        ! kinetic energy of the first progenitor
        Edisp1 = 5.d-1*M1*disc%dV**2.
        !
        ! gas mass of the second progenitor
        M2     = disc_mass(disc2,component='unstr')
        ! kinetic energy of the second progenitor
        Edisp2 = 5.d-1*M2*disc2%dV**2.
        !
        ! We add a fraction of the gravitational interaction energy
        if (disc%rd + disc2%rd .gt. 0.d0) then
            !
            Eint = gravconst_code_unit*(M1*M2)/(1.68d0*(disc%rd + disc2%rd))
        end if
        !
        ! In this case Eturb = 0.
        Eturb = 0.d0
     end if

     if ((M1 + M2) .gt. 0.d0) then
        dV = max(1.d1*sig_star,sqrt(2.d0*(Edisp1 + Edisp2 + Eint + Eturb)/(M1 + M2)))
        disc_update_velocity_dispersion = min(9.9d-1*Vesc,dV)
        !
        if (is_NaN(disc_update_velocity_dispersion)) then
            !
            call IO_print_error_message('new dV is NAN',only_rank=rank,called_by='disc_velocity_dispersion')
            call IO_print_message('used',only_rank=rank,component='disc', &
               param_name=(/'M1                       ','M2                       ', &
                            'Edisp1                   ','Edisp2                   '/), &
               real_param_val=(/M1,M2,Edisp1,Edisp2/))
            stop ! stop the program
        end if
     end if

     return
  end function disc_update_velocity_dispersion

  !*****************************************************************************************************************

  function disc_kappa(r,dm,disc,bulge) ! INTERFACE VERSION

    ! COMPUTE THE EPICYCLIC FREQUENCY OF A ISOLATED DISC
    ! use FREEMAN+1970 eq 16

    implicit none

    real(kind=8), intent(in)     :: r             ! radius (in kpc)
    real(kind=8)                 :: disc_kappa    ! in Gyr^-1

    type(dm_type),intent(in)     :: dm            ! a dark matter component
    type(disc_type), intent(in)  :: disc          ! a disc component
    type(bulge_type), intent(in) :: bulge         ! a bulge component

    disc_kappa = disc_kappa_(r,(/dm%rho_core,dm%r_core,disc%rd,disc_mass(disc,component='no_AGN'),agn_mass(disc%agn),bulge%rb,bulge_mass(bulge)/))

    return
  end function disc_kappa

  !*****************************************************************************************************************

  function disc_kappa_(r,param)

    ! COMPUTE THE EPICYCLIC FREQUENCY OF A ISOLATED DISC
    ! use FREEMAN+1970 eq 16

    implicit none

    real(kind=8),intent(in)   :: r            ! radius (in kpc)
    real(kind=8),intent(in)   :: param(7)     ! parameter array : param(1) = dm%rho_core
                                              !                   param(2) = dm%r_core
                                              !                   param(3) = disc%rd
                                              !                   param(4) = disc_mass(disc,component='no_AGN')
                                              !                   param(5) = agn_mass(disc%agn)
                                              !                   param(6) = bulge%rb
                                              !                   param(7) = bulge_mass(bulge)
    real(kind=8)              :: disc_kappa_  ! in Gyr^-1
    real(kind=8)              :: Omega        ! The angular velocity
    real(kind=8)              :: Omega_u
    real(kind=8)              :: Omega_d
    real(kind=8)              :: dOmega_dr    ! derivative
    real(kind=8), parameter   :: dr = 1.e-5

    Omega       = disc_velocity_(r,param)/r
    Omega_u     = disc_velocity_(r+dr,param)/(r+dr)
    Omega_d     = disc_velocity_(r-dr,param)/(r-dr)
    dOmega_dr   = (Omega_u - Omega_d)/(2.d0*dr)
    disc_kappa_ = 4.d0*pi*Omega**2.*(1.d0+r/(2.d0*Omega)*dOmega_dr)
    disc_kappa_ = sqrt(disc_kappa_)  ! in Gyr^-1

    return
  end function disc_kappa_

  !*****************************************************************************************************************

  function disc_Toomre_parameter(disc)

     ! RETURN THE GLOBAL TOOMRE PARAMETER OF THE DIFFUSE DISC GAS PHASE
     ! WARNING ! This function must be called after the update of disc%dV, disc%kappa, disc%rd

     implicit none

     real(kind=8)                      :: disc_Toomre_parameter
     real(kind=8)                      :: Sigma   ! mean gas mass surface density
     real(kind=8)                      :: Mgas    ! total gas mass in the disc

     type(disc_type),intent(in)        :: disc    ! the disc component

     disc_Toomre_parameter = 1.d0  ! init

     Mgas = disc_mass(disc,r=2.d0*1.68d0*disc%rd,component='unstructured')

     if (Mgas .le. 0.d0) return  ! no disc

     Sigma = Mgas/(pi*(2.d0*1.68d0*disc%rd)**2.)

     disc_Toomre_parameter = disc%kappa*disc%dV/sqrt(3.d0)/(pi*gravconst_code_unit*Sigma)

     return
  end function disc_Toomre_parameter

  !*****************************************************************************************************************

  function kappa_sigma_r(r,param)

    ! RETURN THE PRODUCT OF THE EPICYCLIC FREQUENCY, THE MASS SURFACE DENSITY AND THE GALACTIC RADIUS
    ! this function allows to compute the mass-weighted epiclyclic frequency

    implicit none

    real(kind=8), intent(in)  :: r         ! orbital radius [code unit kpc]
    real(kind=8),intent(in)   :: param(7)  ! parameter array : param(1) = dm%rho_core
                                           !                   param(2) = dm%r_core
                                           !                   param(3) = disc%rd
                                           !                   param(4) = disc_mass(disc,component='no_AGN')
                                           !                   param(5) = agn_mass(disc%agn)
                                           !                   param(6) = bulge%rb
                                           !                   param(7) = bulge_mass(bulge)
    real(kind=8)              :: kappa_sigma_r

    kappa_sigma_r = disc_kappa_(r,param)*disc_mass_surf_density_(r,param(3:4))*r

    return
  end function kappa_sigma_r

  !*****************************************************************************************************************

  function disc_mass_weighted_kappa(r,dm,disc,bulge,error)

    ! RETURN THE MASS WEIGHTED EPICYCLIQUE FREQUENCY
    ! For mass enclose in a radius r

    implicit none

    real(kind=8),intent(in)            :: r         ! orbital radius
    real(kind=8)                       :: kappa     ! mass weighted epicyclic frequency
    real(kind=8),intent(out),optional  :: error     ! integration error
    real(kind=8)                       :: disc_mass_weighted_kappa

    type(dm_type), intent(in)          :: dm        ! a dark matter component
    type(disc_type), intent(in)        :: disc      ! a disc component
    type(bulge_type),intent(in)        :: bulge     ! a bulge component

    kappa = 2.d0*pi/(disc_mass(disc,r=r))
    kappa = kappa*Ronbint(kappa_sigma_r,r_sf,r, &
                            (/dm%rho_core,dm%r_core,disc%rd,disc_mass(disc,component='no_AGN'),agn_mass(disc%agn),bulge%rb,bulge_mass(bulge)/), &
                            error = error, &
                            called_by = 'disc_mass_weighted_kappa')

    disc_mass_weighted_kappa = kappa

    return
  end function

  !*****************************************************************************************************************

  function V_sigma_r(r,param)

    ! RETURN THE PRODUCT OF THE ORBITAL VELOCITY THE MASS SURFACE DENSITY AND THE GALACTIC RADIUS
    ! this function allows to compute the mass-weighted orbital velocity of a disc

    implicit none

    real(kind=8),intent(in)   :: r         ! orbital radius [code unit kpc]
    real(kind=8),intent(in)   :: param(7)  ! parameter array : param(1) = dm%rho_core
                                           !                   param(2) = dm%r_core
                                           !                   param(3) = disc%rd
                                           !                   param(4) = disc_mass(disc,component='no_AGN')
                                           !                   param(5) = agn_mass(disc%agn)
                                           !                   param(6) = bulge%rb
                                           !                   param(7) = bulge_mass(bulge)
    real(kind=8)              :: V_sigma_r

    V_sigma_r = disc_velocity_(r,param)*disc_mass_surf_density_(r,param(3:4))*r

    return
  end function V_sigma_r

  !*****************************************************************************************************************
  !
  ! PRINTING PROCEDURES
  !
  !*****************************************************************************************************************

  subroutine disc_load_disc_data(disc,disc_data)

    ! CREATE THE DISC OUTPUT PROPERTIES LIST

    implicit none

    real(kind=8), intent(inout) :: disc_data(nb_disc_field+nb_stars_field+2*nb_dust_field+nb_agn_field)
    real(kind=8)                :: starsd_data(nb_stars_field)
    real(kind=8)                :: agn_data(nb_agn_field)
    real(kind=8)                :: dust_ISM_data(nb_dust_field)
    real(kind=8)                :: dust_BC_data(nb_dust_field)
    real(kind=8)                :: OH_unstr, OH_str

    type(gas_type)              :: unstr_gas
    type(gas_type)              :: str_gas

    type(disc_type),intent(in)  :: disc   ! disc component

    ! For information
    ! 'disc_stars_age_MW     ','disc_stars_dage_MW    '
    ! 'disc_stars_age_LW     ','disc_stars_dage_LW    '
    ! 'disc_Z_stars_MW       ','disc_Z_stars_LW       '
    call stars_load_stars_data(disc%stars,starsd_data)

    ! For information
    ! 'disc_dust_ISM_mass    ','disc_dust_ISM_f_pah   ','disc_dust_ISM_f_bg    ', &
    ! 'disc_dust_ISM_tau     ','disc_dust_ISM_mZmH    '
    call dust_load_dust_data(disc%dust(2),dust_ISM_data)
    call dust_load_dust_data(disc%dust(1),dust_BC_data)

    ! For information
    ! 'agn_mass              ','agn_L                 ','agn_L_Edd             ',&
    ! 'agn_acc_rate_in       ','agn_ejecta_rate       '
    call agn_load_agn_data(disc%agn,agn_data)

    ! For information
    !'t_since_last_merger   ','disc_unstr_gas        ','disc_unstr_mZ         ',&
    !'disc_unstr_mH         ','disc_unstr_mC         ','disc_unstr_mN         ',&
    !'disc_unstr_mO         ','disc_unstr_mFe        ','disc_unstr_OH_index   ',&
    !'disc_str_gas          ','disc_str_mZ           ','disc_str_mH           ',&
    !'disc_str_mC           ','disc_str_mN           ','disc_str_mO           ',&
    !'disc_str_mFe          ','disc_str_OH_index     ','disc_f_str            ',&
    !'disc_ngc              ','disc_rd               ','disc_h                ',&
    !'disc_V                ','disc_dV               ','disc_Q                ',&
    !'disc_cool_timescale   ','disc_t_dyn            ','disc_t_cool           ',&
    !'disc_t_cascade        ','disc_t_form           ','disc_t_eq             ',&
    !'disc_t_emp            ','disc_gas_acc_rate     ','disc_stripping_rate   ',&
    !'disc_str_rate         ','disc_disrupt_rate     ','disc_ejecta_rate      ',&
    !'disc_sfr              ','disc_sfr_burst        '/)

    unstr_gas = mass_code_unit_in_M_Sun*disc_mass(disc,component='unstr')*disc_gas_signature(disc,component='unstr')
    str_gas   = mass_code_unit_in_M_Sun*(disc_mass(disc,component='str')*disc_gas_signature(disc,component='str'))
    OH_unstr = O_H(unstr_gas)
    OH_str   = O_H(str_gas)
    if (OH_unstr .gt. 0.d0) OH_unstr = 12.d0 + log10(OH_unstr)
    if (OH_str .gt. 0.d0)   OH_str = 12.d0 + log10(OH_str)

    disc_data = (/disc%t_since_last_merger, gas_mass(unstr_gas), &
                  gas_mass(unstr_gas,component='Metals'), &
                  gas_mass(unstr_gas,component='H1'), &
                  gas_mass(unstr_gas,component='C12'), &
                  gas_mass(unstr_gas,component='N14'), &
                  gas_mass(unstr_gas,component='O16'), &
                  gas_mass(unstr_gas,component='Fe56'), &
                  OH_unstr, &
                  gas_mass(str_gas), &
                  gas_mass(str_gas,component='Metals'), &
                  gas_mass(str_gas,component='H1'), &
                  gas_mass(str_gas,component='C12'), &
                  gas_mass(str_gas,component='N14'), &
                  gas_mass(str_gas,component='O16'), &
                  gas_mass(str_gas,component='Fe56'), &
                  OH_str, &
                  disc%f_str, real(disc%ngc,kind=8), disc%rd, disc%h, &
                  disc%V*vel_code_unit_2_kmPerSec, &
                  disc%dV*vel_code_unit_2_kmPerSec, &
                  disc%Q, disc%cooling_timescale, disc%t_dyn, disc%t_cool, &
                  disc%gsh_tab%t_cascade, disc%gsh_tab%t_form, &
                  disc%gsh_tab%t_eq, disc%gsh_tab%t_emp, &
                  mass_rate_code_unit_2_MsunPerYr*gas_mass(disc%fresh_gas_acc_rate), &
                  mass_rate_code_unit_2_MsunPerYr*gas_mass(disc%stripping_rate), &
                  mass_rate_code_unit_2_MsunPerYr*gas_mass(disc%str_rate), &
                  mass_rate_code_unit_2_MsunPerYr*gas_mass(disc%disrupt_rate), &
                  mass_rate_code_unit_2_MsunPerYr*gas_mass(disc%ejecta_rate), &
                  mass_rate_code_unit_2_MsunPerYr*gas_mass(disc%sfr), &
                  mass_rate_code_unit_2_MsunPerYr*disc%sfr_burst, &
                  starsd_data,dust_ISM_data,dust_BC_data,agn_data/)

    return
  end subroutine disc_load_disc_data

  !*****************************************************************************************************************

  subroutine disc_print(unit,form,disc)

    ! PRINT DISC-PROPERTIES IN eGALICS OUTPUT FILE

    implicit none

    integer(kind=4),intent(in)  :: unit  ! file unit
    integer(kind=4)             :: status,i,hdutype

    character(*)                :: form  ! fits or tmp_bin

    real(kind=8)                :: disc_data(nb_disc_field+nb_stars_field+2*nb_dust_field+nb_agn_field)

    type(disc_type),intent(in)  :: disc  ! disc component

    call disc_load_disc_data(disc,disc_data)

    select case (trim(form))
      case ('tmp_bin')
        write(unit) disc_data  ! directly write data in the tmp binary output file
      case ('fits')
        ! move to disc extension
        call ftmahd(unit,hdu_disc,hdutype,status)
        if (status .gt. 0) then
          !
          call IO_print_error_message('ftmahd status', &
                only_rank = rank, called_by = 'disc_print')
          stop ! stop the program
        end if
        !
        ! init
        call ftirow(unit,0,1,status)
        if (status .gt. 0) then
          !
          call IO_print_error_message('ftirow status', &
                only_rank = rank, called_by = 'disc_print')
          stop ! stop the program
        end if
        !
        ! write data in the disc entension
        call ftpcld(unit,1,1,1,1,disc%age_form+disc%life_time,status)
        do i=2, nb_disc_field+nb_stars_field+2*nb_dust_field+nb_agn_field+1
           call ftpcld(unit,i,1,1,1,disc_data(i-1),status)
           if (status .gt. 0) then
              !
              call IO_print_error_message('ftpcld status', &
                          only_rank = rank, called_by = 'disc_print')
              stop ! stop the program
           end if
        end do
      case default
        call IO_print_error_message('Unknwon output data format', &
                only_rank = rank, called_by = 'disc_print')
        stop ! stop the program
    end select

    return
  end subroutine disc_print

  !*****************************************************************************************************************

end module disc
