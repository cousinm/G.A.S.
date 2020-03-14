module global_variables

  !*****************************************************************************************************************
  ! 
  ! OVERVIEW
  !
  ! Here are defined - cosmological, astronomical, physical, mathematical and numerical constants (precision)
  !                  - the model parameters   
  !
  ! G.A.S. uses:     - 1     kpc as lenght scale
  !                  - 1     Gyr as time scale
  !                  - 10^11 Msun as mass scale
  !                  - 1     K    as temperature scale
  !                  - All other units can be deduced from that units
  !
  ! NO SUBROUTINE IN THIS MODULE
  !
  ! NO FUNCTION IN THIS MODULE
  ! 
  !*****************************************************************************************************************
                
  public

  integer(kind=4),parameter :: MAXPATHSIZE                     = 512                                                      ! Maximum length of a path/file name
  !
  ! mathematical constants 
  real(kind=8),parameter    :: pi                              = 3.141592653589793238d0                                   ! Just pi
  ! 
  ! physical constants
  real(kind=8),parameter    :: light_speed_mPerSec             = 2.99792458d8                                             ! speed of light (m/s)
  real(kind=8),parameter    :: light_speed_kmPerSec            = 2.99792458d5                                             ! speed of light (km/s)
  real(kind=8),parameter    :: kb                              = 1.380650424d-23                                          ! Boltzmann constant (Joule/K = kg.m^2/s^2/K)
  real(kind=8),parameter    :: mp                              = 1.67262163783d-27                                        ! mass of a proton (kg)
  real(kind=8),parameter    :: sigma_t                         = 6.65245855827d-29                                        ! Thomson cross-section for electron (in m^2)
  ! astronomical constants 
  real(kind=8),parameter    :: Z_Sun                           = 0.02d0                                                   ! Solar metallicity (mass fraction)
  real(kind=8),parameter    :: L_Sun                           = 3.842714d26                                              ! Sun's bolometric luminosity (W = Joule/s = kg.m^2/s^3)
  real(kind=8),parameter    :: L_Sun_erg_s                     = L_Sun * 1.d7                                             ! Sun's bolometric luminosity (1 Joule = 1.e7 ergs)
  real(kind=8),parameter    :: M_Sun                           = 1.988550d30                                              ! (kg)
  real(kind=8),parameter    :: gravconst                       = 6.6742867d-11                                            ! Gravitational constant (in m^3/kg/s^2)
  !
  ! units conversions & conversion factor 
  real(kind=8),parameter    :: mass_code_unit_in_M_Sun         = 1.d11                                                    ! 1 mass code unit = 10^11 M_sun
  real(kind=8),parameter    :: M_Sun_in_mass_code_unit         = 1.d0/mass_code_unit_in_M_Sun                             ! 1 M_sun in code unit           
  real(kind=8),parameter    :: mass_code_unit_in_kg            = 1.d11*M_Sun                                              ! 1 mass code unit (10^11 Msun) in kg (number of kg in 10^11 Msun) (~ 1.989 10^41 kg)
  real(kind=8),parameter    :: kpc_in_m                        = 3.0856776d19                                             ! number of meters in 1 kpc
  real(kind=8),parameter    :: pc_in_cm                        = 3.0856776d18                                             ! number of centimeter in 1 pc
  real(kind=8),parameter    :: kpc_in_cm                       = pc_in_cm*1.d3
  real(kind=8),parameter    :: kpc_in_km                       = kpc_in_m*1.d-3                                           ! number of kilometer in 1 kpc  
  real(kind=8),parameter    :: Mpc_in_km                       = kpc_in_m                                                 ! number of kilometer in 1 Mpc 
  real(kind=8),parameter    :: Gyr_in_s                        = 3.1556952d16                                             ! number of second in 1 Gyr
  real(kind=8),parameter    :: Gyr_in_yr                       = 1.d9                                                     ! number of yr in 1 Gyr
  real(kind=8),parameter    :: light_speed_code_unit           = light_speed_mPerSec*Gyr_in_s/kpc_in_m                    ! in [kpc/Gyr] (~ 306601)
  real(kind=8),parameter    :: J_in_code_unit                  = Gyr_in_s**2./ mass_code_unit_in_kg/(kpc_in_m**2.)        ! 1 joule in code unit (~ 5.25980 10^-48)
  real(kind=8),parameter    :: E_code_unit_in_J                = 1.d0/J_in_code_unit                                      ! 1 code unit energy (~ 1.90121 10^48 J)
  real(kind=8),parameter    :: L_code_unit_in_W                = mass_code_unit_in_kg*(kpc_in_m)**2./(Gyr_in_s**3.)       ! 1 code unit luminosity [10^11.Msun.kpc^2/Gyr^3] (~ 6.02457 10^30 W)
  real(kind=8),parameter    :: L_Sun_in_code_unit              = L_sun/L_code_unit_in_W                                   ! 1 solar bolometric luminosity in code unit [10^11.Msun.kpc^2/Gyr^3] (~ 6.3822 10^-5)
  real(kind=8),parameter    :: gravconst_code_unit             = gravconst*mass_code_unit_in_kg/kpc_in_m**3.*Gyr_in_s**2. ! Gravitational constant [in code unit] (kpc^3/(10^11 Msun)/Gyr^2) : (~ 4.49862 10^5)   
  real(kind=8),parameter    :: vel_code_unit_2_kmPerSec        = kpc_in_km/Gyr_in_s                                       ! from kpc/Gyr --> km/s
  real(kind=8),parameter    :: vel_code_unit_2_mPerSec         = kpc_in_m/Gyr_in_s                                        ! from kpc/Gyr --> m/s
  real(kind=8),parameter    :: mass_rate_code_unit_2_MsunPerYr = mass_code_unit_in_M_Sun/Gyr_in_yr                        ! from 10^11Msun/Gyr --> Msun/yr
  ! 
  ! cosmological parameters
  real(kind=8)              :: h_0                                        ! [km/s/Mpc]
  real(kind=8)              :: h_0_code_unit                              ! [Gyr]
  real(kind=8)              :: Omega_L                                    ! Dark energy, fixed during initialization (red_cosmology)
  real(kind=8)              :: Omega_m                                    ! Matter, fixed during initialization (red_cosmology)
  real(kind=8)              :: Omega_b                                    ! Baryons, fixed during initialization (red_cosmology)
  real(kind=8)              :: baryon_fraction                            ! Baryonic fraction, fixed during initialization (red_cosmology)
  real(kind=8)              :: rho_crit_code_unit
  !
  ! dm-simulation parameters        
  real(kind=8)              :: L_box                                      ! kpc (code unit)
  real(kind=8)              :: M_tot_min                                  ! minimal halo mass used in the halo finder process (= 20 * dm_particle_mass)   
  real(kind=8)              :: dm_particle_mass                           ! mass of a dark-matter particle in the N-body simulation (depend of the resolution)
  !                                                                       ! fixed during initialization (IO_read_parameter_file)   
  ! model parameters 
  !
  real(kind=8)              :: physical_precision                         ! given in input
  !
  ! cleaning process parameters
  integer(kind=4)           :: nb_of_desc_halo_test                       ! If the cleaning process is turn ON, 
                                                                          ! The cleaning algorithm checks various properties of halos in the next 'nb_of_desc_halo_test' halos.  
  real(kind=4)              :: min_halo_goodness_of_branch                ! goodness minimal value for keep the branch 
  real(kind=4)              :: min_dm_goodness_of_branch                  ! goodness minimal value for keep the branch (only dm-properties are check) 
  real(kind=8)              :: Ecut                                       ! halo is condidered stable (virialize) if Ek < |Ecut*Ep|
  !     
  ! background
  real(kind=8),parameter    :: T_reion                         = 3.d4     ! Valageas & Silk 1999 [K]
  real(kind=8)              :: z_overlap                                  ! The redshift at which the first H II regions begin to overlap (11) = 1. + z_reion
  real(kind=8)              :: z_reion                                    ! The redshift at which most of the medium is re-ionized (10), given in input
  real(kind=8)              :: cold_stream_efficiency                     ! given in input
  real(kind=8)              :: cooling_efficiency                         ! given in input
  real(kind=8)              :: TI_efficiency                              ! given in input
  !
  ! galaxy
  real(kind=8)              :: epsilon_merge                              ! critical mass ratio that separates major and minor mergers
  ! stars
  character(MAXPATHSIZE)    :: IMF                                        ! given in input (kennicutt, scalo, salpeter, TH2, TH4)
  ! disc 
  real(kind=8),parameter    :: Q_crit                          = 1.d0     ! without unit  FIXED HERE
  real(kind=8)              :: disc_stripping_efficiency                  ! given in input
  real(kind=8)              :: disc_ejecta_efficiency                     ! given in input
  ! agn
  real(kind=8)              :: agn_ejecta_efficiency                      ! given in input
  real(kind=8)              :: M_BH_min                                   ! given in input
  real(kind=8)              :: agn_gas_coupling                           ! given in input
  !
  ! END model parameters
  !
  ! numerical parameter
  real(kind=8),parameter    :: x_cut                           = 1.d-4    ! radius parameter ratio (used to gas profil integration)
  real(kind=8),parameter    :: num_precision                   = 1.d-8    ! numerical precision
  real(kind=8),parameter    :: num_accuracy                    = 1.d-12   ! numerical accuracy
  !
  ! other global parameters 
  real(kind=8)              :: dt_min_use
  
end module global_variables
