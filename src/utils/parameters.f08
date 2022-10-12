module parameters

    !*****************************************************************************************************************
    ! 
    ! OVERVIEW
    !
    ! Here are defined - mathematical, physical, astronomical and numerical constants (precision)
    !
    ! G.A.S. uses:     - 1     kpc as lenght scale
    !                  - 1     Gyr as time scale
    !                  - 10^11 Msun as mass scale
    !                  - 1     K    as temperature scale
    !                  - All other units can be deduced from that units
    ! 
    !*****************************************************************************************************************

    public

    ! General parameters
    integer(kind=4), parameter :: MAXPATHSIZE       = 2048
    !
    ! Mathematical constants 
    real(kind=8), parameter    :: pi                = 3.141592653589793238d0                                   ! Just pi
    real(kind=8), parameter    :: num_accuracy      = 1.e-16 
    ! 
    ! Physical constants
    real(kind=8), parameter    :: lightSpeed_m_s    = 2.99792458d8        ! Speed of light (m/s)
    real(kind=8), parameter    :: lightSpeed_km_s   = 2.99792458d5        ! Speed of light (km/s)
    real(kind=8), parameter    :: kb                = 1.380650424d-23     ! Boltzmann constant (Joule/K = kg.m^2/s^2/K)
    real(kind=8), parameter    :: mp                = 1.67262163783d-27   ! Mass of a proton (kg)
    real(kind=8), parameter    :: sigma_t           = 6.65245855827d-29   ! Thomson cross-section for electron (in m^2)
    !
    ! Astronomical constants 
    real(kind=8), parameter    :: ZSun              = 0.02d0              ! Solar metallicity (mass fraction)
    real(kind=8), parameter    :: LSun              = 3.842714d26         ! Sun's bolometric luminosity (W = Joule/s = kg.m^2/s^3)
    real(kind=8), parameter    :: LSun_erg_s        = LSun * 1.d7         ! Sun's bolometric luminosity (1 Joule = 1.e7 ergs)
    real(kind=8), parameter    :: MSun              = 1.988550d30         ! Solar mass (kg)
    real(kind=8), parameter    :: GCst              = 6.6742867d-11       ! Gravitational constant (in m^3/kg/s^2)
    !
    ! Units conversions & conversion factor 
    real(kind=8), parameter    :: Mass2Msun         = 1.d11               ! 1 mass code unit in MSun
    real(kind=8), parameter    :: MSun2mass         = 1.d0/mass2Msun      ! 1 MSun in code unit           
    real(kind=8), parameter    :: Mass_kg           = Mass2Msun*MSun      ! 1 mass code unit in kg (~ 1.989 10^41 kg)
    real(kind=8), parameter    :: kpc2m             = 3.0856776d19        ! 1 kpc in m
    real(kind=8), parameter    :: pc2cm             = 3.0856776d18        ! 1 pc in cm
    real(kind=8), parameter    :: kpc2cm            = kpc2m*1.d2          ! 1 kpc in cm 
    real(kind=8), parameter    :: kpc2km            = kpc2m*1.d-3         ! 1 kpc in km
    real(kind=8), parameter    :: Mpc2km            = kpc2km*1.d3         ! 1 Mpc in km
    real(kind=8), parameter    :: Gyr2s             = 3.1556952d16        ! 1 Gyr in sec
    real(kind=8), parameter    :: Gyr2yr            = 1.d9                ! 1 Gyr in yr
    real(kind=8), parameter    :: lightSpeed_CU     = lightSpeed_m_s*Gyr2s/kpc2m       ! Light Speed in code unit [kpc/Gyr] (~ 306601)
    real(kind=8), parameter    :: Energy_J          = Mass_kg*kpc2m**2./(Gyr2s**2.)    ! 1 CU energy (~ 1.90121 10^48 J)
    real(kind=8), parameter    :: Energy_CU         = 1.d0/Energy_J                    ! 1 Joule in CU [10^11.Msun.kpc^2/Gyr^2] (~ 5.25980 10^-48 )
    real(kind=8), parameter    :: Luminosity_W      = Mass_kg*(kpc2m)**2./(Gyr2s**3.)  ! 1 CU luminosity [10^11Msun.kpc^2/Gyr^3] (~ 6.02457 10^30 W)
    real(kind=8), parameter    :: LSun_CU           = LSun/Luminosity_W                ! 1 solar bolometric luminosity in CU [10^11.Msun.kpc^2/Gyr^3] (~ 6.3822 10^-5)
    real(kind=8), parameter    :: GCst_CU           = GCst*Mass_kg/kpc2m**3.*Gyr2s**2. ! Gravitational constant [in code unit] (kpc^3/(10^11 Msun)/Gyr^2) (~ 4.49862 10^5)
    real(kind=8), parameter    :: Velocity_km_s     = kpc2km/Gyr2s                     ! From CU [kpc/Gyr] to [km/s]
    real(kind=8), parameter    :: Velocity_m_s      = kpc2m/Gyr2s                      ! From CU [kpc/Gyr] to [m/s]
    real(kind=8), parameter    :: MassRate_Msun_Yr  = Mass2Msun/Gyr2yr                 ! From CU [10^11Msun/Gyr] to [Msun/yr]

end module parameters