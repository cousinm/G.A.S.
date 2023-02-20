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
    !                  - All other units can be deduced from these units
    ! 
    !*****************************************************************************************************************

    public

    ! General parameters
    integer, parameter           :: ikd = selected_int_kind(6)
    integer, parameter           :: rkd = selected_real_kind(12)
    integer(kind=ikd), parameter :: MAXPATHSIZE        = 2048
    !
    ! Mathematical constants 
    real(kind=rkd), parameter    :: pi                 = 3.141592653589793238d0   ! Just pi
    real(kind=rkd), parameter    :: num_accuracy       = 1.e-11
    ! 
    ! Physical constants
    real(kind=rkd), parameter    :: lightSpeed_m_s     = 2.99792458d8        ! Speed of light (m/s)
    real(kind=rkd), parameter    :: lightSpeed_km_s    = 2.99792458d5        ! Speed of light (km/s)
    real(kind=rkd), parameter    :: kb                 = 1.380650424d-23     ! Boltzmann constant (Joule/K = kg.m^2/s^2/K)
    real(kind=rkd), parameter    :: mp                 = 1.67262163783d-27   ! Mass of a proton (kg)
    real(kind=rkd), parameter    :: sigma_t            = 6.65245855827d-29   ! Thomson cross-section for electron (in m^2)
    !
    ! Astronomical constants 
    real(kind=rkd), parameter    :: ZSun               = 0.02d0              ! Solar metallicity (mass fraction)
    real(kind=rkd), parameter    :: LSun               = 3.842714d26         ! Sun's bolometric luminosity (W = Joule/s = kg.m^2/s^3)
    real(kind=rkd), parameter    :: LSun_erg_s         = LSun * 1.d7         ! Sun's bolometric luminosity (1 Joule = 1.e7 ergs)
    real(kind=rkd), parameter    :: MSun               = 1.988550d30         ! Solar mass (kg)
    real(kind=rkd), parameter    :: GCst               = 6.6742867d-11       ! Gravitational constant (in m^3/kg/s^2)
    !
    ! Units conversions & conversion factor
    ! Masses
    real(kind=rkd), parameter    :: MassCU2MSun        = 1.d11               ! 1 mass code unit (MassCU) in MSun [Msun/MassCU]
    real(kind=rkd), parameter    :: MSun2MassCU        = 1.d0/MassCU2Msun    ! 1 MSun in MCU [MassCU/Msun]
    real(kind=rkd), parameter    :: MassCU2kg          = MassCU2Msun*MSun    ! 1 MassCU in kg (~ 1.989 10^41 kg) [kg/MassCU]
    real(kind=rkd), parameter    :: kg2MassCU          = 1.d0/MassCU2kg      ! 1 kg in MassCU [MassCU/kg]
    ! Distances
    real(kind=rkd), parameter    :: kpc2m              = 3.0856776d19        ! 1 kpc in m [m/kpc]
    real(kind=rkd), parameter    :: LenCU2m            = kpc2m               ! 1 lenght code unit (LenCU) in m [m/LenCU]
    real(kind=rkd), parameter    :: m2LenCU            = 1.d0/LenCU2m        ! 1 m to LenCU in m [LenCU/m]
    real(kind=rkd), parameter    :: LenCU2cm           = LenCU2m*1.d2        ! 1 kpc in cm [cm/kpc]
    real(kind=rkd), parameter    :: cm2LenCU           = 1.d0/LenCU2cm       ! 1 cm in LenCU [lenCU/cm]
    real(kind=rkd), parameter    :: LenCU2km           = LenCU2m*1.d-3       ! 1 LenCU in km [km/LenCU]
    real(kind=rkd), parameter    :: km2LenCU           = 1.d0/LenCU2km       ! 1 km in LenCU [LenCU/km]
    real(kind=rkd), parameter    :: pc2cm              = 3.0856776d18        ! 1 pc in cm [cm/pc]
    real(kind=rkd), parameter    :: pc2LenCU           = 1.d-3               ! 1 pc in LenCU [LCU/pc]
    real(kind=rkd), parameter    :: Mpc2km             = LenCU2km*1.d3       ! 1 Mpc in km [km/Mpc]
    ! Times
    real(kind=rkd), parameter    :: Gyr2s              = 3.1556952d16        ! 1 Gyr in sec [s/Gyr]
    real(kind=rkd), parameter    :: TimeCU2s           = Gyr2s               ! 1 time code unit (TimeCU) in sec [s/TCU]
    real(kind=rkd), parameter    :: s2TimeCU           = 1.d0/TimeCU2s       ! 1 second in TimeCU [TimeCU/s]
    real(kind=rkd), parameter    :: TimeCU2yr          = 1.d9                ! 1 TimeCU in yr [yr/TimeCU]
    ! Velocity
    real(kind=rkd), parameter    :: km_s2VelCU         = km2LenCU/s2TimeCU                   ! 1 km/s in velocity code unit VCU [VCU/km_s]
    real(kind=rkd), parameter    :: m_s2VelCU          = m2LenCU/s2TimeCU                    ! 1 m in VCU [VCU/m_s]
    real(kind=rkd), parameter    :: lightSpeed_VelCU   = lightSpeed_m_s*m2LenCU/s2TimeCU     ! Light Speed (m/s) in code unit (~ 306601)
    real(kind=rkd), parameter    :: VelCU2km_s         = 1.d0/km_s2VelCU                     ! 1 VCU in km/s [km_s/VCU]
    real(kind=rkd), parameter    :: VelCU2m_s          = 1.d0/m_s2VelCU                      ! 1 VCU in m/s [m_s/VCU]
    ! Ernergy
    real(kind=rkd), parameter    :: J2EnergCU          = kg2MassCU*m_s2VelCU**2              ! 1 J (kg.m^2/s^2) in energy code unit ECU (ECU/kg.m2_s2)
    real(kind=rkd), parameter    :: EnergCU2J          = 1.d0/J2EnergCU                      ! 1 ECU in J
    real(kind=rkd), parameter    :: erg2J              = 1.d-7                               ! 1 erg in Joule [erg/J]
    real(kind=rkd), parameter    :: erg2EnergCU        = erg2J*kg2MassCU*m_s2VelCU**2        ! 1 erg (kg.m^2/s^2) in energy code unit ECU (ECU/kg.m2_s2)
    real(kind=rkd), parameter    :: EnergCU2erg        = 1.d0/erg2EnergCU                    ! 1 ECU in erg
    ! Power / Luminosity
    real(kind=rkd), parameter    :: Watt2LumCU         = kg2MassCU*m2LenCU**2./s2TimeCU**3.  ! 1 Watt in luminosity code unit LumCU [10^11Msun.kpc^2/Gyr^3]
    real(kind=rkd), parameter    :: LSun2LumCU         = LSun*Watt2LumCU                     ! 1 solar bolometric luminosity in LumCU [10^11.Msun.kpc^2/Gyr^3]
    ! Gravitationnal constant
    real(kind=rkd), parameter    :: GCst_CU            = GCst/kg2MassCU*m2LenCU**3./s2TimeCU**2  ! Gravitational constant [in code unit] (kpc^3/(10^11 Msun)/Gyr^2) (~ 4.49862 10^5)
    ! Mass rate
    real(kind=rkd), parameter    :: MassRateCU2Msun_Yr = MassCU2MSun/TimeCU2yr                   ! From CU [10^11Msun/Gyr] to [Msun/yr]
    real(kind=rkd), parameter    :: Msun_Yr2MassRateCU = 1.d0/MassRateCU2Msun_Yr                 ! 1 Msun/yr in CU

end module parameters
