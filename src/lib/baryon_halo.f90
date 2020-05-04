module baryon_halo
  
  
  use dm      ! Contains dm structure definition and dm exploitation functions
  use cooling ! Contains cooling library (e.g. acces to Lambda(T,Z))

  public

 !*****************************************************************************************************************
 ! 
 !  OVERVIEW
 !
 !  baryon_halo_module defines the baryon_halo data structures univ(ts)%halo(ih)%baryon_halo 
 !
 !  This module contains properties and prescriptions used to describe the baryonic phase associated to a halo (without the galaxy)
 !  The baryonic gas phase is composed of a cold and a hot component. The cold (10**4 K) gas can be directly transfered to the galaxy
 !  The hot gas phase has to follow a radiative cooling process, to cool and feed the galaxy disc 
 !  At a given time only a fraction of the hot gas can condensate and feed the galaxy.
 !  This fraction is defined according to the cooling time equation R(t_cool) = T_cool
 !  This fraction is strongly impacted by the thermal instability process that limit the region in which the gas can effectivelly cool down
 !  The hot gas phase is described by an analytical density profile. 
 !  In the header of the module are defined output properties of a baryon_halo component (labels, units and formats)
 !
 !  SUBROUTINES IN THIS MODULE
 !
 !   baryon_halo_set_reference_mass                      : initialize baryon_halo parameters and call subroutines for lower levels
 !
 !   baryon_halo_void                                    : initialize a baryon halo structure
 !      called by : halo_void      
 !  
 !   baryon_halo_copy                                    : copy a baryon halo component in an other
 !
 !   baryon_halo_set_age_form                            : set the formation age of a baryon_halo component
 !      called by : halo_evolve
 !
 !   baryon_halo_add_surrounding_gas                     : add gas to the surrounding halo gas phase
 !      called by : 
 ! 
 !   baryon_halo_sub_surrounding_gas                     : substract gas to the surrounding halo gas phase
 !
 !   baryon_halo_add_transfer_mass                       : add mass to the mass reservoir dedicated to mass transfer between sub and main halo
 !      called by : halo_evolve 
 !
 !   baryon_halo_transfer_hot_mass                       : transfer gas, coming from substructures, to the hot gas atmosphere
 !      called by : halo_evolve
 !
 !   baryon_halo_compute_accretion_rate                  : computed and set baryonic accretion rate (cold and hot)
 !      called by  : halo_evolve
 !      contains   : filtering_mass                      : compute the filtering mass (Okamoto+2008 or Gnedin+2000)
 !                   baryonic_fraction_photoionisation   : compute the effective baryonic fraction 
 ! 
 !   baryon_halo_compute_surrounding_gas_outflow_rate    : compute outflow rate assocaited to the hot surrounding gas phase
 !       called by  : baryon_halo_compute_accretion_rate 
 !
 !   baryon_halo_compute_galaxy_fresh_gas_accretion_rate : computed and set the effective galaxy accretion rate
 !       called by  : 
 !
 !   baryon_halo_set_escape_fraction                     : set the escape hot gas fraction (the fraction of hot mass that leave the hot atmosphere)
 !      called by : halo_evolve
 ! 
 !   baryon_halo_evolve_cold_gas_I                       : first part (PREDICTOR) of the baryon-halo cold gas evolution scheme 
 !      called by : halo_evolve
 !
 !   baryon_halo_compute_cold_outflow_rate               : compute the free-fall rate of the cold halo gas phase   
 !      called by : halo_evolve
 !
 !   baryon_halo_evolve_cold_gas_II                      : evolve the baryon-halo cold gas phase during dt
 !      called by : halo_evolve
 ! 
 !  baryon_halo_evolve_hot_gas_I                         : first part (PREDICTOR) of the baryon-halo hot gas evolution scheme   
 !     called by : halo_evolve
 ! 
 !   baryon_halo_compute_cooling_rate                    : compute the cooling rate of the hot halo gas phase
 !      called by : halo_evolve
 !
 !   baryon_halo_evolve_hot_gas_II                       : evolve the hot gas phase during dt, compute temperatures, escape mass ...
 !      called by : halo_evolve
 !
 !   baryon_halo_update_density_profile                  : update hot atmosphere density profile parameters, r_core, rho_core
 !      called by : baryon_halo_evolve_hot_gas_II
 !                 baryon_halo_merge 
 !
 !   baryon_halo_update_cooling_clock
 !      called by : baryon_halo_evolve_hot_gas           : add the last ddt to the cooling time to take account the last time loop
 !
 !   baryon_halo_update_life_time                        : add the last ddt to the baryon halo time life to take account the last time loop
 !      called by : halo_evolve 
 !   
 !   baryon_halo_merge                                   : merge two baryon halo components
 !      called by : halo_merge 
 ! 
 !  FUNCTIONS IN THIS MODULE
 !
 !   baryon_halo_escape_fraction                         : return the hot atmosphere escape fraction
 !
 !   baryon_halo_hot_gas_fraction                        : return fraction of hot gas in the baryon halo structure, m_hot / (m_hot + m_cold)
 !
 !   baryon_halo_bh_mass                                 : return the baryon halo mass (a selection for cold, hot, igm and hot+cold only are possible) 
 !
 !   baryon_halo_shock_heated_fraction                   : return the effective fraction of hot gas in the global accretion process
 !
 !   f_hot_gas                                           : the geometrical function associated to the hot gas density profile, f_hot_gas is defined by rho = rho(r=0)*f_hot_gas
 !                                                                based on Capelo+10, Mokino+98, Suto+98
 !
 !   f                                                   : a geometrical function used in the f_hot_gas function definition 
 !
 !   b_cont                                              : return the value of the parameter b_cont used in the definition of the hot gas density profile
 ! 
 !   M                                                   : a geometrical function, based on f_hot_gas, used to compute hot mass(r) and the mean hot temperature
 !  
 !   baryon_halo_hot_halo_density_profile                : return the density of the hot gas atmosphere at a given radius
 !
 !   baryon_halo_hot_halo_pressure_profile               : return the pression of the hot gas atmosphere at a given radius
 !
 !   baryon_halo_hot_halo_temperature_profile            : return the temperature of the hot gas atmosphere at a given radius
 !
 !   baryon_halo_integrate_hot_halo_density_profile      : a integration procedure of the hot atmosphere density profile that
 !                                                                returns the mass of the hot halo profile enclosed between r_in and r_rout
 !
 !   T0_crit_polytropic                                  : return the critical(minimum) central temperature of the hot atmosphere
 !
 !   T_mean_polytropic                                   : return the mean temperature of the hot atmosphere, described by a polytropic equation of state
 ! 
 !   T0_polytropic                                       : return the central temperature of the hot atmosphere, described by a polytropic equation of state
 ! 
 !   baryon_halo_bh_mass                                 : return the total mass of a baryon halo component (cold + hot) 
 !
 ! PRINTING PROCEDURES
 !
 !   bh_load_bh_data                                     : create the baryon halo output properties list
 !
 !   baryon_halo_print                                   : print baryon_halo properties in a binary file
 !      called by : halo_print
 ! 
 !*****************************************************************************************************************

 ! DEFINITION OF GLOBAL VARIABLES LINKED TO THE BARYON_HALO STRUCTURE *******************

  type baryon_halo_type  
    ! Evolution
    real(kind=8)             :: life_time                  ! effective evolution time of the baryon halo component
    real(kind=8)             :: age_form                   ! age of the universe when the structure has been formed
    logical                  :: quenched                   ! .TRUE if the halo is quenched = no more cooling
    !
    ! gas reservoirs
    type(gas_type)           :: cold_gas                   ! gas in cold flows
    type(gas_type)           :: hot_gas                    ! gas in the hot atmosphere 
    type(gas_type)           :: surrounding_gas            ! gas stored in the neighborhood of the halo 
                                                           !   (a tampon reservoir for the non accreted photo-ionized gas, 
                                                           !    the gas ejected from the hot halo is also strored into this reservoir)
    ! transfer rates
    type(gas_type)           :: cold_inflow_rate           ! background cold accretion onto the halo
    type(gas_type)           :: hot_inflow_rate            ! background hot accretion onto the halo (shock heated gas)
    type(gas_type)           :: cold_outflow_rate          ! cold accretion rate onto the galaxy
    type(gas_type)           :: cooling_rate               ! cooling rate onto the galaxy
    !  
    ! baryon halo hot atmosphere properties
    real(kind=8)             :: t_cool                     ! effective cooling time of the hot atmosphere (in Gyr)
    real(kind=8)             :: t_TI                       ! thermal instability time scale (in Gyr)
    real(kind=8)             :: r_cool                     ! cooling radius (in Mpc) according to t(r_cool) = t_cool
    real(kind=8)             :: r_TI                       ! thermal instabilitiy radius
    real(kind=8)             :: T0                         ! central temperature of the hot atmosphere
                                                           !   (in the case of a perfect isothermal gas, 
                                                           !    the central temperature is equal to the mean tempeature of the hot gas) 
                                                           ! it is not true in the case of a polytropic gas
    real(kind=8)             :: rho0                       ! central density  [10^11Msun/kpc^3]  
    real(kind=8)             :: escape_fraction            ! escape fraction, fraction (in mass) of gas ejected from the hot atmomsphere 
  end type baryon_halo_type

  ! hdu reference for bh structure
  integer(kind=4)           :: hdu_bh
  ! printable properties for bh structure
  integer(kind=4),parameter :: nb_bh_field = 18 ! number of dm properties saved
  ! Name of each output colomn data
  character(len=ttype_len),dimension(nb_bh_field) :: ttype_bh = (/'cold_gas              ','surrounding_gas       ','hot_gas               ',&
                                                                  'hot_gas_mZ            ','hot_gas_mH            ','hot_gas_mC            ',&
                                                                  'hot_gas_mN            ','hot_gas_mO            ','hot_gas_mFe           ',&
                                                                  'acc_rate              ','ff_rate               ','cooling_rate          ',&
                                                                  'T_hot                 ','t_cool                ','t_TI                  ',&
                                                                  'r_cool                ','r_TI                  ','esc_frac              '/)  
  ! Physical unit of each column data
  character(len=tunit_len),dimension(nb_bh_field) :: tunit_bh = (/'M_sun       ','M_sun       ','M_sun       ','M_sun       ','M_sun       ',&
                                                                  'M_sun       ','M_sun       ','M_sun       ','M_sun       ','M_sun/yr    ',&
                                                                  'M_sun/yr    ','M_sun/yr    ','K           ','Gyr         ','Gyr         ',&
                                                                  'kpc         ','kpc         ','w_o_unit    '/)
  ! Data type of each column data
  character(len=tform_len),dimension(nb_bh_field) :: tform_bh = (/'1E','1E','1E','1E','1E','1E','1E','1E','1E', &
                                                                  '1E','1E','1E','1E','1E','1E','1E','1E','1E'/)

contains

  !*****************************************************************************************************************
  !
  ! SUBROUTINES
  !
  !*****************************************************************************************************************

  subroutine baryon_halo_set_reference_mass 

    ! INITIALIZE REFERENCE PARAMETERS LINKED TO THE DARK-MATTER N-BODY SIMULATION

    implicit none

    call gas_set_reference_mass   ! initialize M_gas_min

    return
  end subroutine

  !*****************************************************************************************************************

  subroutine baryon_halo_void(bh)
    
    ! INIT OR VOID A BARYON_HALO STRUCTURE

    implicit none 

    type(baryon_halo_type),intent(inout) :: bh
    
    ! Evolution 
    bh%life_time          = 0.d0          ! the life time of the baryon_halo composant
    bh%age_form           = 0.d0          ! set in halo_evolve at the begining of the first evolution step
    bh%quenched           = .false.       ! by default the halo can feed the galaxy
    !
    ! gas reservoirs
    call gas_void(bh%cold_gas)            ! void the cold gas component
    call gas_void(bh%hot_gas)             ! void the hot gas component
    call gas_void(bh%surrounding_gas)     ! void the surrounding gas component
    !
    ! transfer rates
    call gas_void(bh%cold_inflow_rate)    ! init the background cold accretion onto the halo
    call gas_void(bh%hot_inflow_rate)     ! init the background hot accretion onto the halo
    call gas_void(bh%cold_outflow_rate)   ! init the cold accretion rate onto the galaxy
    call gas_void(bh%cooling_rate)        ! init the cooling rate onto the galaxy
    !
    ! baryon halo hot atmosphere properties
    bh%t_cool             = 0.d0           
    bh%t_TI               = 0.d0        
    bh%r_cool             = 0.d0
    bh%r_TI               = 0.d0    
    bh%T0                 = -1.d0         ! -1.d0: not already computed, -2.d0: unstable hot gas phase (in the polytropic case) 
    bh%rho0               = -1.d0         ! -1.d0: not already computed (two few gas in the reservoir)
    bh%escape_fraction    = 1.d0
      
    return
  end subroutine baryon_halo_void

  !*****************************************************************************************************************
  
  subroutine baryon_halo_copy(bh1,bh2)
    
    ! COPY PROPERTIES OF A BARYON_HALO STRUCTURE IN AN OTHER (copy bh2 in bh1)

    implicit none

    type(baryon_halo_type),intent(inout) :: bh1       ! a baryon halo component
    type(baryon_halo_type),intent(in)    :: bh2       ! an other one    
   
    ! Evolution 
    bh1%life_time          = bh2%life_time            ! copy the time life
    bh1%age_form           = bh2%age_form             ! copy the formation age (age of the universe when the structure has been formed)
    bh1%quenched           = bh2%quenched             ! copy the quenched state of the baryonic phase
    !
    ! gas reservoirs
    bh1%cold_gas           = bh2%cold_gas             ! copy the cold gas component
    bh1%hot_gas            = bh2%hot_gas              ! copy the hot gas component
    bh1%surrounding_gas    = bh2%surrounding_gas      ! copy the surrounding gas component
    !
    ! transfer rates
    bh1%cold_inflow_rate   = bh2%cold_inflow_rate     ! copy the background cold accretion onto the halo
    bh1%hot_inflow_rate    = bh2%hot_inflow_rate      ! copy the background hot accretion onto the halo
    bh1%cold_outflow_rate  = bh2%cold_outflow_rate    ! copy the cold accretion rate onto the galaxy
    bh1%cooling_rate       = bh2%cooling_rate         ! copy the cooling rate onto the galaxy
    !
    ! baryon halo hot atmosphere properties
    bh1%t_cool              = bh2%t_cool              ! copy the effective cooling time
    bh1%t_TI                = bh2%t_TI                ! copy the thermal_instability time scale
    bh1%r_cool              = bh2%r_cool              ! copy the cooling radius 
    bh1%r_TI                = bh2%r_TI                ! copy the cooling radius 
    bh1%T0                  = bh2%T0                  ! copy the central temperature of the hot atmosphere
    bh1%rho0                = bh2%rho0                ! copy the central density
    bh1%escape_fraction     = bh2%escape_fraction     ! copy the escape hot gas fraction
    
    return   
  end subroutine baryon_halo_copy

  !*****************************************************************************************************************
  
  subroutine baryon_halo_set_age_form(bh,age_univ)
  
    ! SET THE FORMATION AGE OF THE BARYON HALO COMPONENT
  
    implicit none
    
    real(kind=8),intent(in)               :: age_univ ! age of the universe when the galaxy is formed
    
    type(baryon_halo_type),intent(inout)  :: bh       ! a baryon halo component
    
    bh%age_form = age_univ
    
    return
  end subroutine baryon_halo_set_age_form
  
  !*****************************************************************************************************************

  subroutine baryon_halo_add_surrounding_gas(bh,gas)

    ! ADD GAS TO THE SURROUNDING HALO GAS PHASE

    implicit none

    type(baryon_halo_type),intent(inout)  :: bh    ! a baryon halo component
    type(gas_type),intent(in)             :: gas   ! hot gas

    call gas_add(bh%surrounding_gas,gas)
 
    return
  end subroutine baryon_halo_add_surrounding_gas
  
  !*****************************************************************************************************************

  subroutine baryon_halo_sub_surrounding_gas(bh,gas,called_by)

    ! SUBSTRACT GAS TO THE SURROUNDING HALO GAS PHASE

    implicit none
 
    character(MAXPATHSIZE)                :: message     ! a message
    
    real(kind=8)                          :: m1,m2
    
    type(baryon_halo_type),intent(inout)  :: bh             ! a baryon halo component
    type(gas_type),intent(in)             :: gas            ! a gas component
    
    character(*),intent(in),optional      :: called_by      ! name of the function that calls this function

    m1 = gas_mass(bh%surrounding_gas)
    m2 = gas_mass(gas)
    
    if ((m2 .gt. m1) .and. (abs(m1 - m2) .gt. num_accuracy)) then
      !
      if (present(called_by)) then
         !
         write(message,'(a,a)') 'baryon_halo_sub_surrounding_gas, called by ', trim(called_by)
      else
         !
         write(message,'(a)') 'baryon_halo_sub_surrounding_gas'
      endif 
      call IO_print_error_message('Substract too much gas',only_rank=rank,called_by=trim(message))  
      call IO_print_message('used',only_rank=rank,component='bh', &
                      param_name=(/'m1                       ','m2                       '/), &
                      real_param_val=(/m1,m2/))
      stop ! stop the program
    endif
     
    ! substract the gas        
    call gas_sub(bh%surrounding_gas,gas,called_by='baryon_halo_sub_surrounding_gas')

    return
  end subroutine baryon_halo_sub_surrounding_gas
  
  !*****************************************************************************************************************
  
  subroutine baryon_halo_transfer_hot_mass(bh,gas)
  
    ! ADD MASS, COMING FROM SUB-STRUCTURES, TO THE HOT GAS ATMOSPHERE OF THE MAIN HOST HALO

    implicit none

    type(gas_type),intent(in)              :: gas ! the transfered gas 
    type(gas_type)                         :: g   ! a local copy
    type(baryon_halo_type),intent(inout)   :: bh  ! the baryon halo component

    if (gas_mass(gas) .le. 0.d0) return           ! no mass to transfer
    ! 
    call gas_copy(g,gas)
    if (gas_temp(g) .le. 0.d0) then
        !
        ! In the case of a transfer from a halo to its host, the temperature of the transferred gas is not know 
        ! The temperature of the stripped gas is fixed to the temperature of the host hot atmosphere
        call gas_set_component(g,max(gas_temp(bh%hot_gas),diffuse_gas_temp),component='Temp')
    end if  
    !
    ! Add the gas
    call gas_add(bh%hot_gas,g) 
    call gas_void(g)

    return
  end subroutine baryon_halo_transfer_hot_mass
  
  !*****************************************************************************************************************

  subroutine baryon_halo_compute_accretion_rate(bh,dm,z,dt,post_merger)

    ! COMPUTE AND SET THE ACCRETION RATE OF THE BARYONIC HALO PHASE

    implicit none

    logical,intent(in)       :: post_merger                   ! =.true. if the evolution is done onto the remnent structure of a merger
    !
    real(kind=8),intent(in)  :: z                             ! until what redshift the halo must grow 
    real(kind=8),intent(in)  :: dt                            ! time that the halo should grow 
    real(kind=8)             :: M_end                         ! The dark-matter halo mass at the and of the evolution time
    real(kind=8)             :: dm_acc_rate                   ! the dark-matter (smooth) accretion rate
    real(kind=8)             :: f_hot                         ! fraction of the accretion in the hot phase
    real(kind=8)             :: bf                            ! the baryonic fraction (computed by taking into account the photo-ionization process)                    
    real(kind=8)             :: f_sh                          ! the shock heated fraction  
    real(kind=8)             :: M_cold, M_hot                 ! mass of cold gas (target), mass of hot gas
    !
    type(gas_type)           :: acc_rate                      ! accretion rate (hot + cold)
    type(gas_type)           :: dM_hot                        ! a fraction of gas mass that in instantaneously converted in hot gas
    type(gas_type)           :: surrounding_gas_inflow_rate   ! inflow rate (to the surrounding gas phase)
    type(gas_type)           :: surrounding_gas_outflow_rate  ! outflow rate (from the surrounding gas phase)
    !
    type(baryon_halo_type),intent(inout) :: bh                ! a baryon halo component
    type(dm_type),intent(in)             :: dm                ! a dark matter halo component
    
#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('baryon_halo_compute_accretion_rate',component='bh')
! -------------------------------------------------
#endif
  
    ! set the baryon accretion rate
    if (post_merger) then     
       !                           
       dm_acc_rate = dm%acc_rate_left                           
       M_end = dm%M_acc  
    else  
       !
       ! Mass from which the halo should grow 
       ! + accreted mass during dt
       dm_acc_rate = dm%acc_rate_right                     ! constant dark matter accretion rate during the post-merger time-step
       M_end = dm%M_acc + dm_acc_rate * dt                 ! halo mass at the merging time
    end if
    !
    call gas_void(bh%cold_inflow_rate)                     ! init
    call gas_void(bh%hot_inflow_rate)                      ! init
    call gas_void(acc_rate)                                ! init
    call gas_void(surrounding_gas_inflow_rate)             ! init
    call gas_void(surrounding_gas_outflow_rate)            ! init
    !
    ! We assume that the fraction of hot .vs. cold gas follows the Lu et al. prescription
    ! We divide in two parts the baryonic accretion flow
    ! in addition, if at a given time the fraction of cold gas is too high,
    ! we convert (instantaneously) a fraction of the cold gas in hot gas
    ! we assume that, when the hot atmosphere is developping, a fraction of cold stream is disrupted
    !
#ifdef COLD_STREAMS
! -------------------------------------------------   
    !   
    ! COLD FILAMENTARY STRUCTURE 
    !
    ! Respect to Lu et al, a fraction of the cold stream gas is disrupted and converted in hot gas
    call gas_void(dM_hot)
    f_sh   = baryon_halo_shock_heated_fraction(dm%M_acc)
    M_cold = gas_mass(bh%cold_gas) 
    M_hot  = gas_mass(bh%hot_gas) 
    !
    if ((M_cold .gt. 0.d0) .and. (M_hot .gt. 0.d0)) then
        ! 
        f_sh   = max(0.d0,M_cold - max(0.d0,min(1.d0,(1.d0-f_sh)))*(M_cold+M_hot))/M_cold
        ! computed the disrupted mass of cold stream
        dM_hot = f_sh*bh%cold_gas
        ! substract to the cold phase
        ! gas_sub(g1,g2[,therm_mode,called_by])
        call gas_sub(bh%cold_gas,dM_hot,therm_mode='isothermal')
        ! heat the gas, set temperature to bh hot gas temperature
        call gas_set_component(dM_hot,max(diffuse_gas_temp,gas_temp(bh%hot_gas)),component='Temp')
        ! add to the hot atmosphere
        call gas_add(bh%hot_gas,dM_hot)
    end if
! -------------------------------------------------
#endif
! COLD_STREAMS     
    !
    ! Compute the effective baryonic fraction
#ifdef PHOTOIONISATION
! -------------------------------------------------
    bf = baryonic_fraction_photoionisation(M_end,z)   ! use UV photo-ionisation processes 
! ------------
#else
! ------------
    bf = baryon_fraction 
! -------------------------------------------------
#endif
    !
    ! compute baryonic accretion rate
    ! acc_rate corresponds to the full accretion onto both the cold and hot gas phases
    acc_rate = bf*dm_acc_rate*InitAbund(1)   ! cosmological primordial ubundance     
    !
#ifdef REACCRETION
! -------------------------------------------------                               
    if (bf .lt. 9.9d-1*baryon_fraction) then
        !
        ! the residual fraction is added to the surrounding halo gas phase  
        surrounding_gas_inflow_rate = (baryon_fraction - bf)*dm_acc_rate*InitAbund(1) ! cosmological primordial ubundance
        call baryon_halo_add_surrounding_gas(bh,surrounding_gas_inflow_rate*dt)
    else
        !
        ! No more photoionization impact
        ! start to re-accrete surrounding gas
        call baryon_halo_compute_surrounding_gas_outflow_rate(bh,dm,dt,surrounding_gas_outflow_rate)
        acc_rate = acc_rate + surrounding_gas_outflow_rate
        ! the new re-accreted hot surrounding gas has to be remove from the surrounding gas reservoir
        call baryon_halo_sub_surrounding_gas(bh,surrounding_gas_outflow_rate*dt,called_by='baryon_halo_compute_accretion_rate')
    end if      
! -------------------------------------------------
#endif  
! REACCRETION                     
    !
    ! compute f_hot
#ifdef COLD_STREAMS
! -------------------------------------------------
    ! separate in cold and hot phase
    ! (1-f_hot) is the colimated filamentary accretion fraction
    !
    ! compute the fraction of the baryonic halo mass which must be in the isotropic atmosphere 
    if (gas_mass(acc_rate) .gt. 0.d0) then
       !
       ! compute the shock heated fraction
       f_sh  = baryon_halo_shock_heated_fraction(M_end)                                                       
       ! computed the optimal hot gas mass at the end of the time-step (the target mass)
       m_hot = f_sh*(baryon_halo_bh_mass(bh,component='cold+hot')+gas_mass(acc_rate)*dt)                 
       ! compute the hot gas fraction of the new accreted gas
       ! surrounding_gas_outflow_rate is only composed of hot gas --> the minimal hot gas fraction in the accretion has to
       ! be f_hot_min = surrounding_gas_outflow_rate/acc_rate
       f_hot = gas_mass(surrounding_gas_outflow_rate)/gas_mass(acc_rate)
       ! set the optimal hot gas fraction
       f_hot = min(1.d0,max(f_hot,max(0.d0,m_hot-baryon_halo_bh_mass(bh,component='hot'))/(gas_mass(acc_rate)*dt))) 
    end if
! ------------    
#else
! ------------
    ! we consider that the accretion is fully isotropic (no cold collimated gas accretion)
    f_hot = 1.d0 
! -------------------------------------------------
#endif
! COLD_STREAMS
    !    
    if (gas_mass(acc_rate) .ge. 0.d0) then   
       !
       ! set rates 
       bh%hot_inflow_rate  = acc_rate*f_hot        ! in code unit [10^11 Msun / Gyr]
       bh%cold_inflow_rate = acc_rate*(1.d0-f_hot) ! in code unit [10^11 Msun / Gyr]
       !
       ! set temperatures
       ! the hot accreted gas (shock-heated) is set to T_vir 
       call gas_set_component(bh%hot_inflow_rate,dm%T_vir,component='Temp')
       !
       ! the cold accreted gas is setled to ~10^4 K
       call gas_set_component(bh%cold_inflow_rate,diffuse_gas_temp,component='Temp')
       !
       ! evolution scheme
       if (bh%life_time .eq. 0.d0) then
          !
          ! it is the first evolution step of this bh component
          ! set bh%age_form
          call baryon_halo_set_age_form(bh,dm%age_form+dm%life_time)
          ! print starting step of the baryon-halo component
          ! baryon_halo_print(unit,form,bh)
          if (FOLLOW_UP .and. PR_FOLLOW_UP) call baryon_halo_print(follow_up_unit(current_index),'fits',bh)
      end if
    end if
    !
    return

#ifdef PHOTOIONISATION
! -------------------------------------------------
  contains
   
    ! *************************************************
 
    function filtering_mass(z)

      implicit none
      
      real(kind=8),intent(in)      :: z                ! redshift
      real(kind=8)                 :: filtering_mass   ! caracteristic mass (in code unit 10^11 Msun)

#ifdef GNEDIN_2000
! ------------------------------------------------- 
      ! Filtering mass Gnedin+00 and Kravtsov+04 (B2)

      real(kind=8)                 :: a,a0,ar          ! expension factor, expension factor at overlap and expension factor at reionization
      real(kind=8),parameter       :: aa = 6           ! alpha parameter in Kravstov+04 (B2)
      real(kind=8)                 :: f
      a  = 1.d0/(1.d0+z)                               ! expension factor
      a0 = 1.d0/(1.d0+z_overlap)                       ! The epoch where the first HII regions form
      ar = 1.d0/(1.d0+z_reion)                         ! The epoch of the complete reionization

      if (a .le. a0) then
         !
         f = (3.d0*a/((2.d0+aa)*(5.d0+2.d0*aa)))*(a/a0)**(aa)
      else
         !
         if (a .ge. ar) then
            !
            f = a0**2.*(1.d0/(2.d0+aa)-(2.d0*(a/a0)**(-1./2.))/(5.d0+2.d0*aa))
            f = f + (ar**2./1.d1)*(5.d0-4.d0*(a/ar)**(-1./2.)) - (a0**2./1.d1)*(5.d0-4.d0*(a/a0)**(-1./2.))
            f = f + a*ar/3.d0 - ar**2./3.d0*(3.d0-2.d0*(a/ar)**(-1./2.))
            f = 3.d0/a*f
         else
            !
            f = a0**2.*(1.d0/(2.d0+aa)-(2.d0*(a/a0)**(-1./2.))/(5.d0+2.d0*aa))  
            f = f + a**2./1.d1
            f = f - (a0**2./1.d1)*(5.d0-4.d0*(a/a0)**(-1./2.))
            f = 3.d0/a*f
         end if
      end if
      
      filtering_mass = 2.5d11*h_0**(-1.)*Omega_m**(-1./2.)*(mu/mp)**(-3./2.)*M_Sun_in_mass_code_unit*f**(3./2.)
! -------------------------------------------------
#endif
! GNEDIN_2000

#ifdef OKAMOTO_2008
! -------------------------------------------------
      ! Corrected Filtering mass Okamoto+08
      filtering_mass = 1.d10*h_0**(-1.)*exp(-7.d-1*z)*M_Sun_in_mass_code_unit
! -------------------------------------------------
#endif

      return
    end function filtering_mass

    ! *************************************************

    function baryonic_fraction_photoionisation(M_halo,z)

      ! Reionization feedback may cause the halo baryon fraction may be lower than the cosmic baryon fraction
      ! Gnedin+00 Eq. (7), Kravtsov+04 (B3) or Okamoto+08 Eq. (1)

      implicit none

      real(kind=8)                     :: baryonic_fraction_photoionisation  ! baryon fraction
      real(kind=8),intent(in)          :: M_halo                             ! for M_halo (in code unit)
      real(kind=8),intent(in)          :: z                                  ! at redshift z
      real(kind=8)                     :: M_fz,fb                            ! The critical mass below which reionization feedback becomes important
     
      if (z .lt. z_reion) then
         !
         ! apply photo-ionization prescription
         M_fz = filtering_mass(z)
#ifdef GNEDIN_2000
! -------------------------------------------------
         fb   = baryon_fraction / (1.d0 + 2.6d-1*(M_fz/M_halo))**3. 
! -------------------------------------------------
#endif
#ifdef OKAMOTO_2008
! -------------------------------------------------
         fb   = baryon_fraction / (1.d0 + 5.9d-1*(M_fz/M_halo)**2.)**(3./2.) 
! -------------------------------------------------
#endif
      else
         !
         ! apply universal baryonic fraction
         fb = baryon_fraction
      end if

      baryonic_fraction_photoionisation = fb

      return
    end function baryonic_fraction_photoionisation
   
    ! *************************************************
    ! end of contain : baryon_halo_compute_accretion_rate 
#endif
! -------------------------------------------------
! PHOTOIONISATION

  end subroutine baryon_halo_compute_accretion_rate 
  
  !*****************************************************************************************************************   

  subroutine baryon_halo_compute_surrounding_gas_outflow_rate(bh,dm,dt_min,surrounding_gas_outflow_rate)
    
    ! COMPUTE THE OUTFLOW RATE ASSOCIATED TO THE SURROUNDING HOT GAS RESERVOIR

    implicit none
    
    real(kind=8),intent(in)                 :: dt_min                       ! minimal time during wich the outflow_rate have to be maintained 
                                                                            ! At the halo scale, the evolution scheme imposes that tranfer rates are constant 
                                                                            ! during the optimal time step. The surrounding_gas_outflow_rate must be constant 
                                                                            ! during dt_min         
    
    type(gas_type),intent(out)              :: surrounding_gas_outflow_rate 
    type(baryon_halo_type),intent(inout)    :: bh                           ! a baryon halo component
    type(dm_type),intent(in)                :: dm                           ! the dark-matter component
    
    call gas_void(surrounding_gas_outflow_rate)                ! init
    ! 
    if (gas_mass(bh%surrounding_gas) .le. num_accuracy) return ! no mass
    
    surrounding_gas_outflow_rate = (1.d0/max(2.d0*dm%t_dyn,dt_min))*bh%surrounding_gas
    
    return
  end subroutine baryon_halo_compute_surrounding_gas_outflow_rate
  
  !*****************************************************************************************************************
  
  subroutine baryon_halo_compute_galaxy_fresh_gas_accretion_rate(bh,galaxy_fresh_gas_acc_rate)
  
    ! COMPUTE AND RETURN THE EFFECTIVE FRESH GAS ACCRETION RATE ONTO THE GALAXY
    
    implicit none
    
    type(gas_type),intent(out)           :: galaxy_fresh_gas_acc_rate
    type(baryon_halo_type),intent(inout) :: bh                        ! the baryon halo component
    
    call gas_void(galaxy_fresh_gas_acc_rate) ! init
    
    ! The inflow of the galaxy are the outflow of the baryon halo
    galaxy_fresh_gas_acc_rate = bh%cold_outflow_rate + bh%cooling_rate
    
    return
  end subroutine baryon_halo_compute_galaxy_fresh_gas_accretion_rate
  
  !***************************************************************************************************************** 

  subroutine baryon_halo_set_escape_fraction(bh,f_esc)

    ! SET ESCAPE FRACTION

    implicit none

    real(kind=8),intent(in)               :: f_esc  ! the escape fraction (of ejecta)

    type(baryon_halo_type),intent(inout)  :: bh     ! the baryon halo component
    
    bh%escape_fraction = f_esc
       
    return
  end subroutine baryon_halo_set_escape_fraction 
  
  !*****************************************************************************************************************
  
  subroutine baryon_halo_evolve_cold_gas_I(bh,dm,dt_optim)
  
    ! COMPUTE FIRST PART OF BARYON HALO COLD GAS EVOLUTION SCHEME (PREDICTOR)
    ! COMPUTE INPUT, OUTPUT RATES ASSOCIATED TO THE COLD GAS 
    ! These evolution rates are saved into the 'bh' structure
    ! the evolution (in mass) generated by this new transfer rates will be computed and will be saved in the subroutine 'baryon_halo_evolve_cold_gas_II'
    ! compute optimal dt for cold gas evolution
    
    implicit none
    
    real(kind=8),intent(out)              :: dt_optim
    
    type(baryon_halo_type),intent(inout)  :: bh               ! the baryon halo component
    type(dm_type),intent(in)              :: dm               ! the dark matter component

#ifdef PRINTALL 
! -------------------------------------------------
    call IO_print_message('baryon_halo_evolve_cold_gas_I',only_rank=rank,component='galaxy')
! -------------------------------------------------
#endif
      
    ! init
    call gas_void(bh%cold_outflow_rate)               ! init
    !
    ! Compute the outflow rate of the cold reservoir
    call baryon_halo_compute_cold_outflow_rate(bh,dm) 
    ! 
    ! Compute optimal time-step
    dt_optim = gas_dt_optim(bh%cold_gas,bh%cold_inflow_rate,bh%cold_outflow_rate,called_by='baryon_halo_evolve_cold_gas_I')

#ifdef PRINTALL 
! -------------------------------------------------
    call IO_print_message('baryon_halo_evolve_cold_gas_I ... done',only_rank=rank,component='galaxy')
! -------------------------------------------------
#endif
    
    return 
  end subroutine baryon_halo_evolve_cold_gas_I
  
  !*****************************************************************************************************************   

  subroutine baryon_halo_compute_cold_outflow_rate(bh,dm)
    
    ! COMPUTE THE FREE FALL RATE OF THE COLD HALO GAS PHASE

    implicit none
    
    real(kind=8)                         :: M_cold     ! cold gas mass
    real(kind=8)                         :: rate
    real(kind=8),parameter               :: f_str_cold_streams = 1.d0/3.d0 ! Fraction (in mass) of pre-fragmented cold (warm 10^4K) gas in the streams
    
    type(baryon_halo_type),intent(inout) :: bh         ! a baryon halo component
    type(dm_type),intent(in)             :: dm         ! the dark-matter component
    
    call gas_void(bh%cold_outflow_rate)                ! init
    !
    if (cold_stream_efficiency .le. 0.d0) return       ! no cold-streams
    !
    if (bh%quenched) return                            ! the baryonic accretion is quenched
    !
    M_cold = baryon_halo_bh_mass(bh,component='cold')
    ! 
    ! compute rate
    ! definition
    rate =  cold_stream_efficiency*5.d-1*M_cold/dm%t_dyn
    ! gas signature
    bh%cold_outflow_rate = rate*gas_signature(bh%cold_gas,apply_as='rate_builder',called_by='cold_stream_rate')
    !
    ! We assume that the condensed gas has an equilibrium temperature closed to 8000K
    call gas_set_component(bh%cold_outflow_rate,diffuse_gas_temp,component='Temp')
    !
    ! set the structured fraction  
    call gas_set_component(bh%cold_outflow_rate,f_str_cold_streams,component='f_str')
    
    return
  end subroutine baryon_halo_compute_cold_outflow_rate
  
  !*****************************************************************************************************************
  
  subroutine baryon_halo_evolve_cold_gas_II(bh,dt)
  
    ! EVOLVE THE COLD GAS PHASE (CORRECTOR PART)
    
    implicit none
    
    real(kind=8),intent(in)              :: dt    ! time-step duration
    real(kind=8)                         :: U,dU  ! mass check
        
    type(baryon_halo_type),intent(inout) :: bh    ! the baryon halo component

#ifdef PRINTALL 
! -------------------------------------------------
    call IO_print_message('baryon_halo_evolve_cold_gas_II',only_rank=rank,component='galaxy')
! -------------------------------------------------
#endif
        
    U = gas_mass(bh%cold_gas)
    !
    ! ADD MASS
    if (gas_mass(bh%cold_inflow_rate) .gt. 0.d0) then
      !
      ! accretion (filamentary cold stream formation)
      call gas_add(bh%cold_gas,bh%cold_inflow_rate*dt)    
    end if
    !
    ! SUBSTRACT MASS
    if (gas_mass(bh%cold_outflow_rate) .gt. 0.d0) then
      !
      ! freefall onto the galaxy 
      call gas_sub(bh%cold_gas,bh%cold_outflow_rate*dt,therm_mode='iso') 
    end if
    !
    dU = gas_mass(bh%cold_gas) - U
    if (U .gt. M_gas_min) then
       !
       if ((abs(dU)/U) .gt. physical_precision) then
          !
          call IO_print_error_message('No quasi-static state of the cold phase',only_rank=rank,called_by='baryon_halo_evolve_cold_gas')
          call IO_print_message('with',only_rank=rank,component='halo', &
                   param_name = (/'dU/U [%]                 ','U                        '/), &
                   real_param_val  = (/1.d2*abs(dU)/U,U/))
          stop  ! stop the program
       end if
    end if
    !
    if (is_NaN(gas_mass(bh%cold_gas))) then
       ! 
        call IO_print_error_message('cold_gas is NaN ',only_rank=rank,called_by='baryon_halo_evolve_cold_gas')
        stop ! stop the program
    end if

#ifdef PRINTALL 
! -------------------------------------------------
    call IO_print_message('baryon_halo_evolve_cold_gas_II ... done',only_rank=rank,component='galaxy')
! -------------------------------------------------
#endif
    
    return
  end subroutine baryon_halo_evolve_cold_gas_II

  !*****************************************************************************************************************

  subroutine baryon_halo_evolve_hot_gas_I(bh,dm,galaxy_ejecta_rate,f_in, &
                    dt_optim,galaxy_ejecta_rate_Wd,galaxy_ejecta_rate_Wu,f_in_Wd,f_in_Wu)
  
    implicit none
    
    real(kind=8),intent(inout)            :: f_in                     ! fraction of galaxy ejecta that is catched by the hot atmosphere [prediction]
    real(kind=8),intent(out)              :: f_in_Wu                  ! upper value of f_in
    real(kind=8),intent(out)              :: f_in_Wd                  ! lower value of f_in
    real(kind=8)                          :: f_in_crit                ! critical value of f_in to pass from input to output dominated evolution
    real(kind=8),intent(out)              :: dt_optim                 ! optimal time-step
    real(kind=8)                          :: dt_cool
    
    type(gas_type),intent(in)             :: galaxy_ejecta_rate       ! instantaneous galaxy ejecta rate
    type(gas_type),intent(out)            :: galaxy_ejecta_rate_Wd    ! galaxy ejecta barrier (down)
    type(gas_type),intent(out)            :: galaxy_ejecta_rate_Wu    ! galaxy ejecta barrier (up)
    type(gas_type)                        :: in_rate,out_rate         ! global evolution rates
    type(baryon_halo_type),intent(inout)  :: bh                       ! the baryon halo component
    type(dm_type),intent(in)              :: dm                       ! the dark matter component
  
#ifdef PRINTALL 
! -------------------------------------------------
    call IO_print_message('baryon_halo_evolve_hot_gas_I',only_rank=rank,component='galaxy')
! -------------------------------------------------
#endif  
    !
    ! init
    call gas_void(bh%cooling_rate)  
    !
    call baryon_halo_compute_cooling_rate(bh,dm,dt_cool)
    !
    ! save f_esc = 1. - f_in
    call baryon_halo_set_escape_fraction(bh,(1.d0-f_in))  
    !
    ! set f_in barriers
    f_in_Wu = min(1.d0,(1.d0 + physical_precision)*f_in)
    f_in_Wd = (1.d0 - physical_precision)*f_in
    !
    ! update f_in
    if (gas_mass(galaxy_ejecta_rate) .gt. 0.d0) then
        !
        f_in_crit = (gas_mass(bh%cooling_rate)-gas_mass(bh%hot_inflow_rate))/gas_mass(galaxy_ejecta_rate)
        if (f_in .gt. f_in_crit) then
            !
            ! in_rate dominated
            ! to be safe, we assume an increase of f_in
            f_in = f_in_Wu
        else
            !
            ! out_rate dominated
            ! to be safe, we assume a deacrease of f_in
            f_in = f_in_Wd
        end if
    end if
    !
    ! compute in and out rates
    in_rate  = f_in*galaxy_ejecta_rate + bh%hot_inflow_rate
    out_rate = bh%cooling_rate
    ! 
    ! compute optimal evolution time-step 
    ! BARYON-HALO FRAME; takes into account escape fraction of ejecta
    dt_optim = gas_dt_optim(bh%hot_gas,in_rate,out_rate, &
            evolving_rate=f_in*galaxy_ejecta_rate,warning_up=galaxy_ejecta_rate_Wu,warning_down=galaxy_ejecta_rate_Wd, &
            called_by='baryon_halo_evolve_hot_gas_I')
    ! In the hot atmosphere, only a fraction of the mass is already condensed: dM_condensate.
    ! This condensed mass is transfered to the galaxy according to dM_cool/dt = dM_condensate / dt_cool, 
    ! dt_cool = r_cool / v_dm(r=r_cool), r_cool or r_TI if r_TI smaller
    ! In the next time step, the mass transfered into the galaxy cannot be larger than dM_condensate
    ! therefore : dt_optim < dt_cool
    if ((dt_cool > 0.d0) .and. (dt_optim > dt_cool)) then
        !
        dt_optim = 9.9d-1*dt_cool
    end if
    !
    ! GALAXY-FRAME
    ! update galaxy ejecta barriers: WARNING_up,warning_down=WARNING_down
    if (gas_mass(galaxy_ejecta_rate) .gt. 0.d0) then
        !
        if (f_in .gt. 0.d0) then
            !
            galaxy_ejecta_rate_Wu = (1./f_in)*galaxy_ejecta_rate_Wu
            galaxy_ejecta_rate_Wd = (1./f_in)*galaxy_ejecta_rate_Wd
        else
            !
            galaxy_ejecta_rate_Wu = (1.d0 + physical_precision)*galaxy_ejecta_rate
            galaxy_ejecta_rate_Wd = (1.d0 - physical_precision)*galaxy_ejecta_rate
        end if
    end if

#ifdef PRINTALL 
! -------------------------------------------------
    call IO_print_message('baryon_halo_evolve_hot_gas_I ... done',only_rank=rank,component='galaxy')
! -------------------------------------------------
#endif 
    
    return
  end subroutine baryon_halo_evolve_hot_gas_I
  
  !*****************************************************************************************************************
  
  subroutine baryon_halo_compute_cooling_rate(bh,dm,t_ff)
    
    ! COMPUTE THE COOLING RATE OF THE HOT HALO GAS PHASE

    implicit none
    
    integer(kind=4)                      :: loop                        ! loop index
    integer(kind=4)                      :: final_order                 ! final order of the integrator procedure

    real(kind=8)                         :: Mhot                        ! Mass of hot gas
    real(kind=8),intent(out)             :: t_ff                        ! timescale of galaxy feeding
    real(kind=8)                         :: Z_g,T_g,T0,rho_g            ! baryon halo gas properties (temperatures, metalicity and local density)
    real(kind=8)                         :: r_TI
    real(kind=8)                         :: r_cool, lr_max, lr_min, lr  ! local radius variables (mainly in log scale)
    real(kind=8)                         :: t, t_max, t_min             ! cooling time (at a given radius r)
    real(kind=8)                         :: mass                        ! local mass variable
    real(kind=8)                         :: error                       ! error onto integrated hot halo mass (given by the integrator)
    real(kind=8),parameter               :: x_cut_cool = 1.d-4          ! minimal cooling radius                                                                           ! at all redshif R90 (radius which enclose 90% of the galaxy mass) 
                                                                        ! is larger than x_cut_cool*r_core                                                                                                                                      
    type(baryon_halo_type),intent(inout) :: bh                          ! the baryon_halo component
    type(dm_type),intent(in)             :: dm                          ! the dark-matter component
      
#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('baryon_halo_compute_cooling_rate',only_rank=rank,component='bh')
! -------------------------------------------------
#endif
    
    call gas_void(bh%cooling_rate)                   ! init
    r_cool    = bh%r_cool                            ! init
    bh%r_cool = -1.d0                                ! init
    bh%r_TI   = -1.d0                                ! init
    t_ff      = -1.d0                                ! init
    if (bh%quenched) return                          ! the hot atmosphere et quenched, no more cooling 
    if (cooling_efficiency .le. 0.d0) return         ! no cooling 
    if (bh%t_cool .le. 0.d0) return                  ! no cooling
    !
    Mhot = baryon_halo_bh_mass(bh,component='hot')   ! Mass of hot gas
    Z_g = gas_metallicity(bh%hot_gas)                ! compute metalicity (mZ/m)
    T_g = gas_temp(bh%hot_gas)                       ! Temperature of the hot gas
    !
    if (Z_g .lt. 0.d0) then
       !
       call IO_print_error_message('Z_g < 0',only_rank=rank,called_by='baryon_halo_compute_cooling_rate')
       call IO_print_message('with',only_rank=rank,component='bh', &
          param_name=(/'Z_g                      ','hot gas mass             '/), &
          real_param_val=(/Z_g,gas_mass(bh%hot_gas)/))
       stop ! stop the program
    end if
    !
    ! checks temperature
    if ((Mhot .gt. M_gas_crit) .and. (T_g .lt. 9.9d-1*diffuse_gas_temp)) then
        !
        call IO_print_warning_message('T_g < T_cool_gas',only_rank=rank,called_by='baryon_halo_compute_cooling_rate')
#ifdef PRINT_WARNING
! -------------------------------------------------         
            call IO_print_message('with',only_rank=rank,component='bh', &
                    param_name = (/'M_hot                    ','T_g                      '/), &
                    real_param_val  = (/Mhot,T_g/))
#endif
! -------------------------------------------------     
! PRINT_WARNING 
        call gas_set_component(bh%hot_gas,diffuse_gas_temp,component='Temp')
    end if     
    !
    T0    = bh%T0
    !
#ifdef POLYTROPIC
! -------------------------------------------------
    if ((T_g .lt. 0.d0) .or. (T0 .lt. 0.d0)) then
       !
       ! In the case of a polytropic gas profile the hot gas can be unstable (T0 < T0_crit_polytropic)
       if (T0 .eq. -2.d0) then ! -2.d0 is used as a caracteristical value for unstable hot gas phase
          !
          ! compute free fall rate
          bh%cooling_rate = cooling_efficiency/dm%t_dyn*bh%hot_gas  
          bh%r_cool = dm%R_vir ! set cooling radius, in this case, the hot atmosphere is considered as unstable, 
          ! all the gas can cool in a dynamical time and therefore the cooling radius is set to the Virial radius of the dark matter structure. 
          return
       else
          !
          bh%r_cool = 0.d0
          call gas_void(bh%cooling_rate) ! init
          return
       end if
    end if
! ------------    
#else
! ------------
    ! perfect gas case
    if (T_g .le. 0.d0) then
       !
       bh%r_cool = 0.d0  ! set cooling radius
       bh%r_TI   = 0.d0  ! set cooling radius
       call gas_void(bh%cooling_rate) ! init
       return
    end if
! -------------------------------------------------
#endif
! POLYTROPIC
    !
    ! check mass
    ! we integrate the hot atmosphere profile until R_halo and we compare the result with the total hot gas mass
    mass = baryon_halo_integrate_hot_halo_density_profile(dm%R_halo,bh,dm,error, &
         called_by='baryon_halo_compute_cooling_rate')
    error  = abs(gas_mass(bh%hot_gas)-mass)/gas_mass(bh%hot_gas)
    if (error .gt. num_precision) then
       !
       call IO_print_error_message('No corresponding mass',only_rank=rank, &
            called_by='baryon_halo_compute_cooling_rate')
       call IO_print_message('with',only_rank=rank,component='bh', &
          param_name=(/'err (%)                  ','hot gas mass             ','integreted mass          ', &
                       'T_g                      ','rho0                     '/), &
          real_param_val=(/1.d2*error,gas_mass(bh%hot_gas),mass,T_g,bh%rho0/))
       stop
    end if
    !
    ! *************************************************
    ! compute cooling radius 
    ! *************************************************
    !    
    lr_max = log10(dm%R_halo)
    ! deduce properties at r
    rho_g = baryon_halo_hot_halo_density_profile(10.**(lr_max),bh,dm)      ! gas density at radius r_max     
    T_g   = baryon_halo_hot_halo_temperature_profile(10.**(lr_max),bh,dm)  ! temperature at radius r_max
    ! check with the maximal extension --> t is the effective time to cool all the hot gas 
    t_max = t_cool(T_g,rho_g,Z_g)                                          
    if (t_max .gt. bh%t_cool) then
       !
       ! only a fraction of the hot atmosphere can cool during bh%t_cool
       lr_min = log10(x_cut_cool*dm%r_core)
       ! deduce properties at r
       rho_g = baryon_halo_hot_halo_density_profile(10.**(lr_min),bh,dm)      ! gas density at radius r_min   
       T_g   = baryon_halo_hot_halo_temperature_profile(10.**(lr_min),bh,dm)  ! temperature at radius r_min
       ! check with the minimal extension (x_cut_cool*dm%r_core)
       t_min = t_cool(T_g,rho_g,Z_g)  
       ! 
       if (t_min .lt. bh%t_cool) then
          !
          ! *************************************************
          ! effective cooling rate < free-fall rate
          ! *************************************************
          !
          ! @ this point, the cooling radius r_cool in then between r_min = x_cut_cool*dm%r_core and r_max = R_halo
          ! We use a dichotomy algorithm to find the effective cooling radius
          ! The computation is done by using a log scale
          loop = 0  ! loop is an loop index used to limit the number of research loop
          lr = 5.d-1*(lr_max + lr_min)
          ! deduce properties at r
          rho_g = baryon_halo_hot_halo_density_profile(10.**(lr),bh,dm)      ! gas density at radius r     
          T_g   = baryon_halo_hot_halo_temperature_profile(10.**(lr),bh,dm)  ! temperature at radius r
          ! compute cooling timescale
          ! t_cool(T_g,rho_g,Z_g)
          t = t_cool(T_g,rho_g,Z_g) ! t is an increasing function of the radius. Indeed : 
                                    !            - The density of the hot gas deacrese with the radius
                                    !            - t is an deacresing function of the density
          do while ((abs(t-bh%t_cool))/bh%t_cool .gt. num_precision)
             if (loop .gt. 1000) then
                !
                call IO_print_error_message('Too much cycle used in t_cool computation',only_rank=rank, &
                     called_by='baryon_halo_compute_cooling_rate')
                call IO_print_message('Computation used',only_rank=rank,component='bh', &
                     param_name=(/'(dm)r_core               ','R_halo                   ','(dm)rho_core             ', &
                                  'rho0                     ','T_g                      ','T0                       ', &
                                  't_min                    ','t_max                    '/), &
                     real_param_val=(/dm%r_core,dm%R_halo,dm%rho_core,bh%rho0,T_g,T0,t_min,t_max/))
                stop ! stop the program
             end if
             ! t is an increasing function of the radius r
             if (t .gt. bh%t_cool) then  
                !       
                lr_max = lr
             else
                !
                lr_min = lr 
             end if
             loop = loop +1
             lr = 5.d-1*(lr_max + lr_min)
             ! deduce properties at r
             rho_g = baryon_halo_hot_halo_density_profile(10.**(lr),bh,dm)      ! gas density at radius r     
             T_g   = baryon_halo_hot_halo_temperature_profile(10.**(lr),bh,dm)  ! temperature at radius r
             ! compute cooling timescale
             ! t_cool(T_g,rho_g,Z_g)
             t = t_cool(T_g,rho_g,Z_g)
          end do
          !
          ! set r_cool
          r_cool = (1.d1)**(lr)  
          !
          ! test
          if (r_cool .gt. dm%R_halo) then
             !
             call IO_print_warning_message('r_cool > R_halo',only_rank=rank, &
                  called_by='baryon_halo_compute_cooling_rate') 
             r_cool = min(dm%R_halo,(1.d1)**(lr))  ! re-scaling
          end if
          ! 
          if (r_cool .lt. dm%R_halo) then
             !
             ! we must compute a real cooling rate by integreted the hot gas density profil until r_cool
             mass = baryon_halo_integrate_hot_halo_density_profile(r_cool,bh,dm,error, &
                  called_by = 'baryon_halo_compute_cooling_rate',final_order=final_order)   
             if (mass-error .gt. gas_mass(bh%hot_gas)) then
                !
                ! it is possible for numerical integration and density shape profile reasons that the cool mass is greater than the total mass
                ! we check the error and we limit the cooling mass to bh%hot_mass
                call IO_print_warning_message('Integrated hot gas mass (<r_cool) > bh%hot_gas',only_rank=rank, &
                     called_by='baryon_halo_compute_cooling_rate')
                call IO_print_message('used',only_rank=rank,component='disc', &
                param_name = (/'mass                     ','M_hot                    ', &
                               'r_cool                   ','R_halo                   ', &
                               'r_cool/R_halo            '/), &
                real_param_val  = (/mass,Mhot,r_cool,dm%R_halo,r_cool/dm%R_vir/)) 
                mass = gas_mass(bh%hot_gas)
             end if
          else
             !
             ! r_cool = R_halo
             r_cool = dm%R_halo
             mass = gas_mass(bh%hot_gas)
          end if
       else
          !
          ! no mass can cool during bh%t_cool (too short ...)
          r_cool = 0.d0   
          mass   = 0.d0
          !
       end if
    else
       !
       ! free-fall 
       r_cool = dm%R_halo
       mass = gas_mass(bh%hot_gas)
    end if
    !
    ! Compute effective cooling rate
    if ((mass .gt. M_gas_crit) .and. (r_cool .gt. 0.d0)) then
        !
        ! compute dynamical time at the cooling radius
        t_ff = r_cool/dm_velocity(r_cool,dm)
#ifdef THERMAL_INSTABILITIES
! -------------------------------------------------   
        ! Take into account thermal instability
        ! Compute thermal instability time scale (Sharma+2012b, Cornuault+16)
        ! We assume that the condensed mass is unstable, and can be mixed with the hot gas due to thermal instabilities
        ! The effective cooling rate if therefore reduced
        ! compute thermal instabilitie radius
        r_TI = TI_efficiency*sqrt(adiab_ind*r_perfect_gas*T_g/(mu/mp))*bh%t_TI
        !
        if (r_TI .gt. 0.d0) then
            !
            if (r_TI .ge. r_cool) then
                !
                ! no cooling,
                ! The centre of the hot halo phase is dominated by thermal instabilities and gas mixing 
                mass = 0.d0   
                r_TI = 0.d0
                if (dm%M_halo .gt. 1.d1) bh%quenched = .true.
            else
                !
                ! Only a band around the cooling radius is impacted by thermal instabilities
                r_TI = max(0.d0,r_cool - r_TI)
                ! this unstable band cannot cool
                t_ff = r_TI/dm_velocity(r_TI,dm)
                ! we compute the residual cooling mass
                mass = baryon_halo_integrate_hot_halo_density_profile(r_TI,bh,dm,error, &
                        called_by = 'baryon_halo_compute_cooling_rate',final_order=final_order)
            end if  
        else
            !
            r_TI = r_cool
        end if
! -------------------------------------------------     
#endif 
! THERMAL_INSTABILITIES
        bh%cooling_rate = cooling_efficiency*mass/t_ff*gas_signature(bh%hot_gas,apply_as='rate_builder',called_by='cooling_rate')
        !
        ! set warm gas temperature
        ! We assume that the condensed gas has an equilibrium temperature closed to 8000K
        call gas_set_component(bh%cooling_rate,diffuse_gas_temp,component='Temp')
        ! We assume a fix pre-structured fraction
        call gas_set_component(bh%cooling_rate,1.d0/3.d0,component='f_str')
    else
        !
        ! too few condensed mass
        call gas_void(bh%cooling_rate) ! init
        ! update cooling radius
        bh%r_cool = 0.d0
        ! update thermal instability radius
        bh%r_TI = 0.d0
        return
    end if      
    !
    ! update cooling radius
    bh%r_cool = r_cool
    ! update thermal instability radius
    bh%r_TI = r_TI
    !
    if (is_NaN(gas_mass(bh%cooling_rate))) then
        !
        call IO_print_error_message('cooling_rate is NaN ',only_rank=rank,called_by='baryon_halo_compute_cooling_rate')
        write(*,*) mass, r_cool, r_TI, T_g, Z_g
        call gas_void(bh%cooling_rate)
    end if

#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('baryon_halo_compute_cooling_rate ... done',only_rank=rank,component='bh')
! -------------------------------------------------
#endif

    return

  end subroutine baryon_halo_compute_cooling_rate
    
  !*****************************************************************************************************************
  
  subroutine baryon_halo_evolve_hot_gas_II(bh,dm,dt,galaxy_ejecta_rate,Vwind)

    ! EVOLVE THE HOT GAS PHASE (CORRECTOR PART)
    
    implicit none

    real(kind=8),intent(in)               :: dt                   ! time-step duration
    real(kind=8),intent(in)               :: Vwind                ! (time)-mean velocity of galactic wind in [km/s]
    real(kind=8)                          :: Vesc                 ! Escape velocity of the dark matter component
    real(kind=8)                          :: fesc_pred            ! predicted escape fraction (of ejecta)
    real(kind=8)                          :: fesc                 ! escape fraction (of ejecta)
    real(kind=8)                          :: U,dU                 ! mass check
    real(kind=8)                          :: T_ini                ! initial temperature
    real(kind=8)                          :: T_after_infall, T_after_ej, T_after_cooling
    real(kind=8)                          :: M_hot_before_evap
    
    type(gas_type),intent(inout)          :: galaxy_ejecta_rate   ! ejected gas in code unit [10^11.Msun]
    type(gas_type)                        :: m_esc                ! escape mass fraction (wind mass which leaves the hot halo phase, V > Vesc)
    type(gas_type)                        :: m_evap               ! escape mass fraction (hot gas mass which leaves the hot halo phase, V > Vesc)
    type(baryon_halo_type),intent(inout)  :: bh                   ! the baryon halo component
    type(dm_type),intent(in)              :: dm                   ! the dark matter component

#ifdef PRINTALL 
! -------------------------------------------------
  call IO_print_message('baryon_halo_evolve_hot_gas',only_rank=rank,component='bh')
! -------------------------------------------------
#endif

    ! 
    ! Init 
    U = gas_mass(bh%hot_gas)
    T_ini = gas_temp(bh%hot_gas)
    T_after_infall  = -1.d0
    T_after_ej      = -1.d0
    T_after_cooling = -1.d0
    ! 
    call gas_void(m_esc)
    call gas_void(m_evap)
    !
    ! EVOLVE HOT HALO PHASE
    !
    ! INFALL
    if (gas_mass(bh%hot_inflow_rate) .gt. 0.d0) then 
        !
        ! background accretion (onto hot atmosphere)
        ! gas_add(g1,g2[,Vwind,Vesc,m_esc,called_by])
        call gas_add(bh%hot_gas,bh%hot_inflow_rate*dt,called_by='baryon_halo_evolve_hot_gas / accretion')   
        if (gas_mass(bh%hot_gas) .gt. 0.d0) then
            !
            T_after_infall = gas_temp(bh%hot_gas)
        end if
    end if  
    !
    ! EJECTA
    if (gas_mass(galaxy_ejecta_rate) .gt. 0.d0) then
        !
        ! outflow from the galaxy 
        Vesc = dm_escape_velocity(dm%r_core,dm)
        ! gas_add(g1,g2[,Vwind,Vesc,m_esc,called_by])
        call gas_add(bh%hot_gas,galaxy_ejecta_rate*dt,Vwind=Vwind,Vesc=Vesc,m_esc=m_esc,called_by='baryon_halo_evolve_hot_gas / ejecta') 
        fesc_pred = baryon_halo_escape_fraction(bh)                                               
        ! set escape fraction
        fesc = gas_mass(m_esc)/gas_mass(galaxy_ejecta_rate*dt)
        call baryon_halo_set_escape_fraction(bh,fesc)
        if ((fesc .lt. 1.d0) .and. (gas_mass(bh%hot_gas) .gt. 0.d0)) then
            !
            T_after_ej = gas_temp(bh%hot_gas)
        end if  
    end if
    ! 
    ! COOLING
    if (gas_mass(bh%cooling_rate) .gt. 0.d0) then
        !
        ! substract condensed gas du to radiative cooling of the hot phase
        ! gas_sub(g1,g2[,therm_mode,called_by])
        call gas_sub(bh%hot_gas,bh%cooling_rate*dt,therm_mode='cooling',called_by='baryon_halo_evolve_hot_gas / cooling')          
        if (gas_mass(bh%hot_gas) .gt. 0.d0) then
            !
            T_after_cooling = gas_temp(bh%hot_gas)
        end if
    end if
    !
    ! UPDATE THE COOLING CLOCK
    call baryon_halo_update_cooling_clock(bh,dm,dt,(bh%hot_inflow_rate*dt+(galaxy_ejecta_rate*dt-m_esc)))
    !
    ! CHECK I 
    ! compute mass variation of the hot halo phase, crash code if the variation is greater than the physical asked precision
    if (U .gt. M_gas_min) then
        !
        dU = gas_mass(bh%hot_gas) - U
        if ((abs(dU)/U) .gt. physical_precision) then
            !
            call IO_print_error_message('No quasi-static state of the hot phase',only_rank=rank,called_by='baryon_halo_evolve_hot_gas')
            call IO_print_message('with',only_rank=rank,component='bh', &
                    param_name = (/'dt                       ','dU/U [%]                 ','U                        ','hot_inflow               ',&
                                   'galaxy_ejecta            ','cooling                  ','T_ej                     ','Vesc                     ',&
                                   'Vwind                    ','m_esc                    ','fin                      ','fin_pred                 '/), &
                    real_param_val  = (/dt,1.d2*dU/U,U,gas_mass(bh%hot_inflow_rate)*dt, &
                                        gas_mass(galaxy_ejecta_rate)*dt,gas_mass(bh%cooling_rate)*dt,gas_temp(galaxy_ejecta_rate), &
                                        dm_escape_velocity(dm%r_core,dm),Vwind,gas_mass(m_esc),(1.d0-fesc),(1.d0-fesc_pred)/))
            stop  ! stop the program
        end if
    end if
    !
    ! EVAPORATION
    if (gas_mass(bh%hot_gas) .gt. 0.d0) then
        !
        M_hot_before_evap = gas_mass(bh%hot_gas)
        Vesc = dm_escape_velocity(dm%r_core,dm)
        call gas_evap(bh%hot_gas,m_evap,Vesc)
    end if
    !
    ! UPDATE EJECTA_RATE        
    call gas_void(galaxy_ejecta_rate)             ! reset
    galaxy_ejecta_rate = (1.d0/dt)*(m_esc+m_evap) ! update
    !
    ! CHECKS II 
    ! Mass / Temp
    if ((gas_mass(bh%hot_gas) .gt. 0.d0) .and. (gas_temp(bh%hot_gas) .le. 0.d0)) then
        !
        call IO_print_error_message('Massive hot gas phase without setled temperature',only_rank=rank,called_by='baryon_halo_evolve_hot_gas')
        call IO_print_message('with',only_rank=rank,component='bh', &
                    param_name = (/'M_hot                    ','T_hot                    '/), &
                    real_param_val  = (/gas_mass(bh%hot_gas),gas_temp(bh%hot_gas)/))
        call IO_print_message('with',only_rank=rank,component='bh', &
                    param_name = (/'M_hot                    ','M_hot_before_evap        ','T_hot                    ', &
                                   'T_ini                    ','T_after_infall           ', &
                                   'T_after_ej               ','T_after_cooling          '/), &
                    real_param_val  = (/gas_mass(bh%hot_gas),M_hot_before_evap,gas_temp(bh%hot_gas),T_ini,T_after_infall,T_after_ej,T_after_cooling/))          
        stop  ! stop the program
    end if
    !
    if ((gas_mass(bh%hot_gas) .gt. 0.d0) .and. (gas_temp(bh%hot_gas) .lt. 9.9d-1*diffuse_gas_temp)) then
        !
        call IO_print_warning_message('T_hot < T_cool_gas',only_rank=rank,called_by='baryon_halo_evolve_hot_gas')
#ifdef PRINT_WARNING
! -------------------------------------------------         
        call IO_print_message('with',only_rank=rank,component='bh', &
                    param_name = (/'M_hot                    ','M_hot_before_evap        ','T_hot                    ', &
                                   'T_ini                    ','T_after_infall           ', &
                                   'T_after_ej               ','T_after_cooling          '/), &
                    real_param_val  = (/gas_mass(bh%hot_gas),M_hot_before_evap,gas_temp(bh%hot_gas),T_ini,T_after_infall,T_after_ej,T_after_cooling/))
#endif
! -------------------------------------------------     
! PRINT_WARNING 
        call gas_set_component(bh%hot_gas,diffuse_gas_temp,component='Temp')
    end if
    !
    if (is_NaN(gas_mass(bh%hot_gas))) then
        !
        call IO_print_error_message('hot_gas is NaN ',only_rank=rank,called_by='baryon_halo_evolve_hot_gas')
        stop ! stop the program
    end if
    !
    ! UPDATE HOT GAS DENSITY PROFILE PARAMETERS 
    ! baryon_halo_update_density_profile(bh,dm[,called_by])
    call baryon_halo_update_density_profile(bh,dm,called_by='baryon_halo_evolve_hot_gas')
    !
    ! UPDATE TIME  
    ! bh%life_time
    call baryon_halo_update_life_time(bh,dt)
    ! print baryon-halo properties 
    ! baryon_halo_print(unit,form,bh)
    if (FOLLOW_UP .and. PR_FOLLOW_UP) call baryon_halo_print(follow_up_unit(current_index),'fits',bh)
    ! 
     
#ifdef PRINTALL 
! -------------------------------------------------
  call IO_print_message('baryon_halo_evolve_hot_gas ... done ',only_rank=rank,component='bh')
! -------------------------------------------------
#endif

    return
  end subroutine baryon_halo_evolve_hot_gas_II
  
  !*****************************************************************************************************************

  subroutine baryon_halo_update_density_profile(bh,dm,called_by)
 
    ! UPDATE HOT ATMOSPHERE PROFILE PARAMETERS

    implicit none

    character(*),intent(in),optional     :: called_by      ! name of the function which has called this function
    character(MAXPATHSIZE)               :: message        ! a message to display
    
    real(kind=8)                         :: T_hot,T0       ! Temperature (central and mean)
    real(kind=8)                         :: M_hot          ! mass of hot gas
    real(kind=8)                         :: M_test         ! a test variable
    real(kind=8)                         :: rho0           ! core gas density
    real(kind=8)                         :: error          ! integration error
   
    type(baryon_halo_type),intent(inout) :: bh             ! the baryon halo component
  
    type(dm_type),intent(in)             :: dm             ! the dark matter component

#ifdef PRINTALL 
! -------------------------------------------------
  call IO_print_message('baryon_halo_update_density_profile',only_rank=rank,component='bh')
! -------------------------------------------------
#endif
    !
    if (present(called_by)) then
        write(message,'(a,a)') 'called by: ', trim(called_by)
    end if
    !
    T_hot = gas_temp(bh%hot_gas) 
    !
    if (T_hot .gt. 0.d0) then
       !
       ! the hot atmosphere exists 
       ! compute central temperature T0
#ifdef POLYTROPIC
! -------------------------------------------------
       ! if the central temperature is lower than the minimal thershold, return -2.
       T0 = T0_polytropic(bh_tmp,dm) 
! ------------
#else
! ------------
       ! in the ideal gas case the central temperature T0 is equal to the mean temperature T_hot
       T0 = T_hot
#endif   
! -------------------------------------------------       
! POLYTROPIC
       !
       ! UPDATE T0
       bh%T0 = T0    ! update the local copy of bh
       !
       if (T0 .gt. 0.d0) then
          !
          ! Compute rho0 [10^11Msun/kpc^3]
          M_hot = gas_mass(bh%hot_gas)
          rho0 = M_hot / (4.d0*pi*dm%r_core**3.*Ronbint(M,x_cut,dm%concentration,(/b_const(T0,dm),-1.d0/), &
                  called_by = 'baryon_halo_update_density_profile (rho0)'))
          ! UPDATE rho0        
          bh%rho0 = rho0   ! [10^11Msun/kpc^3]                    
          ! check mass by integrating the hot gas density profile
          M_test = baryon_halo_integrate_hot_halo_density_profile(dm%R_halo,bh,dm,error=error, &
               called_by='baryon_halo_update_density_profile')
          if ((abs(M_test-M_hot)/M_hot) .gt. num_precision) then
             !
             call IO_print_error_message('integrated mass (M) != M_hot',only_rank=rank, &
                  called_by='baryon_halo_update_density_profile')
             if (present(called_by)) call IO_print_message(trim(message),only_rank=rank,component='bh')     
             call IO_print_message('Computation used',only_rank=rank,component='bh', &
                     param_name=(/'M                        ','error                    ','M_hot                    ', &
                                  '(dm)r_core               ','R_halo                   ','(dm)rho_core             ', & 
                                  'rho0                     ','T_hot                    ','T0                       '/), &
                     real_param_val=(/M_test,error,baryon_halo_bh_mass(bh,component='hot'),dm%r_core,dm%R_halo,dm%rho_core,rho0,T_hot,T0/))
             stop ! stop the program
          end if
       else
          !
          ! Acording to the polytropic model
          ! the hot atmosphere is unstable 
          bh%rho0 = -1.d0
       end if
    else
       !
       ! the hot atmosphere doesn't exist
       bh%T0   = -1.d0 
       bh%rho0 = -1.d0 
    end if

#ifdef PRINTALL 
! -------------------------------------------------
  call IO_print_message('baryon_halo_update_density_profile ... done ',only_rank=rank,component='bh')
! -------------------------------------------------
#endif

    return
  end subroutine baryon_halo_update_density_profile
  
  !*****************************************************************************************************************
  
  subroutine baryon_halo_update_cooling_clock(bh,dm,dt,gas)
    
    ! UPDATE THE COOLING CLOCK AND THE THERMAL INSTABILITY CLOCK 
    ! compute average cooling time by adding a mass 'm' evolving since 'dt' in the hot atmosphere
    ! The gas is assumed to be continously added to the reservoir,
    ! We take into account the half of the time-step (dt/2.d0), 
    ! WARNING: When this routine is called the mass has already been added to the hot reservoir

    implicit none

    real(kind=8),intent(in)                :: dt       ! the time-step
    real(kind=8)                           :: m1,m2    ! mass of gas
    real(kind=8)                           :: Tg,Zg
    real(kind=8)                           :: old_t_cool 
    real(kind=8)                           :: new_t_TI 
    real(kind=8)                           :: dt_TI

    type(gas_type),intent(in)              :: gas      ! the added gas
    type(dm_type),intent(in)               :: dm       ! the dark matter component
    type(baryon_halo_type),intent(inout)   :: bh       ! A baryon halo component
    
    m1 = gas_mass(bh%hot_gas)                          ! the total hot gas mass
    
    if (m1 .le. M_gas_min) return                           ! no hot atmosphere
    
    m2 = gas_mass(gas)                                 ! the mass of accreted or ejecting (just come in)
    
    ! cooling clock
    ! save previous value 
    old_t_cool = bh%t_cool 
    ! update previous value
    bh%t_cool  = ((old_t_cool+dt)*max(0.d0,(m1-m2)) + (dt/2.d0)*m2)/(m1)  

    if (bh%r_TI .gt. 0.d0) then
        !
        ! thermal instability clock
        Tg = gas_temp(bh%hot_gas)
        Zg = gas_metallicity(bh%hot_gas)
        dt_TI = (1.d0/1.d1)*dlnLambda_dlnT(Tg,Zg)*dt ! 1./10. is an efficiency parameter dt_TI <= dt
        if (dt_TI .gt. 0.d0) then
            !
            ! unstable regime
            new_t_TI = ((bh%t_TI+dt_TI)*max(0.d0,(m1-m2)) + (dt_TI/2.d0)*m2)/(m1)
            if (dm%M_halo .gt. 1.0) then
                !
                if (new_t_TI .gt. bh%t_TI) bh%t_TI  = new_t_TI
            else 
                !
                bh%t_TI  = new_t_TI
            end if
        end if
    end if
        
    return
  end subroutine baryon_halo_update_cooling_clock

  !*****************************************************************************************************************
  
  subroutine baryon_halo_update_life_time(bh,dt)

    ! INCREASE THE TIME LIFE OF A BARYON HALO COMPONENT
    
    implicit none
    
    real(kind=8),intent(in)               :: dt   ! time step

    type(baryon_halo_type),intent(inout)  :: bh   ! the baryon_halo component
    
    if (baryon_halo_bh_mass(bh) .le. 0.d0) return
    
    bh%life_time = bh%life_time + dt              ! increase the life time of the baryon halo component
    
    return
  end subroutine baryon_halo_update_life_time

  !*****************************************************************************************************************

  subroutine baryon_halo_merge(bh1,bh2,dm1)
    
    ! MERGE TWO BARYON HALO COMPONENTS

    implicit none

    real(kind=8)                         :: f_sh                  ! optimal fraction of hot gas in the baryon halo
    real(kind=8)                         :: Mhot_1, Mhot_2        ! hot masses
    real(kind=8)                         :: M_tot, M_cold, M_hot  ! masses
    
    real(kind=8)                         :: Eint                  ! gravitational energy of the merger, injected into the hot gas
    real(kind=8)                         :: mu_gas                ! gas ratio
    real(kind=8)                         :: Vesc                  ! escape velocity of the dark matter host halo
    
    type(gas_type)                       :: dM_hot                ! mass of cold filamentary structure that has been heated during the merger
    type(gas_type)                       :: m_evap                ! evaporated hot mass
    type(dm_type),intent(in)             :: dm1                   ! the remanent dark matter halo
    type(baryon_halo_type),intent(inout) :: bh1                   ! a baryon_halo component (remnent)
    type(baryon_halo_type),intent(in)    :: bh2                   ! an other baryon_halo component
    
#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('baryon_halo_merge',only_rank=rank,component='bh')
! -------------------------------------------------
#endif
    
    !
    ! set age of the formation and life time 
    bh1%age_form  = min(bh1%age_form,bh2%age_form)
    bh1%life_time = max(bh1%life_time,bh2%life_time)
    if (FOLLOW_UP .and. PR_FOLLOW_UP) then
      !
      ! in this case, the bh2 is the followed structure
      ! we have to erase time evolution properties associated to bh1 and replace them by bh2 properties
      bh1%life_time = bh2%life_time  ! the life time of the bh structure
      bh1%age_form  = bh2%age_form   ! the age of the universe when the structure has been formed
    end if
    
    ! merge surounding gas component
    call gas_add(bh1%surrounding_gas,bh2%surrounding_gas)
          
    if (baryon_halo_bh_mass(bh1) .gt. 0.d0) then
      !
      if (baryon_halo_bh_mass(bh2) .gt. 0.d0) then
        ! it is a real merger of baryonic component
        !
        ! HOT ATMOSPHERE I **************************************************************
        !
        Mhot_1 = gas_mass(bh1%hot_gas) 
        Mhot_2 = gas_mass(bh2%hot_gas)
        M_hot  = Mhot_1 + Mhot_2     
        mu_gas = min(Mhot_1,Mhot_2)/max(Mhot_1,Mhot_2)
        !
        ! quenched cases
        ! the status is leaded by the most massive structure
        if (Mhot_1 .gt. Mhot_2) then
           !
           bh1%quenched = bh1%quenched
        else
           !
           bh1%quenched = bh2%quenched
        end if
        !
        ! sum hot gas component
        call gas_add(bh1%hot_gas,bh2%hot_gas)   
        ! inject thermal energy due to the merger
        Eint = gravconst_code_unit*(Mhot_1*Mhot_2)/dm1%R_vir
        call gas_inject_thermal_energy(bh1%hot_gas,Eint)
        !
        ! update cooling clock
        if (M_hot .gt. 0.d0) then 
           !
           ! the cooling timescale is led by the most massive structure
           if (Mhot_1 .gt. Mhot_2) then
              !
              bh1%t_cool = bh1%t_cool
           else
              !
              bh1%t_cool = bh2%t_cool
           end if
           !
#ifdef THERMAL_INSTABILITIES
! -------------------------------------------------   
           ! As for the cooling timescale,
           ! the thermal instability timescale is led by the most massive structure
           if (Mhot_1 .gt. Mhot_2) then
              !
              bh1%t_TI = bh1%t_TI
           else
              !
              bh1%t_TI = bh2%t_TI
           end if    
! -------------------------------------------------
#endif 
        end if   
        !
#ifdef COLD_STREAMS
! -------------------------------------------------   
        !   
        ! COLD FILAMENTARY STRUCTURE ****************************************************
        !
        ! sum cold gas component
        call gas_add(bh1%cold_gas,bh2%cold_gas) 
        ! 
        ! SHOCK HEATED GAS **************************************************************
        !
        ! Respect to Lu et al a fraction of the cold stream gas is disrupted and converted in hot gas
        call gas_void(dM_hot)
        f_sh   = baryon_halo_shock_heated_fraction(dm1%M_acc)
        M_cold = gas_mass(bh1%cold_gas) 
        M_tot  = M_cold + M_hot
        !
        if ((M_cold .gt. 0.d0) .and. (M_hot .gt. 0.d0) .and. (f_sh .lt. 1.d0)) then
            !
            f_sh   = max(0.d0,M_cold - max(0.d0,min(1.d0,(1.d0-f_sh)))*M_tot)/M_cold
            ! computed the disrupted mass of cold stream
            dM_hot = f_sh*bh1%cold_gas
            ! substract to the cold phase
            ! gas_sub(g1,g2[,therm_mode,called_by])
            call gas_sub(bh1%cold_gas,dM_hot,therm_mode='isothermal')
            ! heat the gas, set temperature to bh1 hot gas temperature
            call gas_set_component(dM_hot,max(diffuse_gas_temp,gas_temp(bh1%hot_gas)),component='Temp')
            ! add to the hot atmosphere
            call gas_add(bh1%hot_gas,dM_hot)
        end if
! -------------------------------------------------
#endif
! COLD_STREAMS 
        !
        ! HOT ATMOSPHERE II *************************************************************
        !
        ! checks temperature
        if ((gas_mass(bh1%hot_gas) .gt. M_gas_crit) .and. (gas_temp(bh1%hot_gas) .lt. 9.9d-1*diffuse_gas_temp)) then
            !
            call IO_print_warning_message('T_hot < T_cool_gas',only_rank=rank,called_by='baryon_halo_merge')
#ifdef PRINT_WARNING
! -------------------------------------------------         
            call IO_print_message('with',only_rank=rank,component='bh', &
                    param_name = (/'M_hot                    ','T_hot                    '/), &
                    real_param_val  = (/gas_mass(bh1%hot_gas),gas_temp(bh1%hot_gas)/))
#endif
! -------------------------------------------------     
! PRINT_WARNING 
            call gas_set_component(bh1%hot_gas,diffuse_gas_temp,component='Temp')
        end if        
        !
        ! evaporation
        Vesc = dm_escape_velocity(dm1%r_core,dm1)
        call gas_evap(bh1%hot_gas,m_evap,Vesc)
        ! in case of merger, even if the two progenitors are sub-structures, all the gas ejected is added to the IGM
        call gas_add_igm_gas(m_evap,called_by='baryon_halo_merge / m_evap')  ! feed the IGM
        ! update density profile parameters
        call baryon_halo_update_density_profile(bh1,dm1,called_by='baryon_halo_merge')  
        !
      end if
    else    
       !
       if (baryon_halo_bh_mass(bh2) .gt. 0.d0) then
         ! it the first copy in the descendent merger tree halo
         call baryon_halo_copy(bh1,bh2)
       end if
    end if
    !
    ! CHECKS
    if (is_NaN(gas_mass(bh1%hot_gas))) then
        !
        call IO_print_error_message('hot_gas is NaN ',only_rank=rank,called_by='baryon_halo_merge')
        stop ! stop the program
    end if
    if (is_NaN(gas_mass(bh1%cold_gas))) then
        !
        call IO_print_error_message('cold_gas is NaN ',only_rank=rank,called_by='baryon_halo_merge')
        stop ! stop the program
    end if

#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('baryon_halo_merge ... done',only_rank=rank,component='bh')
! -------------------------------------------------
#endif
    
    return
  end subroutine baryon_halo_merge

  !*****************************************************************************************************************
  !
  ! FUNCTIONS
  !
  !***************************************************************************************************************** 
  
  function baryon_halo_escape_fraction(bh)
    
    ! RETURN THE ESCAPE FRACTION OF THE HOT HALO GAS PHASE

    implicit none

    real(kind=8)                        :: baryon_halo_escape_fraction

    type(baryon_halo_type),intent(in)   :: bh  ! the baryon halo component
    
    baryon_halo_escape_fraction = bh%escape_fraction

    return
  end function baryon_halo_escape_fraction
  
  !***************************************************************************************************************** 
  
  function baryon_halo_hot_gas_fraction(bh)
    
    ! RETURN THE HOT GAS FRACTION

    implicit none

    type(baryon_halo_type),intent(in)   :: bh  ! the baryon halo component
    real(kind=8)                        :: baryon_halo_hot_gas_fraction

    baryon_halo_hot_gas_fraction = 0.d0
 
    if (baryon_halo_bh_mass(bh) .gt. 0.d0) then
       !
       baryon_halo_hot_gas_fraction = gas_mass(bh%hot_gas) / baryon_halo_bh_mass(bh) 
    end if
       
    return
  end function baryon_halo_hot_gas_fraction

  !*****************************************************************************************************************

  function baryon_halo_bh_mass(bh,component)
    
    ! RETURN THE MASS OF A BARYON HALO COMPONENT

    implicit none 

    character(*),intent(in),optional  :: component
    character(MAXPATHSIZE)            :: message        ! a message to display
    
    real(kind=8)                      :: baryon_halo_bh_mass

    type(baryon_halo_type),intent(in) :: bh
    
    if (present(component)) then
       !
       select case (trim(component)) 
       case ('cold')
          baryon_halo_bh_mass = gas_mass(bh%cold_gas)
       case ('hot')
          baryon_halo_bh_mass = gas_mass(bh%hot_gas)
       case('surrounding')
          baryon_halo_bh_mass = gas_mass(bh%surrounding_gas)
       case ('cold+hot','hot+cold')
          baryon_halo_bh_mass = gas_mass(bh%cold_gas + bh%hot_gas)
       case default
        write(message,'(a,a,a)') 'Keyword ', trim(component), ' not defined'
        call IO_print_error_message(message,only_rank=rank,called_by='baryon_halo_bh_mass')    
        stop  ! stop the program
       end select
    else
       !
       baryon_halo_bh_mass = gas_mass(bh%cold_gas + bh%hot_gas)
    end if
    
    return
  end function baryon_halo_bh_mass

  !*****************************************************************************************************************
  
  function baryon_halo_shock_heated_fraction(M_halo)
    
    ! RETURN THE FRACTION OF HOT GAS ACCRETED BY A HALO OF MASS M

    implicit none
    
    real(kind=8)                     :: baryon_halo_shock_heated_fraction  ! fraction of the accreted gas that is shock heated
    real(kind=8),intent(in)          :: M_halo                             ! in a halo of mass M_halo
    real(kind=8),parameter           :: m_tran = 1.d0                      ! in code unit
    real(kind=8),parameter           :: sigma_logm = 1.d0   !0.4d0
    real(kind=8)                     :: x                                  ! x << 0 M<<M_crit, x>0 M>>M_crit
                                                                           
    x  = log(M_halo/m_tran)/sigma_logm 

    baryon_halo_shock_heated_fraction = 5.d-1*(1.0d0+erf_func(x))  
    ! cut for adaptative time step facility
    if (baryon_halo_shock_heated_fraction .lt. 1.d-7) baryon_halo_shock_heated_fraction = 0.d0
    
    return
  end function baryon_halo_shock_heated_fraction

  !*****************************************************************************************************************

  function f_hot_gas(x,b)
    
    ! RETURN THE HOT GAS ADDIMENTIONNAL FUNCTION AT x
    ! All hot gas profile model are extract to : Capelo+10, Mokino+98, Suto+98
    ! f_hot_gas is defined by rho = rho(r=0)*f_hot_gas
    ! we have two different model of hot gas profile : isothermal or polytropic case
    ! in all cases : f_hot_gas is an : 
    !                  - decreasing function of x
    !                  - decresing function of b 
    !                  - increasing function of T

    implicit none
    
    real(kind=8), intent(in)     :: x,b         ! x is the adimentionnal quantity r/dm%r_core
    real(kind=8)                 :: f_hot_gas
    
#ifdef POLYTROPIC
! -------------------------------------------------
    f_hot_gas = (1.d0-b*f(x))**(1.d0/(gamma-1.d0))  ! Eq (42) and (46) Suto+98
#else
    f_hot_gas = exp(-b*f(x))                        ! Eq (8) in Mokino+98
! -------------------------------------------------
#endif 
! POLYTROPIC
    
    return
  end function f_hot_gas

  !*****************************************************************************************************************

  function f(x)

    ! RETURN GEOMETRICAL PROFILE FUNCTION 
    ! Eq (14) suto+98  
    
    implicit none

    real(kind=8),intent(in)    :: x
    real(kind=8)               :: f

    f = 1.d0-log(1.d0+x)/x    

    return
  end function f
  
  !*****************************************************************************************************************

  function b_const(T0,dm)
    
    ! All hot gas profile model are extract to : Capelo+10, Mokino+98, Suto+98
    ! b_const is an : 
    !       - decreasing function of central temperature T0
    !       - increasing function of the core dark-matter radius
    !       - increasing function of the dark-matter core density
    
    implicit none
    
    real(kind=8),intent(in)    :: T0       ! central temperature of the gas
                                           ! in case of perfect gas, the central temperature T0 is equal to the mean temperature T_hot
                                           ! it is not true in case of polytropic gas
    real(kind=8)               :: b_const  ! a constant defined in Eq (8) Suto+98 or Eq. (9) Mokino+98 or Eq (44) Suto+98
    
    type(dm_type),intent(in)   :: dm       ! the dark matter structure
    
    if (T0 .le. 0.d0) then
       !
       call IO_print_error_message('T0 <= 0',only_rank=rank,called_by='b_const')
       stop ! stop the program
    end if
    
    ! perfect gas case
    ! Eq (8) Suto+98 or Eq. (9) Mokino+98
    b_const = 4.d0*pi*gravconst_code_unit*dm%rho_core*dm%r_core**2./(r_perfect_gas*T0)  
#ifdef POLYTROPIC
! -------------------------------------------------
    ! polytropic gas case
    ! gamma is the polytropic index of the gas gamma = 1 + 1/n Eq (44) Suto+98
    b_const = b_const*(gamma-1.d0)/gamma      
! -------------------------------------------------
#endif
! POLYTROPIC
  
    return
  end function b_const

  !*****************************************************************************************************************
  
  function M(x,param)
    
    ! RETURN ADDIMENTIONNAL FUNCTION FOR MASS AND T_MEAN COMPUTATIONS
    ! WARNING: This function depends, through b_const = param(1), of the dark-matter properties (r_core and rho_core) and of the central temperature of the hot atmosphere T0 
    
    implicit none
    
    real(kind=8), intent(in)             :: x
    real(kind=8),intent(in),dimension(2) :: param  ! parameter array
                                                   ! param(1) = b_const
                                                   ! param(2) = gamma (in case of polytropic computation for T_mean computation) or -1.d0 (otherwise)
    real(kind=8)                         :: M
   
    if (param(2) .eq. -1.d0) then
       !
       M = x**2.*f_hot_gas(x,param(1))
    else
       !
       M = x**2.*f_hot_gas(x,param(1))**param(2)
    end if
    
    return
  end function M
  
  !*****************************************************************************************************************
  
  function baryon_halo_hot_halo_density_profile(r,bh,dm)
    
    ! RETURN THE DENSITY OF THE HOT ATMOSPHERE AT RADIUS r
    ! All hot gas profile model are extract to : Capelo+10, Mokino+98, Suto+98
    ! the density profile is define through the function f_hot_gas: rho(x) = rho0*f_hot_gas(x)
    
     implicit none
     
     real(kind=8), intent(in)          :: r  ! radius
     real(kind=8)                      :: baryon_halo_hot_halo_density_profile,b,x
     
     type(baryon_halo_type),intent(in) :: bh ! a baryon halo component     
     type(dm_type),intent(in)          :: dm ! a dark matter component
     
     if (bh%T0 .le. 0.d0) then
        !
        call IO_print_error_message('T0 <= 0',only_rank=rank,called_by='baryon_halo_hot_halo_density_profile')
        stop ! stop the program
     end if
     
     x = r/dm%r_core
     b = b_const(bh%T0,dm)
     
     baryon_halo_hot_halo_density_profile = bh%rho0*f_hot_gas(x,b)
     
     return
   end function baryon_halo_hot_halo_density_profile
   
   !*****************************************************************************************************************
   
   function baryon_halo_hot_halo_pressure_profile(r,bh,dm)
     
     ! RETURN THE PRESSURE OF THE HOT ATMOSPHERE AT RADIUS r
     ! All hot gas profile model are extract to : Capelo+10, Mokino+98, Suto+98
     
     implicit none
     
     real(kind=8),intent(in)           :: r    ! radius
     real(kind=8)                      :: baryon_halo_hot_halo_pressure_profile,b,x
     real(kind=8)                      :: P,P0

     type(baryon_halo_type),intent(in) :: bh   ! a baryon halo component
     type(dm_type),intent(in)          :: dm   ! the dark-matter component
     
     if (bh%T0 .le. 0.d0) then
        !
        call IO_print_error_message('T0 <= 0',only_rank=rank,called_by='baryon_halo_hot_halo_pressure_profile')
        stop ! stop the program
     end if
     
     x = r/dm%r_core
     b = b_const(bh%T0,dm)
     P0 = bh%rho0*r_perfect_gas*bh%T0

#ifdef POLYTROPIC
! -------------------------------------------------
     ! polytropic case
     P  = P0*(f_hot_gas(x,b))**(gamma)   ! Eq (39), (42) et (46) Suto+98
#else
     ! perfect gas case
     ! in the case of a perfect isothermal gas, the central temperature of the hot atmosphere is equal to the mean temperature T_hot
     P = P0*f_hot_gas(x,b)
! -------------------------------------------------
#endif
! POLYTROPIC

     baryon_halo_hot_halo_pressure_profile = P
     
     return
   end function baryon_halo_hot_halo_pressure_profile
      
   !*****************************************************************************************************************
   
   function baryon_halo_hot_halo_temperature_profile(r,bh,dm)
     
     ! RETURN THE TEMPERATURE OF THE HOT HALO PHASE AT RADIUS r
     ! All hot gas profile model are extract to : Capelo+10, Mokino+98, Suto+98
     
     implicit none
     
     real(kind=8), intent(in)          :: r    ! radius
     real(kind=8)                      :: baryon_halo_hot_halo_temperature_profile,b,x,T
     
     type(baryon_halo_type),intent(in) :: bh   ! a baryon halo component
     type(dm_type),intent(in)          :: dm   ! a dark matter component
     
     if (bh%T0 .le. 0.d0) then
        !
        call IO_print_error_message('T0 <= 0',only_rank=rank,called_by='baryon_halo_hot_halo_temperature_profile')
        stop ! stop the program
     end if
     
     x = r/dm%r_core
     b = b_const(bh%T0,dm)
     
#ifdef POLYTROPIC
! -------------------------------------------------
     ! polytropic gas case
     ! Eq (40), (42) et (46) Suto+98
     T = bh%T0*(f_hot_gas(x,b))**(gamma-1.d0)    
#else
     ! perfect isothermal gas case
     T = gas_temp(bh%hot_gas)  ! = bh%T0
! -------------------------------------------------
#endif
! POLYTROPIC
     
     baryon_halo_hot_halo_temperature_profile = T
     
     return
   end function baryon_halo_hot_halo_temperature_profile
   
   !*****************************************************************************************************************
   
  function baryon_halo_integrate_hot_halo_density_profile(r_out,bh,dm,error,final_order,called_by)

    ! RETURN THE INTEGRETED MASS CONTAINS IN THE HOT HALO PHASE between 0 and r_out

     implicit none

     integer(kind=4),intent(out),optional :: final_order   ! final order of the integrator
     
     character(*),intent(in),optional     :: called_by     ! name of the function which has called this function
     character(MAXPATHSIZE)               :: message       ! a message to display
     
     real(kind=8), intent(in)             :: r_out         ! external radius 
     real(kind=8),intent(out),optional    :: error         ! integration error
     real(kind=8)                         :: mass, baryon_halo_integrate_hot_halo_density_profile
     real(kind=8)                         :: x_out         ! adimentional outer radii
     real(kind=8),dimension(2)            :: param         ! a parameter array 

     type(baryon_halo_type),intent(in)    :: bh            ! a baryon halo component
     type(dm_type),intent(in)             :: dm            ! a dark-matter component

     if (present(called_by)) then
        write(message,'(a,a)') 'Function called by : ', trim(called_by)
     end if
     
     baryon_halo_integrate_hot_halo_density_profile = 0.d0  ! init
     
     x_out = r_out/dm%r_core                                ! compute adimential external radius

     if (x_out .le. x_cut) return                           ! too small radius
     
     if (bh%T0 .le. 0.d0) then
        !
        call IO_print_error_message('T0 <= 0',only_rank=rank,called_by='baryon_halo_integrate_hot_halo_density_profile')
        if (present(called_by)) call IO_print_message(trim(message),only_rank=rank,component='bh')
        call IO_print_message('used',only_rank=rank,component='bh', &
             param_name=(/'r_out                    ','x_out                    ','rho0                     ', &
                          'T0                       ','t_cool                   ','r_cool                   '/), &
             real_param_val=(/r_out,x_out,bh%rho0,bh%T0,bh%t_cool,bh%r_cool/))
        stop ! stop the program
     end if
     
     if ((bh%rho0 .lt. 0.d0) .or. (bh%T0 .lt. 0.d0)) then
        !
        call IO_print_error_message('Wrong parameter values',only_rank=rank,called_by='baryon_halo_integrate_hot_halo_density_profile')
        if (present(called_by)) call IO_print_message(trim(message),only_rank=rank,component='bh')
        call IO_print_message('Computation used',only_rank=rank,component='bh', &
             param_name=(/'r_out                    ','x_out                    ','rho0                     ','T0                       '/), &
             real_param_val=(/r_out,x_out,bh%rho0,bh%T0/))
        stop ! stop the program
     end if

     param = (/b_const(bh%T0,dm),-1.d0/) 

     ! compute mass enclose between x_in and x_out and add fixed value in the M(x<x_in) = f_hot_gas(x_in)x**2  
     mass = 4.d0*pi*bh%rho0*dm%r_core**3.*Ronbint(M, x_cut, x_out, param, error = error, final_order = final_order, &
                                             called_by = 'baryon_halo_integrate_hot_halo_density_profile (mass)')

     if (mass .le. 0.d0) then
        !
        call IO_print_error_message('mass <= 0.0',only_rank=rank,called_by='baryon_halo_integrate_hot_halo_density_profile')
        if (present(called_by)) call IO_print_message(trim(message),only_rank=rank,component='bh')
        call IO_print_message('Computation used',only_rank=rank,component='bh', &
             param_name=(/'mass                     ','error                    ','r_out                    ','x_out                    ', &
                          '(dm)r_core               ','rho0                     ','T0                       '/), &
             real_param_val=(/mass,error,r_out,x_out,dm%r_core,bh%rho0,bh%T0/))
        stop
     endif

     baryon_halo_integrate_hot_halo_density_profile = mass

     return
  end function baryon_halo_integrate_hot_halo_density_profile
    
  !*****************************************************************************************************************

#ifdef POLYTROPIC
! -------------------------------------------------
  function T0_crit_polytropic(dm)

    ! RETURN THE CRITICAL CENTRAL TEMPERATURE OF THE POLYTROPIC GEOMETRICAL FUNCTION
    ! The geometrical function f_hot_gas is null at r = +inf for T = T_crit_polytropic

    implicit none
    
    real(kind=8)                :: T0_crit_polytropic  ! minimal temperature
    
    type(dm_type),intent(in)    :: dm                  ! the dark-matter component

    T0_crit_polytropic = gamma/(gamma-1.d0)*4.d0*pi*gravconst_code_unit*dm%rho_core*dm%r_core**2./r_perfect_gas ! Eq (48) Suto+98

    return
  end function T0_crit_polytropic 

  !*****************************************************************************************************************

  function T_mean_polytropic(bh,dm)

    ! RETURN THE MEAN TEMPERATURE OF A HOT HALO 
    ! in the polytropic case T_mean(T0) is a increasing function of T0 for a given set of dm properties

    implicit none
    
    real(kind=8)                         :: err           ! integration process error 
    real(kind=8)                         :: den
    real(kind=8)                         :: T_mean_polytropic

    type(baryon_halo_type),intent(in)    :: bh            ! a baryon halo component
    type(dm_type),intent(in)             :: dm            ! a dark-matter component
   
    T_mean_polytropic = -1.d0 ! init

    if (dm%concentration .lt. x_cut) then
       !
       call IO_print_error_message('dm%concentration < x_cut',only_rank=rank,called_by='T_mean_polytropic')
       stop ! stop the program
    end if

    if (bh%T0 .le. 0.d0) then
       !
       call IO_print_error_message('bh%T0 <= 0',only_rank=rank,called_by='T_mean_polytropic')
       stop ! stop the program
    end if
    
    if (gas_mass(bh%hot_gas) .le. 0.d0) return

    den = Ronbint(M,x_cut,dm%concentration,err,(/b_const(bh%T0,dm),-1.d0/), called_by = 'T_mean_polytropic (den)')
    T_mean_polytropic = bh%T0*Ronbint(M,x_cut,dm%concentration,err,(/b_const(bh%T0,dm),gamma/), called_by = 'T_mean_polytropic (T)')/den
    
    return
  end function T_mean_polytropic

  !*****************************************************************************************************************

  function T0_polytropic(bh,dm)

    ! RETURN THE CENTRAL TEMPERATURE OF A POLYTROPIC GAS PHASE IN FUNCTION OF ITS MEAN TEMPERATURE
    ! In the polytropic case  : 
    !                            - T_hot the mean temperature (Tmean) of the gas is not equal to the central temperature (T0) 
    !                            - The mean temperature (Tmean) is an increasing function of the central temperaure (T0)
    !                            - The mean temperature is always smaller than the central temperature because, the geometrical fucntion, f_hot_gas is an : 
    !                                                             - decreasing function of x
    !                                                             - decresing function of b 
    !                                                             - increasing function of T0
    ! We keep all the hot gas in the virial radius (because we compute escape fraction of the gas), therefore we must have, T(r=r_vir|Tmean) = T(x=c) > 0
    ! For all (previous) reasons if T(x=c|Tmean) > 0, then T(x=c,|T0) > T(x=c|Tmean) > 0, The polytropic profile is consistent and can be computed

    implicit none

    integer(kind=4)                   :: loop                            ! loop index

    real(kind=8)                      :: T0_polytropic                   ! the central temperature of the hot halo phase (-2 if no stable hot atmosphere)
    real(kind=8)                      :: lT0, lT0max, lT0min             ! local temperature variables (mainly log scale defintions)
    real(kind=8)                      :: T, T0, T_hot                    ! other local temperature variables 

    type(baryon_halo_type),intent(in) :: bh                              ! a baryon halo component
    type(baryon_halo_type)            :: bh_tmp                          ! a tmp baryon halo component
    type(dm_type),intent(in)          :: dm                              ! the dark matter component

    T0_polytropic = -1.d0 ! init

    if (baryon_halo_bh_mass(bh,component='hot') .le. 0.d0) return

    ! create a local copy of the baryon halo component bh
    call baryon_halo_copy(bh_tmp,bh)

    ! set the central temperature to the minimal authorized central temperature
    bh_tmp%T0 =  T0_crit_polytropic(dm)    ! central temperature such as T(x=c) = 0 
    !
    T_hot = gas_temp(bh%hot_gas) 
    T     = T_mean_polytropic(bh_tmp,dm)
    
    if (FOLLOW_UP .and. PR_FOLLOW_UP) write(*,*) 'T0, T_hot, T: ', bh_tmp%T0, T_hot, T

    if (T .lt. T_hot) then
       !
       ! OK
       ! the minimal authorized central temperature (T0_crit_polytropic) leads to a mean temperature (T) smaller than the real mean temperature
       ! therefore, as the mean temperature T is an increasing function of the central temperature T0, 
       ! the real central temperature T0 (associated to the real mean temperaure T_hot) is larger than T0_crit_polytropic
       !
       ! search T0 such as T_mean_polytropic(T0) = T_hot, we use log scale
       lT0min    = log10(bh_tmp%T0)             ! init the lawer barrier (in log scale)
       lT0max    = 1.1d1                        ! init the higher barier (in log scale --> 10^11 K in real scale)
       lT0       = 5.d-1*(lT0max + lT0min)      ! init
       bh_tmp%T0 = 1.d1**(lT0)                  ! set 
       T         = T_mean_polytropic(bh_tmp,dm) ! init
       do while ((abs(T-T_hot))/T_hot .gt. num_precision)
          if (loop .gt. 1000) then
             !
             call IO_print_error_message('Too much cycle used in T0 computation',only_rank=rank,called_by='T0_polytropic')
             call IO_print_message('Computation used',only_rank=rank,component='bh', &
                  param_name=(/'T0_min                   ','T0_max                   ','T0                       ', &
                               'T                        ','T_hot                    ','M_hot                    ', &
                               '(dm)rho_core             ','(dm)r_core               ','rho0                     ', &
                               'R_halo                   ','(dm)concentration        '/), &
                  real_param_val=(/1.d1**(lT0min),1.d1**(lT0max),bh_tmp%T0,T,T_hot,gas_mass(bh%hot_gas), &
                              dm%rho_core,dm%r_core,dm%R_halo,dm%concentration/))
             stop ! stop the program
          end if
          ! Tmean is an increasing function of T0
          if (T .gt. T_hot) then
             !
             lT0max = lT0                 
          else
             !
             lT0min = lT0
          end if
          lT0 = 5.d-1*(lT0max + lT0min)    ! update
          bh_tmp%T0 = 1.d1**(lT0)          ! update
          T = T_mean_polytropic(bh_tmp,dm) ! update
          loop = loop +1                   ! increase loop index
       end do
       T0_polytropic = 1.d1**(lT0)
    else
       !
       ! The minimal authorized central temperature (T0_crit_polytropic) leads to a mean temperature (T) 
       ! higher than the real mean temperature (T_hot) 
       ! therefore, as the mean temperature T is an increasing function of the central temperature T0, 
       ! in these conditions the "hot" gas phase is not stable enought  
       T0_polytropic = -2.d0
    end if
    
    return
  end function T0_polytropic
! -------------------------------------------------
#endif  
! POLYTROPIC

  !*****************************************************************************************************************
  !
  ! PRINTING PROCEDURES
  !
  !*****************************************************************************************************************

    subroutine bh_load_bh_data(bh,bh_data)

    ! CREATE THE BH OUTPUT PROPERTIES LIST

    implicit none

    type(baryon_halo_type),intent(in)  :: bh    ! bh component
    real(kind=8), intent(inout)        :: bh_data(nb_bh_field)

    ! For information
    !'cold_gas              ','surrounding_gas       ','hot_gas               ',&
    !'hot_gas_mZ            ','hot_gas_mH            ','hot_gas_mC            ',&
    !'hot_gas_mN            ','hot_gas_mO            ','hot_gas_mFe           ',&
    !'acc_rate              ','ff_rate               ','cooling_rate          ',&
    !'T_hot                 ','t_cool                ','t_TI                  ',&
    !'r_cool                ','r_TI                  ','esc_frac              '
    
    bh_data = (/mass_code_unit_in_M_Sun*gas_mass(bh%cold_gas), &
                mass_code_unit_in_M_Sun*gas_mass(bh%surrounding_gas), &
                mass_code_unit_in_M_Sun*gas_mass(bh%hot_gas), &
                mass_code_unit_in_M_Sun*gas_mass(bh%hot_gas,component='Metals'), &
                mass_code_unit_in_M_Sun*gas_mass(bh%hot_gas,component='H1'), &
                mass_code_unit_in_M_Sun*gas_mass(bh%hot_gas,component='C12'), &
                mass_code_unit_in_M_Sun*gas_mass(bh%hot_gas,component='N14'), &
                mass_code_unit_in_M_Sun*gas_mass(bh%hot_gas,component='O16'), &
                mass_code_unit_in_M_Sun*gas_mass(bh%hot_gas,component='Fe56'), &
                mass_rate_code_unit_2_MsunPerYr*gas_mass(bh%cold_inflow_rate+bh%hot_inflow_rate), mass_rate_code_unit_2_MsunPerYr*gas_mass(bh%cold_outflow_rate), &
                mass_rate_code_unit_2_MsunPerYr*gas_mass(bh%cooling_rate), &
                gas_temp(bh%hot_gas), bh%t_cool, bh%t_TI, bh%r_cool, bh%r_TI, bh%escape_fraction/)

    return
  end subroutine bh_load_bh_data
  
  !*****************************************************************************************************************
  
  subroutine baryon_halo_print(unit,form,bh)
    
    ! PRINT BH-PROPERTIES IN eGALICS OUTPUT FILE
    
    implicit none
    
    integer(kind=4),intent(in)        :: unit                 ! file unit
    integer(kind=4)                   :: status,i,hdutype
     
    character(*)                      :: form                 ! fits or tmp_bin
    
    real(kind=8)                      :: bh_data(nb_bh_field) ! nb_bh_field is given in the header of this bh module
    
    type(baryon_halo_type),intent(in) :: bh                   ! the baryon_halo component
    
    call bh_load_bh_data(bh,bh_data)
    
    select case (trim(form))
      case ('tmp_bin')
        write(unit) bh_data  ! directly write data in the tmp binary output file 
      case ('fits')
        ! move to bh extension
        call ftmahd(unit,hdu_bh,hdutype,status) 
        if (status .gt. 0) then
          !
          call IO_print_error_message('ftmahd status', &
                only_rank = rank, called_by = 'bh_print')
          stop ! stop the program
        end if
        ! init
        call ftirow(unit,0,1,status) 
        if (status .gt. 0) then
          !
          call IO_print_error_message('ftirow status', &
                only_rank = rank, called_by = 'bh_print')
          stop ! stop the program
        end if
        ! write data in the dm entension
        call ftpcld(unit,1,1,1,1,bh%age_form+bh%life_time,status)
        do i=2, nb_bh_field+1   
          call ftpcld(unit,i,1,1,1,bh_data(i-1),status)  
          if (status .gt. 0) then
            !
            call IO_print_error_message('ftpcld status', &
                only_rank = rank, called_by = 'bh_print')
            stop ! stop the program
          end if
        end do
      case default
        call IO_print_error_message('Unknwon output data format', &
                only_rank = rank, called_by = 'bh_print')
        stop ! stop the program
    end select
    
    return
  end subroutine baryon_halo_print
  
  !***************************************************************************************************************`

end module baryon_halo
