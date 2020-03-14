module dm        

  use utilities        ! Acces to integrator Rombint, Bessel function, .... 

  public

  !*****************************************************************************************************************
  ! 
  !  OVERVIEW
  !
  !  dm module defines the dm (dark matter) data structure univ(ts)%halo(ih)%dm
  !
  !  This module contains the complete information about the physical properties of a dark matter halo
  !  A DM halo is described by a NFW model truncated at R_vir
  !  The dark matter halo mass can be followed according to its total instantaneous mass,
  !  its Viral mass and its smoothly accreted mass
  !  In the header of this module are defined output properties of a dm component (labels, units and formats)
  !
  !  SUBROUTINES IN THIS MODULE
  !
  !   dm_set_reference_mass         : initialize some mass parameter linked to a dark-matter component
  !
  !   dm_void                       : init all properties of the dm component
  !      called by : halo_void
  !  
  !   dm_copy                       : copy the dm properties of a halo in an other
  !      called by : halo_save
  !  
  !   dm_set_properties             : set dm the properties directly measured from the N-body simulation
  !      called by : halo_set_properties
  !
  !   dm_set_M_acc                  : set the value of the M_acc (integreted mass onto progenitors) for a dm component
  !      called by : tree_compute_halo_derived_properties
  !
  !   dm_flag_bad_properties        : flag (dm%compute = .false.) halos with bad dark matter properties (total energy, spin, null or negative integreted mass)
  !      called by : tree_compute_halo_derived_properties
  !
  !   dm_set_merger_mass            : set the merger mass 
  !      called by : tree_compute_halo_derived_properties
  ! 
  !   dm_set_merger_rate            : set the merger rate
  !      called by : tree_compute_halo_derived_properties
  !
  !   dm_set_bad_branch_mass        : set the bad branch mass
  !      called by : tree_compute_halo_derived_properties
  !
  !   dm_set_acc_rate_left          : set the accretion rate left (past) of the dark matter content halo
  !      called by : tree_compute_halo_derived_properties
  !
  !   dm_set_acc_rate_right         : set the accretion rate right (future) of the dark matter content halo. This value is used for compute 
  !                                      the evolution of the dark matter mass between the current time-step and the next merging time. 
  !     called by : tree_compute_halo_derived_properties
  !
  !   dm_set_age_form               : set the formation age (age of the universe when the dark matter structure (its oldest progenitor) has been formed)
  !      called by : tree_compute_halo_derived_properties
  !
  !   dm_set_age_merge              : set the merging age (age of the univers when the structure(halo+galaxy) have been generated from its progenitors) 
  !      called by : tree_compute_halo_derived_properties
  !
  !   dm_compute_derived_properties : compute properties of a halo that can't be directly measured in the n-body simulation
  !      called by : tree_compute_halo_derived_properties
  !
  !   dm_update_life_time           : update the evolving time of the dark-matter halo
  !      called by :
  !
  !  FUNCTIONS IN THIS MODULE
  !
  !   dm_density_profile
  !
  !   dm_f_nfw                      : is a function such that 4*pi*rho_core*r_core3**3*dm_f_nfw(x) is the dm mass enclosed within r = x * r_core  
  ! 
  !   dm_concentration              : returns the typical concentration of a halo of mass M_halo at redshift z    
  !  
  !   dm_mass                       : returns the mass encolse in a sphere of radius r (for a NFW profile)
  !
  !   dm_potential                  : returns the gravitationnal potential of the dark matter profile at a radius r
  !
  !   dm_velocity (INTERFACE VERS)  : returns the velocity of the dark matter profile at a radius r
  !
  !   dm_velocity_ (CORE VERS)
  !
  !   dm_escape_velocity            : returns the escape velocity of the dark matter component at a radius r  
  !
  !   dm_orbital_radius             : return the distance between two dark-matter halos
  !
  !   dm_dynamical_friction_time    : return the friction dynamical time associated to a pair of dark-matter structures
  !
  !  PRINTING PROCEDURES
  !
  !   dm_load_dm_data               : load dm data, create the dark-matter output property list 
  !      called by : dm print procedures (e.g. dm_print)
  !    
  !   dm_print                      : print dm component properties of the halo
  !      called by : halo_print
  !
  !*****************************************************************************************************************

  type dm_type 
    ! Evolution
    real(kind=8)    :: life_time            ! evolve time of the dark-matter structure
    real(kind=8)    :: age_form             ! age of the universe when the structure has been formed
    real(kind=8)    :: age_merge            ! age of the univers when the structure(halo+galaxy) have been generated from its progenitors 
    !
    ! dark-matter simulation space and time information
    integer(kind=4) :: ts                   ! timestep
    real(kind=8)    :: pos(3)               ! space coordinates of the dark matter halo in the N-body simulation
    real(kind=8)    :: vel(3)               ! cartesian composante of the velocity of the centre of mass of the dark-matter halo
    real(kind=8)    :: L(3)                 ! angular momentum (Lx, Ly, Lz)
    real(kind=8)    :: Ek_Ep                ! ratio of kinetic to potential energies
    ! cleaning
    logical         :: compute              ! = .true. if the dark matter properties are good
    real(kind=8)    :: goodness_of_branch   ! nb of (dm-)good halos / nb of halos on branch
    ! merger and accretion dynamics
    real(kind=8)    :: dmacc                ! mass of particles that were never part of any halo at any previous tstep (10^11 M_Sun)
    real(kind=8)    :: bad_branch_mass      ! integreted mass of all bad branch (halo%compute = .false.) connected to this halo 
    real(kind=8)    :: merger_mass          ! mass gained by merger along the time evolution of the branch (sum of Mhalo over all non-main progenitor)
    real(kind=8)    :: acc_rate_left        ! accretion rate of the descendent for the evolution between le merger and the next time-step
    real(kind=8)    :: acc_rate_right       ! fraction of the descendent accretion rate for compute the evolution between the last time-step and the merging time
    real(kind=8)    :: merger_rate          ! merger rate
    ! masses
    real(kind=8)    :: M_vir                ! virial mass (10^11 M_Sun) measured from the N-body simulation
    real(kind=8)    :: M_tot                ! total mass (i.e. M_FOF or M_HOP) measured from the N-body simulation
    real(kind=8)    :: M_acc                ! integrated halo mass (sum of dmacc over all progs at all times)
    real(kind=8)    :: M_halo               ! The halo mass used is the baryonic calculation, may be M_tot, M_vir or M_acc  
    ! radius
    real(kind=8)    :: R_vir                ! virial radius (kpc) measured from the N-body simulation
    real(kind=8)    :: R_200                ! virial radius (kpc) computed from M_halo with the Bryan & Norman (1997) fitting formulae
    real(kind=8)    :: R_halo               ! the halo radius used is the baryonic calculation, may be R_vir or R_200 
    ! spin
    real(kind=8)    :: spin                 ! spin parameter (dimensionless) 
    ! halo parameters
    real(kind=8)    :: r_core               ! core radius (kpc)
    real(kind=8)    :: rho_core             ! core density (10^11MSun/kpc**3)
    real(kind=8)    :: V_vir                ! virial velocity (km/s)
    real(kind=8)    :: T_vir                ! virial temperature (K)
    real(kind=8)    :: t_dyn                ! dynamical time (R/V) (Gyr)
    real(kind=8)    :: concentration        ! ratio of virial to core radius

  end type dm_type
  
  ! DEFINE HEADER INFORMATIONS *********************

  ! hdu reference for dm structure
  integer(kind=4)           :: hdu_dm
  ! printable properties for dm structure
  integer(kind=4),parameter :: nb_dm_field = 21 ! number of dm properties saved
  ! Name of each output colomn data
  character(len=ttype_len),dimension(nb_dm_field) :: ttype_dm = (/'x                     ','y                     ','z                     ',&
                                                                  'Vx                    ','Vy                    ','Vz                    ',&
                                                                  'M_vir                 ','M_tot                 ','M_acc                 ',&
                                                                  'R_vir                 ','R_200                 ','spin                  ',&
                                                                  'acc_rate_left         ','acc_rate_right        ','V_vir                 ',&
                                                                  'T_vir                 ','t_dyn                 ','concentration         ',&
                                                                  'r_core                ','dm_rho_core           ','age_form              '/)   
  ! Physical unit of each column data
  character(len=tunit_len),dimension(nb_dm_field) :: tunit_dm = (/'Mpc         ','Mpc         ','Mpc         ','km/s        ','km/s        ','km/s        ', &
                                                                  'M_sun       ','M_sun       ','M_sun       ','kpc         ','kpc         ','w_o_unit    ',&
                                                                  'M_sun/yr    ','M_sun/yr    ','km/s        ','K           ','Gyr         ','w_o_unit    ',&
                                                                  'kpc         ','M_sun/kpc^3 ','Gyr         '/)  
  ! Data type of each column data
  character(len=tform_len),dimension(nb_dm_field) :: tform_dm = (/'1E','1E','1E','1E','1E','1E','1E','1E','1E','1E','1E',&
                                                                  '1E','1E','1E','1E','1E','1E','1E','1E','1E','1E'/)  

contains
  
  !*****************************************************************************************************************
  !
  ! SUBROUTINES
  !
  !*****************************************************************************************************************

  subroutine dm_set_reference_mass
    
    implicit none
 
    M_tot_min = 2.d1*dm_particle_mass
    
    call IO_print_message('use',param_name=(/'M_tot_min [Msun]         '/),real_param_val=(/M_tot_min*mass_code_unit_in_M_Sun/))
    
    return
  end subroutine dm_set_reference_mass
  
  !*****************************************************************************************************************
  
  subroutine dm_void(dm)

    ! INITIALIZE dm COMPONENT TO NULL

    implicit none

    type(dm_type),intent(inout) :: dm
    
    ! Evolution
    dm%life_time          = 0.d0               ! will be update in dm_update_life_time
    dm%age_form           = -1.d0              ! will be set in halo_evolve at the beginning of the first evolution step
    !
    ! dark-matter simulation space and time information
    dm%ts                 = 0                  ! will be set in dm_set_properties
    dm%pos                = (/0.d0,0.d0,0.d0/) ! Will be read in input (dark-matter simulation analysis), will be set in dm_set_properties
    dm%vel                = (/0.d0,0.d0,0.d0/) ! Will be read in input (dark-matter simulation analysis), will be set in dm_set_properties
    dm%L                  = (/0.d0,0.d0,0.d0/) ! Will be read in input (dark-matter simulation analysis), will be set in dm_set_properties
    dm%Ek_Ep              = 0.d0               ! Will be read in input (dark-matter simulation analysis), will be set in dm_set_properties
    !
    ! cleaning
    dm%compute            = .true.             ! will be set in dm_flag_bad_properties
    dm%goodness_of_branch = -1.                ! will be compute in tree_set_dm_goodness_of_branch
    !
    ! merger and accretion dynamics
    dm%dmacc              = 0.d0               ! Will be read in input (dark-matter simulation analysis), will be set in dm_set_properties
    dm%bad_branch_mass    = 0.d0               ! Will be computed in tree_compute_halo_derived_properties and set in dm_set_bad_branch_mass
    dm%merger_mass        = 0.d0               ! will be computed in tree_compute_halo_derived_properties and set in dm_set_merger_mass
    dm%acc_rate_left      = 0.d0               ! will be computed in tree_compute_halo_derived_properties and set in dm_set_acc_rate_left
    dm%acc_rate_right     = 0.d0               ! will be computed in tree_compute_halo_derived_properties and set in dm_set_acc_rate_right
    dm%merger_rate        = 0.d0               ! will be computed in tree_compute_halo_derived_properties and set in dm_set_merger_rate
    !
    ! masses
    dm%M_vir              = 0.d0               ! Will be read in input (dark-matter simulation analysis)
    dm%M_tot              = 0.d0               ! Will be read in input (dark-matter simulation analysis)
    dm%M_acc              = 0.d0               ! Will be computed in tree_compute_halo_derived_properties and set in dm_set_M_acc
    dm%M_halo             = 0.d0               ! Will be set in dm_compute_derived_properties
    !
    ! radius
    dm%R_vir              = 0.d0               ! Will be read in input (dark-matter simulation analysis), will be set in dm_set_properties
    dm%R_200              = 0.d0               ! Will be computed in dm_compute_derived_properties
    dm%R_halo             = 0.d0               ! Will be computed in dm_compute_derived_properties
    !
    ! spin
    dm%spin               = 0.d0               ! ! Will be read in input (dark-matter simulation analysis), will be set in dm_set_properties
    !
    ! halo parameters
    dm%r_core             = 0.d0               ! Will be computed in dm_compute_derived_properties
    dm%rho_core           = 0.d0               ! Will be computed in dm_compute_derived_properties
    dm%V_vir              = 0.d0               ! Will be computed in dm_compute_derived_properties
    dm%T_vir              = 0.d0               ! Will be computed in dm_compute_derived_properties
    dm%t_dyn              = 0.d0               ! Will be computed in dm_compute_derived_properties
    dm%concentration      = 0.d0               ! Will be computed in dm_compute_derived_properties
    
    return
  end subroutine dm_void
  
  !*****************************************************************************************************************
  
  subroutine dm_copy(dm1,dm2)

    implicit none

    type(dm_type),intent(inout) :: dm1
    type(dm_type),intent(in)    :: dm2 

    ! Evolution
    dm1%life_time          = dm2%life_time
    dm1%age_form           = dm2%age_form
    !
    ! dark-matter simulation space and time information
    dm1%ts                 = dm2%ts 
    dm1%pos                = dm2%pos 
    dm1%vel                = dm2%vel
    dm1%L                  = dm2%L    
    dm1%Ek_Ep              = dm2%Ek_Ep   
    ! 
    ! cleaning
    dm1%compute            = dm2%compute                      
    dm1%goodness_of_branch = dm2%goodness_of_branch
    !
    ! merger and accretion dynamics
    dm1%dmacc              = dm2%dmacc 
    dm1%bad_branch_mass    = dm2%bad_branch_mass
    dm1%merger_mass        = dm2%merger_mass
    dm1%acc_rate_left      = dm2%acc_rate_left  
    dm1%acc_rate_right     = dm2%acc_rate_right  
    dm1%merger_rate        = dm2%merger_rate 
    !
    ! masses
    dm1%M_vir              = dm2%M_vir
    dm1%M_tot              = dm2%M_tot
    dm1%M_tot              = dm2%M_tot
    dm1%M_acc              = dm2%M_acc
    dm1%M_halo             = dm2%M_halo
    !
    ! radius
    dm1%R_vir              = dm2%R_vir  
    dm1%R_200              = dm2%R_200
    dm1%R_halo             = dm2%R_halo
    !
    ! spin
    dm1%spin               = dm2%spin                  
    !
    ! halo parameters
    dm1%r_core             = dm2%r_core   
    dm1%rho_core           = dm2%rho_core  
    dm1%V_vir              = dm2%V_vir 
    dm1%T_vir              = dm2%T_vir  
    dm1%t_dyn              = dm2%t_dyn      
    dm1%concentration      = dm2%concentration 
              
    return
  end subroutine dm_copy

  !*****************************************************************************************************************
  
  subroutine dm_set_properties(dm,ts,x,y,z,Vx,Vy,Vz,M_vir,M_tot,R_vir,spin,dmacc,Lx,Ly,Lz,Ek_Ep)

    ! LOAD INTO DM THE PROPERTIES DIRECTLY MEASURED FROM THE N-BOBY SIMULATION

    implicit none

    integer(kind=4),intent(in)  :: ts            ! timestep at which the halo has been identified
    
    real(kind=8),intent(in)     :: x,y,z         ! space coordinates of the dm halo (in the N-body simulation)
    real(kind=8),intent(in)     :: Vx,Vy,Vz      ! Velocity of the centre of mass of the dark-matter halo
    real(kind=8),intent(in)     :: M_vir         ! virial mass (mass in the smallest ellipsoid where the virial theorem is verified to 20% accuracy)
    real(kind=8),intent(in)     :: M_tot         ! total halo mass as define by the halo finder (i.e. FOF or AdaptaHOP) 
    real(kind=8),intent(in)     :: R_vir         ! radius of a sphere with the same volume as the virial ellipsoid
    real(kind=8),intent(in)     :: spin          ! spin parameter
    real(kind=8),intent(in)     :: dmacc         ! mass of the halo particles that are not part of any identified halo at any previous timestep
    real(kind=8),intent(in)     :: Lx,Ly,Lz      ! angular momentum
    real(kind=8),intent(in)     :: Ek_Ep         ! ratio of kinetic to potential energies

    type(dm_type),intent(inout) :: dm

    dm%pos(:)         = (/x*1.e3,y*1.e3,z*1.e3/) ! Position are given in Mpc, convert in [kpc]
    dm%vel(:)         = (/Vx,Vy,Vz/)             ! Velocity of the centre of mass of the dark-matter halo [km/s]
    dm%ts             = ts
    dm%M_vir          = M_vir                    ! in [10^11Msun]
    dm%M_tot          = M_tot                    ! in [10^11Msun]
    dm%R_vir          = R_vir*1.e3               ! the virial radius is given in Mpc, convert in [kpc]
    dm%spin           = spin
    dm%dmacc          = dmacc                    ! in [10^11Msun]
    dm%L              = (/Lx,Ly,Lz/)             ! 
    dm%Ek_Ep          = Ek_Ep

    return
  end subroutine dm_set_properties

  !*****************************************************************************************************************

  subroutine dm_set_M_acc(dm,M_acc)

    ! SET THE DM MASS ACCRETED BY THE HALO BY INTEGRETED OVER ALL PREVIOUS PROGENITORS

    implicit none

    real(kind=8),intent(in)     :: M_acc             ! accreted mass

    type(dm_type),intent(inout) :: dm                ! dark matter halo

    dm%M_acc = M_acc  ! in [10^11Msun]

    return
  end subroutine dm_set_M_acc

  !*****************************************************************************************************************

  subroutine dm_flag_bad_properties(dm)

    ! FLAG (dm%compute = .false.) a dm component in function of energy properties, spin or integrated mass

    implicit none

    type(dm_type),intent(inout) :: dm                ! dark matter halo

#ifdef CLEAN
! -------------------------------------------------
    if ((dm%Ek_Ep .gt. Ecut) .or. (dm%M_acc .eq. 0.d0) .or. (dm%spin .le. 0.d0)) then
       dm%compute = .false.
    end if
#else
    if (dm%M_acc .eq. 0.d0) then
       dm%compute = .false.
    end if
! -------------------------------------------------
#endif

    return
  end subroutine dm_flag_bad_properties

  !*****************************************************************************************************************
  
  subroutine dm_set_merger_mass(dm,merger_mass)
    
    ! SET THE DM MASS GAINED BY MERGER ALONG THE MAIN BRANCH
    
    implicit none
    
    real(kind=8),intent(in)     :: merger_mass       
    
    type(dm_type),intent(inout) :: dm                ! dark matter halo

    dm%merger_mass = merger_mass                     ! in [10^11Msun]

    return
  end subroutine dm_set_merger_mass

  !*****************************************************************************************************************

  subroutine dm_set_merger_rate(dm,rate)

    implicit none

    real(kind=8),intent(in)     :: rate    ! merger_rate
    
    type(dm_type),intent(inout) :: dm      ! dark matter halo

    dm%merger_rate = rate                  ! in [10^11Msun/Gyr]

    return  
  end subroutine dm_set_merger_rate

  !*****************************************************************************************************************

  subroutine dm_set_bad_branch_mass(dm,bad_branch_mass)

    implicit none

    real(kind=8),intent(in)     :: bad_branch_mass  ! the bad_branch_mass connected to this halo
   
    type(dm_type),intent(inout) :: dm               ! dark matter halo

    dm%bad_branch_mass = bad_branch_mass            ! in [10^11Msun]

    return
  end subroutine dm_set_bad_branch_mass

  !*****************************************************************************************************************

  subroutine dm_set_acc_rate_left(dm,acc_rate_left)

    ! SET THE ACCRETION RATE OF THE DARK MATTER CONTENT USED BETWEEN THE MERGING TIME AND THE NEXT TIMESTEP
   
    implicit none

    
    real(kind=8),intent(in)     :: acc_rate_left     ! accretion rate of the halo
    
    type(dm_type),intent(inout) :: dm                ! dark matter halo

    dm%acc_rate_left = acc_rate_left                 ! in [10^11Msun/Gyr]

    return  
  end subroutine dm_set_acc_rate_left

  !*****************************************************************************************************************

  subroutine dm_set_acc_rate_right(dm,acc_rate_right)

    ! SET THE ACCRETION RATE OF THE DARK MATTER CONTENT USED BETWEEN THE LAST TIMESTEP AND THE NEXT MERGING TIME
    !   - estimate with the accretion rate of the descendent and the ratio of the accretion rate of this halo to the sum 
    !   - of all accretion rate of the progenitor of the descendent halos

    implicit none

    real(kind=8),intent(in)     :: acc_rate_right    ! accretion rate of the halo 
    
    type(dm_type),intent(inout) :: dm                ! dark matter halo

    dm%acc_rate_right = acc_rate_right               ! in [10^11Msun/Gyr]

    return  
  end subroutine dm_set_acc_rate_right
  
  !*****************************************************************************************************************
  
  subroutine dm_set_age_form(dm,age_univ)
  
    implicit none
    
    real(kind=8),intent(in)     :: age_univ  ! accretion rate of the halo 
    
    type(dm_type),intent(inout) :: dm        ! dark matter halo
    
    dm%age_form = age_univ
    
    return
  end subroutine
  
  !*****************************************************************************************************************
  
  subroutine dm_set_age_merge(dm,age_univ)
  
    implicit none
    
    real(kind=8),intent(in)     :: age_univ ! accretion rate of the halo 
    
    type(dm_type),intent(inout) :: dm       ! dark matter halo
    
    dm%age_merge = age_univ
    
    return
  end subroutine dm_set_age_merge

  !*****************************************************************************************************************

  subroutine dm_compute_derived_properties(dm,z)

    ! COMPUTE PROPERTIES OF THE HALO THAT CAN'T BE DIRECTLY MEASURED IN THE N-BODY SIMULATION
    !      - Compute R_200, V_vir, T_vir and t_dyn from the integrated mass using the halo density from the fitting formulae in Bryan & Norman 1997
    !      - Core radius and density depend on the concentration, which is computed with the model in dm_concentration

    implicit none

    real(kind=8),intent(in)     :: z                 ! redshift
    real(kind=8)                :: Omega_z,x,Delta   ! variables in the Bryan & Norman 1997 fitting formulae
    real(kind=8)                :: rho_200           ! 10^11 MSun/kpc^3

    type(dm_type),intent(inout) :: dm                ! dark matter halo

    if (dm%M_acc .le. 0.d0) then 
      call IO_print_error_message('halo with M_acc <= 0', &
          only_rank = rank, called_by = 'dm_compute_derived_properties')
      stop
    end if

    ! Compute rho_200 with the fitting formulae in Bryan & Norman (1998) and use it to compute R_200
    Omega_z  = Omega_m*(1.d0+z)**3./(Omega_m*(1.d0+z)**3 + 1.d0-Omega_m)   ! for a flat universe
    x        = Omega_z - 1.d0
    Delta    = 18.d0*pi**2+82.d0*x-39.0d0*x**2                             ! for a flat universe              
    rho_200  = Omega_m*(1.d0+z)**3./Omega_z*Delta*rho_crit_code_unit       ! [10^11Msun/kpc^3]     
    if (dm%M_tot .le. 0.d0) then 
      call IO_print_error_message('halo with M_tot <= 0', &
          only_rank = rank, called_by = 'dm_compute_derived_properties') 
          write(*,*) '> M_tot = ', dm%M_tot
      stop
    end if 
    dm%R_200 = (3.d0*dm%M_tot/(4.d0*pi*rho_200))**(1./3.)                  ! in code unit [kpc]
    if ((dm%M_tot .gt. 0.d0) .and. (dm%R_200 .le. 0.d0)) then
      call IO_print_error_message('Non-zero mass halo with zero radius', &
          only_rank = rank, called_by = 'dm_compute_derived_properties')
      stop
    end if

    ! set M_halo
    dm%M_halo = dm%M_acc
    ! set R_halo
    dm%R_halo = dm%R_200
    ! The concentration is computed with the model specified in the dm_concentration function    
    dm%concentration = dm_concentration(dm%M_halo,z)                              ! without unit
    ! Its value is used to determine the dm core radius
    dm%r_core = dm%R_halo/dm%concentration                                        ! in code unit [kpc]
    ! and core density (assuming that the dm halo is described by an NFW profile)     
    dm%rho_core = dm%M_halo/(4.d0*pi*dm%r_core**3.*dm_f_nfw(dm%concentration))    ! in code unit [10^11 Msun / kpc^3]
    ! V_vir is defined as the orbitral velocity at R_vir
    dm%V_vir = dm_velocity(dm%R_halo,dm) 
    ! Compute T_vir from V_vir
    dm%T_vir = max(1.d4,3.59d1*(dm%V_vir*vel_code_unit_2_kmPerSec)**2.) ! in code unit [K] 
    ! Compute the dynamical time of the dark matter strucure
    dm%t_dyn = (dm%R_halo/dm%V_vir)                          ! in code unit [Gyr]
    if (dm%t_dyn .le. 0.d0) then 
      call IO_print_error_message('halo with t_dyn <= 0', &
          only_rank = rank, called_by = 'dm_compute_derived_properties') 
      stop
    end if
    
    return
  end subroutine dm_compute_derived_properties
  
  !*****************************************************************************************************************

  subroutine dm_update_life_time(dm,dt)

    ! UPDATE (INCREASE) THE TIME LIFE OF A DARK-MATTER STRUCTURE

    implicit none

    real(kind=8),intent(in)      :: dt   ! time step

    type(dm_type),intent(inout)  :: dm   ! the dark matter component

    dm%life_time = dm%life_time + dt

    return
  end subroutine dm_update_life_time

  !*****************************************************************************************************************
  !
  ! FUNCTIONS
  !
  !*****************************************************************************************************************

  function dm_density_profile(x,rho_c)

     implicit none

     real(kind=8), intent(in)     :: x             ! r/r_0
     real(kind=8), intent(in)     :: rho_c         ! density at the core radius of the baryonic hot phase
     real(kind=8)                 :: dm_density_profile

     dm_density_profile = rho_c/(x*(1.d0 + x)**2.)

     return
  end function dm_density_profile

  !*****************************************************************************************************************

  function dm_f_nfw(x)

    ! dm_f_nfw(x) is a function such that 4*pi*rho_core*r_core**3*dm_f_nfw(x) is the dm mass enclosed within r = x * r_core  

    implicit none

    real(kind=8),intent(in) :: x
    real(kind=8)            :: dm_f_nfw

    dm_f_nfw = log(1.0d0 + x) - x/(1.0d0 + x) ! = integral of x/(1+x)**2 from zero to x
                                              ! = -d/dx ln(1+x)/x 

    return
  end function dm_f_nfw

  !*****************************************************************************************************************

  function dm_concentration(M_halo,z)

    ! RETURN THE DM CONCETRATION OF A GIVEN DARK-MATTER SUTRUCTURE
    ! The dm halo concentration is computed with the fitting formulae in Munoz-Cuartas et al 2010
    ! These formulae, established at 0 < z < 2, imply that, for given M_halo, the concentration has a minimum at z = z_turn
    ! This upturn may not be physical (but see Klypin et al 2010)
    ! We follow the conservative approach to use the Munoz-Cuartas et al fitting formulae while setting
    ! concentration(M_halo,z) = concentration(M_halo,z_turn) at z > z_turn

    implicit none

    real(kind=8),intent(in) :: M_halo
    real(kind=8),intent(in) :: z
    real(kind=8)            :: dm_concentration
    real(kind=8),parameter  :: w     = 0.029d0         ! parameters and variables in the Munoz-Cuartas et al fitting formulae
    real(kind=8),parameter  :: m     = 0.097d0
    real(kind=8),parameter  :: alpha = -110.001d0
    real(kind=8),parameter  :: beta  = 2469.720d0
    real(kind=8),parameter  :: gamma = 16.885d0
    real(kind=8)            :: MM,p,q,Delta,a,b,z_turn

    MM = log10(M_halo * 1.d11 * h_0)                   ! 1.0d11 transforms code units into Solar masses
    
    p  = - alpha/w/MM
    q  = - 2.0d0 * beta/w/MM

    Delta = 4.d0/27.d0*p**3+q**2

    ! z_turn is computed by requiring that the derivative of the concentration with respect to z is equal to zero
    ! This requirement is a cubic equation that has one real solution because alpha is negative
    ! We find this solution with the formula of Tartaglia-Cardan

    if ( -q-sqrt(Delta) .ge. 0.0d0) then ! "if" required to make sure that compiler treats cubic root of negative number correctly
       z_turn = (0.5d0*(-q+sqrt(Delta)))**(1./3.) + (0.5d0*(-q-sqrt(Delta)))**(1./3.) - gamma
    else
       z_turn = (0.5d0*(-q+sqrt(Delta)))**(1./3.) - (0.5d0*( q+sqrt(Delta)))**(1./3.) - gamma
    end if

    if (z < z_turn ) then
       a = w * z - m
       b = alpha/(z+gamma)+beta/(z+gamma)**2
    else
       a = w * z_turn - m
       b = alpha/(z_turn+gamma)+beta/(z_turn+gamma)**2
    end if

    dm_concentration = 10.**(a*MM+b)

    return 
  end function dm_concentration

  !*****************************************************************************************************************

  function dm_mass(dm,r)

    ! RETURN THE MASS ENCLOSE IN A SPHERE OF RADIUS r (FOR NFW PROFILE)

    implicit none

    real(kind=8),intent(in),optional :: r          ! radius
    real(kind=8)                     :: dm_mass    ! dark matter mass enclose in radius r, in code unit

    type(dm_type),intent(in)         :: dm         ! dm component

    if (present(r)) then
       dm_mass = 4.d0*pi*dm%rho_core*dm%r_core**3*dm_f_nfw(r/dm%r_core)  ! in code unit
    else
       dm_mass = dm%M_halo                                               ! in code unit
    end if

    return
  end function dm_mass

  !*****************************************************************************************************************

  function dm_potential(r,dm)

    ! RETURN THE DARK MATER HALO GRAVITATIONNAL POTENTIAL (NFW profile)

    implicit none
    
    real(kind=8),intent(in)    :: r             ! radius
    real(kind=8)               :: dm_potential  ! the dark matter gravitationnal potential at radius r (in code unit)
    real(kind=8)               :: x

    type(dm_type),intent(in)   :: dm            ! the dark matter component

    if (r .le. 0.d0) then
      call IO_print_error_message('r <= 0', &
                only_rank = rank, called_by = 'dm_potential')
      stop
    end if

    x = r/dm%r_core
     
    dm_potential = -4.d0*pi*gravconst_code_unit*dm%rho_core*dm%r_core**2*log(1.d0 + x)/x    ! in code unit
    
    return
  end function dm_potential

  !*****************************************************************************************************************

  function dm_velocity(r,dm) ! INTERFACE VERSION

    ! RETURN THE VELOCITY OF THE DARK MATTER PROFILE AT A RADIUS r
    
    implicit none

    real(kind=8), intent(in)     :: r               ! radius
    real(kind=8)                 :: dm_velocity     ! rotational velocity of the dark matter structure at radius r, in code unit

    type(dm_type), intent(in)    :: dm              ! dm component

    if (r .le. 0.d0) then
      call IO_print_error_message('r <= 0', &
              only_rank = rank, called_by = 'dm_velocity')
       stop
    end if
    
    dm_velocity = dm_velocity_(r,(/dm%rho_core,dm%r_core/))

    return
  end function dm_velocity

  !*****************************************************************************************************************

  function dm_velocity_(r,param) ! CORE VERSION

    ! RETURN THE VELOCITY OF THE DARK MATTER PROFILE AT A RADIUS r
    
    implicit none

    real(kind=8), intent(in)  :: r            ! radius
    real(kind=8),intent(in)   :: param(2)     ! parameter array : param(1) = dm%rho_core
                                              !                   param(2) = dm%r_core
    real(kind=8)              :: dm_velocity_ ! rotational velocity of the dark matter structure at radius r, in code unit
    real(kind=8)              :: x            ! addimentionnal radius

    x = r/param(2)

    dm_velocity_ = sqrt(4.d0*pi*gravconst_code_unit*param(1)*param(2)**2*(dm_f_nfw(x)/x)) ! in code unit

    return
  end function dm_velocity_

  !*****************************************************************************************************************

  function dm_escape_velocity(r,dm)

    ! RETURN THE ESCAPE VELOCITY OF THE DARK MATTER AT A RADIUS r
    ! The maximal value of the escape velocity is affected around r~0.1*r_core (it is a very good approximation)
    
    implicit none

    real(kind=8), intent(in)     :: r                    ! radius
    real(kind=8)                 :: dm_escape_velocity   ! escape velocity of the dark matter structure in code unit

    type(dm_type), intent(in)    :: dm                   ! dm component

    if (r .le. 0.d0) then
      call IO_print_error_message('r <= 0', &
                only_rank = rank, called_by = 'dm_escape_velocity')
      stop
    end if

    dm_escape_velocity = sqrt(-2.d0*dm_potential(r,dm))  ! in code unit

    return
  end function dm_escape_velocity
  
  !*****************************************************************************************************************
  
  function dm_orbital_radius(dm1,dm2)
  
    ! RETURN THE DISTANCE BETWEEN TWO dm HALOS
    implicit none
    
    type(dm_type), intent(in)    :: dm1               ! a dm component
    type(dm_type), intent(in)    :: dm2               ! an other dm component
    
    real(kind=8)                 :: dm_orbital_radius ! the distance between the two dark matter halos [kpc]
    
    dm_orbital_radius = sqrt(sum((dm1%pos-dm2%pos)**2.))

    return
  end function dm_orbital_radius
  
    !*****************************************************************************************************************
  
  function dm_dynamical_friction_time(dm1,dm2)
  
    ! RETURN THE DYNAMICAL FRICTION TIME ASSOCIATED TO A COUPLE OF TWO DM HALOS
    ! We consider here only dark-matter component 
    ! the galaxy mass is neglected
    
    implicit none
    
    real(kind=8)              :: dm_dynamical_friction_time
    real(kind=8)              :: M_vir                         ! The virial mass of the main (more massive) dark matter halo
    real(kind=8)              :: V_vir                         ! The virial velovity of the main halo
    real(kind=8)              :: r0_sat                        ! the initial distance between the two halos
    real(kind=8)              :: m_sat                         ! the baryonic mass of the satellite (the low-massive halo)
    real(kind=8)              :: Cl                            ! the Coulomb logarithm
    
    type(dm_type),intent(in)  :: dm1                           ! a given dark_matter structure (the main halo)
    type(dm_type),intent(in)  :: dm2                           ! an other dark_matter structure (the satellite)
    
    ! compute r_sat
    r0_sat = dm_orbital_radius(dm1,dm2)
    
    ! select M_vir, V_vir and m_sat
    if (dm_mass(dm1) .gt. dm_mass(dm2)) then
        ! dm1 is the main structure
        M_vir = dm_mass(dm1)
        V_vir = dm1%V_vir
        m_sat = dm_mass(dm2)
    else
        ! h2 is the main structure
        M_vir = dm_mass(dm2)
        V_vir = dm2%V_vir
        m_sat = dm_mass(dm1)
    end if
    ! compute Cl
    Cl = log(1.d0 + M_vir/m_sat)
    
    dm_dynamical_friction_time = (V_vir*r0_sat**2.)/(2.d0*2.428*gravconst_code_unit*m_sat*Cl)
    
    return
   end function dm_dynamical_friction_time

  !*****************************************************************************************************************
  !
  ! PRINTING PROCEDURES
  !
  !*****************************************************************************************************************

  subroutine dm_load_dm_data(dm,dm_data)

    ! CREATE THE DM OUTPUT PROPERTIES LIST

    implicit none

    real(kind=8), intent(inout) :: dm_data(nb_dm_field)

    type(dm_type),intent(in)    :: dm    ! dm component

    ! For information
    ! 'x                     ','y                     ','z                     ',&
    ! 'Vx                    ','Vy                    ','Vz                    ',&
    ! 'M_vir                 ','M_tot                 ','M_acc                 ',&
    ! 'R_vir                 ','R_200                 ','spin                  ',&
    ! 'acc_rate_left         ','acc_rate_right        ','V_vir                 ',&
    ! 'T_vir                 ','t_dyn                 ','concentration         ',&
    ! 'r_core                ','dm_rho_core           ','age_form              '

    dm_data = (/dm%pos(1)*1.e-3, dm%pos(2)*1.e-3, dm%pos(3)*1.e-3, dm%vel(1), dm%vel(2), dm%vel(3),&
                1.d11*dm%M_vir, 1.d11*dm%M_tot, 1.d11*dm%M_acc, &
                dm%R_vir, dm%R_200, dm%spin, 1.d2*dm%acc_rate_left, 1.d2*dm%acc_rate_right, &
                dm%V_vir, dm%T_vir, dm%t_dyn, dm%concentration, dm%r_core, dm%rho_core, dm%age_form/)

    return
  end subroutine dm_load_dm_data

  !*****************************************************************************************************************

  subroutine dm_print(unit,form,dm)

    ! PRINT DM-PROPERTIES IN eGALICS TMP OUTPUT FILE
     
    implicit none

    integer(kind=4),intent(in)      :: unit  ! file unit
    integer(kind=4)                 :: status,i,hdutype

    character(*)                    :: form  ! fits or tmp_bin 
    
    real(kind=8)                    :: dm_data(nb_dm_field)

    type(dm_type),intent(in)        :: dm    ! dm component

    call dm_load_dm_data(dm,dm_data)

    select case (trim(form))
      case ('tmp_bin')
        write(unit) dm_data  ! directly write data in the tmp binary output file 
      case ('fits')
        ! move to dm extension
        call ftmahd(unit,hdu_dm,hdutype,status) 
        if (status .gt. 0) then
          call IO_print_error_message('ftmahd status', &
                only_rank = rank, called_by = 'dm_print')
          stop ! stop the program
        end if
        ! init
        call ftirow(unit,0,1,status) 
        if (status .gt. 0) then
          call IO_print_error_message('ftirow status', &
                only_rank = rank, called_by = 'dm_print')
          stop ! stop the program
        end if
        ! write data in the dm entension
        call ftpcld(unit,1,1,1,1,dm%age_form+dm%life_time,status)
        do i=2, nb_dm_field+1   
          call ftpcld(unit,i,1,1,1,dm_data(i-1),status)  
          if (status .gt. 0) then
            call IO_print_error_message('ftpcld status', &
                only_rank = rank, called_by = 'dm_print')
            stop ! stop the program
          end if
        end do
      case default
        call IO_print_error_message('Unknwon output data format', &
                only_rank = rank, called_by = 'dm_print')
        stop ! stop the program
    end select 
    
    return
  end subroutine dm_print

  !*****************************************************************************************************************

end module dm
