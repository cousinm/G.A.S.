module gas

  use utilities         ! Acces to integrator Rombint, Bessel function, ....

  public

  !*****************************************************************************************************************
  !
  !  OVERVIEW
  !
  !  The gas module defines all properties and procedures assossiated to a gas object.
  !  In G.A.S. a gas phase is defined by its total mass, metals mass
  !     and also the mass of six main ISM elements (H, He, C12, N14, O16 and Fe56)
  !  This module starts with the definition of this gas type
  !
  !  SUBROUTINES IN THIS MODULE
  !
  !  gas_read_gas_properties           : Read gas properties (Metallicity bins, Initial abundances etc ...)
  !
  !  gas_read_timescales               : Read the inertial gas cascade timescale tables
  !
  !  gas_set_reference_mass            : Initialize gas reference mass : M_gas_min
  !
  !  gas_igm_initialize                : Initialize the IGM phase
  !      called by : G.A.S. main program
  !
  !  gas_void                          : Set to null values all gas properties
  !
  !  gas_void_array                    : Set to null values all cells of a gas array
  !
  !  gas_void_gsh                      : Set to null values, all cells of of fields of a gas structuration history tab
  !
  !  gas_deallocate                    : Deallocate a gas structure (g%elts)
  !
  !  gas_deallocate_array              : Deallocate all gas cells of a gas array
  !
  !  gas_deallocate_gsh                : Deallocate all fields of a gas structuration history tab
  !
  !  gas_copy                          : Copy a gas component in an other one
  !
  !  gas_copy_array                    : Copy an array of gas component in an other one
  !
  !  gas_copy_gsh                      : Copy a gas structuration history tab in an other one
  !
  !  gas_set_component                 : set value of a given gas component
  !
  !  gas_add                           : Add a gas component to an other one
  !
  !  gas_add_igm_gas                   : Add gas to the IGM phase
  !
  !  gas_add_array                     : Add an array of gas components to an other one
  !
  !  gas_sub                           : Substract a gas component to an other one
  !
  !  gas_sub_array                     : Substract an array of gas components to an other one
  !
  !  gas_inject_termal_energy          : Inject thermal internal energy into a given gas phase, and set/update the Temperature
  !
  !  gas_evap                          : Compute the gas evaporation according to the temperature and an escape velocity
  !
  !  gas_finalize                      : Deallocate all data array used by the gas module
  !

  !
  !  FUNCTIONS IN THIS MODULE
  !
  !  gas_igm_Temperature               : Return the IGM temperature at a given redshift
  !
  !  gas_signature                     : Return a normalized gas component
  !                                       (allow to compute gas transfer rates based on a gas component or element mass fractions)
  !
  !  gas_mass                          : return the total gas mass of a gas component (or sub component)
  !
  !  gas_metallicity                   : return the metal mass fraction of a given gas component (Z = mZ/mass)
  !
  !  gas_temp                          : return the gas temperature
  !
  !  gas_dt_optim                      : return the optimal evolution timestep associated to a given gas component
  !
  !  O_H                               : return the Oxygen/Hydrogen mass ratio
  !
  !  O_Fe                              : return the Oxygen/Iron mass ratio
  !
  !  gas_timescales                    : return a specific timescale of the inertial cascade
  !
  !  Bonnor_Ebert_Mass                 : return the Bonnor Ebert mass of the current largest scale of the disc
  !
  !  Saturation_Mass                   : return the saturation mass of the inertial cascade at a given scale k
  !
  !  Maxwell_Boltzman_Vdist            : return F(v) for a Maxwell boltzmann velocity distribution
  !
  !  Maxwell_Boltzman_Vdist_shifted    : return F(v) for a shifted Maxwell velocity boltzmann distribution
  !
  !  Maxwell_Boltzman_Edist            : return F(v) for a Maxwell boltzmann kinetic energy distribution
  !
  !  Maxwell_Boltzman_Edist_shifted    : return F(v) for a shifted Maxwell boltzmann kinetic energy distribution
  !
  !  Maxwell_Boltzman_V_to_T           : return the average temperature T for a Maxwell boltzmann gas with mean velocity V
  !
  !  Maxwell_Boltzman_Etot             : return the internal energy (per perticle) in a Maxwell boltzmann gas with mean temperature T
  !
  !  Maxwell_Boltzman_Vp               : return the most probable velocity in a Maxwell boltzmann gas with mean temperature T
  !
  !  INTERFACE OPERATORS IN THIS MODULE
  !
  !  gas_add_                          : Add a gas component to an other one  --> operator +
  !
  !  gas_add_array_                    : Add an array of gas components to an other one  --> operator +
  !
  !  gas_sub_                          : Substract a gas component to an other one  --> operator -
  !
  !  gas_sub_array_                    : Substract an array of gas components to an other one  --> operator -
  !
  !  gas_scalar_multiply_left          : return gas * a --> operateur *
  !
  !  gas_arr_scalar_multiply_left      : return gas(:) * a --> operateur *
  !
  !  gas_scalar_multiply_right         : return a * gas --> operateur *
  !
  !  gas_arr_scalar_multiply_right     : return a * gas(:) --> operateur *
  !
  !*****************************************************************************************************************

  ! DEFINITION OF GLOBAL VARIABLES LINKED TO THE STARS STRUCTURE *******************

  ! GAS_TYPE DEFINITION *****************************

  type gas_type
     real(kind=8)              :: Temp               ! temperature of the gas
     real(kind=8)              :: mass               ! total mass  (in 10^11 Msun)
     real(kind=8)              :: mZ                 ! mass of metals
     real(kind=8)              :: f_str              ! structured fraction (in mass)
     real(kind=8), allocatable :: elts(:)            ! mass of main ISM elements (e.g H1, He, C12, N14, O16, Fe56)
  end type gas_type

  ! GAS STRUCTURATION HISTORY (GSH_TYPE) DEFINITION *

  type gsh_type
    real(kind=8) :: t_cascade                    ! inertial cascade evolution clock
    real(kind=8) :: t_form                       ! inertial cascade formation timescale
    real(kind=8) :: t_eq                         ! inertial cascade equilibrium timescale
    real(kind=8) :: t_emp                        ! inertial cascade emptying timescale
    real(kind=8) :: t_prod                       ! NON-sfg to SFG production timescale
    !
    type(gas_type), allocatable  :: gas(:)       ! gas mass at a given structuration level
    type(gas_type), allocatable  :: in_rate(:)   ! input mass rate for a given structuration level
    type(gas_type), allocatable  :: out_rate(:)  ! output mass rate for a given structuration level
  end type gsh_type

  ! DEFINE IGM *********************

  type(gas_type)                 :: igm              ! intergalactic medium

  ! DEFINE SOME PARAMETER VALUES ********************

  integer(kind=4)                :: nElts            ! number of Main ISM elements followed
  integer(kind=4)                :: nMetBins         ! nb of metallicity bins in the stellar evolution (gas) model(s)
  integer(kind=4)                :: nAccRateBins     ! nb accretion rate bins use in the turbulent cascade model(s)
  integer(kind=4)                :: nDiscScaleBins   ! nb of disc scale height use in the turbulent cascade model(s)

  character(4),allocatable       :: Elt_names(:)     ! list of main ISM elements followed (e.g H1, He, C12, N14, O16, Fe56)
  character(4),allocatable       :: gas_component(:) ! list of gas component (total mass, metals or elements)

  real(kind=8)                   :: M_gas_null       ! minimal gas mass
  real(kind=8)                   :: M_gas_crit       ! minimal mass to take into account ejection process
  real(kind=8)                   :: M_gas_min        ! the minimal gas mass takes into account in adaptive time step evolution process
                                                     ! = fb * dm_particle_mass
                                                     ! set during initialization (IO_read_parameter_file)
  real(kind=8)                   :: dlAccRate        ! Accretion rate bin size (in log)
  real(kind=8)                   :: dlDiscScale      ! Disc scale height bin size (in log)
  real(kind=8)                   :: AccRatemin       ! minimum value of accretion rate available for the turbulent cascade model
  real(kind=8)                   :: DiscScalemin     ! minimum value for disc scale available for the turbulent cascade model
  real(kind=8),allocatable       :: MetBins(:)       ! mass fraction corresponding to metallicity steps
  real(kind=8),allocatable       :: AccRateBins(:)   ! List of Accretion rates available for turbulent cascade models
  real(kind=8),allocatable       :: DiscScaleBins(:) ! List of disc scale height available for turbulence cascade models
  real(kind=8),allocatable       :: t_str(:,:)       ! gas structuration timescale (first dim: acc rate, second dim: disc scale)
  real(kind=8)                   :: t_str_max        ! maximum value of the gas structuration timescale
  real(kind=8),allocatable       :: t_form(:,:)      ! cascade formation timescale (first dim: acc rate, second dim: disc scale)
  real(kind=8)                   :: t_form_max       ! maximum value of the gas cascade formation timescale
  real(kind=8),allocatable       :: t_eq(:,:)        ! cascade equilibrium timescale (first dim: acc rate, second dim: disc scale)
  real(kind=8)                   :: t_eq_max         ! maximum value of the gas cascade equilibrium timescale  

  type(gas_type),allocatable     :: InitAbund(:)     ! Initial abundances (mass fraction) associated to each metallicity bin
                                                     ! InitAbund is defined as a gas object therefore InitAbund[Z]%mass = 1., InitAbund[Z]%mZ = metBins(Z)
                                                     ! and for each main ISM elements InitAbund[Z]%elts(e) = X0_elt
  real(kind=8),parameter         :: mu                       = 0.62*mp                                  ! mean particle mass for a fully ionised plasma (kg)
  real(kind=8),parameter         :: diffuse_gas_temp         = 8.d3                                     ! temperature of warm diffuse gas from cold-streams and cooling flows
  real(kind=8),parameter         :: r_perfect_gas            = kb/mu*Gyr_in_s**2./(kpc_in_m**2.)        ! specific gas constant [in code unit] (kpc^2/Gyr^2/K)
  real(kind=8),parameter         :: adiab_ind                = 1.4d0
  real(kind=8),parameter         :: large_scale_min_surf_den = 1.0d1/mass_code_unit_in_M_Sun*(1.e3)**2. ! star formation mass surface density threshold
                                                                                                        ! 10 Msun/pc converted in code unit 10^11Msun/kpc**2.
  real(kind=4)                   :: larson_mu_slope          ! column density power law slope
  real(kind=4)                   :: larson_sig_slope         ! velocity dispersion power law slope
  real(kind=8)                   :: l_star                   ! star forming filament size 0.1pc
  real(kind=8)                   :: k_star                   ! associated wavelenght
  real(kind=8)                   :: sig_star                 ! velocity dispersion threshold
  real(kind=8)                   :: mu_star                  ! column density threshold
                                                             ! 150 Msun/pc converted in code unit 10^11Msun/kpc**2.
#ifdef POLYTROPIC
  real(kind=8),parameter         :: gamma                    = 6.d0/5.d0                                ! polytropic index (=1.2)
  real(kind=8),parameter         :: kappa                    = 1.                                       ! without unit FIXED HERE
#endif
  real(kind=8),parameter         :: r_sf                     = 1.d-1/1.e3                               ! 0.1 pc converted in kpc
!
#ifdef NOSFG_2SFG
! -------------------------------------------------
  real(kind=8),parameter         :: input_sfg_frac           = 5.d-3 !1.d-2
! -------------------------------------------------
#endif
! NOSFG_2SFG

  ! INTERFACE OPERATOR DECLARATIONS ****************

   interface assignment (=)    ! allows to copy a gas component by using the symbol '='
     module procedure gas_copy, gas_copy_array, gas_copy_gsh
  end interface assignment (=)

  interface operator (+)    ! allows to sum two gas components by using the symbol '+'
     module procedure gas_add_, gas_add_array_
  end interface operator (+)

  interface operator (-)    ! allows to substract two gas components by using the symbol '-'
     module procedure gas_sub_, gas_sub_array_
  end interface operator (-)

  interface operator (*)    ! allows to ponderate a gas component (a* or *a) by using the symbol '*'
     module procedure gas_scalar_multiply_left,gas_scalar_multiply_right,gas_arr_scalar_multiply_left,gas_arr_scalar_multiply_right
  end interface operator (*)

contains

  !*****************************************************************************************************************
  !
  ! SUBROUTINES
  !
  !*****************************************************************************************************************

  subroutine gas_read_gas_properties

    implicit none

    integer(kind=4)          :: j,e ! loop indexes (metallicity and elements)

    character(MAXPATHSIZE)   :: filename
    character(MAXPATHSIZE)   :: line
    character(MAXPATHSIZE)   :: message

    call IO_print_message('gas_read_gas_properties')

    if (main_process .or. physical_process) then
       write(filename,'(a,a,a,a)') trim(input_path), '/gas_properties.in'
       write(message,'(a,a,a,a)') 'Load data from : ', 'gas_properties.in'
       call IO_print_message(message)
       ! Open the "gas properties" file
       open(unit = gasprop_unit, file = filename, status = 'old')

       do
          read(gasprop_unit, '(a)', end = 2) line
          if (trim(line) .eq. 'START') then
             read(gasprop_unit,*) nElts, nMetBins ! the number main ISM elements followed, the number of Metallicity bin used
             call IO_print_message('use',param_name=(/'nElts                    ','nMetBins                 '/), &
                  int_param_val=(/nElts,nMetBins/))
             ! create the table of the Main ISM element names
             allocate(Elt_names(nElts))
             allocate(gas_component(nElts +2))
             ! create the metallicity bin table
             allocate(MetBins(nMetBins))
             ! create the initial abundance table
             allocate(InitAbund(nMetBins))
             ! read Main ISM elements list
             read(gasprop_unit,*) (Elt_names(e), e=1,nElts)
             ! create a table which allow to run through all gas component (total, metals, and elements)
             gas_component = (/'mass','mZ  ',Elt_names/)
             call IO_print_message('Chemical evolution is followed for:')
             write(message,'(a)') trim(Elt_names(1))
             do e = 2, nElts
                write(message,'(a,a,a)') trim(message), ', ', trim(Elt_names(e))
             end do
             call IO_print_message(message)
             ! read Metallicity bins and initial abundance table
             ! InitAbund is defined as a gas object therefore InitAbund[Z]%mass = 1., InitAbund[Z]%mZ = metBins(Z)
             ! and for each main ISM elements InitAbund[Z]%elts(e) = X0_elt
             do j = 1, nMetBins
                call gas_void(InitAbund(j)) ! init the gas object
                InitAbund(j)%mass = 1.d0
                read(gasprop_unit,*) MetBins(j), (InitAbund(j)%elts(e), e=1,nElts)
                InitAbund(j)%mZ   = MetBins(j)
             end do
             if ((main_process) .and. (.not. physical_process)) then
                if (allocated(InitAbund)) deallocate(InitAbund)
             end if
             exit  ! quit do loop
          end if
          if (line(1:1) .eq. '#') then
             cycle ! header or something like this (skip)
          else
             call IO_print_error_message('Impossible to read the line',only_rank=rank,called_by='gas_read_gas_properties')
             stop ! stop the program
          end if
       end do
2      close(gasprop_unit)
    end if

    return
  end subroutine gas_read_gas_properties

  !*****************************************************************************************************************

  subroutine gas_read_timescales

    ! LOAD GAS INERTIAL CASCADE TIMESCALES

     integer(kind=4)          :: m,l ! loop indexes (accretion rate and disc scaleheight)

    character(MAXPATHSIZE)   :: filename
    character(MAXPATHSIZE)   :: line
    character(MAXPATHSIZE)   :: message

    call IO_print_message('gas_read_timescales')

    write(filename,'(a,a,a,a)') trim(input_path), '/gas_timescales.in'
    write(message,'(a,a,a,a)') 'Load data from : ', 'gas_timescales.in'
    call IO_print_message(message)

    if (main_process .or. physical_process) then
       ! Open the "gas properties" file
       open(unit = gasprop_unit, file = filename, status = 'old')
       !
       do
          read(gasprop_unit, '(a)', end = 2) line
          if (trim(line) .eq. 'START') then
             read(gasprop_unit,*) larson_mu_slope, larson_sig_slope ! read power law indexes
             read(gasprop_unit,*) l_star, sig_star, mu_star         ! read constant parameters
             ! convertion unit
             ! convert from pc to kpc
             l_star = l_star/1.e3
             k_star = 1./l_star
             ! convert from km/s in kpc/Gyr
             sig_star = sig_star/vel_code_unit_2_kmPerSec
             ! convert from Msun/pc^2 to 10^11 Msun/kpc**2
             mu_star = mu_star*M_Sun_in_mass_code_unit*(1.e3)**2
             !
             read(gasprop_unit,*) nAccRateBins, nDiscScaleBins      ! read number of bins
             call IO_print_message('use',param_name=(/'nAccRateBins             ','nDiscScaleBins           '/), &
                  int_param_val=(/nAccRateBins,nDiscScaleBins/))
             if (physical_process) then
                 read(gasprop_unit,*) dlAccRate, dlDiscScale ! read bin size (in log)
                 ! allocate arrays
                 allocate(AccRateBins(nAccRateBins))
                 allocate(DiscScaleBins(nDiscScaleBins))
                 allocate(t_form(nAccRateBins,nDiscScaleBins))
                 allocate(t_eq(nAccRateBins,nDiscScaleBins))
                 allocate(t_str(nAccRateBins,nDiscScaleBins))
                 ! load accretion rates and disc scales
                 read(gasprop_unit,*) AccRateBins
                 read(gasprop_unit,*) DiscScaleBins
                 ! convertion unit
                 ! convert from in Msun/yr to 10^11Msun/Gyr
                 AccRateBins = AccRateBins-log10(mass_rate_code_unit_2_MsunPerYr)
                 ! save minimum values
                 AccRatemin   = minval(AccRateBins)
                 DiscScalemin = minval(DiscScaleBins)
                 !
                 read(gasprop_unit,*) ! a blanck line
                 ! read formation timescale (in log)
                 do m = 1, nAccRateBins
                    read(gasprop_unit,*) (t_form(m,l), l=1,nDiscScaleBins)
                 end do
                 ! set t_form_max
                 t_form_max = 10.**maxval(t_form) ! in Gyr
                 !
                 read(gasprop_unit,*) ! a blanck line
                 ! read equilibrium timescale (in log)
                 do m = 1, nAccRateBins
                    read(gasprop_unit,*) (t_eq(m,l), l=1,nDiscScaleBins)
                 end do
                 ! set t_form_max
                 t_eq_max = 10.**maxval(t_eq) ! in Gyr
                 !
                 read(gasprop_unit,*) ! a blanck line
                 ! read structuration/emptying time (in log)
                 do m = 1, nAccRateBins
                    read(gasprop_unit,*) (t_str(m,l), l=1,nDiscScaleBins)
                 end do
                 ! set t_str_max
                 t_str_max = 10.**maxval(t_str) ! in Gyr
              end if
              exit  ! quit do loop
          end if
          if (line(1:1) .eq. '#') then
            cycle ! header or something like this (skip)
          else
            call IO_print_error_message('Impossible to read the line',only_rank=rank,called_by='gas_read_timescales')
            write(*,*) trim(line)
            stop ! stop the program
          end if
       end do
2      close(gasprop_unit)
    end if

    return
  end subroutine gas_read_timescales
  
  !*****************************************************************************************************************

  subroutine gas_set_reference_mass

    ! INITIALIZE M_gas_min A MASS PARAMETERS LINKED TO THE DARK-MATTER N-BODY SIMULATION

    implicit none

    M_gas_min  = 5.d-1 * baryon_fraction * dm_particle_mass
    M_gas_crit = 1.d-2 * M_gas_min
    M_gas_null = 1.d-11

    call IO_print_message('use',param_name=(/'M_gas_min [Msun]         ','M_gas_crit [Msun]        ','M_gas_null [Msun]        '/), &
             real_param_val=(/M_gas_min*mass_code_unit_in_M_Sun,M_gas_crit*mass_code_unit_in_M_Sun,M_gas_null*mass_code_unit_in_M_Sun/))

    return
  end subroutine gas_set_reference_mass

  !*****************************************************************************************************************

  subroutine gas_igm_initialize

    ! INIT IGM PHASE

    implicit none

    if (physical_process) call gas_void(igm)

    return
  end subroutine gas_igm_initialize

  !*****************************************************************************************************************

  subroutine gas_void(gas)

    ! SET TO NULL VALUE ALL PROPERTIES OF A GAS COMPONENT

    implicit none

    type(gas_type),intent(inout) :: gas  ! the gas component

    gas%Temp  = -1.d0  ! not already setled
    gas%mass  = 0.d0   ! no mass
    gas%mZ    = 0.d0   ! no mass
    gas%f_str = 0.d0   ! fraction of structured gas
    if (.not. allocated(gas%elts)) then
       allocate(gas%elts(nElts))
    end if
    gas%elts(:) = 0.d0

    return
  end subroutine gas_void

  !***********************************************************************************************************

  subroutine gas_void_array(g)

    ! SET TO NULL ALL CELLS OF A GAS ARRAY
    ! WARNING, in case of allocated array, g must be allocated before the call of this subroutine

    implicit none

    integer(kind=4)                          :: N     ! number of element in the gas array
    integer(kind=4)                          :: i,i_l ! loop index
                                                      ! a gsh_tab is indexed from 0 to n_end = nlevel+1 and therefore contains ncells = nlevel+2
    type(gas_type),intent(inout),allocatable :: g(:)

    if (allocated(g)) then
        N = size(g)
        i_l = lbound(g,dim=1)
        if (i_l .eq. 0) N = N -1
        if (N .gt. 0) then
            do i = i_l, N
                call gas_void(g(i))
            end do
        end if
    end if

    return
  end subroutine gas_void_array

  !***********************************************************************************************************

  subroutine gas_void_gsh(gsh)

    ! SET TO NULL ALL CELLS OF ALL FIELDS OF A GAS STRUCTURATION HISTORY TAB

    implicit none

    type(gsh_type), intent(inout)       :: gsh

    ! if the gas structuration history is already allocated, we have to erase it
    if (allocated(gsh%gas)) then
        call gas_deallocate_gsh(gsh)
    endif

    ! Init timescales
    gsh%t_cascade = 0.d0
    gsh%t_eq      = t_eq_max
    gsh%t_form    = 9.9d-1*t_eq_max
    gsh%t_emp     = t_str_max
    gsh%t_prod    = 9.9d-1*t_str_max
    
   
    ! By default the gas structuration history contain the diffuse, the structured/fragmented gas cells
    ! 1 : unstructured / diffuse gas
    ! 2 : structured / fragmented gas 

    if (.not. allocated(gsh%gas)) allocate(gsh%gas(2))            ! create
    call gas_void_array(gsh%gas)                                  ! and init
    if (.not. allocated(gsh%in_rate)) allocate(gsh%in_rate(2))    ! create
    call gas_void_array(gsh%in_rate)                              ! and init
    if (.not. allocated(gsh%out_rate)) allocate(gsh%out_rate(2))  ! create
    call gas_void_array(gsh%out_rate)                             ! and init

    return
  end subroutine gas_void_gsh

  !***********************************************************************************************************

  subroutine gas_deallocate(g)

    ! DEALLOCATE A GAS STRUCTURE

    implicit none

    type(gas_type),intent(inout) :: g  ! the gas component

    if (allocated(g%elts)) deallocate(g%elts)

    return
  end subroutine gas_deallocate

  !***********************************************************************************************************

  subroutine gas_deallocate_array(g)

    ! DEALLOCATE ALL GAS STRUCTURE OF THE ARRAY
    ! WARNING, in case of allocated array, g must be allocated before the call of this subroutine

    implicit none

    integer(kind=4)                          :: N     ! number of element in the gas array
    integer(kind=4)                          :: i,i_l ! loop index
                                                      ! a gsh_tab is indexed from 0 to n_end = nlevel+1 and therefore contains ncells = nlevel+2
    type(gas_type),intent(inout),allocatable :: g(:)

    if (allocated(g)) then
        N = size(g)
        i_l = lbound(g,dim=1)
        if (i_l .eq. 0) N = N -1
        if (N .gt. 0) then
            do i = i_l, N
                call gas_deallocate(g(i))
            end do
        end if
        deallocate(g)
    end if

    return
  end subroutine gas_deallocate_array

  !***********************************************************************************************************

  subroutine gas_deallocate_gsh(gsh)

    ! DEALLOCATE ALL GAS STRUCTURE CELLS IN ALL FIELDS OF A GAS STRUCTURATION HISTORY TAB
    ! WARNING, in case of allocated array, gsh must be allocated before the call of this subroutine

    implicit none

    type(gsh_type),intent(inout) :: gsh

!~     call gas_deallocate_array(gsh%gas)
    call gas_deallocate_array(gsh%in_rate)
    call gas_deallocate_array(gsh%out_rate)

    return
  end subroutine gas_deallocate_gsh

  !***********************************************************************************************************

  subroutine gas_copy(g1,g2)

    ! COPY g2 INTO g1

    implicit none

    type(gas_type),intent(inout) :: g1
    type(gas_type),intent(in)    :: g2

    call gas_void(g1)  ! init

    g1%Temp  = g2%Temp  ! copy the temperature
    g1%mass  = g2%mass  ! copy the total mass
    g1%mZ    = g2%mZ    ! copy the metal mass
    g1%f_str = g2%f_str ! structured fraction

    if (allocated(g2%elts)) then
       if (.not. allocated(g1%elts)) then
          allocate(g1%elts(nElts))  ! create the array
       endif
       g1%elts = g2%elts            ! copy
    else
       call IO_print_error_message('Try to copy an non-existent ISM gas component',only_rank=rank,called_by='gas_copy')
       call IO_print_message('with',only_rank=rank,component='gas', &
            param_name=(/'g2%mass                  ','g2%mZ                    '/),real_param_val=(/g2%mass,g2%mZ/))
       stop ! stop the program
    end if

    return
  end subroutine gas_copy

  !***********************************************************************************************************

  subroutine gas_copy_array(g1,g2)

    ! COPY THE ARRAY OF GAS COMPONENT g2 INTO g1
    ! Warning, in case of allocated array, g2 must be allocated before the call of this subroutine

    implicit none

    integer(kind=4)              :: N1,N2         ! nb elements in the arrays
    integer(kind=4)              :: i,i_l1,i_l2   ! loop index

    type(gas_type),intent(inout),allocatable :: g1(:)
    type(gas_type),intent(in),allocatable    :: g2(:)

    N1 = size(g1)
    N2 = size(g2)
    i_l1 = lbound(g1,dim=1)
    i_l2 = lbound(g2,dim=1)
    if (i_l1 .eq. 0) N1 = N1 -1
    if (i_l2 .eq. 0) N2 = N2 -1

    if (N2 .gt. 0) then
       !
       if (.not. allocated(g1)) then
          allocate(g1(i_l2:N2))
          call gas_void_array(g1)
       else
          if (i_l1 .ne. i_l2) then
             call IO_print_error_message('Try to copy a gas array into an other with a different starting indexes',only_rank=rank,called_by='gas_copy_array')
             call IO_print_message('with',only_rank=rank,component='gas', &
                  param_name=(/'i_l1                     ','i_l2                     '/), int_param_val=(/i_l1,i_l2/))
             stop ! stop the program
          end if
          if (N1 .ne. N2) then
             call IO_print_error_message('Try to copy a gas array into an other with a different size',only_rank=rank,called_by='gas_copy_array')
             call IO_print_message('with',only_rank=rank,component='gas', &
                  param_name=(/'N1                       ','N2                       '/), int_param_val=(/N1,N2/))
             stop ! stop the program
          end if
       end if
       do i = i_l2, N2
          call gas_copy(g1(i),g2(i))
       end do
    else
       call IO_print_error_message('Try to copy a gas arrays without any element',only_rank=rank,called_by='gas_copy_array')
       call IO_print_message('with',only_rank=rank,component='gas', &
            param_name=(/'N2                       '/),int_param_val=(/N2/))
       stop ! stop the program
    end if

    return
  end subroutine gas_copy_array

  ! ***********************************************************************************************************

  subroutine gas_copy_gsh(gsh1,gsh2)

    ! COPY A GAS STRUCTURATION HISTORY TAB IN AN OTHER ONE
    ! WARNING, all fields of gsh2 must be allocated before the call of this subroutine

    implicit none

    type(gsh_type), intent(inout)  :: gsh1  ! a gas structuration history table
    type(gsh_type), intent(in)     :: gsh2  ! an other gas structuration history table

    ! Copy timescales
    gsh1%t_cascade = gsh2%t_cascade  ! Inertial cascade clock
    gsh1%t_form = gsh2%t_form ! Inertial cascade formation timescale
    gsh1%t_eq = gsh2%t_eq ! Inertial cascade equilibrium timescale
    gsh1%t_prod = gsh2%t_prod  ! Inertial cascade production timescale
    gsh1%t_emp = gsh2%t_emp ! Inertial cascade emptying timescale
    
    call gas_copy_array(gsh1%gas,gsh2%gas)
    call gas_copy_array(gsh1%in_rate,gsh2%in_rate)
    call gas_copy_array(gsh1%out_rate,gsh2%out_rate)
    
    return

  end subroutine gas_copy_gsh

  !***********************************************************************************************************

  subroutine gas_set_component(g,v,component)

    ! SET A GAS COMPONENT
    ! (mass, mZ, mH, .... ) to a given value: v

     implicit none

     integer(kind=4)              :: e

     real(kind=8),intent(in)      :: v          ! the value

     character(*),intent(in)      :: component  ! allow to select the gas component

     type(gas_type),intent(inout) :: g          ! the gas component

     select case (trim(component))
     case('Temp','temp','T')
        g%Temp = v
     case ('mass')
        g%mass = v
     case ('mZ')
        g%mZ = v
     case ('f_str')
        g%f_str = v
     case default
        do e = 1, nElts
          if (trim(component) .eq. trim(Elt_names(e))) then
             g%elts(e) = v
             exit
          end if
        end do
     end select

     return
  end subroutine gas_set_component

  !***********************************************************************************************************

  subroutine gas_add(g1,g2,Vwind,Vesc,m_esc,called_by)

    ! STANDARD GAS VERSION
    ! ADD g2 TO THE GAS COMPONENT g1 (g1 = g1 + g2)
    ! Waring: In case of allocated array g1 and g2 must be allocated before the call of the subroutine

    implicit none

    character(*),intent(in),optional    :: called_by    ! name of the function which has called this function

    real(kind=8),intent(in),optional    :: Vwind        ! Macroscopic velocity of g2 in the g1 referential
    real(kind=8),intent(in),optional    :: Vesc         ! Escape velocity of the g1 host structure
    real(kind=8)                        :: E1,Ew        ! Thermal energy
    real(kind=8)                        :: m1,mw        ! Masses
    real(kind=8)                        :: T1,Tw        ! Temeratures
    real(kind=8)                        :: f_in         ! fraction of the wind gas (g2) that is kept into the gas g1
    real(kind=8)                        :: error        ! integration error

    type(gas_type),intent(inout)        :: g1           ! A gas component
    type(gas_type),intent(in)           :: g2           ! An other gas component
    type(gas_type)                      :: g
    type(gas_type),intent(out),optional :: m_esc        ! The escape mass

    ! create a local copy of g2
    call gas_void(g)
    call gas_copy(g,g2)

    if (present(Vwind)) then
        ! init m_esc
        call gas_void(m_esc)
        !
        ! we compute the fraction of particles which have velocity (in probability distribution) larger than the escape velocity
        ! the velocity probability distribution is computed like a bimodal Maxwell-Boltzmann distribution
        ! the first one dedicated to the pre-existing hot gas at temperature T1
        ! the second one dedicated to the wind phase at temperature T2 and shifted of Vwind
        !
        m1 = gas_mass(g1)
        mw = gas_mass(g)
        T1 = gas_temp(g1)
        Tw = gas_temp(g)
        E1 = 0.d0
        if (T1 .gt. 0.d0) E1 = Maxwell_Boltzman_T_to_Eint(T1) ! per mass unit
        !
        if (Tw .gt. 0.d0) then
            ! fraction of gas in the wind which is kept into the halo
            f_in = min(1.d0,max(0.d0,Ronbint(Maxwell_Boltzman_Vdist_shifted,0.d0,Vesc,(/Tw,Vwind/),called_by='gas_add / f_in')))
            !if (f_in .lt. 1.d-2) f_in = 0.d0
            !if (f_in .gt. 9.9d-1) f_in = 1.d0
            !
            if ((f_in .gt. 0.d0) .and. (f_in .le. 1.d0)) then
               ! A fraction of the wind is catched by the hot atmosphere
               ! set escape gas
               m_esc = (1.d0-f_in)*g
               ! update g, the wind phase
               g = f_in*g
               ! Termal energy in the wind phase
               ! We take into account the termal energy of the gas that is catched by the potential well
               Ew = Ronbint(Maxwell_Boltzman_Edist_shifted,0.d0,Vesc,(/Tw,Vwind/),error=error,called_by='gas_add / Internal energy') ! per mass unit
               ! To compute new internal energy of the hot gas phase
               ! We have to add kinetic energy of wind catched particles
               E1 = (E1*m1 + (Ew + 5.d-1*Vwind**2.)*mw*f_in)/(m1 + f_in*mw) ! per mass unit
               ! compute new equilibrium temperature in [K] of the hot gas phase
               T1 = max(diffuse_gas_temp,Maxwell_Boltzman_Eint_to_T(E1))
               !
               if (f_in .lt. 1.d0) then
                    ! Energy in the wind (E1) = total thermal energy in the wind at temperature Tw - total thermal energy in the gas which stay in the halo (Ew; computed previously)
                    E1 = Maxwell_Boltzman_T_to_Eint(Tw)*mw - Ew*mw*f_in ! total thermal energy
                    ! compute equilibrium temperature in [K] of the wind phase
                    Tw = max(diffuse_gas_temp,Maxwell_Boltzman_Eint_to_T(E1/(mw*(1.d0-f_in))))
                    ! set to new temperature
                    if (Tw .le. 0.d0) then
                        call IO_print_error_message('Negative tempertaure',only_rank=rank,called_by='gas_add / Wind temperaure')
                        call IO_print_message('use',only_rank=rank,component='gas', &
                            param_name = (/'Tw                       ','E1                       ','Ew                       ','error                    '/), &
                            real_param_val  = (/Tw,E1,Ew,error/))
                        stop ! stop the program
                    end if
                    call gas_set_component(m_esc,Tw,component='Temp')
                end if
            else
                ! All the wind escape the hot atmosphere
                ! set escape gas
                call gas_copy(m_esc,g)
                ! reset g
                call gas_void(g)
                ! the equilibrium temperature of the gas phase doesn't change
                T1 = gas_temp(g1)
            endif
            if ((m1 .gt. 0.d0) .and. (T1 .le. 0.d0)) then
                call IO_print_error_message('Negative temperature',only_rank=rank,called_by='gas_add / hot phase temperature')
                call IO_print_message('use',only_rank=rank,component='gas', &
                        param_name = (/'T1                       '/),real_param_val  = (/T1/))
                stop ! stop the program
            end if
        else
            call IO_print_error_message('Negative temperature',only_rank=rank,called_by='gas_add / wind temperature')
            call IO_print_message('use',only_rank=rank,component='gas', &
                        param_name = (/'Tw                       '/),real_param_val  = (/Tw/))
            stop ! stop the program
        end if
    else
        if ((g%mass .gt. 0.d0) .and. (g%Temp .gt. 0.d0)) then
            ! g have a mass and a temperature
            if ((g1%mass .gt. 0.d0) .and. (g1%Temp .gt. 0.d0)) then
                ! g1 have a mass and a temperature
                T1 = max(diffuse_gas_temp,(g1%mass*g1%Temp + g%mass*g%Temp)/(g1%mass + g%mass))
            else
                ! no mass in g1, set to g
                T1 = g%Temp
            end if
        else
            ! no mass in g
            ! keep T1 to initial value
            T1 = gas_temp(g1)
        end if
    end if
    !
    ! set to new value
    call gas_set_component(g1,T1,component='Temp')
    ! computed structured fraction
    if ((g1%mass + g%mass) .gt. 0.d0) g1%f_str = (g1%mass*g1%f_str + g%mass*g%f_str)/(g1%mass + g%mass)
    ! add mass
    g1%mass = g1%mass + g%mass
    ! add metal mass
    g1%mZ   = g1%mZ   + g%mZ
    ! sum chemical elements
    if (allocated(g1%elts)) then
       if (allocated(g%elts)) then
          g1%elts = g1%elts + g%elts  ! sum for all main ISM elements
       else
          call IO_print_error_message('Try to add a non-existent ISM gas component',only_rank=rank,called_by=called_by)
          stop ! stop the program
       end if
    else
       if (allocated(g%elts)) then
          call IO_print_error_message('Try to add mass to a non-existent ISM gas component',only_rank=rank,called_by=called_by)
          stop ! stop the program
       end if
    end if

    return
  end subroutine gas_add

  !***********************************************************************************************************

  subroutine gas_add_igm_gas(gas,called_by)

    ! Add gas to the igm

    implicit none

    character(*),intent(in),optional     :: called_by      ! name of the function which has called this function
    character(MAXPATHSIZE)               :: message

    type(gas_type),intent(in)            :: gas

    if (present(called_by)) then
        write(message,'(a,a)') 'gas_add_igm_gas/', trim(called_by)
    else
        message = 'gas_add_igm_gas'
    end if

    call gas_add(igm,gas,called_by=trim(message))

    return
  end subroutine gas_add_igm_gas

  !***********************************************************************************************************

  subroutine gas_add_array(g1,g2)

    ! ADD g2 TO THE ARRAY OF GAS COMPONENTS g1 (g1 = g1 + g2)
    ! Waring: In case of allocated array g1 and g2 must be allocated before the call of the subroutine

    implicit none

    integer(kind=4)              :: N1,N2         ! nb elements in the arrays
    integer(kind=4)              :: i,i_l1,i_l2   ! loop index

    type(gas_type),intent(inout),allocatable :: g1(:)
    type(gas_type),intent(in),allocatable    :: g2(:)

    N1 = size(g1)
    N2 = size(g2)
    i_l1 = lbound(g1,dim=1)
    i_l2 = lbound(g2,dim=1)
    if (i_l1 .eq. 0) N1 = N1 -1
    if (i_l2 .eq. 0) N2 = N2 -1

    if (N2 .gt. 0) then
       !
       if (i_l1 .ne. i_l2) then
          call IO_print_error_message('Try to copy a gas arrays into an other with a different starting indexes',only_rank=rank,called_by='gas_copy_array')
          call IO_print_message('with',only_rank=rank,component='gas', &
                param_name=(/'i_l1                     ','i_l2                     '/), int_param_val=(/i_l1,i_l2/))
          stop ! stop the program
       end if
       if (N1 .ne. N2) then
          call IO_print_error_message('Try to sum two gas arrays with different sizes',only_rank=rank,called_by='gas_copy_array')
          call IO_print_message('with',only_rank=rank,component='gas', &
               param_name=(/'N1                       ','N2                       '/), int_param_val=(/N1,N2/))
          stop ! stop the program
       end if
       do i = i_l2, N2
          call gas_add(g1(i),g2(i))
       end do
    end if

    return
  end subroutine gas_add_array

  !***********************************************************************************************************

  subroutine gas_sub(g1,g2,therm_mode,called_by)

    ! SUBSTRACT g2 TO THE GAS COMPONENT g1 (g1 = g1 - g2)

    implicit none

    integer(kind=4)                  :: e          ! loop index

    character(*),intent(in),optional :: called_by  ! name of the function which has called this function
    character(*),intent(in),optional :: therm_mode ! thermodynamic mode: - isothermal or cooling mode
    character(MAXPATHSIZE)           :: message    ! a message

    real(kind=8)                     :: m1,m2      ! mass of the gas components
    real(kind=8)                     :: T          ! A temperature

    type(gas_type),intent(inout)     :: g1         ! a gas component
    type(gas_type),intent(in)        :: g2         ! an other gas component

    m1 = gas_mass(g1); m2 = gas_mass(g2)
    if ((m2 .gt. m1) .and. (abs(m2 - m1) .gt. num_accuracy)) then
       call IO_print_error_message('Substract too much gas',only_rank=rank,called_by=called_by)
       call IO_print_message('with',only_rank=rank,component='gas', &
            param_name=(/'m1                       ','m2                       '/),real_param_val=(/m1,m2/))
       stop ! stop the program
    end if

    if ((g2%mZ .gt. g1%mZ) .and. (abs(g2%mZ - g1%mZ) .gt. num_accuracy)) then
       call IO_print_error_message('Substract too much metals',only_rank=rank,called_by=called_by)
       call IO_print_message('with',only_rank=rank,component='gas', &
            param_name=(/'g1%mZ                    ','g2%mZ                    ','g1%mass                  ','g2%mass                  '/), &
            real_param_val=(/g1%mZ,g2%mZ,g1%mass,g2%mass/))
       call IO_print_message('with',only_rank=rank,component='gas', &
            param_name=(/'g1%Z                     ','g2%Z                     '/), &
            real_param_val=(/g1%mZ/g1%mass,g2%mZ/g2%mass/))
       stop ! stop the program
    end if

    if (allocated(g1%elts)) then
       if (allocated(g2%elts)) then
          do e = 1, nElts
             if ((g1%elts(e) < g2%elts(e)) .and. (abs(g2%elts(e) - g1%elts(e)) .gt. num_accuracy)) then
                write(message,'(a,a)') 'Substract too much  ', Elt_names(e)
                call IO_print_error_message(message,only_rank=rank,called_by=called_by)
                call IO_print_message('with',only_rank=rank,component='gas', &
                     param_name=(/'g1%elts(e)               ','g2%elts(e)               '/), &
                     real_param_val=(/g1%elts(e),g2%elts(e)/))
                write(errunit,*) 'g1%elts(e) = ', g1%elts(e), ', g2%elts(e) = ', g2%elts(e)
                write(errunit,*) 'g1%mass    = ', g1%mass,    ', g2%mass    = ', g2%mass
                stop ! stop the program
             end if
          end do
       end if
    else
       if (allocated(g2%elts)) then
          call IO_print_error_message('Substract mass to a non-existent ISM gas component',only_rank=rank,called_by=called_by)
          stop ! stop the program
       end if
    end if

    ! @ this point all checks have been done, we can substract g2 to g1

    if (present(therm_mode)) then
        select case(trim(therm_mode))
        case('isothermal','iso')
            ! the original and the substracted gas come from the same reservoir @ T1
            if (abs(g1%Temp - g2%Temp) .gt. num_precision) then
                write(message,'(a)') 'No consistency between temperatures T1 and T2 [isotermal mode]'
                call IO_print_warning_message(message,only_rank=rank,called_by='gas_sub')
                call IO_print_message('used',only_rank=rank,component='bh', &
                        param_name = (/'T1                       ','T2                       ','error                    '/), &
                        real_param_val  = (/g1%Temp,g2%Temp,abs(g1%Temp - g2%Temp)/))
            end if
            T = max(diffuse_gas_temp,gas_temp(g1))
        case('cooler','cooling','cool')
            ! the cool gas is setled to T_cool_gas, the temperature of the hot atmosphere is not modified
            if ((g1%Temp .lt. g2%Temp) .and. (abs(g1%Temp - g2%Temp) .gt. num_precision)) then
                write(message,'(a)') 'No consistency between temperatures T1 and T2 [cooling mode]'
                call IO_print_warning_message(message,only_rank=rank,called_by='gas_sub')
                call IO_print_message('used',only_rank=rank,component='bh', &
                        param_name = (/'T1                       ','T2                       ','error                    '/), &
                        real_param_val  = (/g1%Temp,g2%Temp,abs(g1%Temp - g2%Temp)/))
            end if
            T = max(diffuse_gas_temp,gas_temp(g1))
        case default
            write(message,'(a,a,a)') 'Keyword ', trim(therm_mode), ' not defined'
            call IO_print_error_message(message,only_rank=rank,called_by='gas_sub')
            stop  ! stop the program
        end select
    end if

    ! substract mass
    g1%mass = g1%mass - g2%mass
    if (g1%mass .le. M_gas_null) then
        call gas_void(g1)
    else
        ! substract metal mass
        g1%mZ = g1%mZ - g2%mZ
        if (g1%mZ .lt. M_gas_null) then
            g1 = g1%mass*InitAbund(1)
        else
            g1%elts = g1%elts - g2%elts
            if (minval(g1%elts) .lt. M_gas_null) g1 = g1%mass*InitAbund(1)
        end if
        !
        ! set new gas temperature
        call gas_set_component(g1,T,component='Temp')
    end if

    return
  end subroutine gas_sub

  !***********************************************************************************************************

  subroutine gas_sub_array(g1,g2)

    ! SUBSTRACT THE ARRAY OF GAS COMPONENT g2 TO g1 (g1 = g1 - g2)

    integer(kind=4)              :: N1,N2            ! nb elements in the arrays
    integer(kind=4)              :: i,i_l1,i_l2      ! loop index

    type(gas_type),intent(inout),allocatable :: g1(:)
    type(gas_type),intent(in),allocatable    :: g2(:)

    N1 = size(g1)
    N2 = size(g2)
    i_l1 = lbound(g1,dim=1)
    i_l2 = lbound(g2,dim=1)
    if (i_l1 .eq. 0) N1 = N1 -1
    if (i_l2 .eq. 0) N2 = N2 -1

    if (N2 .gt. 0) then
       !
       if (i_l1 .ne. i_l2) then
          call IO_print_error_message('Try to copy a gas arrays into an other with a different starting indexes',only_rank=rank,called_by='gas_copy_array')
          call IO_print_message('with',only_rank=rank,component='gas', &
                param_name=(/'i_l1                     ','i_l2                     '/), int_param_val=(/i_l1,i_l2/))
          stop ! stop the program
       end if
       if (N1 .ne. N2) then
          call IO_print_error_message('Try to substract two gas arrays with different sizes',only_rank=rank,called_by='gas_copy_array')
          call IO_print_message('with',only_rank=rank,component='gas', &
               param_name=(/'N1                       ','N2                       '/),int_param_val=(/N1,N2/))
          stop ! stop the program
       end if
       do i = il_2, N2
          call gas_sub(g1(i),g2(i))
       end do
    end if

    return
  end subroutine gas_sub_array

  !*****************************************************************************************************************

  subroutine gas_inject_termal_energy(g,Qtherm)

    ! INJECT A THEMAL ENERGY (Qtherm) TO THE GAS (SET/UPDATE THE TEMPERATURE)

    implicit none

    real(kind=8),intent(in)      :: Qtherm    ! thermal energy [code unit]
    real(kind=8)                 :: m,T,E     ! local value, mass, temperature, thermal energy

    type(gas_type),intent(inout) :: g         ! the gas component

    if (Qtherm .le. 0.d0) return              ! no heating
    if (gas_mass(g) .le. 0.d0) return         ! no gas

    m = gas_mass(g)
    T = gas_temp(g)
    E = 0.d0
    if (T .gt. 0) then
        E = Maxwell_Boltzman_T_to_Eint(T)*m ! total thermal energy
    end if
    ! add energy
    E = (E + Qtherm)/m ! per mass unit
    ! compute new equilibrium temperature
    T = Maxwell_Boltzman_Eint_to_T(E)
    call gas_set_component(g,T,component='Temp')

    if (is_NaN(g%Temp)) then
        call IO_print_error_message('T is NaN ',only_rank=rank,called_by='gas_inject_termal_energy')
        stop ! stop the program
    end if

    return
  end subroutine gas_inject_termal_energy

  !*****************************************************************************************************************

  subroutine gas_evap(g,m_evap,Vesc)

    ! APPLY GAS EVAPORATION ONTO A GAS COMPONENT

    implicit none

#ifdef PRINT_WARNING
! -----------------------------------------------
    character(MAXPATHSIZE)       :: message ! a message to display
! -------------------------------------------------
#endif

    real(kind=8),intent(in)      :: Vesc    ! escape velocity of the host gas structure [code unit]
    real(kind=8)                 :: f_in    ! real escape fraction
    real(kind=8)                 :: E1,E2
    real(kind=8)                 :: m,T

    type(gas_type),intent(inout) :: g       ! the gas component
    type(gas_type),intent(out)   :: m_evap  ! the evaporated gas

    call gas_void(m_evap) ! init

#ifndef NO_EVAP
! -------------------------------------------------
    T = gas_temp(g)
    if (T .le. 0.d0) return ! not defined yet
    !
    f_in = min(1.d0,max(0.d0,Ronbint(Maxwell_Boltzman_Vdist,0.d0,Vesc,(/T/),called_by = 'gas_evap / f_in')))
    !
    ! CHECK
    if (f_in .lt. 5.d-1) then
        ! no evap
        return
    end if
    !
    if (f_in .lt. 9.9d-1) then
        !
        ! energies
        ! total mass
        m = gas_mass(g)
        ! total energy
        E1 = Maxwell_Boltzman_T_to_Eint(T)*m
        ! total energy in the enclosed gas
        E2 = Ronbint(Maxwell_Boltzman_Edist,0.d0,Vesc,(/T/),called_by='gas_evap / Internal energy')*m*f_in
        !
        ! compute new equilibrium of enclosed gas
        T = Maxwell_Boltzman_Eint_to_T(E2/(f_in*m))
        ! test new T
        if (T .lt. diffuse_gas_temp) then
            ! no evap
            if (f_in .lt. 8.d-1) then
#ifdef PRINT_WARNING
! -------------------------------------------------
                write(message,'(a,f5.1,a)') 'T_gas < T_cool_gas --> No evap; ', (1.d0-f_in)*1.d2, ' [%]'
                call IO_print_warning_message(message,only_rank=rank,called_by='gas_evap')
! -------------------------------------------------
#endif
            end if
            return
        end if
        call gas_set_component(g,T,component='Temp')
        !
        ! compute escape mass
        m_evap = g*(1.d0-f_in)
        ! compute new equilibrium temperature of ejecta
        T = max(diffuse_gas_temp,Maxwell_Boltzman_Eint_to_T((E1-E2)/((1.d0-f_in)*m)))
        call gas_set_component(m_evap,T,component='Temp')
        !
        ! substract escaped mass
        g = g*f_in
        !
    end if
! -------------------------------------------------
#endif
! NO_EVAP

    return
  end subroutine gas_evap

  !*****************************************************************************************************************

  subroutine gas_finalize

    ! DEALLOCATE ALL DATA ARRAY USED BY GAS MODULE

    implicit none

    if (allocated(Elt_names)) deallocate(Elt_names)
    if (allocated(MetBins))   deallocate(MetBins)
    if (allocated(InitAbund)) deallocate(InitAbund)

    return
  end subroutine gas_finalize

  !*****************************************************************************************************************
  !
  ! FUNCTIONS
  !
  !*****************************************************************************************************************

  function gas_igm_Temperature(z)

    implicit none

    real(kind=8),intent(in)  :: z                      ! current redshift
    real(kind=8)             :: gas_igm_Temperature    ! The igm temperature @ z
    real(kind=8)             :: z_tmp, T_tmp, dz

    if (z .gt. z_reion) then
        gas_igm_Temperature = 3.d2*((1.d0 + z)/1.d2)**2.
    else
        z_tmp = z_reion
        T_tmp = T_reion
        do while (z_tmp .gt. z)
            dz = min(0.01,abs(z_tmp-z))
            T_tmp = T_tmp - 2.d0*T_reion*(1.d0+z)*dz/(1.d0+z_reion)**2.
            z_tmp = z_tmp +dz
        end do
        gas_igm_Temperature = T_tmp
    end if

    return
  end function gas_igm_Temperature

  !*****************************************************************************************************************

  function gas_signature(gas,apply_as,called_by)

    ! RETURN A NORMALIZED GAS COMPONENT (THE FOOTPRINT OF A GAS COMPONENT)
    ! allow to compute gas transfer rate based on a gas component

    implicit none

    integer(kind=4)                   :: comp          ! loop under gas component

    character(MAXPATHSIZE)            :: message       ! a message to display
    character(*),intent(in),optional  :: called_by     ! name of the function which has called this function
    character(*),intent(in),optional  :: apply_as      ! the gas signature function can be use in different cases
                                                       ! when this function is used as an output rate builder, we have to check critical mass
    real(kind=8)                      :: M, Mcomp
    real(kind=8)                      :: signature

    type(gas_type),intent(in)         :: gas           ! a gas component
    type(gas_type)                    :: gas_signature ! a normalized gas component

    call gas_void(gas_signature)

    M = gas_mass(gas)

    if (M .lt. 0.d0) then
        if (present(called_by)) then
            !
            write(message,'(a,a)') 'gas_signature, called by ', trim(called_by)
        else
            !
            write(message,'(a)') 'gas_signature'
        endif
        call IO_print_error_message('gas%mass < 0.',only_rank=rank,called_by=trim(message))
        stop
    end if

    if (present(apply_as)) then
        select case (trim(apply_as))
        case ('rate_builder','rate')
            ! run throup gas components
            do comp = 1, nElts +2
                Mcomp = gas_mass(gas,component=trim(gas_component(comp)))
                if (Mcomp .lt. M_gas_crit) then
                    signature = 0.d0
                else
                    signature = slope_limiter(Mcomp,M_gas_crit)*Mcomp/M
                    if (signature .lt. 1.d-50) signature = 0.d0
                end if
                call gas_set_component(gas_signature,signature,component=trim(gas_component(comp)))
            end do
        case default
            write(message,'(a,a,a)') 'Keyword ', trim(apply_as), ' not defined'
            call IO_print_error_message(message,only_rank=rank,called_by='gas_signature')
            stop  ! stop the program
        end select
    else
      gas_signature = (1.d0/M)*gas
    end if

    return
  end function gas_signature

  !***********************************************************************************************************

  function gas_mass(gas,component)

    ! RETURN THE TOTAL MASS OF THE GAS COMPONENT
    ! The selection between all ISM elements are possible, the case 'Metals' are also takes into account

    implicit none

    integer(kind=4)                   :: e          ! loop index

    character(*),intent(in),optional  :: component  ! allow to select the bulge component (sfg, non-sfg, stars)
    character(MAXPATHSIZE)            :: message    ! a message to display

    logical                           :: found
    real(kind=8)                      :: gas_mass

    type(gas_type),intent(in)         :: gas        ! a gas component

    if (present(component)) then
        found = .false.
        ! select the specific component
        select case (trim(component))
        case ('all','gas','mass')
            gas_mass = gas%mass
            found = .true.
        case ('Metals','metals','mZ')
            gas_mass = gas%mZ
            found = .true.
        case default
        !
        do e = 1, nElts
           if (trim(component) .eq. trim(Elt_names(e))) then
              gas_mass = gas%elts(e)
              found = .true.
              exit
            end if
        end do
        !
        if (.not. found) then
            write(message,'(a,a,a)') 'Keyword ', trim(component), ' not defined'
            call IO_print_error_message(trim(message),only_rank=rank,called_by='gas_mass')
            stop  ! stop the program
        end if
        end select
    else
        ! if no specific component are given, return the total gas mass
        gas_mass = gas%mass
    end if

    return
  end function gas_mass

  !*****************************************************************************************************************

  function gas_metallicity(gas)

    ! RETURN THE METAL MASS FRACTION

    implicit none

    real(kind=8)              :: gas_metallicity

    type(gas_type),intent(in) :: gas            ! a gas component

    if (gas%mass .gt. 0.d0) then
       gas_metallicity = gas%mZ / gas%mass
    else
       gas_metallicity = 0.0d0
    end if

    return
  end function gas_metallicity

  !*****************************************************************************************************************

  function gas_temp(gas)

    ! RETURN THE TEMPERATURE OF THE GAS

    implicit none

    real(kind=8)              :: gas_temp

    type(gas_type),intent(in) :: gas            ! a gas component

    gas_temp = gas%Temp

    return
  end function gas_temp

  !***********************************************************************************************************

  function gas_dt_optim(gas,input,output,pp,evolving_rate,warning_up,warning_down,called_by)

    ! RETURN THE OPTIMAL EVOLUTION TIMESTEP ASSOCIATED TO A GIVEN GAS COMPONENT
    ! input and output rates are the global input and output rates associated to the gas component
    ! pp (optional) is the maximum possible variation in mass fraction
    ! If evolving rate is given in input :
    !       warning_up    is the maximum value of "evolving_rate" compatible with quasi-static evolution
    !       warning_down  is the minimum value of "evolving_rate" compatible with quasi-static evolution

    implicit none

    integer(kind=4)                       :: comp      ! loop index on gas component

    character(MAXPATHSIZE)                :: message   ! a message to display
    character(MAXPATHSIZE)                :: case_applied
    character(*),intent(in),optional      :: called_by ! name of the function which has called this function

    real(kind=8),intent(in),optional      :: pp
    real(kind=8)                          :: warning, dM
    real(kind=8)                          :: gas_dt_optim, gas_dt_optim_prev
    real(kind=8)                          :: in_rate, out_rate, M, rate
    real(kind=8)                          :: e_rate

    type(gas_type),intent(in)             :: gas       ! a gas component
    type(gas_type),intent(in)             :: input     ! the global input rate
    type(gas_type),intent(in)             :: output    ! the global output rate
    type(gas_type),intent(in),optional    :: evolving_rate
    type(gas_type),intent(out),optional   :: warning_up, warning_down

    gas_dt_optim      = 1.d0
    gas_dt_optim_prev = 1.d0

    if ((present(evolving_rate)) .and. (present(warning_up)) .and. (present(warning_down))) then
        call gas_void(warning_up)                              ! init
        call gas_void(warning_down)                            ! init
        warning_up   = (1.d0+physical_precision)*evolving_rate ! set
        warning_down = (1.d0-physical_precision)*evolving_rate ! set
    end if

    ! run throught components
    do comp = 1, nElts +2
       !
       M         = gas_mass(gas,component=trim(gas_component(comp)))     ! mass
       in_rate   = gas_mass(input,component=trim(gas_component(comp)))   ! input rate
       out_rate  = gas_mass(output,component=trim(gas_component(comp)))  ! output rate
       if (present(evolving_rate)) then
         e_rate = gas_mass(evolving_rate,component=trim(gas_component(comp)))
       end if
       !
       ! TEST in_rate and in_rate
       if ((in_rate .lt. 0.d0) .or. (out_rate .lt. 0.d0)) then
          write(message,'(a,a,a)') 'in_rate < 0.  or out_rate < 0. [component: ', trim(gas_component(comp)), ' ]'
          call IO_print_error_message(message,only_rank=rank,called_by=called_by)
#ifdef PRINT_WARNING
! -------------------------------------------------
          call IO_print_message('',only_rank=rank,component='gas', &
                    param_name = (/'M                        ','in_rate                  ','out_rate                 '/), &
                    real_param_val  = (/M,in_rate,out_rate/))
! -------------------------------------------------
#endif
! PRINT_WARNING
          stop
       end if
       !
       rate = in_rate - out_rate  ! global evolution rate
       !
       if (abs(rate) .gt. 0.d0) then
          !
          ! compute optimal time-step
          if (M .gt. M_gas_min) then
             !
             ! adaptative time-step
             !
             if (present(evolving_rate)) then
                 ! take account some evolution (increase or decrease of the evolving rate)
                 if (in_rate .gt. out_rate) then
                    ! INPUT DOMINATED EVOLUTION
                    ! the evolution of this component is leaded by input
                    ! take into account a greater value of input (than the previous time-step) to permit a growing in the next time-step
                    warning = gas_mass(warning_up,component=trim(gas_component(comp)))
                 else
                    ! OUTPUT DOMINATED EVOLUTION
                    ! the evolution of this component is leaded by output
                    ! take into account a lower value of input (than the previous time-step) to permit a decrease in the next time-step
                    warning = gas_mass(warning_down,component=trim(gas_component(comp)))
                 endif
                 ! recompute global evolution rate; i.e replace original evol_rate by warning
                 rate = (in_rate - e_rate + warning) - out_rate
             end if
             !
             if (present(pp)) then
                case_applied = 'adaptive ts [pp%]'
                dM = pp*M
             else
                if (comp .gt. 1) then
                    ! metal and chemical elements
                    dM = max(physical_precision,5.d-1)*M
                    case_applied = 'adaptive ts [50%]'
                else
                    ! total mass
                    dM = physical_precision*M
                    case_applied = 'adaptive ts [standard]'
                end if
             end if
             !
             gas_dt_optim = min(gas_dt_optim, dM/abs(rate))
             if (abs(in_rate)  .gt. 0.d0) gas_dt_optim = min(gas_dt_optim, 5.d-1*M/abs(in_rate))
             if (abs(out_rate) .gt. 0.d0) gas_dt_optim = min(gas_dt_optim, 5.d-1*M/abs(out_rate))
             !
             ! Update the evolving rate window
             if ((present(evolving_rate)) .and. (present(warning_up)) .and. (present(warning_down))) then
                 if (in_rate .gt. out_rate) then
                    ! INPUT DOMINATED EVOLUTION
                    ! warning_up has already been setled
                    ! set warning_down
                    rate = (-dM/gas_dt_optim) - (in_rate-warning) + out_rate
                    rate = max(rate,gas_mass(warning_down,component=trim(gas_component(comp))))
                    call gas_set_component(warning_down,max(0.d0,rate),component=trim(gas_component(comp)))
                 else
                    ! OUTPUT DOMINATED EVOLUTION
                    ! warning_down has already been setled
                    ! set warning_up
                    rate = (dM/gas_dt_optim) - (in_rate-warning) + out_rate
                    rate = min(rate,gas_mass(warning_up,component=trim(gas_component(comp))))
                    call gas_set_component(warning_up,max(0.d0,rate),component=trim(gas_component(comp)))
                 endif
             end if
          else
             !
             ! no adaptive time-step
             !
             if (present(evolving_rate)) then
                ! We take into account the worest case; i.e without e_rate
                ! In this condition, the variation will be lower than predicted (the real abs(rate) < rate)
                rate = (in_rate-e_rate) - out_rate
             end if
             !
             if (abs(rate) .gt. 0.d0) then
                 if (rate .gt. 0.d0) then
                    ! input rate dominated
                    case_applied = 'in rate dominated / reach 2xM_gas_min'
                    dM =  max(2.d0*M_gas_min - M, 3.d-1*M)
                 else
                    ! output rate dominated
                    ! void in the next time step
                    case_applied = 'out rate dominated / void'
                    dM = M - M_gas_null
                 end if
                 !
                 gas_dt_optim = min(gas_dt_optim, dM/abs(rate))
             end if
          end if ! adaptative timestep cases
          !
          ! CHECKS
          if ((present(warning_up)) .and. (present(warning_down))) then
             ! test if warning_down < warning_up
             if (gas_mass(warning_up,component=trim(gas_component(comp))) &
                    .lt. gas_mass(warning_down,component=trim(gas_component(comp)))) then
                write(message,'(a,a,a,a,a)') 'warning_up < warning_down [component: ', trim(gas_component(comp)), ' / ', trim(case_applied),' ]'
                call IO_print_error_message(message,only_rank=rank,called_by='gas_dt_optim')
                stop
             end if
          end if
       end if
    end do

    gas_dt_optim = 9.5d-1*gas_dt_optim

    return
  end function gas_dt_optim

  !***********************************************************************************************************

  function O_H(gas)

    ! Return the Oxygen/Hydrogen Ratio, in number N_O / N_H

    implicit none

    real(kind=8)              :: O_H
    real(kind=8)              :: M_O, M_H

    type(gas_type),intent(in) :: gas

    O_H = -1.0  ! init

    if (gas_mass(gas) .le. 0.d0) return  ! no gas

    M_H = gas_mass(gas,component='H1')
    if (M_H .le. 0.d0) return            ! no H

    M_O  = gas_mass(gas,component='O16')
    if (M_O .le. M_gas_null) return      ! no O

    O_H = M_O / M_H * 1.008d0/15.999d0   ! from mass ratio to number ratio

    return
  end function O_H

  !***********************************************************************************************************

  function O_Fe(gas)

    ! Return the Oxygen/Iron Ratio, in number N_O / N_Fe

    implicit none

    real(kind=8)              :: O_Fe
    real(kind=8)              :: M_O, M_Fe

    type(gas_type),intent(in) :: gas

    O_Fe = -1.0  ! init

    if (gas_mass(gas) .le. 0.d0) return  ! no gas

    M_Fe = gas_mass(gas,component='Fe56')

    if (M_Fe .gt. 0.d0) then
        M_O  = gas_mass(gas,component='O16')
        O_Fe = M_O / M_Fe * 1.008d0/5.585d1  ! from mass ratio to number ratio
    end if

    return
  end function O_Fe

  !*****************************************************************************************************************

   function gas_timescales(acc_rate,h,timescale)

    ! Return one of specific timescale of the inertial cascade 
    ! as a function of the gas accretion rate and the disc scale height
    ! The output is given in Gyr

    implicit none

    character(*),intent(in)     :: timescale       ! selection of the timescale
    character(MAXPATHSIZE)      :: message         ! a message to display
    
    real(kind=8),intent(in)     :: h               ! the disc scale height [kpc]
    real(kind=8)                :: gas_timescales  ! The gas structuration time scale [Gyr]
    real(kind=8)                :: log_acc_rate
    real(kind=8)                :: log_h

    type(gas_type),intent(in)   :: acc_rate        ! the accretion rate onto the largest structure [10^11 Msun/Gyr]

    ! init to the maximaum value
    gas_timescales = t_str_max

    ! Accretion rate
    !
    log_acc_rate = gas_mass(acc_rate) ! in code unit
    if (log_acc_rate .gt. 0.d0) then
        log_acc_rate = log10(log_acc_rate)
        !
        ! Disc scale height
        !
        if (h .gt. 0.d0) then
            log_h = log10(h)
            !
            select case(timescale)
            case ('structuration','str','frag')
                ! Structuration timescale
                gas_timescales = locate2D(log_acc_rate,log_h,t_str, &
                                nAccRateBins,nDiscScaleBins, &
                                AccRatemin,DiscScalemin, &
                                dlAccRate,dlDiscScale)
            case ('formation','form')
                ! Structuration timescale
                gas_timescales = locate2D(log_acc_rate,log_h,t_form, &
                                nAccRateBins,nDiscScaleBins, &
                                AccRatemin,DiscScalemin, &
                                dlAccRate,dlDiscScale)  
            case ('equilibrium','eq')
                ! Structuration timescale
                gas_timescales = locate2D(log_acc_rate,log_h,t_eq, &
                                nAccRateBins,nDiscScaleBins, &
                                AccRatemin,DiscScalemin, &
                                dlAccRate,dlDiscScale)              
            case default
                write(message,'(a,a,a)') 'Keyword ', trim(timescale), ' not defined'
                call IO_print_error_message(message,only_rank=rank,called_by='gas_timescales')
                stop  ! stop the program
            end select                
            gas_timescales = (1.d1)**gas_timescales ! in Gyr
        end if
    end if

    return
  end function gas_timescales

  !*****************************************************************************************************************

  function Bonnor_Ebert_Mass(k)

    ! RETURN THE BONNOR EBERT MASS OF THE SYSTEM

    implicit none

    real(kind=8), intent(in)   :: k                 ! wavenumber [kpc^-1]
    real(kind=8)               :: Bonnor_Ebert_Mass ! [code unit]

    Bonnor_Ebert_Mass = (3.d0/2.d0)*sig_star**4./gravconst_code_unit**2./mu_star*(k/k_star)**(larson_mu_slope-4.*larson_sig_slope) ! in code unit

    return
  end function Bonnor_Ebert_Mass
  
  !*****************************************************************************************************************

  function Saturation_Mass(k, acc_rate)

    ! RETURN THE SATURATION MASS OF THE SYSTEM AT SCALE K

    implicit none

    real(kind=8), intent(in)   :: k               ! wavenumber [kpc^-1]       
    real(kind=8)               :: Saturation_Mass ! [code unit]
    
    type(gas_type),intent(in)  :: acc_rate ! instantaneous accretion rate onto the scale h

    Saturation_Mass = gas_mass(acc_rate) / ((2.d0/27.d0)*(pi*gravconst_code_unit**2.*mu_star**2./sig_star**3./k_star)*(k/k_star)**(6.*larson_sig_slope-larson_mu_slope-3.)) ! in code unit

    return
  end function Saturation_Mass

  !*****************************************************************************************************************

  function Maxwell_Boltzman_Vdist(v,param)

    ! RETURN F(v) FOR A MAXWELL BOLTMAN VELOCITY DISTRIBUTION

    implicit none

    real(kind=8), intent(in)                :: v                      ! velocity [code unit]
    real(kind=8)                            :: vr                     ! velocity [m/s]
    real(kind=8), intent(in), dimension(1)  :: param                  ! parameters array
                                                                      ! param(1) = T
    real(kind=8)                            :: Maxwell_Boltzman_Vdist ! [normalized]

    ! convert velocity input value (code unit) in m/s
    vr = v*vel_code_unit_2_mPerSec
    ! 1.195428d-5 is mu/(2pi kb) and 3.7556d-5 is mu/(2kb)
    ! We have to multiply by 'vel_code_unit_2_mPerSec' because, in the code, integrations are performed on code unit_scale
    ! and that f(v)dv is written and normalised here in the internationnal velocity unit system (m/s)
    Maxwell_Boltzman_Vdist = (4.d0*pi*(1.195428d-5)**(3./2.)*param(1)**(-3./2.)*vr**2.*exp(-3.7556d-5*vr**2./param(1)))*vel_code_unit_2_mPerSec

    return
  end function Maxwell_Boltzman_Vdist

  !*****************************************************************************************************************

  function Maxwell_Boltzman_Vdist_shifted(v,param)

    ! RETURN F(v) FOR A SHIFTED MAXWELL BOLTZMAN VELOCITY DISTRIBUTION

    implicit none

    real(kind=8), intent(in)                :: v                              ! velocity [code unit]
    real(kind=8), intent(in), dimension(2)  :: param                          ! parameters array
                                                                              ! param(1) = T
                                                                              ! param(2) = v_shift [code unit]
    real(kind=8)                            :: Maxwell_Boltzman_Vdist_shifted ! [normalized]! [normalized]
    real(kind=8)                            :: dP

    dP = 0.d0
    if (v-param(2) .gt. 0.d0) then
      dP = Maxwell_Boltzman_Vdist(v-param(2),(/param(1)/))
    end if

    Maxwell_Boltzman_Vdist_shifted = dP ! [normalized]

    return
  end function Maxwell_Boltzman_Vdist_shifted

  !*****************************************************************************************************************

  function Maxwell_Boltzman_Edist(v,param)

    ! RETURN F(v) FOR A MAXWELL BOLTZMAN KINETIC ENERGY DISTRIBUTION

    implicit none

    real(kind=8), intent(in)                :: v                      ! velocity [code unit]
    real(kind=8), intent(in), dimension(1)  :: param                  ! parameters array
                                                                      ! param(1) = T
    real(kind=8)                            :: Maxwell_Boltzman_Edist ! Energy [code unit]

    Maxwell_Boltzman_Edist = 5.d-1*v**2.*Maxwell_Boltzman_Vdist(v,param) ! [code unit]

    return
  end function Maxwell_Boltzman_Edist

  !*****************************************************************************************************************

  function Maxwell_Boltzman_Edist_shifted(v,param)

    ! RETURN F(v) FOR A SHIFTED MAXWELL BOLTZMAN KINETIC ENERGY DISTRIBUTION

    implicit none

    real(kind=8), intent(in)                :: v                              ! velocity [code unit]
    real(kind=8), intent(in), dimension(2)  :: param                          ! parameters array
                                                                              ! param(1) = T
                                                                              ! param(2) = v_shift [code unit]
    real(kind=8)                            :: Maxwell_Boltzman_Edist_shifted ! Energy [code unit]
    real(kind=8)                            :: dP

    dP = 0.d0
    if (v-param(2) .gt. 0.d0) then
      dP = Maxwell_Boltzman_Edist(v-param(2),(/param(1)/))
    end if

    Maxwell_Boltzman_Edist_shifted = dP ! [code unit]

    return
  end function Maxwell_Boltzman_Edist_shifted

  !***********************************************************************************************************

  function Maxwell_Boltzman_V_to_T(v)

    ! RETURN THE AVERAGE TEMPERATURE FOR A MAXWELL BOLTZMAN GAS WITH A MEAN VELOCITY v

    implicit none

    real(kind=8), intent(in)  :: v                       ! the mean velocity [code unit]
    real(kind=8)              :: v_                      ! the mean velocity [m/s]
    real(kind=8)              :: Maxwell_Boltzman_V_to_T ! the mean Temperature corresponding to the mean velocity

    ! convert velocity input value (code unit) in m/s
    v_ = v*vel_code_unit_2_mPerSec

    Maxwell_Boltzman_V_to_T = v_**2.*pi*mu/8.d0/kb  ! in K

    return
  end function Maxwell_Boltzman_V_to_T

  !***********************************************************************************************************

  function Maxwell_Boltzman_T_to_Eint(T)

    ! RETURN THE INTERNAL ENERGY IN A MAXWELL BOLTZMAN GAS WITH A MEAN TEMPERATURE T

    implicit none

    real(kind=8),intent(in)       :: T                          ! temperature of the gas [K]
    real(kind=8)                  :: Maxwell_Boltzman_T_to_Eint ! internal energy [per mass code unit]

    ! 2.0886929086943284d-2 3kbN/2 in code unit with N the number of particle of mass mu*mp in a mass code unit [10^11Msun]
    Maxwell_Boltzman_T_to_Eint = 2.0886929086943284d-2*T

    return
  end function Maxwell_Boltzman_T_to_Eint

  !***********************************************************************************************************

  function Maxwell_Boltzman_Eint_to_T(Eint)

    ! RETURN THE AVERAGE TEMPERATURE FOR A MAXWELL BOLTZMAN GAS WITH AN INTERNAL ENERGY Eint

    implicit none

    real(kind=8),intent(in)  :: Eint                       ! internal energy [per mass code unit]
    real(kind=8)             :: Maxwell_Boltzman_Eint_to_T ! the mean Temperature
    ! 4.787683224457894d1 is 2/(3kbN) in code unit with N the number of particle of mass mu*mp in a mass code unit [10^11Msun]
    Maxwell_Boltzman_Eint_to_T = 4.787683224457894d1*Eint

    return
  end function Maxwell_Boltzman_Eint_to_T

  !***********************************************************************************************************

  function Maxwell_Boltzman_Vp(T)

    ! RETURN THE MOST PROBABLE VELOCITY IN A MAXWELL BOLTZMAN GAS WITH A MEAN TEMPERATURE T

    implicit none

    real(kind=8),intent(in)       :: T                   ! temperature of the gas
    real(kind=8)                  :: Maxwell_Boltzman_Vp ! the most probable speed of the distribution [code unit]

    Maxwell_Boltzman_Vp = sqrt(2.d0*kb*T/mu)
    Maxwell_Boltzman_Vp = Maxwell_Boltzman_Vp/vel_code_unit_2_mPerSec ! [code unit]

    return
  end function Maxwell_Boltzman_Vp

  !***********************************************************************************************************
  !
  ! INTERFACE OPERATOR DEFINITIONS
  !
  !***********************************************************************************************************

  function gas_add_(g1,g2)

    ! RETURN  g1 + g2

    implicit none

    type(gas_type),intent(in)  :: g1,g2
    type(gas_type)             :: gas_add_

    call gas_void(gas_add_)
    call gas_copy(gas_add_,g1)
    call gas_add(gas_add_,g2)

    return
  end function gas_add_

  !***********************************************************************************************************

  function gas_add_array_(g1,g2)

    ! RETURN  g1 + g2

    implicit none

    type(gas_type),intent(in),allocatable   :: g1(:),g2(:)
    type(gas_type),allocatable              :: gas_add_array_(:)

    call gas_copy_array(gas_add_array_,g1)
    call gas_add_array(gas_add_array_,g2)

    return
  end function gas_add_array_

  !***********************************************************************************************************

  function gas_sub_(g1,g2)

    ! RETURN IF POSSIBLE g1 - g2

    type(gas_type),intent(in)    :: g1,g2
    type(gas_type)               :: gas_sub_

    call gas_void(gas_sub_)
    call gas_copy(gas_sub_,g1)
    call gas_sub(gas_sub_,g2)

    return
  end function gas_sub_

  !***********************************************************************************************************

  function gas_sub_array_(g1,g2)

    ! RETURN  g1 + g2

    implicit none

    type(gas_type),intent(in),allocatable  :: g1(:),g2(:)
    type(gas_type),allocatable             :: gas_sub_array_(:)

    call gas_copy_array(gas_sub_array_,g1)
    call gas_sub_array(gas_sub_array_,g2)

    return
  end function gas_sub_array_

  !***********************************************************************************************************

  function gas_scalar_multiply_left(gas,a)

    ! RETURN gas * a

    implicit none

    real(kind=8),intent(in)    :: a

    type(gas_type),intent(in)  :: gas
    type(gas_type)             :: gas_scalar_multiply_left

    call gas_void(gas_scalar_multiply_left)

    gas_scalar_multiply_left%Temp  = gas%Temp
    gas_scalar_multiply_left%f_str = gas%f_str
    gas_scalar_multiply_left%mass  = gas%mass * a
    gas_scalar_multiply_left%mZ    = gas%mZ   * a

    if (allocated(gas%elts)) then
       gas_scalar_multiply_left%elts(:) = gas%elts(:) * a
    else
       call IO_print_error_message('Try to use a non-existent ISM gas component ',only_rank=rank,called_by='gas_scalar_multiply_left/right')
       call IO_print_message('with',only_rank=rank,component='gas', &
            param_name=(/'gas%mass                 ','gas%mZ                   '/), real_param_val=(/gas%mass,gas%mZ/))
       stop ! stop the program
    endif

    return
  end function gas_scalar_multiply_left

  !***********************************************************************************************************

  function gas_arr_scalar_multiply_left(gas,a)

    ! RETURN gas * a

    implicit none

    integer(kind=4)            :: N      ! number of elements in gas
    integer(kind=4)            :: i,i_l  ! loop index

    real(kind=8),intent(in)    :: a

    type(gas_type),intent(in),allocatable  :: gas(:)
    type(gas_type),allocatable             :: gas_arr_scalar_multiply_left(:)

    ! copy the array
    gas_arr_scalar_multiply_left = gas

    ! Temp (not modified)
    gas_arr_scalar_multiply_left(:)%Temp = gas_arr_scalar_multiply_left(:)%Temp

    ! f_str (not modified)
    gas_arr_scalar_multiply_left(:)%f_str = gas_arr_scalar_multiply_left(:)%f_str

    ! apply numerical factor
    gas_arr_scalar_multiply_left(:)%mass = gas_arr_scalar_multiply_left(:)%mass * a
    gas_arr_scalar_multiply_left(:)%mZ   = gas_arr_scalar_multiply_left(:)%mZ  * a

    N   = size(gas)
    i_l = lbound(gas,dim=1)
    if (i_l .eq. 0) N = N -1
    do i_l = 1, N
       if (allocated(gas_arr_scalar_multiply_left(i)%elts)) then
          gas_arr_scalar_multiply_left(i)%elts(:) = gas_arr_scalar_multiply_left(i)%elts(:) * a
       else
          call IO_print_error_message('Try to use a non-existent ISM gas component ',only_rank=rank,called_by='gas_arr_scalar_multiply_left/right')
          call IO_print_message('with',only_rank=rank,component='gas', &
               param_name=(/'gas%mass                 ','gas%mZ                   '/), real_param_val=(/gas%mass,gas%mZ/))
          stop ! stop the program
       endif
    end do

    return
  end function gas_arr_scalar_multiply_left

  !***********************************************************************************************************

  function gas_scalar_multiply_right(a,gas)

    ! RETURN a * gas

    implicit none

    real(kind=8),intent(in)    :: a

    type(gas_type),intent(in)  :: gas
    type(gas_type)             :: gas_scalar_multiply_right

    gas_scalar_multiply_right = gas_scalar_multiply_left(gas,a)

    return
  end function gas_scalar_multiply_right

  !***********************************************************************************************************

  function gas_arr_scalar_multiply_right(a,gas)

    ! RETURN a * gas

    implicit none

    real(kind=8),intent(in)                :: a

    type(gas_type),intent(in),allocatable  :: gas(:)
    type(gas_type),allocatable             :: gas_arr_scalar_multiply_right(:)

    gas_arr_scalar_multiply_right = gas_arr_scalar_multiply_left(gas,a)

    return
  end function gas_arr_scalar_multiply_right

  !***********************************************************************************************************

end module gas
