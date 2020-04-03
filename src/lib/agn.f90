module agn

  use cooling  ! Contains cooling library (acces to Lambda(T,Z))
  use dust     ! Contains dust structure definition and dust exploitation functions
    
  public

  !*****************************************************************************************************************
  ! 
  !  OVERVIEW
  !
  !  agn_module defines the agn data structure univ(ts)%halo(ih)%galaxy%disc%agn
  !
  !  This module defines all properties and procedures assossiated to an agn object. An Active galaxy nuclei is composed
  !  of a Super Massive Black Hole (SMBH) associated with a gas/dust torus. The torus is fed by diffuse gas during mergers.
  !  Gas is then accreted onto the SMBH. A part of this gas allows to increase the mass of the SMBH. The other part is comverted into energy 
  !  The energy is then used to both eject (kinetic energy) and heat the gas (non-kinetic)  
  !  AGN is formed during merger even when the bulge stellar mass is larger than a given threshold
  !  AGN evolve into the centre of the disc. It evolution is closely linked to the diffuse gas content of this disc. 
  !  In the header of the module are defined output properties of an agn component (labels, units and formats)
  !
  !  This module starts with the definition of the agn type
  !
  !  MPI PROCESSES IN THIS MODULE
  !
  !   agn_send_data                               : send specific informations about a agn component
  !
  !   agn_receive_data                            : receive specific informations about a agn component
  !
  !  SUBROUTINES IN THIS MODULE
  !
  !   agn_read_agn_SEDs                           : Read SED templates associated to an active galaxy nuclei (currently Fritz+2006)
  !   
  !   agn_void                                    : set to null values all properties of an agn 
  !  
  !   agn_void_torus                              : void the agn torus (mass = 0.)
  !   
  !   agn_copy                                    : copy an agn component to an other 
  !
  !   agn_create_agn                              : create an SMBH
  !
  !   agn_evolve_I                                : compute the first part of the evolution scheme (predictor) compute optimal evolution time of the AGN component
  !   
  !   agn_evolve_II                               : evolve an AGN component, substract mass to the torus, add mass to the SMBH ...
  !
  !   agn_add_mass                                : add mass to the SMBH
  !
  !   agn_add_torus_mass                          : add mass to the SMBH torus
  !
  !   agn_sub_torus_mass                          : substract mass to the SMBH torus
  !
  !   agn_set_dynamical_time                      : set the value of the accretion dynamical time
  !
  !   agn_reset_t_since_last_merger               : reset (set to 0.) the AGN clock that runs between merger events
  !
  !   agn_spectrum                                : Build the Spectral Energy Distribution associated to the gas/dust torus heating by the AGN
  !
  !  FUNCTIONS IN THIS MODULE
  !
  !   agn_mass                                    : return the mass of the SMBH, the mass of the torus, or the total mass (SMBH + torus)
  !
  !   isolated_agn_velocity_  (CORE VERSION)      : return the circular velocity (at radius r) associated to an isolated SMBH              
  !
  !   agn_merge                                   : return the remnant agn by merging two agns
  !
  !   agn_compute_agn_infall_rate                 : return the total SMBH infall-rate
  !
  !   agn_compute_Bondi_accretion_rate            : return the Bondi accretion rate associated to a given SMBH
  !
  !   agn_compute_Eddington_accretion_rate        : return the Eddington accretion rate limit associated to a given SMBH
  !
  !   agn_compute_agn_accretion_rate              : return the (real) SMBH accretion rate (the SMBH increases in mass through this rate)
  !
  !   agn_compute_agn_ejecta_rate                 : return the gas ejection rate associated to the SMBH activity
  !
  !   agn_compute_agn_velocity_wind               : return the velocity of the wind produced by SMBH activity
  !
  !   agn_compute_agn_turbulent_heating_power     : return the instantaneous turbulent heating power produced by the SMBH actcivity   
  !
  !   agn_compute_agn_non_kinetic_power           : return the instantaneous non-kinetic power produced by SMBH activity
  !
  !   agn_compute_agn_thermal_power               : return the instantaneous thermal power produced by the SMBH actcivity
  !
  !   agn_compute_agn_rad_power                   : return the instantaneous radiative power (non thermal and non kinetic power) produced by the SMBH activity
  !
  !   agn_compute_L                               : return the instantaneous luminosity of the AGN
  !   
  !   agn_compute_L_edd                           : return the Eddington luminosity of a given SMBH
  !
  !  PRINTING PROCEDURES
  !
  !   agn_load_agn_data                           : create the agn output properties list
  !
  !*****************************************************************************************************************
  
  type agn_type
    ! SMBH properties
    real(kind=8)             :: mass                ! the black hole mass
    real(kind=8)             :: L                   ! effective luminosity 
    real(kind=8)             :: L_edd               ! Eddington luminosity 
    real(kind=8)             :: t_dyn               ! dynamical time associated to the accretion coming from the torus
    real(kind=8)             :: t_since_last_merger ! time elapsed since the last merger event
    ! torus
    type(gas_type)           :: torus               ! gas in the torus around to the black hole
    type(dust_type)          :: dust                ! dust component of the torus
    ! transfer rates
    type(gas_type)           :: acc_rate            ! effective accretion onto the SMBH
    type(gas_type)           :: infall_rate         ! infall rate onto the SMBH, coming from the torus to the SMBH
                                                    ! a fraction allows to increase the mass of the SMBH (acc_rate), 
                                                    ! an the other part is ejected, an other is converted in energy
    ! spectrum
    real(kind=4),allocatable :: spectrum(:)         ! agn spectrum Lsun
  end type agn_type
  
  ! parameters used to select the best AGN dust emission spectrum into the library
  integer(kind=4)           :: ntauBins     ! number of equatorial optical depth bins takes into account
  integer(kind=4)           :: nvaBins      ! number of view angles takes into account
  
  real(kind=4),allocatable  :: tauBins(:)      ! list of equatorial optical depth
  real(kind=4),allocatable  :: vaBins(:)       ! list of view angles 
  real(kind=4),allocatable  :: agn_sed(:,:,:)  ! Active Galaxy Nuclei spectrum (nwaves,ntauBins,nvaBins)
  
  real(kind=8),parameter    :: T_agn                = 5.e6    ! in K, fixed temperature of the gas/dust torus
  real(kind=8),parameter    :: AGN_kinetic_fraction = 1.d-3   ! fraction of the total AGN power dedicated to kinetic power 
  real(kind=8),parameter    :: AGN_thermal_fraction = 6.d-1   ! fraction of the non-kinetic AGN power dedicated to thermal power  
  real(kind=8),parameter    :: AGN_mass_energy_conv = 1.d-1   ! fraction of the mass acrreted onto the SMBH that is converted in energy
  real(kind=8),parameter    :: r_torus              = 1.d-2   ! kpc = 10pc, allows to compute the dynamical accretion timescale
  
  ! printable properties for agn structure
  integer(kind=4),parameter :: nb_agn_field = 4 ! Number of agn properties saved in the main FITS output file
  ! Name of each output colomn data
  character(len=ttype_len),dimension(nb_agn_field) :: ttype_agn = (/'agn_mass              ','agn_L                 ','agn_L_Edd             ','agn_acc_rate_in       '/)   
  ! Physical unit of each column data
  character(len=tunit_len),dimension(nb_agn_field) :: tunit_agn = (/'M_sun       ','L_sun       ','L_sun       ','M_sun/yr    '/)   
  ! Data type of each column data
  character(len=tform_len),dimension(nb_agn_field) :: tform_agn = (/'1E','1E','1E','1E'/) 

contains

  !*****************************************************************************************************************
  ! 
  ! MPI PROCESSES
  !
  !*****************************************************************************************************************

#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------  
  subroutine agn_send_data(agn)                                     
  
    ! SEND SPECIFIC INFORMATIONS ABOUT A AGN COMPONENT
  
    implicit none

    logical                     :: go_down
    
    type(agn_type),intent(in)   :: agn         ! agn component
    
    ! data are sent by physical process and receive by luminous process 
    
    go_down = .false. ! by default no AGN 
    
    if (agn_mass(agn) .gt. 0.d0) then
        !
        ! there is a SMBH
        go_down = .true.
        call MPI_SEND(go_down,1,MPI_LOGICAL,rank+1,agn_tag+1,MPI_COMM_WORLD,ierror)
        !
        ! send agn mass
        call MPI_SEND(agn%mass,1,MPI_REAL8,rank+1,agn_tag+2,MPI_COMM_WORLD,ierror)
        !
        ! send agn bolometric luminosity
        call MPI_SEND(agn%L,1,MPI_REAL8,rank+1,agn_tag+3,MPI_COMM_WORLD,ierror)
        !
        ! send dust properties of the torus
        call dust_send_data(agn%dust)
    else
        !
        ! send exit loop order
        call MPI_SEND(go_down,1,MPI_LOGICAL,rank+1,agn_tag+1,MPI_COMM_WORLD,ierror)
    end if
    
    return
  end subroutine agn_send_data 
  
  !*****************************************************************************************************************
  
  subroutine agn_receive_data(agn)                                     
  
    ! RECEIVE SPECIFIC INFORMATIONS ABOUT A AGN COMPONENT
  
    implicit none

    logical                     :: go_down
    
    type(agn_type),intent(out)  :: agn
    
    ! data are sent by physical process and receive by luminous process 
    
    call agn_void(agn) ! init
    
    ! receive go_down order
    call MPI_RECV(go_down,1,MPI_LOGICAL,rank-1,agn_tag+1,MPI_COMM_WORLD,statut,ierror)
    
    if (go_down) then
        !
        ! receive agn informations
        ! receive agn mass
        call MPI_RECV(agn%mass,1,MPI_REAL8,rank-1,agn_tag+2,MPI_COMM_WORLD,statut,ierror)
        !
        ! receive bolometric luminosity
        call MPI_RECV(agn%L,1,MPI_REAL8,rank-1,agn_tag+3,MPI_COMM_WORLD,statut,ierror)
        !
        ! receive dust properties of the agn torus
        call dust_receive_data(agn%dust)
    end if
    
    return
  end subroutine agn_receive_data 
! ------------------------------------------------- 
#endif
! LUMINOUS_PROCESSES

  !*****************************************************************************************************************
  !
  ! SUBROUTINES
  !
  !*****************************************************************************************************************

#ifdef AGN_SPECTRUM
! -------------------------------------------------   
  subroutine agn_read_agn_SEDs
  
    ! READ AGN SED TEMPLATES
    
    implicit none
    
    integer(kind=4)          :: t               ! index loop under equatorial optical depth
    integer(kind=4)          :: l               ! index loop under wavelenghts
    integer(kind=4)          :: nWaves_agn      ! local parameter (number of wavelenght in the agn SED library)
     
    character(MAXPATHSIZE)   :: filename
    character(MAXPATHSIZE)   :: message         ! a message  
    character(MAXPATHSIZE)   :: line
    
    if (main_process .or. luminous_process) then
    
        call IO_print_message('agn_read_agn_SEDs')
        
        ! Build the input filename
        write(filename,'(a,a)') trim(input_path), '/agn_sed.in'

        write(message,'(a,a)') 'Load data from : ', 'agn_sed.in'
        call IO_print_message(message)
        ! Open the library file
        open(unit = agn_sed_unit, file = filename, status = 'old')
           
        do
            read(agn_sed_unit, '(a)', end = 2) line
            if (trim(line) .eq. 'START') then
                !
                ! all lines before the START keywork are header lines: they contain informations about data 
                read(agn_sed_unit,*) nWaves_agn, ntauBins, nvaBins
                !
                if (main_process) then
                    !
                    ! The main process checks if nWaves = nWaves_agn
                    ! If it is not the case the main process kills the run
                    if (nWaves_agn  .ne. nWaves) then
                        !
                        call IO_print_error_message('nWaves_agn /= nWaves',only_rank=rank,called_by='agn_read_agn_SEDs') 
                        stop ! stop the program 
                    end if
                end if
                !
                ! allocate arrays
                ! equatorial optical depth
                allocate(tauBins(ntauBins))
                ! view angles
                allocate(vaBins(nvaBins))
                ! AGN SEDs
                allocate(agn_sed(nWaves,ntauBins,nvaBins))
                !
                read(agn_sed_unit,*) ! skip wavelenght list
                read(agn_sed_unit,*) tauBins ! equatorial optical depth
                read(agn_sed_unit,*) vaBins  ! view angles
                !
                read(agn_sed_unit,*) ! skip a line #
                !
                ! read SEDs
                do t = 1, ntauBins
                    do l = 1, nWaves
                        read(agn_sed_unit,*) agn_sed(l,t,1:nvaBins) ! Lsun/Lsun
                    end do
                end do
                !
                if (main_process) then
                    ! 
                    ! Print some informations 
                    call IO_print_message('use', param_name=(/'ntauBins                 ','nvaBins                  '/),&
                        int_param_val=(/ntauBins,nvaBins/))
                end if
                exit  ! quit do loop 
            end if
            !
            if (line(1:1) .eq. '#') then
                !
                cycle ! header or something like this (skip)
            else
                !
                call IO_print_error_message('Impossible to read the line',only_rank=rank,called_by='agn_read_agn_SEDs') 
                stop ! stop the program
            end if
        end do
    end if  
        
    2 close(abs_sca_unit)

    return
  end subroutine agn_read_agn_SEDs
! ------------------------------------------------- 
#endif
! AGN_SPECTRUM
   
  !*****************************************************************************************************************
    
  subroutine agn_void(agn)

    ! INIT OR VOID AN AGN COMPONENT

    implicit none

    type(agn_type),intent(inout)   :: agn   ! the agn component

    ! SMBH properties
    agn%mass                = 0.d0 
    agn%L                   = 0.d0
    agn%L_edd               = 0.d0 
    agn%t_dyn               = 0.d0
    agn%t_since_last_merger = 0.d0
    ! torus
    call agn_void_torus(agn)
    ! transfer rates
    call gas_void(agn%acc_rate)
    call gas_void(agn%infall_rate)
    
#ifdef AGN_SPECTRUM
! ------------------------------------------------- 
    ! spectrum
    if (allocated(agn%spectrum)) agn%spectrum(:) = 0.d0     
! ------------------------------------------------- 
#endif
! AGN_SPECTRUM

    return
  end subroutine agn_void

  !*****************************************************************************************************************

  subroutine agn_void_torus(agn)

    ! VOID THE AGN ASSOCIATED TORUS
    ! void both the gas and the dust component
    
    implicit none

    type(agn_type),intent(inout)  :: agn

    call gas_void(agn%torus)                     ! void the gas component
    call dust_void(agn%dust,init_geom='torus ')  ! void the dust component

    return
  end subroutine agn_void_torus

  !*****************************************************************************************************************

  subroutine agn_copy(agn1,agn2)

    ! COPY agn2 INTO agn1

    implicit none

    type(agn_type),intent(inout)   :: agn1   ! an agn component
    type(agn_type),intent(in)      :: agn2   ! an other agn component

    ! SMBH properties
    agn1%mass                = agn2%mass
    agn1%L                   = agn2%L
    agn1%L_edd               = agn2%L_edd  
    agn1%t_dyn               = agn2%t_dyn 
    agn1%t_since_last_merger = agn2%t_since_last_merger
    ! torus 
    call gas_copy(agn1%torus,agn2%torus)
    call dust_copy(agn1%dust,agn2%dust)
    ! transfer rates
    call gas_copy(agn1%acc_rate,agn2%acc_rate)
    call gas_copy(agn1%infall_rate,agn2%infall_rate)
    !
#ifdef AGN_SPECTRUM
! -------------------------------------------------     
    if (allocated(agn2%spectrum)) then
        !
        if (.not. allocated(agn1%spectrum)) then
            !
            allocate(agn1%spectrum(nWaves))
            agn1%spectrum = agn2%spectrum
        else
            !
            call IO_print_error_message('Try to erase pre-existing agn%spectrum', &
               only_rank = rank, called_by = 'agn_copy')
            stop ! kill the run
        end if
    end if
! ------------------------------------------------- 
#endif
! AGN_SPECTRUM    
    !
    return
  end subroutine agn_copy

  !*****************************************************************************************************************

  subroutine agn_create_agn(agn,M0)

    ! CREATE AN AGN OF MASS M0

    implicit none

    real(kind=8),intent(in)         :: M0   ! the initial mass

    type(agn_type),intent(inout)    :: agn  ! the agn component
    
    call agn_void(agn)                      ! init
    agn%mass  = M0                          ! set the initial mass of the SMBH
    agn%L_edd = agn_compute_L_edd(agn)      ! compute the Eddington Luminosity

    return
  end subroutine agn_create_agn
    
  !*****************************************************************************************************************
  
  subroutine agn_evolve_I(agn,dt_optim)
  
    ! COMPUTE THE FIRST PART OF THE AGN EVOLUTION SCHEME (PREDICTOR PART)
    ! Compute the overall infall rate and deduce the optimal evolution timescale 
    !   The torus is instantaneously fed during merger events
    !   the secular evolution of an AGN is therfore only linked to the infall onto the SMBH
    
    implicit none
        
    real(kind=8),intent(out)     :: dt_optim            ! the optimal evolution time-scale of the AGN
    real(kind=8)                 :: infall_rate         ! the global infall rate onto the SMBH. In the current scheme, the global infall rate is divided in two parts
                                                        !    - the real accretion rate onto the SMBH and the outflow rate
    type(gas_type)               :: in_rate,out_rate
                                                                                
    type(agn_type),intent(inout) :: agn                 ! the agn component

#ifdef PRINTALL 
! -------------------------------------------------
    call IO_print_message('agn_evolve_I',only_rank=rank,component='agn')
! -------------------------------------------------
#endif
    
    dt_optim = -1.d0        ! init 
    call gas_void(in_rate)  ! init
    call gas_void(out_rate) ! init
    
    if (agn_mass(agn) .gt. 0.d0) then
        !
        ! The SMBH exists and is surrouding by a torus containing some gas
        ! compute global infall rate
        infall_rate = agn_compute_AGN_infall_rate(agn)
        ! set new value of the agn infall rate
        agn%infall_rate = infall_rate*gas_signature(agn%torus,apply_as='rate_builder')
        ! compute and set the real accretion rate onto the SMBH
        agn%acc_rate    = (1.d0-AGN_mass_energy_conv)*agn_compute_AGN_accretion_rate(agn)*gas_signature(agn%torus,apply_as='rate_builder')
        !
        if (infall_rate .gt. 0.d0) then
            !
            ! set out rate, in_rate = 0.
            out_rate = agn%infall_rate
            ! The optimal evolution is computed assuming a max evolution gap of 30%
            dt_optim = gas_dt_optim(agn%torus,in_rate,out_rate,pp=3.d-1,called_by='agn_evolve_I')
        end if
    end if  
    
#ifdef PRINTALL 
! -------------------------------------------------
    call IO_print_message('agn_evolve_I ... done',only_rank=rank,component='agn')
! -------------------------------------------------
#endif  

    return
  end subroutine agn_evolve_I
    
  !*****************************************************************************************************************
  
  subroutine agn_evolve_II(agn,dt,incl)
  
    ! EVOLVE AN AGN COMPONENT DURING dt
    ! Substract mass to the torus, add mass to the SMBH and compute luminosities
  
    implicit none
    
    real(kind=8),intent(in)      :: dt   ! evolution time
    real(kind=8),intent(in)      :: incl ! disc inclination (create a local copy for dust extinction computation)

    type(agn_type),intent(inout) :: agn  ! the agn component

#ifdef PRINTALL 
! -------------------------------------------------
    call IO_print_message('agn_evolve_II',only_rank=rank,component='agn')
! -------------------------------------------------
#endif
    
    if (agn_mass(agn) .gt. 0.d0) then
        ! 
        ! An agn exists and its torus contains gas and dust
        ! update t_since_last_merger
        agn%t_since_last_merger = agn%t_since_last_merger + dt
        ! 
        if (gas_mass(agn%infall_rate) .gt. 0.d0) then
            ! 
            call agn_sub_torus_mass(agn,agn%infall_rate*dt)     ! substract the mass to the torus
            call agn_add_mass(agn,agn%acc_rate*dt)              ! add mass to the SMBH
            ! compute and set the new Eddington luminosity
            agn%L_edd = agn_compute_L_edd(agn) 
            ! compute and set instantaneous luminosity
            agn%L = agn_compute_L(agn)
        else
            !
            agn%L = 0.d0
        end if
        ! dust component
        call dust_evolve(agn%dust,agn%torus,r_torus,incl=incl)
    end if  

#ifdef PRINTALL 
! -------------------------------------------------
    call IO_print_message('agn_evolve_II ... done',only_rank=rank,component='agn')
! -------------------------------------------------
#endif
    
    return
  end subroutine agn_evolve_II 
   
  !*****************************************************************************************************************  

  subroutine agn_add_mass(agn,mass)

    ! ADD MASS TO THE AGN COMPONENT

    implicit none
      
    type(gas_type),intent(in)        :: mass   ! the mass  
    type(agn_type),intent(inout)     :: agn    ! the agn component
     
    if (agn%mass .gt. 0.d0) then
      !
      ! add mass to the SMBH
      agn%mass = agn%mass + gas_mass(mass)
    else  
      !
      if (gas_mass(mass) .gt. 0.d0) then
        !
        call IO_print_error_message('Add mass to an inexisting SMBH',only_rank=rank,called_by='agn_add_mass')
        stop ! stop the program
      end if
    end if

    return
  end subroutine agn_add_mass

  !*****************************************************************************************************************

  subroutine agn_add_torus_mass(agn,mass)

    ! ADD MASS TO THE TORUS RESERVOIR OF THE AGN
    ! During a merger event, gas is transfered into the torus

    implicit none

    type(gas_type),intent(in)       :: mass   ! the gas mass
    type(agn_type),intent(inout)    :: agn    ! the agn component
      
    if (agn_mass(agn) .gt. 0.d0) then
      !
      ! The SMBHG and the torus exist
      ! add mass
      agn%torus = agn%torus + mass
      ! set temperature to T_agn
      call gas_set_component(agn%torus,T_agn,component='Temp')
    else
      !
      call IO_print_error_message('Add mass into a torus without SMBH',only_rank=rank,called_by='agn_add_torus_mass')  
      stop ! stop the program
    end if

    return
  end subroutine agn_add_torus_mass

  !*****************************************************************************************************************

  subroutine agn_sub_torus_mass(agn,gas)
      
    ! SUBSTRACT MASS TO THE TORUS RESERVOIR OF THE AGN

    implicit none

    type(gas_type),intent(in)       :: gas   ! the mass accreted
    type(agn_type),intent(inout)    :: agn   ! the agn component
    
    if (agn_mass(agn) .gt. 0.d0) then
        !
        if (gas_mass(gas) .gt. gas_mass(agn%torus)) then
            !
            call IO_print_error_message('Substract too much gas from the agn torus',only_rank=rank,called_by='agn_sub_torus_mass')
            stop
        end if
        agn%torus = agn%torus - gas
    else
      !
      call IO_print_error_message('Substract gas to a torus without SMBH',only_rank=rank,called_by='agn_sub_torus_mass')  
      stop ! stop the program
    end if

    return
  end subroutine agn_sub_torus_mass
  
  !*****************************************************************************************************************

  subroutine agn_set_dynamical_time(agn,t_dyn)
  
    ! SET THE DYNAMICAL TIME OF THE SMBH
    
    real(kind=8),intent(in)        :: t_dyn    ! the dynamical time [r_torus / Vdisc(r_torus)]
    
    type(agn_type),intent(inout)   :: agn      ! the agn component
  
    agn%t_dyn = t_dyn
    
    return  
  end subroutine agn_set_dynamical_time
  
  !*****************************************************************************************************************

  subroutine agn_reset_t_since_last_merger(agn)
  
    ! RESET THE TIME ELAPSED SINCE THE LAST MERGER EVENT
    
    type(agn_type),intent(inout)   :: agn      ! the agn component
  
    agn%t_since_last_merger = 0.d0
    
    return  
  end subroutine agn_reset_t_since_last_merger
  
  !*****************************************************************************************************************
  
#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------   
#ifdef AGN_SPECTRUM
! -------------------------------------------------    
  subroutine agn_spectrum(agn)
  
    ! BUILD THE SPECTRUM ASSOCIATED TO AN AGN COMPONENT
    
    implicit none
    
    integer(kind=4)              :: itau       ! index of the best AGN SED template (tau)
    integer(kind=4)              :: iva        ! index of the best AGN SED template (view angle)
    
    real(kind=8)                 :: tau_8mic   ! the effective extinction measured at 8 microns
    real(kind=8)                 :: view_angle ! effective view angle of the disc
    real(kind=8)                 :: mu_A, mu_B, mu_C, mu_D
      
    type(agn_type),intent(inout) :: agn        ! An AGN component
    
    ! allocate agn%spectrum
    if (.not. allocated(agn%spectrum)) allocate(agn%spectrum(nWaves))  ! create
    agn%spectrum(:) = 0.d0                                           ! init
    
    if (agn%L .le. 0.d0) return  ! no SMBH activity
    
    ! compute dust extinction in the torus
    ! Fritz+2006 assume a mathis+1987 extinction curve and compute
    ! the equatorial optical depth at 9.7microns
    ! Here we extract the extinction value at 8microns
    call dust_compute_effective_extinction(agn%dust,tau_8mic=tau_8mic)  
    !
    ! select the best SED template as function of the effective extinction
    itau = 1
    do while ((tauBins(itau) < tau_8mic) .and. (itau < ntauBins))
        itau = itau +1
    end do
    ! We select the AGN SED template with the closest tauBins
    itau = max(1,itau-1)
    !
    ! select the best SED template as function of the view angle
    view_angle = max(minval(vaBins),(pi/2.-agn%dust%incl)*1.8d2/pi)
    iva = 1
    do while ((vaBins(iva) < view_angle) .and. (iva < nvaBins))
        iva = iva +1
    end do
    ! We select the AGN SED template with the closest vaBins
    iva = max(1,iva-1)
    !
    ! The effective SED is then interpolated between the 4 available SEDs
    ! compute weighting parameters
    mu_A = sqrt((tauBins(itau)-tau_8mic)**2. + (vaBins(iva)-view_angle)**2.) 
    mu_B = sqrt((tauBins(itau+1)-tau_8mic)**2. + (vaBins(iva)-view_angle)**2.) 
    mu_C = sqrt((tauBins(itau+1)-tau_8mic)**2. + (vaBins(iva+1)-view_angle)**2.) 
    mu_D = sqrt((tauBins(itau)-tau_8mic)**2. + (vaBins(iva+1)-view_angle)**2.) 
    ! build the spectrum
    agn%spectrum = real(mu_A*agn_sed(:,itau,iva) + mu_B*agn_sed(:,itau+1,iva) + &
         mu_C*agn_sed(:,itau+1,iva+1) + mu_D*agn_sed(:,itau,iva+1),4)
    !
    ! The AGN spectrum is normalized to the AGN effective luminosity
    agn%spectrum = agn%spectrum*real(agn%L/L_Sun_in_code_unit/(mu_A + mu_B + mu_C + mu_D),4)  
    
  end subroutine agn_spectrum
! ------------------------------------------------- 
#endif
! AGN_SPECTRUM
! ------------------------------------------------- 
#endif
! LUMINOUS_PROCESSES


  !*****************************************************************************************************************
  !
  ! FUNCTIONS
  !
  !*****************************************************************************************************************

  function agn_mass(agn,component)

    ! RETURN THE AGN MASS

    implicit none

    character(*),intent(in),optional :: component  ! allow to select between different masses (SMBH mass, gas mass in the torus or all)
    character(MAXPATHSIZE)           :: message    ! a message to display

    real(kind=8)                     :: agn_mass   ! The function result

    type(agn_type),intent(in)        :: agn        ! the agn component
    
    if (present(component)) then
      !
      select case (trim(component)) 
      case ('agn','SMBH','black-hole')
        agn_mass = agn%mass
      case ('torus')
        agn_mass = gas_mass(agn%torus)
      case ('all')
        agn_mass = agn%mass + gas_mass(agn%torus)
      case default
        write(message,'(a,a,a)') 'Keyword ', trim(component), ' not defined'
        call IO_print_error_message(message,only_rank=rank,called_by='agn_mass')    
        stop  ! stop the program
      end select
    else
      agn_mass = agn%mass + gas_mass(agn%torus)
    end if
       
    return
  end function agn_mass
  
  !*****************************************************************************************************************
  
  function isolated_agn_velocity_(r,param)
    
    ! COMPUTE THE CIRCULAR VELOCITY ASSOCIATED TO AN ISOLATED SMBH

    implicit none
    
    real(kind=8), intent(in)   :: r                       ! radius (in kpc)
    real(kind=8),intent(in)    :: param(1)                ! parameter array : param(1) = agn%mass

    real(kind=8)               :: isolated_agn_velocity_  ! in kpc/Gyr 
        
    isolated_agn_velocity_ = sqrt(gravconst_code_unit*param(1)/r) ! in kpc/Gyr

    return
  end function isolated_agn_velocity_

  !*****************************************************************************************************************

  function agn_merge(agn1,agn2)

    ! COMPUTE AGN MERGER EVENT

    implicit none
      
    type(agn_type),intent(in)       :: agn1      ! an agn component
    type(agn_type),intent(in)       :: agn2      ! an other agn component
    type(agn_type)                  :: agn_merge ! the remnent agn
     
    call agn_void(agn_merge) ! init
      
    ! reset t_since_last_merger
    agn_merge%t_since_last_merger = 0.d0
    ! add SMBH mass
    agn_merge%mass  = agn1%mass + agn2%mass               
    ! add mass in the torus 
    agn_merge%torus = agn1%torus + agn2%torus             
    ! the dust content will be computed in the next call of agn_evolve_II procedure in disc_evolve_II
    ! update the Eddington luminosity
    agn_merge%L_edd = agn_compute_L_edd(agn_merge)
      
    return
  end function agn_merge
  
  !*****************************************************************************************************************

  function agn_compute_agn_infall_rate(agn)

    ! COMPUTE AGN INFALL RATE

    implicit none
    
    real(kind=8)                 :: torus_acc_rate   ! accretion rate (torus mode)
    real(kind=8)                 :: M_torus, M_BH    ! the torus gas mass and the SMBH mass
               
    real(kind=8)                 :: Bondi_acc_rate   ! accretion rate (Bondi mode)
    real(kind=8)                 :: Edd_acc_rate     ! accretion rate (Eddington limit)
    real(kind=8)                 :: agn_compute_agn_infall_rate

    type(agn_type),intent(in)    :: agn              ! an agn component
    
    agn_compute_agn_infall_rate = 0.d0 ! init
       
    if (agn_mass(agn) .le. 0.d0) return              ! no black hole and therefore no accretion process

    torus_acc_rate = 0.d0 ! init
    Bondi_acc_rate = 0.d0 ! init
    Edd_acc_rate   = 0.d0 ! init
    
    M_torus = agn_mass(agn,component='torus')
    M_BH    = agn_mass(agn,component='SMBH')
    !
    if (M_torus .gt. M_BH_min) then
        !
        ! compute torus accretion rate
        ! dynamical time computed with the rotation speed at r_torus of the SMBH
        torus_acc_rate = 5.d-1*M_torus/agn%t_dyn                   ! in code unit 10^11 M_sun / Gyr
        !
        ! compute Bondi accretion rate
        Bondi_acc_rate = agn_compute_Bondi_accretion_rate(agn)     ! in code unit 10^11 M_sun / Gyr
        !
        ! compute Eddington accretion rate
        Edd_acc_rate   = agn_compute_Eddington_accretion_rate(agn) ! in code unit 10^11 M_sun / Gyr
    end if
    ! compute real accretion rate onto the black hole
    ! the black hole accretion rate cannot be greater than the Eddington luminosity
    ! and cannot be smaller than the Bondi_acc_rate
    agn_compute_agn_infall_rate = min(max(torus_acc_rate,Bondi_acc_rate),Edd_acc_rate)      
    ! 
    ! Check the result and crash the code if:
    !
    ! agn_compute_agn_infall_rate is NAN
    if (is_NaN(agn_compute_agn_infall_rate)) then
      !
      call IO_print_error_message('AGN infall rate is NAN', &
            only_rank = rank, called_by = 'agn_compute_agn_infall_rate')
      call IO_print_message('used',only_rank=rank,component='agn',& 
        param_name=(/'infall rate              ','M_torus                  ','torus_acc_rate           ',&
                     'Bondi_acc_rate           ','Edd_acc_rate             '/), &
        real_param_val =(/agn_compute_agn_infall_rate,M_torus,torus_acc_rate,Bondi_acc_rate,Edd_acc_rate/))
      stop ! stop the program  
    end if
    !
    ! agn_compute_agn_infall_rate is smaller than 0.d0
    if (agn_compute_agn_infall_rate .lt. 0.d0) then   ! check positivity of the ejecta
      !
      call IO_print_error_message('AGN infall rate < 0.', &
            only_rank = rank, called_by = 'agn_compute_agn_infall_rate')
      call IO_print_message('used',only_rank=rank,component='agn',& 
        param_name=(/'infall rate              ','M_torus                  ','torus_acc_rate           ',&
                     'Bondi_acc_rate           ','Edd_acc_rate             '/), &
        real_param_val =(/agn_compute_agn_infall_rate,M_torus,torus_acc_rate,Bondi_acc_rate,Edd_acc_rate/))
      stop ! stop the program  
    end if

    return
  end function agn_compute_agn_infall_rate

  !*****************************************************************************************************************
  
  function agn_compute_Bondi_accretion_rate(agn)

    ! RETURN THE BONDI ACCRETION RATE ONTO A GIVEN SMBH

    implicit none

    real(kind=8)                 :: Zg                                ! gas metalicity
    real(kind=8)                 :: agn_compute_Bondi_accretion_rate  ! in code unit [10^11 Msun / Gyr]
    real(kind=8)                 :: rate

    type(agn_type),intent(in)    :: agn                               ! the agn component

    rate = 0.d0   ! init
    Zg   = gas_metallicity(agn%torus)
      
    if (agn_mass(agn) .gt. 0.d0) then
      ! 
      rate = 3.d0*pi*gravconst*mu/4.d0           ! in m^3/s^2
      rate = rate*kb*T_agn/Lambda(T_agn,Zg)      ! in s^-1  --> WARNING : [cooling_curves_Lambda] = Joule*m^3/s = kg*m^2/s^3)
                                                 !                                           [kb] = Joule/K     = kg.m^2/s^2/K
      rate = rate*Gyr_in_s                       ! in Gyr^-1
      rate = rate*agn_mass(agn)                  ! in 10^11 M_sun/Gyr : code unit for accretion rate
    end if
      
    agn_compute_Bondi_accretion_rate = rate
    
    ! NAN and <0 cases are treated in agn_compute_agn_infall_rate.

    return
  end function agn_compute_Bondi_accretion_rate

  !*****************************************************************************************************************

  function agn_compute_Eddington_accretion_rate(agn)

    ! COMPUTE THE EDDINGTON RATE OF THE AGN COMPONENT

    implicit none

    real(kind=8)                 :: agn_compute_Eddington_accretion_rate  ! in code unit [10^11 Msun / Gyr]
    real(kind=8)                 :: rate

    type(agn_type),intent(in)    :: agn                                   ! the agn component

    rate = 0.d0   ! init
      
    if (agn_mass(agn) .gt. 0.d0) then
      ! 
      rate  = (1.d0 + agn_ejecta_efficiency)*agn_compute_L_edd(agn)/(light_speed_code_unit**2.) ! in code unit 10^11 Msun / Gyr
    end if
      
    agn_compute_Eddington_accretion_rate = rate
    
    ! NAN and <0 cases are treated in agn_compute_agn_infall_rate
      
    return  
  end function agn_compute_Eddington_accretion_rate
    
  !*****************************************************************************************************************

  function agn_compute_agn_accretion_rate(agn)

    ! COMPUTE REAL AGN ACCRETION RATE (real accretion onto the SMBH)

    implicit none

    real(kind=8)              :: f_agn
    real(kind=8)              :: agn_compute_agn_accretion_rate

    type(agn_type),intent(in) :: agn       ! an agn component
    
    agn_compute_agn_accretion_rate = 0.d0  ! init
    !
    if (agn_mass(agn) .le. 0.d0) return              ! no black hole and therefore no accretion process
    !
    f_agn = 1.d0 / (1.d0 + agn_ejecta_efficiency)
    agn_compute_agn_accretion_rate = f_agn*agn_compute_agn_infall_rate(agn)
    
    ! NAN and <0 cases are treated in agn_compute_agn_infall_rate

    return
  end function agn_compute_agn_accretion_rate

  !*****************************************************************************************************************

  function agn_compute_agn_ejecta_rate(agn,Vesc)

    ! COMPUTE AGN EJECTA RATE PRODUCED BY THE AGN
    ! ejecta rate a composed of two parts:
    !  i) a part coming from the torus, the infall rate is divided in two part (acc, eject)
    !  ii) a part coming from the unstructred gas from the disc due to a transfer of momentum

    implicit none

    real(kind=8),intent(in)      :: Vesc ! galaxy escape velocity, must be given in code unit
    real(kind=8)                 :: f_agn
    real(kind=8)                 :: agn_compute_agn_ejecta_rate

    type(agn_type),intent(in)    :: agn  ! an agn component
    
    agn_compute_agn_ejecta_rate = 0.d0   ! init
    
    if (agn_mass(agn) .le. 0.d0) return  ! no black hole and therefore no accretion process
    
    if (Vesc .le. 0.d0) then
      !
      call IO_print_error_message('Vesc < 0.',only_rank=rank,called_by='agn_compute_disc_AGN_ejecta_rate')
      stop ! stop the program
    end if
    !
    f_agn = 2.d0*disc_ejecta_efficiency*AGN_mass_energy_conv*AGN_kinetic_fraction/(agn_ejecta_efficiency + 1.d0)  ! without unit
    !
    if (f_agn .gt. 0.d0) then
        !
        agn_compute_agn_ejecta_rate = f_agn*(light_speed_code_unit/Vesc)**2.*agn_compute_agn_infall_rate(agn)
    end if

    return 
  end function agn_compute_agn_ejecta_rate

  !*****************************************************************************************************************

  function agn_compute_agn_wind_velocity(agn,Vesc)

    ! RETURN THE EFFECTIVE VELOCITY OF THE WIND PRODUCED BY THE AGN
    ! To maximize the ejection process the wind average velocity is set to the escape velocity of the galaxy 
    
    implicit none

    real(kind=8),intent(in)      :: Vesc      ! galaxy escape velocity, must be given in code unit
    real(kind=8)                 :: agn_compute_agn_wind_velocity

    type(agn_type),intent(in)    :: agn       ! an agn component

    agn_compute_agn_wind_velocity = 0.d0      ! init
    
    if (agn_mass(agn) .le. 0.d0) return       ! no black hole and therefore no accretion process

    agn_compute_agn_wind_velocity = Vesc      ! in code unit
    
    return
  end function agn_compute_agn_wind_velocity
  
  !*****************************************************************************************************************

  function agn_compute_agn_turbulent_heating_power(agn)

    ! COMPUTE THE INSTANTANEOUS KINETIC POWER PRODUCED BY THE AGN AND TRANSFERED TO TURBULENCE HEATING/MOTION
    
    implicit none 
    
    real(kind=8)                 :: agn_compute_agn_turbulent_heating_power
    real(kind=8)                 :: acc_rate       ! real accretion rate onto the SMBH
    
    type(agn_type),intent(in)    :: agn            ! an agn component
    
    agn_compute_AGN_turbulent_heating_power = 0.d0 ! init
    
    if (agn_mass(agn) .le. 0.d0) return            ! no black hole and therefore no accretion process

    acc_rate = agn_compute_agn_accretion_rate(agn) ! real infall rate onto the SMBH
    
    if (acc_rate .gt. 0.d0) then
      ! in code unit (10^11.Msun.kpc^2/Gyr^3)
      agn_compute_agn_turbulent_heating_power = (1.d0 - disc_ejecta_efficiency)*AGN_kinetic_fraction &
                                                        *AGN_mass_energy_conv*light_speed_code_unit**2.*acc_rate 
    end if

    return
  end function agn_compute_agn_turbulent_heating_power

  !*****************************************************************************************************************

  function agn_compute_agn_non_kinetic_power(agn)

    ! COMPUTE THE INSTANTANEOUS NON KINETIC POWER PRODUCED BY THE AGN
    
    implicit none 
    
    real(kind=8)                 :: agn_compute_agn_non_kinetic_power
    real(kind=8)                 :: acc_rate       ! real accretion rate onto the SMBH
    
    type(agn_type),intent(in)    :: agn            ! an agn component
    
    agn_compute_agn_non_kinetic_power = 0.d0       ! init
    
    if (agn_mass(agn) .le. 0.d0) return            ! no black hole and therefore no accretion process

    acc_rate = agn_compute_agn_accretion_rate(agn) ! real infall rate onto the SMBH
    
    if (acc_rate .gt. 0.d0) then
      ! in code unit (10^11.Msun.kpc^2/Gyr^3)
      agn_compute_agn_non_kinetic_power = (1.d0 - AGN_kinetic_fraction)*AGN_mass_energy_conv*light_speed_code_unit**2.*acc_rate 
    end if

    return
  end function agn_compute_agn_non_kinetic_power

  !*****************************************************************************************************************

  function agn_compute_agn_thermal_power(agn)

    ! RETURN THE INSTANTANEOUS THERMAL POWER PRODUCED BY THE AGN

    implicit none

    real(kind=8)                 :: agn_compute_agn_thermal_power

    type(agn_type),intent(in)    :: agn       ! an agn component

    agn_compute_agn_thermal_power = AGN_thermal_fraction*agn_compute_agn_non_kinetic_power(agn)

    return
  end function agn_compute_agn_thermal_power

  !*****************************************************************************************************************

  function agn_compute_agn_rad_power(agn)

    ! RETURN THE INSTANTANEOUS RADIATIVE POWER (NON THERMAL AND NON KINETIC POWER) PRODUCED BY THE SMBH ACTIVITY
    
    implicit none
    
    real(kind=8)                 :: agn_compute_agn_rad_power

    type(agn_type),intent(in)    :: agn       ! an agn component

    agn_compute_agn_rad_power = (1.d0 - AGN_thermal_fraction)*agn_compute_agn_non_kinetic_power(agn)

    return
  end function agn_compute_agn_rad_power
  
  !*****************************************************************************************************************
    
  function agn_compute_L(agn)
      
    ! COMPUTE AGN INSTANTANEOUS LUMINOSITY
      
    implicit none
      
    real(kind=8)                 :: agn_compute_L ! The instantaneous luminosity (code unit)

    type(agn_type),intent(in)    :: agn           ! the agn COMPONENT

    agn_compute_L = agn_compute_agn_rad_power(agn)
                                                                                                 
    return
  end function agn_compute_L

  !*****************************************************************************************************************
    
  function agn_compute_L_edd(agn)
      
    ! COMPUTE AGN EDDINGTON LUMINOSITY
      
    implicit none
      
    real(kind=8)                 :: agn_compute_L_edd ! The Eddington luminosity (code unit)

    type(agn_type),intent(in)    :: agn               ! the agn COMPONENT

    agn_compute_L_edd = 4.d0*pi*gravconst_code_unit*agn_mass(agn)*light_speed_code_unit*1.20387d5   ! *1.20387.e5 is (mp/mass_code_unit_in_kg)/(sigma_t/Mpc_in_m**2.)
                                                                                                    ! in code unit
    return
  end function agn_compute_L_edd

  !*****************************************************************************************************************
  !
  ! PRINTING PROCEDURES
  !
  !*****************************************************************************************************************
    
  subroutine agn_load_agn_data(agn,agn_data)

    implicit none

    real(kind=8),intent(inout)  :: agn_data(nb_agn_field)

    type(agn_type),intent(in)   :: agn  ! a agn component
      
    ! for information
    ! 'agn_mass              ','agn_L                 ','agn_L_Edd             ', 'agn_acc_rate_in       '     
    agn_data = (/agn%mass*1.d11, agn%L/L_Sun_in_code_unit, agn%L_edd/L_Sun_in_code_unit,gas_mass(agn%acc_rate)*1.d2/)

    return
  end subroutine agn_load_agn_data
  
  !*****************************************************************************************************************

  end module agn
