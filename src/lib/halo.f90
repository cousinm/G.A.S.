module halo
  
  use baryon_halo ! Contains bh structure and bh exploitation functions
  use galaxy      ! Contains gal structure and gal exploitation functions
  
  public
  
  !*****************************************************************************************************************
  ! 
  !  OVERVIEW
  !
  ! halo_module defines the halo data structures univ(ts)%halo 
  ! This structure contains the entire information about baryonic phase of a halo
  ! The halo strcuture alloxs to link the galaxy component and the baryon halo phase (cold and hot accretion mode)
  ! It also allows to link baryonic phases and the dark-matter component
  ! Constraints on the evolution schemes of baryonic phases (galaxy and baryon halo) are computed at the halo scale 
  !
  !  MPI PROCESSES IN THIS MODULE
  !
  !  halo_send_data                                               : send specific informations about a halo
  !      called by : tree_evolve
  !
  !  halo_receive_data                                            : receive specific informations about a halo
  !      called by : halo_spectrum
  !
  !  SUBROUTINES IN THIS MODULE
  !
  !  halo_set_reference_mass                                      : initialize reference mass
  !      called_by : tree_set_reference_mass
  !
  !  halo_void                                                    : initialize a halo structure
  !  
  !  halo_copy                                                    : copy the baryon_halo and the galaxy component of a halo component in an other
  !  
  !  halo_save                                                    : save the dm, the baryon_halo and the galaxy component of a halo component in a other
  !                                                                 used by tree_evolve to create a complete copy of a halo component  
  ! 
  !  halo_deallocate                                              : erase (deallocate) a halo component
  !     called by : tree_deallocate  
  !
  !  halo_set_halo_evolved                                        : set halo%evolved = .true.
  !     called by : tree_evolve (for progenitor @ ts -1)
  !                 halo_evolve
  !
  !  halo_set_properties                                          : set properties of a halo structure
  !     called by : tree_compute_halo_derived_properties                            
  !
  !  halo_evolve                                                  : compute constraints on the evolution scheme of galaxy and baryon halo components
  !     called by : tree_evolve          
  !
  !  halo_merge                                                   : merge two halo components
  !     called by : tree_evolve           
  !
  !  halo_spectrum                                                : call galaxy scale subroutine to compute galaxy spectra and magnitudes
  !     called by : tree_spectrum    
  !
  !  halo_compute_stripping_rate                                  : Compute the effective instantaneous striping gas rate
  !                                                                    affecting the diffuse disc gas component of a galaxy that 
  !                                                                    progresses into the hot has phase of a host massive halo   
  ! 
  !  FUNCTIONS IN THIS MODULE
  ! 
  !  halo_baryonic_mass                                           : return the baryonic mass of a baryonic halo component
  !
  !  halo_escape_velocity                                         : return the escape velocity of the halo by taking into account dark-matter, baryon halo and galaxy mass
  !
  !  PRINTING PROCEDURES
  !
  !  halo_print                                                   : calls print routines dedicated to all halo components
  !    called by  : tree_print 
  !
  !*****************************************************************************************************************
  
  ! DEFINITION OF GLOBAL VARIABLES LINKED TO THE HALO STRUCTURE *******************

  ! HALO_TYPE DEFINITION **************

  type halo_type                                     
     ! cleaning
     logical                :: compute               ! = .true. if the halo is realy defined
     logical                :: evolved               ! = .true. if the halo has been evolved (baryon)
     logical                :: followed              ! = .true. if the halo is followed
     logical                :: recursively_kill      ! set to .true. if the branch is already kill by the recursive routine (we win some time)
     real(kind=8)           :: goodness_of_branch    ! nb of (halo-)good halos / nb of halos on branch
     ! structure of the halo
     type(dm_type)          :: dm                    ! a dark matter halo
     type(baryon_halo_type) :: baryon_halo           ! a baryonic halo
     type(galaxy_type)      :: galaxy                ! a galaxy
  end type halo_type
  
  ! OTHER DEFINITIONS *******************

  real(kind=8)              :: halo_dt_min_use       ! minimal time-step used at the halo scale

contains

  !*****************************************************************************************************************
  ! 
  ! MPI PROCESSES
  !
  !*****************************************************************************************************************

#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------  
  subroutine halo_send_data(h,z,post_merger)                                     
  
    ! SEND SPECIFIC INFORMATION ABOUT A HALO
  
    implicit none

    logical,intent(in)          :: post_merger
    
    real(kind=8),intent(in)     :: z  ! redshift
    
    type(halo_type),intent(in)  :: h  ! a halo
    
    ! send current redshift
    call MPI_SEND(z,1,MPI_REAL8,rank+1,halo_tag+1,MPI_COMM_WORLD,ierror)
    ! send post merger flag
    call MPI_SEND(post_merger,1,MPI_LOGICAL,rank+1,halo_tag+2,MPI_COMM_WORLD,ierror)
        
    ! go to galaxy scale
    call galaxy_send_data(h%galaxy) 
    
    return
  end subroutine halo_send_data 
  
  !*****************************************************************************************************************
  
  subroutine halo_receive_data(h,z,post_merger)                                     
  
    ! RECEIVE SPECIFIC INFORMATION ABOUT A HALO
  
    implicit none
    
    logical,intent(out)            :: post_merger
    
    real(kind=8),intent(out)       :: z  ! redshift
    
    type(halo_type),intent(inout)  :: h  ! a halo
    
    ! receive redshift
    call MPI_RECV(z,1,MPI_REAL8,rank-1,halo_tag+1,MPI_COMM_WORLD,statut,ierror)
    ! receive post merger flag
    call MPI_RECV(post_merger,1,MPI_LOGICAL,rank-1,halo_tag+2,MPI_COMM_WORLD,statut,ierror)
    
    ! receive properties of the galaxy
    call galaxy_receive_data(h%galaxy)   
  
    return
  end subroutine halo_receive_data
  
! ------------------------------------------------- 
#endif
! LUMINOUS_PROCESSES 
  
  !*****************************************************************************************************************
  !
  ! SUBROUTINES
  !
  !*****************************************************************************************************************

  subroutine halo_set_reference_mass

    ! INITIALIZE REFERENCE PARAMETERS LINKED TO THE DARK-MATTER N-BODY SIMULATION

    implicit none

    call dm_set_reference_mass
    call baryon_halo_set_reference_mass 
    call galaxy_set_reference_mass
    
    return
  end subroutine halo_set_reference_mass

  !*****************************************************************************************************************

  subroutine halo_void(h)

    ! INIT OR VOID A HALO COMPONENT

    implicit none 

    type(halo_type),intent(inout) :: h     ! a halo component

    ! Evolution
    h%evolved             = .false.        ! initially the halo has not been evolved.
    ! follow
    h%followed            = .false.        ! by default a halo is not followed
    ! cleaning
    h%compute             = .false.        ! initially, the halo is not computable. The subroutine tree-clean analyses the branches
    h%recursively_kill    = .false.        ! used and setted by tree_kill
    h%goodness_of_branch  = -1.   
    ! halo contains
    call dm_void(h%dm)                     ! void the dm component
    call baryon_halo_void(h%baryon_halo)   ! void the baryon_halo component
    call galaxy_void(h%galaxy)             ! void the galaxy component
  
    return
  end subroutine halo_void

  !*****************************************************************************************************************

  subroutine halo_copy(h1,h2)

    ! COPY THE BARYON AND THE GALAXY CONTENT (NOT THE DM) OF A HALO h2 IN AN OTHER HALO h1

    implicit none
    
    type(halo_type),intent(inout) :: h1    ! the original halo component
    type(halo_type),intent(in)    :: h2    ! the copie

    ! Evolution
    h1%evolved            = h2%evolved
    ! follow
    h1%followed           = h2%followed
    ! cleaning
    h1%compute            = h2%compute
    h1%recursively_kill   = h2%recursively_kill
    h1%goodness_of_branch = h2%goodness_of_branch 
    ! 
    call baryon_halo_copy(h1%baryon_halo,h2%baryon_halo)     ! Copy the baryon content
    call galaxy_copy(h1%galaxy,h2%galaxy)                    ! Copy the galaxy content  
    
    return
  end subroutine halo_copy

  !*****************************************************************************************************************

  subroutine halo_save(h1,h2)

    ! COPY THE DM, THE BARYON AND THE GALAXY CONTENT OF THE MAIN PROGENITOR h2 INTO ITS DESCENDENT h1
    ! Thi ssubroutine is used to create a intermediate halo in tree_evolve

    implicit none

    type(halo_type),intent(inout)  :: h1   ! the original halo component
    type(halo_type),intent(in)     :: h2   ! the copy

    call dm_copy(h1%dm,h2%dm)              ! copy the dark matter content
    call halo_copy(h1,h2)                  ! copy the halo content

    return
  end subroutine halo_save
  
  !*****************************************************************************************************************

  subroutine halo_deallocate(h)

    ! DEALLOCATE THE HALO VARIABLE (halo%dm doesn't contain allocatable fields)

    implicit none 

    type(halo_type),intent(inout) :: h

    call galaxy_deallocate(h%galaxy)

    return
  end subroutine halo_deallocate

    
  !*****************************************************************************************************************
  
  subroutine halo_set_halo_evolved(h)

    implicit none
    
    type(halo_type),intent(inout) :: h
    
    h%evolved = .true.
    
    return
  end subroutine 
    
  !*****************************************************************************************************************
  
  subroutine halo_set_properties(h,ts,props)

    ! LOAD INTO THE HALO THE PROPERTIES DIRECTLY MEASURED FROM THE N-BOBY SIMULATION

    implicit none

    integer(kind=4),intent(in)    :: ts         ! timestep at which the halo has been identified

    real(kind=4),intent(in)       :: props(:)   ! temporary array in which the properties of dm haloes are read
    real(kind=4)                  :: Ek_Ep      ! ratio of kinetic to potential energies

    type(halo_type),intent(inout) :: h          ! a halo component

#ifdef HORIZON
! -------------------------------------------------
    ! *************************!
    ! props(1:3)   = x,y,z     ! 
    ! props(4:6)   = Vx,Vy,Vz  !
    ! props(7)     = M_tot     !
    ! props(8)     = r         ! 
    ! props(9)     = Spin      !
    ! props(10)    = R_vir     !
    ! props(11)    = M_vir     !
    ! props(12)    = T_vir     !
    ! props(13)    = cvel      !
    ! props(14)    = dmacc     !
    ! props(15)    = frag      !
    ! props(16:18) = Lx,Ly,Lz  !
    ! props(19)    = Ep        !
    ! props(20)    = Ek        !
    ! props(21)    = Et        ! 
    ! *************************!

    if (abs(props(19)) .gt. 0.d0) then
        Ek_Ep = props(20)/abs(props(19))
    else
        Ek_Ep = -1.0
    end if
    ! dm_set_properties(dm,ts,x,y,z,Vx,Vy,Vz,M_vir,M_tot,R_vir,spin,dmacc,Lx,Ly,Lz,Ek_Ep)
    call dm_set_properties(h%dm,ts,real(props(1),8),real(props(2),8),real(props(3),8), &
                           real(props(4),8),real(props(5),8),real(props(6),8),real(props(11),8), &
                           real(props(7),8), real(props(10),8), real(props(9),8), &
                           real(props(14),8),real(props(16),8), real(props(17),8), real(props(18),8), real(Ek_Ep,8))
#endif
! -------------------------------------------------
! HORIZON

#ifdef BOLSHOI
! -------------------------------------------------

    ! *************************!
    ! props(1:3)   = x,y,z     ! 
    ! props(4:6)   = Vx,Vy,Vz  ! 
    ! props(7)     = M_tot     ! 
    ! props(8)     = Spin      ! 
    ! props(9)     = R_vir     ! 
    ! props(10)    = M_vir     ! 
    ! props(11)    = dmacc     ! 
    ! props(12:14) = Lx,Ly,Lz  !
    ! props(15)    = Ep_Ek     !
    ! *************************!
    
    ! dm_set_properties(dm,ts,x,y,z,Vx,Vy,Vz,M_vir,M_tot,R_vir,spin,dmacc,Lx,Ly,Lz,Ek_Ep)
    call dm_set_properties(h%dm,ts,real(props(1),8),real(props(2),8),real(props(3),8), &
                           real(props(4),8),real(props(5),8),real(props(6),8),real(props(10),8), &
                           real(props(7),8), real(props(9),8), real(props(8),8), &
                           real(props(11),8),real(props(12),8), real(props(13),8), real(props(14),8), real(props(15),8))

! -------------------------------------------------
#endif
! BOLSHOI

    return 
  end subroutine halo_set_properties

  !*****************************************************************************************************************
  
  subroutine halo_evolve(h,z,dt,post_merger,host)
    
    ! MAKE EVOLVED THE BARYON AND THE GALXY CONTENT OF THE HALO DURING dt 

    implicit none
    
    logical,intent(in)                     :: post_merger               ! = .true. if the evolution is done onto the remnent structure of a merger

    real(kind=8),intent(in)                :: z                         ! until what redshift the halo must grow 
    real(kind=8),intent(in)                :: dt                        ! time step (in Gyr)
    real(kind=8)                           :: ddt 
    real(kind=8)                           :: time_before_end           ! time until the end of the evolution time-step   

    type(halo_type),intent(inout)          :: h                         ! the halo component
    type(halo_type),intent(inout),optional :: host                      ! the host halo (transfer process or hot ejected mass transfert)
    !
    ! cold phase component
    real(kind=8)                           :: cold_time                 ! evolution time of the cold halo phase
    real(kind=8)                           :: cold_dt_optim             ! optimal time-step of the cold halo phase
    real(kind=8)                           :: cold_next_stop            ! next evolution stop for the cold halo phase
    !
    ! hot phase component
    real(kind=8)                           :: hot_time                  ! evolution time of the hot halo phase
    real(kind=8)                           :: hot_dt_optim              ! optimal time-step of the hot halo phase
    real(kind=8)                           :: hot_next_stop             ! next evolution stop for hot halo phase
    !
    ! galaxy component
    logical                                :: stop_before_end
    real(kind=8)                           :: gal_time                  ! evolution time of the galaxy component (follow the time evolution of the hot halo phase)
    real(kind=8)                           :: gal_dt_optim              ! max value of ddt whitch respect the quasi-static criterion in all components of the halo
    real(kind=8)                           :: Vwind                     ! galaxy wind velocity 
    real(kind=8)                           :: inst_galaxy_stripping_rate! instantaneous value of the stripping rate
    real(kind=8)                           :: inst_galaxy_ejecta_rate   ! the pre galaxy-evolution instantaneous galaxy_ejecta rate used to compute WARNING barriers
    real(kind=8)                           :: f_in                      ! fraction of galaxy ejecta that is catched by the hot atmosphere
    real(kind=8)                           :: f_in_Wd                   ! lower possible value for f_in
    real(kind=8)                           :: f_in_Wu                   ! upper possible value for f_in

    type(gas_type)                         :: galaxy_stripping_rate     ! everage effective stripping rate affecting the diffuse gas contained in the disc
    type(gas_type)                         :: galaxy_ejecta_rate_Wd     ! ejecta barrier (down)
    type(gas_type)                         :: galaxy_ejecta_rate_Wu     ! ejecta barrier (up)
    type(gas_type)                         :: galaxy_ejecta_rate        ! the time(-average) galaxy ejecta rate computed throught ddt and used to evolve the hot gas component
    type(gas_type)                         :: galaxy_fresh_gas_acc_rate ! the net gas accretion rate onto the galaxy
         
    if (dt .le. 0.d0) return 
             
#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('halo_evolve',only_rank=rank,component='halo')
! -------------------------------------------------
#endif

    if (h%dm%life_time .eq. 0.d0) then
        if (FOLLOW_UP .and. PR_FOLLOW_UP) then 
            ! print dark-matter properties at the first step
            ! dm_print(unit,form,dm)
            call dm_print(follow_up_unit(current_index),'fits',h%dm)
        end if
    end if

    ! initialize all gas component used in the procedure
    call gas_void(galaxy_ejecta_rate)        ! init the time(-average) galaxy ejeta rate 
    call gas_void(galaxy_fresh_gas_acc_rate) ! init the net galaxy gas accretion rate

    ! initialize local time-step computation variables
    time_before_end = dt        ! init, time before the end of the time step
    cold_time       = 0.d0      ! init the cold clock
    hot_time        = 0.d0      ! init the hot clock
    cold_dt_optim   = 0.d0      ! init the optimal time-step time of the cold component
    hot_dt_optim    = 0.d0      ! init the optimal time-step time of the hot component
    cold_next_stop  = 0.d0      ! init the next STOP of the cold clock
    hot_next_stop   = 0.d0      ! init the next STOP of the hot clock
    gal_time        = 0.d0      ! init the galaxy clock
    gal_dt_optim    = 0.d0      ! init the optima time-step of the galaxy component
    !
    stop_before_end = .false.
    !
    ! set inflow accretion rate for the baryon halo component
    ! not in the while loop because the halo_cold_inflow_rate and the halo_hot_inflow_rate are constant during 
    ! all the timestep evolution of the halo
    ! baryon_halo_compute_accretion_rate(bh,dm,z,dt,post_merger,cold_inflow_rate,hot_inflow_rate)
    call baryon_halo_compute_accretion_rate(h%baryon_halo,h%dm,z,dt,post_merger)
    !
    ! update hot halo density profile properties
    ! we must update here beacause dm%r_core and dm%rho_core may change between two dm-fixed time-stepes and 
    ! we must taking into account these modifications in the cooling process computation
    call baryon_halo_update_density_profile(h%baryon_halo,h%dm,called_by='halo_evolve_dm_update')
    !
    do while (abs(time_before_end) .gt. num_accuracy)
      !
      ! *********************************************************************************************************
      ! COMPUTE RATE (for cold and hot component of the baryon halo) (PREDICTOR PART) 
      ! *********************************************************************************************************
      !
      ! ****************************************************************
      !
      !
      ! COLD PHASE 
      ! filamentary structure
      !
      !
      ! ****************************************************************
      !
      if (abs(cold_time - cold_next_stop) .lt. num_accuracy) then
        !
        ! compute (and set) input and output rates of the 'bh' cold gas reservoirs
        call baryon_halo_evolve_cold_gas_I(h%baryon_halo,h%dm,cold_dt_optim)
        !
        ! compute the next cold STOP 
        ! cold time does not depend of hot time because galaxy does not affect the cold phase
        cold_next_stop = min(dt,cold_time + cold_dt_optim)    
      end if
      !
      ! ****************************************************************
      !
      !
      ! HOT PHASE
      ! hot gas in hydrostatic equilibrium in the dark-matter potential well
      !
      !
      ! ****************************************************************
      !
      call gas_void(galaxy_ejecta_rate_Wd) ! init
      call gas_void(galaxy_ejecta_rate_Wu) ! init
      !
      if (abs(hot_time - hot_next_stop) .lt. num_accuracy) then
        !       
        ! galaxy_compute_galaxy_feedback_activities(gal,dm,[ejecta_rate,agn_acc_rate,Vwind,Qtherm,Qturb,Qrad])
        call galaxy_compute_galaxy_feedback_activities(h%galaxy,h%dm,ejecta_rate=inst_galaxy_ejecta_rate,f_in=f_in)
        ! convert inst_galaxy_ejecta_rate (real 8) into a gas object
        ! we assume that ejecta have a similar composition than the galaxy disc
        call gas_void(galaxy_ejecta_rate)      ! init
        galaxy_ejecta_rate = inst_galaxy_ejecta_rate*galaxy_gas_signature(h%galaxy,component='disc',apply_as='rate_builder',called_by='halo_evolve')
        !
        call baryon_halo_evolve_hot_gas_I(h%baryon_halo,h%dm,galaxy_ejecta_rate,f_in, &
                    hot_dt_optim,galaxy_ejecta_rate_Wd,galaxy_ejecta_rate_Wu,f_in_Wd,f_in_Wu)                   ! init
        !
        ! compute hot next STOP
        hot_next_stop = min(dt,min(hot_time + hot_dt_optim,cold_next_stop))  ! hot next mark cannot be greater than cold next mark
        !
      end if ! if hot_time
      !   
      ! Reset galaxy ejecta rate
      call gas_void(galaxy_ejecta_rate) ! init
      Vwind = 0.d0   
      ! 
      ! ****************************************************************
      !
      !
      ! GALAXY
      ! Evolve the galaxy component during ddt (if possible) or STOP after gal_dt_optim
      !
      !
      ! ****************************************************************
      !
      ! Compute optimal time-step
      ! the galaxy time is set to the hot gas phase. Indeed, hot gas and galaxy components are closely linked (galaxy winds)  
      ! compute the new galaxy ddt
      ddt = hot_next_stop-gal_time
      ! 
      if ((ddt .lt. num_precision) .and. (ddt .lt. time_before_end)) then        
          call IO_print_warning_message('Very small timestep',only_rank=rank,called_by='halo_evolve')  
      end if           
      if (ddt .le. num_accuracy) then
        call IO_print_error_message('ddt =< num_precision',only_rank=rank,called_by='halo_evolve')  
        stop
      end if ! ddt < num_accuracy
      ! Compute galaxy fresh gas accretion rate
      ! The inflow of the galaxy are the outflow of the baryon halo
      call baryon_halo_compute_galaxy_fresh_gas_accretion_rate(h%baryon_halo,galaxy_fresh_gas_acc_rate)
      !
      ! Compute galaxy stripping 
      inst_galaxy_stripping_rate = 0.d0
#ifdef SUB_STRIPPING
! -------------------------------------------------       
      call halo_compute_stripping_rate(h,inst_galaxy_stripping_rate,host=host)
! -------------------------------------------------
#endif  
! SUB_STRIPPING      
      ! convert inst_galaxy_ejecta_rate (real 8) into a gas object
      ! we assume that the stipped gas has a similar composition than the diffuse gas contained in the disc
      galaxy_stripping_rate = inst_galaxy_stripping_rate*galaxy_gas_signature(h%galaxy,component='disc',subcomponent='diffuse',apply_as='rate_builder',called_by='galaxy_stripping_rate')
      !
      ! compute galaxy evolution during ddt = gal_next_stop-gal_time
      if ((galaxy_mass(h%galaxy) .gt. 0.d0) .or. gas_mass(galaxy_fresh_gas_acc_rate) .gt. 0.d0) then   
        !      
        ! galaxy_evolve(gal,dm,z,dt,fresh_gas_acc_rate,gal_stripping_rate,gal_ejecta_rate_Wd,gal_ejecta_rate_Wu,Qturb_Wd,gal_ejecta_rate,Vwind,dt_optim,stop_before_end)
        call galaxy_evolve(h%galaxy,h%dm,ddt,galaxy_fresh_gas_acc_rate,galaxy_stripping_rate, &
                               galaxy_ejecta_rate_Wd,galaxy_ejecta_rate_Wu,f_in_Wd,f_in_Wu,galaxy_ejecta_rate, &
                               Vwind,gal_dt_optim,stop_before_end) 
      else
        gal_dt_optim    = ddt
        stop_before_end = .false.
      end if
      !   
      ! compute real gal_next_stop
      ! if the galaxy can't be evolve over all the ddt time (ejecta effect) so gal_dt_optim is smaller than ddt
      ! take account this stop and update the galaxy clock
      gal_time = gal_time + gal_dt_optim
      !
      if ((gal_time .gt. hot_next_stop) .and. (abs(gal_time - hot_next_stop) .gt. num_accuracy)) then
        call IO_print_error_message('gal_time > hot_next_stop',only_rank=rank,called_by='halo_evolve')
        call IO_print_message('use',only_rank=rank,component='halo', &
                param_name = (/'gal_time                 ','hot_next_stop            '/), &
                real_param_val  = (/gal_time,hot_next_stop/))
        stop  ! stop the program
      end if
      hot_next_stop = gal_time
      ddt           = gal_dt_optim
      !
      ! update life_time of the dm component
      call dm_update_life_time(h%dm,ddt)
      !
      ! *****************************************************************************************
      !
      !
      ! UPDATE BARYON HALO PROPERTIES (CORRECTOR PART) 
      !
      !
      ! *****************************************************************************************
      !
      ! ******************************************
      !
      ! HOT PHASE 
      !
      ! ******************************************
      !
      ! baryon_halo_evolve_hot_gas(bh,dm,dt,galaxy_ejecta_rate,Vwind)
      call baryon_halo_evolve_hot_gas_II(h%baryon_halo,h%dm,ddt,galaxy_ejecta_rate,Vwind)
      !
      if (present(host)) then
         if (gas_mass((galaxy_stripping_rate+galaxy_ejecta_rate)*ddt) .gt. 0.d0) then
           ! The halo is a sub structured with a "computable" main host halo
           ! Add the ejected gas to the main host halo
           call baryon_halo_transfer_hot_mass(host%baryon_halo,galaxy_ejecta_rate*ddt)
           ! Add the stripped gas to the main host halo
           call baryon_halo_transfer_hot_mass(host%baryon_halo,galaxy_stripping_rate*ddt)
           ! The hot atmosphere of the host halo has been modified
           ! New properties of this hot atmosphere has to be taken in to account at the next ddt (stripping)
           call baryon_halo_update_density_profile(host%baryon_halo,host%dm,called_by='halo_evolve_after_stripping') 
         end if
      else
#ifdef REACCRETION
! ------------------------------------------------- 
          ! Add escapted gas to the surrounding gas phase, will be progressivelly reaccreted by the halo 
          call baryon_halo_add_surrounding_gas(h%baryon_halo,galaxy_ejecta_rate*ddt)
! ---------
#else
! ---------
          ! feed the IGM, the gas is definitivelly lost ! 
          call gas_add_igm_gas(galaxy_ejecta_rate*ddt,called_by='halo_evolve / m_esc --> igm')
! ------------------------------------------------- 
#endif
! REACCRETION         
      end if
      !   
      ! update hot_time
      hot_time = gal_time 
      !
      ! ******************************************
      !
      ! COLD PHASE 
      ! cold phase doesn't depend of the galaxy ejecta. We have to update its properties just when cold_next_stop = gal_time
      ! for machine precision problem we don't ask a strong equality between gal_time and cold_time 
      ! we accept a different of num_accuracy betwenn the two times
      ! check variation of the cold halo phase, crash code if the variation is greater than the physical asked precision
      !
      ! ******************************************  
      !
      if (abs(cold_next_stop - gal_time) .lt. num_accuracy) then
         !
         call baryon_halo_evolve_cold_gas_II(h%baryon_halo,(cold_next_stop-cold_time))
         !
         ! update cold_next_stop
         cold_next_stop = gal_time
         ! update cold_time   
         cold_time = cold_next_stop
         !
         time_before_end = max(0.d0, dt-gal_time)  
      end if                                             
    end do ! end time loop 
    !  
    ! check time evolution of the hot, the cold phase and of the galaxy 
    ! crash code if the difference is greater than num_accuracy
    if (abs(cold_time - gal_time) .gt. num_accuracy) then
      call IO_print_error_message('cold_time != gal_time ',only_rank=rank,called_by='halo_evolve')
      call IO_print_message('with',only_rank=rank,component='halo', &
                param_name = (/'gal_time                 ','cold_time                '/), &
                real_param_val  = (/gal_time,cold_time/))
      stop  ! stop the program
    end if
    if (abs(hot_time - gal_time) .gt. num_accuracy) then
      call IO_print_error_message('hot_time != gal_time ',only_rank=rank,called_by='halo_evolve')
      call IO_print_message('with',only_rank=rank,component='halo', &
                param_name = (/'gal_time                 ','hot_time                 '/), &
                real_param_val  = (/gal_time,hot_time/))
      stop
    end if
    !
    if (FOLLOW_UP .and. PR_FOLLOW_UP) call dm_print(follow_up_unit(current_index),'fits',h%dm)
                
#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('halo_evolve ... done',only_rank=rank,component='halo')
! -------------------------------------------------
#endif
    !   
    return   
  end subroutine halo_evolve
  
  !*****************************************************************************************************************
 
  subroutine halo_merge(h1,h2) 

    ! MERGE THE SMALLER PROGENITOR h2 WITH THE DESCENDENT h1 OF THE MAIN PROGENITOR 
  
    implicit none 

    type(halo_type),intent(inout) :: h1                 ! the descendent halo
    type(halo_type),intent(in)    :: h2                 ! the progenitor
    

#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('halo_merge',only_rank=rank,component='halo')
! -------------------------------------------------
#endif
    !
    ! set (dm) age of the formation and life time 
    h1%dm%age_form  = min(h1%dm%age_form,h2%dm%age_form)
    h1%dm%life_time = max(h1%dm%life_time,h2%dm%life_time)
    if (FOLLOW_UP .and. PR_FOLLOW_UP) then 
       ! the progenitor h2, used here, is a followed structure  
       ! we have to erase properties associated to h1 and replace them by h2 properties
       h1%dm%life_time = h2%dm%life_time
       h1%dm%age_form  = h2%dm%age_form
    end if
    !
    ! 1) Merge the baryonic component (the hot and the cold reservoir)
    call baryon_halo_merge(h1%baryon_halo,h2%baryon_halo,h1%dm)    
    !
    ! 2) Merge the galaxy component
    call galaxy_merge(h1%galaxy,h2%galaxy,h1%dm,h2%dm) 

#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('halo_merge ... done',only_rank=rank,component='halo')
! -------------------------------------------------
#endif

    return
  end subroutine halo_merge
 
#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------   
    
  !*****************************************************************************************************************
  
  subroutine halo_spectrum(h)
  
    implicit none
     
    logical                        :: post_merger
    
    real(kind=8)                   :: z          ! the current redshift
    
    type(halo_type),intent(out)    :: h          ! a halo
    
    call halo_void(h) 
    
    call halo_receive_data(h,z,post_merger)
    !
    ! build the full galaxy spectrum and compute magnitudes associated to each input filters
    call galaxy_spectrum(h%galaxy,z)   
  
    return
  end subroutine halo_spectrum
   
! ------------------------------------------------- 
#endif
! LUMINOUS_PROCESSES 
 
  !*****************************************************************************************************************
  !
  ! FUNCTIONS
  !
  !*****************************************************************************************************************

#ifdef SUB_STRIPPING
! -------------------------------------------------   
  subroutine halo_compute_stripping_rate(h,inst_stripping_rate,host)
  
    ! COMPUTE THE AVERAGE STRIPPING RATE AFFECTING THE SUB-HALO h DURING ITS EVOLUTION IN THE 
    ! HOT GAS PHASE ASSOCIATED WITH THE MAIN HALO host
    ! The stripping efficiency depends on both: the position of the satellite into the hot gas phase according to the density profile.
    !                                           the orientation of the galaxy into the hot gas according to the scalar product dV.L
    ! The stripping rate is proportionnal to the diffuse gas mass. 
    ! The stripping rate is assumed to be constant during a hot gas phase time-step. 
    ! Concequently, we use the mass of the diffuse gas at the begining of the hot gas phase time-step.
    ! The gas coming from the sub-structures is immediatly transfered to the hot gas reservoir of the main host halo 
    ! and its temperature is fixed to the temperature of the host hot halo phase
  
    implicit none
    
    real(kind=8)                        :: Dist                ! the distance between the two halos
    real(kind=8)                        :: dV(3)               ! Differential velocity
    real(kind=8)                        :: L(3)                ! Angular momentum vector
    real(kind=8)                        :: rho_hot, rho_diff   ! gas density
    real(kind=8)                        :: drho, av_rho        ! density contrast   
    real(kind=8)                        :: M_unstr             ! unstructured/diffuse gas mass contained in the disc    
    real(kind=8)                        :: t_fric              ! friction characteristical time
    real(kind=8),intent(out)            :: inst_stripping_rate ! stripping rate
    
    type(halo_type),intent(in)          :: h                   ! the halo containing the stripped galaxy
    type(halo_type),intent(in),optional :: host                ! the host halo containing the hot gas that strip the galaxy

    ! init 
    inst_stripping_rate = 0.d0

    if (present(host)) then
        if (baryon_halo_bh_mass(host%baryon_halo,component='hot') .le. 0.d0) return ! no hot gas in the host halo
        
        ! compute the distance between the two halos
        Dist = dm_orbital_radius(h%dm, host%dm)
        
        if (Dist .gt. host%dm%R_vir) return ! the galaxy is out of the hot gas halo --> No stripping
        
        ! The stripping efficiency is depending on various parameters
        ! 1- The difference between the velocity vector of the host halo and the veloity vector of the galaxy's halo
        dV      = (h%dm%vel - host%dm%vel)
        dV      = dV/sqrt(sum(dV**2.))
        ! 2- The density contrast between the diffuse gas contained in the satellite galaxy and the hot gas of the main halo
        rho_hot  = baryon_halo_hot_halo_density_profile(Dist,host%baryon_halo,host%dm) ! [code unit: 10^11Msun/kpc^3]  
        M_unstr  = disc_mass(h%galaxy%disc,component='unstructured') 
        if ((rho_hot .gt. 0.d0) .and. (M_unstr .gt. M_gas_crit)) then
            ! the diffuse gas evolve in a thick disc 2 times larger than the stellar component 
            ! and with a disc scale height "h"
            ! 99.9% of the mass of an exponential disc with a characteristic radius "rd" in enclosed in a radius r_ext = 11 x rd 
            rho_diff = M_unstr / (pi*(1.1d1*2.d0*h%galaxy%disc%rd)**2.) / h%galaxy%disc%h   ! [code unit: 10^11Msun/kpc^3]  
            av_rho   = 5.d-1*(rho_hot+rho_diff)
            Drho     = max(0.d0, (rho_hot-rho_diff)/av_rho)
            !
            L = h%galaxy%disc%L
            L = L / sqrt(sum(L**2.))
            !
            ! The characteristical time is setled to the min between 
            ! the dm%dynamical time and the dm dynamical friction_time 
            t_fric = dm_dynamical_friction_time(host%dm,h%dm)
            if ((host%dm%t_dyn .gt. 0.d0) .and. (t_fric .gt. 0.d0)) then
                !
                t_fric = min(host%dm%t_dyn,t_fric)
                inst_stripping_rate = disc_stripping_efficiency*abs(sum(L*dV))*Drho*M_unstr/t_fric
            else
                !
                call IO_print_error_message('t_fric <= 0.',only_rank=rank,called_by='halo_compute_stripping_rate')
                stop
            end if
        else
            ! no stripping, no hot gas
            inst_stripping_rate = 0.d0
        end if
        !
        if (is_NaN(inst_stripping_rate)) then
            call IO_print_warning_message('stripping_rate is NaN ',only_rank=rank,called_by='halo_compute_stripping_rate')
            inst_stripping_rate = 0.d0
        end if
        !
    end if

    return
  end subroutine halo_compute_stripping_rate
! -------------------------------------------------
#endif  
! SUB_STRIPPING

  !*****************************************************************************************************************
  
  function halo_baryonic_mass(h)

    ! RETURN THE TOTAL BARYONIC MASS ENCLOSE INTO THE HALO (-> doesn't take into account igm mass)

    implicit none

    real(kind=8)                    :: halo_baryonic_mass

    type(halo_type),intent(in)      :: h                    ! a given halo
    
    halo_baryonic_mass = baryon_halo_bh_mass(h%baryon_halo) + galaxy_mass(h%galaxy)

    return
  end function halo_baryonic_mass

  !*****************************************************************************************************************

  function halo_escape_velocity(h)

    ! RETURN THE ESCAPE VELOCITY OF THE HALO BY TAKING INTO ACCOUNT DARK-MATTER, BARYON HALO AND GALAXY MASS ENCLOSE INTO THE VIRIAL RADIUS

    implicit none

    real(kind=8)                :: halo_escape_velocity
    real(kind=8)                :: mass
 
    type(halo_type),intent(in)  :: h                      ! a given halo
    
    mass = dm_mass(h%dm) + baryon_halo_bh_mass(h%baryon_halo,component='cold+hot') + galaxy_mass(h%galaxy) 

    halo_escape_velocity = sqrt(gravconst_code_unit*(mass)/h%dm%R_vir)

    return
  end function halo_escape_velocity
  
  !*****************************************************************************************************************
  !
  ! PRINTING PROCEDURES
  !
  !*****************************************************************************************************************

  subroutine halo_print(unit,form,data_type,h,go_down)

    ! CALL PRINT SUBROUTINES TO SAVE DARK-MATTER (AND BARYONIC) PROPERTIES OF A GIVEN HALO h

    implicit none

    integer(kind=4),intent(in)    :: unit      ! file unit reference

    logical,optional              :: go_down   ! = .true. if higher level printing procedures have to be called
                                               ! only in merger tree printing procedure
    character(*)                  :: form      ! fits or tmp_bin 
    character(*),intent(in)       :: data_type ! phy or lum data

    type(halo_type),intent(in)    :: h         ! the halo component
    
    select case (trim(data_type))
    case('phy','physical','phy_data')
        !
        ! PHYSICAL PROPERTIES
        !
        if (present(go_down)) then
            if (go_down) then 
                !
                ! dm_print(unit,form,dm)
                call dm_print(unit,form,h%dm) ! call routine to print dm-properties   
                !
#ifdef TREE_EVOLVE
! ------------------------------------------------
                ! If the tree_evolve keyword is turn ON we must also print baryonic properties of galaxies contained in the dark-matter structure 
                ! call routine to print baryon_halo properties (cold and hot gas phase surrounding the galaxy)
                ! baryon_halo_print(unit,form,bh)
                call baryon_halo_print(unit,form,h%baryon_halo)  
                ! call routine to print galaxy properties  
                ! galaxy_print(unit,form,data_type,gal)
                call galaxy_print(unit,form,data_type,h%galaxy,go_down=go_down)  
! ------------------------------------------------
#endif
! TREE_EVOLVE
                !
            end if
        end if
        !
    case('lum','luminous','lum_data')
        !
        ! LUMINOUS PROPERTIES
        !
        ! galaxy_print(unit,form,data_type,gal)
        call galaxy_print(unit,form,data_type,h%galaxy) ! call routine to print galaxy properties
    case default
        call IO_print_error_message('Unknwon data type', &
                only_rank = rank, called_by = 'halo_print')
        stop ! stop the program
    end select  

    return
  end subroutine halo_print

  !*****************************************************************************************************************

end module halo
