module bulge

  use dm       ! Contains dm structure and dm exploitation functions
  use stars    ! Contains stars structure and stars exploitation functions
  use dust     ! Contains dust structure definition and dust exploitation functions
    
  public

  !*****************************************************************************************************************
  ! 
  !  OVERVIEW
  !
  !  bulge_module defines the bulge data structure univ(ts)%halo(ih)%galaxy%bulge
  !
  !  In the header of the module are defined output properties of the bulge component (labels, units and formats).
  !  The bulge module contains the defintion of this structure and all functions and subroutines associated to the bulge structure.
  !  In the current model, the bulge component is considered as a passive component, gas is only considered as a diffuse gas component
  !  The gas contained in the bulge comes from the stellar evolution cycle. 
  !  During a merger event, all the gas produced by stars in the bulge is transfered to the diffuse disc gas component
  !  The bulge evolution procedure is divided in two processes, bulge_evolve_I and bulge_evolve_II 
  !
  !  MPI PROCESSES IN THIS MODULE
  !
  !   bulge_send_data                                 : send specific informations about a bulge in a given galaxy
  !
  !   bulge_receive_data                              : receive specific informations about a bulge in a given galaxy
  !
  !  SUBROUTINES IN THIS MODULE
  !
  !   bulge_void                                      : init all properties of the bulge structure
  !
  !   bulge_copy                                      : copy a bulge component to an other one
  !
  !   bulge_evolve_I                                  : first part (PREDICTOR) of the bulge evolution scheme
  !      called by : galaxy_evolve                         compute input, output and evolution rates associated to the gas structuration history
  !
  !   bulge_evolve_II                                 : second part (CORRECTOR) of the bulge evolution scheme                                                  
  !       called by : galaxy_evolve
  !
  !   bulge_deallocate                                : deallocate all allocated array used in a bulge structure
  !
  !  FUNCTIONS IN THIS MODULE
  !
  !   bulge_mass                                      : returns the mass of the bulge component (or the mass of a sub-bulge component gas or stars)
  !
  !   bulge_gas_signature                             : return the signature of a given gas component
  !
  !   bulge_frac_mass_radius                          : returns the radius which enclose a given fraction of the bulge mass
  !  
  !   isolated_bulge_velocity (INTERFACE VERSION)     : returns the rotational velocity of an isolated bulge (Hernquist profile)
  !
  !   isolated_bulge_velocity_ (CORE VERSION)
  !
  !   stellar_velocity_dispersion                     : returns the stellar velovity dispersion (Herquist 1990 Eq. (10))
  !
  !  PRINTING PROCEDURES
  !
  !   bulge_load_bulge_data                            
  !      called by  : bulge_print                     :  create the output property list of the bulge component
  !
  !   disc_print                                      : print bulge-properties in output files
  !       called by : bulge_evolve_II
  !
  !*****************************************************************************************************************
  
  type bulge_type  
    ! Evolution
    real(kind=8)      :: life_time ! the time evolution of this bulbe component 
    real(kind=8)      :: age_form  ! age of the universe when the first bulge has been formed
    !
    ! main components
    type(gas_type)    :: gas       ! a diffuse gascomponent  (gas coming from stellar wind produced by stars evolving in the bulge)                     
    type(stars_type)  :: stars     ! the stellar component
    type(dust_type)   :: dust      ! the dust component
    ! bulge properties
    real(kind=8)      :: rb        ! the hernquist radius [kpc]
    real(kind=8)      :: t_dyn     ! the dynamical time of the bulge (V(size)/size) in [Gyr]
    real(kind=8)      :: dV        ! velocity dispersion of the stellar component
  end type bulge_type

  ! hdu reference for bulge structure
  integer(kind=4)           :: hdu_bulge
  ! printable properties for bulge structure
  integer(kind=4),parameter :: nb_bulge_field = 12 ! Number of bulge properties saved
  ! Name of each output colomn data
  character(len=ttype_len),dimension(nb_bulge_field) :: ttype_bulge = (/'bulge_age_form        ','bulge_unstr_gas       ','bulge_mZ              ', &
                                                                        'bulge_mH              ','bulge_mC              ','bulge_mN              ', &
                                                                        'bulge_mO              ','bulge_mFe             ','bulge_OH_index        ', &
                                                                        'bulge_rb              ','bulge_t_dyn           ','bulge_dV              '/) 
  ! Physical unit of each column data
  character(len=tunit_len),dimension(nb_bulge_field) :: tunit_bulge = (/'Gyr         ','M_sun       ','M_sun       ','M_sun       ','M_sun       ', &
                                                                        'M_sun       ','M_sun       ','M_sun       ','wounit      ','kpc         ', &
                                                                        'Gyr         ','km/s        '/)   
  ! Data type of each column data
  character(len=tform_len),dimension(nb_bulge_field) :: tform_bulge = (/'1E','1E','1E','1E','1E','1E','1E','1E','1E','1E','1E','1E'/)

  ! *** STARS ***
  ! Name of each output colomn data
  character(len=ttype_len),dimension(nb_stars_field) :: ttype_starsb = (/'bulge_stars_mass      ','bulge_stars_f_young   ', &
                                                                         'bulge_stars_age_MW    ','bulge_stars_dage_MW   ', &
                                                                         'bulge_stars_age_LW    ','bulge_stars_dage_LW   ', &
                                                                         'bulge_stars_Z_MW      ','bulge_stars_Z_LW      '/)
  ! Physical unit of each column data 
  character(len=tunit_len),dimension(nb_stars_field) :: tunit_starsb = (/'M_sun       ','wounit      ','Gyr         ','Gyr         ', &
                                                                         'Gyr         ','Gyr         ','Z_sun       ','Z_sun       '/)
  ! Data type of each column data
  character(len=tform_len),dimension(nb_stars_field) :: tform_starsb = (/'1E','1E','1E','1E','1E','1E','1E','1E'/)
  
  ! *** DUSTS ***
  ! Name of each output colomn data
  character(len=ttype_len),dimension(nb_dust_field) :: ttype_dustb = (/'bulge_dust_f_pah      ','bulge_dust_f_bg       ', &
                                                                       'bulge_dust_tau        ','bulge_dust_mZmH       '/)  
  ! Physical unit of each column data
  character(len=tunit_len),dimension(nb_dust_field) :: tunit_dustb = (/'w_o_unit    ','w_o_unit    ','w_o_unit    ','w_o_unit    '/)
 
  ! Data type of each column data
  character(len=tform_len),dimension(nb_dust_field) :: tform_dustb = (/'1E','1E','1E','1E'/)
  
contains

  !*****************************************************************************************************************
  ! 
  ! MPI PROCESSES
  !
  !*****************************************************************************************************************

#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------  
  subroutine bulge_send_data(bulge)                                     
  
    ! SEND SPECIFIC INFORMATIONS ABOUT A bulge
  
    implicit none

    logical                       :: go_down

    type(bulge_type),intent(in)   :: bulge
    
    ! data are sent by physical process and receive by luminous process
    
    go_down = .false.
    
    if (bulge_mass(bulge,component='stars') .gt. 0.d0) then
        !
        ! the bulge have a stellar population
        go_down = .true. 
        call MPI_SEND(go_down,1,MPI_LOGICAL,rank+1,bulge_tag+1,MPI_COMM_WORLD,ierror)
        ! send bulge typical radius
        call MPI_SEND(bulge%rb,1,MPI_REAL8,rank+1,bulge_tag+2,MPI_COMM_WORLD,ierror)
        ! send stellar population of the bulge
        call stars_send_data(bulge%stars)
        ! send dust component
        call dust_send_data(bulge%dust) 
    else
        !
        ! send exit loop order
        call MPI_SEND(go_down,1,MPI_LOGICAL,rank+1,bulge_tag+1,MPI_COMM_WORLD,ierror)
    end if
    
    return
  end subroutine bulge_send_data 
  
  !*****************************************************************************************************************
  
  subroutine bulge_receive_data(bulge)                                     
  
    ! RECEIVE SPECIFIC INFORMATIONS ABOUT A bulge
  
    implicit none

    logical                       :: go_down
    
    type(bulge_type),intent(out)  :: bulge
    
    ! data are sent by physical process and receive by luminous process
    
    call bulge_void(bulge)  ! init the bulge
    
    ! receive exit loop order
    call MPI_RECV(go_down,1,MPI_LOGICAL,rank-1,bulge_tag+1,MPI_COMM_WORLD,statut,ierror)
    
    if (go_down) then
        !
        ! receive bulge typical radius
        call MPI_RECV(bulge%rb,1,MPI_REAL8,rank-1,bulge_tag+2,MPI_COMM_WORLD,statut,ierror)
        ! receive stellar population of the bulge
        call stars_receive_data(bulge%stars)
        ! receive data from dust component
        call dust_receive_data(bulge%dust)
    end if
    
    return
  end subroutine bulge_receive_data 
! ------------------------------------------------- 
#endif
! LUMINOUS_PROCESSES

  !*****************************************************************************************************************

  subroutine bulge_void(bulge)

    ! INIT OR VOID A BULGE COMPONENT

    implicit none 
    type(bulge_type),intent(inout) :: bulge

    ! Evolution
    bulge%life_time         = 0.d0                ! will be incremented in bulge_evolve_II
    bulge%age_form          = 0.d0                ! will be set in galaxy_merge
    !
    ! main components
    call gas_void(bulge%gas)                      ! void the diffuse gas phase
    call stars_void(bulge%stars)                  ! void the stellar component
    call dust_void(bulge%dust,init_geom='dwek  ') ! void and init the dust component
    ! bulge properties
    bulge%rb                = 0.d0                ! will be set in galaxy_merge during the bulge formation process                      
    bulge%t_dyn             = 0.d0                ! init to null value          
    bulge%dV                = 0.d0                ! init to null value      
    
    return
  end subroutine bulge_void

  !*****************************************************************************************************************

  subroutine bulge_copy(bulge1,bulge2)
    
    ! COPY bulge2 IN bulge1 COMPONENT

    implicit none
    
    type(bulge_type),intent(inout) :: bulge1
    type(bulge_type),intent(in)    :: bulge2

    ! Evolution
    bulge1%life_time         = bulge2%life_time  
    bulge1%age_form          = bulge2%age_form   
    !
    ! main components 
    call gas_copy(bulge1%gas,bulge2%gas) 
    call stars_copy(bulge1%stars,bulge2%stars) 
    call dust_copy(bulge1%dust,bulge2%dust)
    !
    ! bulge properties
    bulge1%rb                = bulge2%rb          
    bulge1%t_dyn             = bulge2%t_dyn                 
    bulge1%dV                = bulge2%dV    
    
    return
  end subroutine bulge_copy
  
  !*****************************************************************************************************************
  
  subroutine bulge_evolve_I(bulge,dt_optim)

    ! FIRST STEP OF THE BULGE EVOLUTION SCHEME
    ! This subroutine compute evolution rate associated to gas and stellar population components.
    
    ! In the same scheme than for a disc structure, the bulge evolution scheme is based on two successive subroutines:
    ! bulge evolve_I and bulge_evolve_II
    ! the first routine computes evolution rates (mass transfer rate) the second one, applies these evolution rates onto gas and/or stellar components
    ! to compute new masses.
    ! The bulge structure is very simple in this model, it is composed of a stellar population and a diffuse gas component.
    ! No stars are formed in the bulge structure. 
    ! The gas accumulated during bulge evolution only comes from the stellar evolution.
    ! During a merger event, the gas is transferred in the disc structure. 
    
    implicit none
    
    real(kind=8),intent(out)       :: dt_optim             ! the global optimal evolution time-step  
    real(kind=8)                   :: dt_gas               ! the optimal time-step for gas component
    real(kind=8)                   :: dt_stars             ! the optimal time-step for the stellar population

    type(gas_type)                 :: void                 ! a null gas component
    type(bulge_type),intent(inout) :: bulge                ! the bulge component

#ifdef PRINTALL 
! -------------------------------------------------
    call IO_print_message('bulge_evolve_I',only_rank=rank,component='bulge')
! -------------------------------------------------
#endif
    
    ! *****************************
    ! GAS
    ! *****************************
    ! 
    ! In the bulge evolution scheme, the gas evolving into the bulge is only produced by stellar evolution, 
    !
    call gas_void(void)
    dt_gas = gas_dt_optim(bulge%gas,bulge%stars%loss_rate,void,called_by='bulge_evolve_I')
    !
    ! *****************************
    ! STARS
    ! *****************************
    !
    ! in the bulge evolution sheme no stars are formed into the bulge. 
    ! In this context the evolution of the the stellar popualtion is only linked to the loss_rate
    dt_stars = -1.d0 ! init
    dt_stars = stars_dt_optim(bulge%stars,void)
    !
    ! *****************************
    ! COMPUTE THE OPTIMAL TIME-STEP
    ! *****************************
    !
    dt_optim = -1.d0 ! init
    if (dt_gas .gt. 0.d0) dt_optim = dt_gas
    !
    if (dt_optim .gt. 0.d0) then
      !
      if (dt_stars .gt. 0.d0) then
        !
        dt_optim = min(dt_stars,dt_optim)
      end if  
    else
      !
      if (dt_stars .gt. 0.d0) then
        !
        dt_optim = dt_stars
      end if 
    end if
    !
    if (dt_optim .eq. 0.d0) then
        !
        call IO_print_warning_message('dt = 0.',only_rank=rank,called_by='bulge_evolve_I') 
#ifdef PRINT_WARNING
! -------------------------------------------------    
        call IO_print_message('used',only_rank=rank,component='bulge', &
                param_name = (/'dt_optim                 ','dt_gas                   ','dt_stars                 '/), &
                real_param_val  = (/dt_optim,dt_gas,dt_stars/)) 
! -------------------------------------------------
#endif
! PRINT_WARNING                
    end if    
    !
    if ((dt_optim .gt. 0.d0) .and. (dt_optim .lt. num_precision)) then
        !
        call IO_print_warning_message('Very small timestep',only_rank=rank,called_by='bulge_evolve_I') 
#ifdef PRINT_WARNING
! -------------------------------------------------    
        call IO_print_message('used',only_rank=rank,component='bulge', &
                param_name = (/'dt_optim                 ','dt_gas                   ','dt_stars                 '/), &
                real_param_val  = (/dt_optim,dt_gas,dt_stars/)) 
! -------------------------------------------------
#endif
! PRINT_WARNING                
    end if
    !
#ifdef PRINTALL 
! -------------------------------------------------
    call IO_print_message('bulge_evolve_I ... done ',only_rank=rank,component='bulge')
! -------------------------------------------------
#endif    

    return
  end subroutine bulge_evolve_I

  !*****************************************************************************************************************

  subroutine bulge_evolve_II(bulge,dt)

    ! SECOND STEP OF THE BULGE EVOLUTION SCHEME (EVOLVE MASS COMPONENT)
    ! take into account evolution rates computed in bulge_evolve_I

    implicit none
        
    real(kind=8),intent(in)           :: dt                    ! the evolution time-step
    real(kind=8)                      :: Ugas, dU            ! mass variation test

    type(gas_type)                    :: loss_rate_prior       ! mass loss rate used as prior value in bulge_evolve_I
    type(gas_type)                    :: loss_rate             ! effective stellar mass loss rate

    type(bulge_type),intent(inout)    :: bulge                 ! the bulge component

#ifdef PRINTALL 
! -------------------------------------------------
    call IO_print_message('bulge_evolve_II',only_rank=rank,component='bulge')
! -------------------------------------------------
#endif
    
    if (dt .le. 0.d0) then
      !
      call IO_print_error_message('dt <= 0.',called_by = 'bulge_evolve_II')
      stop ! stop the program
    end if
     
    if (bulge_mass(bulge) .gt. 0.d0) then
        !
        ! Evolution
        if ((bulge%life_time .eq. 0.d0) .and. (bulge%age_form .gt. 0.d0)) then
            !
            ! print the starting state of the bulge component
            if (FOLLOW_UP .and. PR_FOLLOW_UP) call bulge_print(follow_up_unit(current_index),'fits',bulge)
        end if
        !
        Ugas = bulge_mass(bulge,component='all')
        !
        ! *****************************
        ! STARS
        ! *****************************
        !
        ! Save the previous mass loss rate, used as prior
        loss_rate_prior = bulge%stars%loss_rate
        call stars_evolve(bulge%stars,dt,loss_rate)
        !
        ! *****************************
        ! GAS
        ! *****************************
        ! 
        ! Add gas ejected by stars to the diffuse gas reservoir (gas)
        call gas_add(bulge%gas,loss_rate*dt)
        !
        ! Test nosfg mass variation
        if (Ugas .gt. M_gas_min) then
            !
            dU = abs(bulge_mass(bulge,component='all')-Ugas)
            if ((abs(dU)/Ugas) .gt. physical_precision) then
              !
              call IO_print_error_message('No quasi-static state of the nosfg component',only_rank=rank,called_by='bulge_evolve_II')
              call IO_print_message('with',only_rank=rank,component='bulge', &
                          param_name = (/'dU/U (%)                 ','U                        ','dU                       ',&
                                         'loss_rate_prior          ','loss_rate                '/), &
                          real_param_val = (/1.d2*dU/Ugas,Ugas,dU,gas_mass(loss_rate_prior),gas_mass(loss_rate)/))
              stop  ! stop the program
            end if
        end if
        !
        ! *****************************
        ! DUST
        ! *****************************
        ! evolve dust properties
        ! dust_evolve(d,gas,rc,incl), rc is the half mass radius
        call dust_evolve(bulge%dust,bulge%gas,(1.d0+sqrt(2.d0))*bulge%rb)
        !
        ! *****************************
        ! SET bulge PROPERTIES
        ! *****************************
        !
        ! Add time to the bulge evolution time
        bulge%life_time = bulge%life_time + dt 
        ! Update stellar velocity dispersion
        bulge%dV = stellar_velocity_dispersion(bulge)
        !
        ! Print instantaneous bulge properties of the followed bulge
        if (FOLLOW_UP .and. PR_FOLLOW_UP) call bulge_print(follow_up_unit(current_index),'fits',bulge)
     end if
     
#ifdef PRINTALL 
! -------------------------------------------------
     call IO_print_message('bulge_evolve_II ... done ',only_rank=rank,component='bulge')
! -------------------------------------------------
#endif

     return
  end subroutine bulge_evolve_II
  
  !*****************************************************************************************************************
  
  subroutine bulge_deallocate(bulge)
  
    implicit none
    
    type(bulge_type),intent(inout) :: bulge  ! a bulge component

    call stars_deallocate(bulge%stars)

    return
  end subroutine bulge_deallocate

  !*****************************************************************************************************************
  !
  ! FUNCTIONS
  !
  !*****************************************************************************************************************

  function bulge_mass(bulge,r,component)
    
    ! RETURN THE TOTAL BULGE MASS
    ! if r is present, return mass enclosed in the radius r
    
    implicit none
    
    character(*),intent(in),optional  :: component    ! allow to select the bulge component (sfg, non-sfg, stars)
    character(MAXPATHSIZE)            :: message      ! a message to display

    real(kind=8),intent(in),optional  :: r            ! radius
    real(kind=8)                      :: bulge_mass   ! the bulge mass

    type(bulge_type),intent(in)       :: bulge        ! a bulge component 
    
    bulge_mass = gas_mass(bulge%gas) + stars_mass(bulge%stars)
    
    if (bulge_mass .le. 0.d0) return ! no bulge
    
    if (present(component)) then
       !
       select case (trim(component)) 
       case ('gas','all','diffuse')
          bulge_mass = gas_mass(bulge%gas)
       case ('stars')
          bulge_mass = stars_mass(bulge%stars)
       case default
          write(message,'(a,a,a)') 'Keyword ', trim(component), ' not defined'
          call IO_print_error_message(trim(message),only_rank=rank,called_by='bulge_mass')    
          stop  ! stop the program
       end select
    else
       !
       bulge_mass = gas_mass(bulge%gas) + stars_mass(bulge%stars)
    end if
    !
    ! ADAPT TO THE RADIUS R 
    if (present(r)) then
      !
      if (r .le. 0.d0) then
        !
        call IO_print_error_message('r <= 0.',only_rank=rank,called_by='bulge_mass')
        stop  ! stop the program
      end if
      if (bulge%rb .le. 0.d0) then
        !
        call IO_print_error_message('rb <= 0.',only_rank=rank,called_by='bulge_mass')
        stop  ! stop the program
      endif
      bulge_mass = bulge_mass*r**2./(r+bulge%rb)**2.   ! apply geometrical function
    end if
    !
    ! CHECK the results and crash the code if bulge_mass is NAN
    if (is_NaN(bulge_mass)) then
      !
      call IO_print_error_message('bulge_mass is NaN',only_rank=rank,called_by='bulge_mass')  
      call IO_print_message('with',only_rank=rank,component='bulge', &
               param_name = (/'non-sfg                  ','stars                    '/), &
               real_param_val  = (/gas_mass(bulge%gas),stars_mass(bulge%stars)/))  
      stop  ! stop the program
    end if
      
    return
  end function bulge_mass
  
  !*****************************************************************************************************************
  
  function bulge_gas_signature(bulge,component)
  
    ! RETURN A GAS SIGNATURE (A NORMALIZED GAS OBJECT)
    
    implicit none
    
    character(*),intent(in),optional  :: component           ! allow to select the bulge component (sfg, non-sfg)
    character(MAXPATHSIZE)            :: message             ! a message to display
    
    type(gas_type)                    :: bulge_gas_signature ! the gas signature
    
    type(bulge_type),intent(in)       :: bulge               ! a bulge component
    
    if (present(component)) then
        !
        select case (trim(component)) 
        case ('gas','all','diffuse')
            bulge_gas_signature = gas_signature(bulge%gas)
        case default
            write(message,'(a,a,a)') 'Keyword ', trim(component), ' not defined'
            call IO_print_error_message(message,only_rank=rank,called_by='bulge_gas_signature')    
            stop  ! stop the program
        end select
    else
        !
        ! use all gas as defaut value
        bulge_gas_signature = gas_signature(bulge%gas)
    end if
        
    return
  end function bulge_gas_signature
  
  !*****************************************************************************************************************
  
  function bulge_frac_mass_radius(bulge,frac,component) 
    
    ! RETURN THE RADIUS THAT CONTAINS 100*frac % OF THE bulge MASS
    
    implicit none 
    
    character(*),intent(in),optional  :: component               ! allow to select the bulge component (sfg, non-sfg, stars, agn)

    real(kind=8),intent(in)           :: frac                    ! fraction of the bulge mass
    real(kind=8)                      :: bulge_frac_mass_radius  ! the radius which enclose the given mass fraction 
    real(kind=8)                      :: mass, r
    real(kind=8)                      :: r_min, r_max            ! [r_min : r_max] interval in which a solution for galaxy_frac_mass_radius is searched

    type(bulge_type),intent(in)        :: bulge                  ! a bulge component
    
    bulge_frac_mass_radius = 0.d0
    if (bulge_mass(bulge,component=component) .le. 0.d0) return ! no bulge 
    
    r_min = 0.d0              ! init
    r_max = 1.d2*bulge%rb      ! init
    if (r_max .le. 0.d0) then
      !
      call IO_print_error_message('r_max <= 0', &
            only_rank = rank, called_by = 'bulge_frac_mass_radius')
      call IO_print_message('used',only_rank=rank,component='bulge',& 
        param_name=(/'bulge%rb                 '/), &
        real_param_val =(/bulge%rb/))
      stop ! stop the program
    end if
    !
    ! Test with r_max and crash the code if the optimal radius is larger than this barrier.
    mass = bulge_mass(bulge,r=r_max,component=component)
    if (mass .lt. frac*bulge_mass(bulge,component=component)) then
      !
      call IO_print_error_message('Optimal r > r_max', &
            only_rank = rank, called_by = 'bulge_frac_mass_radius')
      call IO_print_message('used',only_rank=rank,component='bulge',& 
        param_name=(/'bulge%rb                 '/), &
        real_param_val =(/bulge%rb/))
      stop ! stop the program   
    end if
    
    r = r_max  ! init
    do while (((abs(mass-frac*bulge_mass(bulge,component=component))/(frac*bulge_mass(bulge,component=component))) .gt. num_precision) &
      .and. (r_min .ne. r_max))
       !
       r = 5.d-1*(r_min+r_max)
       mass = bulge_mass(bulge,r=r,component=component) 
       ! The bulge mass is an incresing function of the radius
       if (mass .le. frac*bulge_mass(bulge,component=component)) then
          !
          r_min = r
       else
          !
          r_max = r
       end if
    end do
    !
    ! Test the result and crash the code if the optimal radius is smaller than 0.
    if (r .gt. 0.d0) then
      !
      bulge_frac_mass_radius = r
    else
      !
      call IO_print_error_message('Optimal r < 0', &
            only_rank = rank, called_by = 'bulge_frac_mass_radius')
      call IO_print_message('used',only_rank=rank,component='bulge',& 
        param_name=(/'bulge%rb                 '/), &
        real_param_val =(/bulge%rb/))
      stop ! stop the program  
    end if
    
    return
  end function bulge_frac_mass_radius

  !*****************************************************************************************************************

  function isolated_bulge_velocity(r,bulge) ! INTERFACE VERSION

    ! COMPUTE THE ROTATIONNAL VELOCITY OF A ISOLATED BULGE

    implicit none

    real(kind=8),intent(in)       :: r                       ! the radius
    real(kind=8)                  :: isolated_bulge_velocity ! in code unit

    type(bulge_type),intent(in)   :: bulge                   ! the bulge component

    isolated_bulge_velocity = 0.d0
    
    if (bulge_mass(bulge) .le. 0.d0) return
    
    isolated_bulge_velocity = isolated_bulge_velocity_(r,(/bulge%rb,bulge_mass(bulge,r=r)/))
    
    return
  end function isolated_bulge_velocity

  !*****************************************************************************************************************

  function isolated_bulge_velocity_(r,param)

    ! COMPUTE THE ROTATIONNAL VELOCITY OF A ISOLATED BULGE

    implicit none

    real(kind=8),intent(in)   :: r          ! the radius
    real(kind=8),intent(in)   :: param(2)   ! parameter array : param(1) = bulge%rb
                                            !                   param(2) = bulge_mass(bulge,r=r)
    real(kind=8)              :: isolated_bulge_velocity_ ! in code unit
    
    isolated_bulge_velocity_ = 0.d0  ! init
    
    if (param(2) .le. 0.d0) return
    
    isolated_bulge_velocity_ = sqrt(gravconst_code_unit*param(2)/r)
    
    return
  end function isolated_bulge_velocity_
  
  !*****************************************************************************************************************
  
  function stellar_velocity_dispersion(bulge)
  
    ! RETURN THE STELLAR VELOCITY DISPERSION IN A GIVEN BULGE
    ! Based on Herquist+1990 Eq. (10)
    
    implicit none
    
    real(kind=8)                :: r50,r0,rb      ! radius
    real(kind=8)                :: Ms             ! stellar mass
    real(kind=8)                :: dV             ! stellar velocity dispersion
    real(kind=8)                :: stellar_velocity_dispersion
    
    type(bulge_type),intent(in) :: bulge          ! the bulge component
  
    stellar_velocity_dispersion = 0.d0 ! init
    
    ! compute stellar mass
    Ms = bulge_mass(bulge,component='stars')
    
    if (Ms .le. 0.d0) return
    ! @ this point the bulge contains stars
    
    ! compute half-mass radius
    r50 = bulge_frac_mass_radius(bulge,5.d-1)  ! the half_mass radius [~ (1+sqrt(2))*bulge%rb]
    ! save rb
    rb = bulge%rb
    ! compute r0
    r0  = R50/rb
    !                                             
    dV = 1.2d1*r50*(r50+rb)**3./rb**4.*log((r50+rb)/r50)
    dV = dV - (r50/(r50+rb))*(2.5d1+5.2d1*r0+4.2d1*r0**2.+1.2d1*r0**3.)
    dV = dV*gravconst_code_unit*Ms/(1.2d1*rb)
    dV = sqrt(dV)                              ! Code unit kpc/Gyr
    
    stellar_velocity_dispersion = dV
    
    return
  end function stellar_velocity_dispersion

  !*****************************************************************************************************************
  !
  ! PRINTING PROCEDURES
  !
  !*****************************************************************************************************************

  subroutine bulge_load_bulge_data(bulge,bulge_data)

    ! CREATE THE BULGE OUTPUT PROPERTIES LIST

    implicit none

    real(kind=8), intent(inout) :: bulge_data(nb_bulge_field+nb_stars_field+nb_dust_field)
    real(kind=8)                :: starsb_data(nb_stars_field)
    real(kind=8)                :: dust_data(nb_dust_field)

    type(bulge_type),intent(in) :: bulge   ! bulge component

    call stars_load_stars_data(bulge%stars,starsb_data)
       
    call dust_load_dust_data(bulge%dust,dust_data)

    ! For information
    !'bulge_age_form        ','bulge_gas          ','bulge_mZ              ', &
    !'bulge_mH              ','bulge_mC              ','bulge_mN              ', &
    !'bulge_mO              ','bulge_mFe             ','bulge_m_stars         ', &
    !'bulge_rb              ','bulge_t_dyn           ','bulge_dV              '
    
    bulge_data = (/bulge%age_form, mass_code_unit_in_M_Sun*gas_mass(bulge%gas), &
                   mass_code_unit_in_M_Sun*gas_mass(bulge%gas,component='Metals'), &
                   mass_code_unit_in_M_Sun*gas_mass(bulge%gas,component='H1'), &
                   mass_code_unit_in_M_Sun*gas_mass(bulge%gas,component='C12'), &
                   mass_code_unit_in_M_Sun*gas_mass(bulge%gas,component='N14'), &
                   mass_code_unit_in_M_Sun*gas_mass(bulge%gas,component='O16'), &
                   mass_code_unit_in_M_Sun*gas_mass(bulge%gas,component='Fe56'), 12.d0 + log10(O_H(bulge%gas)), &
                   bulge%rb, bulge%t_dyn, bulge%dV*vel_code_unit_2_kmPerSec,starsb_data,dust_data/)

    return
  end subroutine bulge_load_bulge_data

 !*****************************************************************************************************************

  subroutine bulge_print(unit,form,bulge)

    ! PRINT BULGE-PROPERTIES IN eGALICS OUTPUT FILE

    implicit none

    integer(kind=4),intent(in)  :: unit    ! file unit
    integer(kind=4)             :: status,i,hdutype
    
    character(*)                :: form  ! fits or tmp_bin

    real(kind=8)                :: bulge_data(nb_bulge_field+nb_stars_field+nb_dust_field)

    type(bulge_type),intent(in) :: bulge  ! bulge component
    
    call bulge_load_bulge_data(bulge,bulge_data)
    
    select case (trim(form))
      case ('tmp_bin')
        write(unit) bulge_data  ! directly write data in the tmp binary output file 
      case ('fits')
        ! move to dm extension
        call ftmahd(unit,hdu_bulge,hdutype,status) 
        !
        if (status .gt. 0) then
          !
          call IO_print_error_message('ftmahd status', &
                only_rank = rank, called_by = 'bulge_print')
          stop ! stop the program
        end if
        !
        ! init
        call ftirow(unit,0,1,status) 
        if (status .gt. 0) then
          !
          call IO_print_error_message('ftirow status', &
                only_rank = rank, called_by = 'bulge_print')
          stop ! stop the program
        end if
        ! write data in the dm entension
        call ftpcld(unit,1,1,1,1,bulge%age_form+bulge%life_time,status)
        do i=2, nb_bulge_field+nb_stars_field+nb_dust_field+1   
          call ftpcld(unit,i,1,1,1,bulge_data(i-1),status)  
          if (status .gt. 0) then
            !
            call IO_print_error_message('ftpcld status', &
                only_rank = rank, called_by = 'bulge_print')
            stop ! stop the program
          end if
        end do
      case default
        call IO_print_error_message('Unknwon output data format', &
                only_rank = rank, called_by = 'bulge_print')
        stop ! stop the program
    end select

    return
  end subroutine bulge_print

  !*****************************************************************************************************************

end module bulge
