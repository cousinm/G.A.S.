module stars

  use stellar_population_library   ! Contains stellar population library (e.g. acces to loss mass rate, spectra)

  public
  
  !*****************************************************************************************************************
  ! 
  !  OVERVIEW
  !
  !  The stars module define all properties and procedures assossiated to a stellar popualtion.
  !  In G.A.S. a stellar population is mainly defined with a star-formation-history table (sfh-tab). 
  !  This table, composed of a large number of age and metalicity bins, allows to follow the amount of stellar mass
  !  in each age and metallicity bin.
  !  This module starts with the definition of the stars type
  !  In the header of the module are defined output properties of a star component (labels, units and formats)
  !
  !  MPI PROCESSES IN THIS MODULE
  !
  !  stars_send_data                   : send a stars structure from a MPI process to an other
  ! 
  !  stars_receive_data                : receive a stars structure from a given MPI process
  ! 
  !  SUBROUTINES IN THIS MODULE
  !
  !  stars_set_reference_mass_and_age  : Initialize stars reference mass : M_stars_min and young_stars_iMaxAge
  !
  !  stars_void                        : init all properties of a stellar popualtion
  !     called by : disc_void
  !               : bulge_void
  !  
  !  stars_reset_properties            : reset all sfh_tab integrated properties (Age, Z, Lum, SN_rate, inst_loss_rate ...)
  !     called by : stars_add
  !               : stars_sub
  !               : stars_evolve
  !
  !  stars_copy                        : copy the properties of a stellar population in an other
  !     called by : disc_copy
  !               : bulge_copy
  !
  !  stars_compute_properties          : compute sfh_tab integrated parameters (Age, Z, SN_rate, inst_loss_rate ...)
  !
  !  stars_transfer_young_stars        : transfer a young stellar population from a stellar population to an other one
  !     called by: galaxy merge
  !
  !  stars_add                         : add a stellar population to an other one
  !
  !  stars_sub                         : substract a stellar population to an other one
  ! 
  !  stars_add_new_stars               : add new formed stars in a (new or pre-existing) stellar population 
  !
  !  stars_evolve                      : evolve a given stellar population during dt
  !
  !  stars_build_stellar_spectra       : build the mass-weighted stellar spectra of a given stellar population
  !
  !  stars_deallocate                  : deallocate the sfh_tab array of a given stellar population
  !
  !  FUNCTIONS IN THIS MODULE 
  !
  !  stars_mass                        : return the stellar mass of a given stellar population
  !
  !  stars_dt_optim                    : return the optimal evolution timestep associated to a given stars component
  !
  !  PRINTING PROCEDURES
  !
  !  stars_load_stars_data             : create the stars output properties list
  !
  !  INTERFACE OPERATORS IN THIS MODULE
  !
  !  stars_add_                        : return stars1 + stars2
  !
  !  stars_sub_                        : return if possible stars1 - stars2
  !
  !  stars_scalar_multiply_left        : return a * stars (all masses multiplied)
  !
  !  stars_scalar_multiply_right       : return stars * a (all masses multiplied)
  !
  !*****************************************************************************************************************

  ! DEFINITION OF GLOBAL VARIABLES LINKED TO THE STARS STRUCTURE *******************

  ! STARS_TYPE DEFINITION **************

  type stars_type
     ! main properties
     real(kind=8)             :: mass               ! stellar mass (1.d11 M_Sun)
     real(kind=8)             :: Z(2)               ! mass/lum-weighted metalicity of the stellar population Z(1): mass-weighted, Z(2): lum-weighted
     real(kind=8)             :: dZ(2)              ! error on the mass/lum-weighted metalicity of the stellar population
     real(kind=8)             :: Age(2)             ! mass/lum-weighted age of the stellar population (Gyr) age(1): mass-weighted, age(2): lum-weighted
     real(kind=8)             :: dAge(2)            ! error on the mass/lum-weighted age of the stellar population (Gyr)
     real(kind=8)             :: f_young            ! fraction of young stars in the population
     ! star formation history
     logical                  :: exist              ! true if sfh_tab is allocated, false otherwise. 
     integer(kind=4)          :: iMaxAge            ! index in AgeBins(:)
     real(kind=8),allocatable :: sfh_tab(:,:,:)     ! sfh_tab(:,:,1): stellar mass in age and metallicity bin (1.d11 M_Sun) 
                                                    ! sfh_tab(:,:,2) clock (Gyr)
     ! other properties
     type(gas_type)           :: loss_rate          ! instantaneous mass loss rate rate of the stellar population
     real(kind=8)             :: SN_rate            ! instantaneous SN event rate [nb/Gyr]
     real(kind=8)             :: Lum                ! bolometric stellar luminosity
     real(kind=8)             :: ISRF(2)            ! mass weighted Interstellar radiation field (computed under stellar population ages 1: young stars only / 2: old stars only)
  end type stars_type

  ! INTERFACE OPERATOR DECLARATIONS ****************

  interface assignment (=)  ! allows to assign stars1 = stars2 (based on stars copy)
     module procedure stars_copy
  end interface assignment (=)

  interface operator (+)    ! allows to sum two star components by using the symbol '+'
     module procedure stars_add_
  end interface

  interface operator (-)    ! allows to substract two star components by using the symbol '-'
     module procedure stars_sub_
  end interface
   
  interface operator (*) 
     module procedure stars_scalar_multiply_left,stars_scalar_multiply_right
  end interface operator (*)

  ! OTHER DEFINITIONS *******************

  integer(kind=4)           :: young_stars_iMaxAge         ! the index of the AgeBin corresponding to the young stars maximum age
  
  ! minimal mass in a sfh_tab bin to compute stellar evolution in this bin
  real(kind=8)              :: M_stars_min                 ! Mminimal stellar mass takes into account in adaptive time step evolution process
                                                           ! set during initialization (IO_read_parameter_file)
  ! printable properties for stars structure
  integer(kind=4),parameter :: nb_stars_field = 8          ! Number of stars properties saved

contains
  
  !*****************************************************************************************************************
  ! 
  ! MPI PROCESSES
  !
  !*****************************************************************************************************************
 
#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------    
  subroutine stars_send_data(s)
  
    implicit none
    
    integer(kind=4)             :: Ncells ! number of cells in the sfh_tab with mass
    
    real(kind=8)                :: Ms     ! The stellar mass
    
    type(stars_type),intent(in) :: s      ! the stellar population
    
    Ms = stars_mass(s)
    
    ! send the existing criteron s%exist
    call MPI_SEND(s%exist,1,MPI_LOGICAL,rank+1,stars_tag+1,MPI_COMM_WORLD,ierror)
    !
    if (s%exist) then
        !
        ! The stellar population exists --> send data
        !
        ! send mass
        call MPI_SEND(Ms,1,MPI_REAL8,rank+1,stars_tag+2,MPI_COMM_WORLD,ierror)
        !
        ! send iMaxAge
        call MPI_SEND(s%iMaxAge,1,MPI_INTEGER4,rank+1,stars_tag+3,MPI_COMM_WORLD,ierror)
        !   
        ! send sfh_tab
        ncells = s%iMaxAge*nMetBins*2
        call MPI_SEND(s%sfh_tab(1:s%iMaxAge,:,:),ncells,MPI_REAL8,rank+1,stars_tag+4,MPI_COMM_WORLD,ierror)
    else
        !
        if (Ms .gt. 0.d0) then
            !
            call IO_print_error_message('No existing stellar population with no null mass', &
                        only_rank = rank, called_by = 'stars_send_data')
            call IO_print_message('with',only_rank=rank,component='stars', &
                        param_name = (/'Ms                       '/), real_param_val  = (/Ms/))
            stop  ! stop the program
        end if
    end if
    
    return
  end subroutine stars_send_data
  
  !*****************************************************************************************************************
  
  subroutine stars_receive_data(s)
  
    implicit none
    
    integer(kind=4)              :: Ncells ! number of cells in the sfh_tab with mass
    
    type(stars_type),intent(out) :: s      ! the stellar population
    
    call stars_void(s)  ! init the stellar population
    
    ! receive the existing criterion
    call MPI_RECV(s%exist,1,MPI_LOGICAL,rank-1,stars_tag+1,MPI_COMM_WORLD,statut,ierror)
    !
    if (s%exist) then
        !
        ! receive the mass
        call MPI_RECV(s%mass,1,MPI_REAL8,rank-1,stars_tag+2,MPI_COMM_WORLD,statut,ierror)
        !
        ! receive iMaxAge
        call MPI_RECV(s%iMaxAge,1,MPI_INTEGER4,rank-1,stars_tag+3,MPI_COMM_WORLD,statut,ierror)
        !
        Ncells = s%iMaxAge*nMetBins*2
        ! create the sfh_tab
        allocate(s%sfh_tab(nAgeBins,nMetBins,2))
        ! init the sfh_tab
        s%sfh_tab(:,:,:) = 0.d0
        ! receive sfh_tab
        call MPI_RECV(s%sfh_tab(1:s%iMaxAge,:,:),Ncells,MPI_REAL8,rank-1,stars_tag+4,MPI_COMM_WORLD,statut,ierror)
    end if
    
    return
  end subroutine stars_receive_data
  
! ------------------------------------------------- 
#endif
! LUMINOUS_PROCESSES

  !*****************************************************************************************************************
  !
  ! SUBROUTINES
  !
  !*****************************************************************************************************************

  subroutine stars_set_reference_mass_and_age

    ! INITIALIZE - M_stars_min: A MASS PARAMETERS LINKED TO THE DARK-MATTER N-BODY SIMULATION
    !            - young_stars_iMaxAge
    ! WARNING: Must be call after stellar_population_library_read_stellar_spectra

    implicit none

    integer(kind=4)  :: iAge  ! index for age loop
    
    M_stars_min = physical_precision * baryon_fraction * dm_particle_mass
    
    do iAge = 2, nAgeBins
       !
       if (AgeBins(iAge) .gt. young_stars_MaxAge) then
          !
          young_stars_iMaxAge = iAge -1
          exit
       end if
    end do
    
    call IO_print_message('use',param_name=(/'M_stars_min  [Msun]      ','young_stars_MaxAge [Gyr] '/), &
            real_param_val=(/M_stars_min*mass_code_unit_in_M_Sun,young_stars_MaxAge/))
    call IO_print_message('use',param_name=(/'young_stars_iMaxAge      '/),int_param_val=(/young_stars_iMaxAge/))
                            
    return
  end subroutine stars_set_reference_mass_and_age

  !*****************************************************************************************************************

  subroutine stars_void(s)

    ! SET TO NULL VALUES ALL COMPONENTS OF A STARS OBJECT

    implicit none 

    type(stars_type),intent(inout) :: s
    
    s%mass    = 0.d0       ! initialy the stellar popualtin does not exsit and therfore does not contain mass
    s%iMaxAge = -1         ! index of the maximal occuped age bin
    if (allocated(s%sfh_tab)) deallocate(s%sfh_tab)
    if (s%exist) s%exist = .false.
    call stars_reset_properties(s)

    return
  end subroutine stars_void

  !*****************************************************************************************************************
  
  subroutine stars_reset_properties(stars)   
  
    ! RESET SFH_TAB INTEGRATED PROERTIES   
    
    implicit none
        
    type(stars_type),intent(inout) :: stars
  
    stars%f_young = 0.d0   ! erase
    stars%Age     = 0.d0   ! erase
    stars%dAge    = 0.d0   ! erase
    stars%Z       = 0.d0   ! erase
    stars%dZ      = 0.d0   ! erase
    stars%Lum     = 0.d0   ! erase
    stars%sn_rate = 0.d0   ! erase
    stars%ISRF    = 0.d0   ! erase
        
    call gas_void(stars%loss_rate)  ! erase
    
    return
  end subroutine stars_reset_properties
  
  !*****************************************************************************************************************

  subroutine stars_copy(stars1,stars2) 

    ! SET STARS1 = STARS2 (if possible : stars1 not already exist !!)
    
    implicit none

    type(stars_type),intent(inout) :: stars1  ! a stellar population
    type(stars_type),intent(in)    :: stars2  ! an other one
    
    if (.not. stars1%exist) then ! stars1 no exist
       !
       ! the stellar population stars1 doesn't exists yet
       if (stars_mass(stars1) .gt. 0.d0) then
          !
          ! but it has got a mass of stars
          call IO_print_error_message('No existing stellar population with no null mass (stars1)', &
               only_rank = rank, called_by = 'stars_copy')
          call IO_print_message('with',only_rank=rank,component='stars', &
               param_name = (/'stars_mass               '/), real_param_val  = (/stars_mass(stars1)/))
          stop  ! stop the program
       end if
    end if
    
    ! at this point, the stellar population stars1 doesn't exist

    if (stars2%exist) then
       !
       ! the stellar population stars2 already exists
       stars1%mass    = stars2%mass    ! copy the stellar mass
       stars1%Z       = stars2%Z       ! copy the stellar metalicity
       stars1%dZ      = stars2%dZ      ! copy the error on the stellar metalicity 
       stars1%Age     = stars2%Age     ! copy the stellar age
       stars1%dAge    = stars2%dAge    ! copy the error on the stellar age
       stars1%iMaxAge = stars2%iMaxAge ! copy max occuped age index in the sfh_tab 
       !
       if (.not. allocated(stars1%sfh_tab)) then
          !
          ! the sfh_tab of the stellar population stars1 doesn't exist yet (not allocated)
          allocate(stars1%sfh_tab(nAgeBins,nMetBins,2)) ! create the sfh_tab
          stars1%sfh_tab(:,:,:)     = 0.d0              ! set to null values         
       end if
       !
       stars1%sfh_tab            = stars2%sfh_tab    ! copy the sfh_tab
       stars1%exist              = .true.            ! Now the stellar population exists
       ! copy the stellar wind rate
       call gas_copy(stars1%loss_rate,stars2%loss_rate) 
       stars1%SN_rate = stars2%SN_rate  ! copy the instantaneous SN event rate
       stars1%Lum     = stars2%Lum      ! copy the bolometric stellar luminosity
       stars1%ISRF    = stars2%ISRF     ! Interstellar radiation field
       ! at this point stars1 = stars2
    end if

    return
  end subroutine stars_copy

  !*****************************************************************************************************************

  subroutine stars_compute_properties(stars,iAge,iMet)

    ! COMPUTE: - THE MASS/LUM-WEIGHTED AGE 
    !          - THE MASS/LUM-WEIGHTED STELLAR METALLICITY 
    !          - THE BOLOMETRIC STELLAR LUMINOSITY 
    !          - THE INSTANTANEOUS SN RATE
    !          - THE INSTANTANEOUS MASS LOSS RATE
    ! OF A GIVEN STELLAR POPUALTION

    implicit none

    integer(kind=4),optional       :: iAge,iMet  ! indexes
    
    real(kind=8)                   :: Ms,MsL     ! total stellar mass, Ms*luminosity
    real(kind=8)                   :: m_bin      ! mass in a bin
    real(kind=8)                   :: ml_bin     ! mass*luminosity in a bin
    real(kind=8)                   :: t_bin      ! clock linked to the bin

    type(stars_type),intent(inout) :: stars      ! a stellar population     
    
    if (stars%iMaxAge .le. 0) return 
    
    if (present(iAge) .and. present(iMet)) then
      ! the routine has been called in a integration loop
      !
      m_bin  = stars%sfh_tab(iAge,iMet,1)   
      ml_bin = stars%sfh_tab(iAge,iMet,1)*LumBins(iAge,iMet)  ! LumBins(iAge,iMet) is in [Lsun/(10^11 Msun)]
      t_bin  = stars%sfh_tab(iAge,iMet,2)
      !
      if (m_bin .gt. 0.d0) then  
        !
        ! STELLAR AGE
        !
        ! mass-weighted
        stars%Age(1)  = stars%Age(1)  + m_bin*t_bin
        ! luminosity-weighted
        stars%Age(2)  = stars%Age(2)  + ml_bin*t_bin
        ! errors
        ! we use **2. to compute the quadradic error
        !
        if (iAge .lt. nAgeBins) then
           !
           stars%dAge(1) = stars%dAge(1) + (m_bin*5.d-1*(AgeBins(iAge+1) - AgeBins(iAge)))**2. 
           stars%dage(2) = stars%dage(2) + (ml_bin*5.d-1*(AgeBins(iAge+1) - AgeBins(iAge)))**2. 
        else
           !
           stars%dAge(1) = stars%dAge(1) + (m_bin*AgeBins(iAge))**2. 
           stars%dAge(2) = stars%dAge(2) + (ml_bin*AgeBins(iAge))**2. 
        end if
        !
        ! STELLAR METALLICITY
        !
        ! mass-weighted
        stars%Z(1)  = stars%Z(1) + m_bin*MetBins(iMet) 
        ! luminosity-weighted
        stars%Z(2)  = stars%Z(2) + ml_bin*MetBins(iMet)
        ! errors
        ! we use **2. to compute the quadradic error
        !
        if (iMet .lt. nMetBins) then
           !
           stars%dZ(1) = stars%dZ(1) + (m_bin*5.d-1*(MetBins(iMet+1)-MetBins(iMet)))**2.
           stars%dZ(2) = stars%dZ(2) + (ml_bin*5.d-1*(MetBins(iMet+1)-MetBins(iMet)))**2.
        else
           !
           stars%dZ(1) = stars%dZ(1) + (m_bin*MetBins(iMet))**2.
           stars%dZ(2) = stars%dZ(2) + (ml_bin*MetBins(iMet))**2.
        end if
        !
        ! STELLAR LUMINOSITY
        !
        stars%Lum = stars%Lum + ml_bin
        !
        ! INSTANTANEOUS MASS LOSS RATE
        !
        stars%loss_rate = stars%loss_rate + m_bin*mass_loss_rates(iAge,iMet) 
        !
        ! INSTANTANEOUS SN EVENT RATE
        !
        stars%SN_rate = stars%SN_rate + m_bin*SN_rates(iAge,iMet) 
      end if
    else
        !
        ! the routine has been called at the end of a integration loop
        ! normalized intagration value to stellar mass
        !
        ! STELLAR MASS (fraction of young stars in the population)
        Ms  = stars_mass(stars)
        stars%f_young = 0.d0
        if ((Ms .gt. 0.d0) .and. (allocated(stars%sfh_tab))) then
            !
            stars%f_young = max(0.d0,min(1.d0,sum(stars%sfh_tab(1:young_stars_iMaxAge,:,1))/Ms))
            !
            if ((stars%f_young .gt. 1.d0) .or. (stars%f_young .lt. 0.d0)) then
                !
                call IO_print_error_message('f_young < 0. or f_young > 1.',only_rank=rank,called_by='stars_compute_properties')
                write(*,*) 'f_young: ', stars%f_young
                stop ! stop the program
            end if 
        end if
        !
        ! luminosity weighted factor
        MsL = sum(stars%sfh_tab(:,:,1)*LumBins(:,:))
        !
        ! STELLAR AGE
        !
        ! AgeBin(1) = 0, we set the minimal age of the stellar population to 0.5*AgeBins(2)
        ! For the new stellar population, the error onto the age (dage) is equal to its age
        ! mass-weighted
        stars%Age(1)  = max(5.d-1*AgeBins(2),stars%Age(1)/Ms)   
        ! we use here a mass weighted quadratic error
        stars%dAge(1) = max(5.d-1*AgeBins(2),sqrt(stars%dAge(1))/Ms) 
        ! luminosity-weighted
        stars%Age(2)  = max(5.d-1*AgeBins(2),stars%Age(2)/MsL)   
        ! we use here a luminosity weighted quadratic error
        stars%dAge(2) = max(5.d-1*AgeBins(2),sqrt(stars%dAge(2))/MsL) 
        !
        ! STELLAR METALLICITY
        !
        ! we set the minimal metallicity of the stellar population to MetBins(1)
        ! For a new stellar population, the error onto the metallicity is equal to difference between the two first bins
        ! mass-weighted
        stars%Z(1)  = max(MetBins(1),stars%Z(1)/Ms)
        ! we use here a mass weighted quadratic error
        stars%dZ(1) =  max(MetBins(2)-MetBins(1),sqrt(stars%dZ(1)/Ms)) 
        ! luminosity-weighted
        stars%Z(2)  = max(MetBins(1),stars%Z(2)/MsL)
        ! we use here a mass weighted quadratic error
        stars%dZ(2) =  max(MetBins(2)-MetBins(1),sqrt(stars%dZ(2)/MsL)) 
        !
        ! TESTS
        if (gas_metallicity(stars%loss_rate) .lt. 0.d0) then
            !
            call IO_print_error_message('Metallicity of the instantaneous stellar ejecta is negative (i.e Z < 0)', &
                only_rank = rank, called_by = 'stars_compute_properties')
            call IO_print_message('with',only_rank=rank,component='stars', &
                param_name = (/'loss_rate%mass           ','loss_rate%mZ             '/), real_param_val  = (/stars%loss_rate%mass,stars%loss_rate%mZ/))
            stop ! stop the program
        end if
        !
        if (stars%SN_rate .lt. 0.d0) then
            !
            call IO_print_error_message('SN event rate is negative)', &
                only_rank = rank, called_by = 'stars_compute_properties')
            call IO_print_message('with',only_rank=rank,component='stars', &
                param_name = (/'stars%mass               ','stars%SN_rate            '/), real_param_val  = (/stars%mass,stars%SN_rate/))
            stop ! stop the program
        end if
    end if
    
    return
  end subroutine stars_compute_properties
  
  !*****************************************************************************************************************
  
  subroutine stars_transfer_young_stars(s1,s2)
  
    ! TRANSFER YOUNG STARS FROM s1 TO s2
    
    implicit none
    
    type(stars_type),intent(inout) :: s1    ! a stellar population
    type(stars_type),intent(inout) :: s2    ! an other one
    type(stars_type)               :: s
    
    if (s1%exist) then
      
        if (s1%iMaxAge .le. young_stars_iMaxAge) then
            ! The original stellar population s1 contains only young stars
            ! add directly s1 to s2
            call stars_add(s2,s1,called_by='stars_transfer_young_stars')
            ! erase s1
            call stars_void(s1)
        else
            ! s1 contains young and old stars
            ! young stars are stored in SFH_TAB until young_stars_iMaxAge
            ! create a local copy of the young stellar population 
            call stars_void(s)
            ! allocate the sfh_tab
            allocate(s%sfh_tab(nAgeBins,nMetBins,2))
            ! initialization of the sfh_tab to null values 
            s%sfh_tab(:,:,:) = 0.d0     
            ! copy young stars
            s%sfh_tab(1:young_stars_iMaxAge,:,:) = s1%sfh_tab(1:young_stars_iMaxAge,:,:)
            ! update mass
            s%mass = sum(s%sfh_tab(1:young_stars_iMaxAge,:,1))
            ! set properties
            s%exist   = .true.   ! now the stellar population exists 
            s%iMaxAge = min(young_stars_iMaxAge,s1%iMaxAge)     
        
            ! substract s to s1
            call stars_sub(s1,s,called_by='stars_transfer_young_stars')
            ! @ this point s1 does not contained young stars
        
            if (.not. s2%exist) then
                ! In case of major merger the disc remnent structure is not already created
                call stars_void(s2)
                ! allocate the sfh_tab
                allocate(s2%sfh_tab(nAgeBins,nMetBins,2))
                ! initialization of the sfh_tab to null values 
                s2%sfh_tab(:,:,:) = 0.d0  
                s2%exist = .true.
            end if
            
            ! add s to s2
            call stars_add(s2,s,called_by='stars_transfer_young_stars')
                
            ! erase s
            call stars_void(s)
            call stars_deallocate(s)
        end if ! iMaxAge < young_stars_iMaxAge
    end if
    
    return
  
  end subroutine stars_transfer_young_stars
  
  !*****************************************************************************************************************

  subroutine stars_add(stars1,stars2,called_by) 

    ! ADD A STELLAR POPULATION TO AN OTHER ONE
    
    implicit none
    
    integer(kind=4)                :: iAge, iMet
    
    character(*),optional          :: called_by      ! name of the function which has called this function
    character(MAXPATHSIZE)         :: mes

    type(stars_type),intent(inout) :: stars1         ! a stellar population
    type(stars_type),intent(in)    :: stars2         ! an other one
    
    if (.not.(stars2%exist)) then
       if (stars_mass(stars2) .gt. 0.d0) then
          ! the stellar popualtion stars2 doesn't exist but it has a stellar mass
          call IO_print_error_message('No existing stellar population with no null mass (stars2)', &
               only_rank = rank, called_by = 'stars_add')
          call IO_print_message('with',only_rank=rank,component='stars', &
               param_name = (/'stars2%mass              '/), real_param_val  = (/stars_mass(stars2)/))
          if (present(called_by)) then
             write(mes,'(a,a)') 'called by : ', trim(called_by)
             call IO_print_message(trim(mes),only_rank=rank,component='stars')
          end if
          stop ! stop the program
       end if
       return  ! the stellar population stars2 doesn't exist, no stellar population to add --> exit
    end if

    ! at this point the stellar population stars2 exists and can be added to the stellar popualtion stars1

    if (stars1%exist) then
       ! the stellar population stars1 exists
       stars1%mass     = stars1%mass + stars2%mass                              ! add stellar mass
       stars1%iMaxAge  = max(stars1%iMaxAge,stars2%iMaxAge)                     ! compute max index of sfh_tab 
       where ((stars1%sfh_tab(:,:,1) .gt. 0.d0) .or. (stars2%sfh_tab(:,:,1) .gt. 0.d0)) 
            stars1%sfh_tab(:,:,2) = (stars1%sfh_tab(:,:,1)*stars1%sfh_tab(:,:,2) + & ! mass-weighted evolution time
                                        stars2%sfh_tab(:,:,1)*stars2%sfh_tab(:,:,2)) / &
                                        (stars1%sfh_tab(:,:,1)+stars2%sfh_tab(:,:,1))
       end where
       stars1%sfh_tab(:,:,1) = stars1%sfh_tab(:,:,1) + stars2%sfh_tab(:,:,1)    ! sum sfh_tab array
       if (physical_process) then
          call stars_reset_properties(stars1)  
          do iAge = 1, stars1%iMaxAge
            do iMet = 1, nMetBins
                call stars_compute_properties(stars1,iAge=iAge,iMet=iMet) ! integration
            end do
          end do
          call stars_compute_properties(stars1) ! normalisation
       end if
    else  
       ! stars1 doesn't exist
       if (stars_mass(stars1) .gt. 0.d0) then
          ! but it has a stellar mass
          call IO_print_error_message('No existing stellar population with no null mass (stars1)', &
               only_rank = rank, called_by = 'stars_add')
          call IO_print_message('with',only_rank=rank,component='stars', &
               param_name = (/'stars1%mass              '/), real_param_val  = (/stars_mass(stars1)/))
          if (present(called_by)) then
             write(mes,'(a,a)') 'called by : ', trim(called_by)
             call IO_print_message(trim(mes),only_rank=rank,component='stars')
          end if
          stop ! stop the program
       end if
       ! the stellar population stars1 doesn't exist, the stars_add is therfore equivalent to a stars_copy
       call stars_void(stars1)
       call stars_copy(stars1,stars2)
    end if

    return   
  end subroutine stars_add

  !*****************************************************************************************************************

  subroutine stars_sub(stars1,stars2,called_by) 
    
    ! SUBSTRACT STARS2 to STARS1
    
    implicit none

    integer(kind=4)                :: iAge, iMet
    integer(kind=4),dimension(2)   :: min_pos

    character(*),optional          :: called_by  ! name of the function which has called this function
    character(MAXPATHSIZE)         :: message
    
    type(stars_type),intent(inout) :: stars1     ! a stellar population
    type(stars_type),intent(in)    :: stars2     ! an other one
    
    write(message,'(a,a)') 'stars_sub / ', trim(called_by)
    
    if (.not.(stars2%exist)) then
       ! the stellar population stars2 doesn't exist
       if (stars_mass(stars2) .gt. 0.d0) then
          ! but it has a stellar mass
          call IO_print_error_message('No existing stellar population with no null mass', &
               only_rank = rank, called_by = trim(message))
          call IO_print_message('with',only_rank=rank,component='stars', &
               param_name = (/'stars2%mass              '/), real_param_val  = (/stars_mass(stars2)/))
          stop ! stop the program
       end if
       return  ! there is no stellar popualtion to substract --> exit
    end if

    if (stars1%exist) then
       ! the stellar population exists
       if (stars1%mass .ge. stars2%mass) then
          ! stars1 is more massive than stars2 --> OK
          stars1%mass = stars1%mass - stars2%mass  ! substarct stellar mass
          if (stars1%mass .lt. M_gas_null) then  ! the all stars mass in lower than 0.1 solar mass, which is negligible
             call stars_void(stars1)             ! erase the stellar population and --> exit
             return 
          else
             ! the remanent stellar mass is not negligable 
             stars1%iMaxAge = max(stars1%iMaxAge,stars2%iMaxAge)                     ! compute max index of sfh_tab 
             if (minval(stars1%sfh_tab-stars2%sfh_tab) .lt. 0.d0) then    ! bin with negative mass
                min_pos = minloc(stars1%sfh_tab(:,:,1)-stars2%sfh_tab(:,:,1))
                call IO_print_error_message('min(sfh_tab) < 0.', &
                     only_rank = rank, called_by = trim(message))
                call IO_print_message('with',only_rank=rank,component='stars', &
                     param_name = (/'stars1%sfh_tab(minloc)   ','stars1%sfh_tab(minloc)   ' /), &
                     real_param_val  = (/stars1%sfh_tab(min_pos(1),min_pos(2),1),stars2%sfh_tab(min_pos(1),min_pos(2),1)/))
                stop ! stop the program
             end if
             ! no negative value, we can substract bin mass
             stars1%sfh_tab(:,:,1) = stars1%sfh_tab(:,:,1) - stars2%sfh_tab(:,:,1) ! the bin clock doesn't change
             if (physical_process) then
                call stars_reset_properties(stars1)  
                do iAge = 1, stars1%iMaxAge
                    do iMet = 1, nMetBins
                        call stars_compute_properties(stars1,iAge=iAge,iMet=iMet) ! integration
                    end do
                end do
                call stars_compute_properties(stars1) ! normalisation
            end if
          end if
       else
          call IO_print_error_message('stars1%mass < stars2%mass', &
               only_rank = rank, called_by = trim(message))
          call IO_print_message('with',only_rank=rank,component='stars', &
               param_name = (/'stars1%mass              ','stars2%mass              ' /), real_param_val  = (/stars1%mass,stars2%mass/))
          stop ! stop the program
       end if
    else 
       ! the stellar popualtion stars1 doesn't exist
       if (stars_mass(stars1) .gt. 0.d0) then
          ! but it has a stellar mass
          call IO_print_error_message('No existing stellar population with no null mass', &
               only_rank = rank, called_by = trim(message))
          call IO_print_message('with',only_rank=rank,component='stars', &
               param_name = (/'stars1%mass              '/), real_param_val  = (/stars_mass(stars1)/))
          stop ! stop the program
       end if
       call IO_print_error_message('Attempt to substract a stellar popualtion from an inexistent component', &
               only_rank = rank, called_by = trim(message))
       stop ! stop the program
    end if

    return   
  end subroutine stars_sub

  !*****************************************************************************************************************

  subroutine stars_add_new_stars(stars,gas,gas_return)
    
    ! ADD NEW STARS IN THE STELLAR POPULATION (IN THE FIRST AGE BIN AND IN THE APPROPRIATED METALICITY BIN)
    
    implicit none
    
    integer(kind=4)                :: j,e      ! loop indexes

    real(kind=8)                   :: M_gas    ! mass of gas
    real(kind=8)                   :: Z_gas    ! gas metallicity
    real(kind=8)                   :: M_formed(nMetBins)

    type(stars_type),intent(inout) :: stars      ! the stellar population
    type(gas_type),intent(in)      :: gas        ! the amount of gas which are converted in stars
    type(gas_type),intent(out)     :: gas_return ! the amount of gas that is return to the ISM 
    
#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('stars_add_new_stars',component='stars')
! -------------------------------------------------
#endif
      
    call gas_void(gas_return)         ! initialize the output
    
    ! check the mass of strat forming gas
    ! if there is no gas --> exit
    M_gas = gas_mass(gas)  
    if (M_gas .le. 0.d0) return    

    if (.not. stars%exist) then
       ! the stellar population stars doesn't exist
       if (stars_mass(stars) .gt. 0.d0) then
          ! but it has already stellar mass
          call IO_print_error_message('No existent stellar population with mass', &
               only_rank = rank, called_by = 'stars_add_new_stars')
          call IO_print_message('with',only_rank=rank,component='stars', &
               param_name = (/'stars%mass               '/), real_param_val  = (/stars_mass(stars)/))
          stop ! stop the program
       end if
       if (allocated(stars%sfh_tab)) then
          ! but its sfh_tab is already allocated
          call IO_print_error_message('No existent stellar population with already allocated sfh_tab', &
               only_rank = rank, called_by = 'stars_add_new_stars')
          call IO_print_message('with',only_rank=rank,component='stars', &
               param_name = (/'stars%mass               '/), real_param_val  = (/stars_mass(stars)/))
          call stars_void(stars)       ! we erase the pre-existing stellar popualtion
       end if
       ! allocate the sfh_tab
       allocate(stars%sfh_tab(nAgeBins,nMetBins,2))
       stars%sfh_tab(:,:,:) = 0.d0     ! initialization of the sfh_tab to null values 
       stars%exist          = .true.   ! now the stellar population exists 
       stars%iMaxAge        = 1        ! init
    end if
    !
    ! initialize gas_return and M_formed
    gas_return  = gas
    Z_gas       = gas_metallicity(gas_return)
    M_formed(:) = 0.d0
    !
    if (Z_gas .lt. MetBins(2)) then
        M_formed(1) = M_gas
        call gas_void(gas_return)
    else
        ! run onto metallicity bins to build the optimal combination of stellar population
        do j = nMetBins, 1, -1
           !
           ! compute the mass of available gas 
           M_gas = gas_mass(gas_return) 
           ! compute the metallicity of available gas 
           Z_gas = gas_metallicity(gas_return)
           !
           if (Z_gas .ge. MetBins(j)) then
              ! test the metal mass
              M_formed(j) = min(M_gas, gas_return%mZ/MetBins(j))   
              do e = 1, nElts
                 if (InitAbund(j)%elts(e) .gt. 0.d0) then
                    M_formed(j) = min(M_formed(j), gas_return%elts(e)/InitAbund(j)%elts(e))
                    if (abs(M_formed(j)) .lt. M_gas_null) M_formed(j) = 0.d0
                 end if
              end do  
              ! Remove the gas used to formed the stellar popualtion
              call gas_sub(gas_return,InitAbund(j)*M_formed(j),called_by='gas_return/M_formed')
              if (gas_mass(gas_return) .le. M_gas_null) then
                call gas_void(gas_return)
                exit
              end if
           end if
       end do
    end if  
    !
    ! increment the stellar mass, the gas is converted in stars
    stars%mass = stars%mass + sum(M_formed) 
    ! Create a stellar popualtion 
    stars%sfh_tab(1,:,1) = stars%sfh_tab(1,:,1) + M_formed
    
#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('stars_add_new_stars ... done',component='stars')
! -------------------------------------------------
#endif
    
    return    
  end subroutine stars_add_new_stars

  !*****************************************************************************************************************

  subroutine stars_evolve(stars,dt,loss_rate,sfr,gas_return)

    ! EVOLVE A STELLAR POPULATION : 
    ! TAKES INTO ACCOUNT : - LOSS RATES
    !                      - STAR FORMATION RATE 

    implicit none

    integer(kind=4)                     :: iAge,iAge2,iMaxAge   ! dumb index to loop over age
    integer(kind=4)                     :: iMet                 ! dumb index to loop over metallicity

    real(kind=8),intent(in)             :: dt                   ! time interval (Gyr)
    real(kind=8)                        :: ddt1, ddt2           ! sub time-step ddt1 + ddt2 = dt
    real(kind=8)                        :: m_bin, dm_bin        ! mass, mass loss in an age and metallicity bin
    real(kind=8)                        :: t_bin                ! clock in a gien bin
    real(kind=8)                        :: mass_save            ! mass of the stellar population before the evolution loop
    real(kind=8)                        :: loss_rate_prior      ! mass loss rate used as prior in the ddt computation
    real(kind=8)                        :: SN_rate_save         ! previous SN rate
    real(kind=8)                        :: err                  ! relative error between mass and sum(sfh_tab)
    real(kind=8)                        :: U, dU                ! mass evolution test

    type(stars_type),intent(inout)      :: stars                ! stellar population to evolve
    type(gas_type),intent(in),optional  :: sfr                  ! star formation rate
    type(gas_type),intent(out),optional :: gas_return           ! 
    type(gas_type),intent(out)          :: loss_rate            ! stellar instantaneous loss rate
    type(gas_type)                      :: dloss_rate

#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('stars_evolve',only_rank=rank,component='stars')
! -------------------------------------------------
#endif
    
    mass_save = stars_mass(stars)               ! save 
    loss_rate_prior = gas_mass(stars%loss_rate) ! save
    SN_rate_save = stars%SN_rate
                             
    if (present(gas_return)) then
        call gas_void(gas_return)      ! init
    end if
    call gas_void(dloss_rate)          ! init
    call gas_void(loss_rate)           ! init
    call stars_reset_properties(stars) ! init

    if (stars%exist) then
       !
       !*********************************************************
       ! STARS EXIST 
       !*********************************************************
       !
       iMaxAge = min(stars%iMaxAge,nAgeBins) ! fix the end of the age loop 
       !
       ! evolution loop over Age and metalicity
       do iAge = iMaxAge, 0, -1  
          do iMet = nMetBins, 1, -1
             !         
             if (iAge .ge. 1) then    
                !       
                ! take into account stellar evolution
                m_bin = stars%sfh_tab(iAge,iMet,1)
                !
                if (m_bin .gt. 0.d0) then ! if there are stars in the bin
                      !
                      t_bin = stars%sfh_tab(iAge,iMet,2)  
                      !     
                      if (t_bin + dt .gt. AgeBins(iAge+1)) then
                         ! The mass of the current AgeBin must be move to the next AgeBin
                         ddt1 = AgeBins(iAge+1) - t_bin         ! time until the end of the current MetBin
                         ddt2 = (t_bin + dt) - AgeBins(iAge+1)  ! time in the next AgeBin 
                         ! stellar ejecta rate
                         ! compute dloss_rate (time-average)
                         !dloss_rate = m_bin*(mass_loss_rates(iAge,iMet)*ddt1 + mass_loss_rates(iAge+1,iMet)*ddt2)*(1.d0/(ddt1 + ddt2))
                         dloss_rate = m_bin*mass_loss_rates(iAge,iMet)
                         ! update loss rate
                         loss_rate = loss_rate + dloss_rate   
                         ! compute dm_bin
                         dm_bin = gas_mass(dloss_rate)*dt    
                         ! Remove the ejected mass from the current AgeBin
                         if (dm_bin .gt. m_bin) then
                            call IO_print_error_message('Substract to much mass to a sfh_tab bin ! ', &
                                only_rank = rank, called_by = 'stars_evolve')
                            call IO_print_message('with',only_rank=rank,component='stars', &
                                param_name = (/'iAge                     ','iMet                     '/), int_param_val = (/iAge,iMet/))
                            call IO_print_message('with',only_rank=rank,component='stars', &
                                param_name = (/'dm_bin                   ','m_bin                    '/), real_param_val = (/dm_bin, m_bin/))
                            stop ! stop the program
                         end if
                         !
                         ! Remove ejected mass from the bin
                         stars%sfh_tab(iAge,iMet,1) = stars%sfh_tab(iAge,iMet,1) - dm_bin   ! remove from the age/metalicity bin
                         stars%mass                 = stars%mass                 - dm_bin   ! remove from the total stellar mass
                         ! clock
                         stars%sfh_tab(iAge+1,iMet,2) = (stars%sfh_tab(iAge+1,iMet,1)*stars%sfh_tab(iAge+1,iMet,2) + stars%sfh_tab(iAge,iMet,1)*(t_bin + dt)) / &
                              (stars%sfh_tab(iAge+1,iMet,1) + stars%sfh_tab(iAge,iMet,1))
                         ! add mass to the next AgeBin
                         stars%sfh_tab(iAge+1,iMet,1) = stars%sfh_tab(iAge+1,iMet,1) + stars%sfh_tab(iAge,iMet,1)
                         ! remove mass and reset the clock from the current AgeBins
                         stars%sfh_tab(iAge,iMet,:) = 0.d0
                         !
                         ! update iMaxAge
                         stars%iMaxAge = max(stars%iMaxAge,min(nAgeBins,iAge+1))
                      else
                         ! no mass transfer have to be done 
                         ! stellar ejecta rate
                         dloss_rate = m_bin*mass_loss_rates(iAge,iMet)
                         ! update loss rate
                         loss_rate = loss_rate + dloss_rate   
                         ! compute dm_bin
                         dm_bin = gas_mass(dloss_rate)*dt  
                         ! Remove the ejected mass from the current AgeBin
                         if (dm_bin .gt. m_bin) then
                            call IO_print_error_message('Substract to much mass to a sfh_tab bin ! ', &
                                only_rank = rank, called_by = 'stars_evolve')
                            call IO_print_message('with',only_rank=rank,component='stars', &
                                param_name = (/'iAge                     ','iMet                     '/), int_param_val = (/iAge,iMet/))
                            call IO_print_message('with',only_rank=rank,component='stars', &
                                param_name = (/'dm_bin                   ','m_bin                    '/), real_param_val = (/dm_bin, m_bin/))
                            stop ! stop the program
                         end if
                         !
                         ! Remove ejected mass from the bin
                         stars%sfh_tab(iAge,iMet,1) = stars%sfh_tab(iAge,iMet,1) - dm_bin   ! remove from the age/metalicity bin
                         stars%mass                 = stars%mass                 - dm_bin   ! remove from the total stellar mass
                         ! update clock
                         stars%sfh_tab(iAge,iMet,2) = stars%sfh_tab(iAge,iMet,2) + dt
                      endif
                      !
                      ! check very low mass bin
                      ! recompute m_bin and t_bin for the current AgeBin
                      m_bin = stars%sfh_tab(iAge,iMet,1)
                      t_bin = stars%sfh_tab(iAge,iMet,2)
                      if ((m_bin .lt. M_gas_null) .and. (m_bin .ne. 0.d0)) then
                         if (abs(m_bin) .lt. M_gas_null) then
                            stars%sfh_tab(iAge,iMet,1) = 0.d0
                            stars%sfh_tab(iAge,iMet,2) = 0.d0
                            ! remove from the total mass
                            stars%mass = stars%mass - max(0.d0,m_bin)
                         else
                            ! If the next larger AgeBin contains mass, we transfert the residual mass of the current AgeBin
                            ! otherwise, this residual mass is kept in the current mass bin 
                            if (stars%sfh_tab(iAge+1,iMet,1) .gt. 0.d0) then
                                ! clock
                                stars%sfh_tab(iAge+1,iMet,2) = (stars%sfh_tab(iAge+1,iMet,1)*stars%sfh_tab(iAge+1,iMet,2) + m_bin*t_bin) / &
                                    (stars%sfh_tab(iAge+1,iMet,1) + m_bin)
                                ! transfert mass
                                stars%sfh_tab(iAge+1,iMet,1) = stars%sfh_tab(iAge+1,iMet,1) + m_bin
                                ! reset mass and clock
                                stars%sfh_tab(iAge,iMet,:) = 0.d0
                            end if
                         endif
                      end if
                end if ! if m_bin > 0
             end if  ! iAge .ge. 1
             ! 
             ! a run over all sfh_tab need time, to optimize the computation time we recompute stellar properties just after the uptade of the sfh_tab. 
             ! we use the same loop, but with +1 in iAge loop for taking into account stellar evolution, we win time. 
             !
             iAge2 = iAge + 1
             if ((iAge2 .ge. 1) .and. (iAge2 .le. nAgeBins)) then
                ! we cannot recompute ejecta of the first bin age before add new stars
                if (iAge2 .eq. 1) then
                   if (iMet .eq. nMetBins) then
                      ! it's the end of the age loop and the begining of the metal loop
                      ! evolution of all bins has already been done
                      ! we must add new stars and compute the new end life lost rate of the first age bin
                      if (present(sfr) .and. present(gas_return)) then
                         call stars_add_new_stars(stars,sfr*dt,gas_return)  ! create new stars 
                      end if
                   end if
                end if
                call stars_compute_properties(stars,iAge=iAge2,iMet=iMet)   ! integration loop
             end if
          end do ! iMet
       end do ! iAge
       !
       ! end of Age and metallicity loops
       call stars_compute_properties(stars) ! normalization
       if (SN_rate_save .gt. 0.d0) then
           stars%SN_rate = 5.d-1*(SN_rate_save+stars%SN_rate)
       end if
    else
       !
       ! *************************************************************************
       ! STARS DOESN'T EXIST YET
       ! *************************************************************************
       !
       ! add new stars
       if (present(sfr) .and. present(gas_return)) then
          call stars_add_new_stars(stars,sfr*dt,gas_return)  ! create new stars 
          do iMet = 1, nMetBins
            call stars_compute_properties(stars,iAge=1,iMet=iMet)
          end do
          call stars_compute_properties(stars)
       end if
    end if ! stars%exist
    !
    ! **********************************************
    ! TESTS
    ! **********************************************
    !
    if (stars%exist) then
       if (stars_mass(stars) .gt. 0.d0) then
          !
          ! test consistency between stars%mass and sum(sfh_tab)
          err = stars_mass(stars) - sum(stars%sfh_tab(:,:,1))
          if (err .gt. 1.d1*M_gas_null) then
             call IO_print_warning_message('sum(sfh_tab) /= stars%mass', &
                  only_rank = rank, called_by = 'stars_evolve')
             !     
#ifdef PRINT_WARNING
! -------------------------------------------------                  
             call IO_print_message('with',only_rank=rank,component='stars', &
                  param_name = (/'stars%mass               ','sum(sfh_tab)             ', &
                  'dt                       ','err                      ','mass_save                '/), &
                  real_param_val  = (/stars_mass(stars),sum(stars%sfh_tab(:,:,1)),dt,err,mass_save/))
! -------------------------------------------------
#endif
! PRINT_WARNING
             !
             stars%mass = sum(stars%sfh_tab(:,:,1)) ! reset
          endif
          !
          err = gas_mass(loss_rate)-loss_rate_prior
          if (abs(err) .gt. num_accuracy) then
             call IO_print_warning_message('loss_rate /= loss_rate_prior', &
                  only_rank = rank, called_by = 'stars_evolve')
             !     
#ifdef PRINT_WARNING
! -------------------------------------------------                   
            call IO_print_message('with',only_rank=rank,component='stars', &
                  param_name = (/'loss_rate                ','loss_rate_prior          ','err                      ',&
                                 'stars%mass               ','mass_save                '/), &
                  real_param_val  = (/gas_mass(loss_rate),loss_rate_prior,gas_mass(loss_rate)-loss_rate_prior, &
                                      stars_mass(stars),mass_save/))
! -------------------------------------------------
#endif
! PRINT_WARNING
             !
          end if
          !       
          ! check positivity of the ejecta mass
          if (gas_mass(stars%loss_rate) .lt. 0.d0) then
              call IO_print_error_message('stars%loss_rate < 0', &
                    only_rank = rank, called_by = 'stars_evolve')
              call IO_print_message('with',only_rank=rank,component='stars', param_name = (/'loss_rate                '/),real_param_val  = (/gas_mass(stars%loss_rate)/))
              stop ! stop the program
          end if
          !
          ! test stellar mass variation
          if (mass_save .gt. M_stars_min) then
            U  = mass_save
            dU = abs(stars_mass(stars)-mass_save)
            if ((abs(dU)/U) .gt. physical_precision) then
                call IO_print_error_message('No quasi-static state of the stellar component',only_rank=rank,called_by='stars_evolve')
                call IO_print_message('with',only_rank=rank,component='stars', &
                      param_name = (/'dt                       ','dU/U (%)                 ','U                        ','sfr                      ',&
                                     'loss_rate_prior          ','loss_rate                '/), &
                      real_param_val = (/dt,1.d2*dU/U,U,gas_mass(sfr),loss_rate_prior,gas_mass(loss_rate)/))
                stop  ! stop the program
            end if
          end if
       end if
    end if
    
#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('stars_evolve ... done',only_rank=rank,component='stars')
! -------------------------------------------------
#endif
    
    return
  end subroutine stars_evolve

  !*****************************************************************************************************************

  subroutine stars_build_stellar_spectra(stars,spectra)

    ! BUILD THE MASS-WEIGHTED STELLAR SPECTRA OF A GIVEN STELLAR POPULATION
    ! We also compte the mass weighted ISRF

    implicit none

    integer(kind=4)                        :: i,j            ! loop indexes over Ages and Metalicity bins

    real(kind=8)                           :: m_bin, t_bin   ! mass in a bin, mass average age of the bin
    real(kind=8)                           :: dt1, dt2
    real(kind=8)                           :: min_test
    real(kind=4),intent(out),allocatable   :: spectra(:,:)   ! young stars and old stars spectra
    real(kind=4),allocatable               :: tmp_spectra(:) ! the age-average stellar spectra
    
    type(stars_type),intent(inout)         :: stars          ! a stellar population

    ! create spectra
    allocate(spectra(nWaves,2))   ! spectra[:,1] spectrum of the young stars
                                  ! spectra[:,2] spectrum of the old stars
    spectra(:,:) = 0.d0           ! init spectra
    ! create the tmp_spectra
    allocate(tmp_spectra(nWaves))
    tmp_spectra(:) = 0.d0 ! init
    
    stars%ISRF = 0.d0   ! init
     
    if (stars%exist) then
      !
      do i = 1, stars%iMaxAge
        ! Age loop
        do j = 1, nMetBins
          ! Met loop
          ! Compute the mass of the stars with age AgeBins(i) and metallicity MetBins(j)
          m_bin = stars%sfh_tab(i,j,1)
          if (m_bin .gt. 0.d0) then
              t_bin = stars%sfh_tab(i,j,2) ! mass average age of the current bin
              ! 
              if (t_bin .gt. AgeBins(nAgeBins)) then
                 ! t_bin > AgeBins(nAgeBins)
                 ! sb_sed(:,i,j) is in [Lsun] and normalized to M*_liv + M_rem = 1Msun
                 tmp_spectra = sb_sed(:,nAgeBins,j)  
              else
                 ! AgeBins(i) < t_bin < AgeBins(i+1)
                 dt1 = AgeBins(i+1) - t_bin
                 dt2 = t_bin - AgeBins(i)
                 ! compute age-average spectra 
                 ! sb_sed(:,i,j) is in [Lsun/Msun] and normalized to M*_liv + M_rem = 1Msun
                 tmp_spectra = (real(dt1,4)*sb_sed(:,i,j)+real(dt2,4)*sb_sed(:,min(nAgeBins,i+1),j))/real((dt1 + dt2),4) 
              end if
              !
              ! Add their contribution to the stellar spectrum
              if (i .le. young_stars_iMaxAge) then
                !
                ! young stars
                spectra(:,1) = spectra(:,1) + real(m_bin*mass_code_unit_in_M_Sun,4)*tmp_spectra  ! in [Lsun]
              else
                !
                ! old stars              
                spectra(:,2) = spectra(:,2) + real(m_bin*mass_code_unit_in_M_Sun,4)*tmp_spectra  ! in [Lsun]
              end if
              !
           end if
        end do
      end do
      !
      min_test = minval(spectra(:,1))
      if (min_test .lt. 0.d0) then
            call IO_print_warning_message('Stellar spectra, negative value [young stars]',only_rank = rank, called_by = 'stars_build_stellar_spectra')
            call IO_print_message('with',only_rank=rank,component='gal', &
                      param_name = (/'min_test                 '/), real_param_val = (/min_test/))
            where (spectra(:,1) .lt. 0.d0)  
                spectra(:,1) = 0.d0
            endwhere      
      end if
      !
      min_test = minval(spectra(:,2))
      if (min_test .lt. 0.d0) then
            call IO_print_warning_message('Stellar spectra, negative value [old stars]',only_rank = rank, called_by = 'stars_build_stellar_spectra')
            call IO_print_message('with',only_rank=rank,component='gal', &
                      param_name = (/'min_test                 '/), real_param_val = (/min_test/))
            where (spectra(:,2) .lt. 0.d0)  
                spectra(:,2) = 0.d0
            endwhere          
      end if
      !   
      ! set ISRF [Lsun] 
      stars%ISRF(1) = trap(spectra(i_14eV:i_6eV,1)/Waves(i_14eV:i_6eV),Waves(i_14eV:i_6eV))  ! young stars [Lsun]
      stars%ISRF(2) = trap(spectra(i_14eV:i_6eV,2)/Waves(i_14eV:i_6eV),Waves(i_14eV:i_6eV))  ! old stars   [Lsun]
      !
    endif

    return
  end subroutine stars_build_stellar_spectra

  !*****************************************************************************************************************

  subroutine stars_deallocate(s)

    ! DEALLOCATE THE SFH_TAB ARRAY OF A GIVEN STELLAR POPULATION

    implicit none 

    type(stars_type),intent(inout) :: s
    
    if (allocated(s%sfh_tab)) deallocate(s%sfh_tab)
    if (s%exist) s%exist = .false.

    return
  end subroutine stars_deallocate

  !*****************************************************************************************************************
  !
  ! FUNCTIONS
  !
  !*****************************************************************************************************************

  function stars_mass(s)

    ! RETURN THE STELLAR MASS OF A STARS OBJECT

    implicit none 

    real(kind=8)                :: stars_mass
    
    type(stars_type),intent(in) :: s
    
    stars_mass = s%mass
    
    return
  end function stars_mass
  
  !*****************************************************************************************************************
  
  function stars_dt_optim(s,input)
  
    ! RETURN THE OPTIMAL EVOLUTION TIMESTEP ASSOCIATED TO A GIVEN STELLAR POPULATION
    
    implicit none
    
    real(kind=8)                :: stars_dt_optim
    real(kind=8)                :: M, rate
    
    type(gas_type),intent(in)   :: input    ! the global input rate (here SFR)
    
    type(stars_type),intent(in) :: s
    
    stars_dt_optim = -1.d0 ! init
    
    if (s%exist) then
        M    = stars_mass(s)
        rate = abs(gas_mass(input) - gas_mass(s%loss_rate))  
        if ((M .gt. M_stars_min) .and. (rate .gt. 0.d0)) then
            ! compute optimal time-step for stellar population
            stars_dt_optim = min(StellarTimeStep,5.d-1*physical_precision*M/abs(rate))
        else
            stars_dt_optim = StellarTimeStep
        end if
    end if
    
    return
  end function stars_dt_optim

  !***********************************************************************************************************
  !
  ! PRINTING PROCEDURES
  !
  !***********************************************************************************************************

  subroutine stars_load_stars_data(stars,stars_data)

    ! CREATE THE STARS OUTPUT PROPERTIES LIST

    implicit none

    real(kind=8),intent(inout)   :: stars_data(nb_stars_field)

    type(stars_type),intent(in)  :: stars  ! disc component
    
    ! For information
    !'stars_mass       ','stars_f_young    ', &
    !'stars_age_MW     ','stars_dage_MW    ', &
    !'stars_age_LW     ','stars_dage_LW    ', &
    !'stars_Z_MW       ','stars_Z_LW       '/)

    stars_data = (/mass_code_unit_in_M_Sun*stars%mass, stars%f_young, stars%Age(1), stars%dAge(1), stars%Age(2), stars%dAge(2), &
        stars%Z(1)/Z_sun, stars%Z(2)/Z_sun/)
 
    return
  end subroutine stars_load_stars_data

  !***********************************************************************************************************
  !
  ! INTERFACE OPERATOR DEFINITIONS
  !
  !***********************************************************************************************************

  function stars_add_(s1,s2)
    
    ! RETURN g1 + g2

    implicit none 

    type(stars_type),intent(in)  :: s1, s2
    type(stars_type)             :: stars_add_

    call stars_void(stars_add_)
    call stars_copy(stars_add_,s1)
    call stars_add(stars_add_,s2)
    
    return
  end function stars_add_

  !***********************************************************************************************************

  function stars_sub_(s1,s2)
    
    ! RETURN if possible g1 - g2

    implicit none

    type(stars_type),intent(in)  :: s1, s2
    type(stars_type)             :: stars_sub_

    call stars_void(stars_sub_)
    call stars_copy(stars_sub_,s1)
    call stars_sub(stars_sub_,s2)
    
    return
  end function stars_sub_

  !***********************************************************************************************************

  function stars_scalar_multiply_left(stars,a) 

    ! RETURN a * STARS (ALL MASSES MULTIPLIED)

    implicit none

    type(stars_type),intent(in)  :: stars                      ! a stellar population
    real(kind=8),intent(in)      :: a                          ! scallar coefficient
    type(stars_type)             :: stars_scalar_multiply_left ! output stellar population = a * stars
   
    call stars_void(stars_scalar_multiply_left)                ! init with null values
    
    stars_scalar_multiply_left%mass      = stars%mass * a          
    stars_scalar_multiply_left%Z         = stars%Z   
    stars_scalar_multiply_left%dZ        = stars%dZ
    stars_scalar_multiply_left%Age       = stars%Age 
    stars_scalar_multiply_left%dAge      = stars%dAge 
    stars_scalar_multiply_left%iMaxAge   = stars%iMaxAge  
    stars_scalar_multiply_left%loss_rate = stars%loss_rate*a

    if (stars%exist) then
       allocate(stars_scalar_multiply_left%sfh_tab(nAgeBins,nMetBins,2))    ! create the sfh_tab
       stars_scalar_multiply_left%sfh_tab(:,:,:) = 0.d0                     ! initialization of the sfh_tab 
       stars_scalar_multiply_left%sfh_tab(:,:,1) = stars%sfh_tab(:,:,1) * a ! copy sfh_tab
       stars_scalar_multiply_left%sfh_tab(:,:,2) = stars%sfh_tab(:,:,2)     ! copy sfh_tab 
       if (minval(stars_scalar_multiply_left%sfh_tab(:,:,1)) .lt. 0.d0) then
          call IO_print_error_message('min(sfh_tab) < 0',only_rank=rank,called_by='stars_scalar_multiply_left')
          call IO_print_message('with',only_rank=rank,component='stars', &
                  param_name = (/'stars%mass               ','sum(sfh_tab)             ','min(sfh_tab(:,:,1)       '/), &
                  real_param_val  = (/stars_mass(stars),sum(stars%sfh_tab(:,:,1)),minval(stars%sfh_tab(:,:,1))/))
          stop ! stop the program
       end if
       stars_scalar_multiply_left%exist = .true.  ! the stellar population exists     
    else  
       stars_scalar_multiply_left%exist = .false. ! the stellar population doesn't exist
    end if

    return
  end function stars_scalar_multiply_left

  !***********************************************************************************************************

  function stars_scalar_multiply_right(a,stars) 
    
    ! RETURN a * STARS (ALL MASSES MULTIPLIED)

    implicit none

    real(kind=8),intent(in)      :: a                            ! scallar coefficient

    type(stars_type),intent(in)  :: stars                        ! a stellar popualtion     
    type(stars_type)             :: stars_scalar_multiply_right  ! output stellar population = a * stars
   
    call stars_void(stars_scalar_multiply_right)
    call stars_copy(stars_scalar_multiply_right,stars_scalar_multiply_left(stars,a)) 

    return
  end function stars_scalar_multiply_right

  !***********************************************************************************************************

end module stars
