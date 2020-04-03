module tree

  use halo             ! Contains halo structure and halo exploitation functions
 
  public
  
  !*****************************************************************************************************************
  ! 
  ! OVERVIEW
  !
  ! tree_module defines the univ and the tree data structures: univ(ts) and univ(ts)%tree
  !
  ! univ(ts) contains the informations about merger trees evolution at a given timestep (ts)
  !
  ! univ(ts)%tree contains all links existing between interacting halos through their evolution
  !
  ! The tree module contains also the definitions for: 
  !  - tree_type
  !  - treeIDs_type
  !  - timestep_type
  ! 
  ! The tree module stores all evolution procedures allowing to compute halo evolutions through the cosmic time
  ! In the header of the module are defined all tree output properties (labels, units and formats)
  !
  !  MPI PROCESSES IN THIS MODULE
  !
  !  tree_send_z_data                                    : send a list of global (cosmic time) average quantity associated 
  !                                                           to galaxy and halo scale (average accretion rate, SFR, ... )
  !
  !  tree_send_data                                      : send specific informations about a halo in a tree (HID, host_HID, level ...)
  !
  !  tree_receive_data                                   : receive specific informations about a halo in a tree 
  !
  ! SUBROUTINES IN THIS MODULE
  !
  ! tree_set_reference_mass                              : Call procedures to set reference masses
  !
  ! tree_read                                            : Read dark-matter (dm) properties and halo links extracted from the Nbody dark-matter simulation
  !     called by main program
  !     contains : tree_read_timestep_properties         : Read timestep properties of the tree in tree_file
  !                tree_read_tree_structure              : Read tree structure in tree_file
  !                tree_read_halo_properties             : Read dm properties of halo contained in the tree in props_files 
  !
  !  tree_compute_halo_derived_properties                : Compute derived properties of dm component which 
  !     called by : main program                             can't be computed directly in the Nbody simulation like (integreted mass, accretion rate ...)
  !      
  !  tree_clean                                          : Flag haloes in function of their dm-properties and their position in the tree (good or bad) branch
  !     called by : tree_compute_halo_derived_properties
  !     contains :  tree_compute_dm_goodness_of_branch   
  !                 tree_set_dm_goodness_of_branch      
  !                 tree_halo_become_good                
  !                 tree_compute_halo_goodness_of_branch 
  !                 tree_set_halo_goodness_of_branch                                       
  !                 tree_interpolate_bad_properties  
  !
  !  tree_kill_branch                                    : Set halo%compute of all propenitors of a given halo to .false.
  ! 
  !  tree_followed_branch                                : Select followed haloes in a merger tree
  !
  !  tree_deallocate                                     : erase a merger tree structure
  !    called by   : tree_evolve
  !
  !  tree_evolve                                         : evolve baryonic phase of the tree
  !    called by   : main program                              each halo evolve and merge 
  !    contains    : tree_evolve_halo                    : evolve a given halo during dt
  !                  last_computable_progenitor_of_a_main_branch 
  !
  !  tree_spectrum                                       : receive univ, tree and halo informations (HID, HID_host, level, aexp) 
  !    called by   : main program                              and calls higher protocols dedicated to spectrum builder
  !
  !  FUNCTIONS IN THIS MODULE
  ! 
  !   tree_most_massive_prog_index                       : Return the index of the most massive progenitor of a halo in a merger tree
  !
  !  PRINTING PROCEDURES
  !
  !  tree_load_tree_data                                 : load tree data
  !    called by : tree print procedure (e.g. tree_print, tree_follow_up_print)
  !
  !  tree_print                                          : print tree component properties   
  !    called by   : tree_evolve                       
  ! 
  !  tree_print_tree                                     : print halo properties over a tree (taking account all computable progenitors)
  !    called by   : tree_evolve
  ! 
  !*****************************************************************************************************************
  
  ! DEFINITION OF GLOBAL VARIABLES LINKED TO THE TREE STRUCTURE *******************
  
  ! TREE_TYPE DEFINITION **************
  ! A tree_type object stores the indexes that link a halo to its descendents, progenitors, hosts and subhaloes
  ! these indexes work locally within a tree file and differ from the unique IDs that identify haloes in a simulation
  
  type tree_type
    integer(kind=8) :: HID_orig        ! The original Halo identification number (TreeMaker, Rockstar)
    integer(kind=8) :: HID             ! Halo identification number (G.A.S.)
    integer(kind=8) :: host_HID        ! Holo identification number of the host halo (in case of substructure)
    integer(kind=4) :: me              ! my number (index at timestep ts)
    integer(kind=4) :: descendent      ! index of my descendent (index at timestep ts +1)
    integer(kind=4) :: ndads           ! number of progenitors
    integer(kind=4) :: n_computed_dads ! number of computable progenitors
    integer(kind=4) :: firstProg       ! index of my first progenitor (firstProg = -1 when there are no progenitors) (index at timestep ts -1)
    integer(kind=4) :: nextProg        ! index of the next progenitor of my descendent (my brother next on the line, -1 if not) 
    integer(kind=4) :: hosthalo        ! index of main halo
    integer(kind=4) :: hostsub         ! index of host subhalo
    integer(kind=4) :: nextsub         ! index of next substructure in main halo
    integer(kind=4) :: level           ! level of substructure (level = 1 for a main halo)
  end type tree_type
    
  ! TIMESTEP_TYPE DEFINITION ********** 
  ! A timestep_type object is a data structure that contains the entire information on the universe at a given timestep
  ! There is a one-to-one link between the tree structure and the halo structure 
  ! i.e. tree(ih) contains the tree information for halo(ih)
  
  type timestep_type
    integer(kind=4)                :: ts                      ! the timestep index
    integer(kind=4)                :: nb_of_halos             ! nb of haloes in a tree file at a given ts (includes subhaloes)
    integer(kind=4)                :: nb_of_computable_halos  ! nb of computable (compute = .true.) haloes in a tree file at a given ts (includes subhaloes)
                                                              ! In case of the cleaning algorithm is applied, only the number of computable halos (and sub-halos) is given 
    integer(kind=4)                :: nb_of_followed_up_halos ! In follow mode only a selection of halos are computed                       
    real(kind=4)                   :: age_univ                ! age of the universe (in Gyr) at this time-step
    real(kind=4)                   :: aexp                    ! expansion factor at this time-step (=1 at z=0) 
    
    type(tree_type),allocatable    :: tree(:)                 ! tree table: contains the list of tree-type structure; one for halo identified at the given time-step
    type(halo_type),allocatable    :: halo(:)                 ! halo table: contains the list of halo-type structures, a halo type structure groups : 
                                                              ! - the dm-structure, 
                                                              ! - the baryon-halo-structure and the galay-structure for each halo identified in the dark-matter simulation
  end type timestep_type

  ! CREATE univ(:) *******************
  type(timestep_type),allocatable  :: univ(:)                 ! univ(ts) contains the information (list of timestep_type structures) on the entire simulated universe at timestep ts
                                                              ! ts goes from 0 to nsteps but ts = 0 has no haloes
  type(timestep_type)              :: u                       ! A local tiny copy a the univ structure, used to build galaxy spectra                                                            

  ! DEFINE HEADER INFORMATIONS *********************
                                                             
  ! hdu reference for tree structure
  integer(kind=4)                  :: hdu_tree                                                            
  ! printable properties for tree structure
  integer(kind=4),parameter :: nb_tree_field = 6              ! Number of tree properties saved
  ! Name of each output colomn data
  character(len=ttype_len),dimension(nb_tree_field) :: ttype_tree = (/'redshift              ','ts                    ','HID                   ',&
                                                                      'host_HID              ','level                 ','n_dads                '/)  
  ! Physical unit of each column data
  character(len=tunit_len),dimension(nb_tree_field) :: tunit_tree = (/'w_o_unit    ','#           ','#           ',&
                                                                      '#           ','#           ','#           '/)
  ! Data type of each column data
  character(len=tform_len),dimension(nb_tree_field) :: tform_tree = (/'1E','1J','1K','1K','1I','1J'/)
  
  !
  ! z table
  integer(kind=4),parameter                    :: n_z_data = 12     
  ! legend for the z table 
  character(len=ttype_len),dimension(n_z_data) :: z_table_legend = (/'ts                    ','redshift              ','Age_Universe          ',&
                                                                     'merger rate           ','dm-accretion rate     ','Halo virial mass      ',&
                                                                     'bar-accretion rate    ','cold-accretion rate   ','cooling rate          ',&
                                                                     'SFR                   ','Stellar mass (disc)   ','Stellar mass (bulge)  '/)  
  character(len=ttype_len),dimension(n_z_data) :: z_table_legend_unit   = (/'#                     ','wo-u                  ','Gyr                   ',&
                                                                            'Msun/yr/Mpc^3         ','Msun/yr/Mpc^3         ','Msun/Mpc^3            ',&
                                                                            'Msun/yr/Mpc^3         ','Msun/yr/Mpc^3         ','Msun/yr/Mpc^3         ',&
                                                                            'Msun/yr/Mpc^3         ','Msun/Mpc^3            ','Msun/Mpc^3            '/)  

contains

  !*****************************************************************************************************************
  ! 
  ! MPI PROCESSES
  !
  !*****************************************************************************************************************
  
  subroutine tree_send_z_data(u)
  
    implicit none
    
    real(kind=8),allocatable       :: buffer(:)
    
    type(timestep_type),intent(in) :: u         ! a local copy of the universe
    
    ! For information
    !'ts                    ','redshift              ','Age_Universe          '
    !'merger rate           ','dm-accretion rate     ','Halo virial mass      '
    !'bar-accretion rate    ','cold-accretion rate   ','cooling rate          '
    !'SFR                   ','Stellar mass (disc)   ','Stellar mass (bulge)  '

    !
    ! send the ts
    call MPI_SEND(u%ts,1,MPI_INTEGER4,main_process_rank,ztable_tag,MPI_COMM_WORLD,ierror)
    !
    ! send the number of halos
    call MPI_SEND(u%nb_of_computable_halos,1,MPI_INTEGER4,main_process_rank,ztable_tag+1,MPI_COMM_WORLD,ierror)
    !
    if (u%nb_of_computable_halos .gt. 0) then
        !
        allocate(buffer(n_z_data-1))
        !   
        buffer = (/1.d0/u%aexp - 1.d0, real(u%age_univ,8),& 
                mass_rate_code_unit_2_MsunPerYr*sum(u%halo(:)%dm%merger_rate,MASK=(u%halo(:)%compute)),&
                mass_rate_code_unit_2_MsunPerYr*sum(u%halo(:)%dm%acc_rate_left,MASK=(u%halo(:)%compute)),&
                mass_code_unit_in_M_Sun*sum(u%halo(:)%dm%M_vir,MASK=(u%halo(:)%compute)),&
                mass_rate_code_unit_2_MsunPerYr*(sum(u%halo(:)%baryon_halo%cold_inflow_rate%mass,MASK=(u%halo(:)%compute))+sum(u%halo(:)%baryon_halo%hot_inflow_rate%mass,MASK=(u%halo(:)%compute))),&
                mass_rate_code_unit_2_MsunPerYr*sum(u%halo(:)%baryon_halo%cold_outflow_rate%mass,MASK=(u%halo(:)%compute)),&
                mass_rate_code_unit_2_MsunPerYr*sum(u%halo(:)%baryon_halo%cooling_rate%mass,MASK=(u%halo(:)%compute)),&
                mass_rate_code_unit_2_MsunPerYr*sum(u%halo(:)%galaxy%disc%sfr%mass,MASK=(u%halo(:)%compute)),&
                mass_code_unit_in_M_Sun*sum(u%halo(:)%galaxy%disc%stars%mass,MASK=(u%halo(:)%compute)),&
                mass_code_unit_in_M_Sun*sum(u%halo(:)%galaxy%bulge%stars%mass,MASK=(u%halo(:)%compute))/)
        !
        ! send data 
        call MPI_SEND(buffer,n_z_data-1,MPI_REAL8,main_process_rank,ztable_tag+2,MPI_COMM_WORLD,ierror)
    end if
    !
    if(allocated(buffer)) deallocate(buffer)
    
    return
  end subroutine tree_send_z_data
  
  !*****************************************************************************************************************

#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------  
  subroutine tree_send_data(ts,ih)                                     
  
    ! SEND SPECIFIC INFORMATIONS ABOUT A HALO IN A TREE (HID, host_HID, level ...)
  
    implicit none
    
    integer(kind=4),intent(in),optional     :: ts       ! timestep index 
    integer(kind=4),intent(in),optional     :: ih       ! halos index 
    
    logical(kind=4)                         :: go_down
    
    ! data are sent by physical process and receive by luminous process 
    
    go_down = .false. ! init
    
    if (present(ts) .and. present(ih)) go_down = .true.
    !
    ! send exit loop order
    call MPI_SEND(go_down,1,MPI_LOGICAL,rank+1,tree_tag+1,MPI_COMM_WORLD,ierror)
    !
    if (go_down) then
        !
        ! send the timestep index (use to print in the good tmp output file)
        call MPI_SEND(ts,1,MPI_INTEGER4,rank+1,tree_tag+2,MPI_COMM_WORLD,ierror)
        !
        ! send the number of computable halo at the current timestep
        call MPI_SEND(univ(ts)%nb_of_computable_halos,1,MPI_INTEGER4,rank+1,tree_tag+3,MPI_COMM_WORLD,ierror)
        !
        ! send the expansion factor
        call MPI_SEND(univ(ts)%aexp,1,MPI_REAL4,rank+1,tree_tag+4,MPI_COMM_WORLD,ierror)
        !
        ! send the halo identification parameter
        call MPI_SEND(univ(ts)%tree(ih)%HID,1,MPI_INTEGER8,rank+1,tree_tag+5,MPI_COMM_WORLD,ierror)
        !
        ! send the host halo identification parameter
        call MPI_SEND(univ(ts)%tree(ih)%host_HID,1,MPI_INTEGER8,rank+1,tree_tag+6,MPI_COMM_WORLD,ierror)
        !
        ! send the halo level
        call MPI_SEND(univ(ts)%tree(ih)%level,1,MPI_INTEGER4,rank+1,tree_tag+7,MPI_COMM_WORLD,ierror)
        !
        ! send the number of computed dads
        call MPI_SEND(univ(ts)%tree(ih)%n_computed_dads,1,MPI_INTEGER4,rank+1,tree_tag+8,MPI_COMM_WORLD,ierror)
        !
    end if

    return
  end subroutine tree_send_data 
  
  !*****************************************************************************************************************
  
  subroutine tree_receive_data(u,go_down)                                     
  
    ! RECEIVE SPECIFIC INFORMATIONS ABOUT A HALO IN A TREE (HID, host_HID, level ...)
  
    implicit none
    
    logical(kind=4),intent(out)         :: go_down
    
    type(timestep_type), intent(inout)  :: u        ! the local copy of the universe structure
        
    ! receive the exit loop order
    call MPI_RECV(go_down,1,MPI_LOGICAL,rank-1,tree_tag+1,MPI_COMM_WORLD,statut,ierror)
    !
    if (go_down) then 
        !
        ! receive the timestep index
        call MPI_RECV(u%ts,1,MPI_INTEGER4,rank-1,tree_tag+2,MPI_COMM_WORLD,statut,ierror)
        !
        ! receive the number of computable halo at the current timestep
        call MPI_RECV(u%nb_of_computable_halos,1,MPI_INTEGER4,rank-1,tree_tag+3,MPI_COMM_WORLD,statut,ierror)
        !
        ! receive the expansion factor
        call MPI_RECV(u%aexp,1,MPI_REAL4,rank-1,tree_tag+4,MPI_COMM_WORLD,statut,ierror)
        !
        ! receive the halo identification parameter
        call MPI_RECV(u%tree(1)%HID,1,MPI_INTEGER8,rank-1,tree_tag+5,MPI_COMM_WORLD,statut,ierror)
        !
        ! receive the host halo identification parameter
        call MPI_RECV(u%tree(1)%host_HID,1,MPI_INTEGER8,rank-1,tree_tag+6,MPI_COMM_WORLD,statut,ierror)
        !
        ! receive the halo level
        call MPI_RECV(u%tree(1)%level,1,MPI_INTEGER4,rank-1,tree_tag+7,MPI_COMM_WORLD,statut,ierror)
        !
        ! receive the number of computed dads
        call MPI_RECV(u%tree(1)%n_computed_dads,1,MPI_INTEGER4,rank-1,tree_tag+8,MPI_COMM_WORLD,statut,ierror)
        !
    end if
            
    return
  end subroutine tree_receive_data 
! ------------------------------------------------- 
#endif
! LUMINOUS_PROCESSES
  
  !*****************************************************************************************************************
  !
  ! SUBROUTINES
  !
  !*****************************************************************************************************************
  
  subroutine tree_set_reference_mass
    
    ! INITIALIZE REFERENCE MASS
    
    implicit none
    
    
    call IO_print_message('Resolution scale / Reference mass')
    
    call halo_set_reference_mass
    
    return
  end subroutine tree_set_reference_mass
  
  !*****************************************************************************************************************    
    
  subroutine tree_read(i)

    ! READ DARK-MATTER PROPERTIES COMPUTED IN THE DARK-MATTER SIMULATION
    ! - time step properties
    ! - tree structures
    ! - halo properties

    implicit none 

    ! name of step-file, tree-file, prop_files (input_path+filename)

    integer(kind=4),intent(in) :: i
    
    character(MAXPATHSIZE)     :: stepsfile
    character(MAXPATHSIZE)     :: treefile
    character(MAXPATHSIZE)     :: propsfile 
    character(MAXPATHSIZE)     :: message    

    write(message,'(a,a)') 'Load tree file set : ', trim(list_of_computed_files(i))
    call IO_print_message(message,only_rank=rank,component='tree')
    !
    ! input_path is defined in IO module
    ! time_step_filename, tree_struct_filename and tree_props_filename are given in input and saved in the IO module
    ! list_of_computed_files is defined in the IO module
    ! stepfile : contains nb of haloes, cosmic expansion factor and age of the universe at each ts
    write(stepsfile,'(a,a,a,a)') trim(tree_input_path), '/tstep/', trim(timestep_filename), trim(list_of_computed_files(i))
    ! tree-file : contains the information on the tree structure
    write(treefile,'(a,a,a,a)')  trim(tree_input_path), '/tree/', trim(tree_struct_filename), trim(list_of_computed_files(i))
    ! propsfile : contains the dm halo properties measured from the N-body simulation
    write(propsfile,'(a,a,a,a)') trim(tree_input_path), '/props/', trim(tree_props_filename), trim(list_of_computed_files(i))

    ! first read the properties of the trees 
    call tree_read_timestep_properties(stepsfile)
    ! then read the tree structure
    call tree_read_tree_structure(treefile,IO_return_ifile_from_file_index(list_of_computed_files(i)))
    ! finally read properties of each dark-matter halo in the tree
    call tree_read_halo_properties(propsfile)

    return

  contains

    ! *************************************************

    subroutine tree_read_timestep_properties(timestep_property_file)

      ! READ TIMESTEP PROPERTIES OF TREES BUILT FROM THE DARK-MATTER SIMULATION

      implicit none

      integer(kind=4)                   :: i                        ! loop index (to read data in a input file)
      integer(kind=4)                   :: stat                     ! error status for allocate process

      character(MAXPATHSIZE),intent(in) :: timestep_property_file   ! the global file path
      
      real(kind=4)                      :: dage

#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('tree_read_timestep_properties',only_rank=rank,component='tree')
! -------------------------------------------------
#endif
      
      ! Open the file with the properties of the trees at each timestep
      open(unit = tstep_file_unit, file = timestep_property_file, status = "old", form = "unformatted")
        
      ! Skip the number of timestep, already done
      read(tstep_file_unit)

      ! Dimension univ and read into univ the properties of the trees at each timestep
      if (.not. allocated(univ)) then
        !
        allocate(univ(0:nsteps),stat=stat)  ! allocate memory
        if (stat .ne. 0) then    
          !
          call IO_print_error_message('Cannot allocate univ array',only_rank = rank,called_by = 'tree_read_timestep_properties')           
          stop
        end if
      else
        !
        call IO_print_error_message('Cannot allocate univ array',only_rank = rank,called_by = 'tree_read_timestep_properties')           
        stop
      end if
      !
      ! Read informations associated to time-steps
      univ(:)%nb_of_halos = 0
      read(tstep_file_unit) (univ(i)%nb_of_halos,i=1,nsteps)   ! The number of haloes 
      read(tstep_file_unit) (univ(i)%aexp,i=1,nsteps)          ! The cosmic expansion factor 
      read(tstep_file_unit) (univ(i)%age_univ,i=1,nsteps)      ! The age of the universe in Gyr
      
      ! compute dage
      dage = univ(2)%age_univ - univ(1)%age_univ
      ! afect age at ts0
      univ(0)%age_univ = univ(1)%age_univ - dage

      ! Close the file
      close(tstep_file_unit)

#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('tree_read_timestep_properties ... done',only_rank=rank,component='tree')
! -------------------------------------------------
#endif

      return
    end subroutine tree_read_timestep_properties

    ! *************************************************

    subroutine tree_read_tree_structure(tree_structure_file,ifile)

      ! READ TREE STRUCTURE COMPUTE IN THE DARK-MATTER SIMULATION

      implicit none

      integer(kind=4)                   :: ifile                ! index of the intput tree properties file (allow to build the HID identification number)
      integer(kind=4)                   :: nIndexesPerHalo      ! must be .ge. 9
      integer(kind=4)                   :: nIDsPerHalo          ! nb of fields in treeID type 
      integer(kind=4),allocatable       :: indexes(:,:)         ! temporary array in which the information on the tree structure is read 
                                                                ! before it is stored into univ(:)%tree(:)
      integer(kind=4),allocatable       :: IDs(:,:)             ! tremporary array in which the information on the tree structure is read
                                                                ! before it is stored into univ(:)%treeID(:)
      integer(kind=4)                   :: ts                   ! timestep index loop
      integer(kind=4)                   :: ih                   ! halos index loop
      integer(kind=4)                   :: i,j                  ! loop indexes
      integer(kind=4)                   :: host                 ! host halo index
      integer(kind=4)                   :: level                ! structure level
      integer(kind=4)                   :: stat                 ! error status for allocate process

      character(MAXPATHSIZE),intent(in) :: tree_structure_file  ! the global file path

#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('tree_read_tree_structure',only_rank=rank,component='tree')
! -------------------------------------------------
#endif

      ! Open the file with the indexes that link a halo to its descendents, progenitors, hosts and subhaloes
      open(unit = tree_file_unit, file = tree_structure_file, status = 'old', form = 'unformatted')

      read(tree_file_unit) nts, nIDsPerHalo, nIndexesPerHalo 
      read(tree_file_unit) ! skip nb_of_halos(1:nsteps) because already read from stepsfile
      
      do ts = 0,nsteps
        ! 
        if (univ(ts)%nb_of_halos .gt. 0) then ! for all timesteps where there are haloes
            !
            ! Allocate univ(:)%tree(:)
            if (.not. allocated(univ(ts)%tree)) then
                !
                allocate(univ(ts)%tree(univ(ts)%nb_of_halos), stat = stat)     ! allocate memory
                if (stat .ne. 0) then     
                  !
                  call IO_print_error_message('Cannot allocate univ(ts)%tree array', only_rank = rank, called_by = 'tree_read_tree_structure')           
                  stop          
                end if
            else
                !
                call IO_print_error_message('Cannot allocate univ(ts)%tree array. This structure is already allocated', &
                only_rank = rank, called_by = 'tree_read_tree_structure')           
                stop 
            end if  
            !
            ! Dimension the temporary array into which we read the input data before we load them into univ
            if (.not. allocated(IDs)) then
                !
                allocate(IDs(nIDsPerHalo,univ(ts)%nb_of_halos), stat = stat)       ! allocate memory
                if (stat .ne. 0) then               
                    !
                    call IO_print_error_message('Cannot allocate IDs array',only_rank = rank,called_by = 'tree_read_tree_structure')  
                    call IO_print_message('with',only_rank=rank,component='tree',param_name=(/'stat                     '/),int_param_val=(/stat/))     
                    stop
                end if
            else
                !
                call IO_print_error_message('Cannot allocate IDs array. This structure is already allocated', & 
                only_rank = rank, called_by = 'tree_read_tree_structure')           
                stop
            end if
            !
            if (.not. allocated(indexes)) then
                !
                allocate(indexes(nIndexesPerHalo,univ(ts)%nb_of_halos), stat = stat)       ! allocate memory
                if (stat .ne. 0) then               
                    !
                    call IO_print_error_message('Cannot allocate indexes array',only_rank = rank,called_by = 'tree_read_tree_structure')  
                    call IO_print_message('with',only_rank=rank,component='tree',param_name=(/'stat                     '/),int_param_val=(/stat/))     
                    stop
                end if
            else
                !
                call IO_print_error_message('Cannot allocate indexes array. This structure is already allocated', &
                only_rank = rank, called_by = 'tree_read_tree_structure')           
                stop
            end if
            !
#ifdef HORIZON
! -------------------------------------------------
            ! 
            ! Read the tree structure into the temporary array      
            ! Read the ID of halos into the temporary array
            read(tree_file_unit) ((IDs(i,j),i=1,nIDsPerHalo),j=1,univ(ts)%nb_of_halos) 
            !
            ! *********************************!        
            ! IDs(1,ih)  = BushID              !
            ! IDs(2,ih)  = TreeID              !
            ! IDs(3,ih)  = HaloID              !
            ! IDs(4,ih)  = haloNum             !
            ! IDs(5,ih)  = haloTimestep        !                                      
            ! IDs(6,ih)  = FirstProgenitorID   !
            ! IDs(7,ih)  = NextProgenitorID    !
            ! IDs(8,ih)  = DescendentID        !
            ! IDs(9,ih)  = LastProgenitorID    !                                             
            ! IDs(10,ih) = HostHaloID          ! 
            ! IDs(11,ih) = HostSubID           !
            ! IDs(12,ih) = NextSubID           !
            ! *********************************! 
            !
            ! ****************************!
            ! IDs(2,ih) = my_number -> allows match with HaloMaker outputs (or halos_results.xxx files) 
            ! ****************************!      
            !
            ! Read the indexes into the temporary array
            read(tree_file_unit) ((indexes(i,j),i=1,nIndexesPerHalo),j=1,univ(ts)%nb_of_halos)  
            !
            ! ****************************!
            ! indexes(1,ih) = me          !
            ! indexes(2,ih) = descendent  !
            ! indexes(3,ih) = ndads       !        
            ! indexes(4,ih) = firstProg   !         
            ! indexes(5,ih) = nextProg    !       
            ! indexes(6,ih) = hosthalo    !        
            ! indexes(7,ih) = hostsub     !       
            ! indexes(8,ih) = nextsub     !        
            ! indexes(9,ih) = level       !        
            ! ****************************!
            !
            do ih = 1,univ(ts)%nb_of_halos
                !
                univ(ts)%tree(ih)%HID        = IO_generate_HID(ifile,ts,ih)
                univ(ts)%tree(ih)%HID_orig   = IDs(3,ih)       ! == my_number   -> allows match with HaloMaker outputs (or halos_results.xxx files) 
                univ(ts)%tree(ih)%me         = indexes(1,ih)   ! my number (index at timestep ts)
                univ(ts)%tree(ih)%descendent = indexes(2,ih)   ! index of my descendent (index at timestep ts +1)
                univ(ts)%tree(ih)%ndads      = indexes(3,ih)   ! number of progenitors
                univ(ts)%tree(ih)%firstProg  = indexes(4,ih)   ! index of my first progenitor (firstProg = -1 when there are no progenitors) (index at timestep ts -1)
                univ(ts)%tree(ih)%nextProg   = indexes(5,ih)   ! index of the next progenitor of my descendent (my brother next on the line, -1 if not) 
                univ(ts)%tree(ih)%hosthalo   = indexes(6,ih)   ! index of main halo
                univ(ts)%tree(ih)%hostsub    = indexes(7,ih)   ! index of host subhalo
                univ(ts)%tree(ih)%nextsub    = indexes(8,ih)   ! index of next substructure in main halo
                univ(ts)%tree(ih)%level      = indexes(9,ih)   ! level of substructure (level = 1 for a main halo)
            end do
            
! -------------------------------------------------
#endif
! HORIZON

#ifdef BOLSHOI
! -------------------------------------------------         
            
            do ih = 1,univ(ts)%nb_of_halos
                ! 
                ! Read the tree structure into the temporary array
                read(tree_file_unit) (IDs(i,ih),i=1,nIDsPerHalo), (indexes(j,ih),j=1,nIndexesPerHalo)
                
                ! ****************************!
                ! IDs(1,ih) = my_number -> allows match with Rockstar outputs
                ! ****************************!  
            
                ! ****************************!
                ! indexes(1,ih) = me          !
                ! indexes(2,ih) = descendent  !
                ! indexes(3,ih) = ndads       !        
                ! indexes(4,ih) = firstProg   !         
                ! indexes(5,ih) = nextProg    !       
                ! indexes(6,ih) = hosthalo    !        
                ! indexes(7,ih) = hostsub     !       
                ! indexes(8,ih) = nextsub     !        
                ! indexes(9,ih) = level       !        
                ! ****************************!
                
                univ(ts)%tree(ih)%HID        = IO_generate_HID(ifile,ts,ih)
                univ(ts)%tree(ih)%HID_orig   = IDs(1,ih)       ! == my_number   -> allows match with RockStar outputs
                univ(ts)%tree(ih)%me         = indexes(1,ih)   ! my number (index at timestep ts)
                univ(ts)%tree(ih)%descendent = indexes(2,ih)   ! index of my descendent (index at timestep ts +1)
                univ(ts)%tree(ih)%ndads      = indexes(3,ih)   ! number of progenitors
                univ(ts)%tree(ih)%firstProg  = indexes(4,ih)   ! index of my first progenitor (firstProg = -1 when there are no progenitors) (index at timestep ts -1)
                univ(ts)%tree(ih)%nextProg   = indexes(5,ih)   ! index of the next progenitor of my descendent (my brother next on the line, -1 if not) 
                univ(ts)%tree(ih)%hosthalo   = indexes(6,ih)   ! index of main halo
                univ(ts)%tree(ih)%hostsub    = indexes(7,ih)   ! index of host subhalo
                univ(ts)%tree(ih)%nextsub    = indexes(8,ih)   ! index of next substructure in main halo
                univ(ts)%tree(ih)%level      = indexes(9,ih)   ! level of substructure (level = 1 for a main halo)
            end do
            
! -------------------------------------------------
#endif
! BOLSHOI
            !
#ifndef UNLINKED_TREES  
! -------------------------------------------------     
            !
            ! set host_HID
            ! If the halo is identified as a sub-structure (satellite) of an other halo
            ! host_HID is the halo identification number of the host (main) halo
            ! If the halo is a main structure host_HID = HID
            do ih = 1,univ(ts)%nb_of_halos
                !
                level = univ(ts)%tree(ih)%level
                if (level .gt. 1) then
                    !
                    ! The halo is a substructure
                    host = univ(ts)%tree(ih)%hosthalo
                    if (host .gt. 0) then
                        !
                        if (univ(ts)%tree(host)%level .eq. 1) then
                            !
                            univ(ts)%tree(ih)%host_HID = univ(ts)%tree(host)%HID
                        else
                            !
                            call IO_print_error_message('Host halo with level != 1',only_rank=rank,called_by='tree_read_tree_structure')
                            call IO_print_message('used',only_rank=rank,component='tree', &
                                    param_name=(/'tree(ih)%level           ','tree(ih)%hosthalo        ','tree(host)%level         '/), &
                                    int_param_val=(/level,host,univ(ts)%tree(host)%level/))  
                            stop ! stop the program
                        end if
                    else
                        !
                        call IO_print_error_message('Substructure without host halo',only_rank=rank,called_by='tree_read_tree_structure')
                        call IO_print_message('used',only_rank=rank,component='tree', &
                                    param_name=(/'tree(ih)%level           ','tree(ih)%hosthalo        '/), &
                                    int_param_val=(/level,host/))  
                        stop ! stop the program
                    end if
                else
                    !
                    ! the halo is a main-structure
                    univ(ts)%tree(ih)%host_HID = univ(ts)%tree(ih)%HID
                end if
            end do
! -------------------------------------------------   
#endif         
! UNLINKED_TREE            
            !
        end if
        !
        ! Deallocate the temporary array
        if (allocated(indexes)) deallocate(indexes) 
        if (allocated(IDs))     deallocate(IDs)
      end do

      ! Close the file with indexes    
      close(tree_file_unit)

#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('tree_read_tree_structure ... done ',only_rank=rank,component='tree')
! -------------------------------------------------
#endif
      
      return   
    end subroutine tree_read_tree_structure
    
    ! *************************************************

    subroutine tree_read_halo_properties(halo_property_file)

      ! READ HALO PROPERTIES COMPUTED IN THE DARK-MATTER SIMULATION

      implicit none

      integer(kind=4)                   :: nPropsPerHalo
      integer(kind=4)                   :: ts                 ! loop index (timesteps)
      integer(kind=4)                   :: ih                 ! loop index (halos)
      integer(kind=4)                   :: i                  ! loop indexes (for read data in the input file)
      integer(kind=4)                   :: stat               ! error status for allocate process
      
      character(MAXPATHSIZE),intent(in) :: halo_property_file ! the global file path

      real(kind=4),allocatable          :: props(:,:)         ! temporary array in which the properties of dm haloes are read 
                                                              ! before they are stored into univ(:)%halo(:)
#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('tree_read_halo_properties',only_rank=rank,component='tree')
! -------------------------------------------------
#endif

      ! Open the file with the dm halo properties 
      open(unit = props_file_unit, file = halo_property_file, status = 'old', form = 'unformatted')

      read(props_file_unit) nts,nPropsPerHalo
      read(props_file_unit) ! skip nb_of_halos(1:nsteps) because already read from stepsfile
       
      do ts = 0,nsteps
        !
        univ(ts)%ts = ts ! save the timestep index
        !
        if (univ(ts)%nb_of_halos > 0) then
            !
            ! Dimension univ(ts)%halo(:) (same dimensions as univ(:)%tree(:))
            if (.not. allocated(univ(ts)%halo)) then
                !
                allocate(univ(ts)%halo(univ(ts)%nb_of_halos), stat = stat)
                if (stat .ne. 0) then               
                  !
                  call IO_print_error_message('Cannot allocate univ(ts)%halo array', only_rank = rank, called_by = 'tree_read_halo_properties') 
                  stop
                end if
            else
                !
                call IO_print_error_message('Cannot allocate univ(ts)%halo array. This structure is already allocated', &
                only_rank = rank, called_by = 'tree_read_halo_properties') 
                stop
            end if
            !
            ! Dimension the temporary array props appropriately before reading dm properties into it 
            if (.not. allocated(props)) then 
                !
                allocate(props(nPropsPerHalo,univ(ts)%nb_of_halos), stat = stat)
                if (stat .ne. 0) then       
                    !         
                    call IO_print_error_message('Cannot allocate props array', only_rank = rank, called_by = 'tree_read_halo_properties') 
                    stop
                end if
            else
                ! 
                call IO_print_error_message('Cannot allocate props array. This structure is already allocated', &
                only_rank = rank, called_by = 'tree_read_halo_properties') 
                stop
            end if
            !
#ifdef HORIZON
! -------------------------------------------------
            !
            ! Read the DM properties into the temporary array
            read(props_file_unit) ((props(i,ih),i=1,nPropsPerHalo),ih=1,univ(ts)%nb_of_halos) 
            !
            ! ****************************!
            ! props(1:3,ih)   = x,y,z     ! 
            ! props(4:6,ih)   = Vx,Vy,Vz  !
            ! props(7,ih)     = M_tot     !
            ! props(8,ih)     = r         ! 
            ! props(9,ih)     = Spin      !
            ! props(10,ih)    = R_vir     !
            ! props(11,ih)    = M_vir     !
            ! props(12,ih)    = T_vir     !
            ! props(13,ih)    = cvel      !
            ! props(14,ih)    = dmacc     !
            ! props(15,ih)    = frag      !
            ! props(16:18,ih) = Lx,Ly,Lz  !
            ! props(19,ih)    = Ep        !
            ! props(20,ih)    = Ek        !
            ! props(21,ih)    = Et        ! 
            ! ****************************!
            !
            do ih = 1,univ(ts)%nb_of_halos
              !
              ! Initialize the halo strcuture univ(ts)%halo
              call halo_void(univ(ts)%halo(ih))
              ! Transfer them into univ(:)%halo(:) by calling subroutine halo_set_props(h,ts,props(:,ih))
              call halo_set_properties(univ(ts)%halo(ih),ts,props(:,ih))
            end do
        
! -------------------------------------------------
#endif
! HORIZON

#ifdef BOLSHOI
! -------------------------------------------------

            do ih = 1,univ(ts)%nb_of_halos
                ! 
                ! Read the DM properties into the temporary array
                read(props_file_unit) (props(i,ih),i=1,nPropsPerHalo)
                !
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
                !
                ! Initialize the halo structure univ(ts)%halo
                call halo_void(univ(ts)%halo(ih))
                ! Transfer them into univ(:)%halo(:) by calling subroutine halo_set_props(h,ts,props(:,ih))
                call halo_set_properties(univ(ts)%halo(ih),ts,props(:,ih))
            end do

! -------------------------------------------------
#endif
! BOLSHOI
            
            ! Deallocate the temporary array
            if (allocated(props)) deallocate(props)
        end if ! nb_halo > 0
      end do
      
      ! Close the file
      close(props_file_unit)

#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('tree_read_halo_properties ... done ',only_rank=rank,component='tree')
! -------------------------------------------------
#endif
    
      return
    end subroutine tree_read_halo_properties
   
    ! *************************************************
    ! end contains for tree read

  end subroutine tree_read

  !*****************************************************************************************************************
  
  subroutine tree_compute_halo_derived_properties

    ! COMPUTE DARK-MATTER (dm) PROPERTIES WHICH CANNOT BE COMPUTED DIRECTLY IN THE DARK-MATTER SIMULATION
    ! integreted mass Macc
    ! dark matter accretion rate (left and right)
    ! take account dmacc part and (bad-branch) accretion

    implicit none

    integer(kind=4)    :: ts              ! loop index (timesteps)
    integer(kind=4)    :: ih              ! loop index (halos)
    integer(kind=4)    :: ih2
    integer(kind=4)    :: ip,idesc        ! progenitor and descendent indexes 
    
    real(kind=8)       :: Macc            ! integreted accreted mass
    real(kind=8)       :: z               ! redshift
    real(kind=8)       :: tot, me         ! tmp variable, the total mass (over all progenitors), the halo mass
    real(kind=8)       :: dt_left         ! time-step 
    real(kind=8)       :: acc_rate_right  ! dark matter accretion rate (during the last half time-step)   
    real(kind=8)       :: acc_rate_left   ! dark matter accretion rate (during the next half time-step) 
    real(kind=8)       :: merger_mass     ! Mass acquires by merger 
    real(kind=8)       :: bad_branch_mass ! integreted mass contains in dead branches connected to a halo
    real(kind=8)       :: t1, t2          ! age of the universe at two consecutive time-steps
    real(kind=8)       :: age_form        ! age of formation of a given dm-structure (age of the universe) 
    real(kind=8)       :: age_merge       ! age of the univers when the structure(halo+galaxy) have been generated from its progenitors 
    real(kind=8)       :: r               ! a random number

#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('tree_compute_halo_derived_properties',only_rank=rank,component='tree')
! -------------------------------------------------
#endif
    !
    ! ***********************************
    ! TREE LINKS
    ! ***********************************
    ! @ ts = nsteps: set desc = -1
    if (univ(nsteps)%nb_of_halos .gt. 0) then       ! if the timestep contains some halos
        !
        do ih = 1,univ(nsteps)%nb_of_halos
            !
            univ(nsteps)%tree(ih)%descendent  = -1  ! stop the run throught the descendents
        end do
    end if
    !
    ! ***********************************
    ! ACCRETED MASS Macc
    ! ***********************************
    !
    ! Compute the DM mass accreted by the halo by integrating over all its progenitors
    ! The loop is computed from the first halos (ts = 1) to the last identified halo (ts = nsteps)
    ! Next flag bad halos (Positive total energy, negatives spin, or null integrated mass)
    !
    do ts = 1,nsteps
      !
      if (univ(ts)%nb_of_halos .gt. 0) then                ! if the timestep contains some halos
        !
        do ih = 1,univ(ts)%nb_of_halos
          !
          ! Sum the accreted mass over all progenitors, dmacc is the new smooth accreted mass  
          Macc = univ(ts)%halo(ih)%dm%dmacc                ! start from the halo itself  
          ip   = univ(ts)%tree(ih)%firstProg               ! move to the main progenitor
          do while (ip .gt. 0)                             ! continue until there is no progenitor left
            !
            Macc = Macc + univ(ts-1)%halo(ip)%dm%M_acc     ! add progenitors
            ip = univ(ts-1)%tree(ip)%nextProg              ! jump to the next progenitor
          end do
          call dm_set_M_acc(univ(ts)%halo(ih)%dm,Macc)      ! set M_acc
          ! Flag bad dm properties (Positive total energy, negatives spin, or null integrated mass)
          call dm_flag_bad_properties(univ(ts)%halo(ih)%dm)
        end do ! end ih loop
      end if
    end do ! end ts loop
    !
    ! ***********************************
    ! CLEAN PROCESS
    ! ***********************************
    !
    ! Compute halo%compute flag taking into account dm-flag and fixed goodness of branch
    call tree_clean 
    !
    ! ***********************************
    ! HALO SELECTION (FOLLOW MODE)
    ! ***********************************
    !
    ! select followed halos
    ! first run, select halos
    if (nb_follow_up_halos .gt. 0.d0) then
       !
       do ts = 1,nsteps 
            !
            if (univ(ts)%nb_of_halos .gt. 0) then
                !
                do ih = 1,univ(ts)%nb_of_halos
                    !
                    do ih2 = 1, nb_follow_up_halos 
                        !
                        if (univ(ts)%tree(ih)%HID .eq. list_of_follow_up_halo(ih2)) then
                            !
                            call tree_followed_branch(ih,ts)
                            exit
                        end if
                    end do
                end do
            end if
        end do
    end if
#ifdef FOLLOW 
! -------------------------------------------------   
    ! in the follow_up mode only a selection of halos are evolved
    ! these halos are the progenitors of a given list of haloes   
    ! second run, kill not followed halo
    do ts = 1,nsteps 
        !
        if (univ(ts)%nb_of_halos .gt. 0) then
            !
            do ih = 1,univ(ts)%nb_of_halos
                !
                if (.not. univ(ts)%halo(ih)%followed) univ(ts)%halo(ih)%compute = .false.
            end do
        end if
    end do       
! -------------------------------------------------    
#endif
! FOLLOW  
    !
    ! ***********************************
    ! DM-ACCRETION RATES
    ! ***********************************
    !
    ! ACCRETION RATE LEFT
    ! Compute the accretion rate "left" of the dark matter component of the halo
    ! Take account the bad branch part. A bad branch is considered as an exceptional accretion event onto the dark matter halo structure
    ! This kind of branch is computed 
    ! - like a smooth accretion on to the halo
    ! - don't have baryonic gas evolution before the accretion onto the halo
    !   
    do ts = 1,nsteps
      !   
      dt_left = univ(ts)%age_univ - univ(ts-1)%age_univ
      if (dt_left .le. 0.d0) then
        !
        call IO_print_error_message('Bad time-step value: dt_left < 0', &
          only_rank = rank, called_by = 'tree_compute_halo_derived_properties') 
        call IO_print_message('',only_rank=rank,component='tree',param_name=(/'dt_left                  '/),real_param_val=(/dt_left/))
        stop
      end if
      !
      univ(ts)%nb_of_computable_halos  = 0 ! initialize the number of computable halos at a given timestep
      univ(ts)%nb_of_followed_up_halos = 0 ! initialize the number of follow_up halos at a given timestep
      ! if the current time-step contains some halos
      if (univ(ts)%nb_of_halos .gt. 0) then
        !
        do ih = 1,univ(ts)%nb_of_halos
          ! 
          if (univ(ts)%halo(ih)%followed) univ(ts)%nb_of_followed_up_halos = univ(ts)%nb_of_followed_up_halos +1
          acc_rate_left = 0.d0   ! init
          if (univ(ts)%halo(ih)%compute) then                                     ! only if the halo is on a good branch
            !
            univ(ts)%nb_of_computable_halos = univ(ts)%nb_of_computable_halos +1  ! compute number of computable haloes
            ip = univ(ts)%tree(ih)%firstProg                                      ! move to the main progenitor
            merger_mass     = 0.d0   ! initialize the merger mass (the dark-matter mass acquires by merger event)
            bad_branch_mass = 0.d0   ! initialize the bad_branch mass (the mass accreted by all dead branch connected to the halo)                      
            do while (ip .gt. 0)     ! continue until there is no progenitor left
              !
              if (univ(ts-1)%halo(ip)%compute) then
                !
                ! compute merger mass
                if (ip .eq. univ(ts)%tree(ih)%firstProg) then                     ! we are on the main branch
                  !
                  merger_mass = univ(ts-1)%halo(ip)%dm%merger_mass                ! copy the previous value at the previous time-step save in the main progenitor
                else
                  !
                  merger_mass = merger_mass + univ(ts-1)%halo(ip)%dm%M_acc        ! add mass of the halo the merger mass
                end if
              else
                !
                ! this progenitor halo is not computable
                ! add its mass to the bad_branch_mass
                bad_branch_mass = bad_branch_mass + univ(ts-1)%halo(ip)%dm%M_acc 
              end if  
              ip = univ(ts-1)%tree(ip)%nextProg                                  ! jump to the next progenitor
            end do
            !
            ! the accretion is the sum of real smooth accretion (dmacc) and the exceptional accretion due to bad-branch mass
            acc_rate_left = (univ(ts)%halo(ih)%dm%dmacc + bad_branch_mass) / dt_left  ! dt_left is the time between time-step ts and ts-1
            !
            if (acc_rate_left .lt. 0.d0) then
              !
              call IO_print_error_message('Bad acc_rate_left value: acc_rate_left < 0', &
                only_rank = rank, called_by = 'tree_compute_halo_derived_properties') 
              write(errunit,*) ' --> acc_rate_left :', acc_rate_left 
              write(errunit,*) ' --> With : dmacc  :', univ(ts)%halo(ih)%dm%dmacc, &
                ' bad_branch_mass : ', bad_branch_mass, ' dt_left : ', dt_left
              stop  ! stop the program
            end if
            !
            ! save some dark-matter properties in the dm-strcture
            call dm_set_acc_rate_left(univ(ts)%halo(ih)%dm,acc_rate_left)     ! set accretion rate left
            call dm_set_merger_mass(univ(ts)%halo(ih)%dm,merger_mass)         ! set merger mass
            call dm_set_merger_rate(univ(ts)%halo(ih)%dm,merger_mass/dt_left) ! set merger rate
            call dm_set_bad_branch_mass(univ(ts)%halo(ih)%dm,bad_branch_mass) ! set bad branch mass
            !
          end if ! end if compute
        end do ! end ih loop
      end if ! nb_of_halos > 0
    end do ! end ts loop
    !
    ! ACCRETION RATE RIGHT 
    ! Compute the accretion rate right of the dark matter component
    ! It's a redistribution of the accretion rate of the descendent halo onto its progenitor
    ! This redistribution is done by taking into account 
    !     - the accretion rate through the progenitors 
    !     - if all progenitors haven't accretion the redistribution is done by taking into account the integrated mass Macc
    do ts = 1, nsteps-1 ! there is no halos at the time-step 0, so we start the loop at ts = 1
      !
      if (univ(ts)%nb_of_halos .gt. 0) then
        !
        do ih = 1,univ(ts)%nb_of_halos
          acc_rate_right = 0.d0                                   ! init
          if (univ(ts)%halo(ih)%compute) then                     ! only if the halo is on a good branch
            !
            idesc = univ(ts)%tree(ih)%descendent                  ! index of my descendent
            if (idesc .gt. 0) then
              !
              if (univ(ts+1)%halo(idesc)%dm%acc_rate_left .le. 0.d0) then
                !
                ! the descendent has not accretion, don't report on progenitors
                acc_rate_right = 0.d0                       
              else
                !
                ! the descendent halo has got posivive accretion rate (acc_rate_left > 0.)
                if (univ(ts+1)%tree(idesc)%n_computed_dads .eq. 1) then
                  !
                  ! the descendent has just one progenitor and i'am this progenitor because univ(ts)%halo(ih)%compute = .true.
                  acc_rate_right = univ(ts+1)%halo(idesc)%dm%acc_rate_left
                else
                  !
                  ! the descendent halo has got more than one progenitors
                  ! separate accretion rate of the descendent onto all its progenitors 
                  ! take into account progenitor accretion rate weight
                  me  = univ(ts)%halo(ih)%dm%acc_rate_left           ! my accretion rate
                  tot = 0.d0                                         ! initialize to 0. the total accretion rate through progenitors 
                  ip  = univ(ts+1)%tree(idesc)%firstProg             ! move to the main progenitor of my descendent
                  do while (ip .gt. 0)                               ! continue until there is no progenitor left
                    !
                    if (univ(ts)%halo(ip)%compute) then
                      tot = tot + univ(ts)%halo(ip)%dm%acc_rate_left ! only if is a computable halo
                    end if
                    ip  = univ(ts)%tree(ip)%nextProg                 ! jump to the next progenitor
                  end do
                  if (tot .eq. 0.d0) then
                    !
                    ! take into account progenitor mass weight
                    me  = univ(ts)%halo(ih)%dm%M_tot          ! my halo finder mass
                    tot = 0.d0                                ! initialize to 0. the total accreted mass through progenitors
                    ip  = univ(ts+1)%tree(idesc)%firstProg    ! move to the main progenitor of my descendent
                    do while (ip .gt. 0)                      ! continue until there is no progenitor left
                      !
                      if (univ(ts)%halo(ip)%compute) then
                        ! 
                        tot = tot + univ(ts)%halo(ip)%dm%M_tot ! only if is a computable halo
                      end if
                      ip  = univ(ts)%tree(ip)%nextProg        ! jump to the next progenitor
                    end do
                    if (tot .eq. 0.d0) then
                      !
                      call IO_print_error_message('Bad accretion rate distribution: tot = 0', &
                        only_rank = rank, called_by = 'tree_compute_halo_derived_properties') 
                      stop ! stop the program
                    end if
                  end if
                  ! compute accretion rate left associated to the currrent ih halo
                  ! fraction of the descendent accretion rate
                  if (me/tot .gt. 1.d0) then
                    !
                    call IO_print_error_message('Bad accretion rate distribution: (me/tot) > 1', &
                      only_rank = rank, called_by = 'tree_compute_halo_derived_properties') 
                    stop ! stop the program
                  end if
                  acc_rate_right = (me/tot) * univ(ts+1)%halo(idesc)%dm%acc_rate_left
                end if ! n_computed_dads > 1
              end if ! acc_rate_left
              if (acc_rate_right .lt. 0.d0) then
                !
                call IO_print_error_message('Wrong value for acc_rate_right', &
                  only_rank = rank, called_by = 'tree_compute_halo_derived_properties') 
                stop ! stop the program
              end if
              !
              if (acc_rate_right .gt. univ(ts+1)%halo(idesc)%dm%acc_rate_left) then
                !
                call IO_print_error_message('Wrong value for acc_rate_right. acc_rate_right(prog) > acc_rate_left(desc)', &
                  only_rank = rank, called_by = 'tree_compute_halo_derived_properties') 
                write(errunit,*) ' --> With acc_rate_right : ', acc_rate_right, ', acc_rate_left : ', univ(ts+1)%halo(idesc)%dm%acc_rate_left
                stop ! stop the program
              end if
              !
              ! save acc_rate_right value in the dm structure
              call dm_set_acc_rate_right(univ(ts)%halo(ih)%dm,acc_rate_right)  ! set accretion rate right
              !
            end if ! idesc > 0
          end if ! end if compute
        end do ! end ih loop
      end if ! there are haloes
    end do ! end ts loop
    !
    ! ***********************************
    ! OTHER DM-PROPERTIES
    ! ***********************************
    !
    ! Compute properties of the halo that can't be directly measured in the N-body simulation
    ! - Concentration, virial Temperature, dynamical time .... 
    ! - A formation age (age of the universe) is also compute for each dm-structure
    do ts = 1,nsteps
      !
      z =  1./univ(ts)%aexp - 1.  ! compute redshift
      if (univ(ts)%nb_of_halos .gt. 0) then
        !
        do ih = 1,univ(ts)%nb_of_halos
          !
          if (univ(ts)%halo(ih)%compute) then
            !
            ! set formation age
            age_form = -1.d0 ! init
            if (univ(ts)%tree(ih)%n_computed_dads .gt. 0) then
                !
                ! run throught progenitors
                ip = univ(ts)%tree(ih)%firstProg     ! move to the main progenitor                   
                do while (ip .gt. 0)                 ! continue until there is no progenitor left
                   !
                   if (univ(ts-1)%halo(ip)%compute) then
                      !
                      ! the current progenitor is a computable halo
                       if (age_form .lt. 0.d0) then
                          !
                          age_form = univ(ts-1)%halo(ip)%dm%age_form ! init
                       else
                          !
                          ! set to minimum value
                          age_form = min(age_form,univ(ts-1)%halo(ip)%dm%age_form)
                       end if
                    end if   
                    ip = univ(ts-1)%tree(ip)%nextProg  ! jump to the next progenitor
                end do
            else
                !
                ! if a halo doesn't have any computable progenitor,
                ! this new halo has been formed between the last and the current time-step
                ! the formation age (age of the universe) is randomly set between these two main time-steps
                call random_number(r)
                t1        = univ(ts-1)%age_univ   ! age of the universe at ts-1 
                t2        = univ(ts)%age_univ     ! age of the universe at ts, t2 > t1
                ! the formation age is randomly set between the two time-steps
                age_form = t1 + num_precision+r*(t2-t1-2.d0*num_precision)  ! > t1 and < t2
                !
                if ((age_form .le. t1) .or. (age_form .ge. t2)) then
                   !
                   call IO_print_error_message('Wrong value for age_form', &
                        only_rank = rank, called_by = 'tree_compute_halo_derived_properties') 
                   call IO_print_message('',only_rank=rank,component='halo', &
                            param_name = (/'t1                       ','tf                       ','t2                       '/), &
                            real_param_val  = (/t1,age_form,t2/))   
                   stop ! stop the program
                end if
            end if
            call dm_set_age_form(univ(ts)%halo(ih)%dm,age_form)
            !
            ! set merging age 
            age_merge = -1.d0 ! init
            if (univ(ts)%tree(ih)%n_computed_dads .gt. 0) then  
                !
                t1 = univ(ts-1)%age_univ  ! age of the universe at ts-1 
                t2 = univ(ts)%age_univ    ! age of the universe a ts, t2 > t1
                !
                if (univ(ts)%tree(ih)%n_computed_dads .gt. 1) then
                    !
                    ! the merging age (age of the universe) is randomly set between these two main time-steps
                    call random_number(r)
                    age_merge = t1 + num_precision+r*(t2-t1-2.d0*num_precision)   
                else
                    !
                    ! In case of secular evolution, the "merging" process occurs at the half time-step
                    age_merge = 5.d-1*(t1+t2)
                end if
                !
                if ((age_merge .le. t1) .or. (age_merge .ge. t2))  then
                    !
                    call IO_print_error_message('Wrong value for age_merge', &
                        only_rank = rank, called_by = 'tree_compute_halo_derived_properties') 
                    call IO_print_message('',only_rank=rank,component='halo', &
                        param_name = (/'t1                       ','tm                       ','t2                       '/), &
                        real_param_val  = (/t1,age_merge,t2/))  
                    stop ! stop the program
                end if
            end if
            call dm_set_age_merge(univ(ts)%halo(ih)%dm,age_merge)
            !
            call dm_compute_derived_properties(univ(ts)%halo(ih)%dm,z)  
          end if
        end do ! end ih loop
      end if
    end do ! end ts loop
    
#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('tree_compute_halo_derived_properties .. done ',only_rank=rank,component='tree')
! -------------------------------------------------
#endif

    return 
  end subroutine tree_compute_halo_derived_properties
   
  ! *****************************************************************************************************************

  subroutine tree_clean

    ! CLEAN UP THE TREE ACCORDING TO ITS PROPERTIES (total energy, spin, integreted mass)

    implicit none

    integer(kind=4)                :: ts,ih        ! dumb indexes      
    
    real(kind=8)                   :: goodness     ! tmp variable

#ifdef PRINTALL
! -------------------------------------------------
    if (nbproc .eq. 1) call IO_print_message('tree_clean',only_rank=rank,component='tree')
! -------------------------------------------------
#endif

    ! 1/ scan all halos and define goodness of tree branches as nb of (dm-)good halos / nb of halos on branch
    ! -> save goodness into dm%goodness_of_branch (initialized to -1 elsewhere)
    ! NB: a branch connects main progs and descendent.
    ! the tree is run through ts = 1 to ts = nsteps, a branch is checked until its first connection with an other branch
    do ts = 1, nsteps 
      do ih = 1, univ(ts)%nb_of_halos
        if (univ(ts)%halo(ih)%dm%goodness_of_branch .lt. 0) then  ! haven't done this branch yet ...
          goodness = tree_compute_dm_goodness_of_branch(ih,ts)    ! get fraction of (dm-)good halos in branch 
          call tree_set_dm_goodness_of_branch(ih,ts,goodness)     ! apply this value for all halos in the branch
        end if
      end do
    end do
    !
    ! 2/ initialization of halo%compute flag -> set them to dm%compute
    do ts = 1, nsteps 
      do ih = 1, univ(ts)%nb_of_halos
        if (univ(ts)%tree(ih)%level .eq. 1) then
          univ(ts)%halo(ih)%compute = univ(ts)%halo(ih)%dm%compute
        else
          ! if it is a sub-structure, we consider (a priori) that is a good structure (in term of halo compute flag)
          univ(ts)%halo(ih)%compute = .true.
        end if
      end do
    end do
    !
    ! 3/ RESCUE STEP ! set main-halos (level = 1) compute flag to true if they become a sub-halo or a good halo shortly (i.e. within "nb_of_desc_halo_test" timesteps)
    ! NB: this is done following main branch (i.e. does not affect halos which have a sub as descendent but are not the main prog of that sub)
    do ts = 1, nsteps 
      do ih = 1, univ(ts)%nb_of_halos
        ! we apply the rescue step only on :
        !  - bad main-structure 
        !  - with an integrated accretd mass Macc > 0 
        !  - only on a branch with a dm%goodness_of_branch greater than the treshold "min_dm_goodness_of_branch"
        if ((.not. univ(ts)%halo(ih)%compute) .and. (univ(ts)%tree(ih)%level .eq. 1) &
          .and. (univ(ts)%halo(ih)%dm%M_acc .gt. 0.) &
          .and. (univ(ts)%halo(ih)%dm%goodness_of_branch .gt. min_dm_goodness_of_branch)) then
          ! Here the current halo have halo%compute = .false.
          univ(ts)%halo(ih)%compute = tree_halo_become_good(ih,ts) ! The function return .true. if the current halo respect at least one rescue condition
        end if
      end do
    end do
    !
    ! 4/ compute halo_goodness_of_branch as the number of main-halos (level 1) with halo%compute = true over the number of main-halos (level 1)
    ! this ratio is computed only through the main branch   
    do ts = 1, nsteps 
      do ih = 1, univ(ts)%nb_of_halos
        if (univ(ts)%halo(ih)%goodness_of_branch .lt. 0) then    ! haven't done this branch yet ...
          goodness = tree_compute_halo_goodness_of_branch(ih,ts) ! get fraction of (halo-)good halos in branch 
          call tree_set_halo_goodness_of_branch(ih,ts,goodness)  ! apply this value for all halos in the branch
        end if
      end do
    end do
    !
    ! 5/ re-define halo%compute as follows:
    ! set to .true. if the halo :
    ! - has an integrated accretd mass Macc > 0 
    ! - is in a branch with halo_goodness_of_branch >= "min_goodness_of_branch"
    ! all other halos (in branches with h_g_o_b < min_goodness_of_branch) are set to halo%compute = .false.
    do ts = 1, nsteps 
      do ih = 1, univ(ts)%nb_of_halos
        if ((univ(ts)%halo(ih)%dm%M_acc .gt. 0.) .and. &
          (univ(ts)%halo(ih)%goodness_of_branch .ge. min_halo_goodness_of_branch)) then  
          univ(ts)%halo(ih)%compute = .true.    ! we hope that is a good branch
        else
          univ(ts)%halo(ih)%compute = .false.   ! it is definitively not a good halo                                          
        end if
      end do
    end do
    !
    ! 6/ set halo%compute to false for:
    ! - computable sub-structures 
    ! - without progenitor
    do ts = 1, nsteps 
      do ih = 1, univ(ts)%nb_of_halos
        if (univ(ts)%tree(ih)%level .gt. 1) then            ! The halo is a sub halo
          if (univ(ts)%halo(ih)%compute) then               ! computable
            if (univ(ts)%tree(ih)%firstProg .lt. 0) then    ! but without any progenitor                                
              univ(ts)%halo(ih)%compute = .false.           ! it is definitively not a good halo
            else
              ! if the sub-halo has a progenitors, we propagate the halo%compute state of the first progenitor
              univ(ts)%halo(ih)%compute = univ(ts-1)%halo(univ(ts)%tree(ih)%firstProg)%compute
            end if
          end if
        endif
      end do
    end do
    !
    ! 7/ kill branches connected to a pathological halo (halo%compute = .fasle.)
    do ts = nsteps, 1, -1
      do ih = 1, univ(ts)%nb_of_halos
        if (.not. univ(ts)%halo(ih)%compute) then ! remove full branch (and those attached)
          call tree_kill_branch(ih,ts)
        end if
      end do
    end do
    !
    ! 8/ computed and set n_computed_dads
    do ts = 1,nsteps
       do ih = 1, univ(ts)%nb_of_halos
          call tree_computed_and_set_n_computed_dads(ih,ts)
       end do
    end do
    !
    ! 9/ interpolate properties of (sub-)halo that will be computed but which have bad dm properties (spin)
    do ts = 1,nsteps
       do ih = 1, univ(ts)%nb_of_halos
          if (univ(ts)%halo(ih)%compute .and. (.not. univ(ts)%halo(ih)%dm%compute)) then
             call tree_interpolate_bad_properties(ih,ts)
          end if
       end do
    end do

#ifdef PRINTALL
! -------------------------------------------------
    if (nbproc .eq. 1) call IO_print_message('tree_clean .. done ',only_rank=rank,component='tree')
! -------------------------------------------------
#endif

    return
      
  contains
        
    ! *************************************************

    function tree_compute_dm_goodness_of_branch(ih,ts)

      ! COMPUTE THE NUMBER OF (dm-)GOOD MAIN-HALO / NUMBER OF MAIN HALO IN THE CURRENT BRANCH

      implicit none

      integer(kind=4),intent(in)  :: ts                                 ! the current time-step
      integer(kind=4),intent(in)  :: ih                                 ! the current halo
      integer(kind=4)             :: me,idesc,ifirst                    ! some indexes for: a halo, a descendent, a first prog
      integer(kind=4)             :: ts2                                ! a local time-step  
      integer(kind=4)             :: n_halos_in_branch                  ! number of halos in the current branch
      integer(kind=4)             :: n_good_halos                       ! number_of good halos in the current branch
      
      real(kind=8)                :: tree_compute_dm_goodness_of_branch ! number of (dm-)good main-halo/number of main halo in the current branch
      
      n_good_halos      = 0  ! init 
      n_halos_in_branch = 0  ! init 

      ! check for the current halo (the halo ih at ts)
      ! a main halo is caracerised by level = 1
      if (univ(ts)%tree(ih)%level .eq. 1) then
        !
        n_halos_in_branch = 1
        ! if the current halo is a (dm-)computable halo
        ! this flag has been computed previously by: dm_flag_bad_properties in tree_compute_halo_derived_properties
        if (univ(ts)%halo(ih)%dm%compute) n_good_halos = 1
      end if
      !
      ! scan descendents
      idesc = univ(ts)%tree(ih)%descendent                 ! index of my (halo ih at ts) descendent at the next timestep
      ts2   = ts + 1                                       ! time-step of my descendent
      me    = ih                                           ! my index
      do while (idesc .gt. 0)
        !
        if (univ(ts2)%tree(idesc)%firstProg .eq. me) then  ! the descendent is still on main branch (i.e. main prog of desc)
          !
          if (univ(ts2)%tree(idesc)%level .eq. 1) then     ! the halo is a main halo
            n_halos_in_branch = n_halos_in_branch + 1      ! update the number of main-halo in the current branch (good and bad)  
            ! update the number of good identified main-halo
            if (univ(ts2)%halo(idesc)%dm%compute) n_good_halos = n_good_halos + 1                
          end if
            me    = idesc
            idesc = univ(ts2)%tree(me)%descendent          ! go to the next descendent (on the same branch)
            ts2   = ts2 + 1                                ! go in time
        else 
          exit                                             ! the next halo is not on the same branch, stop the run
        end if
      end do
      !
      ! scan progenitors
      ifirst = univ(ts)%tree(ih)%firstProg                 ! index of the first progenitor of the halo ih at ts
      ts2    = ts -1                                       ! time-step of the main progenitor
      do while (ifirst .gt. 0)
        !
        if (univ(ts2)%tree(ifirst)%level .eq. 1) then      ! only main-halos 
          n_halos_in_branch = n_halos_in_branch + 1        ! update the number of main halo (good and bad)
          ! update the number of good identified main-halo
          if (univ(ts2)%halo(ifirst)%dm%compute) n_good_halos = n_good_halos + 1                
        end if
        me     = ifirst
        ifirst = univ(ts2)%tree(me)%firstProg              ! go to the next progenitor (on the same branch)
        ts2    = ts2 - 1                                   ! go back in time
      end do
      !
      if (n_halos_in_branch .ne. 0) then
         !
         tree_compute_dm_goodness_of_branch = real(n_good_halos,8)/real(n_halos_in_branch,8)
      else
         !
         tree_compute_dm_goodness_of_branch = 0.
      end if
      
      return
    end function tree_compute_dm_goodness_of_branch

    ! *************************************************

    subroutine tree_set_dm_goodness_of_branch(ih,ts,goodness)

      ! SET (dm%goodness_of_branch) FOR EACH HALO IN THE CURENT BRANCH THE (dm-)GOODNESS OF ITS CURRENT BRANCH 

      implicit  none

      integer(kind=4),intent(in)  :: ts              ! the current time-step
      integer(kind=4),intent(in)  :: ih              ! the current halo
      integer(kind=4)             :: me,idesc,ifirst ! some indexes for: a halo, a descendent, a first prog
      integer(kind=4)             :: ts2             ! a local time-step 

      real(kind=8),intent(in)     :: goodness        ! number of (dm-)good main-halo/number of main halo in the current branch

      ! start with the current halo (ih at ts)
      univ(ts)%halo(ih)%dm%goodness_of_branch = goodness
      ! 
      ! scan descendents  
      idesc = univ(ts)%tree(ih)%descendent           ! index of my descendent
      ts2   = ts + 1                                 ! time-step of my descendent
      me    = ih                                     ! my index
      
      do while (idesc > 0)
        ! 
        if (univ(ts2)%tree(idesc)%firstProg .eq. me) then ! the descendent is still on the main branch (i.e. main prog of desc)
          !
          univ(ts2)%halo(idesc)%dm%goodness_of_branch = goodness
          me    = idesc                              ! go to my next descendent (on the same branch)
          idesc = univ(ts2)%tree(me)%descendent          
          ts2   = ts2 + 1                            ! go in time
        else 
          exit  ! th descendent is not still on the same branch
        end if
      end do
      ! 
      ! scan progenitors
      ifirst = univ(ts)%tree(ih)%firstProg           ! index of my first progenitor
      ts2    = ts -1                                 ! time-step of my fist progenitor
      
      do while (ifirst .gt. 0)
        !
        univ(ts2)%halo(ifirst)%dm%goodness_of_branch = goodness
        me     = ifirst                              ! go to my next progenitor (on the same branch)
        ifirst = univ(ts2)%tree(me)%firstProg             
        ts2    = ts2 - 1                             ! go back in time
      end do

      return
    end subroutine tree_set_dm_goodness_of_branch

    ! *************************************************

    function tree_halo_become_good(ih,ts) 

      ! CHECK RESCUE CONDITION FOR A GIVEN HALO
      ! return .true. is the halo must be saved

      implicit none

      integer(kind=4),intent(in)  :: ts                           ! the current time-step
      integer(kind=4),intent(in)  :: ih                           ! the current halo
      integer(kind=4)             :: me,idesc,ifirst              ! some indexes for: a halo, a descendent, a first prog
      integer(kind=4)             :: ts2                          ! a local time-step
      integer(kind=4)             :: nb_of_desc_halo_already_test ! number of descendent halos already tested 

      logical                     :: tree_halo_become_good

      tree_halo_become_good = .false.                  ! in input, the the halo is a BAD halo
      !
      ! if there aren't enouth future time step for test halos, we keep the current halo
      if (ts .gt. (nsteps-nb_of_desc_halo_test)) then
        !
        tree_halo_become_good = .true.
        return
      endif 
      !
      idesc = univ(ts)%tree(ih)%descendent             ! the first descendent of the current halo (halo ih at ts)
      ts2   = ts +1                                    ! the time-step of the firts descendent
      me    = ih              
      !
      nb_of_desc_halo_already_test  = 0                ! init the number of descendent halo already tested
      !
      do while ((idesc .gt. 0) .and. (ts2 .le. nsteps) .and. (nb_of_desc_halo_already_test .lt. nb_of_desc_halo_test))
        !
        nb_of_desc_halo_already_test = nb_of_desc_halo_already_test +1
        ifirst = univ(ts2)%tree(idesc)%firstProg   
        if (ifirst .ne. me) return                     ! if i'm not the main prog of my descendent, screw me ! 
        !
        if (univ(ts2)%tree(idesc)%level .gt. 1) then   ! if my desc is a sub halo, keep me !
          !
          tree_halo_become_good = .true.
          return
        else
          !
          ! my descendent is a main-halo
          if (univ(ts2)%halo(idesc)%dm%compute) then      
            !
            tree_halo_become_good = .true.             ! my descendent is a halo with good dark matter properties, keep me !
            return 
          end if
        end if
        me    = idesc                                  ! go to the next descendent (on the same branch)
        idesc = univ(ts2)%tree(idesc)%descendent         
        ts2   = ts2 +1                                 ! go in time
      end do
                                                                              
      return    
    end function tree_halo_become_good

    ! *************************************************

    function tree_compute_halo_goodness_of_branch(ih,ts)

      ! COMPUTE THE NUMBER OF (halo-)GOOD MAIN-HALO / NUMBER OF MAIN HALO IN THE CURRENT BRANCH

      implicit none

      integer(kind=4),intent(in)  :: ts                                   ! the current time-step
      integer(kind=4),intent(in)  :: ih                                   ! the current halo
      integer(kind=4)             :: me,idesc,ifirst                      ! some indexes for: a halo, a descendent, a first prog
      integer(kind=4)             :: ts2                                  ! a local time-step  
      integer(kind=4)             :: n_halos_in_branch                    ! number of halos in the current branch
      integer(kind=4)             :: n_good_halos                         ! number_of good halos in the current branch
      
      real(kind=8)                :: tree_compute_halo_goodness_of_branch ! number of (dm-)good main-halo/number of main halo in the current branch
      
      n_good_halos      = 0  ! init 
      n_halos_in_branch = 0  ! init 

      ! check for the current halo (the halo ih at ts)
      ! a main halo is caracerised by level = 1
      if (univ(ts)%tree(ih)%level .eq. 1) then
        !
        n_halos_in_branch = 1
        ! if the current halo is a (halo-)computable halo
        ! this floag has been computed previously by: dm_flag_bad_properties in tree_compute_halo_derived_properties
        if (univ(ts)%halo(ih)%compute) n_good_halos = 1
      end if
      !
      ! scan descendents
      idesc = univ(ts)%tree(ih)%descendent                 ! index of my (halo ih at ts) descendent at the next timestep
      ts2   = ts + 1                                       ! time-step of my descendent
      me    = ih                                           ! my index
      
      do while (idesc .gt. 0)
        !
        if (univ(ts2)%tree(idesc)%firstProg .eq. me) then  ! the descendent is still on main branch (i.e. main prog of desc)
          !
          if (univ(ts2)%tree(idesc)%level .eq. 1) then     ! the halo is a main halo
            !
            n_halos_in_branch = n_halos_in_branch + 1      ! update the number of main-halo in the current branch (good and bad)  
            ! update the number of good identified main-halo
            if (univ(ts2)%halo(idesc)%compute) n_good_halos = n_good_halos + 1                
          end if
          me    = idesc
          idesc = univ(ts2)%tree(me)%descendent            ! go to the next descendent (on the same branch)
          ts2   = ts2 + 1                                  ! go in time
        else 
          exit                                             ! the next halo is not on the same branch, stop the run
        end if
      end do
      !
      ! scan progenitors
      ifirst = univ(ts)%tree(ih)%firstProg                 ! index of the first progenitor of the halo ih at ts
      ts2    = ts -1                                       ! time-step of the main progenitor
      do while (ifirst .gt. 0)
        !
        if (univ(ts2)%tree(ifirst)%level .eq. 1) then      ! only main-halos 
          !
          n_halos_in_branch = n_halos_in_branch + 1        ! update the number of main halo (good and bad)
          ! update the number of good identified main-halo
          if (univ(ts2)%halo(ifirst)%compute) n_good_halos = n_good_halos + 1                
        end if
        me     = ifirst
        ifirst = univ(ts2)%tree(me)%firstProg              ! go to the next progenitor (on the same branch)
        ts2    = ts2 - 1                                   ! go back in time
      end do
      !
      if (n_halos_in_branch .ne. 0) then
         !
         tree_compute_halo_goodness_of_branch = real(n_good_halos,8)/real(n_halos_in_branch,8)
      else
         !
         tree_compute_halo_goodness_of_branch = 0.
      end if

      return
    end function tree_compute_halo_goodness_of_branch

    ! *************************************************

    subroutine tree_set_halo_goodness_of_branch(ih,ts,goodness)

      ! SET (halo%goodness_of_branch) FOR EACH HALO IN THE CURENT BRANCH THE (halo-)GOODNESS OF ITS CURRENT BRANCH 

      implicit  none

      integer(kind=4),intent(in)  :: ts              ! the current time-step
      integer(kind=4),intent(in)  :: ih              ! the current halo
      integer(kind=4)             :: me,idesc,ifirst ! some indexes for: a halo, a descendent, a first prog
      integer(kind=4)             :: ts2             ! a local time-step 

      real(kind=8),intent(in)     :: goodness        ! number of (dm-)good main-halo/number of main halo in the current branch

      ! start with the current halo (ih at ts)
      univ(ts)%halo(ih)%goodness_of_branch = goodness
      ! 
      ! scan descendents  
      idesc = univ(ts)%tree(ih)%descendent           ! index of my descendent
      ts2   = ts + 1                                 ! time-step of my descendent
      me    = ih                                     ! my index
      do while (idesc > 0)
        !
        if (univ(ts2)%tree(idesc)%firstProg .eq. me) then ! the descendent is still on the main branch (i.e. main prog of desc)
          !
          univ(ts2)%halo(idesc)%goodness_of_branch = goodness
          me    = idesc                              ! go to my next descendent (on the same branch)
          idesc = univ(ts2)%tree(me)%descendent          
          ts2   = ts2 + 1                            ! go in time
        else 
          exit  ! th descendent is not still on the same branch
        end if
      end do
      ! 
      ! scan progenitors
      ifirst = univ(ts)%tree(ih)%firstProg           ! index of my first progenitor
      ts2    = ts -1                                 ! time-step of my fist progenitor
      
      do while (ifirst .gt. 0)
        !
        univ(ts2)%halo(ifirst)%goodness_of_branch = goodness
        me     = ifirst                              ! go to my next progenitor (on the same branch)
        ifirst = univ(ts2)%tree(me)%firstProg             
        ts2    = ts2 - 1                             ! go back in time
      end do

      return
    end subroutine tree_set_halo_goodness_of_branch
    
    ! *************************************************

    subroutine tree_computed_and_set_n_computed_dads(ih,ts)

      ! COMPUTED, FOR A GIVEN HALO, THE NUMBER OF COMPUTABLE PROGENITORS

      implicit none
      
      integer(kind=4),intent(in)  :: ts     ! the current time-step
      integer(kind=4),intent(in)  :: ih     ! the current halo
      integer(kind=4)             :: iprog  ! index of a progenitor
      
      univ(ts)%tree(ih)%n_computed_dads = 0 ! The number of computable halos is initially set to 0

      iprog = univ(ts)%tree(ih)%firstProg   ! Initialize to the first progenitor

      do while (iprog .gt. 0) 
        !
        if (univ(ts-1)%halo(iprog)%compute) then 
          !
          ! only computable halo are take into account
          univ(ts)%tree(ih)%n_computed_dads = univ(ts)%tree(ih)%n_computed_dads + 1
        end if
        iprog = univ(ts-1)%tree(iprog)%nextProg  ! go the next progenitor (at the same time-step)
      end do

      return
    end subroutine tree_computed_and_set_n_computed_dads
    
    ! *************************************************

    subroutine tree_interpolate_bad_properties(ih,ts)

      ! INTERPOLATE dm-PROPERTIES (here only spin) BETWEEN TWO TIME_STEPS OF COMPUTED HALO IF THIS ONE HASN'T WELL DEFINED dm-PROPERTIES

      implicit none
      
      integer(kind=4),intent(in)   :: ts,ih                ! the time-step and the halo id
      integer(kind=4)              :: idesc,iprog,ts2,me   ! local id

      logical                      :: find_prog, find_desc ! research flags

      real(kind=8)                 :: spin_prog,spin_desc  ! spin of the first good progenitor and first good descendent.
      real(kind=8)                 :: z_desc,z_prog        ! redshift of the first good progenitor and first good descendent
      real(kind=8)                 :: a,b                  ! interpolation coefficients
      
      ! down in the tree along the main prog for find a halo with well defined DM props.
      find_prog = .false.
      iprog = univ(ts)%tree(ih)%firstProg
      ts2 = ts -1
      
      if (iprog .gt. 0) then
         !
         spin_prog = log10(univ(ts2)%halo(iprog)%dm%spin)
         z_prog    = 1./univ(ts2)%aexp - 1.
         find_prog = .true.
      end if
      
      ! up in the tree along the desc for find a halo with well defined DM props.
      find_desc = .false.
      me        = ih 
      idesc     = univ(ts)%tree(me)%descendent
      ts2       = ts + 1
      do while (idesc .gt. 0)
         !
         ! check that halo is main prog of its descendent to be used (or not) for interpolation
         if (univ(ts2)%tree(idesc)%firstProg .ne. me) exit
         !
         ! if descendent is in the same branch, use it if its properties are well defined
         if (univ(ts2)%halo(idesc)%dm%compute) then
            !
            spin_desc = log10(univ(ts2)%halo(idesc)%dm%spin)
            z_desc    = 1./univ(ts2)%aexp - 1.
            find_desc = .true.
            exit
         else 
            !
            ! if desc not well defined move to next
            me    = idesc
            idesc = univ(ts2)%tree(me)%descendent 
            ts2   = ts2 +1
         end if
      end do
      !
      if (find_prog) then
         !
         if (find_desc) then
            !
            a = (spin_desc - spin_prog)/(z_desc - z_prog)
            b = spin_prog - a*z_prog
            univ(ts)%halo(ih)%dm%spin = 10.**(a*(1./univ(ts)%aexp - 1.d0) + b)  ! erase the previous value
         else
            !
            univ(ts)%halo(ih)%dm%spin = 10.**(spin_prog)
         end if
      else
         !
         if (find_desc) then
            !
            univ(ts)%halo(ih)%dm%spin = 10.**(spin_desc)
         else
            !
            ! neither valid prog nor desc in the whole tree -> assign mean value.
            if (univ(ts)%tree(ih)%level .eq. 1) then
               !
               univ(ts)%halo(ih)%dm%spin = 0.056d0      ! for main halos
            else
               !
               univ(ts)%halo(ih)%dm%spin = 0.051d0      ! for sub-halos
            endif
         end if
      end if
      
      if (univ(ts)%halo(ih)%dm%spin .le. 0.d0) then
         !
         call IO_print_error_message('spin = 0',only_rank=rank,called_by='tree_interpolate_bad_properties')
         stop ! stop the program 
      end if
      
      return
    end subroutine tree_interpolate_bad_properties
       
    ! *************************************************
    ! end contains tree_clean

  end subroutine tree_clean
  
  !*****************************************************************************************************************

  recursive subroutine tree_kill_branch(ih,ts)

    ! KILL (FLAG halo%compute = .false.) ALL HALOS IN A GIVEN BRANCH

    implicit none

    integer(kind=4),intent(in)  :: ts        ! the current time-step
    integer(kind=4),intent(in)  :: ih        ! the current halo
    integer(kind=4)             :: ifirst    ! index of the first prog
      
    univ(ts)%halo(ih)%compute          = .false.
    univ(ts)%halo(ih)%followed         = .false.
    univ(ts)%halo(ih)%recursively_kill = .true.
    !
    ! go to te first prog
    ifirst  = univ(ts)%tree(ih)%firstProg
    do while (ifirst .gt. 0)
      !
      if (univ(ts-1)%halo(ifirst)%compute .and. (.not. univ(ts-1)%halo(ifirst)%recursively_kill)) then
        !
        call tree_kill_branch(ifirst,ts-1)
      end if
      ifirst  = univ(ts-1)%tree(ifirst)%nextProg           ! go to the next progenitors
    end do
    
    return
  end subroutine tree_kill_branch
  
  !*****************************************************************************************************************

  recursive subroutine tree_followed_branch(ih,ts)

    ! SELECT ALL HALOS IN A GIVEN BRANCH WITH FLAG (halo%followed = .true.) 

    implicit none

    integer(kind=4),intent(in)  :: ts        ! the current time-step
    integer(kind=4),intent(in)  :: ih        ! the current halo
    integer(kind=4)             :: ifirst    ! index of the first prog
      
    if (univ(ts)%halo(ih)%compute) univ(ts)%halo(ih)%followed = .true.
    !
    ! go to the first prog
    ifirst  = univ(ts)%tree(ih)%firstProg
    do while (ifirst .gt. 0)
      !
      if (univ(ts-1)%halo(ifirst)%compute) then
        !
        call tree_followed_branch(ifirst,ts-1)
      end if
      ifirst  = univ(ts-1)%tree(ifirst)%nextProg  ! go to the next progenitors
    end do

    return
  end subroutine tree_followed_branch
  
  !*****************************************************************************************************************
  
  subroutine tree_deallocate
    
    ! ERASE ALL INFORMATION ABOUT THE TREES AND THEIR CONTENT TO FREE MEMORY
    
    implicit none 
    
    integer(kind=4) :: ts,ih
    
    do ts = 1,nsteps
       !
       if (univ(ts)%nb_of_halos > 0) then                           ! if there are any haloes at the timestep ts
          !
          if (allocated(univ(ts)%tree)) deallocate(univ(ts)%tree)   ! deallocate the tree information
          !
          do ih = 1,univ(ts)%nb_of_halos
             call halo_deallocate(univ(ts)%halo(ih))                ! delete the information inside each halo
          end do
          !
          if (allocated(univ(ts)%halo)) deallocate(univ(ts)%halo)   ! deallocate the halo array
          if (allocated(univ(ts)%tree)) deallocate(univ(ts)%tree)   ! deallocate the tree array
       end if
    end do
    if (allocated(univ)) deallocate(univ)                           ! deallocate the container data structure
    
    return  
  end subroutine tree_deallocate

  !*****************************************************************************************************************

  subroutine tree_evolve
    
    ! COMPUTE THE EVOLUTION OF THE BARYONIC COMPONENT OF ALL HALOS IN A GIVEN TREE

    implicit none
    
    integer(kind=4)              :: i                                     ! a loop index
    integer(kind=4)              :: ts                                    ! the current time-step
    integer(kind=4)              :: ots                                   ! index in the output timestep list
    integer(kind=4)              :: ih                                    ! the current halo
    integer(kind=4)              :: prog                                  ! index for a progenitor
    integer(kind=4)              :: idesc                                 ! index for a progenitor
    integer(kind=4)              :: upper_level_of_substruct              ! the maximal level of substructure in the main halo loop (at ts)
    integer(kind=4)              :: upper_level_of_substruct2             ! the maximal level of substructure in the halo progenitor loop (at ts-1)
    integer(kind=4)              :: l,l2                                  ! level (main/sub) indexes loop 
    integer(kind=4)              :: n_computable_prog_computed            ! number of computable progenitors already compute   
    integer(kind=4)              :: nb_current_follow_up_halos            ! number of HID associated to the current computed ifile
    integer(kind=4)              :: current_follow_up_halo_index          ! current index in the follow_up_halo_list
    integer(kind=4),allocatable  :: ih_follow_up(:), prog_follow_up(:)    ! list of halo index for follow_up halos
    integer(kind=4),allocatable  :: ts_follow_up(:)
    integer(kind=4),allocatable  :: ts_stop_follow_up(:)                  ! list of last time-step of the current follow_up halo
    integer(kind=4),allocatable  :: current_unit_list(:) 

    logical                      :: OUTPUT_TS                             ! = .true. if the current ts is in the output timestep list 

    character(MAXPATHSIZE)       :: message                               ! a message 
    
    real(kind=8)                 :: zu, zd, z, zref, dz                   ! redshift zu = z(ts) and zd = z(ts -1), z = z(t)
    real(kind=8)                 :: m1
    real(kind=8)                 :: dt1,dt2                               ! timestep duration before and after the merger process
    real(kind=8)                 :: tu, td                                ! age of the univers td = t(ts-1) and tu = t(ts) 
    
    type(halo_type)              :: h                                     ! a tmp halo
        
#ifdef PRINTALL
! -------------------------------------------------
    call IO_print_message('tree_evolve',only_rank=rank,component='tree')
! -------------------------------------------------
#endif
    
    halo_dt_min_use  = StellarTimeStep                                      ! init
    disc_dt_min_use  = StellarTimeStep                                      ! init
    bulge_dt_min_use = StellarTimeStep                                      ! init
    dt_min_use = min(halo_dt_min_use,min(disc_dt_min_use,bulge_dt_min_use)) ! init
    ! 
    ! ********************************
    ! SELECT CURRENT FOLLOWED HALOS
    ! ********************************
    !
    ! count the number oh HID associated to the current computed ifile
    nb_current_follow_up_halos = 0
    
    do ih = 1, nb_follow_up_halos
        !
        if (IO_return_ifile_from_HID(list_of_follow_up_halo(ih)) .eq. IO_return_ifile_from_HID(minval(univ(nsteps)%tree(:)%HID))) then 
            !
            nb_current_follow_up_halos = nb_current_follow_up_halos +1
        end if
    end do
    
    if (nb_current_follow_up_halos .gt. 0) then
        !
        ! create follow_up halo index list
        allocate(ih_follow_up(nb_current_follow_up_halos))             ! the current halo index at ts_follow_up
        allocate(prog_follow_up(nb_current_follow_up_halos))           ! the index of the main progenitor for the halo ih_follow_up
        allocate(ts_follow_up(nb_current_follow_up_halos))             ! the current ts of the halo ih_follow_up
        allocate(ts_stop_follow_up(nb_current_follow_up_halos))        ! the last time step follow_up evolution for the halo ih_follow_up
        allocate(current_unit_list(nb_current_follow_up_halos))        ! list of indexes available in the global follow_up list
        i = 0
        
        do ih = 1, nb_follow_up_halos
            !
            if (IO_return_ifile_from_HID(list_of_follow_up_halo(ih)) .eq. IO_return_ifile_from_HID(minval(univ(nsteps)%tree(:)%HID))) then 
                !
                i = i+1
                ! associate the current index with the real index in the global follow_up list
                current_unit_list(i) = ih
                ! generate the ts_stop_follow_up list
                ts_stop_follow_up(i) = IO_return_ts_from_HID(list_of_follow_up_halo(ih))
                ! generate the ih_follow_up list (indexes of the last computable progenitor of the main)
                call last_computable_progenitor_of_a_main_branch(ts_stop_follow_up(i),IO_return_ih_from_HID(list_of_follow_up_halo(ih)), &
                    ts_follow_up(i),ih_follow_up(i))
                ! generate prog_follow_up list
                ! @ the fist time step, the follow_up halo has no progenitor
                prog_follow_up(i) = -1  ! init
            end if
        end do
    end if
    !
    ! init ots
    ots = 1
    !
    ! ********************************
    ! EVOLUTION
    ! ********************************
    !
    call IO_print_message('Tree[Evolve]: Start',only_rank=rank,component='tree')
    !
    do ts = 1, ts_STOP
        !
        write(message,'(a,i3.3,a,3(i6.6,a))') 'ts = ', ts, ', [computed, followed, tot] = [', &
                univ(ts)%nb_of_computable_halos, ', ', univ(ts)%nb_of_followed_up_halos, ', ', univ(ts)%nb_of_halos, ']' 
        call IO_print_message(message,component='tree',only_rank=rank)
        !
        ! set halo%evolve keyword of progenitor halos to false
        ! allow to control if a main structure has been evolved before a sub-halo mass transfert
        if (univ(ts-1)%nb_of_halos .gt. 0) then
            !
            univ(ts-1)%halo(:)%evolved = .false.
        end if    
        !
        ! init output keyword
        OUTPUT_TS = .false.
        !
        ! Check if the current ts is in the output timestep list
        if (ts .eq. tsout(ots)) then
            !
            ! the current time-step is an output timestep
            ! turn on keyword
            OUTPUT_TS = .true.
            ! print the number of computable halo registered in this tree-file at this time-step
            write(IO_get_temporary_output_files_unit(ts)) univ(ts)%nb_of_computable_halos
        end if        
        !
        ! redshifts   
        if (ts .gt. 1) zd = 1./univ(ts-1)%aexp - 1.    ! redshift of snapshot ts -1 
        zu = 1./univ(ts)%aexp - 1.                     ! redshift of snapshot ts
        if (ts .gt. 1) td  = univ(ts-1)%age_univ       ! age of the universe at ts-1 
        tu  = univ(ts)%age_univ                        ! age of the universe at ts 
        !
        if (univ(ts)%nb_of_computable_halos .gt. 0) then
            ! 
            ! run throught structure level
            upper_level_of_substruct = maxval(univ(ts)%tree(:)%level) ! 
            ! If gas is transfered from higher structure level to main structure, 
            ! evolution loop has to run throught sub-structures before computre evolution of main halos
            do l = upper_level_of_substruct, 1, -1
                ! 
                ! run throught halo at the current timestep  
                do ih = univ(ts)%nb_of_halos, 1, -1
                    !
                    if (univ(ts)%halo(ih)%compute .and. (univ(ts)%tree(ih)%level .eq. l)) then
#ifdef TREE_EVOLVE
! -------------------------------------------------                         
                        ! @ this point the halo is a computable halo and it level is equal to the current computed level
                        !
                        FOLLOW_UP    = .false.  ! init, by default the evolution of the current halo is not followed
                        PR_FOLLOW_UP = .false.  ! init  by default the evolution of progenitors of the current halo is not followed
                        if (nb_current_follow_up_halos .gt. 0) then
                            !
                            ! loop under the follow_up list
                            do i = 1, nb_current_follow_up_halos
                                !
                                if ((ts .eq. ts_follow_up(i)) .and. (ih .eq. ih_follow_up(i))) then
                                    current_follow_up_halo_index = i
                                    current_index                = current_unit_list(i)
                                    ! for print in the file, FOLLOW_UP AND PR_FOLLOW_UP must be .true.
                                    FOLLOW_UP                    = .true.  ! turn ON
                                    PR_FOLLOW_UP                 = .true.  ! turn ON
                                    exit
                                end if
                            end do
                        end if                    
                        !
#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------     
                        if (OUTPUT_TS) then       
                            !
                            ! send tree(ih) information (redshift, HID, host_HID ...) to the associated luminous process
                            call tree_send_data(ts=ts,ih=ih)    
                        end if
! ------------------------------------------------- 
#endif
! LUMINOUS_PROCESSES
                        !                     
                        if (univ(ts)%tree(ih)%n_computed_dads .eq. 0) then                  ! no computable progenitor  
                            !
                            ! *****************************************************************
                            ! NO PROGENITOR
                            ! *****************************************************************
                            !
                            ! Evolve the structure passively for the entire timestep dt2
                            ! dt2 is the time laps between the formation age of the dm-halo and the current main time-step
                            dt2 = tu - univ(ts)%halo(ih)%dm%age_form 
                            !         
                        else
                            !
                            ! *****************************************************************
                            ! MERGER OR ONE PROGENITOR
                            ! *****************************************************************
                            !
                            dt1 = univ(ts)%halo(ih)%dm%age_merge  - td   ! compute time laps before the merging time
                            dt2 = tu - univ(ts)%halo(ih)%dm%age_merge    ! compute time laps after the merging time
                            !
                            ! Evolve all progenitors until the merging time 
                            m1 = (1.d0 - dt1/(abs(tu-td)))
                            zref = m1*zd + (1.d0-m1)*zu
                            dz = dt2dz(dt1,zref)      ! negative
                            z = zd + dz               ! until what redshift the halo must grow
                            !
                            if (dz .ge. 0.) then
                                !
                                call IO_print_error_message('dz >= 0',only_rank = rank, called_by = 'tree_evolve') 
                                call IO_print_message('  ',only_rank=rank,component='tree', &
                                    param_name = (/'zd                       ','z                        ','zu                       ','dz                       ', &
                                                   'dt1                      ','dt2                      ','td                       ', &
                                                   'tm                       ','tu                       '/), &
                                    real_param_val  = (/zd,z,zu,dz,dt1,dt2,td,univ(ts)%halo(ih)%dm%age_merge,tu/))
                                stop
                            end if
                            !
                            n_computable_prog_computed = 0                                    ! we set to 0 the number of computable progenitors
                            ! 
                            ! run throught progenitors levels
                            upper_level_of_substruct2  = maxval(univ(ts-1)%tree(:)%level)     ! we set the max level of structure (in the progenitor timestep ts-1)
                            !
                            do l2 = upper_level_of_substruct2, 1, -1
                                !   
                                ! run throught progenitors
                                prog = univ(ts)%tree(ih)%firstProg                              ! we star the evolution by the main progenitor
                                do while (prog .gt. 0)   
                                    !
                                    ! As for main evolution loop through halos, the progenitors are evolved by following substructure level
                                    ! First, the highest substructure level is evolved  
                                    if (univ(ts-1)%halo(prog)%compute .and. (univ(ts-1)%tree(prog)%level .eq. l2)) then  
                                        !
                                        ! @ this point, the current progenitor (prog) is computable and has a level equal to the current computed level: l2 
                                        !
                                        ! reset to false the progenitor keyword, by default the current progenitor is not followed                                                                                                    
                                        PR_FOLLOW_UP = .false.                                      
                                        if (current_follow_up_halo_index .gt. 0) then
                                            !
                                            if ((FOLLOW_UP) .and. (prog_follow_up(current_follow_up_halo_index) .gt. 0)) then
                                                !
                                                ! the FOLLOW_UP keyword is turn ON
                                                ! for print in file "FOLLOW_UP" AND "PR_FOLLOW_UP" must be .true.
                                                if (prog .eq. prog_follow_up(current_follow_up_halo_index)) PR_FOLLOW_UP = .true. 
                                            end if
                                        end if  
                                        !
                                        ! HALO EVOLVE
                                        ! ***********************
                                        ! Evolve the progenitor until the merging time
                                        ! create a local copy of the progenitor
                                        call halo_void(h)
                                        call halo_save(h,univ(ts-1)%halo(prog))
                                        call tree_evolve_halo(h,ts-1,prog,z,dt1,post_merger=.false.)
                                        ! ***********************
                                        !
                                        ! HALO MERGE
                                        ! ***********************
                                        ! merge baryonic component in the descendent
                                        call halo_merge(univ(ts)%halo(ih),h)  
                                        call halo_void(h)                        ! erase the tmp halo 
                                        call halo_deallocate(h)                  ! deallocate the tmp halo 
                                        ! ***********************
                                        !
                                        ! add a progenitor to the computed list
                                        n_computable_prog_computed = n_computable_prog_computed +1  
                                        !
                                        if (n_computable_prog_computed .eq. univ(ts)%tree(ih)%n_computed_dads) exit ! all computable progenitors have been computed
                                    end if
                                    ! go to the next prog
                                    prog = univ(ts-1)%tree(prog)%nextProg
                                end do ! while loop under progenitor
                                !
                                if (n_computable_prog_computed .eq. univ(ts)%tree(ih)%n_computed_dads) exit ! all computable progenitors have been computed
                                !
                            end do ! do loop under sub-structure level
                            !
                            if (n_computable_prog_computed .ne. univ(ts)%tree(ih)%n_computed_dads) then
                                !
                                call IO_print_error_message('Evolution scheme does not evolve the right number of progenitors !', &
                                    only_rank = rank, called_by = 'tree_evolve') 
                                write(errunit,*) ' --> With prog_computed : ', n_computable_prog_computed, &
                                        ', n_computed_dads : ', univ(ts)%tree(ih)%n_computed_dads       
                                stop
                            end if
                            !
                            if (FOLLOW_UP) PR_FOLLOW_UP = .true.  ! for print in file FOLLOW_UP AND PR_FOLLOW_UP must be true
                            !
                        end if ! MERGER OR ONE PROGENITOR
                        !
                        ! **********************
                        ! EVOLVE THE REMNANT HALO
                        ! **********************
                        !
                        z = zu
                        !
                        ! HALO EVOLVE
                        ! ***********************
                        ! Evolve the merger remnant until the end of the timestep 
                        call tree_evolve_halo(univ(ts)%halo(ih),ts,ih,z,dt2,post_merger=.true.)
                        ! ***********************
                        !
#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------               
                        ! send halo, galaxy, stars data to the associated luminous process
                        if (OUTPUT_TS) call halo_send_data(univ(ts)%halo(ih),z,post_merger=.true.)  
! ------------------------------------------------- 
#endif
! LUMINOUS_PROCESSES 
                        !              
                        if (FOLLOW_UP) then
                            !
                            ! if FOLLOW_UP ts = ts_follow_up(current_follow_up_halo_index)
                            !
                            ! print tree information
                            call tree_print(follow_up_unit(current_index),'fits','phy',univ(ts),ih=ih)      
                            !
                            ! for merger tree building
                            if (ts_follow_up(current_follow_up_halo_index) .eq. ts_stop_follow_up(current_follow_up_halo_index)) then
                                !
                                ! the tree evolution is complete, we must print the merger-tree in the corresponding output-file
                                call tree_print_tree(ts_stop_follow_up(current_follow_up_halo_index), &
                                        ih_follow_up(current_follow_up_halo_index),first_call=.true.)
                            else
                                !
                                ! The evolution of the current followed halo is not complete,
                                ! update prog_follow_up, ih_follow_up and ts_follow_up
                                ! the current followed halo, become the progenitor
                                prog_follow_up(current_follow_up_halo_index) = ih_follow_up(current_follow_up_halo_index)
                                !
                                ! the next follup_up halo, is obviously the descendent 
                                idesc = univ(ts)%tree(ih_follow_up(current_follow_up_halo_index))%descendent
                                if (idesc .gt. 0) then
                                    !
                                    ! halo has a descendent !
                                    if (univ(ts+1)%halo(idesc)%compute) then
                                        !
                                        ! the descendent halo is computable
                                        ih_follow_up(current_follow_up_halo_index) = idesc
                                    else
                                        !
                                        ! the descendent is not computable
                                        ih_follow_up(current_follow_up_halo_index) = -1 ! set to null value
                                    end if
                                else  
                                    !
                                    ! the ih halo has not descendent
                                    ih_follow_up(current_follow_up_halo_index) = -1   ! set to null value
                                end if 
                                ts_follow_up(current_follow_up_halo_index) = ts +1  ! go to next timestep
                                !
                            endif
                        end if 
                        !
! -------------------------------------------------
#endif
! TREE_EVOLVE    
                        if (OUTPUT_TS) then             
                            !
                            ! print tree, halo, galaxy data in tmp file
                            call tree_print(IO_get_temporary_output_files_unit(ts),'tmp_bin','phy',univ(ts),ih=ih,go_down=.true.)                                                     
                        end if ! OUTPUT_TS                   
                    end if ! univ(ts)%halo(ih)%compute .and. univ(ts)%tree(ih)%level .eq. l
                end do !end ih loop
            end do ! end level loop 
            !
            ! Erase past tree to save memory... 
            if (ts .gt. 1) then 
                !
                if (univ(ts-1)%nb_of_halos .gt. 0) then 
                    !
                    do ih = 1,univ(ts-1)%nb_of_halos
                        !
                        call halo_deallocate(univ(ts-1)%halo(ih))   ! All the progenitors are deallocated (one by one)
                    end do
                end if
            end if
        end if ! univ(ts)%nb_of_computable_halos .gt. 0
        !       
        ! send z data
        call tree_send_z_data(univ(ts))
        !
        ! update ots
        if (OUTPUT_TS) ots = min(ots +1,nout)
    end do ! ts loop
    !
    call IO_print_message('Tree[Evolve]: End',only_rank=rank,component='tree')
    !   
    if (allocated(ih_follow_up))      deallocate(ih_follow_up)      ! the current halo index at ts_follow_up
    if (allocated(prog_follow_up))    deallocate(prog_follow_up)    ! the index of the main progenitor for the halo ih_follow_up
    if (allocated(ts_follow_up))      deallocate(ts_follow_up)      ! the current ts of the halo ih_follow_up
    if (allocated(ts_stop_follow_up)) deallocate(ts_stop_follow_up) ! the last time step follow_up evolution for the halo ih_follow_up
    !
#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------         
    call tree_send_data ! a last call --> stop the loop   
    call IO_print_message('tree[send]: LAST CALL !!',only_rank=rank,component='tree')   
! ------------------------------------------------- 
#endif
! LUMINOUS_PROCESSES                
    !
    !
    ! deallocate arrays (stars ...) to keep memory space
    call tree_deallocate 
    !
    return
    
    contains 

      ! *************************************************
      
      subroutine tree_evolve_halo(h,ts,ih,z,dt,post_merger)
      
        ! EVOLVE THE HALO "ih" evolving @ ts during dt
        ! this routine allows to take into account all model option (quenching, gas transfert)
        
        implicit none
      
        integer(kind=4),intent(in)    :: ts            ! the time-step
        integer(kind=4),intent(in)    :: ih            ! the halo
        integer(kind=4)               :: host          ! host halo index
        
        logical,intent(in)            :: post_merger   ! .true. if the evolution is done onto the remanent structure of a merger
        
        real(kind=8),intent(in)       :: dt            ! the time-step
        real(kind=8),intent(in)       :: z             ! the redshift at the end of the timestep
        
        type(halo_type),intent(inout) :: h             ! the halo

        ! by default the quenching process is turned OFF
        h%baryon_halo%quenched = .false.   
                                      
#ifdef SUB_QUENCHING
! -------------------------------------------------
        if (univ(ts)%tree(ih)%level .gt. 1) then
          !
          ! the halo is a sub-structure
#ifdef UNLINKED_TREES
! -------------------------------------------------
          ! the model runs on a set of unlinked merger trees 
          ! i.e no link between sub-ahlo and main halo
          h%baryon_halo%quenched = .true.
! ------------
#else
! ------------  
          ! th emodel runs on a set of linked-merger-trees, 
          ! i.e. connexions between halo and sub-halo is available.
          ! find the host halo
          host  = univ(ts)%tree(ih)%hosthalo                             
          if (host .gt. 0) then
            !
            ! The sub-halo is linked to a main halo
            ! Is the main halo "computable" ?
            if (univ(ts)%halo(host)%compute) then
                !
                ! The quenching process is turned ON
                h%baryon_halo%quenched = .true.
            end if
          end if
! -------------------------------------------------
#endif
! UNLINKED_TREES        
        end if   
! -------------------------------------------------                
#endif 
! SUB_QUENCHING 
        !
#ifdef SUB_TRANSFER
! -------------------------------------------------
        ! For sub-structure, the gas ejected of the gravitational may be transferred to the host main halo
        if (univ(ts)%tree(ih)%level .gt. 1) then
           !
           ! the halo is a sub-structure
           ! find the host halo
           host  = univ(ts)%tree(ih)%hosthalo  ! set to the main host halo
           !
           if (host .gt. 0) then
              !
              if (univ(ts)%halo(host)%compute) then 
                !
                call halo_evolve(h,z,dt,post_merger,host=univ(ts)%halo(host)) 
                call halo_set_halo_evolved(univ(ts)%halo(ih)) 
                return    
             end if
           endif
        end if
! -------------------------------------------------                
#endif 
! SUB_TRANSFER 
        !
        ! in all other cases, the ejected gas is directly transfer to the IGM
        call halo_evolve(h,z,dt,post_merger) 
        ! 
        ! The halo has been evolved
        call halo_set_halo_evolved(univ(ts)%halo(ih)) 
        
        return
      end subroutine tree_evolve_halo
      
      ! *************************************************
      
      subroutine last_computable_progenitor_of_a_main_branch(ts,ih,last_ts,last_ih)

        ! RETURN "last_ts" AND "last_ih" THE TIME-STEP AND THE HALO INDEX OF THE LAST COMPUTABLE PROGENITOR 
        ! OF THE BRANCH WHICH CONTAINS THE HALO "ih" EVOLVING AT "ts" 

        implicit none

        integer(kind=4),intent(in)  :: ts       ! the current time-step
        integer(kind=4),intent(in)  :: ih       ! the current halo
        integer(kind=4),intent(out) :: last_ts  ! time-step of the last halo of the main branch
        integer(kind=4),intent(out) :: last_ih  ! index of the last computable halo of the main branch
        integer(kind=4)             :: iprog    ! the first progenitor index

        last_ts = ts                                       ! init
        last_ih = ih                                       ! init
        iprog   = univ(last_ts)%tree(last_ih)%firstProg    ! init
        
        do while (iprog .gt. 0)  
          !
          if (univ(last_ts-1)%halo(iprog)%compute) then
              ! 
              last_ts = last_ts -1
              last_ih = iprog
              iprog   = univ(last_ts)%tree(last_ih)%firstProg
          else
            !
            iprog = -1
          endif     
        end do
        
        return
      end subroutine last_computable_progenitor_of_a_main_branch

      ! *************************************************
      !end contains tree_evolve

  end subroutine tree_evolve
  
  !*****************************************************************************************************************

#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------  
  subroutine tree_spectrum
  
    ! RECEIVE HALO INFORMATION (HID, HID_host, level, aexp) AND CALL HIGHER SUBROUTINE DEDICATED TO SPECTRUM BUILER (halo_spectrum)
    
    implicit none

    integer(kind=4) :: ots                         ! index of the output timestep list
    integer(kind=4) :: current_ts                  ! the current timestap
    integer(kind=4) :: current_nb_computable_halo  ! the number of computable halo at the current timestep
    integer(kind=4) :: nb_gal_spectum              ! number of galaxy spectrum computed
         
    logical         :: go_down                     ! if go_down = .false. --> exit the loop 
    
    ! the local copy of the univers contains 1 halo at a given age of the univers
    if (.not. allocated(u%tree)) then
        !
        u%nb_of_halos = 1
        allocate(u%tree(u%nb_of_halos))
        allocate(u%halo(u%nb_of_halos))
    end if
    
    current_ts = -1
    ots = 1

    call IO_print_message('Tree[Spectrum]: Start',only_rank=rank,component='tree')
    do 
       !
       ! receive exit loop order, univ, halo and tree informations (HID, HID_host, level and aexp)
       call tree_receive_data(u,go_down)
       !
       if ((u%ts .gt. current_ts) .or. (.not. go_down)) then
          !
          ! Start a new time step or end of evolution
          if (current_ts .eq. tsout(ots)) then
              !
              ! The previous timestep was an output time step
              ! compare the number of galaxy spectrum build during the last time step with the number of computable halo
              ! associated to this timestep
              ! crash the code if these two values are diffrents
              if (nb_gal_spectum .ne. current_nb_computable_halo) then
                 !
                 call IO_print_error_message('nb_gal_spectrum /= current_nb_computable_halo',only_rank=rank,called_by='tree_spectrum')
                 call IO_print_message('used',only_rank=rank,component='tree', &
                          param_name=(/'nb_gal_spectum           ','nb_computable_halo       '/), &
                          int_param_val=(/nb_gal_spectum,current_nb_computable_halo/))  
                 stop
              end if
              ots = ots +1
          end if
          ! update current_ts
          current_ts = u%ts
          ! update current_nb_computable_halo
          current_nb_computable_halo = u%nb_of_computable_halos
          ! reset nb_gal_spectum
          nb_gal_spectum = 0
          !
          ! print the number of computable halo in this block
          if (go_down) then
             !
             write(IO_get_temporary_output_files_unit(current_ts)) current_nb_computable_halo
          end if    
       end if
       !
       ! test loop exit
       if (.not. go_down) then
          !
          call IO_print_message('Tree[Spectrum]: End',only_rank=rank,component='tree')
          exit
       else
          !
          ! Go to the halo scale
          ! At the halo scale data coming from tree and galaxy scale will be merged
          ! The stellar spectra, and associated data, will be save in the dedicated tmp file
          call halo_spectrum(u%halo(1))
          nb_gal_spectum = nb_gal_spectum +1
          !
          ! Save stellar spectra and tree informations in the dedicated tmp files
           call tree_print(IO_get_temporary_output_files_unit(current_ts),'tmp_bin','lum',u)
       end if
       !   
       ! deallocate the halo structure (galaxy)
       call halo_deallocate(u%halo(1))   ! the halo is deallocated
    end do
    
    return
  end subroutine tree_spectrum
! ------------------------------------------------- 
#endif
! LUMINOUS_PROCESSES  
  
  !*****************************************************************************************************************
  !
  ! FUNCTIONS
  !
  !*****************************************************************************************************************
  
  function tree_most_massive_galaxy_prog_index(ts,ih)
    
    ! SEND THE INDEX OF THE MOST MASSIVE PROGENITOR OF THE ih HALO AT THE TIME STEP ts
    ! In fact firstprog isn't always most massive prog
    
    implicit none
    
    integer(kind=4),intent(in) :: ts,ih
    integer(kind=4)            :: tree_most_massive_galaxy_prog_index
    integer(kind=4)            :: p
    
    real(kind=8)               :: m
    
    tree_most_massive_galaxy_prog_index = -1
    m                            = -1.0d0
    
    p = univ(ts)%tree(ih)%firstprog 
    do while (p .gt. 0)
       !
       if (univ(ts-1)%halo(p)%compute) then
          !
          if (galaxy_mass(univ(ts-1)%halo(p)%galaxy) .gt. m) then
             !
             m = galaxy_mass(univ(ts-1)%halo(p)%galaxy)
             tree_most_massive_galaxy_prog_index = p
          end if
       end if
       p = univ(ts-1)%tree(p)%nextprog
    end do
    
    return
  end function tree_most_massive_galaxy_prog_index
    
  !*****************************************************************************************************************
  !
  ! PRINTING PROCEDURES
  !
  !*****************************************************************************************************************

  subroutine tree_load_tree_data(u,ih,tree_data)

    ! CREATE THE TREE OUTPUT PROPERTIES LIST

    implicit none

    integer(kind=4),intent(in)     :: ih ! halo index
    
    real(kind=8), intent(inout)    :: tree_data(nb_tree_field)
    real(kind=8)                   :: halo_compute
 
    type(timestep_type),intent(in) :: u  ! a local copy of the universe (only with one halo)
    
    halo_compute = 0.d0
    if (u%halo(ih)%compute) halo_compute = 1.d0

    ! For information 
    !'z                     ','ts                    ','HID                   '
    !'host_HID              ','level                 ','n_dads                '
                                                                     
    tree_data = (/1.d0/u%aexp - 1.d0, real(u%ts,8), real(u%tree(ih)%HID,8), &
                  real(u%tree(ih)%host_HID,8),real(u%tree(ih)%level,8), real(u%tree(ih)%n_computed_dads,8)/)

    return
  end subroutine tree_load_tree_data

  !****************************************************************************************************************
    
  subroutine tree_print(unit,form,data_type,u,ih,go_down)
     
    ! PRINT TREE PROPERTIES IN THE TMP FILE unit 
      
    integer(kind=4),intent(in)          :: unit                       ! file unit
    integer(kind=4),intent(in),optional :: ih ! halo index
    integer(kind=4)                     :: ih_loc
    
    logical,optional                    :: go_down                    ! = .true. if higher level printing procedures have to be called
                                                                      ! only in 'fits' merger tree printing procedure and in tmp_bin mode
                                               
    character(*),intent(in)             :: form                       ! fits or tmp_bin 
    character(*),intent(in)             :: data_type                  ! phy or lum data

    real(kind=8)                        :: tree_data(nb_tree_field)   ! nb_tree_fields is give nin the header of this tree module
    real(kind=8)                        :: age
    
    type(timestep_type),intent(in)      :: u                          ! a local copy of the universe (only with one halo)
   
    select case (trim(data_type))
    case('phy','physical','phy_data')
        !
        ! PHYSICAL PROPERTIES
        !
        if (.not. present(ih)) then
            !
            call IO_print_error_message('No ih index in phy data print mode', &
                only_rank = rank, called_by = 'tree_print')
            stop ! stop the program
        end if
        !
        ! @ this point ih is present
        ih_loc = ih
        call tree_load_tree_data(u,ih,tree_data)
        ! select the output format
        select case (trim(form))
        case ('tmp_bin')
            !
            write(unit) tree_data ! print tree properties  
        case ('fits')
            !
            ! move to dm extension
            call ftmahd(unit,hdu_tree,hdutype,status) 
            !
            if (status .gt. 0) then
                !
                call IO_print_error_message('ftmahd status', &
                    only_rank = rank, called_by = 'tree_print')
                stop ! stop the program
            end if
            !
            ! init
            call ftirow(unit,0,1,status) 
            if (status .gt. 0) then
                !
                call IO_print_error_message('ftirow status', &
                    only_rank = rank, called_by = 'tree_print')
                stop ! stop the program
            end if
            !
            ! write data in the tree entension
            age = u%age_univ
            call ftpcld(unit,1,1,1,1,age,status)
            do i=2, nb_tree_field+1   
                !
                call ftpcld(unit,i,1,1,1,tree_data(i-1),status)  
                if (status .gt. 0) then
                    !
                    call IO_print_error_message('ftpcld status', &
                    only_rank = rank, called_by = 'tree_print')
                    stop ! stop the program
                end if
            end do
        case default
            !
            call IO_print_error_message('Unknwon output data format', &
                only_rank = rank, called_by = 'tree_print')
            stop ! stop the program
        end select
        !
    case('lum','luminous','lum_data')
        !
        ! LUMINOUS PROPERTIES
        !
        ih_loc = 1
        call tree_load_tree_data(u,ih_loc,tree_data)
        write(unit) tree_data ! print tree properties  
    case default
        !
        call IO_print_error_message('Unknwon data type', &
                only_rank = rank, called_by = 'tree_print')
        stop ! stop the program
    end select
    !
    ! call routine for print halo properties
    ! halo_print(unit,form,data_type,h)
    call halo_print(unit,form,data_type,u%halo(ih_loc),go_down=go_down)  
                  
    return
  end subroutine tree_print
      
  !****************************************************************************************************************
 
  recursive subroutine tree_print_tree(ts,ih,first_call)

    ! FOLLOW A HALO MERGER TREE AND PRINT TREE INFORMATIONS, BARYONIC AND DARK-MATTER PROPERTIES

    implicit none
      
    integer(kind=4),intent(in)  :: ts, ih ! timestep and halo index
    integer(kind=4)             :: p

    logical,intent(in),optional :: first_call
    
#ifdef PRINTALL
! -------------------------------------------------
    if ((nbproc .eq. 1) .and. first_call) call IO_print_message('tree_print_tree',component='tree')
! -------------------------------------------------
#endif

    ! tree_print(unit,form,data_type,u,ih=ih)
    call tree_print(merger_tree_unit(current_index),'fits','phy',univ(ts),ih=ih,go_down=.true.)
    if (univ(ts)%tree(ih)%n_computed_dads .gt. 0) then
        !
        p = univ(ts)%tree(ih)%firstprog
        do while (p .gt. 0) 
            !
            if (univ(ts-1)%halo(p)%compute) call tree_print_tree(ts-1,p,first_call=.false.)
            p = univ(ts-1)%tree(p)%nextprog 
        end do
    end if

    return
  end subroutine tree_print_tree
 
  !****************************************************************************************************************
  
end module tree
