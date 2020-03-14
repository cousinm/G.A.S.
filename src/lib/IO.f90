module IO
  
  use PrDi              ! Managment of MPI processes
  
  public
  
  !*****************************************************************************************************************
  ! 
  ! OVERVIEW
  !
  ! This module defined global parameters and procedures in link to input processes of merger-tree set files
  ! Message printing procedures are also defined in this module
  !
  ! SUBROUTINES IN THIS MODULE
  !
  !  IO_print_message(message,only_rank)            : Print a message 
  !    called by : all functions or subroutines which have to print something
  !   
  !  IO_print_error_message                         : Print a error message 
  !    called by : all functions or subroutines which have to print an error message
  !
  !  IO_print_warning_message                       : Print a error message (on screen)
  !    called by : all functions or subroutines which have to print a warning message
  ! 
  !  IO_read_parameter_file                         : Read the simualtion input parameter file. 
  !                                                      This parameter file allows to defined input/output path
  !    calley by : main program
  !    contains  : read_cosmology
  !                read_dmsim_param
  !                read_output_timestep_list
  !                read_redshift_list           
  !
  !  IO_distribute_tree_files                       : Allow to distribute all input merger-tree set files onto the differents MPI processes used by the program
  !    calley by : main program
  !    contains  : read_treefiles_index_list
  !                create_follow_up_list       
  !
  ! FUNCTION IN THIS MODULE
  !
  !  IO_get_tmp_filelist_unit                       : Return file unit for tmp filelist files 
  !
  !  IO_get_log_process_file_unit                   : Return file unit for log process files
  !
  !  IO_get_temporary_output_files_unit             : Return file unit for temporary output files
  !
  !  IO_return_ifile_from_file_index                : Return the ifile index associated to an input tree-fil index: 'nsteps.ifile'
  !
  !  IO_generate_HID                                : generate the halo identification number (HID) from 'ifile', 'ts' and 'ih'
  !
  !  IO_return_ifile_from_HID                       : return the 'ifile' index from the halo identification number
  !
  !  IO_return_ts_from_HID                          : return the last 'ts' from the halo identification number
  !
  !  IO_return_ih_from_HID                          : return the 'ih' index from the halo identification number 
  !
  ! FUNCTIONS IN THIS MODULE
  !
  !   No function in this module
  !
  ! PRINTING PROCEDURES
  !
  !   No printing procedure
  !
  !*****************************************************************************************************************
  !
  ! global parameters
  integer(kind=4)                      :: res                              ! resolution of the dark-matter simulation used
  integer(kind=4)                      :: nout                             ! number of ouput timestep (depend of the resolution of the simulation used)
  integer(kind=4)                      :: nts                              ! number of Nbody simulation output time-steps (depend of the resolution of the simulation used)
  integer(kind=4)                      :: nsteps                           ! number of Nbody simulation output time-steps (depend of the resolution of the simulation used)
                                                                           ! it is an overkill with nts (nts is used before have read nstep)
  integer(kind=4)                      :: ts_STOP                          ! STOP baryonic evolution time-step (given in the parameter file) 
  integer(kind=4),allocatable          :: tsout(:)                         ! list of output timestep (depend of the resolution of the simulation used)
  real(kind=4),allocatable             :: zout(:)                          ! list of output redshif (depend of the resolution of the simulation used)
  !
  ! Path and reference filename for input files
  character(MAXPATHSIZE)               :: input_path                       ! path for input library
  character(MAXPATHSIZE)               :: stellar_population_input_path    ! path to stellar population libraries
  character(MAXPATHSIZE)               :: tree_input_path                  ! path for input tree files
  character(MAXPATHSIZE)               :: output_path                      ! path for the output files
  character(MAXPATHSIZE)               :: timestep_filename                ! name of the timestep input file, given in input
  character(MAXPATHSIZE)               :: tree_struct_filename             ! name of the tree structure input file, given in input
  character(MAXPATHSIZE)               :: tree_props_filename              ! name of the tree properties input file, givne in input
  !
  ! units for input files and display protocole
  integer(kind=4),parameter            :: errunit                   = 6    ! stdout unit, for screen print message
  integer(kind=4),parameter            :: input_parameter_file_unit = 110  ! used to read input parameter file (IO_read_parameter_file)
  integer(kind=4),parameter            :: simu_param_file_unit      = 111  ! used to read dm-simulation parameter file
  integer(kind=4),parameter            :: cosmology_file_unit       = 112  ! used to read cosmology parameter file
  integer(kind=4),parameter            :: treefile_list_unit        = 113
  integer(kind=4),parameter            :: output_timestep_file_unit = 114
  integer(kind=4),parameter            :: redshift_list_unit        = 115  ! used to read redshift list (IO_read_parameter_file)
  integer(kind=4),parameter            :: follow_up_list_unit       = 116
  integer(kind=4),parameter            :: gasprop_unit              = 117  ! used to read stellar ejecta files (stellar_population_read_mass_loss_rates)
  integer(kind=4),parameter            :: massloss_unit             = 118  ! used to read stellar ejecta files (stellar_population_read_mass_loss_rates)
  integer(kind=4),parameter            :: snrate_unit               = 119  ! used to read SN event rates (stellar_population_read_SN_rates)
  integer(kind=4),parameter            :: stellarspectra_unit       = 120  ! used to read stellar spectra files (stellar_population_read_stellar_spectra)
  integer(kind=4),parameter            :: abs_sca_unit              = 121  ! used to read absorbtion, scattering dust properties (dust_read_abs_sca_properties)
  integer(kind=4),parameter            :: dust_sed_unit             = 122  ! used to read dust SED library
  integer(kind=4),parameter            :: cooling_unit              = 123  ! used to read cooling tables
  integer(kind=4),parameter            :: tree_file_unit            = 124  ! used to read merger-tree file extracted to the dark-matter simulation
  integer(kind=4),parameter            :: tstep_file_unit           = 125  ! used to read merger-tree file extracted to the dark-matter simulation
  integer(kind=4),parameter            :: props_file_unit           = 126  ! used to read merger-tree file extracted to the dark-matter simulation
  integer(kind=4),parameter            :: filter_list_unit          = 127  ! used to read the filelist of filter (filters_read_filters)
  integer(kind=4),parameter            :: filter_unit               = 128  ! used to read filter files
  integer(kind=4),parameter            :: agn_sed_unit              = 129  ! used to read agn SED library
  !
  ! unit for processes log
  integer(kind=4),parameter            :: log_process_unit          = 2000  
  !
  ! units for outputs files
  ! fits outputs parameter
  integer(kind=4),parameter            :: tmp_filelist_unit         = 3000 ! base of tmp_filelist units
  integer(kind=4),parameter            :: z_table_unit              = 211  ! used to open z table output file
  integer(kind=4),parameter            :: log_file_unit             = 212  ! unit used for the log file
  integer(kind=4)                      :: main_lum_FITS_file_unit        
  integer(kind=4)                      :: nb_main_output_files_to_print    ! number of main output files (phy or lum) that I have to print
  integer(kind=4),allocatable          :: list_of_ts_to_print(:)           ! list of timestep taht the process have to print
  integer(kind=4)                      :: nb_field                         ! number of field in the output FITS file, (its value depends of the TREE_EVOLVE paramerter)
  integer(kind=4),parameter            :: ttype_len                 = 22
  integer(kind=4),parameter            :: tunit_len                 = 12
  integer(kind=4),parameter            :: tform_len                 = 2
  !
  character(len=ttype_len),allocatable :: ttype(:)                         ! Name of each output colomn data
  character(len=tunit_len),allocatable :: tunit(:)                         ! Physical unit of each column data
  character(len=tform_len),allocatable :: tform(:)                         ! Data type of each column data
  integer(kind=4),parameter            :: varidat                   = 0    ! FITSIO library parameter. Fixed here definitively
  integer(kind=4)                      :: status                           ! if > 0 -> some problems in FITSIO library process                                          
  integer(kind=4)                      :: blocksize                        ! FITSIO parameter (output of routines)

  !
  ! follow_up and merger_tree file
  integer(kind=4)                      :: nb_follow_up_halos               ! number of follow_up halos
  integer(kind=4),target,allocatable   :: follow_up_unit(:)                ! unit used for follow_up file 
  integer(kind=4),target,allocatable   :: merger_tree_unit(:)              ! unit used for tree-file 
  integer(kind=4)                      :: current_index
  integer(kind=8),allocatable          :: list_of_follow_up_halo(:)        ! list of halo identification number associated to follow_up halos
  !
  ! keyword
  logical                              :: FOLLOW_UP, PR_FOLLOW_UP          ! =.true. and =.true. if halo is followed and must be printed.
  logical                              :: BARYON                           ! =.true. if baryon properties are computed and must be save.

  contains

  !*****************************************************************************************************************

  ! DISPLAY PROTOCOLE
  ! DEFINE SCREEN PRINTING PROCEDURE

  !*****************************************************************************************************************

  ! For standard message

  subroutine IO_print_message(message,only_rank,component,param_name,real_param_val,int_param_val)

    ! PRINT A MESSAGE (ON SCREEN)

    implicit none

    integer(kind=4),intent(in),optional :: only_rank        ! use for print a message for a given process
    integer(kind=4)                     :: i                ! loop over physical parameters that should be printed
    integer(kind=4)                     :: u
    integer(kind=4),intent(in),optional :: int_param_val(:) ! value of a physical parameters (integer format)

    character(*),intent(in)             :: message          ! the message  
    character(*),intent(in), optional   :: component        ! component (tree or halo or dm or galaxy or disc or bulge .... ) 
    character(30)                       :: space            ! indentation
    character(30)                       :: current_time     ! current_time
    character(25),intent(in),optional   :: param_name(:)    ! name of a physical parameters

    real(kind=8),intent(in),optional    :: real_param_val(:)  ! value of a physical parameters (real format)

    call cpu_time(current)  ! current time
    hour = floor((current-start)/60./60.)            ! number of hours
    mn   = floor(((current-start)-hour*60.*60.)/60.) ! number of min
    sec  = (current-start)-hour*60.*60-mn*60.        ! number of sec
    write(current_time,'(i2.2,a,i2.2,a,f5.2,a)') int(hour), 'h ', int(mn), 'mn ', sec, 'sec'
    
    if (present(component)) then
       select case (trim(component))
          case ('tree') 
             space = '|-'
          case ('halo') 
             space = '|--'
          case ('dm') 
             space = '|---'
          case ('baryon_halo','bh') 
             space = '|---'
          case ('galaxy','gal') 
             space = '|---'
          case ('disc') 
             space = '|----'
          case ('bulge') 
             space = '|----'
          case ('clumps') 
             space = '|-----'
          case ('agn') 
             space = '|-----'
          case ('gas') 
             space = '|------'
          case ('stars') 
             space = '|------'
          case ('dust') 
             space = '|------'
          case default 
             space = '|-'
       end select
    else
       space = '|-'
    end if
   
    if (present(only_rank)) then
       !    
#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------  
    if ((nbproc .gt. 2) .and. (.not. main_process)) then!
       u = IO_get_log_process_file_unit(rank)
    else
       u = errunit
    end if
! ------------               
#else
! ------------
    if ((nbproc .gt. 2) .and. (.not. main_process)) then
        u = IO_get_log_process_file_unit(rank)
    else
        u = errunit
    end if
! ------------------------------------------------- 
#endif
! LUMINOUS_PROCESSES     
       ! if only_rank is given in input of the subroutine
       if (rank .eq. only_rank) then 
          if (present(param_name)) then
            do i = 1, size(param_name)
              if (present(real_param_val)) then
                write(u,'(a,a,i3.3,a,a,a,a,a,a,a,es15.8)') trim(space), '-> <Pr ', rank, ' / ', trim(current_time), '> ', trim(message), ': ', param_name(i), ' = ', real_param_val(i)
              else
                 if (present(int_param_val)) then
                    write(u,'(a,a,i3.3,a,a,a,a,a,a,a,i4)') trim(space), '-> <Pr ', rank, ' / ', trim(current_time), '> ', trim(message), ': ', param_name(i), ' = ', int_param_val(i)
                 end if
              end if
            end do  
          else

             write(u,'(a,a,i3.3,a,a,a,a)') trim(space), '-> <Pr ', rank, ' / ', trim(current_time), '> ', trim(message)
          end if
       endif
    else
       if (main_process) then
          if (present(param_name)) then
            do i = 1, size(param_name) 
              if (present(real_param_val)) then
                write(errunit,'(a,a,a,a,a,a,es15.8)') trim(space), '-> ', trim(message), ': ', param_name(i), ' = ', real_param_val(i)
              else
                 if (present(int_param_val)) then
                    write(errunit,'(a,a,a,a,a,a,i4)') trim(space), '-> ', trim(message), ': ', param_name(i), ' = ', int_param_val(i)
                 else
                    write(errunit,'(a,a,a,a,a)') trim(space), '-> ', trim(message), ': ', param_name(i)
                 end if
              end if  
            end do
          else
             write(errunit,'(a,a,a)') trim(space),'-> ', trim(message) 
          end if
       endif
    end if
    !
    return
  end subroutine IO_print_message
  
  !*****************************************************************************************************************

  ! For error message

  subroutine IO_print_error_message(error_message,only_rank,called_by)  
      
    ! PRINT AN ERROR MESSAGE

    implicit none
    
    integer(kind=4),intent(in),optional :: only_rank     ! use for print an error message for a given process

    character(30)                       :: current_time  ! current_time
    character(*),intent(in)             :: error_message ! the error message
    character(*),intent(in),optional    :: called_by     ! name of the subroutine in which the error occurs 

    call cpu_time(current)  ! current time
    hour = floor((current-start)/60./60.)            ! number of hours
    mn   = floor(((current-start)-hour*60.*60.)/60.) ! number of min
    sec  = (current-start)-hour*60.*60-mn*60.        ! number of sec
    write(current_time,'(i2.2,a,i2.2,a,f5.2,a)') int(hour), 'h ', int(mn), 'mn ', sec, 'sec'
    
    if (present(only_rank)) then
      ! if only_rank is given in input of the subroutine
!~       if (rank .eq. only_rank) then
        if (present(called_by)) write(errunit, '(a,i3.3,a,a)') '|--> <Pr ', rank, '> !!! ERROR !!! | ', trim(called_by)
        write(errunit,'(a,i3.3,a,a,a,a)') ,'|--> <Pr ', rank, ' / ', trim(current_time), '> ', trim(error_message)
!~       endif
    else
      if (main_process) then
        if (present(called_by)) write(errunit, '(a,a)') '!!! ERROR in : ', trim(called_by)
        write(errunit,'(a,a)') '|--> ', trim(error_message)
      endif
    endif
    !         
    return
  end subroutine IO_print_error_message
   
  !*****************************************************************************************************************

  ! Forwarning message

  subroutine IO_print_warning_message(warning_message,only_rank,called_by)  
      
    ! PRINT AN WARNING MESSAGE

    implicit none
    
    integer(kind=4),intent(in),optional :: only_rank       ! use for print an warning message for a given process
    integer(kind=4)                     :: u

    character(30)                       :: current_time    ! current_time
    character(*),intent(in)             :: warning_message ! the error message
    character(*),intent(in),optional    :: called_by       ! name of the subroutine in which the error occurs 
    
    call cpu_time(current)  ! current time
    hour = floor((current-start)/60./60.)            ! number of hours
    mn   = floor(((current-start)-hour*60.*60.)/60.) ! number of min
    sec  = (current-start)-hour*60.*60-mn*60.        ! number of sec
    write(current_time,'(i2.2,a,i2.2,a,f5.2,a)') int(hour), 'h ', int(mn), 'mn ', sec, 'sec'

#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------  
    if (nbproc .gt. 3) then
        u = IO_get_log_process_file_unit(rank)
    else
        u = errunit
    end if
! ------------               
#else
! ------------
    if (nbproc .gt. 2) then
        u = IO_get_log_process_file_unit(rank)
    else
        u = errunit
    end if
! ------------------------------------------------- 
#endif
! LUMINOUS_PROCESSES 

    if (present(only_rank)) then
      ! if only_rank is given in input of the subroutine
      if (rank .eq. only_rank) then
        if (present(called_by)) write(errunit, '(a,i3.3,a,a)') '|--> <Pr ', rank, '> !!! WARNING !!! | ', trim(called_by)
        write(u,'(a,i3.3,a,a,a,a)') ,'|--> <Pr ', rank, ' / ', trim(current_time), '> ', trim(warning_message)
      endif
    else
      if (main_process) then
        if (present(called_by)) write(errunit, '(a,a)') '!!! WARNING !!! | ', trim(called_by)
        write(errunit,'(a,a)') '|--> ', trim(warning_message)
      endif
    endif
    !         
    return
  end subroutine IO_print_warning_message

  !*****************************************************************************************************************

  ! INPUT PROTOCOLE
  ! reding procedure for input libraries and files distribution through MPI processes

  !*****************************************************************************************************************

  subroutine IO_read_parameter_file
    
    ! READ SIMULATION USER PARAMETER FILE
    ! input / output path
    ! tree cleaning parameters
    ! physical parameters

    implicit none

    integer(kind = 4)       :: loop

    character(MAXPATHSIZE)  :: filename
    character(MAXPATHSIZE)  :: line
    character(MAXPATHSIZE)  :: name
    character(MAXPATHSIZE)  :: val

    !           
    ! read path acces to the parameters file
    call getarg(1,filename)
    !
    call IO_print_message('Read parameters file :')
    call IO_print_message(trim(filename))
    !
    close(input_parameter_file_unit)
    open(input_parameter_file_unit, file = filename, status = 'old')
    do
      read(input_parameter_file_unit, '(a)', end = 2) line
      loop = scan(line, '=')
      if ((loop .eq. 0) .or. (line(1:1) .eq. '#')) cycle ! '#' is the symbol for header comments in all input data files
      name = trim(adjustl(line(:loop-1)))
      val  = trim(adjustl(line(loop+1:)))
      ! '!' is the symbol for comments in all input data files
      loop = scan(val, '!')
      if (loop .ne. 0) val = trim(adjustl(val(:loop-1))) ! cut the comment part of the readed line
      !
      ! List of all keyword
      select case (trim(name)) 
      case ('input', 'input_path')        
        read(val, '(a)') input_path                 ! library input path input path
        ! set stellar_population_input_path
#ifdef TWOMYRS
! -------------------------------------------------
        write(stellar_population_input_path,'(a,a)') trim(input_path), '/2Myrs'
! ------------
#else
! ------------
        write(stellar_population_input_path,'(a,a)') trim(input_path), '/1Myr'
! -------------------------------------------------
#endif
! TWOMYRS     
        !
      case ('tree_input_path','tree_input','simu_input')
        read(val, '(a)') tree_input_path            ! tree files input path input path
        ! read cosmology
        call read_cosmology
        ! read simu_param.in
        call read_dmsim_param
      case ('output', 'output_path')               
        read(val, '(a)') output_path                ! output_path
      case ('nb_of_desc_halo_test')        
        read(val,*) nb_of_desc_halo_test            ! nb_of_desc_halo_test
      case ('min_dm_goodness_of_branch')   
        read(val,*) min_dm_goodness_of_branch       ! min_dm_goodness_of_branch
      case ('min_halo_goodness_of_branch') 
        read(val,*) min_halo_goodness_of_branch     ! min_halo_goodness_of_branch
      case ('Energy_threshold') 
        read(val,*) Ecut                            ! Energy_threshold
#ifdef TREE_EVOLVE
! -------------------------------------------------
      case ('z_reion')
        read(val,*) z_reion                         ! redshift of reionization 
        z_overlap = z_reion + 1.d0
      case ('physical_precision')          
        read(val,*) physical_precision              ! physical_precision
      case ('cold_stream_efficiency')          
        read(val,*) cold_stream_efficiency          ! take account cold stream or not
      case ('cooling_efficiency')          
        read(val,*) cooling_efficiency              ! take account cooling or not
      case ('TI_efficiency')          
        read(val,*) TI_efficiency                   ! propagation efficiency of Thermal instabilities
      case ('epsilon_merge')
        read(val,*) epsilon_merge                   ! critical mass ratio that separates major and minor mergers  
      case ('IMF')                                 
        read(val,*) IMF                             ! Initial mass function name (Salpeter, Kennicutt, Scalo)
      case ('disc_stripping_efficiency','stripping_efficiency','strip_efficiency')                
        read(val,*) disc_stripping_efficiency       ! gas stripping efficiency
      case ('disc_ejecta_efficiency','ejecta_efficiency','ej_efficiency')          
        read(val,*) disc_ejecta_efficiency          ! efficency of the ejecta process in the disc       
      case ('agn_ejecta_efficiency')
        read(val,*) agn_ejecta_efficiency           ! Mout/Macc
      case ('M_BH_min')
        read(val,*) M_BH_min                        ! minimal black hole mass (in M_sun) 
        ! conversion in code unit
        M_BH_min = M_BH_min*M_Sun_in_mass_code_unit 
      case ('agn_gas_coupling')                        
        read(val,*) agn_gas_coupling                ! gas coupling efficieny agn jet / gas
! -------------------------------------------------
#endif
! TREE_EVOLVE
      case default
        call IO_print_error_message('Input keyword not allowed',called_by = 'IO_read_parameter_file')  
        write(line,'(a,a,a)') 'Keyword: ', trim(name), ' unknown !'
        call IO_print_message(trim(line))
        stop ! stop the program
      end select
    end do  
    !
2   close(11)    
    !
    if (main_process) then  
      write(errunit,'(a)')         '| COSMOLOGY --------------------------------------------------------------' 
      write(errunit,'(a,f4.2)')    '|--> h                             : ', h_0
      write(errunit,'(a,f6.4)')    '|--> Omega_L                       : ', Omega_L
      write(errunit,'(a,f6.4)')    '|--> Omega_m                       : ', Omega_m
      write(errunit,'(a,f6.4)')    '|--> Omega_b                       : ', Omega_b
      write(errunit,'(a)')         '| DM-SIMULATION ----------------------------------------------------------' 
      write(errunit,'(a,i4.4,a)')  '|--> resolution                    : ', res, '^3'
      write(errunit,'(a,f5.1)')    '|--> L_box [Mpc/h]                 : ', L_box
      write(errunit,'(a)')         '| PATH -------------------------------------------------------------------'
      write(errunit,'(a,a)')       '|--> input_path                    : ', trim(input_path)
      write(errunit,'(a,a)')       '|--> tree_input_path               : ', trim(tree_input_path)
      write(errunit,'(a,a)')       '|--> output_path                   : ', trim(output_path) 
      write(errunit,'(a)')         '| DM CLEAN ---------------------------------------------------------------'
      write(errunit,'(a,i1.1)')    '|--> nb_of_desc_halo_test          : ', nb_of_desc_halo_test
      write(errunit,'(a,f4.2)')    '|--> min_dm_goodness_of_branch     : ', min_dm_goodness_of_branch
      write(errunit,'(a,f4.2)')    '|--> min_halo_goodness_of_branch   : ', min_halo_goodness_of_branch
      write(errunit,'(a,f4.2,a)')  '|--> Energy threshold              : ', Ecut
#ifdef TREE_EVOLVE
! -------------------------------------------------
      write(errunit,'(a)')         '| BACKGROUND -------------------------------------------------------------'
      write(errunit,'(a,e10.4)')   '|--> z_reion                       : ', z_reion
      write(errunit,'(a)')         '| ACCRETION --------------------------------------------------------------'
      write(errunit,'(a,e10.4)')   '|--> cold_stream_efficiency        : ', cold_stream_efficiency
      write(errunit,'(a,e10.4)')   '|--> cooling_efficiency            : ', cooling_efficiency
#ifdef THERMAL_INSTABILITIES 
! ------------------------------------------------- 
      if (TI_efficiency .gt. 0.d0)    write(errunit,'(a,e10.4)')   '|--> TI_efficiency                 : ', TI_efficiency
! ------------------------------------------------- 
#endif     
      write(errunit,'(a)')         '| GALAXY -----------------------------------------------------------------'
      if (epsilon_merge .gt. 0.d0)    write(errunit,'(a,e10.4)')   '|--> epsilon_merge                 : ', epsilon_merge
      write(errunit,'(a)')         '| DISC -------------------------------------------------------------------'
#ifdef SUB_STRIPPING
! ------------------------------------------------- 
      if (disc_stripping_efficiency .gt. 0.d0) write(errunit,'(a,e10.4)')   '|--> disc_stripping_efficiency     : ', disc_stripping_efficiency
! ------------------------------------------------- 
#endif         
      if (disc_ejecta_efficiency .gt. 0.d0)    write(errunit,'(a,e10.4)')   '|--> disc_ejecta_efficiency        : ', disc_ejecta_efficiency
      write(errunit,'(a)')         '| STARS ------------------------------------------------------------------'
      write(errunit,'(a,a)')       '|--> IMF                           : ', trim(IMF)
      write(errunit,'(a)')         '| AGN --------------------------------------------------------------------'
      if (agn_ejecta_efficiency .gt. 0.d0) write(errunit,'(a,e10.4)')   '|--> agn_ejecta_efficiency         : ', agn_ejecta_efficiency
      if (M_BH_min .gt. 0.d0)              write(errunit,'(a,e10.4)')   '|--> M_BH_min                      : ', M_BH_min
      if (agn_gas_coupling .gt. 0.d0)      write(errunit,'(a,e10.4)')   '|--> agn_gas_coupling              : ', agn_gas_coupling
! -------------------------------------------------
#endif
! TREE_EVOLVE
      write(errunit,'(a)')         '|------------------------------------------------------------------------- ' 
      !
    end if  
    !
    ! read the output timestep list
    call read_output_timestep_list
    !
    ! read redshift list
    call read_redshift_list
    
    return

    contains
    
   !**********************************************************
    
    subroutine read_cosmology
    
        ! READ COSMOLOGICAL PARAMETERS ASSOCIATED TO THE INPUT DM-SIMULATION
        
        implicit none
        
        integer(kind = 4)       :: loop

        character(MAXPATHSIZE)  :: filename
        character(MAXPATHSIZE)  :: line
        character(MAXPATHSIZE)  :: name
        character(MAXPATHSIZE)  :: val
        
        ! build the filename
        call IO_print_message('Load data from : cosmology.in')
        write(filename,'(a,a,a)') trim(tree_input_path), '/cosmology.in'
        open(unit = cosmology_file_unit, file = trim(filename), status = 'old')
    
        do
            read(cosmology_file_unit, '(a)', end = 2) line
            loop = scan(line, '=')
            if ((loop .eq. 0) .or. (line(1:1) .eq. '#')) cycle ! '#' is the symbol for header comments in all input data files
            name = trim(adjustl(line(:loop-1)))
            val  = trim(adjustl(line(loop+1:)))
            loop = scan(val, '!')  ! '!' is the symbol for comments in all input data files
            if (loop .ne. 0) val = trim(adjustl(val(:loop-1))) ! cut the comment part of the readed line
            !
            ! List of keyword
            select case (trim(name)) 
            case ('h','h0','hubble')           
                read(val,*) h_0                   ! reduced hubble parameter
            case ('Omega_L','OmegaL','Omega_lambda')
                read(val,*) Omega_L               ! Dark energy
            case ('Omega_0','Omega_m','Omegam','Omega_matter')
                read(val,*) Omega_m               ! Matter
            case ('Omega_b','Omegab','Omega_baryon','Omega_bar')
                read(val,*) Omega_b               ! Baryon
                baryon_fraction = Omega_b/Omega_m ! compute baryonic fraction   
            case ('n','slope')
            
            case ('sig8','sigma8')  
            
            case default
                call IO_print_error_message('Input keyword not allowed',called_by = 'IO_read_dmsim_param')  
                write(line,'(a,a,a)') 'Keyword: ', trim(name), ' unknown !'
                call IO_print_message(trim(line))
                stop ! stop the program
            end select
        end do  
        !
2       close(cosmology_file_unit)
        
        h_0_code_unit      = h_0/Mpc_in_km*Gyr_in_s                                      ! from km/s/Mpc --> Gyr  
        rho_crit_code_unit = 3.d0*(1.d2*h_0_code_unit)**2./(8.d0*pi*gravconst_code_unit) ! in 10^11Msun / kpc^3 (~)
        
        return
    end subroutine read_cosmology
    
    !**********************************************************
    
    subroutine read_dmsim_param
    
        ! READ PARAMETERS OF THE INPUT DM-SIMULATION
        
        implicit none
        
        integer(kind = 4)       :: loop

        character(MAXPATHSIZE)  :: filename
        character(MAXPATHSIZE)  :: line
        character(MAXPATHSIZE)  :: name
        character(MAXPATHSIZE)  :: val
        
        ! build the filename
        call IO_print_message('Load data from : simu_param.in')
        write(filename,'(a,a,a)') trim(tree_input_path), '/simu_param.in'
        open(unit = simu_param_file_unit, file = trim(filename), status = 'old')
    
        do
            read(simu_param_file_unit, '(a)', end = 2) line
            loop = scan(line, '=')
            if ((loop .eq. 0) .or. (line(1:1) .eq. '#')) cycle ! '#' is the symbol for header comments in all input data files
            name = trim(adjustl(line(:loop-1)))
            val  = trim(adjustl(line(loop+1:)))
            loop = scan(val, '!')  ! '!' is the symbol for comments in all input data files
            if (loop .ne. 0) val = trim(adjustl(val(:loop-1))) ! cut the comment part of the readed line
            !
            ! List of keyword
            select case (trim(name)) 
            case ('L_box','Lbox','L')           
                read(val,*) L_box                   ! box size [Mpc/h]
            case ('N_particles','Nparticles','Npart','N')        
                read(val,*) res                ! resolution (number of particule in the simulation)
            case ('dm_particle_mass','dm_mass','dm_part_mass')
                read(val,*) dm_particle_mass   ! mass of one dark-matter particle [Msun/h]
                dm_particle_mass = dm_particle_mass * M_Sun_in_mass_code_unit ! convert in code unit
            case default
                call IO_print_error_message('Input keyword not allowed',called_by = 'IO_read_dmsim_param')  
                write(line,'(a,a,a)') 'Keyword: ', trim(name), ' unknown !'
                call IO_print_message(trim(line))
                stop ! stop the program
            end select
        end do  
        !
2       close(simu_param_file_unit)
        
        return
    end subroutine read_dmsim_param
    
    !**********************************************************

    subroutine read_output_timestep_list

      ! READ THE OUTPUT TIMESTEP LIST 
      ! correspondig to a selection in the timestep list of the dark-matter simulation used (depends of the resolution)

      implicit none

      integer(kind=4)            :: ots        ! timestep

      character(MAXPATHSIZE)     :: filename
      character(MAXPATHSIZE)     :: message    ! a message  
      character(MAXPATHSIZE)     :: line       ! a lilne read in the file

      ! build the filename
      call IO_print_message('Load data from : output_timestep_list.in')
      write(filename,'(a,a,a)') trim(tree_input_path), '/output_timestep_list.in'
      open(unit = output_timestep_file_unit, file = trim(filename), status = "old")
      
      do
          read(output_timestep_file_unit, '(a)', end = 2) line
          if (trim(line) .eq. 'START') then
             ! all lines before the START keywork are header lines: they contain informations about data 
             !
             read(output_timestep_file_unit,*) nout
             !
             if (nout .gt. 0) then
                ! create tsout(:)
                allocate(tsout(nout))
                !
                do ots = 1, nout
                    read(output_timestep_file_unit,*) tsout(ots)
                end do
                ! set ts_STOP and nsteps
                ts_STOP  = maxval(tsout)
                !
                if (main_process) then
                    ! print some information
                    if (nout .le. 8) then
                        call IO_print_message('Data will be saved for timesteps:') 
                        write(message,'(i3.3)') tsout(1)
                        do ots = 2, nout
                            write(message,'(a,a,(2x,i3.3))') trim(message), ',', tsout(ots)
                        end do
                    else
                        write(message,'(a,i3.3,a)') 'Data will be saved for ', nout, ' timesteps'
                    end if
                    call IO_print_message(message)      
                    if(ts_STOP .gt. 0) then
                        write(message,'(a,i3.3)') 'ts_STOP : ', ts_STOP
                        call IO_print_message(message)   
                    end if
                end if 
             end if
             exit  ! quit do loop
          end if
          if (line(1:1) .eq. '#') then
             cycle ! header or something like this (skip)
          else
             call IO_print_error_message('Impossible to read the line',only_rank=rank,called_by='IO_read_output_timestep_list') 
             stop ! stop the program
          end if
      end do
2     close(output_timestep_file_unit)

    end subroutine read_output_timestep_list

    !**********************************************************

    subroutine read_redshift_list

      ! READ THE REDSHIFT LIST 
      ! correspondig to the timestep list of the dark-matter simulation used (depends of the resolution)

      implicit none

      integer(kind=4)            :: ts         ! timestep index loop
      integer(kind=4)            :: its        ! timestep
      integer(kind=4)            :: ots        ! output timestep index loop

      character(MAXPATHSIZE)     :: filename
      character(MAXPATHSIZE)     :: message    ! a message  
      character(MAXPATHSIZE)     :: line       ! a lilne read in the file
      
      real(kind=4)               :: z          ! redshift value

      ! build the filename
      call IO_print_message('Load data from : redshift_list.in')
      write(filename,'(a,a)') trim(tree_input_path), '/redshift_list.in' 
      open(unit = redshift_list_unit, file = trim(filename), status = "old")

      do
          read(redshift_list_unit, '(a)', end = 2) line
          if (trim(line) .eq. 'START') then
             ! all lines before the START keywork are header lines: they contain informations about data 
             !
             read(redshift_list_unit,*) nts
             ! set nsteps
#ifdef UNLINKED_TREES
! -------------------------------------------------
             nsteps = nts
! ------------
#else
! ------------
             nsteps = min(ts_STOP + nb_of_desc_halo_test, nts)
! -------------------------------------------------
#endif
! UNLINKED_TREES
             !
             if (nout .gt. 0) then
                ! create zout(:)
                allocate(zout(nout))
                !
                ots = 1
                do ts = 1, ts_STOP
                    read(redshift_list_unit,*) its, z
                    if (its .eq. tsout(ots)) then
                        zout(ots) = z
                        ots = ots + 1
                    end if
                end do
                if (main_process) then
                    ! print some informations
                    if (nout .le. 8) then
                        call IO_print_message('Output timesteps corresponds to redshifts:') 
                        write(message,'(f5.3)') zout(1)
                        do ots = 2, nout
                            write(message,'(a,a,(2x,f5.3))') trim(message), ',', zout(ots)
                        end do
                        call IO_print_message(message) 
                    end if
                end if
            end if
            exit  ! quit do loop
          end if
          if (line(1:1) .eq. '#') then
             cycle ! header or something like this (skip)
          else
             call IO_print_error_message('Impossible to read the line',only_rank=rank,called_by='IO_read_redshift_list') 
             stop ! stop the program
          end if
      end do
2     close(redshift_list_unit)

    end subroutine read_redshift_list

    !**********************************************************
    ! end contains
  
  end subroutine IO_read_parameter_file

  !*****************************************************************************************************************

  subroutine IO_distribute_tree_files

    ! DISTRIBUTE ALL TREE-FILE BETWEEN ALL MPI PROCESSES

    implicit none

    integer(kind=4)                    :: ifile
    integer(kind=4)                    :: n                                    ! Index for loop throught MPI process list
    integer(kind=4)                    :: ts                                   ! timestep index loop for read data in file
    integer(kind=8)                    :: max_nb_halos_in_a_tree_file          ! Maximal number of halos in a input tree-file  
    integer(kind=4)                    :: nb_files_takes_into_account          ! Number of tree-files take into account into the MPI distribution
    integer(kind=4)                    :: ifile_of_max                         ! index of the (not already computed) tree-file that contains the largest number of halos
    integer(kind=4)                    :: mloc
    integer(kind=8)                    :: all_halos
    integer(kind=8),allocatable        :: nb_halos_in_each_tree_file(:)        ! Total number of halos in a given tree-file
    integer(kind=4),allocatable        :: tmp(:)                               ! A temporary arry which save the number of halos identified at a given time-step
    integer(kind=8),allocatable        :: nb_halos_computed_per_proc(:)        ! Number of haloes computed per each MPI process
    integer(kind=4),allocatable        :: nb_tree_files_computed_per_proc(:)   ! Number of tree-files read by each MPI process
                                                                               ! ONLY FOR RANK 0                                                                     
    character(len=10),allocatable      :: list_of_computed_files_per_proc(:,:) ! Contains indexes of tree-file computed by each MPI procces
                                                                               ! ONLY FOR RANK 0  
    character(MAXPATHSIZE)             :: stepsfile                            ! filename
    character(MAXPATHSIZE)             :: format

    ! This routine distribute all the input tree-files into the diffrents processes used
    !      
    !
    ! read the tree-file index list
    call read_treefiles_index_list(list_of_computed_files)
    !  
    ! create compute file list for different processus 
    if (main_process) then ! The Main process compute the number of haloes in the different input data files and separate them into the differents processes
                           ! All physical process must compute (in mean) the same number of haloes.

      write(errunit,'(a,i2.2,a)') '|--> Distribute tree files onto ',nb_physical_processes,' (physical) processes' 

      allocate(nb_halos_in_each_tree_file(nb_tree_files_computed))                      ! Contains the number of haloes (sum over all time-step) 
                                                                                        ! for each tree-file
      do ifile = 1, nb_tree_files_computed
        ! Open the file where properties of the trees at each timestep is saved
        ! Build the path throught the stepfile: ifile
        write(stepsfile,'(a,a,a,a)') trim(tree_input_path), '/tstep/', trim(timestep_filename), trim(list_of_computed_files(ifile))  
        open(unit = tree_file_unit, file = stepsfile, status = "old", form = "unformatted")
        read(tree_file_unit) ! nsteps                                                   ! Skip the number of timestep, already read previously
        if (ifile .eq. 1) then
           allocate(tmp(nsteps))                                                        ! In the first loop allocate a temporary array which contains,  
                                                                                        ! the number of halos at a given time-stpep in a give tree-file: ifile
        end if
        read(tree_file_unit) (tmp(ts),ts=1,nsteps)                                      ! Read the number of halos at each timestep for the file: ifile
        nb_halos_in_each_tree_file(ifile) = sum(tmp(:))                                 ! Total number of halos (sun over all ts) in the file: ifile
        close(tree_file_unit)                                                           ! Close the tree-file: ifile
      end do
      
      ! deallocate the temporary array
      deallocate(tmp)
      !
      ! Print some stats
      all_halos = sum(nb_halos_in_each_tree_file,dim=1)
      max_nb_halos_in_a_tree_file = maxval(nb_halos_in_each_tree_file,dim=1) ! the maximal number of halos
      mloc = maxloc(nb_halos_in_each_tree_file, dim=1)
      write(errunit,'(a,i9.9,a,i3.3,a)') '|---> Simulation contains : ', all_halos, ' halos in ', nb_tree_files_computed, ' differents input tree-files'
      write(errunit,'(a,i3.3,a)') '|---> All merger-trees evolve onto : ', nsteps, ' main time-steps '
      write(errunit, '(a,i3.3,a,i8.8,a)') '|---> The largest number of halos is contained in the tree-file : ', mloc, ' which contain : ', max_nb_halos_in_a_tree_file , ' halos'  
      !
      ! Distribution throught the diferent processus 
      ! the processus 0 doesn't compute halo evolution
      ! the biggest tree-file is computed by process 1)
      ! allocate some arrays
      allocate(list_of_computed_files_per_proc(nb_tree_files_computed,nb_physical_processes))  
      list_of_computed_files_per_proc(:,:) = ''   ! Init the list         
      ! 
      allocate(nb_halos_computed_per_proc(nb_physical_processes))       ! Number of halos computed by each MPI process
      nb_halos_computed_per_proc(:)        = 0                          ! Init the table
      ! 
      allocate(nb_tree_files_computed_per_proc(nb_physical_processes))  ! Number of files read by each MPI process
      nb_tree_files_computed_per_proc(:)   = 0                          ! Init the table
      !
      nb_files_takes_into_account          = 0                          ! Init the number of file taking account by the repartition algorithm
      !
      ! First run, build the list of tree-file computed by each MPI process  
      do n = 1, nb_physical_processes
        do while ((nb_halos_computed_per_proc(n) + maxval(nb_halos_in_each_tree_file)) .le. max_nb_halos_in_a_tree_file)
          if (nb_files_takes_into_account .eq. nb_tree_files_computed) exit 
          ! Search the index of the tree-file which contains the largest number of halos
          ifile_of_max                                                   = maxloc(nb_halos_in_each_tree_file,dim=1) 
          nb_halos_computed_per_proc(n)                                  = nb_halos_computed_per_proc(n) + nb_halos_in_each_tree_file(ifile_of_max)   
          nb_halos_in_each_tree_file(ifile_of_max)                       = 0  ! This tree-file is take into account, reset its number of halos
          ! The proc "n" must compute a other file
          nb_tree_files_computed_per_proc(n)                             = nb_tree_files_computed_per_proc(n) + 1 
          ! Add the index of the tree-file which contains the largest number of halo into the list of tree-files computed by the proc "n"
          list_of_computed_files_per_proc(nb_tree_files_computed_per_proc(n),n) = list_of_computed_files(ifile_of_max)  ! add the index of this file 
          ! A other tree-file has been take into account in the distribution process
          nb_files_takes_into_account                                    = nb_files_takes_into_account + 1        
          if (nb_files_takes_into_account .eq. nb_tree_files_computed) exit 
        end do
      end do
      !
      if (nb_files_takes_into_account .ne. nb_tree_files_computed) then
        do while (nb_files_takes_into_account .lt. nb_tree_files_computed)
          n = minloc(nb_halos_computed_per_proc, dim=1, mask=(nb_halos_computed_per_proc .gt. 0))  
          ! Search the index of the tree-file which contains the largest number of halos
          ifile_of_max                                                          = maxloc(nb_halos_in_each_tree_file,dim=1)
          nb_halos_computed_per_proc(n)                                         = nb_halos_computed_per_proc(n) + nb_halos_in_each_tree_file(ifile_of_max)
          nb_halos_in_each_tree_file(ifile_of_max)                              = 0   ! This tree-file is take into account, reset its number of halos
          nb_tree_files_computed_per_proc(n)                                    = nb_tree_files_computed_per_proc(n) + 1
          ! Add the index of the tree-file which contains the largest number of halo into the list of tree-files computed by the proc "n"
          list_of_computed_files_per_proc(nb_tree_files_computed_per_proc(n),n) = list_of_computed_files(ifile_of_max)  ! add the index of this file
          ! A other tree-file has been take into account in the distribution process
          nb_files_takes_into_account                                           = nb_files_takes_into_account + 1
        end do
      endif  
      !
      ! Second run (send and reveive)
      ! The process of rank 0 which know the global list, send to each other process its list of tree-file 
      do n = 1, nb_physical_processes 
        !
        ! The process "list_of_physical_processes(n)" must compute some halos
        ! Send the number of tree-file that the process "list_of_physical_processes(n)" must take into account
        call MPI_SEND(nb_tree_files_computed_per_proc(n),1,MPI_INTEGER4,list_of_physical_processes(n),io_tag+1,MPI_COMM_WORLD,ierror)
        !
#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------            
        ! Send the number of tree-file that the luminous process associated to the physical process "list_of_physical_processes(n)" must take into account
        call MPI_SEND(nb_tree_files_computed_per_proc(n),1,MPI_INTEGER4,list_of_luminous_processes(n),io_tag+1,MPI_COMM_WORLD,ierror)
! ------------------------------------------------- 
#endif
! LUMINOUS_PROCESSES         
        !
        if (nb_halos_computed_per_proc(n) .gt. 0) then
          ! print some informations
          if (nb_tree_files_computed_per_proc(n) .le. 10) then 
            ! print all reference numbers
            write(format,'(a,i3.3,a)') '(a,i3.3,a,', nb_tree_files_computed_per_proc(n), '(a,2x))'
            write(errunit,format) '|---> MPI process ', list_of_physical_processes(n), ' will compute tree-files: ', &
                list_of_computed_files_per_proc(1:nb_tree_files_computed_per_proc(n),n)
          else
            ! print only the two first and the two last reference numbers
            write(errunit,'(a,i3.3,a,i3.3,a,2(a,2x),a,2(a,2x),a)') '|---> MPI process ', list_of_physical_processes(n), &
                  ' will compute ', nb_tree_files_computed_per_proc(n), ' tree-files: (  ', list_of_computed_files_per_proc(1:2,n), '...', &
                  list_of_computed_files_per_proc(nb_tree_files_computed_per_proc(n)-1:nb_tree_files_computed_per_proc(n),n), ')'
          end if
          !  
          ! I am the process of rank 0 and I have to send the list into the other MPI processes
          ! send to the process "n" the list of tree-files that it must compute    
          call MPI_SEND(nb_halos_computed_per_proc(n),1,MPI_INTEGER8,list_of_physical_processes(n),io_tag+2,MPI_COMM_WORLD,ierror)           
          !
          ! Send the explicit list off tree-files that the process "list_of_physical_processes(n)" must take into account
          call MPI_SEND(list_of_computed_files_per_proc(1:nb_tree_files_computed_per_proc(n),n), &
            nb_tree_files_computed_per_proc(n)*10,MPI_CHARACTER,list_of_physical_processes(n),io_tag+3,MPI_COMM_WORLD,ierror)
        end if
      end do 
      ! 
      ! deallocate array used to distributed files throught the differents MPI processes
      if (allocated(nb_halos_in_each_tree_file))      deallocate(nb_halos_in_each_tree_file)   ! Total number of halos in a given tree-file 
      if (allocated(tmp))                             deallocate(tmp)                          ! A temporary arry which save the number of halos identified at a given time-step
      if (allocated(nb_halos_computed_per_proc))      deallocate(nb_halos_computed_per_proc)   ! Number of haloes computed per each MPI process
      !         
    else  ! end if rank = main_process_rank
      ! I am not the main process
      nb_files_computed = 0 ! init
      if (physical_process) then
        ! I am a physical process
        call MPI_RECV(nb_files_computed,1,MPI_INTEGER4,main_process_rank,io_tag+1,MPI_COMM_WORLD,statut,ierror)  ! Receive the number of tree-files
        if (nb_files_computed .gt. 0) then        
            call MPI_RECV(nb_halos_computed,1,MPI_INTEGER8,main_process_rank,io_tag+2,MPI_COMM_WORLD,statut,ierror)      ! Receive the number of halos 
            allocate(list_of_computed_files(nb_files_computed))                                    ! Allocate the array which contain the explicit list of tree-files indexes
            ! Receive the explicit list of tree-files 
            call MPI_RECV(list_of_computed_files,nb_files_computed*10,MPI_CHARACTER,main_process_rank,io_tag+3,MPI_COMM_WORLD,statut,ierror) 
        end if
      endif
#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------        
      if (luminous_process) then
        ! I am a luminous process
        ! receive the number of tree-file computed by my associated physical process
        ! this info is usefull for compute the number of blocks in tmp output files
        call MPI_RECV(nb_files_computed,1,MPI_INTEGER4,main_process_rank,io_tag+1,MPI_COMM_WORLD,statut,ierror) ! Receive the number of tree-files
      end if 
! ------------------------------------------------- 
#endif
! LUMINOUS_PROCESSES        
    endif ! end if rank = 0
    
    call create_follow_up_list(nb_tree_files_computed_per_proc,list_of_computed_files_per_proc)
    
    if (allocated(list_of_computed_files_per_proc))  deallocate(list_of_computed_files_per_proc)   ! Contains index of tree-file computed by each MPI procces 
    if (allocated(nb_tree_files_computed_per_proc))  deallocate(nb_tree_files_computed_per_proc)   ! Number of tree-files read by each MPI process   
       
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    !
    return

    contains
  
    !**********************************************************
    
    subroutine read_treefiles_index_list(list_of_computed_files)
    
        ! READ THE LIST OF INDEXES ASSOCIATED TO TREE FILES USED IN INPUT OF THE CURRENT RUN 
         
        implicit none
        
        integer(kind=4)                             :: ifile      ! loop index 
        
        character(len=10),allocatable,intent(inout) :: list_of_computed_files(:)
        character(MAXPATHSIZE)                      :: filename
        character(MAXPATHSIZE)                      :: line       ! a lilne read in the file

        ! build the filename
        call IO_print_message('Load data from : treefile_list.in')
        write(filename,'(a,a,a)') trim(tree_input_path), '/treefile_list.in'
        open(unit = treefile_list_unit, file = trim(filename), status = "old")
      
        do
          read(treefile_list_unit, '(a)', end = 2) line
          if (trim(line) .eq. 'START') then
             ! all lines before the START keywork are header lines: they contain informations about data 
             !
             read(treefile_list_unit,*) timestep_filename          ! name of the timestep file
             read(treefile_list_unit,*) tree_struct_filename       ! name of the tree structure file
             read(treefile_list_unit,*) tree_props_filename        ! name of the tree properties file
             if (main_process) then
                read(treefile_list_unit,*) nb_tree_files_computed     ! number of tree files
                !
                if (nb_tree_files_computed .gt. 0) then
                    ! create list_of_computed_files(:)
                    allocate(list_of_computed_files(nb_tree_files_computed))
                    !
                    do ifile  = 1, nb_tree_files_computed
                        read(treefile_list_unit,*) list_of_computed_files(ifile)
                    end do
                end if
             end if
             exit  ! quit do loop
          end if
          if (line(1:1) .eq. '#') then
             cycle ! header or something like this (skip)
          else
             call IO_print_error_message('Impossible to read the line',only_rank=rank,called_by='IO_read_treefiles_index_list') 
             stop ! stop the program
          end if
        end do
2       close(treefile_list_unit)
    
        return
    end subroutine read_treefiles_index_list
  
    ! *************************************************
  
    subroutine create_follow_up_list(nb_tree_files_computed_per_proc,list_of_computed_files_per_proc)
  
        ! DISTRIBUTE "FOLLOW_UP" HALOS ONTO THE DIFFERENT PROCESSES ACCORDING TO THEIR HALO IDENTIFICATION NUMBER
    
        implicit none

        integer(kind=4)                          :: nb_halos_followed                    ! the number of "follow_up" halos 
        integer(kind=4)                          :: ih,n,ifile                           ! loop indexes
        integer(kind=4),intent(in),allocatable   :: nb_tree_files_computed_per_proc(:)   ! Number of tree-files read by each MPI process
                                                                                       ! ONLY FOR RANK 0
        character(len=10),intent(in),allocatable :: list_of_computed_files_per_proc(:,:) ! Contains index of tree-file computed by each MPI procces
                                                                                       ! ONLY FOR RANK 0    
        integer(kind=8),allocatable              :: list_of_follow_up_files(:)
        integer(kind=8),allocatable              :: list_of_follow_up_files_per_proc(:,:)
        integer(kind=4),allocatable              :: nb_of_follow_up_files_per_proc(:)
    
        character(MAXPATHSIZE)                   :: filename           ! the filename of the halo identification number list
        character(MAXPATHSIZE)                   :: message            ! a message  
        character(MAXPATHSIZE)                   :: line               ! a line read in the file
        character(MAXPATHSIZE)                   :: format
        
        if (main_process) then
            !
            call IO_print_message('IO_create_follow_up_list')
        
            ! Build the complete filename
            call IO_print_message('Load data from : HID_list.in')
            write(filename,'(a,a)') trim(tree_input_path), '/HID_list.in'
    
            ! Read the follow_up_halo list
            open(unit = follow_up_list_unit, file = filename, status = "old")

            do
                read(follow_up_list_unit, '(a)', end = 2) line
                if (trim(line) .eq. 'START') then
                    !
                    ! Read the number of halo identification number given in the list
                    read(follow_up_list_unit,*) nb_halos_followed   ! read the number of "follow_up" halos 
                    !
                    ! allocate nb_of_follow_up_files_per_proc
                    allocate(nb_of_follow_up_files_per_proc(nb_physical_processes))
                    nb_of_follow_up_files_per_proc(:) = 0
                    !
                    if (nb_halos_followed .gt. 0) then
                        !   
                        write(message,'(a,i3.3,a)') ' The current run will generated ',  nb_halos_followed, ' time-lines and merger-tree files'
                        call IO_print_message(message)
                        !
                        ! allocate list_of_follow_up_files_per_proc
                        allocate(list_of_follow_up_files_per_proc(nb_halos_followed,nb_physical_processes))
                        !
                        ! allocate list_of_follow_up_files
                        allocate(list_of_follow_up_files(nb_halos_followed))
                        !
                        do ih = 1, nb_halos_followed
                            read(follow_up_list_unit,*) list_of_follow_up_files(ih)
                        end do  
                        do n = 1, nb_physical_processes
                            do ifile = 1, nb_tree_files_computed_per_proc(n)
                                do ih = 1, nb_halos_followed
                                    if (IO_return_ifile_from_HID(list_of_follow_up_files(ih)) .eq. &
                                                IO_return_ifile_from_file_index(list_of_computed_files_per_proc(ifile,n))) then
                                        nb_of_follow_up_files_per_proc(n) = nb_of_follow_up_files_per_proc(n) +1
                                        list_of_follow_up_files_per_proc(nb_of_follow_up_files_per_proc(n),n) = list_of_follow_up_files(ih)
                                    end if
                                end do
                            end do
                        end do
                    else
                        call IO_print_warning_message('The "follow up" list is empty',called_by='IO_create_follow_up_list') 
                    end if
                    !
                    exit  ! quit do loop
                end if
                if (line(1:1) .eq. '#') then
                    cycle ! header or something like this (skip)
                else
                    call IO_print_error_message('Impossible to read the line',only_rank=rank,called_by='IO_create_follow_up_list') 
                    stop ! stop the program
                end if  
            end do
2           close(follow_up_list_unit)
            !
            do n = 1, nb_physical_processes
                !
                ! I am the main process and I have to send the list into the other MPI processes
                ! send to the process "list_of_physical_processes(n)" the list of halos identification number that it must follow
                call MPI_SEND(nb_of_follow_up_files_per_proc(n),1,MPI_INTEGER4,list_of_physical_processes(n),io_tag,MPI_COMM_WORLD,ierror)
                if (nb_of_follow_up_files_per_proc(n) .gt. 0) then
                    ! 
                    ! print some informations
                    if (nb_of_follow_up_files_per_proc(n) .gt. 5) then
                        write(errunit,'(a,i3.3,a,i3.3,a,2(x,i11.11),a,2(x,i11.11),a)') '|---> MPI process ', list_of_physical_processes(n), ' will follow ', &
                                nb_of_follow_up_files_per_proc(n), ' halos: ( ', &
                                list_of_follow_up_files_per_proc(1:2,n), ' ... ', &
                                list_of_follow_up_files_per_proc(nb_of_follow_up_files_per_proc(n)-1:nb_of_follow_up_files_per_proc(n),n), ' )'
                    else
                        write(format,'(a,i3.3,a)') '(a,i3.3,a,', nb_of_follow_up_files_per_proc(n), '(i11.11,2x))'
                        write(errunit,format) '|---> MPI process ', list_of_physical_processes(n), ' will follow halos: ', &
                                list_of_follow_up_files_per_proc(1:nb_of_follow_up_files_per_proc(n),n)
                    end if
                    ! The process "list_of_physical_processes(n)" must follow some halos
                    ! Send the list
                    call MPI_SEND(list_of_follow_up_files_per_proc(1:nb_of_follow_up_files_per_proc(n),n), &
                            nb_of_follow_up_files_per_proc(n),MPI_INTEGER8,list_of_physical_processes(n),io_tag,MPI_COMM_WORLD,ierror)
                end if
            end do
        end if
        !
        if (physical_process) then
            ! I am a physical process
            call MPI_RECV(nb_follow_up_halos,1,MPI_INTEGER4,main_process_rank,io_tag,MPI_COMM_WORLD,statut,ierror)                          ! Receive the follow_up halos
            if (nb_follow_up_halos .gt. 0) then        
                allocate(list_of_follow_up_halo(nb_follow_up_halos))                                                                        ! Allocate the array which contain the explicit list
                call MPI_RECV(list_of_follow_up_halo,nb_follow_up_halos,MPI_INTEGER8,main_process_rank,io_tag,MPI_COMM_WORLD,statut,ierror) ! Receive the explicit list
            endif 
        end if
        
        if (allocated(nb_of_follow_up_files_per_proc))   deallocate(nb_of_follow_up_files_per_proc)
        if (allocated(list_of_follow_up_files_per_proc)) deallocate(list_of_follow_up_files_per_proc)
    
        return
    end subroutine create_follow_up_list
  
    ! *************************************************
    ! end contains for IO_distribute_tree_files
  
  end subroutine IO_distribute_tree_files
    
  !*****************************************************************************************************************
  !
  ! FUNCTIONS
  !
  !*****************************************************************************************************************
  
  function IO_get_tmp_filelist_unit(ts,r)
  
    implicit none 
    
    integer(kind=4),intent(in)  :: ts  ! timestep
    integer(kind=4),intent(in)  :: r   ! rank
    integer(kind=4)             :: IO_get_tmp_filelist_unit
    
    IO_get_tmp_filelist_unit = ts * 10000 + tmp_filelist_unit + r
    
    return
  end function IO_get_tmp_filelist_unit
  
  !*****************************************************************************************************************
  
  
  function IO_get_log_process_file_unit(r)
  
    implicit none 
    
    integer(kind=4),intent(in)  :: r  ! rank of the process
    integer(kind=4)             :: IO_get_log_process_file_unit
    
    IO_get_log_process_file_unit = log_process_unit + r 
    
    return
  end function IO_get_log_process_file_unit
  
  !*****************************************************************************************************************

  function IO_get_temporary_output_files_unit(ts)

    ! RETURN THE FILE UNIT OF A GIVEN TMP OUTPUT FILE

    implicit none

    integer(kind=4),intent(in)  :: ts  ! the timestep index
    integer(kind=4)             :: IO_get_temporary_output_files_unit

    ! rank is the label of the MPI process (it is saved in the MPI_process module)
    ! the file unit "u" is build with the label of the MPI process and the timestep 
    ! this file unit will be re-build (with the same formulae) at this end of the program to print and to close this same file.
    IO_get_temporary_output_files_unit = 1000*(rank+1)+ts
    ! 1001 for the MPI process 0 and the firts output timestep to
    ! nb_physical_processes * 1000 + ts for the last MPI process used and the last output timestep 

    return
  end function IO_get_temporary_output_files_unit 
  
  !*****************************************************************************************************************
  
  function IO_return_ifile_from_file_index(ind)
  
    ! RETURN THE IFILE INDEX ASSOCIATED TO AN INPUT TREE-FILE INDEX
    ! The tree file index is build as followed: 'nsteps.ifile'
  
    implicit none
    
    integer(kind=4)                :: IO_return_ifile_from_file_index
    integer(kind=4)                :: loc_point
    
    character(len=10),intent(in)   :: ind ! index of the input tree-file
    character(len=10)              :: ifile
    
    loc_point = scan(ind, '.')
    ifile = ind(loc_point+1:)
    read(ifile,'(I10)') IO_return_ifile_from_file_index
    
    return
  end function IO_return_ifile_from_file_index
  
  !*****************************************************************************************************************
 
  function IO_generate_HID(ifile,ts,ih)
  
    ! RETURN THE HALO ID ASSOCIATED TO A HALO 'ih' EVOLVING AT 'ts' AND SAVED IN THE 'ifile' TREE-FILE
    
    implicit none
    
    integer(kind=4),intent(in)   :: ifile  ! the index of the tree-file
    integer(kind=4),intent(in)   :: ih     ! the index of the halo at time step ts
    integer(kind=4),intent(in)   :: ts     ! the time-step of the halo
    integer(kind=8)              :: IO_generate_HID
    
    IO_generate_HID = int(ifile,8)*int(100000000,8)+int(ts,8)*int(100000,8)+int(ih,8)
    
    return
    
  end function IO_generate_HID
  
  !*****************************************************************************************************************
  
  function IO_return_ifile_from_HID(HID)
  
    ! RETURN THE IFILE INDEX ASSOCIATED TO A GIVEN HID
    
    implicit none
    
    integer(kind=8),intent(in)   :: HID ! the halo identification number
    integer(kind=4)              :: IO_return_ifile_from_HID
    
    IO_return_ifile_from_HID = int(HID/int(100000000,8),4)
    
    return
    
  end function IO_return_ifile_from_HID
  
  !*****************************************************************************************************************
  
  function IO_return_ts_from_HID(HID)
  
    ! RETURN THE TIMESTEP INDEX ASSOCIATED TO A GIVEN HID 
  
    implicit none
    
    integer(kind=8),intent(in)   :: HID ! the halo identification number
    integer(kind=4)              :: IO_return_ts_from_HID
    
    IO_return_ts_from_HID = int((HID - int(100000000,8)*IO_return_ifile_from_HID(HID))/int(100000,8),4)
    
    return
    
  end function IO_return_ts_from_HID

  !*****************************************************************************************************************
  
  function IO_return_ih_from_HID(HID)
  
    ! RETURN THE HALO INDEX ASSOCIATED TO A GIVEN HID 
  
    implicit none
    
    integer(kind=8),intent(in)   :: HID ! the halo identification number
    integer(kind=4)              :: IO_return_ih_from_HID
    
    IO_return_ih_from_HID = int(HID - int(100000,8)*int(HID/int(100000,8)),4)
    
    return
    
  end function IO_return_ih_from_HID
  
  !*****************************************************************************************************************  

end module IO
