program G_A_S

  ! main program

  use IO_next  ! IO interface

  implicit none 
 
  integer(kind=4)         :: i        ! loop under ifiles.
  
  character(MAXPATHSIZE)  :: message  ! a message
  !
  ! Initialize 
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nbproc, ierror)  ! return the number of processes used in the computation
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)   
  !
  ! Print welcome message
  if (rank .eq. 0) then
    write(*,'(a)') ''
    write(*,'(a)') ''
    write(*,'(a)') ''
    write(*,'(a)') ''
    write(*,'(a)') '       GGGGGGG        AAAAAAAA       SSSSSSSS      '
    write(*,'(a)') '      G       G      A        A     S              '
    write(*,'(a)') '      G              A        A     S              '
    write(*,'(a)') '      G   GGGGG      AAAAAAAAAA      SSSSSSSS      '
    write(*,'(a)') '      G       G      A        A              S     '
    write(*,'(a)') '      G       G  ..  A        A  ..          S  .. '
    write(*,'(a)') '       GGGGGGG   ..  A        A  ..  SSSSSSSS   .. '
    write(*,'(a)') ''
    write(*,'(a)') ''
    write(*,'(a)') ''
    write(*,'(a)') '****************************************************'
    write(*,'(a)') '* The Galaxy Assembler from dark-matter Simulation *' 
    write(*,'(a)') '*                                                  *'
    write(*,'(a)') '*                       LAM                        *' 
    write(*,'(a)') '*       UMR7326 Technopole de Chateau Gombert      *'
    write(*,'(a)') '*                    MARSEILLE                     *'
    write(*,'(a)') '*                                                  *'
    write(*,'(a)') '*                 COUSIN Morgane                   *'
    write(*,'(a)') '*           morgane.cousin86@gmail.com             *'
    write(*,'(a)') '*                                                  *'
    write(*,'(a)') '****************************************************'
    write(*,'(a)') ''
    write(*,'(a)') ''
    !
  end if ! MPI main_process
  !             
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  !
  ! Atribute process labels (main, physical or luminous) 
  call PrDi_attribute_process_function   
  !
  if (main_process) then 
    call date_and_time(date=date,time=time)
    write(message,'(13(a))') 'Run start: ', date(7:8), '/', date(5:6), '/', date(1:4), ' at ', time(1:2), 'h ', time(3:4), 'mn ', time(5:), 'sec' 
    call IO_print_message(message)  
  end if      
  !
  ! Read parameter file 
  call IO_read_parameter_file
  !
  ! Open log process files
#ifndef RECOVERY_MODE 
! ------------------------------------------------- 
  call IO_open_process_log_files 
  !
  ! Distribute all tree-file into the differents MPI process
  call IO_distribute_tree_files
! -------------------------------------------------
#endif  
  !
  ! Open temporary output files. For a given list of time-step, all properties of halos "computed" are saved 
  ! These tmp output files are then used to build the single output FITS file. 
  call IO_open_and_distribute_temporary_output_files
  ! 
#ifdef TREE_EVOLVE
! -------------------------------------------------
  ! Read input library
  call IO_print_message('Read input libraries')
  !
  ! gas
  ! Read gas properties 
  call gas_read_gas_properties
  ! gas 
  call gas_read_emptying_timescale
  call gas_read_ngc_table 
  ! Set Inter-Galactic Medium
  call gas_igm_initialize  
  !
  ! Cooling process
  ! Read cooling efficiency table
  call cooling_read_cooling_efficiency     
  ! Real thermal instability parameter
  call cooling_read_thermal_instability                 
  !
  ! Stellar population
  ! Read library given the stellar mass loss rate and some other parmeter such as nAgeBins, StellarTimeStep etc ... 
  call stellar_population_read_mass_loss_rates
  ! Read library given the instantaneous SN event rate
  call stellar_population_read_SN_rates
  ! Read stellar spectra library
  call stellar_population_read_stellar_spectra
  !
#ifdef LUMINOUS_PROCESSES
! ------------------------------------------------- 
  ! Dust
  ! Read absorbtion and scattering properties of dusts
  call dust_read_abs_sca_properties
  ! Read spectral energy distributions
  call dust_read_dust_SEDs
  !
  ! AGN
#ifdef AGN_SPECTRUM
! -------------------------------------------------      
  ! Read AGN spectral energy distribution 
  call agn_read_agn_SEDs
! -------------------------------------------------  
#endif
! AGN_SPECTRUM
  !
  ! Filters
  call filters_load_filters
! -------------------------------------------------  
#endif
! LUMINOUS_PROCESSES
  !
  ! Initialize some global parameters linked to the dark-matter n-body simulation
  call tree_set_reference_mass
! -------------------------------------------------  
#endif
! TREE_EVOLVE  
  !
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  !
  call cpu_time(start)  ! for each MPI process we save the starting computation time 
  if (main_process) then 
    !
    call IO_print_message('Now compute tree evolution') 
    !
    ! build the z_table by collecting data from all processes at all ts
#ifndef RECOVERY_MODE 
! ------------------------------------------------- 
    call IO_create_build_and_print_z_table
! -------------------------------------------------  
#endif    
    !
  end if 
  !
#ifndef RECOVERY_MODE 
! -------------------------------------------------   
  if (physical_process) then
    !
    ! PHYSICAL PROCESSES (compute physical properties évolution) 
    !
    if (nb_files_computed .gt. 0) then
        ! RUN MAIN LOOP on input files ****************
        do i = 1, nb_files_computed  ! nb_files_computed are defined in the PrDi module
            !
            ! The "follow" mode allows to follow a structure through cosmic time with:
            !   - a time-line and a merger-tree structure    
            ! Create and open "follow_up", "merger_tree" output files 
            call IO_create_and_open_follow_up_and_merger_tree_output_files(i)
            !
            ! read input files: stepsfile,treefile and propsfile 
            call tree_read(i)
            ! The G.A.S. model uses results from a pure dark-matter N-body simulation
            ! Properties of dark-matter structures are extracted from that simulation
            ! In addition to these first order properties, the G.A.S. model need some other properties
            call tree_compute_halo_derived_properties
            !
            ! ************* TREE EVOLVE *****************
            ! the routine "tree_evolve" computes the evolution processes through merger-trees
            call tree_evolve
            ! ************* TREE EVOLVE *****************
            !
            ! Close "merger_tree" and "follow_up" files
            call IO_close_follow_up_and_merger_tree_output_files(i)
            !
        end do
        ! END OF THE MAIN LOOP ************************
        !
    end if
  end if
  !
#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------  
  if (luminous_process) then
    !
    ! LUMINOUS PROCESSES (compute stellar spectra) 
    !
#ifdef TREE_EVOLVE
! -------------------------------------------------
    if (nb_files_computed .gt. 0) then
        ! RUN MAIN LOOP on input files ****************
        do i = 1, nb_files_computed  ! nb_files_computed are defined in the PrDi module
            !
            ! ************* TREE SPECTRUM *****************
            ! the routine "tree_spectrum" allows to compute galaxies spectrum trhought merger-trees
            call tree_spectrum
            ! ************* TREE SPECTRUM *****************
            !
        end do
    end if
! -------------------------------------------------
#endif
! TREE_EVOLVE
    !
  end if
! ------------------------------------------------- 
#endif
! LUMINOUS_PROCESSES
!
! -------------------------------------------------
#endif
!RECOVERY_MODE
  !  
  call cpu_time(finish)
  !
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  !
  if (.not. main_process) then
    hour = floor((finish-start)/60./60.)            ! number of hours
    mn   = floor(((finish-start)-hour*60.*60.)/60.) ! number of min
    sec  = (finish-start)-hour*60.*60-mn*60.        ! number of sec
    write(message,'(a,i2.2,a,i2.2,a,f5.2,a)') 'Tree evolution completed in ', int(hour), 'h ', int(mn), 'mn ', sec, 'sec'
    call IO_print_message(message,only_rank=rank)
    !
    ! close tmp output file
    call IO_close_temporary_output_files
    !
#ifndef NO_OUTPUTS
! -------------------------------------------------
    ! create FITS headers
    call IO_create_FITS_header
    !
    ! Create open and print in main output FITS files
    call IO_create_open_and_print_in_main_FITS_output_files
! -------------------------------------------------    
#endif
! NO_OUTPUTS
    !     
  end if
  !
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  !
  if (main_process) then
    !
#ifndef RECOVERY_MODE 
! -------------------------------------------------      
    ! Print log file
    call IO_print_log_file
#ifndef NOT_REMOVE_TMP_FILES
! ------------------------------------------------- 
    ! Delete all tmp output files
    call IO_remove_all_tmp_output_files
! ------------------------------------------------- 
#endif
! NOT_REMOVE_TMP_FILES
! -------------------------------------------------    
#endif
! RECOVERY_MODE 
    !
    call date_and_time(date=date,time=time)
    write(message,'(13(a))') 'Run done: ', date(7:8), '/', date(5:6), '/', date(1:4), ' at ', time(1:2), 'h ', time(3:4), 'mn ', time(5:), 'sec' 
    call IO_print_message(message)
    !
    ! print bye bye message
    write(*,'(a)') ''
    write(*,'(a)') ''
    write(*,'(a)') '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
    write(*,'(a)') '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@$@@@@@@@@@@@@@@@@@@@@@@@@@@'
    write(*,'(a)') '@@@@@@@@@@@@@@@@@@@@@@EB@BE@@@@@@@@@@@@@Q@@@@@@@@@@@@@@@'
    write(*,'(a)') '@@@@@@@@@@@@@@@@BEEEEE5E5555SEBB@@@@@E@@@E@E@@@@@@@@@@@@'
    write(*,'(a)') '@@@@@@@@@@@@@@EEE5EEEEtztZ35EEEF5SB@EEEBEEE@Q@@@@@@@@@@@'
    write(*,'(a)') '@@@@@@@@@@@@@E5zE@EEEE3tt=tz33BE5E3EEEE@E@E@@@@@@@@@@@@@'
    write(*,'(a)') '@@@@@@@@@@@@EE35EEE@EE1:.;;;@EEEE@E3EEEK25@@B@@@@@@@@@@@'
    write(*,'(a)') '@@@@@@@@@@@EEEt3@@@@@@@Ek233EE3@ESEEE23@E@@@@@@@@@@@@@@@'
    write(*,'(a)') '@@@@@@@@@@B@@Et!3@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
    write(*,'(a)') '@@@@@@@@@@@@@EEzJ3E@@@@@@@@@@@@@@@@@@@@@@@@@@@@@B@@@@@@@'
    write(*,'(a)') '@@@@@@@@@@@@@@@@kz33E@@@@@@@@@@@@@@@@@@@@@@@@BE@@@@@@@@@'
    write(*,'(a)') '@@@@@@@@@@@@@@@@@@@z235SE$B$B@@@@@@@@@@@@@@@@@@@@@@@@@@@'
    write(*,'(a)') '@@@@@@@@@@@@@@@@@@@@@@@@EEE2EEESEEEEE@@@@@@@@@@@@@@@@@@@'
    write(*,'(a)') '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
    write(*,'(a)') '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
    write(*,'(a)') '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
    write(*,'(a)') '@@@@@@@@@@@@@@@@B@@@@@@@@@@@@@@@B@@@@@@@@@@@@@@@@@@@@@@@'
    write(*,'(a)') '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@BE5@@@@@@@@@@@@@@@@@@@@@@@'
    write(*,'(a)') '@@@@@@@@@@@@@@@@@@@@@@@@@@BEEEEg@@@@@@@@@@@@@@@@@@@@@@@@'
    write(*,'(a)') '@@@@@@@@@@@@@@@@@@@@@@BE5=!7;@@@@@@@@@@@@@@@@@@@@@@@@@@@'
    write(*,'(a)') '@@@@@@@@@@@@@@@@@@@@EF::t;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
    write(*,'(a)') '@@@@@@@@@@@@@@@@@B5j3@@EE@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
    write(*,'(a)') '@@@@@@@@@@@@@@@@E5@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
    write(*,'(a)') '@@@@@@@@@@@@@@@E@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
    write(*,'(a)') '@@@@@@$@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
    write(*,'(a)') '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
    write(*,'(a)') '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@' 
    write(*,'(a)') ''
    write(*,'(a)') ''
    write(*,'(a)') '********************************************************'
    write(*,'(a)') '*                     END OF G.A.S.                    *'
    write(*,'(a)') '*   The Galaxy Assembler from dark-matter Simulation   *'  
    write(*,'(a)') '*                                                      *'
    write(*,'(a)') '*                   Thanks for use !                   *'
    write(*,'(a)') '*                                                      *'
    write(*,'(a)') '********************************************************'
    write(*,'(a)') ''
    write(*,'(a)') ''  
  endif ! main_process
  !
  ! Close log process files
  call IO_close_process_log_files
  !
  ! deallocate data
  call stellar_population_finalize
  call cooling_function_finalize
  call gas_finalize
  !
  call PrDi_finalize
  !
end program G_A_S
