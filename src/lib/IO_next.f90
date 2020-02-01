module IO_next

  use tree    ! acces to tree structure 

  public
  
  !*****************************************************************************************************************
  ! 
  ! OVERVIEW
  !
  !  This module defined global parameters and procedures in link to data output processes
  !  parameters and procedure linked to log files, tmp files and main FITS output files are defined in this module
  !
  ! SUBROUTINES IN THIS MODULE
  !
  !  IO_create_build_and_print_z_table                          : Allows to create a z table. 
  !    called by  : main program                                     This "z" table saved some global halos/galaxies properties through cosmic time
  !
  !  IO_open_process_log_files                                  : Create and open log files for each used processes (physical and luminous)
  !    called by  : main program
  !
  !  IO_close_process_log_files                                 : Close log files for each used processes (physical and luminous)
  !    called by  : main program
  !
  !  IO_open_and_distribute_temporary_output_files              : Create and open temporary output files. These files allow to saved halos/galaxies properties.
  !    called by  : main program                                      These tmp files will be re-readed at the end of the program and re-formated in a single FITS output file  
  !
  !  IO_close_temporary_output_files                            : Close temporary output files
  !    called by  : main program
  !
  !  IO_remove_all_tmp_output_files                             : Remove all tmp output files opened and used by all MPI process
  !    call by    : main program
  !
  !  IO_create_FITS_header                                      : Build Header for outputs FITS files
  !    call by    : main program
  !
  !  IO_create_and_open_follow_up_and_merger_tree_output_files  : Create and open a "follow_up"/"merger-tree" file.
  !    call by    : main program                                      In the context of a evolution model,
  !                                                                  it is possible to follow properties of a halos through cosmic time (it is a time-line)
  !
  !  IO_close_follow_up_and_merger_tree_output_files            : Close a given "follow-up"/"merger-tree" file 
  !    call by    : tree_evolve                                       
  !
  !  IO_print_log_file                                          : Print in a "log" file all informations about the run
  !    call by    : main program                                      input parameter values, option used ...
  !
  !  IO_create_open_and_print_in_main_FITS_output_files         : Allows to re-format tmp output files and save results in a single FITS file
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
  
  contains

  !*****************************************************************************************************************
  !
  ! SUBROUTINES
  !
  !*****************************************************************************************************************

  subroutine IO_create_build_and_print_z_table

    ! CREATE, BUILD AND PRINT THE "Z"_TABLE TO SAVED GLOBAL PROPERTIES THROUGH COSMIC TIME 
    ! The G.A.S. model can follow some halo/galaxy properties through cosmic time. 
    ! These properties are saved in a z_table for all dark-matter timesteps. 

    implicit none

    integer(kind=4)          :: nb_block_receive
    integer(kind=4)          :: current_ts, l, nb_halos
    
    character(MAXPATHSIZE)   :: filename
    character(MAXPATHSIZE)   :: legend
    !
    real(kind=8),allocatable :: z_table(:,:)
    real(kind=8),allocatable :: buffer(:)
    
    if (main_process) then
      call IO_print_message('IO_create_build_and_print_z_table')    
      call IO_print_message('Create and open redshift evolution properties file')
      ! build the filename using the resolution and the output _path 
      write(filename, '(a,a,i4.4,a)') trim(output_path), '/GAS_res_',res,'_z_table.dat'
      ! Open the file
      open(unit = z_table_unit, file = filename, status = 'new')
      ! Allocate z table
      allocate(z_table(nsteps,n_z_data-1)) 
      z_table(:,:) = 0.d0                ! init
      ! allocate the buffer
      allocate(buffer(n_z_data-1))          
      buffer(:) = 0.d0                   ! init
      !
      nb_block_receive = 0
      do while (nb_block_receive .lt. (nb_tree_files_computed*ts_STOP))
        ! information from each tree file must be receive, each tree-file give the evolution from ts = 1, to ts = nsteps
        ! 
        ! receive information about the current ts
        call MPI_RECV(current_ts,1,MPI_INTEGER4,MPI_ANY_SOURCE,ztable_tag,MPI_COMM_WORLD,statut,ierror)  
        ! 
        ! receive information about the number of computable halo at this ts
        call MPI_RECV(nb_halos,1,MPI_INTEGER4,MPI_ANY_SOURCE,ztable_tag+1,MPI_COMM_WORLD,statut,ierror)    
        !
        if (nb_halos .gt. 0) then
            !
            ! receive data 
            call MPI_RECV(buffer,n_z_data-1,MPI_REAL8,MPI_ANY_SOURCE,ztable_tag+2,MPI_COMM_WORLD,statut,ierror) 
            ! sum data
            z_table(current_ts,3:) = z_table(current_ts,3:) + buffer(3:)/real(nb_halos,8)
            ! earse data for redshift and universe age
            z_table(current_ts,1:2) = buffer(1:2)
        end if
        !
        ! a new block has been received
        nb_block_receive = nb_block_receive +1
      end do
      !
      ! the table is completed
      ! apply volume factor
      z_table(current_ts,3:) = z_table(current_ts,3:)/(L_box**3.)
      ! print it
      ! write the legend
      ! init the legend
      legend = '# '
      do l = 1, n_z_data
        write(legend,'(a,i2.2,a,a,a)') trim(legend), l, ') ', z_table_legend(l), ' | '
      end do
      write(z_table_unit,'(a)') trim(legend)
      legend = '# '
      do l = 1, n_z_data
        write(legend,'(a,a,a,a)') trim(legend), '    ', z_table_legend_unit(l), ' | '
      end do
      write(z_table_unit,'(a)') trim(legend)
      write(legend,'(a,i2.2,a)') '(i2.2,',n_z_data, '(2x,E14.6))' 
      do current_ts = 1, nsteps
        write(z_table_unit,legend) current_ts, z_table(current_ts,:)
      end do
    end if

    call IO_print_message('IO_create_build_and_print_z_table',only_rank=rank)
    
    return
  end subroutine IO_create_build_and_print_z_table

  !***************************************************************************************************************** 
  
  subroutine IO_open_process_log_files
  
    implicit none
    
    integer(kind=4)              :: u
    
    character(len=MAXPATHSIZE)   :: filename
    
    if (.not. main_process) then
      !
      ! get file unit
      u = IO_get_log_process_file_unit(rank)
      !
      ! build the filename
      write(filename, '(a,a,i3.3,a)') trim(output_path), '/GAS_process_',rank,'.log'
      !
      ! open the file
      open(unit = u, file = trim(filename), status = 'new')
      !     
      write(u,'(a,i3.3)') '# G.A.S. log file associated to process: ', rank
      write(u,'(a,a)') '# Process runing on: ', trim(nodename)
    end if
    
    return
  end subroutine IO_open_process_log_files
  
  !***************************************************************************************************************** 
  
  subroutine IO_close_process_log_files
  
    implicit none
    
    integer(kind=4)              :: u
    
    if (.not. main_process) then
      !
      ! get file unit
      u = IO_get_log_process_file_unit(rank)
      !
      ! open the file
      close(u)
      !     
    end if
    
    return
  end subroutine IO_close_process_log_files

  !*****************************************************************************************************************
    
  subroutine IO_open_and_distribute_temporary_output_files

    ! OPEN AND DITRIBUTE TEMPORARY OUTPUT FILE
    ! we print result in temporary files (one per output timestep and per process in MPI mode)
    ! after the computation each phy (resp. lum) process print a set of main phy (resp. lum) files

    implicit none

    integer(kind=4)              :: ots  ! index loop under output timestep list   
    integer(kind=4)              :: r,r2 ! index loop under rank   
    integer(kind=4)              :: u    ! index loop under file unit
    integer(kind=4),allocatable  :: printer_process_list(:,:)
    integer(kind=4),allocatable  :: nb_files_per_printer_process(:)
    
    character(len=MAXPATHSIZE)   :: filename
    character(len=MAXPATHSIZE)   :: path
    character(MAXPATHSIZE)       :: format

    call IO_print_message('Open and distribute tmp output files')
    
    ! each physical (luminous) process print result in nout files (1 per output timestep)
    ! each files is splitted in some blocks (1 block per tree_file computed by this process), 
    ! A tmp output file is therfore composed of "nb_files_computed" blocks of data
    !
    if (main_process) then
      !
      ! all tmp filenames are listed in differents files (one per ts)
      !
      allocate(printer_process_list(0:nbproc-1,nout))
      printer_process_list = 0
      allocate(nb_files_per_printer_process(0:nbproc-1))
      nb_files_per_printer_process = 0
      !
      ots = 1  ! init the output timestep list index
      !
      do while (ots .le. nout) ! tsout contains the list of output timesteps
        !
        do r = 1, max(nb_physical_processes, nb_luminous_processes)
          !
#ifndef ONLY_LUM 
! -------------------------------------------------           
          r2 = list_of_physical_processes(r)
          !
          nb_files_per_printer_process(r2) = nb_files_per_printer_process(r2) +1 
          printer_process_list(r2,nb_files_per_printer_process(r2)) = tsout(ots) 
          !
#ifndef RECOVERY_MODE 
! -------------------------------------------------           
          u = IO_get_tmp_filelist_unit(tsout(ots),r2)
          write(filename, '(a,a,i3.3,a)') trim(output_path),'/GAS_phy_tmp_filename_list_ts_',tsout(ots),'.tmp'
          open(unit = u, file = trim(filename), status = 'new')
          ! print the number of path listed after
          write(u,'(i3.3)') nb_physical_processes
          ! print the explicit list of phy_tmp filename in the list
          do r2 = 1, nb_physical_processes
            write(path, '(a,a,i4.4,a,i3.3,a,i3.3,a)') trim(output_path), &
                '/GAS_phy_res_',res,'_ts_',tsout(ots),'_rank_',list_of_physical_processes(r2),'.tmp'  ! res is the resolution of the N-body simualtion used
            write(u,'(a)') trim(path)    
          end do 
          close(u)
! ------------------------------------------------- 
#endif
! RECOVERY_MODE  
#endif
! ONLY_LUM         
          !
#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------  
          !
          r2 = list_of_luminous_processes(r)
          !
          nb_files_per_printer_process(r2) = nb_files_per_printer_process(r2) +1 
          printer_process_list(r2,nb_files_per_printer_process(r2)) = tsout(ots) 
          !
#ifndef RECOVERY_MODE 
! -------------------------------------------------           
          u = IO_get_tmp_filelist_unit(tsout(ots),r2)
          write(filename, '(a,a,i3.3,a)') trim(output_path),'/GAS_lum_tmp_filename_list_ts_',tsout(ots),'.tmp'
          open(unit = u, file = trim(filename), status = 'new')
          ! print the number of path listed after
          write(u,'(i3.3)') nb_luminous_processes
          ! print the explicit list of lum_tmp filename in the list
          do r2 = 1, nb_luminous_processes
            write(path, '(a,a,i4.4,a,i3.3,a,i3.3,a)') trim(output_path), &
                '/GAS_lum_res_',res,'_ts_',tsout(ots),'_rank_',list_of_luminous_processes(r2),'.tmp'  ! res is the resolution of the N-body simualtion used
            write(u,'(a)') trim(path)    
          end do
          close(u)
! ------------------------------------------------- 
#endif
! RECOVERY_MODE           
          !
! ------------------------------------------------- 
#endif
! LUMINOUS_PROCESSES   
          !
          ots = ots +1  
          if (ots .gt. nout) exit
          !
        end do
      end do 
      !
      ! send list to phy and lum porcesses
      do r = 0, nbproc-2
        ! send to process r, the number of main output files that it have to print 
        call MPI_SEND(nb_files_per_printer_process(r),1,MPI_INTEGER4,r,io_tag+1,MPI_COMM_WORLD,ierror)
        if (nb_files_per_printer_process(r) .gt. 0) then
           ! Send the explicit list ts
           call MPI_SEND(printer_process_list(r,1:nb_files_per_printer_process(r)), &
              nb_files_per_printer_process(r),MPI_INTEGER4,r,io_tag+2,MPI_COMM_WORLD,ierror)
           !
           if (nb_files_per_printer_process(r) .le. 10) then 
            ! print ts index of files that the current process have to print
            write(format,'(a,i3.3,a)') '(a,i3.3,a,', nb_files_per_printer_process(r), '(i3.3,2x))'
            write(errunit,format) '|---> MPI process ', r, ' will print main output files for ts: ', printer_process_list(r,1:nb_files_per_printer_process(r))
          else
            ! print only the two first and the two last ts index
            write(errunit,'(a,i3.3,a,i3.3,a,2(i3.3,2x),a,2(i3.3,2x),a)') '|---> MPI process ', r, &
                  ' will print ', nb_files_per_printer_process(r), ' main output files : ( ', printer_process_list(r,1:2), '...  ', &
                  printer_process_list(r,nb_files_per_printer_process(r)-1:nb_files_per_printer_process(r)), ' )'
          end if
        end if
      end do
      !
      deallocate(printer_process_list)
      deallocate(nb_files_per_printer_process)
      !
    else
      ! I am a physical or a luminous process
      ! receive the number of main output files that I have to print 
      call MPI_RECV(nb_main_output_files_to_print,1,MPI_INTEGER4,main_process_rank,io_tag+1,MPI_COMM_WORLD,statut,ierror)  
      if (nb_main_output_files_to_print .gt. 0) then
        allocate(list_of_ts_to_print(nb_main_output_files_to_print))
        ! Receive the explicit list
        call MPI_RECV(list_of_ts_to_print,nb_main_output_files_to_print,MPI_INTEGER4,main_process_rank,io_tag+2,MPI_COMM_WORLD,statut,ierror)
      end if
      !
      ! open tmp files associated to the current process
      do ots = 1, nout
        if (physical_process) then  
          write(filename, '(a,a,i4.4,a,i3.3,a,i3.3,a)') trim(output_path), &
              '/GAS_phy_res_',res,'_ts_',tsout(ots),'_rank_',rank,'.tmp'  ! res is the resolution of the N-body simualtion used  
        endif
        !
#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------          
        if (luminous_process) then   
          write(filename, '(a,a,i4.4,a,i3.3,a,i3.3,a)') trim(output_path), &
              '/GAS_lum_res_',res,'_ts_',tsout(ots),'_rank_',rank,'.tmp'   ! res is the resolution of the N-body simualtion used  
        endif                                   
! ------------------------------------------------- 
#endif
! LUMINOUS_PROCESSES  
        !
#ifndef RECOVERY_MODE 
! -------------------------------------------------         
        ! open tmp file
        if ((physical_process) .or. (luminous_process)) then 
            u =  IO_get_temporary_output_files_unit(tsout(ots))
            open(unit = u, file = trim(filename), status = 'new',form = 'unformatted')
            write(u) nb_files_computed, tsout(ots), ots  ! print : - nb_files_computed corresponding to the number of blocks print in this tmp output file,
                                                        !         - ts, the timestep value of the N-body simulation used 
                                                        !         - ots (index over the output timestep list) 
        end if
! ------------------------------------------------- 
#endif
! RECOVERY_MODE         
      end do
    end if
    
    return
  end subroutine IO_open_and_distribute_temporary_output_files
    
  !*****************************************************************************************************************

  subroutine IO_close_temporary_output_files

    ! CLOSE TMP OUTPUT FILE

    integer(kind=4)        :: ots  ! index loop under output timestep list
    integer(kind=4)        :: ts   ! index loop under current timestep    
    integer(kind=4)        :: u    ! index loop under file unit

    call IO_print_message('Close tmp output files')

    if (physical_process .or. luminous_process) then    
      ots = 1          ! init the output timestep list index

      do ts = 1, nts   ! nts is the number of output timesteps of the N-body simulation
                       ! nts is a global variable save in the IO header module.  
        if(ts .eq. tsout(ots)) then ! tsout save the list of output timesteps
          ! each processus print result in nout files (1 per output timestep)
          u = IO_get_temporary_output_files_unit(ts)
          close(unit = u)
          ! update output timestep loop "ots"
          if (ots .eq. nout) then  ! nout is the number of output timestep
            ots = 1
          else
            ots = ots + 1  
          end if
        end if ! ts = tsout
      end do ! end ts 
    end if  

    return
  end subroutine IO_close_temporary_output_files

  !*****************************************************************************************************************  

  subroutine IO_remove_all_tmp_output_files

    ! REMOVE ALL TMP OUTPUT FILES
    ! The code used various tmp output files. 
    ! These file are used, at the end of the computation to generate main FITS files

    implicit none

    character(MAXPATHSIZE) :: command

    call IO_print_message('Remove tmp files')
    !
    write(command, '(a,a,a)') 'rm -f ', trim(output_path), '/*.tmp'
    call system(command)

    return
  end subroutine IO_remove_all_tmp_output_files  

  !***************************************************************************************************************** 

  ! OUTPUT FITS FILE
  ! The G.A.S. model saves halos/galaxies properties in differents output FITS files.
  ! It is possible (FOLLOW option) to follow time evolution of some specific halos/galaxies. 
  ! These results are printed in two different files: time-line or merger-tree formats

  !***************************************************************************************************************** 

  subroutine IO_create_FITS_header

    implicit none

#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------    
    integer(kind=4)            :: ifilt ! loop index onto filter list
! ------------------------------------------------- 
#endif
! LUMINOUS_PROCESSES      

    ! CREATE ALL FITS HEADERS
    ! G.A.S. results (physical properties of all galaxies at various timesteps) are printed in FITS files.
    ! The main FITS file is build at the end of the program. We used all tmp output files.
    ! If the FOLLOW_UP keyword is turned ON, its is possible to follow the evolution of individual galaxies through cosmic time. 
    ! This kind of "followed" results are printed in a different type of FITS file.
    ! In addition, if the FOLLOW_UP keyword is turned ON, it is possible to follow evolution of a complete merger tree linked to a halo structure
    ! This kind of "merger-tree" results are printed in a different type of FITS file.

    call IO_print_message('Create FITS headers')

    ! Create header informations
    !
#ifdef TREE_EVOLVE    
! -------------------------------------------------               
     BARYON = .true.   ! If the TREE_EVOLVE keyword is turned ON, the model must compute baryonic process.
                       ! The output FITS file contains a large number of physical properties   
! ------------                       
#else
! ------------
     BARYON = .false.  ! If the TREE_EVOLVE keyword is turn OFF, only the dard matter structure evolution is done
                       ! The output FITS file contains a limited number of field (only dark matter halos properties)
! -------------------------------------------------                       
#endif
! endif TREE_EVOLVE
    ! 
    ! Model with baryonic processes followed
    if (BARYON) then
      !  
      ! we must print all data (tree/dm/bh/gal/disc/starsd/dustISMd/dustBCd/agn/bulge/starsb/dustb)
      nb_field = nb_tree_field + nb_dm_field + nb_bh_field + nb_gal_field + &
                nb_disc_field + nb_stars_field + 2*nb_dust_field + nb_agn_field + &
                nb_bulge_field + nb_stars_field + nb_dust_field

      ! allocate ttype, tform and tunit array. These arrays contains FITS header information about physical properpies saved by the G.A.S. model.  
      allocate(ttype(nb_field))  ! Name of each output colomn data
      allocate(tform(nb_field))  ! Data type of each column data
      allocate(tunit(nb_field))  ! Physical unit of each column data
      !
      ! build the ttype list (the list of field name). Properties are group by topics (tree structure, dark matter properties then halo and galaxies properties)
      ttype = (/ttype_tree,      & ! tree properties
                ttype_dm,        & ! dm properties
                ttype_bh,        & ! baryon_halo properties
                ttype_gal,       & ! galaxy properties
                ttype_disc,      & ! disc properties
                ttype_starsd,    & ! disc stars properties
                ttype_dust_ISMd, & ! dust properties, ISM
                ttype_dust_BCd,  & ! dust properties, BC
                ttype_agn,       & ! agn properties
                ttype_bulge,     & ! bulge properties
                ttype_starsb,    & ! bulge stars properties
                ttype_dustb/)      ! dust properties
                    
      !
      ! build the list of physical unit, in general distance are given in kpc, times are given in Gyr, masses are given in Msun and velocities in km/s
      tunit = (/tunit_tree,   & ! tree properties
                tunit_dm,     & ! dm properties
                tunit_bh,     & ! baryon_halo properties
                tunit_gal,    & ! galaxy properties
                tunit_disc,   & ! disc properties
                tunit_starsd, & ! disc stars properties
                tunit_dustd,  & ! dust properties, ISM  
                tunit_dustd,  & ! dust properties, BC 
                tunit_agn,    & ! agn properties  
                tunit_bulge,  & ! bulge properties
                tunit_starsb, & ! bulge stars properties
                tunit_dustb/)   ! dust properties
      !
      ! build the tform array (the list of data format type).
      tform = (/tform_tree,   & ! tree properties
                tform_dm,     & ! dm properties
                tform_bh,     & ! baryon_halo properties
                tform_gal,    & ! galaxy properties
                tform_disc,   & ! disc properties
                tform_starsd, & ! disc stars properties
                tform_dustd,  & ! dust properties, ISM 
                tform_dustd,  & ! dust properties, BC 
                tform_agn,    & ! agn properties 
                tform_bulge,  & ! bulge properties
                tform_starsb, & ! bulge stars properties
                tform_dustb/)   ! dust properties
      !
#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------  
#ifdef GAL_SPECTRA
! ------------------------------------------------- 
      !
      ! build ttype, tunit and tform arrays for spectra 
      ttype_spec = 'spectra '
      tunit_spec = 'micron      ' 
      write(tform_spec,'(i4.4,a)') nWaves, 'E'
      !
! ------------------------------------------------- 
#endif
! GAL_SPECTRA
      !
      ! build ttype, tunit and tform arrays for magnitudes
      ! we use the list of spectra to define ttype_mag
      !
      allocate(ttype_mag(nfilters))
      allocate(tunit_mag(nfilters))
      allocate(tform_mag(nfilters))
      !
      do ifilt = 1, nfilters
        ttype_mag(ifilt) = trim(filters_tab(ifilt)%name)
        tunit_mag(ifilt) = 'AB Mag      '
        tform_mag(ifilt) = '1E'
      end do
      !
! ------------------------------------------------- 
#endif
! LUMINOUS_PROCESSES                            
      !                                                  
    else
      ! In this second case, the baryonic processes are not followed. Only the dark-matter properties are given
      ! We must just print dm data
      nb_field = nb_tree_field + nb_dm_field   ! The number of printed field
      ! allocate ttype, tform and tunit array. These arrays contains FITS header information about physical properpies saved by the G.A.S. model.  
      allocate(ttype(nb_field))  ! Name of each output colomn data
      allocate(tform(nb_field))  ! Data type of each column data
      allocate(tunit(nb_field))  ! Physical unit of each column data
      !
      ! build the ttype list
      ttype = (/ttype_tree, & ! tree properties
                ttype_dm/)    ! dm properties
      !
      ! build the tunit list
      tunit = (/tunit_tree, & ! tree properties
                tunit_dm/)    ! dm properties
            !
      ! build the tform list
      tform = (/tform_tree, & ! tree properties
                tform_dm/)    ! dm properties
    !
    end if ! end BARYON  

  end subroutine IO_create_FITS_header

  !*****************************************************************************************************************

  ! "FOLLOW UP" AND "MERGER_TREE" FILES (DEFINITION, OPEN PROCEDURE, PRINTING PROCEDURE ECT ... )

  subroutine IO_create_and_open_follow_up_and_merger_tree_output_files(ifile)

    ! CREATE ANS OPEN A GIVEN "FOLLOW_UP" FILE OR A GIVEN "MERGER_TREE" FILE
     
    implicit none

    integer(kind=4)                            :: ifile       ! index of the input mergertree file
    integer(kind=4)                            :: typ,ih      ! loop index (under output file type and follow_up halos)
    integer(kind=4), pointer                   :: pu          ! a pointer throught the file unit

    character(MAXPATHSIZE)                     :: filename    ! a filename
    character(MAXPATHSIZE)                     :: file        ! a filename
    character(MAXPATHSIZE)                     :: message     ! a message 
    character(len=11),dimension(2),parameter   :: output_type = (/'follow_up  ','merger_tree'/) ! this subroutine is used to create and open two different files, 
                                                                                                ! a follow_up file and a merger tree_file

    if (physical_process) then
        if (nb_follow_up_halos .gt. 0) then
            !
            ! create the follow_up_unit list
            if (.not. allocated(follow_up_unit)) allocate(follow_up_unit(nb_follow_up_halos))
            ! create the merger_tree_unit_list
            if (.not. allocated(merger_tree_unit)) allocate(merger_tree_unit(nb_follow_up_halos))
            !
            do ih = 1, nb_follow_up_halos
                !
                if (IO_return_ifile_from_HID(list_of_follow_up_halo(ih)) .eq. &
                        IO_return_ifile_from_file_index(list_of_computed_files(ifile))) then
                    !
                    ! open only follow_up and merger tree used in the current merger tree 
                    !
                    do typ = 1, 2 ! loop under the two different output files, follow_up and merger_tree
                        !
                        ! reset the FITSIO status (= 0.)
                        status = 0
                        !
                        select case (trim(output_type(typ)))
                        case ('follow_up')
                            ! create the file with unit "follow_up_unit" 
                            call ftgiou(follow_up_unit(ih),status)
                            pu => follow_up_unit(ih)
                            ! build the filename  
                            write(filename, '(a,i11.11,a)') 'GAS_follow_up_',list_of_follow_up_halo(ih),'.fits' ! list_of_computed_files is defined in the MPI_process module
                            write(message, '(a,a,a,i2.2,a)') 'Open file (u) : ', trim(filename), ' (', pu, ')'
                            call IO_print_message(message,only_rank = rank)
                            write(file, '(a,a,a)') trim(output_path), '/', trim(filename)
                            !
                        case ('merger_tree')
                            ! create the file with unit "merger_tree_unit" 
                            call ftgiou(merger_tree_unit(ih),status)
                            pu => merger_tree_unit(ih)
                            ! build the filename  
                            write(filename, '(a,i11.11,a)') 'GAS_merger_tree_',list_of_follow_up_halo(ih),'.fits' ! list_of_computed_files is defined in the MPI_process module
                            write(message, '(a,a,a,i2.2,a)') 'Open file (u) : ', trim(filename), ' (', pu, ')'
                            call IO_print_message(message,only_rank = rank)
                            write(file, '(a,a,a)') trim(output_path), '/', trim(filename)
                            !
                        case default 
                            write(message,'(a,a)') 'Unkwnon keyword : ', trim(output_type(typ))
                            call IO_print_error_message(trim(message),called_by='IO_create_and_open_follow_up_and_merger_tree_output_files')  
                        end select
                        !
                        ! check the status     
                        if (status .gt. 0) then  
                            call IO_print_error_message('Cannot create follow_up/merger_tree file',only_rank=rank,called_by='IO_create_and_open_follow_up_and_merger_tree_output_files')  
                            write(message,'(a,i4,a,a)') 'ftgiou status: ', status, ' for file: ', trim(file)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        !
                        ! open a new FITS file
                        call ftinit(pu,trim(file),blocksize,status)  
                        ! check the status
                        if (status .gt. 0) then  
                            call IO_print_error_message('Cannot open follow_up/merger_tree file',only_rank=rank,called_by='IO_create_and_open_follow_up_and_merger_tree_output_files')
                            write(message,'(a,i4,a,a)') 'ftinit status: ', status, ' for file: ', trim(file)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        !
                        ! create extensions
                        ! The follow_up file contains, 6 extensions dedicated to tree structure, dark-matter, baryon-halo (global) galaxy disc and bulge component
                        !
                        ! create extension for tree structure 
                        call ftibin(pu,0,nb_tree_field+1,&
                            (/'Age_Universe          ',ttype_tree/),&
                            (/'1E'                    ,tform_tree/),&
                            (/'Gyr         '          ,tunit_tree/),&
                            'tree',varidat,status) ! here nrowsll is set to 0.
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot create tree extension for follow_up/merger_tree file',only_rank=rank,called_by='IO_create_and_open_follow_up_and_merger_tree_output_files')  
                            write(message,'(a,i4,a,a)') 'ftibin status: ', status, ' for file: ', trim(file)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        ! attribute hdu parameter: "hdu_tree"
                        call ftghdn(pu,hdu_tree) 
                        !
                        ! create extension for dm data 
                        call ftibin(pu,0,nb_dm_field+1,&
                            (/'life_time             ',ttype_dm/), &
                            (/'1E'                    ,tform_dm/), &
                            (/'Gyr         '          ,tunit_dm/), &
                            'dm',varidat,status) ! here nrowsll is set to 0.
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot create dm extension for follow_up/merger_tree file',only_rank=rank,called_by='IO_create_and_open_follow_up_and_merger_tree_output_files')  
                            write(message,'(a,i4,a,a)') 'ftibin status: ', status, ' for file: ', trim(file)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        ! attribute hdu parameter: "hdu_dm"
                        call ftghdn(pu,hdu_dm) 
                        !
#ifdef TREE_EVOLVE
! -------------------------------------------------
                        ! create extension for baryon_halo data 
                        call ftibin(pu,0,nb_bh_field+1, &
                            (/'life_time             ',ttype_bh/), &
                            (/'1E'                    ,tform_bh/), &
                            (/'Gyr         '          ,tunit_bh/), &
                            'bh',varidat,status) ! here nrowsll is set to 0.
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot create bh extension for follow_up/merger_tree file',only_rank=rank,called_by='IO_create_and_open_follow_up_and_merger_tree_output_files')  
                            write(message,'(a,i4,a,a)') 'ftibin status: ', status, ' for file: ', trim(file)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        ! attribute hdu parameter: "hdu_bh"
                        call ftghdn(pu,hdu_bh)
                        !
                        ! create extension for galaxy data
                        call ftibin(pu,0,nb_gal_field+1, &
                            (/'life_time             ',ttype_gal/), &
                            (/'1E'                    ,tform_gal/), &
                            (/'Gyr         '          ,tunit_gal/), &
                            'galaxy',varidat,status) ! here nrowsll is set to 0.
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot create gal extension for follow_up/merger_tree file',only_rank=rank,called_by='IO_create_and_open_follow_up_and_merger_tree_output_files')  
                            write(message,'(a,i4,a,a)') 'ftibin status: ', status, ' for file: ', trim(file)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        ! attribute hdu parameter: "hdu_gal"
                        call ftghdn(pu,hdu_gal)
                        !
                        ! create extension for disc data 
                        ! header informations are build with disc disc_starsd and agn structures
                        call ftibin(pu,0,nb_disc_field+nb_stars_field+2*nb_dust_field+nb_agn_field+1, &
                            (/'life_time             ',ttype_disc,ttype_starsd,ttype_dust_ISMd,ttype_dust_BCd,ttype_agn/), &
                            (/'1E'                    ,tform_disc,tform_starsd,tform_dustd,tform_dustd,tform_agn/), &
                            (/'Gyr         '          ,tunit_disc,tunit_starsd,tunit_dustd,tunit_dustd,tunit_agn/), &
                            'disc',varidat,status) ! here nrowsll is set to 0.
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot create disc extension for follow_up/merger_tree file',only_rank=rank,called_by='IO_create_and_open_follow_up_and_merger_tree_output_files')  
                            write(message,'(a,i4,a,a)') 'ftibin status: ', status, ' for file: ', trim(file)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        ! attribute hdu parameter: "hdu_disc"
                        call ftghdn(pu,hdu_disc) 
                        !
                        ! create extension for bulge data 
                        ! header informations are build with bulge bulge_starsb and agn structures
                        call ftibin(pu,0,nb_bulge_field+nb_stars_field+nb_dust_field+1, &
                            (/'life_time             ',ttype_bulge,ttype_starsb,ttype_dustb/), &
                            (/'1E'                    ,tform_bulge,tform_starsb,tform_dustb/), &
                            (/'Gyr         '          ,tunit_bulge,tunit_starsb,tunit_dustb/), &
                            'bulge',varidat,status) ! here nrowsll is set to 0.
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot create bulge extension for follow_up/merger_tree file',only_rank=rank,called_by='IO_create_and_open_follow_up_and_merger_tree_output_files')  
                            write(message,'(a,i4,a,a)') 'ftibin status: ', status, ' for file: ', trim(file)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        ! attribute hdu parameter: "hdu_bulge"
                        call ftghdn(pu,hdu_bulge) 
#endif
! -------------------------------------------------
! endif TREE_EVOLVE  
                    end do ! typ
                end if
            end do ! ih
        end if ! nb_follow_up_halos > 0
    end if ! physical_process

    return
  end subroutine IO_create_and_open_follow_up_and_merger_tree_output_files

  !*****************************************************************************************************************    

  subroutine IO_close_follow_up_and_merger_tree_output_files(ifile)

    ! CLOSE A "FOLLOW_UP" AND "MERGER_TREE" FILE

    implicit none

    integer(kind=4)                            :: ifile       ! index of the input mergertree file
    integer(kind=4)                            :: typ, ih
    integer(kind=4), pointer                   :: pu 
    
    character(len=11),dimension(2),parameter   :: output_type = (/'follow_up  ','merger_tree'/)
    character(MAXPATHSIZE)                     :: message

    if (physical_process) then
        if (nb_follow_up_halos .gt. 0) then
            do ih = 1, nb_follow_up_halos
                !
                if (IO_return_ifile_from_HID(list_of_follow_up_halo(ih)) .eq. &
                        IO_return_ifile_from_file_index(list_of_computed_files(ifile))) then
                    !
                    ! close only follow_up and merger tree used in the current merger tree 
                    !
                    do typ = 1, 2
                        !
                        select case (trim(output_type(typ)))
                        case ('follow_up')
                            pu => follow_up_unit(ih)
                        case ('merger_tree')
                            pu => merger_tree_unit(ih)
                        case default 
                            write(message,'(a,a)') 'Unkwnon keyword : ', trim(output_type(typ))
                            call IO_print_error_message(trim(message),called_by='IO_close_follow_up_and_merger_tree_output_files')  
                            stop
                        end select
                        !
                        ! close the file referenced by the unit
                        write(message,'(a,i2.2)') 'Close file unit: ', pu
                        call IO_print_message(message,only_rank = rank)
                        call ftclos(pu, status)
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot close follow_up/merger_tree file',only_rank=rank,called_by='IO_close_follow_up_and_merger_tree_output_files')  
                            write(message,'(a,i4)') 'ftclos status: ', status
                            call IO_print_message(message,only_rank=rank)  
                            stop   ! stop the program
                        end if
                        ! deallocate the unit
                        call ftfiou(pu, status)
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot deallocate follow_up_unit/merger_tree_unit',only_rank=rank,called_by='IO_close_follow_up_and_merger_tree_output_files')  
                            write(message,'(a,i4)') 'ftfiou status: ', status
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                    end do ! typ
                end if
            end do ! ih
        end if
    end if   ! physical_process
    
    return
  end subroutine IO_close_follow_up_and_merger_tree_output_files

  !*****************************************************************************************************************

  subroutine IO_print_log_file

    ! WRITE THE LOG FILE

    implicit none
    
    integer(kind=4)              :: i
    integer(kind=4)              :: ni

    character(MAXPATHSIZE)       :: filename

    ! build the filename (with output_path)
    write(filename, '(a,a)') trim(output_path), '/GAS.log'
    ! open the corresponding file
    open(unit = log_file_unit, file = filename, status = 'unknown')
    !
    ! print
    write(log_file_unit,'(a)')      ''
    write(log_file_unit,'(a)')      '       GGGGGGG        AAAAAAAA       SSSSSSSS      '
    write(log_file_unit,'(a)')      '      G       G      A        A     S              '
    write(log_file_unit,'(a)')      '      G              A        A     S              '
    write(log_file_unit,'(a)')      '      G   GGGGG      AAAAAAAAAA      SSSSSSSS      '
    write(log_file_unit,'(a)')      '      G       G      A        A              S     '
    write(log_file_unit,'(a)')      '      G       G  ..  A        A  ..          S  .. '
    write(log_file_unit,'(a)')      '       GGGGGGG   ..  A        A  ..  SSSSSSSS   .. '
    write(log_file_unit,'(a)')      ''
    write(log_file_unit,'(a)')      ''
    write(log_file_unit,'(a)')      ''
    write(log_file_unit,'(a)')      '****************************************************'
    write(log_file_unit,'(a)')      '* The Galaxy Assembler from dark-matter Simulation *'  
    write(log_file_unit,'(a)')      '*                                                  *'
    write(log_file_unit,'(a)')      '*                       LAM                        *' 
    write(log_file_unit,'(a)')      '*       UMR7326 Technopole de Chateau Gombert      *'
    write(log_file_unit,'(a)')      '*                    MARSEILLE                     *'
    write(log_file_unit,'(a)')      '*                                                  *'
    write(log_file_unit,'(a)')      '*                 COUSIN Morgane                   *'
    write(log_file_unit,'(a)')      '*              morgane.cousin@lam.fr               *'
    write(log_file_unit,'(a)')      '*                                                  *'
    write(log_file_unit,'(a)')      '****************************************************'
    write(log_file_unit,'(a)')      ''
    write(log_file_unit,'(a)')      ''
    write(log_file_unit,'(a)')      ''
    write(log_file_unit,'(a)')      '---------------Simulation parameter---------------'
    write(log_file_unit,'(a,i4.4)') 'resolution                  = ', res
    write(log_file_unit,'(a,a)')    'input_path                  = ', trim(input_path)
    write(log_file_unit,'(a,a)')    'output_path                 = ', trim(output_path)
    write(log_file_unit,'(a,i3.3)') 'nb_tree_files_computed      = ', nb_tree_files_computed
    write(log_file_unit,'(a,i1.1)') 'nb_of_desc_halo_test        = ', nb_of_desc_halo_test
    write(log_file_unit,'(a,f4.2)') 'min_dm_goodness_of_branch   = ', min_dm_goodness_of_branch
    write(log_file_unit,'(a,f4.2)') 'min_halo_goodness_of_branch = ', min_halo_goodness_of_branch
    write(log_file_unit,'(a,f4.2)') 'energy threshold            = ', Ecut
    write(log_file_unit,'(a,f4.2)') 'physical_precision          = ', physical_precision
    write(log_file_unit,'(a)')      '------------------Model Options-------------------'
    !
#ifdef HORIZON
! -------------------------------------------------
    write(log_file_unit,'(a)')      ' HORIZON'
! -------------------------------------------------
#endif 
! endif HORIZON     
    !
#ifdef BOLSHOI
! -------------------------------------------------
    write(log_file_unit,'(a)')      ' BOLSHOI'
! -------------------------------------------------
#endif 
! endif BOLSHOI
    !
#ifdef UNLINKED_TREES
! -------------------------------------------------
    write(log_file_unit,'(a)')      ' UNLINKED_TREES'
! -------------------------------------------------
#endif 
! endif UNLINKED_TREES  
    !    
#ifdef CLEAN
! -------------------------------------------------
    write(log_file_unit,'(a)')      ' CLEAN'
! -------------------------------------------------
#endif 
! endif CLEAN
    !
#ifdef FOLLOW
! -------------------------------------------------
    write(log_file_unit,'(a)')      ' FOLLOW'
! -------------------------------------------------
#endif 
! endif FOLLOW
    !
#ifdef TREE_EVOLVE
! -------------------------------------------------
    write(log_file_unit,'(a)')      ' TREE_EVOLVE'
    !
#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------
    write(log_file_unit,'(a)')      ' LUMINOUS_PROCESSES'
! -------------------------------------------------
#endif 
! endif LUMINOUS_PROCESSES
    !
#ifdef GAL_SPECTRA
! -------------------------------------------------
    write(log_file_unit,'(a)')      ' GAL_SPECTRA'
! -------------------------------------------------
#endif 
! endif GAL_SPECTRA
    !     
#ifdef DAGN_SPECTRUM
! -------------------------------------------------
    write(log_file_unit,'(a)')      ' DAGN_SPECTRUM'
! -------------------------------------------------
#endif 
! endif DAGN_SPECTRUM
    !
#ifdef PHOTOIONISATION
! -------------------------------------------------
    write(log_file_unit,'(a)')      ' PHOTOIONISATION'
    !
#ifdef GNEDIN_2000
! -------------------------------------------------
    write(log_file_unit,'(a)')      ' - GNEDIN_2000'
! -------------------------------------------------
#endif
! endif GNEDIN
    !
#ifdef OKAMOTO_2008
! -------------------------------------------------
    write(log_file_unit,'(a)')      ' - OKAMOTO_2008'
! -------------------------------------------------
#endif
! endif OKAMOTO
    !
! -------------------------------------------------
#endif
! endif PHOTOIONISATION
    !  
#ifdef NO_EVAP
! -------------------------------------------------
    write(log_file_unit,'(a)')      ' NO_EVAP'
! -------------------------------------------------
#endif 
! endif NO_EVAP
    !  
#ifdef COLD_STREAMS
! -------------------------------------------------
    write(log_file_unit,'(a)')      ' COLD_STREAMS'
! -------------------------------------------------
#endif
! endif COLD_STREAMS
    !  
#ifdef REACCRETION
! -------------------------------------------------
    write(log_file_unit,'(a)')      ' REACCRETION'
! -------------------------------------------------
#endif
! endif REACCRETION
    !
#ifdef SN_FEEDBACK_PROP_TO_SFR
! -------------------------------------------------
    write(log_file_unit,'(a)')      ' SN_FEEDBACK_PROP_TO_SFR'
! -------------------------------------------------
#endif 
! endif SN_FEEDBACK_PROP_TO_SFR
    !       
#ifdef SUB_QUENCHING
! -------------------------------------------------
    write(log_file_unit,'(a)')      ' SUB_QUENCHING'
! -------------------------------------------------
#endif
! endif SUB_QUENCHING
!
#ifdef SUB_STRIPPING
! -------------------------------------------------
    write(log_file_unit,'(a)')      ' SUB_STRIPPING'
! -------------------------------------------------    
#endif
! endif SUB_STRIPPING
    !
#ifdef SUB_TRANSFER
! -------------------------------------------------
    write(log_file_unit,'(a)')      ' SUB_TRANSFER'
! -------------------------------------------------    
#endif
! endif SUB_TRANSFER
    !    
#ifdef POLYTROPIC
! -------------------------------------------------
    write(log_file_unit,'(a)')      ' POLYTROPIC'
! -------------------------------------------------    
#endif
! endif POLYTROPIC
    !
#ifdef THERMAL_INSTABILITIES
! -------------------------------------------------
    write(log_file_unit,'(a)')      ' THERMAL_INSTABILITIES'
! -------------------------------------------------    
#endif
! endif THERMAL_INSTABILITIES
    !
    write(log_file_unit,'(a)') '--------------output timestep list---------------'
    write(log_file_unit,'(a,i2.2)') 'n_output_timestep = ', nout
    do i = 1, nout
      write(log_file_unit,'(a,i2.2,a,f6.3)') ' - ', tsout(i), ' | z = ', zout(i)
    end do
    write(log_file_unit,'(a)')         '-----------------Model parameters----------------'
    write(log_file_unit,'(a)')         ' ACCRETION --------------------------------------'
    write(log_file_unit,'(a,e10.4)')   ' cold_stream_efficiency        : ', cold_stream_efficiency
    write(log_file_unit,'(a,e10.4)')   ' cooling_efficiency            : ', cooling_efficiency
#ifdef THERMAL_INSTABILITIES 
! ------------------------------------------------- 
    write(log_file_unit,'(a,e10.4)')   ' TI_efficiency                 : ', TI_efficiency
! ------------------------------------------------- 
#endif      
    write(log_file_unit,'(a)')         ' STARS ------------------------------------------'
    write(log_file_unit,'(a,a)')       ' IMF                           : ', trim(IMF)
    write(log_file_unit,'(a)')         ' GALAXY -----------------------------------------'
    write(log_file_unit,'(a,e10.4)')   ' epsilon_merge                 : ', epsilon_merge
    write(log_file_unit,'(a)')         ' DISC -------------------------------------------'
#ifdef SUB_STRIPPING
! ------------------------------------------------- 
    write(log_file_unit,'(a,e10.4)')   ' disc_stripping_efficiency     : ', disc_stripping_efficiency
! ------------------------------------------------- 
#endif      
    write(log_file_unit,'(a,e10.4)')   ' disc_ejecta_efficiency        : ', disc_ejecta_efficiency    
    write(log_file_unit,'(a)')         ' AGN --------------------------------------------'
    write(log_file_unit,'(a,e10.4)')   ' agn_ejecta_efficiency         : ', agn_ejecta_efficiency
    write(log_file_unit,'(a,e10.4)')   ' M_BH_min [Msun]               : ', M_BH_min*mass_code_unit_in_M_Sun
    write(log_file_unit,'(a,e10.4)')   ' agn_gas_coupling              : ', agn_gas_coupling
    write(log_file_unit,'(a)') ' ------------------------------------------------'
    !
! -------------------------------------------------
#endif
! endif TREE_EVOLVE
    write(log_file_unit,'(a)') ' -------------output parameter order-------------'
    ni = 0 ! init the number of field labels printed
    do i = 1, nb_tree_field
      ni = ni +1
      write(log_file_unit,'(i3.3,a,a,a)') ni, ') ',ttype_tree(i), tunit_tree(i)  
    end do
    write(log_file_unit,'(a)') '  # dm properties -------------------------------'
    do i = 1, nb_dm_field
      ni = ni +1
      write(log_file_unit,'(i3.3,a,a,a)') ni, ') ',ttype_dm(i), tunit_dm(i)  
    end do
#ifdef TREE_EVOLVE
! -------------------------------------------------
    write(log_file_unit,'(a)') ' # baryon halo properties -----------------------'
    do i = 1, nb_bh_field
      ni = ni +1
      write(log_file_unit,'(i3.3,a,a,a)') ni, ') ', ttype_bh(i), tunit_bh(i)  
    end do
    write(log_file_unit,'(a)') ' # galaxy global properties ---------------------'
    do i = 1, nb_gal_field
      ni = ni +1
      write(log_file_unit,'(i3.3,a,a,a)') ni, ') ', ttype_gal(i), tunit_gal(i)  
    end do
    write(log_file_unit,'(a)') ' # disc properties ------------------------------'
    do i = 1, nb_disc_field
      ni = ni +1
      write(log_file_unit,'(i3.3,a,a,a)') ni, ') ', ttype_disc(i), tunit_disc(i)  
    end do
    write(log_file_unit,'(a)') ' # disc stellar population properties -----------'
    do i = 1, nb_stars_field
      ni = ni +1
      write(log_file_unit,'(i3.3,a,a,a)') ni, ') ', ttype_starsd(i), tunit_starsd(i)  
    end do
    write(log_file_unit,'(a)') ' # disc dust properties -------------------------'
    do i = 1, nb_dust_field
      ni = ni +1
      write(log_file_unit,'(i3.3,a,a,a)') ni, ') ', ttype_dust_ISMd(i), tunit_dustd(i)  
    end do
    do i = 1, nb_dust_field
      ni = ni +1
      write(log_file_unit,'(i3.3,a,a,a)') ni, ') ', ttype_dust_BCd(i), tunit_dustd(i)  
    end do
    write(log_file_unit,'(a)') ' # agn properties -------------------------------'
    do i = 1, nb_agn_field
      ni = ni +1
      write(log_file_unit,'(i3.3,a,a,a)') ni, ') ', ttype_agn(i), tunit_agn(i)  
    end do
    write(log_file_unit,'(a)') ' # bulge properties -----------------------------'
    do i = 1, nb_bulge_field
      ni = ni +1
      write(log_file_unit,'(i3.3,a,a,a)') ni, ') ', ttype_bulge(i), tunit_bulge(i)  
    end do
    write(log_file_unit,'(a)') ' # bulge stellar population properties ----------'
    do i = 1, nb_stars_field
      ni = ni +1
      write(log_file_unit,'(i3.3,a,a,a)') ni, ') ', ttype_starsb(i), tunit_starsb(i)  
    end do
    write(log_file_unit,'(a)') ' # bulge dust properties ------------------------'
    do i = 1, nb_dust_field
      ni = ni +1
      write(log_file_unit,'(i3.3,a,a,a)') ni, ') ', ttype_dustb(i), tunit_dustb(i)  
    end do
! -------------------------------------------------    
#endif
! endif TREE_EVOLVE 
    !
#ifdef LUMINOUS_PROCESSES
! -----------------------------------------------
    write(log_file_unit,'(a)') ' '
    write(log_file_unit,'(a)') ' ----------------- filter list ------------------'
    do i = 1, nfilters
        write(log_file_unit,'(i3.3,a,a,a,e9.3,a)') i, ') ', trim(filters_tab(i)%name), '  [', filters_tab(i)%clambda, ' microns]'
    end do
! ------------------------------------------------- 
#endif
! LUMINOUS_PROCESSES    
    
    write(log_file_unit,'(a)') ''
    write(log_file_unit,'(a)') '*************************************************'
    write(log_file_unit,'(a)') '*                    END OF                     *'
    write(log_file_unit,'(a)') '*               G.A.S. log file                 *'
    write(log_file_unit,'(a)') '*************************************************'
    !
    close(log_file_unit)
    !  
  end subroutine IO_print_log_file 

  !*****************************************************************************************************************

  subroutine IO_create_open_and_print_in_main_FITS_output_files

    ! CREATE, OPEN AND PRINT IN MAIN FITS OUTPUT FILE
    ! All results are saved in various tmp output files (one per MPI process and per ts) 
    ! This routine create, open and print in data in the main output FITS file
     
    implicit none

    integer(kind=4)             :: nfiles               ! number of tmp output files
    integer(kind=4)             :: its                  ! loop index
    integer(kind=4)             :: i,f                  ! loop indexes
    integer(kind=4)             :: n                    ! loop index (1 --> n_blocks)
    integer(kind=4)             :: ih                   ! loop index (halos)
    integer(kind=4)             :: main_phy_FITS_file_unit
    integer(kind=4)             :: n_blocks             ! number of data blocks in a tmp output file
    integer(kind=4)             :: ts                   ! the timestep index (dark matter simulation reference)
    integer(kind=4)             :: ots                  ! the timestep index 
    integer(kind=4)             :: u,u2                 ! tmp filelist unit
    integer(kind=4)             :: nb_halos             ! number of halos in a data block
    integer(kind=4)             :: nrows                ! numbers of rows in the current FITS extension
    integer(kind=4)             :: nrows_in_buffer
    integer(kind=4)             :: fits_nrows
    integer(kind=4)             :: slash
    integer(kind=4)             :: nb_current
    integer(kind=4),parameter   :: buffer_size = 50000

    character(MAXPATHSIZE)      :: new_file
    character(MAXPATHSIZE)      :: filename
    character(MAXPATHSIZE)      :: adress
    character(MAXPATHSIZE)      :: extname              ! extension name
    character(MAXPATHSIZE)      :: message              ! a message

    real(kind=8),allocatable    :: buff_data(:,:)       ! buffer for physical data
    
    ! Additionnal variables linked to other compilation options ///////////////////////
        
#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------
    integer(kind=4)             :: j                    ! loop indexes
    integer(kind=4)             :: hdutype
    integer(kind=4)             :: main_lum_FITS_file_unit
#ifdef GAL_SPECTRA
! -------------------------------------------------      
    real(kind=4),allocatable    :: buff_spec1(:,:)      ! buffer for the first spectrum
    real(kind=4),allocatable    :: buff_spec2(:,:)      ! buffer for the second spectrum
    real(kind=4),allocatable    :: buff_spec3(:,:)      ! buffer for the third spectrum
! -------------------------------------------------
#endif
!     
    real(kind=4),allocatable    :: buff_mags(:,:)       ! buffer for magnitudes
! ------------------------------------------------- 
#endif
! LUMINOUS_PROCESSES  
    
    ! /////////////////////////////////////////////////////////////////////////////////
    !
    call IO_print_message('Create open and print in main FITS output files')
    !
    ! read temporary output files and convert data
    !
    ! PHYSICAL DATA
    !
    if (physical_process) then 
        !
        do f = 1, nb_main_output_files_to_print
            !
            ! Open tmp filename list
            write(filename, '(a,a,i3.3,a)') trim(output_path), '/GAS_phy_tmp_filename_list_ts_', list_of_ts_to_print(f), '.tmp'
            u = IO_get_tmp_filelist_unit(list_of_ts_to_print(f),rank)
            open(unit = u, file = filename, status = 'old')
            !
            ! Each MPI physical (or luminous) process create and used nout tmp output files (one per output timestep)
            ! read the number of file listed after
            read(u,'(i3.3)') nfiles  
            write(message,'(a,i3.3,a)') 'RE-format output data from ', nfiles, ' phy_tmp files'
            call IO_print_message(trim(message),only_rank=rank)
            !                              
            ! init table
            nrows = 0              ! init                                     
            nb_halos_computed = 0  ! init
            !
            ! Main loop under tmp output files                                  
            do its = 1, nfiles
              ! read the tmp output filename
              read(u,'(a)') adress
              slash = scan(adress, '/', back = .TRUE.)
              write(filename, '(a,a)') trim(output_path), trim(adjustl(adress(slash:)))
              !
              write(message,'(a,a)') 'Load : ', trim(filename)
              call IO_print_message(trim(message),only_rank=rank)
              !
              ! open the tmp output file
              u2 = IO_get_temporary_output_files_unit(list_of_ts_to_print(f))
              open(unit = u2, file = trim(filename), status = 'old',form = 'unformatted')
              ! read : - nb_files_computed corresponding to the number of blocks print in this tmp output file,
              !        - ts, the timestep value of the N-body simulation used 
              !        - ots (index over the output timestep list)
              read(u2) n_blocks, ts, ots
              !
              if (n_blocks .gt. 0) then ! the file contains n_blocks of data
                ! the temporary files are read in ts order ts = 1 --> ts = nts
                do n = 1, n_blocks
                  ! read number of galaxies in the following block 
                  read(u2) nb_halos
                  if (nb_halos .gt. 0) then
                    if (nrows .eq. 0) then
                        ! build the filename of the new file
                        write(new_file, '(a,i3.3,a)') '/GAS_phy_ts_',ts,'.fits'
                        write(filename, '(a,a)') trim(output_path), trim(new_file)
                        !
                        ! Open the FITS file
                        status = 0  
                        ! create the file with unit "main_phy_FITS_file_unit" 
                        call ftgiou(main_phy_FITS_file_unit,status)
                        ! check the status     
                        if (status .gt. 0) then  
                            call IO_print_error_message('Cannot create (physical properties) FITS output file',only_rank=rank, &
                                    called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                            write(message,'(a,i4,a,a)') 'ftgiou status: ', status, ' for file: ', trim(new_file)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        ! open the new FITS file
                        call ftinit(main_phy_FITS_file_unit,trim(filename),blocksize,status)  
                        ! check the status
                        if (status .gt. 0) then  
                            call IO_print_error_message('Cannot open (physical properties) FITS output file',only_rank=rank, &
                                        called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                            write(message,'(a,i4,a,a)') 'ftinit status: ', status, ' for file: ', trim(new_file)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        ! create a new extension with name ts_i
                        write(extname, '(a,i3.3)') 'ts_', ts
                        call ftibin(main_phy_FITS_file_unit,0,nb_field,ttype,tform,tunit,extname,varidat,status)  ! here nrowsll = 0
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot create the main extension for (physical properties) FITS output file', &
                                    only_rank=rank, called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                            write(message,'(a,i4,a,a)') 'ftibin status: ', status, ' for file: ', trim(filename)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        ! @ this point the main extension is created
                        call ftpkye(main_phy_FITS_file_unit, 'z', zout(ots), 3,'redshift', status)
                        call ftpkyj(main_phy_FITS_file_unit, 'RES', res, 'simulation resolution (RES^3 particules)', status)
                        write(message,'(a,a,i3.3,a,f6.3)') trim(new_file), ' ts : ', ts, ' / z = ', zout(ots) 
                        call IO_print_message(message,only_rank=rank)
                        !
                        ! allocate buffer
                        allocate(buff_data(nb_field,buffer_size))
                    end if
                    !
                    nb_halos_computed = nb_halos_computed + nb_halos 
                    nrows_in_buffer = 0   ! init 
                    do ih = 1, nb_halos
                      if ((nrows_in_buffer .lt. buffer_size) .and. (ih .le. nb_halos)) then    ! load buffer
                        nrows_in_buffer = nrows_in_buffer +1
                        if (BARYON) then
                          ! read properties
                          ! tree properties *****************************************
                          read(u2)  (buff_data(i,nrows_in_buffer),i=1,nb_tree_field)
                          nb_current = nb_tree_field
                          ! dm properties *******************************************
                          read(u2)  (buff_data(i,nrows_in_buffer),i=nb_current +1,nb_current+nb_dm_field)
                          nb_current = nb_current + nb_dm_field
                          ! baryon halo properties **********************************
                          read(u2)  (buff_data(i,nrows_in_buffer),i=nb_current +1,nb_current+nb_bh_field)
                          nb_current = nb_current +nb_bh_field
                          ! galaxy global properties ********************************
                          read(u2)  (buff_data(i,nrows_in_buffer),i=nb_current +1,nb_current+nb_gal_field)
                          nb_current = nb_current +nb_gal_field
                          ! disc properties *****************************************
                          read(u2)  (buff_data(i,nrows_in_buffer),i=nb_current +1, &
                            nb_current+nb_disc_field+nb_stars_field+2*nb_dust_field+nb_agn_field)
                          nb_current = nb_current +nb_disc_field +nb_stars_field +2*nb_dust_field+nb_agn_field
                          ! bulge properties ****************************************
                          read(u2)  (buff_data(i,nrows_in_buffer),i=nb_current +1, &
                            nb_current+nb_bulge_field+nb_stars_field+nb_dust_field)
                          nb_current = nb_current +nb_bulge_field +nb_stars_field +nb_dust_field
                        else                 
                          ! read only tree *******************************************
                          read(u2)  (buff_data(i,nrows_in_buffer),i=1,nb_tree_field)
                          nb_current = nb_tree_field
                          ! and dm properties ****************************************
                          read(u2)  (buff_data(i,nrows_in_buffer),i=nb_current +1,nb_current+nb_dm_field)
                          nb_current = nb_current +nb_dm_field
                        end if ! end if BARYON
                      end if
                      !
                      if ((nrows_in_buffer .eq. buffer_size) .or. (ih .eq. nb_halos)) then  
                        ! full buffer or end of the galaxy list, print into the fits file
                        ! create nrows_in_buffer new rows
                        call ftirow(main_phy_FITS_file_unit,0,nrows_in_buffer,status)
                        ! print data
                        do i = 1, nb_field
                          call ftpcld(main_phy_FITS_file_unit,i,1,1,nrows_in_buffer,buff_data(i,1:nrows_in_buffer),status)    
                        end do
                        nrows           = nrows + nrows_in_buffer
                        nrows_in_buffer = 0      ! reset
                        buff_data(:,:)  = 0.d0   ! reset
                      end if
                    end do ! end ih loop  
                    if (nrows .ne. nb_halos_computed) then
                      call IO_print_error_message('nrows != nb_halos_computed [phy properties]',called_by='IO_create_open_and_print_in_main_FITS_output_files')
                      stop
                    end if
                    call ftgnrw(main_phy_FITS_file_unit,fits_nrows,status)
                    if (fits_nrows .ne. nrows) then
                      call ftdrow(main_phy_FITS_file_unit,nrows+1,fits_nrows-nrows,status)
                    end if
                  end if
                end do ! end n_blocks loop 
              end if ! n_blocks
              !
              close(u2)
              !
            end do ! end its loop
            !
            ! dellocate buffer
            if (allocated(buff_data)) deallocate(buff_data) 
            !
            ! close tmp filelist file
            close(u)   
            !
            ! Close FITS file   
            call ftclos(main_phy_FITS_file_unit, status)
            ! check the status     
            if (status .gt. 0) then  
                call IO_print_error_message('Cannot close (physical properties) FITS output file',only_rank=rank, &
                    called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                write(message,'(a,i4,a,a)') 'ftclos status: ', status, ' for file: ', trim(new_file)
                call IO_print_message(message,only_rank=rank)
                stop   ! stop the program
            end if
            call ftfiou(main_phy_FITS_file_unit, status)
            ! check the status     
            if (status .gt. 0) then  
                call IO_print_error_message('Cannot deallocate unit index for (physical properties) FITS output file',only_rank=rank, &
                    called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                write(message,'(a,i4,a,a)') 'ftfiou status: ', status, ' for file: ', trim(new_file)
                call IO_print_message(message,only_rank=rank)
                stop   ! stop the program
            end if
        end do
    end if
    !
#ifdef TREE_EVOLVE
! -------------------------------------------------  
    !
#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------  
    !
    ! LUMINOUS DATA
    !
    if (luminous_process) then 
        !
        do f = 1, nb_main_output_files_to_print
            !
            ! Open tmp filename list
            write(filename, '(a,a,i3.3,a)') trim(output_path),'/GAS_lum_tmp_filename_list_ts_',list_of_ts_to_print(f),'.tmp'
            u = IO_get_tmp_filelist_unit(list_of_ts_to_print(f),rank)
            open(unit = u, file = filename, status = 'old')
            !
            ! Each MPI physical (or luminous) process create and used nout tmp output files (one per output timestep)
            ! read the number of file listed after
            read(u,'(i3.3)') nfiles  
            write(message,'(a,i3.3,a)') 'RE-format output data from ', nfiles, ' lum_tmp files'
            call IO_print_message(trim(message),only_rank=rank)
            !
            ! init table 
            nrows             = 0          ! init
            nb_halos_computed = 0          ! init
            !
            ! Main loop under tmp output files                                  
            do its = 1, nfiles
              ! read the tmp output filename
              read(u,'(a)') adress
              slash = scan(adress, '/', back = .TRUE.)
              write(filename, '(a,a)') trim(output_path), trim(adjustl(adress(slash:)))
              !
              write(message,'(a,a)') 'Load : ', trim(filename)
              call IO_print_message(trim(message),only_rank=rank)
              !
              ! open the tmp output file
              u2 = IO_get_temporary_output_files_unit(list_of_ts_to_print(f))
              open(u2, file = trim(filename), status = 'old', form = 'unformatted')
              ! read : - nb_files_computed corresponding to the number of blocks print in this tmp output file,
              !        - ts, the timestep value of the N-body simulation used 
              !        - ots (index over the output timestep list)
              read(u2) n_blocks, ts, ots
              !
              if (n_blocks .gt. 0) then ! the file contains n_blocks of data
                ! the temporary files are read in ts order ts = 1 --> ts = nts we create a fits extension for each ts 
                do n = 1, n_blocks
                  ! read number of galaxies in the following block 
                  read(u2) nb_halos
                  if (nb_halos .gt. 0) then
                      if (nrows .eq. 0) then
                        ! build the filename of the new file
                        write(new_file, '(a,i3.3,a)') '/GAS_lum_ts_',ts,'.fits'
                        write(filename, '(a,a)') trim(output_path), trim(new_file)
                        !
                        ! Open the FITS file
                        status = 0  
                        ! create the file with unit "main_lum_FITS_file_unit" 
                        call ftgiou(main_lum_FITS_file_unit,status)
                        ! check the status     
                        if (status .gt. 0) then  
                            call IO_print_error_message('Cannot create (luminous properties) FITS output file', &
                                    only_rank=rank, called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                            write(message,'(a,i4,a,a)') ' --> ftgiou status: ', status, ' for file: ', trim(new_file)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        ! open the new FITS file
                        call ftinit(main_lum_FITS_file_unit,trim(filename),blocksize,status)  
                        ! check the status
                        if (status .gt. 0) then  
                            call IO_print_error_message('Cannot open (luminous properties) FITS output file', &
                                        called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                            write(message,'(a,i4,a,a)') ' --> ftinit status: ', status, ' for file: ', trim(new_file)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        ! create extension for tree structure 
                        call ftibin(main_lum_FITS_file_unit,0,nb_tree_field,ttype_tree,tform_tree,tunit_tree,'tree',varidat,status) ! here nrowsll is set to 0.
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot create tree extension',only_rank=rank, &
                                    called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                            write(message,'(a,i4,a,a)') 'ftibin status: ', status, ' for file: ', trim(filename)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        ! @ this point the main (tree) extension is created
                        call ftpkye(main_lum_FITS_file_unit, 'z', zout(ots), 3,'redshift', status)
                        call ftpkyj(main_lum_FITS_file_unit, 'RES', res, 'simulation resolution (RES^3 particules)', status)
                        ! attribute hdu parameter: "hdu_tree"
                        call ftghdn(main_lum_FITS_file_unit,hdu_tree) 
                        !
                        ! create extension for NON extincted observer-frame magnitudes
                        call ftibin(main_lum_FITS_file_unit,0,nfilters,ttype_mag,tform_mag,tunit_mag,'Non_Ext_OBS_Mags',varidat,status) ! here nrowsll is set to 0.
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot create extension for Non extincted OBS AB magnitudes',only_rank=rank, &
                                        called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                            write(message,'(a,i4,a,a)') 'ftibin status: ', status, ' for file: ', trim(filename)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        ! @ this point the extension dedicated to NON extincted observer-frame magnitudes
                        ! attribute hdu parameter: "hdu_NE_Mag_obs"
                        call ftghdn(main_lum_FITS_file_unit,hdu_NE_Mag_obs) 
                        !
                        ! create extension for NON extinced rest-frame magnitudes
                        call ftibin(main_lum_FITS_file_unit,0,nfilters,ttype_mag,tform_mag,tunit_mag,'Non_Ext_GAL_Mags',varidat,status) ! here nrowsll is set to 0.
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot create extension for Non extincted GAL AB magnitudes', &
                                    only_rank=rank, called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                            write(message,'(a,i4,a,a)') 'ftibin status: ', status, ' for file: ', trim(filename)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        ! @ this point the extension dedicated to NON extinced rest-frame magnitudes
                        ! attribute hdu parameter: "hdu_NE_Mag_gal"
                        call ftghdn(main_lum_FITS_file_unit,hdu_NE_Mag_gal) 
                        !
                        ! create extension for observer-frame AB magnitudes
                        call ftibin(main_lum_FITS_file_unit,0,nfilters,ttype_mag,tform_mag,tunit_mag,'OBS_Mags',varidat,status) ! here nrowsll is set to 0.
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot create extension for OBS AB magnitudes',only_rank=rank, &
                                        called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                            write(message,'(a,i4,a,a)') 'ftibin status: ', status, ' for file: ', trim(filename)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        ! @ this point the extension dedicated to OBS AB magnitudes is created
                        ! attribute hdu parameter: "hdu_Mag_obs"
                        call ftghdn(main_lum_FITS_file_unit,hdu_Mag_obs) 
                        !
                        ! create extension for observer-frame AB magnitudes (for disc only)
                        call ftibin(main_lum_FITS_file_unit,0,nfilters,ttype_mag,tform_mag,tunit_mag,'OBS_Mags_disc',varidat,status) ! here nrowsll is set to 0.
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot create extension for OBS AB disc magnitudes', &
                                        only_rank=rank, called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                            write(message,'(a,i4,a,a)') 'ftibin status: ', status, ' for file: ', trim(filename)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        ! @ this point the extension dedicated to OBS AB magnitudes is created
                        ! attribute hdu parameter: "hdu_Mag_disc_obs"
                        call ftghdn(main_lum_FITS_file_unit,hdu_Mag_disc_obs) 
                        !
                        ! create extension for rest-frame AB magnitudes
                        call ftibin(main_lum_FITS_file_unit,0,nfilters,ttype_mag,tform_mag,tunit_mag,'GAL_Mags',varidat,status) ! here nrowsll is set to 0.
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot create extension for GAL AB magnitudes',only_rank=rank, &
                                        called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                            write(message,'(a,i4,a,a)') 'ftibin status: ', status, ' for file: ', trim(filename)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        ! @ this point the extension dedicated to GAL AB magnitudes is created
                        ! attribute hdu parameter: "hdu_Mag_gal"
                        call ftghdn(main_lum_FITS_file_unit,hdu_Mag_gal) 
                        !
                        ! create extension for observer-frame AB first order derivative magnitudes
                        call ftibin(main_lum_FITS_file_unit,0,nfilters,ttype_mag,tform_mag,tunit_mag,'OBS_dMags',varidat,status) ! here nrowsll is set to 0.
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot create extension for OBS first order derivative magnitudes', &
                                        only_rank=rank,called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                            write(message,'(a,i4,a,a)') 'ftibin status: ', status, ' for file: ', trim(filename)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        ! @ this point the extension dedicated to observer-frame AB first order derivative magnitudes
                        ! attribute hdu parameter: "hdu_dMag_obs"
                        call ftghdn(main_lum_FITS_file_unit,hdu_dMag_obs) 
                        !
                        ! create extension for galaxy luminous properties
                        call ftibin(main_lum_FITS_file_unit,0,nlumprops,ttype_lumprops,tform_lumprops,tunit_lumprops,'GAL_Props',varidat,status) ! here nrowsll is set to 0.
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot create extension for GAL luminous properties',only_rank=rank, &
                                        called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                            write(message,'(a,i4,a,a)') 'ftibin status: ', status, ' for file: ', trim(filename)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        ! @ this point the extension dedicated to luminous properties is created
                        ! attribute hdu parameter: "hdu_Props"
                        call ftghdn(main_lum_FITS_file_unit,hdu_Props) 
                        !
#ifdef GAL_SPECTRA
! -------------------------------------------------                
                        ! create extension for young star spectra
                        call ftibin(main_lum_FITS_file_unit,0,nb_spec_field,ttype_spec,tform_spec,tunit_spec,'young_stars',varidat,status) ! here nrowsll is set to 0.
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot create extension for young stars spectra',only_rank=rank, &
                                        called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                            write(message,'(a,i4,a,a)') 'ftibin status: ', status, ' for file: ', trim(filename)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        ! @ this point the extension dedicated to young stars spectrum is created
                        ! attribute hdu parameter: "hdu_ystars"
                        call ftghdn(main_lum_FITS_file_unit,hdu_ystars) 
                        !
                        ! create extension for old star spectra
                        call ftibin(main_lum_FITS_file_unit,0,nb_spec_field,ttype_spec,tform_spec,tunit_spec,'old_stars',varidat,status) ! here nrowsll is set to 0.
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot create extension for old stars spectra',only_rank=rank, &
                                        called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                            write(message,'(a,i4,a,a)') 'ftibin status: ', status, ' for file: ', trim(filename)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        ! @ this point the extension dedicated to young stars spectrum is created
                        ! attribute hdu parameter: "hdu_ostars"
                        call ftghdn(main_lum_FITS_file_unit,hdu_ostars) 
                        !
                        ! create extension for full galaxy spectrum
                        call ftibin(main_lum_FITS_file_unit,0,nb_spec_field,ttype_spec,tform_spec,tunit_spec,'full_gal_spectrum',varidat,status) ! here nrowsll is set to 0.
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot create extension for full galaxy spectrum',only_rank=rank, &
                                        called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                            write(message,'(a,i4,a,a)') 'ftibin status: ', status, ' for file: ', trim(filename)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        ! @ this point the last extension dedicated to the full galaxy spectrum is created
                        ! attribute hdu parameter: "hdu_fullspec"
                        call ftghdn(main_lum_FITS_file_unit,hdu_fullspec) 
! -------------------------------------------------- 
#endif
! GAL_SPECTRA            
                        ! print information about the new file
                        write(message,'(a,a,i3.3,a,f6.3)') trim(new_file), ' ts : ', ts, ' / z = ', zout(ots) 
                        call IO_print_message(message,only_rank=rank)
                        !
                        ! allocate buffer
                        allocate(buff_data(nb_tree_field,buffer_size))
                        allocate(buff_mags(6*nfilters+nlumprops,buffer_size))
#ifdef GAL_SPECTRA
! -------------------------------------------------         
                        allocate(buff_spec1(nWaves,buffer_size))
                        allocate(buff_spec2(nWaves,buffer_size))
                        allocate(buff_spec3(nWaves,buffer_size))
! ------------------------------------------------- 
#endif
! GAL_SPECTRA
                    end if
                    !
                    nb_halos_computed = nb_halos_computed + nb_halos 
                    nrows_in_buffer = 0   ! init 
                    do ih = 1, nb_halos
                      if ((nrows_in_buffer .lt. buffer_size) .and. (ih .le. nb_halos)) then    ! load buffer
                        nrows_in_buffer = nrows_in_buffer +1
                        if (BARYON) then
                          ! read properties
                          ! tree properties *****************************************
                          read(u2) (buff_data(i,nrows_in_buffer),i=1,nb_tree_field)
                          read(u2) (buff_mags(i,nrows_in_buffer),i=1,6*nfilters+nlumprops)
#ifdef GAL_SPECTRA
! -------------------------------------------------                   
                          ! young stars spectrum ************************************
                          read(u2) (buff_spec1(i,nrows_in_buffer),i=1,nWaves)
                          ! old stars spectrum **************************************
                          read(u2) (buff_spec2(i,nrows_in_buffer),i=1,nWaves)
                          ! full galaxy spectrum ************************************
                          read(u2) (buff_spec3(i,nrows_in_buffer),i=1,nWaves)
! ------------------------------------------------- 
#endif
! GAL_SPECTRA
                        end if ! end if BARYON
                      end if
                      !
                      if ((nrows_in_buffer .eq. buffer_size) .or. (ih .eq. nb_halos)) then  
                        ! full buffer or end of the galaxy list, print into the fits file
                        !
                        ! move to tree extension
                        call ftmahd(main_lum_FITS_file_unit,hdu_tree,hdutype,status) 
                        ! create nrows_in_buffer new rows
                        call ftirow(main_lum_FITS_file_unit,0,nrows_in_buffer,status)
                        ! print data
                        do i = 1, nb_tree_field
                          call ftpcld(main_lum_FITS_file_unit,i,1,1,nrows_in_buffer,buff_data(i,1:nrows_in_buffer),status)    
                        end do
                        !
                        ! move to Non Extincted OBS Magnitudes extension
                        call ftmahd(main_lum_FITS_file_unit,hdu_NE_Mag_obs,hdutype,status) 
                        ! create nrows_in_buffer new rows
                        call ftirow(main_lum_FITS_file_unit,0,nrows_in_buffer,status)
                        ! print data
                        do i = 1, nfilters
                            call ftpcle(main_lum_FITS_file_unit,i,1,1,nrows_in_buffer,buff_mags(i,1:nrows_in_buffer),status)  
                        end do  
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot print new set of Non Extincetd OBS Mags',only_rank=rank, &
                                    called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                            write(message,'(a,i4,a,a)') 'ftpcle status: ', status, ' for file: ', trim(filename)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        !
                        ! move to Non Extincted GAL Magnitudes
                        call ftmahd(main_lum_FITS_file_unit,hdu_NE_Mag_gal,hdutype,status) 
                        ! create nrows_in_buffer new rows
                        call ftirow(main_lum_FITS_file_unit,0,nrows_in_buffer,status)
                        ! print data
                        do i = 1, nfilters
                            j = i + nfilters
                            call ftpcle(main_lum_FITS_file_unit,i,1,1,nrows_in_buffer,buff_mags(j,1:nrows_in_buffer),status)  
                        end do  
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot print new set of Non Extincted Gal Mags',only_rank=rank, &
                                    called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                            write(message,'(a,i4,a,a)') 'ftpcle status: ', status, ' for file: ', trim(filename)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        !
                        ! move to OBS Magnitudes extension
                        call ftmahd(main_lum_FITS_file_unit,hdu_Mag_obs,hdutype,status) 
                        ! create nrows_in_buffer new rows
                        call ftirow(main_lum_FITS_file_unit,0,nrows_in_buffer,status)
                        ! print data
                        do i = 1, nfilters
                            j = i + 2*nfilters
                            call ftpcle(main_lum_FITS_file_unit,i,1,1,nrows_in_buffer,buff_mags(j,1:nrows_in_buffer),status)  
                        end do  
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot print new set of OBS Mags',only_rank=rank, &
                                    called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                            write(message,'(a,i4,a,a)') 'ftpcle status: ', status, ' for file: ', trim(filename)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        !
                        ! move to OBS Magnitudes disc only extension
                        call ftmahd(main_lum_FITS_file_unit,hdu_Mag_disc_obs,hdutype,status) 
                        ! create nrows_in_buffer new rows
                        call ftirow(main_lum_FITS_file_unit,0,nrows_in_buffer,status)
                        ! print data
                        do i = 1, nfilters
                            j = i + 3*nfilters
                            call ftpcle(main_lum_FITS_file_unit,i,1,1,nrows_in_buffer,buff_mags(j,1:nrows_in_buffer),status)  
                        end do  
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot print new set of OBS Mags',only_rank=rank, &
                                    called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                            write(message,'(a,i4,a,a)') 'ftpcle status: ', status, ' for file: ', trim(filename)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        !
                        ! move to GAL Magnitudes
                        call ftmahd(main_lum_FITS_file_unit,hdu_Mag_gal,hdutype,status) 
                        ! create nrows_in_buffer new rows
                        call ftirow(main_lum_FITS_file_unit,0,nrows_in_buffer,status)
                        ! print data
                        do i = 1, nfilters
                            j = i + 4*nfilters
                            call ftpcle(main_lum_FITS_file_unit,i,1,1,nrows_in_buffer,buff_mags(j,1:nrows_in_buffer),status)  
                        end do  
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot print new set of Gal Mags',only_rank=rank, &
                                    called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                            write(message,'(a,i4,a,a)') 'ftpcle status: ', status, ' for file: ', trim(filename)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        !
                        ! move to OBS first derivative magnitudes 
                        call ftmahd(main_lum_FITS_file_unit,hdu_dMag_obs,hdutype,status) 
                        ! create nrows_in_buffer new rows
                        call ftirow(main_lum_FITS_file_unit,0,nrows_in_buffer,status)
                        ! print data
                        do i = 1, nfilters
                            j = i + 5*nfilters
                            call ftpcle(main_lum_FITS_file_unit,i,1,1,nrows_in_buffer,buff_mags(j,1:nrows_in_buffer),status)  
                        end do  
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot print new set of OBS dMags',only_rank=rank, &
                                    called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                            write(message,'(a,i4,a,a)') 'ftpcle status: ', status, ' for file: ', trim(filename)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        !
                        ! move to Luminous properties
                        call ftmahd(main_lum_FITS_file_unit,hdu_Props,hdutype,status) 
                        ! create nrows_in_buffer new rows
                        call ftirow(main_lum_FITS_file_unit,0,nrows_in_buffer,status)
                        ! print data
                        do i = 1, nlumprops
                            j = i + 6*nfilters
                            call ftpcle(main_lum_FITS_file_unit,i,1,1,nrows_in_buffer,buff_mags(j,1:nrows_in_buffer),status)  
                        end do  
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot print new set of Luminous Properties',only_rank=rank, &
                                    called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                            write(message,'(a,i4,a,a)') 'ftpcle status: ', status, ' for file: ', trim(filename)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        !
#ifdef GAL_SPECTRA
! -------------------------------------------------                 
                        ! move to young stars spectrum extension
                        call ftmahd(main_lum_FITS_file_unit,hdu_ystars,hdutype,status) 
                        ! create nrows_in_buffer new rows
                        call ftirow(main_lum_FITS_file_unit,0,nrows_in_buffer,status)
                        ! print data
                        call ftpcle(main_lum_FITS_file_unit,1,1,1,nrows_in_buffer*nWaves,buff_spec1(:,1:nrows_in_buffer),status)    
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot print new young stars spectra',only_rank=rank, &
                                    called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                            write(message,'(a,i4,a,a)') 'ftpcle status: ', status, ' for file: ', trim(filename)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if  
                        !
                        ! move to old stars spectrum extension
                        call ftmahd(main_lum_FITS_file_unit,hdu_ostars,hdutype,status) 
                        ! create nrows_in_buffer new rows
                        call ftirow(main_lum_FITS_file_unit,0,nrows_in_buffer,status)
                        ! print data
                        call ftpcle(main_lum_FITS_file_unit,1,1,1,nrows_in_buffer*nWaves,buff_spec2(:,1:nrows_in_buffer),status)    
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot print new old stars spectra',only_rank=rank, &
                                    called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                            write(message,'(a,i4,a,a)') 'ftpcle status: ', status, ' for file: ', trim(filename)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
                        !
                        ! move to full galaxy spectum extension
                        call ftmahd(main_lum_FITS_file_unit,hdu_fullspec,hdutype,status) 
                        ! create nrows_in_buffer new rows
                        call ftirow(main_lum_FITS_file_unit,0,nrows_in_buffer,status)
                        ! print data
                        call ftpcle(main_lum_FITS_file_unit,1,1,1,nrows_in_buffer*nWaves,buff_spec3(:,1:nrows_in_buffer),status)    
                        ! check the status
                        if (status .gt. 0) then
                            call IO_print_error_message('Cannot print new full galaxy spectra',only_rank=rank, &
                                    called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                            write(message,'(a,i4,a,a)') 'ftpcle status: ', status, ' for file: ', trim(filename)
                            call IO_print_message(message,only_rank=rank)
                            stop   ! stop the program
                        end if
! ------------------------------------------------- 
#endif
! GAL_SPECTRA    
                        nrows           = nrows + nrows_in_buffer
                        nrows_in_buffer = 0    ! reset
                        buff_data(:,:)  = 0.d0 ! reset
                        buff_mags(:,:)  = 0.d0 ! reset
#ifdef GAL_SPECTRA
! -------------------------------------------------                  
                        buff_spec1(:,:) = 0.   ! reset
                        buff_spec2(:,:) = 0.   ! reset
                        buff_spec3(:,:) = 0.   ! reset
! ------------------------------------------------- 
#endif
! GAL_SPECTRA                 
                      end if
                      !
                    end do ! end ih loop  
                    !
                    if (nrows .ne. nb_halos_computed) then
                      call IO_print_error_message('nrows != nb_halos_computed(ots) [lum properties]', &
                                only_rank=rank, called_by='IO_create_open_and_print_in_main_FITS_output_files') 
                      stop
                    end if
                    !
                    ! move to tree extension
                    call ftmahd(main_lum_FITS_file_unit,hdu_tree,hdutype,status)
                    ! extract the number of rows 
                    call ftgnrw(main_lum_FITS_file_unit,fits_nrows,status)
                    if (fits_nrows .ne. nrows) then
                      call ftdrow(main_lum_FITS_file_unit,nrows+1,fits_nrows-nrows,status)
                    end if
                    !
                    ! move to Non extincted OBS Mags extension
                    call ftmahd(main_lum_FITS_file_unit,hdu_NE_Mag_obs,hdutype,status)
                    ! extract the number of rows 
                    call ftgnrw(main_lum_FITS_file_unit,fits_nrows,status)
                    if (fits_nrows .ne. nrows) then
                      call ftdrow(main_lum_FITS_file_unit,nrows+1,fits_nrows-nrows,status)
                    end if
                    !
                    ! move to Non extincted GAL Mags extension
                    call ftmahd(main_lum_FITS_file_unit,hdu_NE_Mag_gal,hdutype,status)
                    ! extract the number of rows 
                    call ftgnrw(main_lum_FITS_file_unit,fits_nrows,status)
                    if (fits_nrows .ne. nrows) then
                      call ftdrow(main_lum_FITS_file_unit,nrows+1,fits_nrows-nrows,status)
                    end if
                    !
                    ! move to OBS Mags extension
                    call ftmahd(main_lum_FITS_file_unit,hdu_Mag_obs,hdutype,status)
                    ! extract the number of rows 
                    call ftgnrw(main_lum_FITS_file_unit,fits_nrows,status)
                    if (fits_nrows .ne. nrows) then
                      call ftdrow(main_lum_FITS_file_unit,nrows+1,fits_nrows-nrows,status)
                    end if
                    !
                    ! move to OBS Mags (disc only) extension
                    call ftmahd(main_lum_FITS_file_unit,hdu_Mag_disc_obs,hdutype,status)
                    ! extract the number of rows 
                    call ftgnrw(main_lum_FITS_file_unit,fits_nrows,status)
                    if (fits_nrows .ne. nrows) then
                      call ftdrow(main_lum_FITS_file_unit,nrows+1,fits_nrows-nrows,status)
                    end if
                    !
                    ! move to GAL Mags extension
                    call ftmahd(main_lum_FITS_file_unit,hdu_Mag_gal,hdutype,status)
                    ! extract the number of rows 
                    call ftgnrw(main_lum_FITS_file_unit,fits_nrows,status)
                    if (fits_nrows .ne. nrows) then
                      call ftdrow(main_lum_FITS_file_unit,nrows+1,fits_nrows-nrows,status)
                    end if
                    !
                    ! move to OBS dMags extension
                    call ftmahd(main_lum_FITS_file_unit,hdu_dMag_obs,hdutype,status)
                    ! extract the number of rows 
                    call ftgnrw(main_lum_FITS_file_unit,fits_nrows,status)
                    if (fits_nrows .ne. nrows) then
                      call ftdrow(main_lum_FITS_file_unit,nrows+1,fits_nrows-nrows,status)
                    end if
                    !
                    ! move to Luminous properties extension
                    call ftmahd(main_lum_FITS_file_unit,hdu_Props,hdutype,status)
                    ! extract the number of rows 
                    call ftgnrw(main_lum_FITS_file_unit,fits_nrows,status)
                    if (fits_nrows .ne. nrows) then
                      call ftdrow(main_lum_FITS_file_unit,nrows+1,fits_nrows-nrows,status)
                    end if
                    !
#ifdef GAL_SPECTRA
! -------------------------------------------------             
                    ! move to young stars spectra extension
                    call ftmahd(main_lum_FITS_file_unit,hdu_ystars,hdutype,status)
                    ! extract the number of rows 
                    call ftgnrw(main_lum_FITS_file_unit,fits_nrows,status)
                    if (fits_nrows .ne. nrows) then
                      call ftdrow(main_lum_FITS_file_unit,nrows+1,fits_nrows-nrows,status)
                    end if
                    !
                    ! move to old stras spectra extension
                    call ftmahd(main_lum_FITS_file_unit,hdu_ostars,hdutype,status)
                    ! extract the number of rows 
                    call ftgnrw(main_lum_FITS_file_unit,fits_nrows,status)
                    if (fits_nrows .ne. nrows) then
                      call ftdrow(main_lum_FITS_file_unit,nrows+1,fits_nrows-nrows,status)
                    end if
                    !
                    ! move to full galaxy spectra extension
                    call ftmahd(main_lum_FITS_file_unit,hdu_fullspec,hdutype,status)
                    ! extract the number of rows 
                    call ftgnrw(main_lum_FITS_file_unit,fits_nrows,status)
                    if (fits_nrows .ne. nrows) then
                      call ftdrow(main_lum_FITS_file_unit,nrows+1,fits_nrows-nrows,status)
                    end if
! ------------------------------------------------- 
#endif
! GAL_SPECTRA
                    !
                  end if
                end do ! end n_blocks loop 
              end if ! nblock
              !
              close(u2)
              !
            end do ! end ifile loop
            !
            ! close tmp filelist file
            close(u)   
            !
            ! dellocate buffer
            if (allocated(buff_data)) deallocate(buff_data) 
            if (allocated(buff_mags)) deallocate(buff_mags) 
#ifdef GAL_SPECTRA
! -------------------------------------------------        
            if (allocated(buff_spec1)) deallocate(buff_spec1) 
            if (allocated(buff_spec2)) deallocate(buff_spec2) 
            if (allocated(buff_spec3)) deallocate(buff_spec3) 
! ------------------------------------------------- 
#endif
! GAL_SPECTRA
            !
            ! Close FITS file      
            call ftclos(main_lum_FITS_file_unit, status)
            ! check the status     
            if (status .gt. 0) then  
                call IO_print_error_message('Cannot close (luminous properties) FITS output file',only_rank=rank, &
                    called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                write(message,'(a,i4,a,a)') 'ftclos status: ', status, ' for file: ', trim(new_file)
                call IO_print_message(message,only_rank=rank)
                stop   ! stop the program
            end if
            call ftfiou(main_lum_FITS_file_unit, status)
            ! check the status     
            if (status .gt. 0) then  
                call IO_print_error_message('Cannot deallocate unit index for (luminous properties) FITS output file',only_rank=rank, &
                    called_by='IO_create_open_and_print_in_main_FITS_output_files')  
                write(message,'(a,i4,a,a)') 'ftfiou status: ', status, ' for file: ', trim(new_file)
                call IO_print_message(message,only_rank=rank)
                stop   ! stop the program
            end if
            !
        end do
    end if
    !
! ------------------------------------------------- 
#endif
! LUMINOUS_PROCESSES  
    !
! ------------------------------------------------- 
#endif
! TREE_EVOLVE

    return
  end subroutine IO_create_open_and_print_in_main_FITS_output_files
  
  !*****************************************************************************************************************
  !
  ! FUNCTIONS
  !     
  !*****************************************************************************************************************  

end module IO_next
