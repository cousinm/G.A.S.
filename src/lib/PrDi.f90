module PrDi

  use global_variables  ! acces to global defintion and properties
  use mpi               ! The message passing interface module

  public

  !*****************************************************************************************************************
  !
  ! OVERVIEW
  !
  ! Here are defined global parameter used for process distribution.
  ! G.A.S. uses a set of physical (and luminous) processes that allow to compute in parallel
  ! physical evolution and luminous propreties of all modeled galaxies.
  ! The G.A.S model is fully compatible with MPI cluster
  ! the standard scheme uses 2N+1 process
  ! All specific TAG associated to each module (tree/halo/galaxy/disc ...) are also defined in this module
  ! Those TAGs allow to control and follow exchanges between the different galaxy components trhough the different MPI processes
  !
  ! SUBROUTINES IN THIS MODULE
  !
  !  PrDi_attribute_process_function                   : Attribute process labels (main, physical or luminous)
  !      called by :    main program                          Generate the namelist of computer-node (MAP)
  !
  !  PrDi_finalize                                     : Finalize MPI processes
  !      called by :    main program
  !
  ! FUNCTIONS IN THIS MODULE
  !
  !*****************************************************************************************************************

  ! DEFINITION OF GLOBAL VARIABLE LINKED TO PROCESS DISTRIBUTION *******************
  !
  ! MPI process variable, rank, number of processes
  integer(kind=4)                      :: rank                             ! Rank of the processus
  integer(kind=4)                      :: main_process_rank                ! rank of the main process
  integer(kind=4)                      :: nbproc                           ! Number of processus available for the computation
  integer(kind=4)                      :: nbproc_used                      ! Number of processus used for the computation
                                                                           ! Can be smaller (-1) than nbproc
                                                                           ! case of luminous process computation (2N +1 computation sheme)
  integer(kind=4),parameter            :: map_process_unit  = 109          ! map process file
  integer(kind=4)                      :: ierror                           ! Error index (for information about MPI process fail)
  integer(kind=4),parameter            :: prdi_tag          = 2100         ! When messages are transfered from one process to an other
  integer(kind=4),parameter            :: io_tag            = 2200         ! these messages have to be formally identified
  integer(kind=4),parameter            :: tree_tag          = 2300         ! G.A.S. uses a set of tag associated with the different levels
  integer(kind=4),parameter            :: halo_tag          = 2400         ! tree > halo > galaxy > disc > gas ...
  integer(kind=4),parameter            :: gal_tag           = 2500
  integer(kind=4),parameter            :: disc_tag          = 2600
  integer(kind=4),parameter            :: bulge_tag         = 2700
  integer(kind=4),parameter            :: stars_tag         = 2800
  integer(kind=4),parameter            :: dust_tag          = 2900
  integer(kind=4),parameter            :: agn_tag           = 3000
  integer(kind=4),parameter            :: ztable_tag        = 3100
  integer(kind=4)                      :: statut(MPI_STATUS_SIZE)           ! MPI process statut, 0 if all good or > 0 if some problem in the execution of the MPI process
  !
  ! Distribute processes (physical, luminous and main processes)
  logical                              :: main_process                     ! The main process compute tree_files distribution and is used to write output data
  logical                              :: physical_process                 ! Physical processes compute baryonic mechanism following merger trees
  logical                              :: luminous_process                 ! Luminous processes compute luminous properties such as spectra and magnitudes
  integer(kind=4)                      :: nb_physical_processes            ! The number nodes affected to physical process
  integer(kind=4)                      :: nb_luminous_processes            ! The number nodes affected to luminous process
  integer(kind=4),allocatable          :: list_of_physical_processes(:)    ! The complete list of nodes affected to physical process
  integer(kind=4),allocatable          :: list_of_luminous_processes(:)    ! The complete list of nodes affected to luminous process
  character(MPI_MAX_PROCESSOR_NAME)    :: nodename                         ! name of nodes
  !
  ! MPI files distribution parameters, All tree-files are analysed and,
  ! each process receive a list of tree-files in function of the quantity of haloes contains in each files
  ! At the end of the distribution, each process may rougly compute the same number of haloes
  !
  integer(kind=4)                      :: nb_tree_files_computed           ! total number of tree-files sets computed.
  integer(kind=8)                      :: nb_halos_computed                ! Number of halos computed by the processus 'rank'
  integer(kind=4)                      :: nb_files_computed                ! Number of tree-file sets computed by the processus 'rank'

  character(len=10),allocatable        :: list_of_computed_files(:)        ! Indexes of tree-file computed by the processus 'rank'
  !
  ! computationnal time
  character(len=8)                     :: date                             ! the date
  character(len=10)                    :: time                             ! the local time
  character(len=5)                     :: zone                             ! zone
  real(kind=4)                         :: start,current,finish             ! (for each MPI process) starting, current and ending time of computation
  real(kind=4)                         :: hour,mn,sec


  contains

  !*****************************************************************************************************************
  !
  ! SUBROUTINES
  !
  !*****************************************************************************************************************

  subroutine PrDi_attribute_process_function

    ! ATTRIBUTE PROCESS LABEL

    implicit none

    integer(kind=4)         :: i         ! loop index
    integer(kind=4)         :: ierr,res
#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------
    integer(kind=4)         :: j,k    ! loop index
! -------------------------------------------------
#endif
! LUMINOUS_PROCESSES

    character(MAXPATHSIZE)                         :: format
    character(MAXPATHSIZE)                         :: filename
    character(MPI_MAX_PROCESSOR_NAME),allocatable  :: hostname(:)

    ! init nbproc_used to nb_proc,
    ! this initial value can be modified in case of 2N +1 computation scheme
    ! used for parallel computation of luminous and physical properties
    nbproc_used = nbproc

    ! set to null values
    main_process      = .false.
    physical_process  = .false.
    luminous_process  = .false.

    ! by default the main process is the last process
    main_process_rank = nbproc-1

    if (rank .eq. main_process_rank) write(*,'(a)') '|--> PrDi_attribute_process_function'

#ifdef ONLY_LUM
! -------------------------------------------------
    if (rank .eq. main_process_rank) then
        ! the process (nbproc-1) is the main process
        main_process = .true.
    else
        ! other processes compute luminous properties
        luminous_process = .true.
    end if
    !
    if (main_process) then
        ! create list of luminous processes
        nb_physical_processes = 0
        nb_luminous_processes = nbproc_used-1
        allocate(list_of_luminous_processes(nb_luminous_processes))
        do i = 0, nbproc_used-2
            list_of_luminous_processes(i+1) = i
        end do
        write(format,'(a,i3.3,a)') '(a,', nb_luminous_processes, '(i3.3,2x))'
        write(*,format) '|---> List of luminous processes: ', list_of_luminous_processes
    end if
    !
! ------------------------------------------------
#else
! ------------------------------------------------
#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------
    if (nbproc .ge. 3) then
        ! we are in a real multi-processes run
        ! The standard scheme uses 2N +1 processes
        if (rank .eq. main_process_rank) then
            write(*,'(a)') '|---> Take physical & luminous processes into account'
            write(*,'(a)') '|---> Apply 2N+1 process distribution '
        end if
        ! If nbproc is a even number, we have to remove one process
        if (mod(nbproc,2) .eq. 0) then
            nbproc_used = nbproc -1
            if (rank .eq. 0) then
                write(*,'(a)') '|---> 1 process has been disable '
            end if
        end if
        if (rank .eq. main_process_rank) then
            ! the process (nbproc-1) is the main process
            main_process = .true.
        else
          if (rank .le. nbproc_used-2) then
             if (mod(rank,2) .eq. 0) then
                ! even processes compute physical evolution of halos
                physical_process = .true.
             else
                ! odd processes compute luminous properties
                luminous_process = .true.
             end if
          end if
       end if
       !
       if (main_process) then
            ! create list of physical and luminous processes
            nb_physical_processes = (nbproc_used-1)/2
            nb_luminous_processes = nb_physical_processes
            allocate(list_of_physical_processes(nb_physical_processes))
            allocate(list_of_luminous_processes(nb_luminous_processes))
            j = 1
            k = 1
            do i = 0, nbproc_used-2
                if (mod(i,2) .eq. 0) then
                    ! even processes compute physical evolution of halos
                    list_of_physical_processes(j) = i
                    j = j +1
                else
                    ! odd processes compute luminous properties
                    list_of_luminous_processes(k) = i
                    k = k +1
                end if
            end do
            write(format,'(a,i3.3,a)') '(a,', nb_physical_processes, '(i3.3,2x))'
            write(*,format) '|---> List of physical processes: ', list_of_physical_processes
            write(*,format) '|---> List of luminous processes: ', list_of_luminous_processes
        end if
    else
        if (rank .eq. main_process_rank) then
            ! no sufficient processes to distribute physical and luminous processes
            write(*,'(a)') '!!! ERROR in PrDi_attribute_process_function'
            write(*,'(a)') '|---> No enought cpus for parallel computation of physical & luminous properties'
            write(*,'(a)') '|--->   Please used, at least, 2N +1 = 3 processes'
        end if
        call PrDi_finalize
        stop
    end if
!
! ------------------------------------------------
#else
! ------------------------------------------------
! NO LUMINOUS PROCESS
    !
    if (rank .eq. main_process_rank) write(*,'(a)') '|---> Take only physical processes into account'
    if (nbproc .gt. 1) then
       ! we are in a real multi-processes run
       if (rank .eq. main_process_rank) then
          ! the process 'nbproc-1' is the main process
          main_process = .true.
       else
          physical_process = .true.
       end if
    else
       if (rank .eq. main_process_rank) write(*,'(a,i1,a)') '|---> Only 1 process used (rank = ', rank, ')'
       ! only one process is used
       main_process      = .true.
       physical_process  = .true.
    end if
    !
    if (main_process) then
      nb_physical_processes = max(1,nbproc-1)
      nb_luminous_processes = 0
      allocate(list_of_physical_processes(nb_physical_processes))
      !
      do i = 0, nb_physical_processes -1
        list_of_physical_processes(i+1) = i
      end do
      write(format,'(a,i3.3,a)') '(a,', nb_physical_processes, '(i3.3,2x))'
      write(*,format) '|---> List of physical processes: ', list_of_physical_processes
    end if
! -------------------------------------------------
#endif
! LUMINOUS_PROCESSES
! -------------------------------------------------
#endif
! ONLY_LUM
    !
    ! build process MAP
    allocate(hostname(nbproc))
    !
    ! get hostname
    call MPI_Get_processor_name(nodename,res,ierr)
    !
    if (nbproc .gt. 1) then
        if (.not. main_process) then
            ! send my name to the main process
            call MPI_SEND(nodename,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,main_process_rank,prdi_tag+rank,MPI_COMM_WORLD,ierror)
        else
            ! the main process receives names of nodes
            hostname(nbproc) = nodename
            do i = 0, nbproc-2
                call MPI_RECV(hostname(i+1),MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,i,prdi_tag+i,MPI_COMM_WORLD,statut,ierror)
            end do
            !
            ! write process MAP file
            ! build the filename
            write(filename, '(a)') 'MAP_process.log'
            !
            ! open the file
            open(unit = map_process_unit, file = trim(filename), status = 'new')
            write(map_process_unit,'(a)') '# MAP Process log file'
            write(map_process_unit,'(a)') '# rank | hostname'
            do i = 1, nbproc
                if (i-1 .eq. main_process_rank) then
                    write(map_process_unit,'(i3.3,a,a,a,a)') i-1, ' | ', trim(hostname(i)), '[main process]'
                else
                    write(map_process_unit,'(i3.3,a,a,a)') i-1, ' | ', trim(hostname(i))
                end if
            end do
            !
            close(map_process_unit)
            deallocate(hostname)
        end if
    end if

    return
  end subroutine PrDi_attribute_process_function

  !*****************************************************************************************************************

  subroutine PrDi_finalize

    implicit none

    ! Finalize MPI processes
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call MPI_FINALIZE(ierror)

    return
  end subroutine PrDi_finalize

  !*****************************************************************************************************************

  end module PrDi
