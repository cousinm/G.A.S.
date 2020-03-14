module filters
  
  use stellar_population_library   ! Contains stellar population library (e.g. acces to loss mass rate)

  public

  !*****************************************************************************************************************
  ! 
  !  OVERVIEW
  !
  !  filters defines the filter type and the filter_table_list. These filters are used to compute AB magnitudes
  !
  !  SUBROUTINES IN THIS MODULE
  !
  !   filters_void                       : init all properties of a filter
  !
  !   filters_copy                       : copy a filter structure into an other one
  !
  !   filters_deallocate                 : deallocate a filter (set parameters to null values and deallocate dynamical arrays)
  !
  !   filters_load_filters               : Read the filter input list and load corresponding filters
  !  
  !   filters_precompute_interp          : Define a_i, b_i so that for x in [lambda(i);lambda(i+1)]: trans(x) = a_i * x + b_i 
  !
  !   filters_precompute_norm_and_cc     : pre-compute normalization and correction factors 
  !      
  !  FUNCTIONS IN THIS MODULE
  !
  !   No function in this module
  !
  !  PRINTING PROCEDURES
  !
  !   No printing procedure
  !
  !*****************************************************************************************************************

  type filter_type
     character(25)             :: name            ! name of filter 
     integer(kind=4)           :: nlambda         ! number of wavelengths in the filter
     real(kind=4)              :: clambda         ! central wavelenght of the filter
     real(kind=8)              :: norm            ! normalisation factor
     real(kind=8)              :: cc              ! color correction (lamIlam = cte)
     real(kind=4),allocatable  :: lambda(:)       ! wavelength table
     real(kind=4),allocatable  :: a(:),b(:)       ! linear interpolation coefficients between lambda points. (dimension f%nlambda)
  end type filter_type

  ! DEFINE THE FILTER TABLE *************

  integer(kind=4)                :: nfilters       ! Total number of filters through which we compute Mags
  type(filter_type),allocatable  :: filters_tab(:) ! filter table (dimension: nfilters)
  
contains 
  
  !*****************************************************************************************************************
  !
  ! SUBROUTINES
  !
  !*****************************************************************************************************************

  subroutine filters_void(f)
    
    ! INIT A FILTER
    
    implicit none 

    type(filter_type),intent(inout) :: f
      
    f%name      = ''    ! name of filter 
    f%nlambda   = 0     ! number of wavelengths in the filter
    f%clambda   = 0.d0  ! center wavelenght of the filter
    f%norm      = 0.d0  ! normalisation factor
    f%cc        = 0.d0  ! color correction

    return
  end subroutine filters_void

  !*****************************************************************************************************************
  
  subroutine filters_copy(f1,f2)
    
    ! COPY FILTER f2 IN f1 
    
    implicit none 

    type(filter_type),intent(out) :: f1
    type(filter_type),intent(in)  :: f2
      
    f1%name      = f2%name      ! name of filter 
    f1%nlambda   = f2%nlambda   ! number of wavelengths in the filter
    f1%clambda   = f2%clambda   ! center wavelenght of the filter
    f1%norm      = f2%norm      ! normalisation factor
    f1%cc        = f2%cc        ! color correction
    !
    ! arrays
    if (allocated(f2%lambda)) then
        if (.not. allocated(f1%lambda)) then
            allocate(f1%lambda(f1%nlambda))
            f1%lambda = f2%lambda
        else
            call IO_print_error_message('Try to erase pre-existing f1%lambda', &
               only_rank = rank, called_by = 'filters_copy')
        end if
    end if
    !
    if (allocated(f2%a)) then
        if (.not. allocated(f1%a)) then
            allocate(f1%a(f1%nlambda))
            f1%a = f2%a
        else
            call IO_print_error_message('Try to erase pre-existing f1%a', &
               only_rank = rank, called_by = 'filters_copy')
        end if
    end if
    !
    if (allocated(f2%b)) then
        if (.not. allocated(f1%b)) then
            allocate(f1%b(f1%nlambda))
            f1%b = f2%b
        else
            call IO_print_error_message('Try to erase pre-existing f1%b', &
               only_rank = rank, called_by = 'filters_copy')
        end if
    end if

    return
  end subroutine filters_copy

  !*****************************************************************************************************************

  subroutine filters_deallocate(f)
    
    ! INIT A FILTER
    
    implicit none 

    type(filter_type),intent(inout) :: f
      
    call filters_void(f)                          ! void simple fields
    if (allocated(f%lambda)) deallocate(f%lambda) ! wavelength table
    if (allocated(f%a))      deallocate(f%a)      ! first interpolation coefficient
    if (allocated(f%b))      deallocate(f%b)      ! second interpolation coefficient

    return
  end subroutine filters_deallocate

  !*****************************************************************************************************************

  subroutine filters_load_filters
    
    ! READ FILTER LIST AND ASSOCIATED FILTERS
    ! Performs some pre-computations on them to ease magnitude computations later on.

    implicit none
    
    integer(kind=4)                  :: ifilt,ilam,nl,ind

    character(MAXPATHSIZE)           :: filename           ! filename of the input filter list 
    character(MAXPATHSIZE)           :: line               ! a file line
    character(MAXPATHSIZE)           :: filter_path        ! local path for filters 
    character(MAXPATHSIZE)           :: message            ! a message
    
    real(kind=4),allocatable         :: ltrans(:)          ! a local array for filter transmission
    
    type(filter_type),allocatable    :: tmp_filters_tab(:) ! tmp filter table (dimension: nfilters)
                                                           ! filters are read in tmp_filter_tab 
                                                           ! and the organized (in increasing clambda) in filters_tab

#ifdef PRINTALL
! -------------------------------------------------
    if (nbproc .eq. 1) call IO_print_message('filters_load_filters')
! -------------------------------------------------
#endif

    if (main_process .or. luminous_process) then
      ! build the filename with the input path
      write(filename,'(a,a)') trim(input_path), '/filters/filters_list.in' 
      ! open the file, using filter_unit
      open(unit = filter_list_unit,file = trim(filename),status='old')
      do
        read(filter_list_unit, '(a)') line
        if (trim(line) .eq. 'START') then
          !
          read(filter_list_unit,*) nfilters  ! read the number of filters used 
          !
          ! allocate filter_table_list(:)
          allocate(tmp_filters_tab(nfilters))         
          !
          do ifilt = 1, nfilters
            ! read the local filter path (in the input/filters input path)
            read(filter_list_unit,'(a)') filter_path
            write(filename,'(a,a,a)') trim(input_path), '/filters/', trim(filter_path) 
            open(unit = filter_unit,file = trim(filename),status='old')
            !
            ! init the filter 
            call filters_void(tmp_filters_tab(ifilt))
            !
            do
              read(filter_unit, '(a)') line
              if (trim(line) .eq. 'START') then
                !
                read(filter_unit,*) tmp_filters_tab(ifilt)%name     ! read the name of the filter
                read(filter_unit,*) nl                              ! read the number of wavelenght
                !
                ! The filter band pass must have null transmission at both end of the wavelenght range
                ! Add two null cells
                tmp_filters_tab(ifilt)%nlambda = nl +2
                nl = tmp_filters_tab(ifilt)%nlambda                 ! create a local copy
                !
                read(filter_unit,*) tmp_filters_tab(ifilt)%clambda  ! read the value of the central wavelenght
                if (luminous_process) then
                  !
                  ! allocate arrays
                  allocate(tmp_filters_tab(ifilt)%lambda(nl))         ! wavelenght table
                  tmp_filters_tab(ifilt)%lambda(:) = 0.d0             ! init
                  allocate(ltrans(nl))                                ! local transmission table
                  ltrans(:) = 0.d0                                    ! init
                  !
                  ! read lambda table and transmission
                  do ilam = 2, nl-1 
                    read(filter_unit,*) tmp_filters_tab(ifilt)%lambda(ilam),ltrans(ilam)  ! lambda is in microns, trans is without unit 
                  end do
                  !
                  ! complete wavelenght table
                  tmp_filters_tab(ifilt)%lambda(1) = tmp_filters_tab(ifilt)%lambda(2)*0.999
                  tmp_filters_tab(ifilt)%lambda(nl) = tmp_filters_tab(ifilt)%lambda(nl -1)*1.001
                  !
                  ! define interpolation coefficients (a and b) to interpolate filter curve onto any lambda grid later on
                  call filters_precompute_interp(tmp_filters_tab(ifilt),ltrans)
                  ! 
                  ! compute normalisation factor and color correction factor
                  call filters_precompute_norm_and_cc(tmp_filters_tab(ifilt),ltrans)
                  deallocate(ltrans)
                  !
                end if ! luminous process 
                exit  ! quit do loop
              end if  ! START  
              if (line(1:1) .eq. '#') then
                cycle ! header or something like this (skip)
              else
                write(message,'(a,a)') 'Impossible to read a line in the filter file: ', trim(filter_path)
                call IO_print_error_message(message,only_rank=rank,called_by='cooling_function_read_Lambda_table') 
                stop ! stop the program
              end if
            end do ! reading loop (in the filter file)   
            ! close the current filter file
            close(filter_unit)
          end do  ! filters loop
          exit  ! quit do loop
        end if ! START 
        if (line(1:1) .eq. '#') then
            cycle ! header or something like this (skip)
        else
          call IO_print_error_message('Impossible to read a line in a filter file',only_rank=rank,called_by='filters_load_filters') 
          stop ! stop the program
        end if
      end do ! reading loop (in the filter list file)
      ! close filter list file
      close(filter_list_unit)  
      !
      ! organize filters_tab (filters are listed in increasing order of clambda) 
      allocate(filters_tab(nfilters))
      do ifilt = 1, nfilters
        ind = maxloc(tmp_filters_tab%clambda,dim=1)                            ! locate the next filter
        call filters_copy(filters_tab(nfilters-ifilt+1),tmp_filters_tab(ind))  ! copy
        call filters_deallocate(tmp_filters_tab(ind))                          ! erase
      end do
      !
      ! deallocate tmp_filters_tab
      if (allocated(tmp_filters_tab)) deallocate(tmp_filters_tab)
      ! 
      ! Print some informations
      if (main_process) then
        write(message,'(a,i2.2,a)') 'Load ', nfilters, ' filters'
        call IO_print_message(message)
        if (nfilters .lt. 5) then
            write(message,'(a,a)') '( ', trim(filters_tab(1)%name)
            do ifilt = 2, nfilters
                write(message,'(3(a))') trim(message), ', ', trim(filters_tab(ifilt)%name)
            end do
            write(message,'(a,a)') trim(message), ' )'
            call IO_print_message(message)
        else
            write(message,'(9(a))') '( ', trim(filters_tab(1)%name), ', ', trim(filters_tab(2)%name), &
                ' ... ', trim(filters_tab(nfilters-1)%name), ', ', trim(filters_tab(nfilters)%name), ' )'
            call IO_print_message(message)
        end if
      end if 
      !
    end if

    return
  end subroutine filters_load_filters

  !*****************************************************************************************************************

  subroutine filters_precompute_interp(f,ltrans)
    
    ! DEFINE a_i, b_i SO THAT FOR X IN [lambda(i);lambda(i+1)]: trans(x) = a_i * x + b_i 

    implicit none
    
    integer(kind=4)                 :: ilam      ! loop index

    real(kind=4),intent(in)         :: ltrans(:) ! the transmission function

    type(filter_type),intent(inout) :: f         ! a given filter
    
    if (f%nlambda .eq. 0) then
      call IO_print_error_message('f%nlambda = 0',only_rank=rank,called_by='filters_precompute_interp') 
      stop ! stop the program
    end if  

    if (f%nlambda .ne. size(ltrans)) then
      call IO_print_error_message('f%nlambda /= size(ltrans)',only_rank=rank,called_by='filters_precompute_interp') 
      stop ! stop the program
    end if     
    
    if ((allocated(f%a)) .or. (allocated(f%b))) then
      call IO_print_error_message('f%a or f%b already allocated',only_rank=rank,called_by='filters_precompute_interp') 
      stop ! stop the program
    end if

    ! allocate interpolation coefficient arrays
    allocate(f%a(f%nlambda))
    allocate(f%b(f%nlambda))

    do ilam = 1,f%nlambda -1
      ! trans(x) = a_i * x + b_i 
      f%a(ilam) = (ltrans(ilam+1) - ltrans(ilam)) / (f%lambda(ilam+1) - f%lambda(ilam))
      f%b(ilam) = ltrans(ilam) - f%a(ilam) * f%lambda(ilam)
    end do
    ! filter's last points is at zero
    ilam = f%nlambda 
    f%a(ilam) = 0.0d0
    f%b(ilam) = 0.0d0

    return
  end subroutine filters_precompute_interp

  !*****************************************************************************************************************

  subroutine filters_precompute_norm_and_cc(f,trans)

    ! DEFINE THE COLOR CORRECTION LINKED TO A GIVEN FILTER
    ! Indeed, we assume that all monochromatic fluxes are given in the case of a lam Ilam = cte spectrum

    implicit none

    real(kind=8)                    :: cte      ! constant : lam Ilam = cte (the reference spectrum)
    real(kind=4),intent(in)         :: trans(:) ! the transmission function

    type(filter_type),intent(inout) :: f        ! a given filter

    ! The transmission function must be normalized
    
    f%norm = trap(trans,f%lambda)   
    
    ! compute cte --> 1/cte = int(T(lam)/lam dlam)
!~     cte = trap(trans/f%lambda,f%lambda)/f%norm
    cte = 1.d0
    
    ! compute the color correction coefficient
    ! integration with filter is done with Flamb
    ! magnitude are computed by using Fnu
    ! f%cc allow to compute Fv from Flamb (nu*Fnu = lamb*Flamb --> Fnu = lamb**2*Ilamb/c)
    ! lambda is in microns, light_speed have to be in micron/s
    ! spectrum is given in lamb*Flamb --> Lsun, convert to erg/s/Hz
    f%cc = L_Sun_erg_s*f%clambda**2./(light_speed_mPerSec*1.e6*cte)

    return
  end subroutine filters_precompute_norm_and_cc
  
  !*****************************************************************************************************************
  
  subroutine filters_spec2mag(spectrum,z,Mags,frame)
  
    ! COMPUTE MAGNITUDES BY APPLYING A LIST OF OBSERVATORY FILTERS ONTO A GALAXY SPECTRUM
  
    implicit none
    
    integer(kind=4)                        :: ifilt                ! loop index onto filter list
    integer(kind=4)                        :: l, l_inf, l_sup      ! indexes
    integer(kind=4)                        :: l2, l3               ! loop indexes onto wavelenght
    integer(kind=4)                        :: nwaves_in_filter
    
    character(*),intent(in)                :: frame                ! allows to select galaxy-frame or observer-frame
    character(MAXPATHSIZE)                 :: message              ! a message to display
    
    real(kind=8),intent(in)                :: z                    ! redshift of the galaxy        
    real(kind=4),intent(in)                :: spectrum(:)          ! spectrum of the galaxy
    real(kind=4),intent(inout)             :: Mags(:)              ! Magnitudes
    real(kind=4),allocatable               :: tmp_Waves(:)
    real(kind=4),allocatable               :: tmp_spectrum(:)
    real(kind=4)                           :: wave                 ! current wavelenght
    real(kind=4)                           :: Wave_min, Wave_max   ! miminum and maximum wavelenght of the spectra table
    real(kind=4)                           :: fWave_min, fWave_max ! miminum and maximum wavelenght of the filter table
    real(kind=4),allocatable               :: filter_trans(:)      ! tranmission of the current filter interpolated onto the galaxy spectra wavelenght table
    real(kind=8)                           :: Fnu                  ! Flux in erg/s/Hz
    real(kind=8)                           :: dil_factor           ! In case of Observer frame computation, The spectrum are diluted due to 1+z factor
    
    ! create a local copy of Waves
    allocate(tmp_Waves(nWaves))
    !
    ! in input the spectrum is given in lamb*Flamb
    ! flux in filters are computed by using int(Flam*T(lamb)*dlamb)
    ! we have to convert spectum in Flamb
    ! create a local copy of spectrum
    allocate(tmp_spectrum(nWaves))
    tmp_spectrum = spectrum/Waves  ! Flamb  Lsun/micron
    !
    select case (trim(frame)) 
    case ('Observer','obs','observer')
        ! Observer-frame
        ! Apply redshift correction onto wavelenght scale
        tmp_Waves  = real((1.d0 + z),4)*Waves
        dil_factor = 2.5d0*log10((1.d0 + z)**2.)
    case ('Galaxy','galaxy','gal') 
        ! Galaxy-frame
        ! Copy Waves
        tmp_Waves  = Waves
        dil_factor = 0.d0
    case default
        write(message,'(a,a,a)') 'Keyword ', trim(frame), ' not defined'
        call IO_print_error_message(message,only_rank=rank,called_by='filters_spec2mag')    
        stop  ! stop the program   
    end select
    ! 
    Wave_min = minval(tmp_Waves)
    Wave_max = maxval(tmp_Waves)
    !   
    l = 1
    do ifilt = 1, nfilters
        !
        ! test if the filter is applicable
        fWave_min = minval(filters_tab(ifilt)%lambda) 
        fWave_max = maxval(filters_tab(ifilt)%lambda) 
        if ((fWave_min .ge. Wave_min) .and. (fWave_max .le. Wave_max)) then
            ! the filter is inluded in the whole wavelenght range
            do while(tmp_Waves(l) .gt. fWave_min)
                l = l -1
            end do
            ! @ this point tmp_Waves(l) < fWave_min
            l_inf = nWaves
            do while(tmp_Waves(l) .le. fWave_max)
                if ((tmp_Waves(l) .ge. fWave_min) .and. (tmp_Waves(l) .lt. tmp_Waves(l_inf))) l_inf = l
                l = l +1
            end do
            l_sup = l -1
            ! @ this point tmp_Waves(l) > fWave_max
            !
            ! define the filter transmission
            nwaves_in_filter = l_sup - l_inf +1 +2
            allocate(filter_trans(nwaves_in_filter))
            filter_trans = 0.
            !
            l3 = 1
            do l2 = 2, nwaves_in_filter-1
                wave = tmp_Waves(l_inf+l2-2)
                do while (filters_tab(ifilt)%lambda(l3) .le. wave)
                    l3 = l3 +1
                end do
                l3 = max(l3 -1,1)
                filter_trans(l2) = filters_tab(ifilt)%a(l3)*wave + filters_tab(ifilt)%b(l3)
            end do
            ! @ this point the filter transmission is interpolated onto the spectrum wavelenght table 
            !
            ! compute flux
            ! spectrum is in Flamb
            Fnu = filters_tab(ifilt)%cc*real(trap(tmp_spectrum(l_inf:l_sup)*filter_trans,tmp_Waves(l_inf:l_sup)),8)/filters_tab(ifilt)%norm             
            ! @ this point Fnu is in erg/s/Hz
            !
            ! compute magnitude
            ! 100.19478 is 2.5*Log10(4 pi [10pc in cm]^2)
            Mags(ifilt) = real(-2.5d0*log10(Fnu) - 48.6d0 + 100.19478d0 + dil_factor,4)
            !
            ! deallocate filter_trans
            if (allocated(filter_trans)) deallocate(filter_trans)
        end if
    end do ! filters
    
    return
  end subroutine filters_spec2mag
 
  !*****************************************************************************************************************

end module filters
