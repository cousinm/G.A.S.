module dust

  use stellar_population_library  ! Contains stellar population library (e.g. acces to loss mass rate, spectra)

  public
  
  !*****************************************************************************************************************
  ! 
  !  OVERVIEW
  ! 
  !  The dust module defines all properties and all procedures assossiated to a dust object. 
  !  A dust object is linked to a gas object
  !  A dust component is associated to each gas component used in the model: 
  !    - Diffuse gas in disc and bulge
  !    - dense/structured/fragmented gas in the disc
  !    - gas torus associated to the AGN component
  !    - ...
  !  It groups the dust mass, extinction geometyrical functions, attenuation, etc ...
  !  In the header of the module are defined output properties of a dust component (labels, units and formats)
  !
  !  MPI PROCESSES IN THIS MODULE
  !
  !   dust_send_data                      : send specific informations about a dust component
  !
  !   dust_receive_data                   : receive specific informations about a dust component
  !
  !  SUBROUTINES IN THIS MODULE
  !
  !   dust_void                          : initialize a dust component
  !
  !   dust_copy                          : copy a dust component in a other one
  !      called by: disc_copy
  !
  !   dust_mass_fraction                 : compute the effective fraction (in amss) of each component constituing the dust
  !
  !   dust_evolve                        : evolve dust properties (mass, extinction, temperature)
  !      called by : disc_evolve
  !
  !   dust_read_abs_sca_properties       : read absorbtion and scattering properties of dust 
  !      called by : main program 
  !
  !   dust_read_dust_SEDs                : read dust ir spectral energy distribution
  !      called by : main program 
  !
  !   dust_compute_effective_extinction  : compute effective extinction due to the dust component
  !
  !   dust_build_dust_spectrum           : build the spectral energy distribution of a given dust component
  !
  !  FUNCTIONS IN THIS MODULE
  !
  !   dust_mass                          : return the mass of the dust component
  !
  !   dust_extinction                    : return the FUV extinction (based on Boquien+13)  
  !
  !   dust_AV                            : return A(V) from the effective extinction profile of a dust component
  !
  !   dust_E_BV                          : return E(B-V) from the effective extinction profile of a dust component
  !
  !  PRINTING PROCEDURES
  !
  !   dust_load_dust_data                : load dust data, create the dust output properties list 
  !      called by : disc_load_disc_data
  !
  !*****************************************************************************************************************

  ! DEFINITION OF GLOBAL VARIABLES LINKED TO THE DUST STRUCTURE *******************

  ! DUST_TYPE DEFINITION **********************************************************

  type dust_type
     character(6)             :: geom         ! geometry distribution 'slab(disc) or dwek(bulge)'
     real(kind=8)             :: mass         ! mass of dust
     real(kind=8)             :: MZ_MH        ! metal mass over hydrogen mass (use to convert M_pah/M_H in M_pah/M_dust)
     real(kind=8)             :: incl         ! inclination (a local copy of disc%incl) 
     real(kind=8)             :: f_pah        ! fraction of pah
     real(kind=8)             :: f_bg         ! fraction of big grain
     real(kind=8)             :: tau          ! ISM extinction in FUV or V band 
     real(kind=4),allocatable :: eff_ext(:)   ! effective extinction
     real(kind=4),allocatable :: spectrum(:)  ! dust spectrum erg/s/NH
  end type dust_type

  ! PARAMETERS **************************************
  
  integer(kind=4)           :: ndtypes            ! number of grain types
  integer(kind=4)           :: nISRFBins          ! number of ISRF bins in dust SEDs
  
;  character(4),allocatable  :: dust_types(:)      ! names of dust types
  
  real(kind=8),parameter    :: f_pah_ref    = 4.58d-2 
  real(kind=8),parameter    :: f_vsg_ref    = 16.0d-2 
  real(kind=8),parameter    :: O_over_H_min = 6.5d0
  real(kind=8),parameter    :: MZ_MH_solar  = 2.d-2
  
  ! SPECTRAL ENERGY DISTRIBUTION OF DUST ************

  real(kind=4),allocatable  :: dust_sed(:,:,:,:) ! Infrared spectrum for different grain types and ISRF intensities 
                                                 ! The library contains SED for young and old stellar populations sp_age = 1, 2 
                                                 ! (nWaves, nISRFbins, sp_age, ndtypes)

  ! ABSORBTION & SCATTERING PROPERTIES **************

  real(kind=8),allocatable  :: dust_Abs(:,:)     ! dust absorbtion (gt,lambda)
  real(kind=8),allocatable  :: dust_Sca(:,:)     ! dust scattering (gt,lambda)
  
  ! DEFINE HEADER INFORMATIONS **********************
  
  integer(kind=4),parameter :: nb_dust_field = 4  ! Number of dust properties saved

contains

  !*****************************************************************************************************************
  ! 
  ! MPI PROCESSES
  !
  !*****************************************************************************************************************

#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------  
  subroutine dust_send_data(dust)                                     
  
    ! SEND SPECIFIC INFORMATIONS ABOUT A DUST COMPONENT
  
    implicit none

    logical                      :: go_down

    type(dust_type),intent(in)   :: dust         ! dust component
    
    ! data are sent by physical process (odd) and receive by luminous process (even)
    
    go_down = .false.
    
    if (dust%tau .gt. 0.d0) then
        ! there is some dust
        !
        go_down = .true.
        call MPI_SEND(go_down,1,MPI_LOGICAL,rank+1,dust_tag+1,MPI_COMM_WORLD,ierror)
        !
        ! send dust geom
        call MPI_SEND(dust%geom,8,MPI_CHAR,rank+1,dust_tag+2,MPI_COMM_WORLD,ierror)
        ! send dust mass
        call MPI_SEND(dust%mass,1,MPI_REAL8,rank+1,dust_tag+3,MPI_COMM_WORLD,ierror)
        ! send Mz_MH
        call MPI_SEND(dust%MZ_MH,1,MPI_REAL8,rank+1,dust_tag+5,MPI_COMM_WORLD,ierror)
        ! send fraction of pah
        call MPI_SEND(dust%f_pah,1,MPI_REAL8,rank+1,dust_tag+6,MPI_COMM_WORLD,ierror)
        ! send fraction of big grain
        call MPI_SEND(dust%f_bg,1,MPI_REAL8,rank+1,dust_tag+7,MPI_COMM_WORLD,ierror)
        ! send dust inclination
        call MPI_SEND(dust%incl,1,MPI_REAL8,rank+1,dust_tag+8,MPI_COMM_WORLD,ierror)
        ! send dust extinction
        call MPI_SEND(dust%tau,1,MPI_REAL8,rank+1,dust_tag+9,MPI_COMM_WORLD,ierror)
    else
        ! send exit loop order
        call MPI_SEND(go_down,1,MPI_LOGICAL,rank+1,dust_tag+1,MPI_COMM_WORLD,ierror)
    end if
    
    return
  end subroutine dust_send_data 
  
  !*****************************************************************************************************************
  
  subroutine dust_receive_data(dust)                                     
  
    ! RECEIVE SPECIFIC INFORMATIONS ABOUT A DUST COMPONENT
  
    implicit none

    logical                      :: go_down
    
    type(dust_type),intent(out)  :: dust
    
    ! data are sent by physical process (odd) and receive by luminous process (even)
    
    call dust_void(dust) ! init
    
    ! receive exit loop order
    call MPI_RECV(go_down,1,MPI_LOGICAL,rank-1,dust_tag+1,MPI_COMM_WORLD,statut,ierror)
    
    if (go_down) then
        ! receive dust information
        ! receive dust geom
        call MPI_RECV(dust%geom,8,MPI_CHAR,rank-1,dust_tag+2,MPI_COMM_WORLD,statut,ierror)
        ! receive dust mass
        call MPI_RECV(dust%mass,1,MPI_REAL8,rank-1,dust_tag+3,MPI_COMM_WORLD,statut,ierror)
        ! receive Mz_MH
        call MPI_RECV(dust%MZ_MH,1,MPI_REAL8,rank-1,dust_tag+5,MPI_COMM_WORLD,statut,ierror)
        ! receive fraction of pah
        call MPI_RECV(dust%f_pah,1,MPI_REAL8,rank-1,dust_tag+6,MPI_COMM_WORLD,statut,ierror)
        ! receive fraction of big grain
        call MPI_RECV(dust%f_bg,1,MPI_REAL8,rank-1,dust_tag+7,MPI_COMM_WORLD,statut,ierror)
        ! receive dust inclination
        call MPI_RECV(dust%incl,1,MPI_REAL8,rank-1,dust_tag+8,MPI_COMM_WORLD,statut,ierror)
        ! receive dust extinction in ISM
        call MPI_RECV(dust%tau,1,MPI_REAL8,rank-1,dust_tag+9,MPI_COMM_WORLD,statut,ierror)
    end if
    
    return
  end subroutine dust_receive_data 
! ------------------------------------------------- 
#endif
! LUMINOUS_PROCESSES

  !*****************************************************************************************************************
  !
  ! SUBROUTINES
  !
  !*****************************************************************************************************************

  subroutine dust_void(dust,init_geom)
    
    ! SET TO NULL VALUES ALL DUST PROPERTIES

    implicit none

    character(6),intent(in),optional  :: init_geom ! the initial geometry of dust (slab for disc and dwek for bulge)
    type(dust_type),intent(inout)     :: dust      ! a dust component

    if (present(init_geom)) then
        dust%geom = init_geom
    else
        dust%geom = ''
    end if
    dust%mass   = 0.d0
    dust%MZ_MH  = 0.d0
    dust%incl   = 0.d0
    dust%f_pah  = -1.d0
    dust%f_bg   = -1.d0
    dust%tau    = -1.d0
    if (allocated(dust%eff_ext)) dust%eff_ext(:)   = 0.d0
    if (allocated(dust%spectrum)) dust%spectrum(:) = 0.d0
    
    return
  end subroutine dust_void

  !*****************************************************************************************************************
  
  subroutine dust_copy(d1,d2)
    
    ! COPY A DUST COMPONENT IN AN OTHER 

    implicit none

    type(dust_type),intent(inout)  :: d1 ! a dust component
    type(dust_type),intent(in)     :: d2 ! a dust component

    d1%geom   = d2%geom
    d1%mass   = d2%mass
    d1%MZ_MH  = d2%MZ_MH
    d1%incl   = d2%incl
    d1%f_pah  = d2%f_pah
    d1%f_bg   = d2%f_bg
    d1%tau    = d2%tau
    !
    if (allocated(d2%eff_ext)) then
        if (.not. allocated(d1%eff_ext)) then
            allocate(d1%eff_ext(nWaves))
            d1%eff_ext = d2%eff_ext
        else
            call IO_print_error_message('Try to erase pre-existing dust%eff_ext', &
               only_rank = rank, called_by = 'dust_copy')
        end if
    end if
    !
    if (allocated(d2%spectrum)) then
        if (.not. allocated(d1%spectrum)) then
            allocate(d1%spectrum(nWaves))
            d1%spectrum = d2%spectrum
        else
            call IO_print_error_message('Try to erase pre-existing dust%spectrum', &
               only_rank = rank, called_by = 'dust_copy')
        end if
    end if
    !
    return
  end subroutine dust_copy
  
  !*****************************************************************************************************************
  
  subroutine dust_mass_fraction(dust,gas)
    
    ! COMPUTE THE FRACTION OF PAH, VSG and BG IN THE DUST COMPONENT ! M_pah/M_dust
    ! use Remy-Ruyer+15 Eq 5 (Metallicity)
    ! The relation measured by Remy-Ruyer+15 shows a scatter of 0.35 dex. 
    ! We assume a dispersion for the log normal distribution of 3sig = 0.35 --> sig = 0.12  

    implicit none
   
    real(kind=8)                   :: M_H, M_Z
    real(kind=8)                   :: lf
    real(kind=8)                   :: O_over_H   ! oxygen to hydrogen relative abundance
    real(kind=8)                   :: r          ! a random number
    real(kind=8)                   :: f_min_bg = 4.d0/5.d0
    real(kind=8)                   :: new_f_pah
    real(kind=8)                   :: new_f_bg
    
    type(dust_type),intent(inout)  :: dust       ! the dust component
    
    type(gas_type),intent(in)      :: gas        ! the gas component of the remnant disc

    dust%f_pah = -1.d0; dust%f_bg = -1.d0
    
#ifdef PRINTALL 
    ! -------------------------------------------------
    call IO_print_message('dust_mass_fraction',only_rank=rank,component='dust')
    ! -------------------------------------------------
#endif
 
#ifndef NO_EXTINCTION
! ------------------------------------------------- 
    !   
    ! Oxygen relative abundance (compare to Hydrogen, in number)
    O_over_H = O_H(gas)
    
    if (O_over_H .le. 0.d0) return ! no gas enough 

    if (12.d0 + log10(O_over_H) .lt. O_over_H_min) return 

    M_H = gas_mass(gas,component='H1')
    M_Z = gas_mass(gas,component='metals')
    !
    dust%MZ_MH = M_Z / M_H
    ! We apply a shift of 0.2 dex
    ! To take into account variations due to the different prescriptions (from PP04 O3N2 to Tex measurments)
    lf = -11.d0 + 1.3d0*(12.d0 + log10(O_H(gas)) -2.d-1) 
    lf = max(-6.d0, min(5.d-1,normal_distribution(lf,2.d-1))) 
    new_f_pah = (10.**lf)*f_pah_ref
    ! 
    call random_number(r)
    !
    if (dust%f_pah .lt. 0.d0) then
        ! init
        ! M_pah/M_dust [0% : 10%]
        dust%f_pah = new_f_pah 
    else
        ! The new value is randomly taken between the previous and the new reference value: new_f_pah
        dust%f_pah = max(0.d0,dust%f_pah + r*(new_f_pah-dust%f_pah)) 
    end if  
    ! 
    new_f_bg = f_min_bg*(1.d0-dust%f_pah)
    if (dust%f_bg .lt. 0.d0) then
        ! init
        dust%f_bg = new_f_bg
    else
        ! We assume that the target mass fraction of big grain is almost "f_min_bg" of the residual dust (1.d0-dust%f_pah)
        dust%f_bg = max(0.d0,dust%f_bg + r*(new_f_bg - dust%f_bg)) ! M_bg/M_dust
    end if
    !
    if ((dust%f_bg .le. 0.d0) .or. (dust%f_pah .le. 0.d0)) then
        call IO_print_error_message('grain fraction < 0. !',only_rank = rank, called_by = 'dust_mass_fraction')
        call IO_print_message('with',only_rank=rank,component='gal', &
                   param_name = (/'f_bg                     ','f_pah                    ','O_over_H                 '/), &
                   real_param_val = (/dust%f_bg,dust%f_pah,O_over_H/))
        stop
    end if
    
! -------------------------------------------------
#endif 

#ifdef PRINTALL 
    ! -------------------------------------------------
    call IO_print_message('dust_mass_fraction ... done',only_rank=rank,component='dust')
    ! -------------------------------------------------
#endif

    return
  end subroutine dust_mass_fraction
   
  !*****************************************************************************************************************  
  
  subroutine dust_evolve(d,gas,rc,incl)
  
    ! COMPUTE DUST PROPERTIES
    
    implicit none

    real(kind=8),intent(in)           :: rc    ! the half mass radius of the host component (from rd or rb)
    real(kind=8),intent(in),optional  :: incl  ! disc inclination (create a local copy for dust extinction computation)
    
    type(gas_type),intent(in)         :: gas   ! the gas component of the host component
    
    type(dust_type),intent(inout)     :: d     ! the dust component
    
#ifdef PRINTALL 
! -------------------------------------------------
    call IO_print_message('dust_evolve',only_rank=rank,component='dust')
! -------------------------------------------------
#endif
    
    ! dust inclination
    d%incl = 0.d0
    if (present(incl)) d%incl = incl
    ! compute the fraction of pah and bg
    call dust_mass_fraction(d,gas)
    ! compute the dust FUV extinction
    d%tau = dust_extinction(d,gas,rc)

#ifdef PRINTALL 
! -------------------------------------------------
    call IO_print_message('dust_evolve ... done',only_rank=rank,component='dust')
! -------------------------------------------------
#endif
    
    return
  end subroutine dust_evolve

  !*****************************************************************************************************************
 
#ifdef LUMINOUS_PROCESSES
! -------------------------------------------------    
  subroutine dust_read_abs_sca_properties
  
    ! READ ABSORBTION AND SCATTERING PROPERTIES OF DUST
    
    implicit none
    
    integer(kind=4)          :: gt,l         ! index loop under dust types and wavelenght
    integer(kind=4)          :: nWaves_dust  ! local variable (process 0, checks if nWaves = nWaves_dust)
     
    character(MAXPATHSIZE)   :: filename
    character(MAXPATHSIZE)   :: message      ! a message  
    character(MAXPATHSIZE)   :: line
    
    if (main_process .or. luminous_process) then
    
        call IO_print_message('dust_read_abs_sca_properties')

        ! Build the input filename
        write(filename,'(a,a,a,a)') trim(input_path), '/dust_abs_sca.in'

        write(message,'(a,a,a,a)') 'Load data from : ', 'dust_abs_sca.in'
        call IO_print_message(message)
        ! Open the library file
        open(unit = abs_sca_unit, file = filename, status = 'old')
       
        do
            read(abs_sca_unit, '(a)', end = 2) line
            if (trim(line) .eq. 'START') then
                ! all lines before the START keywork are header lines: they contain informations about data 
                !
                read(abs_sca_unit,*) nWaves_dust, ndtypes
                !
                if (main_process) then
                    ! The main process checks if nWaves = nWaves_dust
                    ! If it is not the case the main process kills the run
                    if (nWaves_dust  .ne. nWaves) then
                        call IO_print_error_message('nWaves_dust /= nWaves',only_rank=rank,called_by='dust_read_abs_sca_properties') 
                        stop ! stop the program 
                    end if
                end if
                !
                ! allocate arrays
                ! dust type
                allocate(dust_types(ndtypes))
                ! for Absorption
                allocate(dust_Abs(ndtypes,nWaves))
                ! for scattering
                allocate(dust_Sca(ndtypes,nWaves))
                !
                read(abs_sca_unit,*) dust_types(1:ndtypes)
                read(abs_sca_unit,*) ! skip wavelenght list
                !
                ! read abs
                do gt = 1, ndtypes
                    read(abs_sca_unit,*) dust_Abs(gt,1:nWaves)
                end do
                !
                do gt = 1, ndtypes
                    read(abs_sca_unit,*) dust_Sca(gt,1:nWaves)
                end do
                !
                if (main_process) then
                    ! Print some informations 
                    write(message,'(i1.1,a)') ndtypes, ' dust types are taken into account'
                    call IO_print_message(message)
                    write(message,'(a)') trim(dust_types(1))
                    do gt = 2, ndtypes
                        write(message,'(a,a,a)') trim(message), ', ', trim(dust_types(gt))
                    end do
                    call IO_print_message(message)
                    call IO_print_message('use', param_name=(/'nWaves                   '/),&
                        int_param_val=(/nWaves_dust/))
                end if
                exit  ! quit do loop 
            end if
            !
            if (line(1:1) .eq. '#') then
                cycle ! header or something like this (skip)
            else
                call IO_print_error_message('Impossible to read the line',only_rank=rank,called_by='dust_read_abs_sca_properties') 
                stop ! stop the program
            end if
        end do
        !
        l = 1
        ! Initialize i_14eV, i_FUV, i_6eV, i_B and i_V
        ! 6 eV    => 0.20664 mu m
        ! 13.6 eV => 0.09116 mu m
        do while (Waves(l) < 1.0)
            if ((Waves(l) < 0.0912) .and. (Waves(l+1) > 0.0912)) then
                i_14eV = l
                write(message,'(a,i3.3,a,f6.4,a)') 'i_14eV: ', i_14eV, ' l_14eV: ', Waves(i_14eV), ' microns'
                call IO_print_message(message)
            end if
            if ((Waves(l) < 0.15) .and. (Waves(l+1) > 0.15)) then
                i_FUV = l
                write(message,'(a,i3.3,a,f6.4,a)') 'i_FUV : ', i_FUV, ' l_FUV : ', Waves(i_FUV), ' microns'
                call IO_print_message(message)
            end if
            if ((Waves(l) < 0.207) .and. (Waves(l+1) > 0.207)) then
                i_6eV = l
                write(message,'(a,i3.3,a,f6.4,a)') 'i_6eV : ', i_6eV, ' l_6eV : ', Waves(i_6eV), ' microns'
                call IO_print_message(message)
            end if
            if ((Waves(l) < 0.43) .and. (Waves(l+1) > 0.43)) then
                i_B = l
                write(message,'(a,i3.3,a,f6.4,a)') 'i_B   : ', i_B, ' l_B   : ', Waves(i_B), ' microns'
                call IO_print_message(message)
            end if
            if ((Waves(l) < 0.55) .and. (Waves(l+1) > 0.55)) then
                i_V = l
                write(message,'(a,i3.3,a,f6.4,a)') 'i_V   : ', i_V, ' l_V   : ', Waves(i_V), ' microns'
                call IO_print_message(message)
                exit
            end if
            l = l +1
        end do
        !
    end if
    
2   close(abs_sca_unit)
    
    return  
  end subroutine dust_read_abs_sca_properties
  
  !*****************************************************************************************************************
  
  subroutine dust_read_dust_SEDs
  
    ! READ DUST SEDs
    ! In this new version, dust SEDs are parametrized by the InterStellar Radiation Field (ISRF)
    ! Dust SEDs are available for ndtypes differents types of dusts
    
    implicit none
    
    integer(kind=4)          :: i,l          ! loop indexes under ISRF and wavelenght
    integer(kind=4)          :: nWaves_dust  ! local variable (process 0, checks if nWaves  = nWaves_dust)
    integer(kind=4)          :: ndtypes_seds ! local variable (process 0, checks if ndtypes = ndtypes_seds)
     
    character(MAXPATHSIZE)   :: filename
    character(MAXPATHSIZE)   :: message      ! a message  
    character(MAXPATHSIZE)   :: line
    
    if (main_process .or. luminous_process) then
    
        call IO_print_message('dust_read_dust_SEDs')

        ! Build the input filename
        write(filename,'(a,a)') trim(input_path), '/dust_sed.in'

        write(message,'(a,a)') 'Load data from : ', 'dust_sed.in'
        call IO_print_message(message)
        ! Open the library file
        open(unit = dust_sed_unit, file = filename, status = 'old')
       
        do
            read(dust_sed_unit, '(a)', end = 2) line
            if (trim(line) .eq. 'START') then
                ! all lines before the START keywork are header lines: they contain informations about data 
                !
                read(dust_sed_unit,*) nWaves_dust, nISRFBins, ndtypes_seds
                !
                if (main_process) then
                    ! The main process checks if nWaves = nWaves_dust
                    ! If it is not the case the main process kills the run
                    if (nWaves_dust  .ne. nWaves) then
                        call IO_print_error_message('nWaves_dust /= nWaves',only_rank=rank,called_by='dust_read_dust_SEDs') 
                        stop ! stop the program 
                    end if
                    ! The main process checks if ndtypes = ndtypes_seds
                    ! If it is not the case the main process kills the run
                    if (ndtypes_seds  .ne. ndtypes) then
                        call IO_print_error_message('ndtypes_seds /= ndtypes',only_rank=rank,called_by='dust_read_dust_SEDs') 
                        stop ! stop the program 
                    end if
                end if
                !
                ! allocate arrays
                ! ISRF scalling factor
                allocate(ISRFBins(nISRFBins))
                ! IR SEDs
                allocate(dust_sed(nWaves,nISRFBins,2,ndtypes))
                !
                read(dust_sed_unit,*) ! skip wavelenght list
                read(dust_sed_unit,*) ISRFBins(1:nISRFBins)
                read(dust_sed_unit,*) ISRF_ref    ! ISRF reference value (G0 = 1) [erg/s/cm**2]
                ! convert in code unit [Lsun/kpc**2]
                ISRF_ref = ISRF_ref/L_Sun_erg_s*(kpc_in_cm)**2   
                !
                read(dust_sed_unit,*) ! skip a line #
                !
                ! read SEDs
                do i = 1, nISRFBins
                    do l = 1, nWaves
                        read(dust_sed_unit,*) dust_sed(l,i,1,1:ndtypes), dust_sed(l,i,2,1:ndtypes) ! Lsun/Msun
                    end do
                end do
                !
                if (main_process) then
                    ! Print some informations 
                    call IO_print_message('use', param_name=(/'nISRFBins                ','nWaves                   '/),&
                        int_param_val=(/nISRFBins,nWaves/))
                end if
                exit  ! quit do loop 
            end if
            !
            if (line(1:1) .eq. '#') then
                cycle ! header or something like this (skip)
            else
                call IO_print_error_message('Impossible to read the line',only_rank=rank,called_by='dust_read_dust_SEDs') 
                stop ! stop the program
            end if
        end do
        !
        ! set i_8mic and i_1000mic
        i = 1
        i_8mic = -1
        do while (Waves(i) .le. 1000.)
            if ((Waves(i) .gt. 8.0) .and. (i_8mic .lt. 0)) i_8mic = i-1
            i = i +1
        end do
        i_1000mic = i -1
        
    end if
    
2   close(abs_sca_unit)
    
    return  
  end subroutine dust_read_dust_SEDs
 
  !*****************************************************************************************************************
  
  subroutine dust_compute_effective_extinction(dust,tau_8mic)
  
    ! COMPUTE EFFECTIVE EXTINCTION 
    ! tau(lambda) and A(lambda) for different geometries
    
    implicit none

    integer(kind=4)                   :: gt                      ! dust type
    character(MAXPATHSIZE)            :: message                 ! a message to display
 
    real(kind=8)                      :: f_vsg                   ! fraction of vsg                     ! 
    real(kind=8)                      :: f_dust(3)
    real(kind=8)                      :: aspect
    real(kind=8),allocatable          :: eff_Ext(:), eff_Alb(:)  ! tmp arrays (extinction and albedo)
    real(kind=8),allocatable          :: alamb(:)
    real(kind=8),intent(out),optional :: tau_8mic                ! the extinction measured at 8 microns (used to select AGN SED template)

    type(dust_type),intent(inout)     :: dust                    ! a dust component

    !
    ! allocate dust%eff_ext
    if (.not. allocated(dust%eff_ext)) allocate(dust%eff_ext(nWaves))  ! create
    dust%eff_ext = 1.d0                                                ! init 
    
    if (dust%tau .le. 0.d0) return ! no extinction
    
    ! ***********************************************
    ! 1st STEP: BUILD THE EFFECTIVE EXTINCTION CURVE
    ! ***********************************************
    !
    ! take into accound dust composition (f_pah, f_vsg and f_bg,  [M_/M_dust])
    f_vsg = 1.d0 - (dust%f_pah + dust%f_bg)
    f_dust = (/dust%f_pah,f_vsg,dust%f_bg/)
    !
    ! allocate eff_Ext and eff_Alb
    allocate(eff_Ext(nWaves))  ! create
    eff_Ext(:) = 0.d0          ! init
    allocate(eff_Alb(nWaves))  ! create
    eff_Alb(:) = 0.d0          ! init
    !
    ! compute eff_Ext and eff_Alb
    ! Effective extinction is the sum of Abs and Sca for each dust component
    do gt = 1, ndtypes
        eff_Ext = eff_Ext + f_dust(gt)*(dust_Abs(gt,:) + dust_Sca(gt,:))
    end do
    ! check 
    if (minval(eff_Ext) .lt. 0.d0) then
        write(message,'(a,a,a,e10.3)') 'eff_Ext < 0 [ geom: ', trim(dust%geom), ' ] : ', minval(eff_Ext)
        call IO_print_error_message(message,only_rank=rank,called_by='dust_compute_effective_extinction') 
#ifdef PRINT_WARNING
! -------------------------------------------------        
        call IO_print_message('use',only_rank=rank,component='dust', &
               param_name = (/'tau                      ','dust%f_pah               ','dust%f_bg                ', &
                              'f_vsg                    '/), &
               real_param_val  = (/dust%tau,dust%f_pah,dust%f_bg,f_vsg/))                  
! -------------------------------------------------
#endif
! PRINT_WARNING
        stop ! stop the program
    end if  
    ! Effective albedo in the ratio of Sca onto the effective extinction
    do gt = 1, ndtypes
        eff_Alb = eff_Alb + f_dust(gt)*dust_Sca(gt,:)/eff_Ext
    end do
    if (minval(eff_Alb) .lt. 0.d0) then
        write(message,'(a,a,a,e10.3)') 'eff_Alb < 0 [ geom: ', trim(dust%geom), ' ] : ', minval(eff_Alb)
        call IO_print_error_message(message,only_rank=rank,called_by='dust_compute_effective_extinction') 
#ifdef PRINT_WARNING
! -------------------------------------------------        
        call IO_print_message('use',only_rank=rank,component='dust', &
               param_name = (/'tau                      ','f_pah                    ', &
                              'f_vsg                    ','f_bg                     '/), &
               real_param_val  = (/dust%tau,dust%f_pah,f_vsg,dust%f_bg/))                  
! -------------------------------------------------
#endif
! PRINT_WARNING
        stop ! stop the program
    end if  
    !
    ! ***********************************************
    ! 2d STEP: APPLY GEOMETRY
    ! ***********************************************
    !
    ! 
    allocate(alamb(nWaves))  ! create
    alamb(:) = 0.d0          ! init
    !
    select case (trim(dust%geom))
    case('slab','screen','BC','sdwich','disc','disk','clumps','clumpy','torus')
       !
       ! In all these cases, we use FUV depth
       ! We have to normilized eff_Ext to i_FUV
       eff_Ext = dust%tau*eff_Ext/eff_Ext(i_FUV)
       !
       ! treat the extreme extinction case
       aspect = max(1.d-4, abs(cos(dust%incl))) 
       !
       alamb = sqrt(1.d0 - eff_Alb)*eff_Ext/aspect
       !
       select case (trim(dust%geom))
       case ('slab','disc','disk')
            ! compute effective extinction
            dust%eff_ext = real((1.d0 - exp(-1.d0*alamb))/alamb,4) 
            ! 
       case ('screen','BC')
            ! compute effetive extinction
            dust%eff_ext = real(exp(-1.d0*alamb),4)
            !
       case ('sdwich')
            ! compute effetive extinction
            dust%eff_ext = real(2.5d-1 + 5.d-1*(1.d0 - exp(-1.d0*alamb))/alamb + 2.5d-1*exp(-1.d0*alamb),4)
            !
       case ('clumps','clumpy')
            ! compute effetive extinction
            ! we assume nclumps = 4.
            dust%eff_ext = real(exp(-4.d0*(1.d0 - exp(-1.d0*alamb))),4)
            !
       case ('torus')
            eff_Ext = real(5.d0*eff_Ext,4)
            !
            if (present(tau_8mic)) then
                tau_8mic = eff_Ext(i_8mic)
                return
            end if
       end select
        !
    case ('dwek','bulge')
        !
        ! In this case we use optical depth
        ! We have to normilized eff_Ext to end select(i_V)
        ! we also apply correction factor (Devriendt+99)
        eff_Ext = 2.619d0*eff_Ext/eff_Ext(i_V)
        !
        alamb = 1.d0
        where (eff_Ext .gt. 1.e-4)
            alamb = (3.d0/(4.d0*eff_Ext))*(1.d0 - (1.d0/(2.d0*eff_Ext**2.)) + (1.d0/eff_Ext + (1.d0/(2.d0*eff_Ext**2.)))*exp(-2.d0*eff_Ext))
        endwhere
        !
        if (minval(alamb) .le. 0.d0) then
            write(message,'(a,a,a,e10.3)') 'alamb <= 0 [ geom: ', trim(dust%geom), ' ] : ', minval(alamb)
            call IO_print_warning_message(message,only_rank=rank,called_by='dust_compute_effective_extinction') 
#ifdef PRINT_WARNING
! -------------------------------------------------        
            call IO_print_message('use',only_rank=rank,component='dust', &
                param_name = (/'tau                      ','f_pah                    ', &
                              'f_vsg                    ','f_bg                     '/), &
                real_param_val  = (/dust%tau,dust%f_pah,f_vsg,dust%f_bg/))                  
! -------------------------------------------------
#endif
! PRINT_WARNING
        end if
        !
        ! compute effective extinction
        dust%eff_ext = real(alamb/(1.d0 - eff_Alb + eff_Alb*alamb),4)
        !
    case default
        write(message,'(a,a,a)') 'Dust geometry distribution', trim(dust%geom), ' not defined'
        call IO_print_error_message(message,only_rank=rank,called_by='dust_compute_effective_extinction')    
        stop  ! stop the program
    end select
    !
    if (minval(dust%eff_Ext) .lt. 0.d0) then
        write(message,'(a,a,a,e10.3,a)') 'un-physical value of dust%eff_Ext [ geom: ', trim(dust%geom), ' ], [ min: ', minval(eff_Ext), ' ]'
        call IO_print_error_message(message,only_rank=rank,called_by='dust_compute_effective_extinction') 
#ifdef PRINT_WARNING
! -------------------------------------------------        
        call IO_print_message('use',only_rank=rank,component='dust', &
               param_name = (/'tau                      ','f_pah                    ', &
                              'f_vsg                    ','f_bg                     '/), &
               real_param_val  = (/dust%tau,dust%f_pah,f_vsg,dust%f_bg/))                  
! -------------------------------------------------
#endif
! PRINT_WARNING
        stop ! stop the program
    end if
      
    if (allocated(eff_Ext)) deallocate(eff_Ext)  ! erase
    if (allocated(eff_Alb)) deallocate(eff_Alb)  ! erase
    if (allocated(alamb))   deallocate(alamb)    ! erase

    return
  end subroutine dust_compute_effective_extinction
  
  !*****************************************************************************************************************
  
  subroutine dust_build_dust_spectrum(dust,ISRF,dspt)
  
    ! BUILD THE INFRARED SPECTRUM ASSOCIATED TO THE DUST COMPONENT
    
    implicit none
    
    integer(kind=4)                   :: iISRF,gt,iage    ! ISRF,dust type and dspt loop indexes
    
    character(*),intent(in),optional  :: dspt             ! dominated stellar population type (old young stars)
    character(MAXPATHSIZE)            :: message          ! a message to display
    
    real(kind=8),intent(in)           :: ISRF             ! ISRF scalling factor
    real(kind=8)                      :: f_vsg            ! fraction of vsg
    real(kind=8)                      :: f_dust(3)

    type(dust_type),intent(inout)     :: dust             ! a dust component
    
    !
    ! allocate dust%spectrum
    if (.not. allocated(dust%spectrum)) allocate(dust%spectrum(nWaves))  ! create
    dust%spectrum(:) = 0.d0                                              ! init
    
    if (dust%tau .le. 0.d0) return ! no dust
    !
    ! select the best SED template as function of the ISRF
    iISRF = 1
    do while ((ISRFBins(iISRF) < ISRF) .and. (iISRF < nISRFBins))
        iISRF = iISRF +1
    end do
    ! We select the SED template with the closest ISRFBins  but always lower than the current galaxy ISRF
    iISRF = max(1,iISRF -2)
    !
    ! take into accound dust composition (f_pah, f_vsg and f_bg,  [M_/M_dust])
    f_vsg = 1.d0 - (dust%f_pah + dust%f_bg)
    f_dust = (/dust%f_pah,f_vsg,dust%f_bg/)
    !
    select case (trim(dspt)) 
        case ('young','ystars','young stars')
            iage = 1
        case ('old','ostars','old stars')
            iage = 2
        case default
            write(message,'(a,a,a)') 'Keyword ', trim(dspt), ' not defined'
            call IO_print_error_message(message,only_rank=rank,called_by='dust_build_dust_spectrum')    
            stop  ! stop the program
    end select
     
    do gt = 1, ndtypes
        dust%spectrum = dust%spectrum + real(f_dust(gt)*dust_sed(:,iISRF,iage,gt),4)    ! Lsun/Msun/Md
    end do
    
    return
  end subroutine dust_build_dust_spectrum
  
! ------------------------------------------------- 
#endif
! LUMINOUS_PROCESSES

  !*****************************************************************************************************************
  !
  ! FUNCTIONS
  !
  !*****************************************************************************************************************
 
  function dust_mass(dust)
  
    ! RETURN THE MASS OF DUST
    
    implicit none
    
    real(kind=8)                :: dust_mass
    
    type(dust_type),intent(in)  :: dust
    
    dust_mass = dust%mass
    
    return
  end function
  
  !*****************************************************************************************************************  
  
  function dust_extinction(dust,gas,rc)
  
    ! RETURN THE FACE-ON OPTICAL DEPTH 
    ! Use a metallicity dependant prescription based on Boquien+13 Eq. 13
    ! This FUV-band extinction (for disc) is a function of both the hydrogen surface density and the oxygen abundance O/H
    ! The V-band extinction (for bulges) is a function of both the Hydrogen column density and the metallicity 
    
    implicit none
    
    character(MAXPATHSIZE)     :: message          ! a message 
    
    real(kind=8),intent(in)    :: rc               ! radius
    real(kind=8)               :: rc_in_cm         ! rc converted in cm
    real(kind=8)               :: M_H, N_H         ! hydrogen mass, column density
    real(kind=8)               :: Sig_H            ! Hydrogen surface density
    real(kind=8)               :: O_over_H         ! Oxygen relative abundance (compare to Hydrogen)
    real(kind=8)               :: dust_extinction  ! The FUV extinction tating into account metallicity dependance 
    
    type(dust_type),intent(in) :: dust             ! the dusc component
    type(gas_type),intent(in)  :: gas              ! the gas component
 
#ifdef PRINTALL 
    ! -------------------------------------------------
    call IO_print_message('dust_extinction',only_rank=rank,component='dust')
    ! -------------------------------------------------
#endif   
    
    dust_extinction = -1.d0 ! init

    if ((dust%f_pah .le. 0.d0) .or. (dust%f_bg .le. 0.d0)) return
    if (rc .le. 0.d0) return
    
    M_H = gas_mass(gas,component='H1')
    
    ! Oxygen relative abundance (compare to Hydrogen, in number)
    O_over_H = O_H(gas)
 
    select case (trim(dust%geom)) 
    case('slab','sdwich','disc','disk')
        !
        ! compute hydrogen column density
        ! the diffuse gas evolve in a thick disc 2 times larger than the stellar component 
        ! 50.0% of the mass of an exponential disc with a characteristic radius "rd" in enclosed in a radius r_50 = 1.68 x rd 
        ! If stars are distributed in the mid-plan of the disc, the light produced by stars only 'see' 
        ! in average half of the gas mass before to leave the disc
        Sig_H = 5.d-1*5.d-1*M_H*mass_code_unit_in_M_Sun/(pi*(2.d0*1.68d0*rc*1.d3)**2.)    ! in Msun/pc**2
        !
        ! tau_FUV
        dust_extinction = 1.926d0 + 5.1d-2*Sig_H
        dust_extinction = dust_extinction*1.d1**(9.47d-1*(12.d0+log10(O_over_H)-9.d0))
        !
    case('clumps','clumpy','BC','screen')
        !
        ! compute hydrogen column density
        Sig_H = M_H*mass_code_unit_in_M_Sun/(pi*(5.d-1*rc*1.d3)**2.)     ! in Msun/pc**2
        !
        ! tau_FUV
        dust_extinction = 1.364d0 + 2.9d-2*Sig_H
        dust_extinction = dust_extinction*1.d1**(4.27d-1*(12.d0+log10(O_over_H)-9.d0))
        !
    case('dwek','bulge','spheroid','torus')
        !
        ! compute hydrogen density
        rc_in_cm = rc*kpc_in_m*1.e2
        N_H = (5.d-1*M_H*mass_code_unit_in_kg)/(2.d0*mp*pi*(rc_in_cm)**2.) ! at/cm**2
        !
        ! WARNING tau_V
        dust_extinction = (gas_metallicity(gas)/Z_sun)**(1.6)*(N_H/2.1d21)
        !
    case default
        call IO_print_error_message('geometry not allowed',called_by = 'dust_extinction')  
        write(message,'(a,a,a)') 'Keyword: ', trim(dust%geom), ' unknown !'
        call IO_print_message(trim(message))
        stop ! stop the program
    end select  
    
    ! In all cases the reference depth is limited to tau = [0.01 ; 1000].
    dust_extinction = max(1.d-2,min(dust_extinction,1.d3))   

#ifdef PRINTALL 
    ! -------------------------------------------------
    call IO_print_message('dust_extinction ... done',only_rank=rank,component='dust')
    ! -------------------------------------------------
#endif  

    return
  end function dust_extinction
  
  !*****************************************************************************************************************
  
  function dust_AV(d)                    
  
    ! RETURN A(V) FROM THE EFFECTIVE EXTINCTION PROFILE OF A DUST COMPONENT 
    
    implicit none
    
    real(kind=4)                 :: dust_AV   ! the attenuation in the Visible band
     
    type(dust_type),intent(in)   :: d         ! A dust component
    
    dust_AV = -2.5*log10(d%eff_Ext(i_V))
    
    return
  end function dust_AV
  
  !*****************************************************************************************************************
  
  function dust_E_BV(d)    
  
    ! RETURN E(B-V) FROM THE EFFECTIVE EXTINCTION PROFILE OF A DUST COMPONENT 
    
    implicit none
    
    real(kind=4)                 :: dust_E_BV   ! the attenuation
     
    type(dust_type),intent(in)   :: d           ! A dust component
    
    dust_E_BV = -2.5*log10(d%eff_Ext(i_B)) + 2.5*log10(d%eff_Ext(i_V))
    
    return
  end function dust_E_BV
  
  !*****************************************************************************************************************
  !
  ! PRINTING PROCEDURES
  !
  !*****************************************************************************************************************

  subroutine dust_load_dust_data(dust,dust_data)

    ! CREATE THE DUST OUTPUT PROPERTIES LIST

    implicit none
    
    real(kind=8), intent(inout) :: dust_data(nb_dust_field)

    type(dust_type),intent(in)  :: dust    ! dust component

    ! For information
    ! 'f_pah   ','f_bg    ', 'tau     ','mZmH    '       

    dust_data = (/dust%f_pah,dust%f_bg,dust%tau,dust%mZ_mH/)

    return
  end subroutine dust_load_dust_data
  
  !***********************************************************************************************************

end module dust
