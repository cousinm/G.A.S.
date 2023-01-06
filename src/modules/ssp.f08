module ssp_mod

    !*****************************************************************************************************************
    !
    !  OVERVIEW
    !
    ! The ssp module defines all properties and procedures
    ! associated to single stellar populations
    ! A ssp is a stellar population for specific age/metalicity bin 
    !
    !*****************************************************************************************************************

    use parameters  ! Acces to global defintions and properties
    use model_mod   ! Acces to model parameters and integration scheme configuration
    use gas_mod     ! Acces to gas properties and procedures
    use status_mod  ! Acces to transfer rates status
    use config_mod  ! Acces to configurations parameters (path)

    implicit none

    public

    ! Single stellar population
    type ssp
        integer(kind=ikd)  :: iAge   ! Index of the Age bin in the sfh
        integer(kind=ikd)  :: jMet   ! Index of the Metalicity bin in the sfh
        real(kind=rkd)     :: mass   ! Mass of the ssp
        real(kind=rkd)     :: tform  ! Formation time
        real(kind=rkd)     :: avgAge ! Average age of the sp
    contains
        procedure    :: create => ssp_create      ! Create/Init a ssp
        procedure    :: mlr => ssp_mlr         ! Return the current sn ejecta rate
        procedure    :: transfer => ssp_transfer  ! Return the transfer rate from a age bin to the next age bin
        procedure    :: stevolve => ssp_stevolve  ! Evolve, for dt, the next step of the integration scheme
        procedure    :: evolve => ssp_evolve      !         according to internal/external input/output
        procedure    :: update => ssp_update
        procedure    :: solve => ssp_solve        ! Solve, for dt, the next step of the integration scheme
    end type ssp
    !
    ! Define ssp specific parameters
    integer(kind=ikd)           :: nAgeBins            ! Number of stellar age bins
    integer(kind=ikd)           :: nWaves              ! Number of wavelenghts used to defined a SSP SED

    real(kind=rkd)              :: spTimeStep          ! Minimal stellar evolution time step

    real(kind=rkd), allocatable :: ageBins(:)          ! Age bin values
    real(kind=rkd), allocatable :: lumBins(:, :)       ! Luminosity of ssp SEDs for a given age
                                                       !  (allows luminosity weight computation)
    real(kind=rkd), allocatable :: waves(:)            ! SED wavelengths
    real(kind=rkd), allocatable :: SNR(:, :)           ! SN Rates according to age and metalicity
    real(kind=rkd), allocatable :: spSED(:, :, :)      ! SSP SED

    type(gas), allocatable      :: MLR(:, :)           ! Gas ejection rates according to age and metlicity
    type(status), target        :: sspStatus           ! The current ssp status

    ! INTERFACE OPERATOR DECLARATIONS

    interface assignment (=)  ! allows to copy a ssp component by using the symbol '='
        module procedure ssp_copy
    end interface assignment (=)

contains

    !
    ! SUBROUTINES
    !
    ! **********************************
    subroutine ssp_create(this, iAge, jMet)

        ! Initialize a single stellar population

        implicit none

        integer(kind=ikd), intent(in)  :: iAge, jMet

        class(ssp)   :: this

        this%iAge   = iAge
        this%jMet   = jMet
        this%mass   = 0.d0
        this%avgAge = 0.d0
        this%tform  = 0.d0

    end subroutine ssp_create

    ! **********************************
    subroutine ssp_copy(ssp1, ssp2)

        ! Copy a ssp object ssp2 into the ssp object ssp1

        implicit none

        type(ssp), intent(in)      :: ssp2

        class(ssp), intent(inout)  :: ssp1

        ssp1%iAge   = ssp2%iAge
        ssp1%jMet   = ssp2%jMet
        ssp1%mass   = ssp2%mass
        ssp1%avgAge = ssp2%avgAge
        ssp1%tform  = ssp2%tform

    end subroutine ssp_copy

    ! **********************************
    subroutine ssp_init()

        ! Initialize the ssp module and dependancies
        !  Read stellar population properties (MLR, SNR, spSED)

        implicit none

        ! Read mass loss rates
        call ssp_read_MLR()

        ! Read SN Rates
        call ssp_read_SNR()

        ! Read SSP SED
        call ssp_read_SED()

    end subroutine ssp_init

    ! **********************************
    subroutine ssp_read_MLR

        ! Read the mass loss rates for a stellar population of given age and metallicity
        ! The model used is linked to the initial mass function used
    
        implicit none
        
        integer(kind=ikd)        :: nMetBins_    ! tmp nMetBins, 
                                                 !   allows to compared values read here with 
                                                 !   nMetBins data used/saved in the gas module
        integer(kind=ikd)        :: nElts_       ! tmp nElts
                                                 !   allows to compared values read here with 
                                                 !   nElts data used/saved in the gas module
        integer(kind=ikd)        :: i            ! loop index for stellar age
        integer(kind=ikd)        :: j            ! loop index for metalicity
        integer(kind=ikd)        :: e            ! loop index for elements
        
        character(MAXPATHSIZE)   :: filename
        character(MAXPATHSIZE)   :: message
        character(MAXPATHSIZE)   :: line
    
        ! Build the input filename, according to the initial mass function (IMF)
        write(filename,'(a,a,a,a)') trim(stellarPopPath), '/sp_MLR_[BC03]_', trim(IMF), '.in'
        write(message,'(a,a)') 'Load data from: ', trim(filename)
        call log_message(message)
        !
        ! Open the file
        open(unit=massLoss_unit, file=trim(filename), status='old')
        ! Read and load data
        do
            read(massLoss_unit, '(a)', end=2) line
            if (trim(line) .eq. '----') then
                !
                ! After the header lines:
                ! Read nAgeBins  : Number of stellar age bin used in the stellar population model
                ! Read nMetBins  : Number of metallicity bins used in the stellar population model
                ! Read nElts     : Number of elements followed by the stellar population model (C,N,O,Fe)
                read(massloss_unit,*) nAgeBins, nMetBins_, nElts_
                !
                ! Test nMetBins and nElts already saved in the gas module
                if (nMetBins_ .ne. nMetBins) then
                    call log_message('Inconsistency between metalicity bin counts ssp .vs. gas', &
                                     logLevel=LOG_ERROR, calledBy='sp_read_MLR')
                endif
                !
                if (nElts_ .ne. nElts) then
                    call log_message('Inconsistency between elements counts ssp .vs. gas', &
                                     logLevel=LOG_ERROR, calledBy='sp_read_MLR')
                end if
                !
                read(massLoss_unit, *) ! Skip the line with Metallicity bins, these data are already saved in the gas module
                !
                ! Skip the blank line
                read(massLoss_unit, *)
                !
                ! Allocate arrays
                ! Metalicity bins and element names are already saved in the gas module
                ! Allocate stellar age (Gyr) table corresponding to each step of the stellar model used
                allocate(ageBins(nAgeBins))
                ! 
                ! Allocate the massLossRate table,
                ! Storing ejecta rates associated to elements (H, He, C, N, O and Fe) [CU]
                allocate(MLR(nAgeBins, nMetBins))
                !
                do i = 1, nAgeBins
                    do j = 1, nMetBins
                        call MLR(i, j)%create() ! Init
                    end do                 
                    ! For each stellar age bin read:
                    !  - The age of the population,
                    !  - The mass loss rate 
                    !  - The yield of each elements
                    read(massLoss_unit,*) AgeBins(i), &
                                          (MLR(i, j)%mass, j=1, nMetBins), &
                                          (MLR(i, j)%mZ, j=1, nMetBins), &
                                          ((MLR(i, j)%elts(e), j=1, nMetBins), e=1, nElts)
                end do
                !       
                ! Define spTimeStep
                ! as the minimal step in the stellar population evolution model
                spTimeStep= minval(AgeBins(2 : nAgeBins) - AgeBins(: nAgeBins - 1))
                !
                exit
            end if
            if (line(1:1) .eq. '#') then
                cycle ! Header or something like this (skip)
            else
                write(message, '(a, a)') 'Impossible to read a line in ', trim(filename)
                call log_message(message, &
                                 logLevel=LOG_ERROR, &
                                 calledBy='sp_read_MLR')
            end if
        end do
2       close(massLoss_unit)
    
        return
    end subroutine ssp_read_MLR

    ! **********************************
    subroutine ssp_read_SNR()

        ! Read and load SN rates (SNII + SNIa)

        implicit none

        integer(kind=ikd)        :: nMetBins_    ! tmp nMetBins, 
                                                 !   allows to compared values read here with 
                                                 !   nMetBins data used/saved in the gas module
        integer(kind=ikd)        :: nAgeBins_    ! tmp nAgeBins
                                                 !   allows to compared values read here with 
                                                 !   nAgeBins data used previously in ssp_read_MLR
        integer(kind=ikd)        :: i            ! loop index under stellar age
        integer(kind=ikd)        :: j            ! loop index under metalicity

        character(MAXPATHSIZE)   :: filename
        character(MAXPATHSIZE)   :: message
        character(MAXPATHSIZE)   :: line

        real(kind=rkd)           :: aAgeBin      ! A local variable used to read the stellar age

        ! Build the input filename, according to the initial mass function (IMF)
        write(filename,'(a,a,a,a)') trim(stellarPopPath), '/sp_SNR_[BC03]_', trim(IMF), '.in'
        write(message,'(a,a)') 'Load data from: ', trim(filename)
        call log_message(message)
        !
        ! Open the library file
        open(unit=snRate_unit, file=trim(filename), status='old')
        ! Read and load data
        do
            read(snrate_unit, '(a)', end = 2) line
            if (trim(line) .eq. '----') then
                !
                ! After the header lines:
                ! Read nAgeBins  : Number of stellar age bin used in the stellar population model
                ! Read nMetBins  : Number of metallicity bins used in the stellar population model
                read(snrate_unit,*) nAgeBins_, nMetBins_
                !
                ! Test nMetBins, nAgeBins and nElts already saved previously
                if (nMetBins_ .ne. nMetBins) then
                    call log_message('Inconsistency between metalicity bin counts sp .vs. gas', &
                                     logLevel=LOG_ERROR, calledBy='sp_read_SNR')
                endif
                !
                if (nAgeBins_ .ne. nAgeBins) then
                    call log_message('Inconsistency between age counts', &
                                     logLevel=LOG_ERROR, calledBy='sp_read_SNR')
                end if
                read(snrate_unit,*) ! Skip the line with Metallicity bins, these data are already saved in the gas module
                !
                ! Skip the blank line
                read(snrate_unit, *)
                !
                ! Allocate arrays
                ! SNRates: Number of SN events per Gyr for a 1 Msun SSP
                allocate(SNR(nAgeBins, nMetBins))
                do i = 1, nAgeBins
                    ! For each stellar age bin read:
                    ! - The age of the stellar population
                    ! - The SN event rate
                    read(snrate_unit, *) aAgebin, (SNR(i, j), j=1, nMetBins)
                end do
                !
                ! Convert in code unit
                ! Initially SN rates are given [nb/Gyr] for a 1Msun SSP convert in Mass CU
                SNR = SNR / MSun2mass
                !
                exit
            end if
            if (line(1:1) .eq. '#') then
                cycle ! Header or something like this (skip)
            else
                write(message, '(a, a)') 'Impossible to read a line in ', trim(filename)
                call log_message(message, &
                                 logLevel=LOG_ERROR, &
                                 calledBy='sp_read_SNR')
            end if
        end do
2       close(snrate_unit)

    end subroutine ssp_read_SNR

    ! **********************************
    subroutine ssp_read_SED
    
        ! Read and load ssp spectral energy distribution.
        ! Each spectrum is given for a ssp stellar mass of 1 Msun
        ! Intensity is given in Lsun lamb * Ilamb
        
        implicit none
        
        integer(kind=4)          :: i            ! index loop under stellar age
        integer(kind=4)          :: j, l         ! index loop under wavelenght
        integer(kind=ikd)        :: nMetBins_    ! tmp nMetBins, 
                                                 !   allows to compared values read here with 
                                                 !   nMetBins data used/saved in the gas module
        integer(kind=ikd)        :: nAgeBins_    ! tmp nAgeBins
                                                 !   allows to compared values read here with 
                                                 !   nAgeBins data used previously in ssp_read_MLR
         
        character(MAXPATHSIZE)   :: filename
        character(MAXPATHSIZE)   :: message         ! a message  
        character(MAXPATHSIZE)   :: line

        ! Build the input filename, according to the initial mass function (IMF)
        write(filename,'(a,a,a,a)') trim(stellarPopPath), '/sp_SED_[BC03]_', trim(IMF), '.in'
        write(message,'(a,a)') 'Load data from: ', trim(filename)
        call log_message(message)

        ! Open the library file
        open(unit=spsed_unit, file=filename, status='old')
        ! Read and load data
        do
            read(spsed_unit, '(a)', end = 2) line
            if (trim(line) .eq. '----') then
                !
                ! After the header lines:
                ! Read nAgeBins  : Number of stellar age bin used in the stellar population model
                ! Read nMetBins  : Number of metallicity bins used in the stellar population model
                ! Read nwaves    : Number of wavelenght used to defined the SED of each SSP
                read(spsed_unit, *) nAgeBins_, nMetBins_, nWaves
                !
                ! Test nMetBins and nAgeBins already saved previously
                if (nMetBins_ .ne. nMetBins) then
                    call log_message('Inconsistency between metalicity bin counts ssp .vs. gas', &
                                     logLevel=LOG_ERROR, calledBy='sp_read_SED')
                endif
                !
                if (nAgeBins_ .ne. nAgeBins) then
                    call log_message('Inconsistency between age counts', &
                                     logLevel=LOG_ERROR, calledBy='sp_read_SED')
                end if
                !
                ! Allocate arrays
                allocate(waves(nWaves))
                allocate(lumBins(nAgeBins, nMetBins))
                allocate(spSED(nAgeBins, nMetBins, nWaves))
                !
                read(spsed_unit, *) ! Skip the line with ages bins, these data are already saved previously
                read(spsed_unit, *) ! Skip the line with metallicity bins, these data are already saved previously
                ! Read wavelenght
                read(spsed_unit, *) waves(1:nWaves)
                ! Read elementary luminosities
                read(spsed_unit, *) ((lumBins(i, j), i=1, nAgeBins), j=1, nMetBins)  ! [Lsun] for a 1 Msun SSP
                ! Convert in code unit
                ! Initially luminosities are given for a 1 Msun SSP convert in Mass CU
                LumBins = LumBins / MSun2mass
                !
                ! Skip the blank line
                read(spsed_unit, *)
                !
                ! Read stellar SED
                do i = 1, nAgeBins
                    do l = 1, nWaves
                        read(spsed_unit, *) spSED(i, 1:nMetBins, l) ! [Lsun]
                    end do
                end do
                !
                exit
            end if
            if (line(1:1) .eq. '#') then
                cycle ! header or something like this (skip)
            else
                write(message, '(a, a)') 'Impossible to read a line in ', trim(filename)
                call log_message(message, &
                                 logLevel=LOG_ERROR, &
                                 calledBy='ssp_read_SED')
            end if
        end do
        !
    2   close(spsed_unit)
        return
    end subroutine ssp_read_SED

    ! **********************************
    subroutine ssp_stevolve(this, dt, st, inRate, outRate, pStatus, pSsp)

        ! Compute current status of the ssp structure
        ! according to :
        !   - Internal evolution processes (SN gas ejection)
        !   - A constant input rate (due to external process: SFR) "inRate" and
        ! Then
        ! Apply the next step "st" of the complete integration scheme
        !
        ! In output: 
        !   - inRate is setled to the current transfer rate to the next age bin
        !   - outRate is setled to the current ejection rate (SN gas ejection)

        implicit none

        integer(kind=ikd), intent(in)  :: st       ! Step index of evolution scheme

        real(kind=rkd), intent(in)     :: dt       ! Evolution time-step

        type(gas), intent(inout)       :: inRate   ! The input rate
        type(gas), intent(out)         :: outRate  ! The output rate
        type(gas)                      :: trRate   ! The tranfer rate

        type(ssp), pointer             :: pSsp     ! Pointer to the curent evolved state
        type(status), pointer          :: pStatus  ! Pointer to the final evolution status
        type(status)                   :: aStatus   ! Current status

        class(ssp)                     :: this     ! The current ssp

        ! Compute current global status due to
        ! Internal process + external input
        ! For ssp, no external ouput are taken into account,
        !
        ! Set new current status from the latest stage
        ! Tranfer rate
        trRate = pSsp%transfer()
        ! Output rate
        outRate = pSsp%mlr()
        call aStatus%set(inRate, trRate, outRate)
        !
        ! Apply one more solver step and update the current intermediate stage
        pSsp = this%solve(st, dt, aStatus, pStatus)
        !
        ! Update (for the next age bin) the input rate
        ! to the transfert rate of the current bin
        inRate = trRate
        return

    end subroutine ssp_stevolve

    ! **********************************
    subroutine ssp_evolve(this, dt, inRate, pOutRate)

        ! Evolve, during dt, the current ssp structure 
        ! according to :
        ! - Internal evolution processes
        ! - A constant input rate (SFR) "inRate" and
        !
        ! In output of the subroutine, "outRate" is the wind/sn output rate

        implicit none

        integer(kind=ikd)           :: st         ! Step index of evolution scheme
        real(kind=rkd), intent(in)  :: dt         ! The ssp is evolve during dt

        type(gas), intent(inout)    :: inRate     ! The (dt-)constant input rate
        type(gas), pointer          :: pOutRate   ! The wind/sn output rate

        type(ssp), target           :: aSsp       ! a ssp matching the current state
        type(ssp), pointer          :: pSsp       ! Pointer to the current state
        type(status), pointer       :: pStatus    ! Pointer to the final status

        class(ssp)                  :: this       ! The current ssp

        ! Init intermediate status with the current status
        aSsp = this
        pSsp => aSsp
        ! Define complete corrected status
        call sspStatus%reset()
        pStatus => sspStatus
        do st = 1, nSolverStep
            call this%stevolve(dt, st, inRate, pOutRate, pStatus, pSsp)
        end do
        ! Save complete evolved state
        this = aSsp

    end subroutine ssp_evolve

    !
    ! FUNCTIONS
    !

    ! **********************************
    function ssp_mlr(this) result(rate)

        ! Return the current status of the ssp
        ! according to its internal evolution only

        type(gas)      :: rate     ! output rate

        class(ssp)     :: this

        ! Compute wind/sn ejecta
        rate = this%mass * MLR(this%iAge, this%jMet)

        return

    end function ssp_mlr

    ! **********************************
    function ssp_transfer(this) result(rate)

        ! Return the transfer rate from a bin to the next age bin

        real(kind=rkd) :: t_form
        real(kind=rkd) :: vRate

        type(gas)      :: rate     ! Transfer rate

        class(ssp)     :: this

        ! A transfert is done if the formation time associated to the bin
        ! is higher than the next age bin
        
        if (this%iAge < nAgeBins) then
            t_form = this%tform + ageBins(this%iAge)
            if (t_form > ageBins(this%iAge + 1)) then
                vRate = this%mass/t_form
            end if
        end if
        !
        ! The transfer rate gas object is build according 
        ! to the initial abundance associated to the current metallicity bin
        ! This is just to create a gas object.
        ! In the context of sp evolution only the total mass (gas%mass) is used
        ! to process evolution
        rate = vRate * initAbund(this%jMet)
        return

    end function ssp_transfer

    ! **********************************
    function ssp_update(this, dt, aStatus) result(aSsp)

        ! Update the ssp during dt according to its current status

        implicit none

        real(kind=rkd)              :: dt                ! The ssp is evolve during dt
        real(kind=rkd)              :: dmOut, dmIn, dmTr ! Output, input, and transfered masses
        real(kind=rkd)              :: mass

        type(ssp)                   :: aSsp              ! The evolved ssp
        type(status), intent(in)    :: aStatus           ! The current ssp status (input and output rate)

        class(ssp)                  :: this              ! The current ssp

        ! Copy of the current state of the ssp
        aSsp = this
        !
        ! Save initial mass
        mass = this%mass
        ! Compute mass variation
        dmIn = dt * aStatus%in%mass
        dmTr = dt * aStatus%tr%mass
        dmOut = dt * aStatus%out%mass
        !
        ! Test mass evolution
        if (mass - dmTr - dmOut + dmIn < real(0.d0, kind=rkd)) then
            call log_message('Current evolution step leads to negative mass', &
                             logLevel=LOG_ERROR, calledBy='ssp_update')
        end if
        !
        ! Update the average age of the ssp
        ! Incoming mass is homogeneously received during dt
        ! its average age is therefore ageBin(i) + dt/2.
        if (mass - dmTr - dmOut + dmIn > num_accuracy) then
            aSsp%avgAge = ((aSsp%avgAge + dt) * max(real(0.d0, kind=rkd), mass - dmOut - dmTr) &
                        + dmIn * (ageBins(aSsp%iAge) + dt / real(2.d0, kind=rkd))) &
                        / (mass - dmTr - dmOut + dmIn)
        else
            aSsp%avgAge = real(0.d0, kind=rkd)
        end if
        !
        ! Evolution
        aSsp%mass = aSsp%mass - dmTr - dmOut + dmIn
        !
        ! Update formation time of this single stellar population
        if (aSsp%mass > 0.d0) then
            aSsp%tform = aSsp%tform + dt
        end if
        return

    end function ssp_update

    ! **********************************
    function ssp_solve(this, it, dt, aStatus, pStatus) result(aSsp)

        ! Apply one solver step

        implicit none

        integer(kind=ikd), intent(in)     :: it       ! Current solver iteration

        real(kind=rkd)                    :: w
        real(kind=rkd), intent(in)        :: dt       ! Full time step

        type(ssp)                         :: aSsp     ! The new intermediate evolution stage

        type(status), intent(in)          :: aStatus  ! Current status
        type(status), intent(in), pointer :: pStatus  ! Pointer to the final status

        class(ssp)                        :: this     ! The current ssp

        select case (trim(solver))
        case ('RK4')
            ! Range-Kutta 4th order
            select case (it)
            case (1)
                ! Reset corrected status with status 1
                w = real(1.d0, kind=rkd)/real(6.d0, kind=rkd)
                call pStatus%resetFrom(aStatus, w)
                ! Update by dt/2. using status 1
                aSsp = this%update(dt/real(2.d0, kind=rkd), aStatus)
            case (2, 3)
                ! Update corrected status with status 2 or 3
                w = real(2.d0, kind=rkd)/real(6.d0, kind=rkd)
                call pStatus%updateFrom(aStatus, w)
                ! Update by dt/2. using status 2
                aSsp = this%update(dt/real(2.d0, kind=rkd), aStatus)
            case (4)
                ! Final step
                ! Update corrected status with status 4
                w = real(1.d0, kind=rkd)/real(6.d0, kind=rkd)
                call pStatus%updateFrom(aStatus, w)
                ! Update by dt using final status (pStatus)
                aSsp = this%update(dt, pStatus)
            end select
        case default
            ! Range-Kutta 2d order
            select case (it)
            case (1)
                ! Reset corrected status with status 1
                w = real(1.d0, kind=rkd)
                call pStatus%resetFrom(aStatus, w)
                ! Update by dt/2. using status 1
                aSsp = this%update(dt/real(2.d0, kind=rkd), aStatus)
            case (2)
                ! Reset corrected status with status 2
                w = real(1.d0, kind=rkd)
                call pStatus%resetFrom(aStatus, w)
                ! Update by dt using complete corrected pStatus
                aSsp = this%update(dt, pStatus)
            end select
        end select

    end function ssp_solve

end module ssp_mod