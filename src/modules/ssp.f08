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
        type(ssp), pointer :: pSsp   ! Pointer to the reference
    contains
        procedure    :: create => ssp_create      ! Create/Init a ssp
        procedure    :: copy => ssp_copy          ! Copy a spp
        procedure    :: setAsRef => ssp_setAsRef  ! Set pointer to reference
        procedure    :: mlr => ssp_mlr            ! Return the insta,taneous SN ejecta rate
        procedure    :: snp => ssp_snp            ! Return the instantaneous SN power
        procedure    :: transfer => ssp_transfer  ! Return the transfer rate from a age bin to the next age bin
        procedure    :: ststatus => ssp_ststatus  ! Get the current status of the ssp
        procedure    :: stevolve => ssp_stevolve  ! Apply a integrator scheme step
        procedure    :: dtoptim => ssp_dtoptim    ! Get the optimal time-step
        procedure    :: update => ssp_update      ! Update a ssp
        procedure    :: evolve => ssp_evolve      ! Evolve the scale, for dt 
                                                  ! according to internal/external input/output
    end type ssp
    !
    ! Define specific constants
    real(kind=rkd), parameter   :: Esn = 1.d51 * erg2J * J2EnergCU ! From erg to J to CU
    real(kind=rkd), parameter   :: ssp_mass_accuracy = real(1.d-19, kind=rkd)
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
    
    ! STATUS
    type(status), target        :: sspStatus           ! The current status
    type(status), target        :: sspFinalStatus      ! The final status

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

        ! Init or create fields
        this%iAge   = iAge
        this%jMet   = jMet
        this%mass   = 0.d0
        this%avgAge = 0.d0
        this%tform  = 0.d0
        ! Set pointer
        this%pSsp => null()

    end subroutine ssp_create

    ! **********************************
    subroutine ssp_copy(this, aSsp)

        ! Copy in "this" the ssp object aSsp
        ! Pointer to reference IS NOT COPIED

        implicit none

        type(ssp), intent(in)      :: aSsp

        class(ssp), intent(inout)  :: this

        this%iAge   = aSsp%iAge
        this%jMet   = aSsp%jMet
        this%mass   = aSsp%mass
        this%avgAge = aSsp%avgAge
        this%tform  = aSsp%tform

    end subroutine ssp_copy

    ! **********************************
    subroutine ssp_setAsRef(this, aSsp)

        ! Set pointer to the new reference
    
        implicit none

        type(ssp), intent(in), target  :: aSsp

        class(ssp), intent(inout)      :: this

        ! Pointer to new refence
        this%pSsp => aSsp

    end subroutine ssp_setAsRef

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
        write(filename,'(a, a, a, a)') trim(stellarPopPath), '/sp_MLR_[BC03]_', trim(IMF), '.in'
        write(message,'(a, a)') 'Load data from: ', trim(filename)
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
                SNR = SNR / MSun2massCU
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
                LumBins = LumBins / MSun2massCU
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
    subroutine ssp_evolve(this, dt, inRate)

        ! Evolve, during dt the ssp structure
        ! according to :
        !   - Internal evolution processes (mass transfert between age bins)
        !   - A constant input rate (due to external process: SFR) "inRate" and
        !
        ! In output "outRate" is the current ejection rate (SN gas ejection)

        implicit none

        integer(kind=ikd)             :: st         ! Step index of evolution scheme
        real(kind=rkd), intent(inout) :: dt         ! The ssp is evolve during dt

        type(gas), intent(in)         :: inRate     ! The (dt-)constant input rate

        type(ssp), target             :: aSsp       ! a ssp matching the current state
        type(ssp), pointer            :: pSsp       ! Pointer to the current state
        
        type(status), pointer         :: pStatus       ! Pointer to the "current" status
        type(status), pointer         :: pFinalStatus  ! Pointer to the "final" status

        class(ssp)                    :: this       ! The current ssp

        ! Init intermediate status with the current status
        aSsp = this
        ! Pointer to status
        pFinalStatus => sspFinalStatus
        ! Pointer to current state
        pSsp => aSsp
        !
        ! Loop over integration steps
        do st = 1, nSolverStep
            !
            ! Get the current status of the ssp and update "final" status
            pStatus => sspStatus
            call pSsp%ststatus(inRate, pStatus)
            !
            ! Update final status
            call pFinalStatus%update(st, pStatus)
            !
            if (solver_isFinalStep(st)) then
                ! Switch to final status
                pStatus => sspFinalStatus
            end if
            !
            ! Adaptative time-step
            dt = this%dtoptim(st, dt, pStatus)
            !
            ! Get new state
            call this%stevolve(st, dt, pStatus, pSsp)
        end do
        ! Save complete evolved state
        this = aSsp

    end subroutine ssp_evolve

    ! **********************************
    subroutine ssp_ststatus(this, inRate, pStatus)

        ! Compute current status of the ssp structure
        ! according to :
        ! - Internal evolution processes "trRate"
        ! - A constant input rate (due to external process) "inRate"

        implicit none

        type(gas), intent(in)         :: inRate     ! The (dt-)constant input rate
        type(gas)                     :: outRate    ! The output rate (ejecta)
        type(gas)                     :: trRate     ! Internal transfer rate of the ssp

        type(status)                  :: pStatus    ! Pointer to the "current" status

        class(ssp)                    :: this       ! The current ssp

        ! Internal transfer rate
        trRate = this%transfer()
        !
        ! Ejecta rate
        outRate = this%mlr()
        !
        ! Set status
        call pStatus%set(inRate, trRate, outRate)

    end subroutine ssp_ststatus

    ! **********************************
    subroutine ssp_stevolve(this, st, dt, pStatus, pSsp)

        ! Apply a new integrator step according to "aStatus"

        implicit none

        integer(kind=ikd), intent(in)  :: st       ! Step index of evolution scheme

        real(kind=rkd), intent(in)     :: dt       ! Evolution time-step
        real(kind=rkd)                 :: wdt

        type(ssp), pointer             :: pSsp     ! Pointer to the intermediate state

        type(status), pointer          :: pStatus  ! Pointer to the final evolution status

        class(ssp)                     :: this     ! The current ssp

        ! Compute new state
        wdt = solver_wdt(st) 
        call this%update(wdt * dt, pStatus, pSsp)

    end subroutine ssp_stevolve

    ! **********************************
    subroutine ssp_update(this, dt, pStatus, pSsp)

        ! Update "this" during dt to the next intermediate state (given by: pSsp) 
        ! according to a current status "pStatus

        implicit none

        real(kind=rkd)              :: dt                ! The ssp is evolve during dt
        real(kind=rkd)              :: dmOut, dmIn, dmTr ! Output, input, and transfered masses
        real(kind=rkd)              :: mass, finalMass, sspMass

        type(status), pointer       :: pStatus           ! The current ssp status (input and output rate)

        type(ssp), pointer          :: pSsp              ! Pointer to the intermediate state

        class(ssp)                  :: this              ! The current ssp

        ! Save current mass
        sspMass = pSsp%mass
        ! Copy of the current state of the ssp
        pSsp = this
        !
        ! Save ref mass
        mass = this%mass
        ! Compute mass variation
        dmIn = dt * pStatus%in%mass
        dmTr = dt * pStatus%tr%mass
        dmOut = dt * pStatus%out%mass
        !
        ! Test very low final mass (lead to numerical precision)
        finalMass = mass - dmTr - dmOut + dmIn
        if (mass > real(0.d0, kind=rkd) & 
            .and. abs(finalMass) >= real(0.d0, kind=rkd) &
            .and. abs(finalMass) < ssp_mass_accuracy) then
            ! The final will be too low,
            ! Reset the ssp, set pointer to reference and return
            call pSsp%create(this%iAge, this%jMet)
            call pSsp%setAsRef(this)
            return
        end if
        !
        ! From this point, mass is sufficient
        ! Test mass evolution
        if (finalMass < real(0.d0, kind=rkd)) then
            call log_message('Current evolution step leads to negative mass', &
                             logLevel=LOG_ERROR, calledBy='ssp_update')
        end if
        !
        ! Update the average age of the ssp
        ! Incoming mass is homogeneously received during dt
        ! its average age is therefore ageBin(i) + dt/2.
        if (finalMass > real(0.d0, kind=rkd)) then
            pSsp%avgAge = ((pSsp%avgAge + dt) * max(real(0.d0, kind=rkd), mass - dmOut - dmTr) &
                            + dmIn * (ageBins(pSsp%iAge) + dt / real(2.d0, kind=rkd))) &
                            / (mass - dmTr - dmOut + dmIn)
            if (this%iAge < nAgeBins) pSsp%avgAge = min(pSsp%avgAge, ageBins(this%iAge + 1))
        end if
        !
        ! Evolution
        pSsp%mass = pSsp%mass - dmTr - dmOut + dmIn
        !
        ! Update formation time of this single stellar population
        if (pSsp%mass > 0.d0) then
            pSsp%tform = pSsp%tform + dt
        end if

    end subroutine ssp_update

    !
    ! FUNCTIONS
    !

    ! **********************************
    function ssp_mlr(this) result(rate)

        ! Return the instantaneous ejecta rate

        type(gas)      :: rate     ! output rate

        class(ssp)     :: this

        ! Compute wind/sn ejecta
        ! call rate%create()
        rate = this%mass * MLR(this%iAge, this%jMet)

    end function ssp_mlr

    ! **********************************
    function ssp_snp(this) result(power)

        ! Return the instantaneous SN disruption power

        real(kind=rkd)  :: power     ! output rate

        class(ssp)      :: this

        ! Compute SN power
        power = this%mass * SNR(this%iAge, this%jMet)*Esn   ! [CU]

    end function ssp_snp

    ! **********************************
    function ssp_transfer(this) result(rate)

        ! Return the transfer rate from a bin to the next age bin

        real(kind=rkd) :: t_tr, t_avg
        real(kind=rkd) :: vRate

        type(gas)      :: rate     ! Transfer rate

        class(ssp)     :: this

        ! A transfert is done if the formation time associated to the bin
        ! is higher than the next age bin
        vRate = real(0.d0, kind=rkd)
        if (this%iAge < nAgeBins) then
            if (ageBins(this%iAge) + this%tform > ageBins(this%iAge + 1)) then
                t_avg = (ageBins(this%iAge + 1) + ageBins(this%iAge)) / real(2.d0, kind=rkd)
                t_tr = abs((ageBins(this%iAge + 1) - ageBins(this%iAge)) - real(2.d0, kind=rkd)*abs(this%avgAge - t_avg))
                vRate = this%mass/max(t_tr, real(2./3., kind=rkd)*solver_dt)
            end if
        end if
        !
        ! The transfer rate gas object is build according 
        ! to the initial abundance associated to the current metallicity bin
        ! This is just to create a gas object.
        ! In the context of a ssp evolution only the total mass (gas%mass) is used
        ! to process evolution
        rate = vRate * initAbund(this%jMet)
        return

    end function ssp_transfer

    ! **********************************
    function ssp_dtoptim(this, st, dt, pStatus) result(adt)

        ! Return the optimal time-step accoring to current status
        ! and mass

        implicit none

        integer(kind=ikd), intent(in) :: st

        real(kind=rkd), intent(in)    :: dt       ! time-step
        real(kind=rkd)                :: adt

        type(status), pointer         :: pStatus  ! Pointer to the current status

        class(ssp)                    :: this     ! The current state

        adt = pStatus%dtMax(st, dt, this%mass)

    end function ssp_dtoptim

end module ssp_mod