module sp_mod

    !*****************************************************************************************************************
    !
    !  OVERVIEW
    !
    ! The sp module defines all properties and procedures
    ! associated to stellar populations
    ! This structure allows to "follow" stellar evoltion process
    !    age/metalicity and build stellar population spectra
    !
    !*****************************************************************************************************************

    use parameters  ! Acces to global defintions and properties
    use model_mod   ! Acces to model parameters and integration scheme configuration
    use gas_mod     ! Acces to gas properties and procedures
    use config_mod  ! Acces to configurations parameters (path)

    implicit none

    public

    ! A complete stellar population is built from a set of single stellar population
    !   A ssp evolved in a set of age and metalicity bin
    !
    ! Single stellar population
    type ssp
        real(kind=rkd)     :: mass   ! Mass of the ssp
        real(kind=rkd)     :: tform  ! Formation time
        real(kind=rkd)     :: avgAge ! Average age of the sp
    contains
        procedure    :: create => ssp_create  ! Create/Init a ssp
        procedure    :: evolve => ssp_evolve  ! Evolve a sap during dt
    end type ssp

    ! Complete stellar population
    type sp
        real(kind=rkd)           :: mass          ! Total mass of the stellar population
        real(kind=rkd)           :: mAge          ! Mass-weighted average age
        real(kind=rkd)           :: lAge          ! Luminosity-weighted average age
        real(kind=rkd)           :: mZ            ! Mass-weighted average metelicity
        real(kind=rkd)           :: lZ            ! Luminosity-weighted average metalicity
        type(ssp), allocatable   :: sfh(:, :)     ! Star Formation History
    contains
        procedure  :: create => sp_create         ! Create a complete stellar population
        procedure  :: delete => sp_delete         ! Delete a stellar population structure
        procedure  :: copy => sp_copy             ! Copy a stellar population structure
        procedure  :: isCreated => sp_isCreated   ! Test if sp object is already created
        procedure  :: evolve => sp_evolve         ! Evolve a stellar population structure by dt
        procedure  :: status => sp_status         ! Return the current status (outRates) of a stellar population
        procedure  :: transfer => sp_transfer     ! Transfer mass from an ageBin to the next one
        procedure  :: solve => sp_solve           ! Solve the input output system during dt
    end type sp

    ! Define ssp and sp specific parameters
    integer(kind=ikd)           :: nAgeBins            ! Number of stellar are bins

    real(kind=rkd), allocatable :: ageBins(:)          ! Age bin values
    real(kind=rkd)              :: spTimeStep          ! Minimal stellar evolution time step
    real(kind=rkd), allocatable :: SNRates(:, :)       ! SN Rates according to age and metalicity
    real(kind=rkd), allocatable :: gas2sp(:, :)        ! gas to sp matrix conversion, based on abundancies

    type(gas), allocatable    :: massLossRates(:, :) ! Ejection rates according to age and metlicity

    ! INTERFACE OPERATOR DECLARATIONS

    interface assignment (=)  ! allows to copy a sp component by using the symbol '='
        module procedure sp_copy, ssp_copy
    end interface assignment (=)

contains

    !
    ! SUBROUTINES
    !

    ! **********************************
    subroutine ssp_create(this)

        ! Initialize a single stellar population

        implicit none

        class(ssp)   :: this

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

        ssp1%mass   = ssp2%mass
        ssp1%avgAge = ssp2%avgAge
        ssp1%tform  = ssp1%tform

    end subroutine ssp_copy

    ! **********************************
    subroutine sp_init()

        ! Initialize the sp module and dependancies
        !  Read stellar population properties (mass loss rates, SN rates, spectrum)

        implicit none

        ! Read mass loss rates
        call sp_read_mass_loss_rates()

        ! Read SN Rates
        call sp_read_SN_rates()

        ! Build gas to stellar population conversion matrix
        call sp_build_gas2sp_matrix()

    end subroutine sp_init

    ! **********************************
    subroutine sp_read_mass_loss_rates

        ! Read the mass loss rates for a stellar population of given age and metallicity
        ! The model used is linked to the initial mass function used
    
        implicit none
        
        integer(kind=ikd)          :: nMetBins_    ! tmp nMetBins, 
                                                 !   allows to compared values read here with 
                                                 !   nMetBins data used/saved in the gas module
        integer(kind=ikd)          :: nElts_       ! tmp nElts
                                                 !   allows to compared values read here with 
                                                 !   nElts data used/saved in the gas module
        integer(kind=ikd)          :: i            ! loop index for stellar age
        integer(kind=ikd)          :: j            ! loop index for metalicity
        integer(kind=ikd)          :: e            ! loop index for elements
        
        character(MAXPATHSIZE)   :: filename
        character(MAXPATHSIZE)   :: message
        character(MAXPATHSIZE)   :: line
    
        ! Build the input filename, according to the initial mass function (IMF)
        write(filename,'(a,a,a,a)') trim(stellarPopPath), '/sp_mass_loss_rates[BC03]_', trim(IMF), '.in'
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
                    call log_message('Inconsistency between metalicity bin counts st .vs. gas', &
                                     logLevel=LOG_ERROR, calledBy='sp_read_mass_loss_rates')
                endif
                !
                if (nElts_ .ne. nElts) then
                    call log_message('Inconsistency between elements counts st .vs. gas', &
                                     logLevel=LOG_ERROR, calledBy='sp_read_mass_loss_rates')
                end if
                !
                read(massLoss_unit,*) ! Skip the line with Metallicity bins, these data are already saved in the gas module
                !
                ! Allocate arrays
                ! Metalicity bins and element names are already saved in the gas module
                ! Allocate stellar age (Gyr) table corresponding to each step of the stellar model used
                allocate(ageBins(nAgeBins))
                ! 
                ! Allocate the massLossRate table,
                ! Storing ejecta rates associated to elements (H, He, C, N, O and Fe) [CU]
                allocate(massLossRates(nAgeBins, nMetBins))
                !
                do i = 1, nAgeBins
                    do j = 1, nMetBins
                        call massLossRates(i, j)%create() ! Init
                    end do                 
                    ! For each stellar age bin read:
                    !  - The age of the population,
                    !  - The mass loss rate 
                    !  - The yield of each elements
                    read(massLoss_unit,*) AgeBins(i), &
                                          (massLossRates(i, j)%mass, j=1, nMetBins), &
                                          (massLossRates(i, j)%mZ, j=1, nMetBins), &
                                          ((massLossRates(i, j)%elts(e), j=1, nMetBins), e=1, nElts)
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
                                 calledBy='sp_read_mass_loss_rates')
            end if
        end do
2       close(massLoss_unit)
    
        return
    end subroutine sp_read_mass_loss_rates

    ! **********************************
    subroutine sp_read_SN_rates()

        implicit none

        integer(kind=ikd)          :: nMetBins_    ! tmp nMetBins, 
                                                 !   allows to compared values read here with 
                                                 !   nMetBins data used/saved in the gas module
        integer(kind=ikd)          :: nAgeBins_    ! tmp nAgeBins
                                                 !   allows to compared values read here with 
                                                 !   nAgeBins data used/saved in the gas module
        integer(kind=ikd)          :: i            ! loop index under stellar age
        integer(kind=ikd)          :: j            ! loop index under metalicity

        character(MAXPATHSIZE)   :: filename
        character(MAXPATHSIZE)   :: message
        character(MAXPATHSIZE)   :: line

        real(kind=rkd)             :: aAgeBin      ! a local variable used to read the stellar age

        ! Build the input filename, according to the initial mass function (IMF)
        write(filename,'(a,a,a,a)') trim(stellarPopPath), '/sp_SN_rates[BC03]_', trim(IMF), '.in'
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
                    call log_message('Inconsistency between metalicity bin counts st .vs. gas', &
                                     logLevel=LOG_ERROR, calledBy='sp_read_SN_rates')
                endif
                !
                if (nAgeBins_ .ne. nAgeBins) then
                    call log_message('Inconsistency between age counts st .vs. gas', &
                                     logLevel=LOG_ERROR, calledBy='sp_read_SN_rates')
                end if
                read(snrate_unit,*) ! skip the line with Metallicity bins, these data are already saved in the gas module
                !
                ! allocate arrays  
                ! SNRates: Number of SN events per Gyr and per Msun [nb/Gyr/Msun]
                allocate(SNRates(nAgeBins, nMetBins))
                do i = 1, nAgeBins
                    ! For each stellar age bin read:
                    ! - The age of the stellar population
                    ! - The SN event rate
                    read(snrate_unit,*) aAgebin, (SNRates(i, j), j=1, nMetBins)
                end do
                !
                ! Convert in code unit
                ! Initially SN rates are given [nb/Gyr/Msun] convert in Mass CU
                SNRates(:,:) = SNRates(:, :)/MSun2mass
                !
                exit
            end if
            if (line(1:1) .eq. '#') then
                cycle ! Header or something like this (skip)
            else
                write(message, '(a, a)') 'Impossible to read a line in ', trim(filename)
                call log_message(message, &
                                 logLevel=LOG_ERROR, &
                                 calledBy='sp_read_SN_rates')
            end if
        end do
2       close(snrate_unit)

    end subroutine sp_read_SN_rates

    ! **********************************
    subroutine sp_build_gas2sp_matrix()

        ! Build the conversion matrix from gas component to stellar population

        implicit none

        integer(kind=ikd)     :: i, j

        ! Allocate the matrix
        allocate(gas2sp(nMetBins, nMetBins))

        do i = 1, nMetBins
            do j = 1, nMetBins
                if (i == 1) then
                    gas2sp(i, j) = 1.
                else
                    gas2sp(i, j) = initAbund(j)%elts(i-1)
                end if
            end do
        end do

    end subroutine

    ! **********************************
    subroutine sp_create(this)

        ! Create a complete stellar population

        implicit none

        integer(kind=ikd)       :: i ! Age loop index
        integer(kind=ikd)       :: j ! Metalicity loop index

        class(sp)             :: this

        this%mass = 0.d0   ! Total mass of the stellar population
        this%mAge = 0.d0   ! Mass-weighted average age
        this%lAge = 0.d0   ! Luminosity-weighted average age
        this%mZ = 0.d0     ! Mass-weighted average metelicity
        this%lZ = 0.d0     ! Luminosity-weighted average metalicity

        allocate(this%sfh(nAgeBins, nMetBins)) ! Star Formation History
        do i = 1, nAgeBins
            do j = 1, nMetBins
                call this%sfh(i, j)%create()
            end do
        end do
    end subroutine sp_create

    ! **********************************
    subroutine sp_delete(this)

        ! Delete a complete stellar population

        implicit none

        class(sp)     :: this

        this%mass = 0.d0   ! Total mass of the stellar population
        this%mAge = 0.d0   ! Mass-weighted average age
        this%lAge = 0.d0   ! Luminosity-weighted average age
        this%mZ = 0.d0     ! Mass-weighted average metelicity
        this%lZ = 0.d0     ! Luminosity-weighted average metalicity

        if (allocated(this%sfh)) deallocate(this%sfh)
    end subroutine sp_delete

    ! **********************************
    subroutine sp_copy(sp1, sp2)

        ! Copy the sp object sp2 into the sp object sp1

        implicit none

        integer(kind=ikd)            i, j  ! loop indexes

        type(sp), intent(in)       :: sp2

        class(sp), intent(inout)   :: sp1

        if (.not. sp1%isCreated()) then
            ! Create g1
            call sp1%create()
        end if

        ! Copy fields
        sp1%mass = sp2%mass   ! Total mass of the stellar population
        sp1%mAge = sp2%mAge   ! Mass-weighted average age
        sp1%lAge = sp2%lAge   ! Luminosity-weighted average age
        sp1%mZ = sp2%mZ       ! Mass-weighted average metelicity
        sp1%lZ = sp2%lZ       ! Luminosity-weighted average metalicity

        ! Copy star formation history
        do i = 1, nAgeBins
            do j = 1, nMetBins
                sp1%sfh(i, j) = sp2%sfh(i, j)
            end do
        end do

    end subroutine sp_copy

    ! **********************************
    subroutine sp_transfer(this, iAge, iMet)

        ! Transfert
        implicit none

        integer(kind=ikd), intent(in)  :: iAge  ! Age bin index of the ssp
        integer(kind=ikd), intent(in)  :: iMet  ! Metalicity bin index of the ssp

        real(kind=rkd)                 :: m     ! Mass in the bin
        real(kind=rkd)                 :: dm    ! Mass transfert form a bin to an other
        real(kind=rkd)                 :: dt    ! time spend over the current age bin

        class(sp)                      :: this  ! The current ssp

        if (this%sfh(iAge, iMet)%tform > ageBins(iAge) .and. &
            iAge < nAgeBins) then
            !
            ! Transfet is needed and possible
            dt = this%sfh(iAge, iMet)%tform - ageBins(iAge)
            ! Assume homogeneous distribution of the mass
            ! TO DO, mass distribution scaling by "avgAge"
            m = this%sfh(iAge, iMet)%mass
            dm = dt / this%sfh(iAge, iMet)%tform * m
            !
            ! Mass is substract from the current bin
            ! Average age of the current age bin is impacted
            this%sfh(iAge, iMet)%avgAge = (this%sfh(iAge, iMet)%avgAge * m - &
                                           dm * (ageBins(iAge) + dt/2.d0)) / (m - dm)
            this%sfh(iAge, iMet)%mass = this%sfh(iAge, iMet)%mass - dm
            !
            ! And add to the next one
            ! Average age of the next age bin is also impacted
            m = this%sfh(iAge + 1, iMet)%mass
            this%sfh(iAge + 1, iMet)%avgAge = (this%sfh(iAge + 1, iMet)%avgAge * m + &
                                               dm * (ageBins(iAge) + dt/2.d0)) / (m + dm)
            this%sfh(iAge + 1, iMet)%mass = this%sfh(iAge + 1, iMet)%mass + dm
        end if
    end subroutine sp_transfer

    ! **********************************
    subroutine sp_solve(this, dt, inRate)

        ! Solve the evolution of the current stellar population
        !  according to a constant external input rate "inRate"

        implicit none

        integer(kind=ikd)          :: i, j

        real(kind=rkd), intent(in) :: dt            ! The gsh is evolved during dt

        type(gas), intent(in)    :: inRate          ! The (dt-)constant input rate
        type(gas), allocatable   :: outRates(:, :)  ! The corrected set of output rates
        type(gas), allocatable   :: outRates1(:, :) ! Intermediates set of output rates
        type(gas), allocatable   :: outRates2(:, :)
        type(gas), allocatable   :: outRates3(:, :)
        type(gas), allocatable   :: outRates4(:, :)

        type(sp)                 :: sp_tmp   ! Intermediate states of the stellar population

        class(sp)                :: this     ! The current stellar population

        ! Init OutRates
        allocate(outRates(nAgeBins, nMetBins))
        !
        ! Get curent status
        outRates1 = this%status()
        ! Performed evolution for dt/2
        sp_tmp = this%evolve(dt/2.d0, inRate, outRates1)
        !
        select case (trim(intScheme))
            !
            case ('RK4')
                ! Range-Kutta 4th order
                ! Get curent status
                outRates2 = sp_tmp%status()
                ! Get second intermediate state
                sp_tmp = this%evolve(dt/2.d0, inRate, outRates2)
                outRates3 = sp_tmp%status()
                ! Get last intermediate state
                sp_tmp = this%evolve(dt, inRate, outRates3)
                outRates4 = sp_tmp%status()
                ! Perform complete (corrected) evolution
                do i = 1, nAgeBins
                    do j = 1, nMetBins
                        outRates(i, j) = real(1.d0/6.d0,kind=rkd)*(outRates1(i, j) + &
                                            real(2.d0, kind=rkd)*outRates2(i, j) + &
                                            real(2.d0, kind=rkd)*outRates3(i, j) + &
                                            outRates4(i, j))
                    end do
                end do
                ! Final evolution with corrected output rates
                this = this%evolve(dt, inRate, outRates)
            case default
                ! Range-Kutta 2d order
                ! Get intermedite status
                outRates = sp_tmp%status()
                ! Perform complete evolution
                this = this%evolve(dt, inRate, outRates)
        end select
        !
        ! Delete temporary gas object
        do i = 1, nAgeBins
            do j = 1, nMetBins
                call outRates(i, j)%delete()
                call outRates1(i, j)%delete()
                if (allocated(outRates2)) call outRates2(i, j)%delete()
                if (allocated(outRates3)) call outRates3(i, j)%delete()
                if (allocated(outRates4)) call outRates4(i, j)%delete()
            end do
        end do
        ! Deallocate arrays
        deallocate(outRates, outRates1)
        if (allocated(outRates2)) deallocate(outRates2)
        if (allocated(outRates3)) deallocate(outRates3)
        if (allocated(outRates4)) deallocate(outRates4)

    end subroutine sp_solve

    !
    ! FUNCTIONS
    !

    ! **********************************
    function sp_isCreated(this) result(isValid)

        ! Test the validity of a sp object

        logical     :: isValid

        class(sp)  :: this

        isValid = allocated(this%sfh)
    end function sp_isCreated

    ! **********************************
    function sp_status(this) result(rates)

        ! Return the set of output rates
        !  e.i one per age and metalicity bin

        implicit none

        integer(kind=ikd)       :: i, j

        real(kind=rkd)          :: mbin

        type(gas), allocatable  :: rates(:, :) ! The set of output rates
        
        class(sp)               :: this
        !
        ! Create and set values of output rates
        allocate(rates(nAgeBins, nMetBins))
        do i = 1, nAgeBins
            do j = 1, nMetBins
                mbin = this%sfh(i, j)%mass
                rates(i ,j) = mbin * massLossRates(i, j)
            end do
        end do
    end function sp_status

    ! **********************************
    function g2s(g) result(B)

        ! Convert gas to stellar population

        implicit none

        external                       :: dgesv     ! Solver

        integer(kind=ikd)              :: i
        integer(kind=ikd)              :: rc
        integer(kind=ikd), allocatable :: pvt(:)    ! Pivot indices (list of swap operations).

        real(kind=rkd)                 :: mass
        real(kind=rkd)                 :: max, r
        real(kind=rkd), allocatable    :: B(:)      ! Mass in each bin
        real(kind=rkd), allocatable    :: A(:, :)   ! Conversion matrix

        type(gas), intent(in)          :: g         ! Gas to convert in sp

        ! Allocate
        allocate(A(nMetBins, nMetBins))
        allocate(B(nMetBins))
        allocate(pvt(nMetBins))

        ! Set masses
        do i = 1, nMetBins
            if (i == 1) then
                B(i) = g%mass
            else
                B(i) = g%elts(i-1)
            end if
        end do
        !
        ! Solve the system
        ! The matrix is used and modified in "dgesv", a copy is needed
        A = gas2sp
        call dgesv(nMetBins, 1, A, nMetBins, pvt, B, nMetBins, rc)
        !
        ! Clean up, very low values
        max = maxval(B)
        do i = 1, nMetBins
            r = abs(B(i)/max)
            if (r < 1.d-5) B(i) = 0.d0
        end do
        ! Normalisation
        mass = sum(B)
        B = g%mass/mass * B

        ! Deallocate
        deallocate(A)
        deallocate(pvt)

    end function g2s

    ! **********************************
    function ssp_evolve(this, dt, inRate, outRate) result(essp)

        ! Evolve a single stellar population by dt
        ! according to an input rate (SFR) "inRate"
        ! and an output rate (wind and sn ejecta) "outRate"

        implicit none

        real(kind=rkd), intent(in) :: dt

        type(gas), intent(in)      :: inRate  ! The input rate
        type(gas), intent(in)      :: outRate ! The output rate
        type(gas)                  :: newStar

        type(ssp)                 :: essp    ! The evolved ssp

        class(ssp)                 :: this    ! The current ssp

        ! Init the evolved ssp with the current one
        essp = this

        ! Substract ejected gas
        essp%mass = essp%mass - outRate%mass * dt
        ! Ejected mass do not impact the average age

        ! Add new formed stars
        newStar = inRate * dt
        if (newStar%mass > 0.d0) then
            ! New stars are formed continuously during dt
            ! The average age is therfore impacted
            essp%avgAge = (essp%avgAge * essp%mass + dt/2.d0 * newStar%mass) /&
                          (essp%mass + newStar%mass)
            ! Update mass
            essp%mass = essp%mass + newStar%mass
        end if

        ! Check mass in bin and reset if the stored mass is too low
        if ((essp%mass > 0.d0) .and. &
            (essp%mass < num_accuracy)) call this%create()

        ! Update formation time of this single stellar population
        if (essp%mass > 0.d0) essp%tform = essp%tform + dt

    end function ssp_evolve

    ! **********************************
    function sp_evolve(this, dt, inRate, outRates) result(esp)

        ! Evolve a stellar population by dt
        ! according to an input rate (SFR) "inRate"
        ! and a set of output rates (wind and sn ejecta) "outRates"

        implicit none
    
        integer(kind=ikd)                   :: i, j

        real(kind=rkd), intent(in)          :: dt             ! The time-step
        real(kind=rkd)                      :: mass           ! The updated total mass of the sp
        real(kind=rkd), dimension(nMetBins) :: inRateBins     ! Bin distribution of the gas turned into new stars

        type(gas), intent(in)               :: inRate         ! The input rate
        type(gas)                           :: sInRate
        type(gas), intent(in), allocatable  :: outRates(:, :) ! The set of output rates

        type(sp)                            :: esp            ! The evolved sp

        class(sp)                           :: this           ! The current sp

        ! Init gs structure from current
        esp = this
        !
        ! Distribute gas input stellar metalicity bins
        inRateBins = g2s(inRate * dt)
        !
        ! Loop over age and metalicity bins
        mass = 0. ! Init the total mass
        do i = 1, nAgeBins
            do j = 1, nMetBins
                ! Init the ssp input rate, null by default
                call sInRate%create()
                if (i == 1) then
                    ! First age bin,
                    ! Add input rate according to metalicity bin distribution
                    sInRate = inRateBins(j) * inRate
                end if
                this%sfh(i, j) = this%sfh(i, j)%evolve(dt, sInRate, outRates(i, j))
                call sInRate%delete()
                !
                ! Try mass transfer over age bins
                call this%transfer(i, j)
                ! 
                ! Update total mass
                mass = mass + this%sfh(i, j)%mass
            end do
        end do
        !
        ! Set new total gas mass
        this%mass = mass
        !
        ! sInRate is delete
        call sInRate%delete()

    end function sp_evolve

end module sp_mod