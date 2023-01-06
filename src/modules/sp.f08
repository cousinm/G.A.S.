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
    use ssp_mod     ! Acess to single stellar populations methods ans properties
    use config_mod  ! Acces to configurations parameters (path)

    implicit none

    public

    ! A complete stellar population is built from a set of single stellar population
    !   A ssp evolved in a set of age and metalicity bin
    !
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
        procedure  :: stevolve => sp_stevolve     ! 
        procedure  :: evolve => sp_evolve         ! Evolve a stellar population structure by dt
    end type sp

    ! Define sp specific parameters
    real(kind=rkd), allocatable :: gas2sp(:, :)   ! gas to sp matrix conversion, based on abundancies

    ! TEMPORARY INTERMEDIATE STATUS
    ! Targets
    type(status), allocatable, target :: spStatus(:, :)

    ! INTERFACE OPERATOR DECLARATIONS

    interface assignment (=)  ! allows to copy a sp component by using the symbol '='
        module procedure sp_copy
    end interface assignment (=)

contains

    !
    ! SUBROUTINES
    !

    ! **********************************
    subroutine sp_init()

        ! Initialize the sp module and dependancies
        !  Read stellar population properties (mass loss rates, SN rates, spectrum)

        implicit none

        integer(kind=ikd)     :: i, j

        ! Init ssp
        call ssp_init()

        ! Init intermediate status
        allocate(spStatus(nAgeBins, nMetBins))

        ! Build gas to stellar population conversion matrix
        ! Allocate the matrix
        allocate(gas2sp(nMetBins, nMetBins))
        ! Init it
        do i = 1, nMetBins
            do j = 1, nMetBins
                if (i == 1) then
                    gas2sp(i, j) = 1.
                else
                    gas2sp(i, j) = initAbund(j)%elts(i - 1)
                end if
            end do
        end do

    end subroutine sp_init

    ! **********************************
    subroutine sp_finalize()

        ! Delete and deallocate all sp objects

        implicit none

        integer(kind=ikd)     :: i, j

        ! Deallocate spStatus
        if (allocated(spStatus)) then
            do i = 1, nAgeBins
                do j = 1, nMetBins
                    call spStatus(i, j)%delete()
                end do
            end do
            deallocate(spStatus)
        end if

        ! Build gas to stellar population conversion matrix
        ! Allocate the matrix
        if (allocated(gas2sp)) deallocate(gas2sp)

    end subroutine sp_finalize

    ! **********************************
    subroutine sp_create(this)

        ! Create a complete stellar population

        implicit none

        integer(kind=ikd)    :: iAge ! Age loop index
        integer(kind=ikd)    :: jMet ! Metalicity loop index

        class(sp)            :: this

        this%mass = 0.d0   ! Total mass of the stellar population
        this%mAge = 0.d0   ! Mass-weighted average age
        this%lAge = 0.d0   ! Luminosity-weighted average age
        this%mZ = 0.d0     ! Mass-weighted average metelicity
        this%lZ = 0.d0     ! Luminosity-weighted average metalicity

        allocate(this%sfh(nAgeBins, nMetBins)) ! Star Formation History
        do iAge = 1, nAgeBins
            do jMet = 1, nMetBins
                call this%sfh(iAge, jMet)%create(iAge, jMet)
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

    subroutine sp_stevolve(this, dt, st, inRate, outRate, pStatus, pSp)

        ! Compute current status of the stellar population
        ! according to :
        ! - Internal evolution processes
        ! - A constant input rate (due to external process) "inRate" and
        !
        ! No ejection from any external processes can be applied to sp
        ! Output rate is an output parameter storing wind/sn ejection from the current sp
        !
        ! Then
        ! Apply the next step "st" of the complete integration scheme

        implicit none
    
        integer(kind=ikd), intent(in)   :: st
        integer(kind=ikd)               :: iAge, jMet

        real(kind=rkd), intent(in)      :: dt             ! The time-step
        real(kind=rkd)                  :: mass           ! Total stellar mass
        real(kind=rkd)                  :: mAge           ! mass weighted age
        real(kind=rkd)                  :: mZ             ! average metalicity

        type(gas), intent(in)           :: inRate         ! The external input rate (SFR)
        type(gas), intent(out)          :: outRate        ! Global output rate
        type(gas)                       :: spOutRate      ! global wind/sn gas ejection rate
        type(gas)                       :: sspInRate      ! Input rate of the current ssp
        type(gas)                       :: sspOutRate     ! Wind/sn gas ejection rate of the current ssp
        type(gas), allocatable          :: inRateBins(:)  ! Gas distribution into metalicity bins

        type(ssp), pointer              :: pSsp           ! Pointer to the current ssp stage
        type(sp), pointer               :: pSp            ! pointer to the current stage
        
        type(status), pointer           :: pStatus(:, :)  ! Pointer to the final status of the current sp
        type(status), pointer           :: pSspStatus     ! Pointer to the current status of the ssp

        class(sp)                       :: this           ! The current stellar population

        ! Init sOutRate, sInRate and outRate
        call sspInRate%create()
        call sspOutRate%create()
        call spOutRate%create()
        if (.not. outRate%isCreated()) then
            call outRate%create()
        end if

        ! Init total mass
        mass = real(0.d0, kind=rkd)
        do jMet = 1, nMetBins
            !
            do iAge = 1, nAgeBins
                !
                ! Input rate (SFR) is only applied to the youngest ssp
                !
                if (iAge == 1) then
                    ! Distribute gas input stellar metalicity bins
                    inRateBins = g2s(inRate)
                    sspInRate = sspInRate + inRateBins(jMet)
                end if
                !
                
                ! Get pointers
                if (st == 1) then
                    ! Reset "final" status to null
                    call pStatus(iAge, jMet)%reset()
                end if
                pSspStatus => pStatus(iAge, jMet)
                pSsp => pSp%sfh(iAge, jMet)
                !
                if (pSsp%mass > 0.d0 .or. sspInRate%mass > 0.d0) then
                    !
                    ! Apply the new solver step
                    ! In output of stevolve, sspInRate contains the transfer rate
                    ! from the current age bin to the next one
                    call this%sfh(iAge, jMet)%stevolve(dt, st, sspInRate, sspOutRate, pSspStatus, pSsp)
                    !
                    ! Update global wind/sn ejection rate
                    spOutRate = spOutRate + pSspStatus%out
                    !
                    ! Update total stellar mass
                    mass = mass + pSp%sfh(iAge, jMet)%mass
                    !
                    ! Update mass weighted age
                    mAge = mAge + pSp%sfh(iAge, jMet)%mass * pSp%sfh(iAge, jMet)%avgAge
                    !
                    ! Update average metalicity
                    mZ = mZ + pSp%sfh(iAge, jMet)%mass * metBins(jMet)
                end if
            end do ! Met loop
        end do ! Age loop
        !
        ! Set total mass
        pSp%mass = mass
        !
        ! Set mass weighted age
        pSp%mAge = mAge / mass
        !
        ! Set average metalicity
        pSp%mZ = mZ / mass
        !
        ! Set current output rate and delete local output rate object
        outRate = spOutRate
        call spOutRate%delete()
        !
        ! Delete all other gas object created here
        call sspInRate%delete()
        call sspOutRate%delete()

    end subroutine sp_stevolve

    ! **********************************
    subroutine sp_evolve(this, dt, inRate, outRate)

        ! Evolve, during dt, the current complete stellar population
        ! according to :
        ! - Internal evolution processes
        ! - A constant input rate (due to external process) "inRate"
        !
        ! No ejection from any external processes can be applied to sp
        ! Output rate is an output parameter storing wind/sn ejection from the current sp

        implicit none

        integer(kind=ikd)           :: st             ! Step index of evolution scheme

        real(kind=rkd), intent(in)  :: dt             ! The scale is evolve during dt

        type(gas), intent(in)       :: inRate         ! The (dt-)constant input rate
        type(gas), intent(out)      :: outRate        ! Pointer to the wind/sn output rate

        type(sp), target            :: aSp            ! A sp matching the current state
        type(sp), pointer           :: pSp            ! Pointer to the current state
        
        type(status), pointer       :: pStatus(:, :)  ! Pointer to the final status of the current sp

        class(sp)                   :: this           ! The current stellar population

        ! Init intermediate status with the current status
        aSp = this
        pSp => aSp
        ! Define complete corrected status
        pStatus => spStatus
        do st = 1, nSolverStep
            call this%stevolve(dt, st, inRate, outRate, pStatus, pSp)
        end do
        ! Save complete evolved state
        this = aSp

    end subroutine sp_evolve

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
    function g2s(g) result(gasBins)

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
        type(gas), allocatable         :: gasBins(:)

        ! Allocate B and gasBins
        allocate(B(nMetBins))
        allocate(gasBins(nMetBins))
        !
        ! Set masses
        do i = 1, nMetBins
            if (i == 1) then
                B(i) = g%mass
            else
                B(i) = g%elts(i-1)
            end if
        end do
        !
        ! Test null values
        if (sum(B) <= num_accuracy) then
            do i = 1, nMetBins
                gasBins(i) = real(0.d0, kind=rkd) * initAbund(i)
            end do
            deallocate(B)
            return
        end if
        !
        ! Allocate other structures
        allocate(A(nMetBins, nMetBins))
        allocate(pvt(nMetBins))
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
        !
        ! Convert in gas object
        do i = 1, nMetBins
            gasBins(i) = B(i) * initAbund(i)
        end do

        ! Deallocate
        deallocate(A)
        deallocate(B)
        deallocate(pvt)

    end function g2s

end module sp_mod