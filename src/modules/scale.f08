module scale_mod

    !*****************************************************************************************************************
    !
    !  OVERVIEW
    !
    !  The scale module defines all properties and procedures
    !  assossiated to a gas scale structuration.
    !  A set of scale obeject is then used to describe the gas structuration cascade
    !
    !*****************************************************************************************************************

    use parameters  ! Acces to global defintions and properties
    use log_mod     ! Acces to logging procedures
    use gas_mod     ! Acces to gas properties and procedures
    use status_mod  ! Acces to transfer rates status
    use model_mod   ! Acces to model parameters and integration scheme configuration

    implicit none
    
    public

    ! SCALE TYPE DEFINITION

    type scale
        real(kind=rkd)        :: l      ! The current scale
        type(gas)             :: gas    ! The gas stored at this scale > l-1 and < l
        type(scale), pointer  :: pScl   ! Pointer to the reference
    contains
        procedure   :: create => scale_create               ! Create a scale object
        procedure   :: delete => scale_delete               ! Delete a scale object
        procedure   :: copy => scale_copy                   ! Copy a scale object
        procedure   :: setAsRef => scale_setAsRef           ! Set pointer to reference
        procedure   :: add => scale_add                     ! Add gas to the scale
        procedure   :: sub => scale_sub                     ! Substract gas to the scale
        procedure   :: isValid => scale_isValid             ! Test the validity of the scale
        procedure   :: nClouds => scale_get_N_clouds        ! Get the number of active clouds at this scale
        procedure   :: BE_mass => scale_get_BE_mass         ! Get the Bonnor Ebert Mass associated to the scale
        procedure   :: volume => scale_volume               ! Get the volume of a sphere at the scale
        procedure   :: sV => scale_get_velocity_dispersion  ! Get the scaled velocity dispersion
        procedure   :: mu => scale_get_mass_surface_density ! Get the scaled mass surface density
        procedure   :: transfer => scale_transfer           ! Get the transfer rate to the lowest scale
        procedure   :: ststatus => scale_ststatus           ! Get the current status of the scale
        procedure   :: stevolve => scale_stevolve           ! Apply a integrator scheme step
        procedure   :: dtoptim => scale_dtoptim             ! Get the optimal time-step
        procedure   :: evolve => scale_evolve               ! Evolve the scale, for dt 
                                                            ! according to internal/external input/output
    end type scale
    !
    ! Define scale specific parameters
    real(kind=rkd), parameter :: scale_mass_accuracy = real(1.d-17, kind=rkd)

    ! STATUS
    type(status), target      :: sclStatus        ! The current status
    type(status), target      :: sclFinalStatus   ! The final status

    ! INTERFACE OPERATOR DECLARATIONS

    interface assignment (=)  ! allows to copy a gas component by using the symbol '='
        module procedure scale_copy
    end interface assignment (=)

contains

    !
    ! SUBROUTINES
    !
    ! **********************************
    subroutine scale_init()

        ! Initialize the scale module
        !  compute the slope index of the velocity disperion scalling relation

        implicit none

        ! Init dependencies
        call gas_init()

        ! Unit conversion
        ! Initialy lStar is given in pc
        lStar = lStar * pc2LenCU
        ! Initialy muStar is given in Msun/pc^2
        muStar = muStar * MSun2MassCU / pc2LenCU**2.
        ! Initialy sigmaStar is given in km/sec
        sigmaStar = sigmaStar * km2LenCU / s2TimeCU
        !
        ! The slope index of the velocity dispersion scalling relation is deduced from
        ! energy transfert conservation
        sigma_slope = real(1./3., kind=rkd)*(2.0 - mu_slope)
        !
        ! Compute the energy transfert rate per unit of volume
        ETRV = muStar*sigmaStar**3./lStar**2.

    end subroutine scale_init

    ! **********************************
    subroutine scale_create(this, l)

        ! Create a scale object

        implicit none

        real(kind=rkd)  :: l

        class(scale)  :: this

        ! Init or create fields
        this%l = l
        call this%gas%create()
        ! Set pointer
        this%pScl => null()

    end subroutine

    ! **********************************
    subroutine scale_delete(this)

        ! Delete a scale object

        implicit none

        class(scale)   :: this

        ! Reset or delete fields
        this%l = 0.d0
        call this%gas%delete()
        ! Unset pointer
        this%pScl => null()

    end subroutine

    ! **********************************
    subroutine scale_copy(this, aScl)

        ! Copy in "this" the scale object aScl
        ! Pointer to reference IS NOT COPIED
    
        implicit none

        type(scale), intent(in), target  :: aScl

        class(scale), intent(inout)      :: this
        
        ! Copy fields
        this%l = aScl%l
        this%gas = aScl%gas

    end subroutine scale_copy

    ! **********************************
    subroutine scale_setAsRef(this, aScl)

        ! Set pointer to the new reference
    
        implicit none

        type(scale), intent(in), target  :: aScl

        class(scale), intent(inout)      :: this

        ! Pointer to new refence
        this%pScl => aScl

    end subroutine scale_setAsRef

    ! **********************************
    subroutine scale_add(this, g)

        ! Add gas to the current scale

        implicit none

        type(gas), intent(in)            :: g

        class(scale)                     :: this

        call this%gas%add(g)

    end subroutine

    ! **********************************
    subroutine scale_sub(this, g)

        ! Substract gas to the current scale object

        implicit none

        type(gas), intent(in)            :: g

        class(scale)                     :: this

        call this%gas%sub(g)

    end subroutine

    ! **********************************
    subroutine scale_isValid(this, calledBy)

        ! Test a scale component
        !    After evolution, scale%l should be > 0
        !    The gas component should also be valid

        implicit none

        character(MAXPATHSIZE)              :: message
        character(MAXPATHSIZE), intent(in)  :: calledBy

        class(scale)                        :: this

        ! Test scale value
        if (this%l < 0.) then
            write(message, '(a)') 'Current scale is not valid.'
            call log_message(message, &
                             logLevel=LOG_ERROR, &
                             paramNames=(/'this%l                   '/), &
                             realParams=(/this%l/), &
                             calledBy=log_calledBy('scale_test', calledBy))
        end if
        !
        ! Test gas component of the scale
        call this%gas%isValid(log_calledBy('scale_test', calledBy))

    end subroutine scale_isValid

    ! **********************************
    subroutine scale_evolve(this, dt, inRate, outRate)

        ! Evolve, during dt, the current scale structure 
        ! according to :
        ! - Internal evolution processes
        ! - A constant input rate (due to external process) "inRate" and
        ! - A constant output rate (due to external process) "outRate"

        implicit none

        integer(kind=ikd)              :: st                ! Step index of evolution scheme

        logical                        :: hasBeenCorrected  ! False is status has been corrected
                                                            ! Only use for gsh status computation
        real(kind=rkd), intent(inout)  :: dt                ! The scale is evolve during dt

        type(gas), intent(in)          :: inRate            ! The (dt-)constant input rate
        type(gas), intent(in)          :: outRate           ! The (dt-)constant output rate

        type(scale), target            :: aScl              ! The current intermediate evolved scale
        type(scale), pointer           :: pScl              ! pointer to the intermediate state

        type(status), pointer          :: pFinalStatus      ! Pointer to the final status
        type(status), pointer          :: pStatus           ! Pointer to the current status

        class(scale), target           :: this              ! The current scale

        ! Init the first intermediate status to the current status
        aScl = this
        ! Set "this" a the reference
        call aScl%setAsRef(this)
        ! Pointer to current state
        pScl => aScl
        ! Pointer to status
        pFinalStatus => sclFinalStatus
        ! Loop over integration steps
        do st = 1, nSolverStep
            !
            ! Get the current status of the scale
            pStatus => sclStatus
            hasBeenCorrected = aScl%ststatus(inRate, outRate, pStatus)
            !
            ! Update final status
            call pFinalStatus%update(st, pStatus)
            !
            if (solver_isFinalStep(st)) then
                ! Switch to the final status
                pStatus => pFinalStatus
            end if
            !
            ! Adaptative time-step
            dt = this%dtoptim(st, dt, pStatus)
            !
            ! Get new state
            call this%stevolve(st, dt, pStatus, pScl)
        end do
        ! Save complete evolved state
        this = aScl
        !
        ! Delete tmp state
        call aScl%delete()

    end subroutine scale_evolve

    ! **********************************
    subroutine scale_stevolve(this, st, dt, pStatus, pScl)

        ! Apply on the intermediate scale structure (through its pointer: pScl)
        ! a new integrator step according to the current status (given through its pointer pStatus)

        implicit none
 
        integer(kind=ikd), intent(in) :: st                ! Integartion scheme index

        character(MAXPATHSIZE)        :: calledBy

        real(kind=rkd), intent(in)    :: dt                ! time-step
        real(kind=rkd)                :: wdt
        real(kind=rkd)                :: adt               ! The scale is evolve during dt
        real(kind=rkd)                :: dmOut, dmIn, dmTr ! Output, input, and transfered masses
        real(kind=rkd)                :: mass, finalMass

        type(status), pointer         :: pStatus  ! Pointer to the current status

        type(scale), pointer          :: pScl     ! Pointer to the intermediate state

        class(scale), target          :: this     ! The current state

        ! Compute new state
        wdt = solver_wdt(st)
        adt = wdt * dt
        !
        ! Copy of the current scale
        pScl = this
        !
        ! Save initial mass
        mass = this%gas%mass
        ! Compute mass variation
        dmIn = adt * pStatus%in%mass
        dmTr = adt * pStatus%tr%mass
        dmOut = adt * pStatus%out%mass
        !
        ! Test very low final mass (lead to numerical precision)
        finalMass = mass - dmTr - dmOut + dmIn
        if (solver_isFinalStep(st) .and. &                      ! In the last step od the integration scheme
            (dmTr > real(0.d0, kind=rkd) .or. dmOut > real(0.d0, kind=rkd)) .and. &   ! with mass transfer
            abs(finalMass) < 1.d-2 * scale_mass_accuracy) then  ! too low
            ! The final will be too low,
            ! Reset the scale, associated with its reference and return
            call pScl%create(this%l)
            call pScl%setAsRef(this)
            return
        end if
        !
        ! From this point, mass is sufficient
        ! Test mass evolution
        if (finalMass < real(0.d0, kind=rkd)) then
            call log_message('Current evolution step leads to negative mass', &
                            logLevel=LOG_ERROR, calledBy='scale_stevolve')
        end if
        !
        ! Evolution
        ! Add input mass
        call pScl%add(adt * pStatus%in)
        ! Substract external output mass
        call pScl%sub(adt * pStatus%out)
        ! Substract transfered mass
        call pScl%sub(adt * pStatus%tr)
        !
        ! Test the validity of the current scale
        write(calledBy, '(a)') 'scale_stevolve'
        call pScl%isValid(calledBy)

    end subroutine scale_stevolve

    !
    ! FUNCTIONS
    !

    ! **********************************
    function scale_get_N_clouds(this) result(N)

        ! Return the number of active clouds at this scale

        implicit none

        real(kind=rkd) :: N

        class(scale)   :: this

        N = real(floor(this%gas%mass/this%BE_mass(), kind=rkd), kind=rkd)
    
    end function scale_get_N_clouds

    ! **********************************
    function scale_get_BE_mass(this) result(M_BE)

        ! Return the Bonnot Ebert mass associated to the scale

        implicit none

        real(kind=rkd)  :: M_BE

        class(scale)    :: this

        M_BE = real(3./2., kind=rkd) * sigmaStar**4./GCst_CU**2./muStar*(this%l/lStar)**(4.*sigma_slope - mu_slope)
    
    end function scale_get_BE_mass

    ! **********************************
    function scale_volume(this) result(V)

        ! Return the volume of a sphere of scale l

        real(kind=rkd)   :: V

        class(scale)     :: this

        V = real(4./3., kind=rkd) * pi * (this%l/real(2., kind=rkd))**3.

    end function scale_volume

    ! **********************************
    function scale_get_velocity_dispersion(this) result(sV)

        ! Return the mean 1D-velocity dispersion at this scale

        implicit none

        real(kind=rkd)    :: sV

        class(scale)      :: this

        sV = sigmaStar * (this%l / lStar)**sigma_slope

    end function scale_get_velocity_dispersion

    ! **********************************
    function scale_get_mass_surface_density(this) result(mu)

        ! Return the mean mass surface density at this scale

        real(kind=rkd)   :: mu

        class(scale)     :: this

        mu = muStar * (this%l / lStar)**mu_slope

    end function scale_get_mass_surface_density

    ! **********************************
    function scale_transfer(this) result(rate)

        ! Return the current transfer rate of the scale

        real(kind=rkd)         :: vRate

        type(gas)              :: rate  ! output rate

        class(scale)           :: this

        ! Transfer mass to the lower scale is only possible
        ! if the current scale is unstable (this%nClouds() > 0)
        ! If the ref value is null, no transfert
        vRate = real(3.d0/2.d0, kind=rkd) * ETRV * this%volume() / this%sV()**2. * this%nClouds()
        !
        ! The transfer rate gas object is build according to the gas signature
        rate = vRate * this%gas%signature()

    end function scale_transfer

    ! **********************************
    function scale_ststatus(this, inRate, outRate, pStatus) result(isOK)

        ! Compute current status of the scale structure
        ! according to :
        ! - Internal evolution processes "trRate"
        ! - A constant input rate (due to external process) "inRate"
        ! - A constant output rate (due to external process) "outRate"

        implicit none

        logical                       :: isOK

        type(gas), intent(in)         :: inRate       ! The (dt-)constant input rate
        type(gas), intent(in)         :: outRate      ! The (dt-)constant output rate
        type(gas)                     :: trRate       ! Internal transfer rate of the scale
        type(gas)                     :: corrOutRate  ! Corrected output rate (if needed)

        type(status), pointer         :: pStatus      ! Pointer to the "current" status

        class(scale)                  :: this         ! The current scale

        ! Init
        isOK = .True.

        ! Internal transfer rate
        trRate = this%transfer()
        !
        ! Init corrected outRate
        corrOutRate = outRate
        !
        ! Test consistency with ref value
        if (.not. associated(this%pScl)) then
            call log_message('pScl pointer not associated ! ', logLevel=LOG_ERROR)
            return
        end if
        if (this%pScl%gas%mass == 0.d0) then
            ! Reference mass is null
            if (outRate%mass + trRate%mass > inRate%mass) then
                ! For this state,
                ! external and internal processes lead to mass depletion
                ! Set trRate = 0 and outRate = inRate to ensure no mass variation and keep mass = 0.
                call trRate%create()  ! no transfer
                corrOutRate = inRate  ! out = in
                isOK = .False.
            end if
        end if
        !
        ! Set status
        call pStatus%set(inRate, trRate, corrOutRate)
        !
        ! Delete tmp objects
        call corrOutRate%delete()

    end function scale_ststatus

    ! **********************************
    function scale_dtoptim(this, st, dt, pStatus) result(adt)

        ! Return the optimal time-step accoring to current status
        ! and mass

        implicit none

        integer(kind=ikd), intent(in) :: st

        real(kind=rkd), intent(in)    :: dt       ! time-step
        real(kind=rkd)                :: adt

        type(status), pointer         :: pStatus  ! Pointer to the current status

        class(scale)                  :: this     ! The current state

        adt = pStatus%dtMax(st, dt, this%gas%mass)

    end function scale_dtoptim

end module scale_mod

