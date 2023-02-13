module gsh_mod

    !*****************************************************************************************************************
    !
    !  OVERVIEW
    !
    ! The gsh module defines all properties and procedures
    ! associated to gas structuration history structure
    ! This structure allows to "follow" gas structuration process
    ! to the different binned scales
    !
    !*****************************************************************************************************************

    use parameters  ! Acces to global defintions and properties
    use model_mod   ! Acces to model parameters and integration scheme configuration
    use status_mod  ! Acces to transfer rates status
    use scale_mod   ! Acces to scale properties and procedures

    implicit none

    public

    type gsh
        integer(kind=ikd)          :: ih             ! Current largest scale (at least one Cloud)
        real(kind=rkd)             :: mass           ! Total mass of the structure
        type(status)               :: myStatus       ! Global status of the gsh, intput/transfer/output
        type(scale), allocatable   :: cascade(:)     ! Gas structuration cascade, set of scale
    contains
        procedure  :: create => gsh_create           ! Create the gas structuration (cascade)
        procedure  :: delete => gsh_delete           ! Delete the gas structuration (cascade)
        procedure  :: copy => gsh_copy               ! Copy a gsh object
        procedure  :: ststatus => gsh_ststatus       !
        procedure  :: stevolve => gsh_stevolve       ! Evolve, for dt, the next step of the integration scheme
        procedure  :: dtoptim => gsh_dtoptim         ! Get the optimal time-step
        procedure  :: evolve => gsh_evolve           !         according to internal/external input/output
        procedure  :: Vesc => gsh_escape_velocity    ! Return the escape velocity of the gsh
        procedure  :: iSFR => gsh_instantaneous_SFR  ! Return the instantaneous SFR
                                                     ! i.e the SFR computed according to the current state
                                                     ! the value is different of the RK solver value
                                                     ! that use a combination of instantaneous SFR
        procedure  :: disrupt => gsh_disrupt         ! Compute instantaneous gas disrupt rate
                                                     ! according to the instantaneous injected power
    end type gsh

    ! STATUS
    type(status), allocatable, target     :: gshStatus(:)
    type(status), allocatable, target     :: gshFinalStatus(:)

    ! INTERFACE OPERATOR DECLARATIONS

    interface assignment (=)  ! allows to copy a gsh component by using the symbol '='
        module procedure gsh_copy
    end interface assignment (=)

contains

    !
    ! SUBROUTINES
    !

    ! **********************************
    subroutine gsh_init()

        ! Initialize the gsh module and dependancies

        implicit none

        ! Init scale
        call scale_init()

        ! Init status
        if (.not. allocated(gshStatus)) then
            allocate(gshStatus(nScales))
        end if
        if (.not. allocated(gshFinalStatus)) then
            allocate(gshFinalStatus(nScales))
        end if

    end subroutine gsh_init

    ! **********************************
    subroutine gsh_create(this)

        ! Create (allocate) a gas history structure

        implicit none

        integer(kind=ikd)  :: i

        real(kind=rkd)     :: l

        class(gsh)       :: this

        this%ih = 1                     ! Init the one-cloud scale index
        this%mass = 0.d0                ! Total mass of the structure
        call this%myStatus%create()     ! Create and init status
        !
        ! Create the cascade with nScales cells
        allocate(this%cascade(nScales)) 
        ! Create gas element in each cells
        do i = 1, nScales
            l = gsh_i2l(i)
            call this%cascade(i)%create(l)
        end do

    end subroutine gsh_create

    ! **********************************
    subroutine gsh_delete(this)

        ! Delete a gas structuration history structure

        implicit none

        integer(kind=ikd)  :: il

        class(gsh)         :: this

        this%ih   = 1               ! cell index of the largest scale (at least one cloud)
        this%mass = 0.d0            ! Total mass of the structure
        call this%myStatus%delete() !

        ! Delete gas element in each cells
        do il = 1, nScales
            call this%cascade(il)%delete()
        end do

        deallocate(this%cascade) ! Delete the cascade

    end subroutine gsh_delete

    ! **********************************
    subroutine gsh_copy(this, aGsh)

        ! Copy in "this" the gsh object aGsh

        implicit none

        integer(kind=ikd)           :: s

        type(gsh), intent(in)       :: aGsh

        class(gsh), intent(inout)   :: this

        ! Copy
        this%ih       = aGsh%ih
        this%mass     = aGsh%mass
        this%myStatus = aGsh%myStatus
        ! 
        ! Allocate cascade
        if (.not. allocated(this%cascade)) then
            ! Create the cascade with nScales cells
            allocate(this%cascade(nScales))
        end if
        ! Copy each scale
        do s = 1, nScales
            this%cascade(s) = aGsh%cascade(s)
            call this%cascade(s)%setAsRef(aGsh%cascade(s))
        end do
    
    end subroutine gsh_copy

    ! **********************************
    subroutine gsh_evolve(this, dt, l, Q, inRate)

        ! Evolve, during dt, the current gas structuration
        ! according to :
        ! - Internal evolution processes
        ! - A gas input rate (due to external fresh accretion): "inRate" and
        ! - A power input (due to disruptive external processes): Q
        ! - fresh gas accreation is injected at scale "l"
        !
        ! In output
        ! - "outRate" is the global output rate produced by disruption processes on the cascade

        implicit none

        integer(kind=ikd)                :: st               ! Step index of evolution scheme

        real(kind=rkd), intent(inout)    :: dt               ! The scale is evolve during dt
        real(kind=rkd), intent(in)       :: l                ! The injection scale
        real(kind=rkd), intent(in)       :: Q                ! Instantaneous gas injected power

        type(gas), intent(in)            :: inRate           ! The (dt-)constant input rate

        type(gsh), target                :: aGsh             ! The current state
        type(gsh), pointer               :: pGsh             ! Pointer to the current state

        type(status), pointer            :: pFinalStatus(:)  ! Pointer to the final status
        type(status), pointer            :: pStatus(:)       ! Pointer to the current status

        class(gsh)                       :: this             ! The current gas structuration

        ! Init the first intermediate status to the current status
        aGsh = this
        ! Pointer to the current state
        pGsh => aGsh
        ! Pointer to status
        pFinalStatus => gshFinalStatus
        ! Loop over integration steps
        do st = 1, nSolverStep
            !
            ! Get the current status of the scale
            pStatus => gshStatus
            call pGsh%ststatus(st, l, Q, inRate, pStatus, pFinalStatus)
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
            call this%stevolve(st, dt, pStatus, pGsh)
        end do
        ! Save complete evolved state
        this = aGsh
        ! Delete the tmp gsh structure
        call aGsh%delete()

    end subroutine gsh_evolve

    ! **********************************
    subroutine gsh_ststatus(this, st, l, Q, inRate, pStatus, pFinalStatus)

        ! Compute current status of the gas structuration history
        ! according to :
        ! - Internal evolution processes
        ! - A gas input rate (due to external fresh accretion): "inRate" and
        ! - A power input (due to disruptive external processes): Q
        ! - fresh gas accreation is injected at scale "l"

        implicit none

        integer(kind=ikd), intent(in)   :: st
        integer(kind=ikd)               :: il
        integer(kind=ikd)               :: s

        logical                         :: sclIsOK
        logical                         :: isOK

        real(kind=rkd), intent(in)      :: l                 ! The injection scale
        real(kind=rkd), intent(in)      :: Q                 ! The gas injected powerde

        type(gas), intent(in)           :: inRate            ! The external input rate
        type(gas)                       :: outRate           ! Large scale ejection rate

        type(gas)                       :: sclInRate         ! Input rate
        type(gas)                       :: sclOutRate        ! Output rate
        type(gas), allocatable          :: disruptRates(:)   ! Disruption rates
                                                             ! Gas transfered from lower to larger scale
        type(scale), pointer            :: pScl              ! Pointer to the current scale

        type(status), pointer           :: pFinalStatus(:)   ! Pointer to the "final" status
        type(status), pointer           :: pStatus(:)        ! Pointer to the "current" status
        type(status), pointer           :: pSclStatus        ! Pointer to the current state of a scale

        class(gsh), target              :: this              ! The current gsh
        !
        ! Compute disruption rates according to power injection
        disruptRates = this%disrupt(Q)
        !
        ! Get the scale index associated to fresh gas injection scale
        il = gsh_l2i(l)
        !
        ! Init sclInRate
        call sclInRate%create()
        ! Init large scale output
        call outRate%create()
        ! Reset global status
        call this%myStatus%create()
        !
        ! Init corrected status
        sclIsOK = .True.
        isOK = .True.
        !
        ! Get the current status of the gas structuration
        ! Run through the different scales,
        ! from the largest one to the lowest one
        ! FIRST loop (predictor)
        do s = nScales, 1, -1
            !
            ! Get pointers
            pSclStatus => pStatus(s)
            pScl => this%cascade(s)
            !
            ! Compute current global status due to
            ! external input/output
            ! The largest scale can only be fed by external input
            !
            ! At injection scale, add fresh input rate
            if (s == il) then
                sclInRate = sclInRate + inRate
            end if
            ! Take into account mass transfer due to disruption process
            call sclOutRate%create()  ! Reset local output rate
            ! In all cases, disrupt rate should be takes into account
            if (disruptRates(s)%isCreated()) sclOutRate = disruptRates(s)
            ! But a disrupt rate is an injection rate at the upper scale only below the (one-cloud scale)
            if (s > 1 .and. s <= this%ih) then
                ! Scale 's' receives disrution rate from scale s-1
                if (disruptRates(s - 1)%isCreated()) sclInRate = sclInRate + disruptRates(s - 1)
            end if
            !
            ! Update scale status, take into account internal transfer rate
            ! In output indicates if a correction has been done on the status
            sclIsOK = pScl%ststatus(sclInRate, sclOutRate, pSclStatus)
            if (.not. sclIsOK .and. isOK) then
                ! A correction has been done on this scl,
                ! output and transfer status has been impacted
                isOK = .False.
            end if
            !
            ! The transfer rate of the scale "s" become an input rate the next scale
            sclInRate = pSclStatus%tr
        end do
        !
        ! Reset inRate
        call sclInRate%create()
        !
        ! SECOND loop (corrector if need and set final status)
        do s = 1, nScales
            !
            ! Get pointers
            pSclStatus => pStatus(s)
            pScl => this%cascade(s)
            !
            ! Corrector ?
            if (.not. isOK) then
                ! Apply corrections
                !
                if (s < nScales) then
                    sclInRate = sclInRate + pStatus(s + 1)%tr
                end if
                ! At injection scale, add fresh input rate
                if (s == il) then
                    sclInRate = sclInRate + inRate
                end if
                ! Update scale status
                sclIsOK = pScl%ststatus(sclInRate, pSclStatus%out, pSclStatus)
                !
                ! The output rate of the scale "s" become an input rate the next scale
                if (s < this%ih) then
                    sclInRate = pSclStatus%out
                else
                    call sclInRate%create()
                end if
            end if
            !
            ! Update final status
            call pFinalStatus(s)%update(st, pSclStatus)
            !
            ! Update large scale ejecta rate
            if (s >= this%ih .and. pFinalStatus(s)%out%isCreated()) then
                outRate = outRate + pFinalStatus(s)%out
            end if
        end do
        !
        ! Set the average input rate applied
        ! Set the average transfert rate  at the lowest scale
        ! Set the average large scale output applied
        call this%myStatus%set(pFinalStatus(il)%in, pFinalStatus(1)%tr, outRate)
        !
        ! Deallocate local structure
        deallocate(disruptRates)

    end subroutine gsh_ststatus

    ! **********************************
    subroutine gsh_stevolve(this, st, dt, pStatus, pGsh)

        ! Apply the next step "st" of the complete integration scheme
        ! according to aStatus
        !
        ! In output
        ! - "outRate" is the global output rate produced by disruption processes on the cascade

        implicit none
    
        integer(kind=ikd), intent(in)   :: st
        integer(kind=ikd)               :: ih
        integer(kind=ikd)               :: s

        real(kind=rkd), intent(in)      :: dt              ! The time-step
        real(kind=rkd)                  :: mass            ! mass contains in the complete cascade

        type(scale), pointer            :: pScl            ! Pointer to the current scale

        type(status), pointer           :: pStatus(:)      ! Pointer to the current status
        type(status), pointer           :: pSclStatus

        type(gsh), pointer              :: pGsh            ! Pointer to the intermediate state

        class(gsh)                      :: this            ! The current gsh

        ! Init
        ! Index of the one-cloud scale (or largest scale with mass)
        ! This index is assumed to be constant for a complete solver step evaluation (st = 1, 2, 3, 4)
        ! and it is therefore up to date only for the final step of the integration scheme
        ih = -1
        ! Total mass
        mass = real(0.d0, kind=rkd)
        ! Get current status
        this%myStatus = pGsh%myStatus
        ! And reset the new state to the current one (status is kept)
        pGsh = this
        !
        ! Compute new state
        do s = nScales, 1, -1
            !
            ! Get pointers
            pSclStatus => pStatus(s)
            pScl => pGsh%cascade(s)
            !
            if (pScl%gas%mass > 0.d0 .or. pSclStatus%in%mass > 0.d0) then
                !
                ! Apply the new solver step for the scale
                call this%cascade(s)%stevolve(st, dt, pSclStatus, pScl)
                !
                ! Update total cascade mass
                mass = mass + pScl%gas%mass
                !
                ! Update the one-cloud scale index
                if (ih == -1) then
                    ih = s  ! Init to the first scale containing mass
                end if
                if (s < nScales) then
                    if (pScl%nClouds() > 0 .and. this%cascade(s + 1)%nClouds() == 0) then
                        ! Set to the first scale with at least one Cloud
                        ih = s
                    end if
                end if
            end if
        end do
        !
        ! Set new ih
        ! This index is assumed to be constant for a complete solver step evaluation (st = 1, 2, 3, 4)
        ! and it is therefore up to date only for the final step of the integration scheme
        if (solver_isFinalStep(st)) pGsh%ih = ih
        !
        ! Set total mass
        pGsh%mass = mass

    end subroutine gsh_stevolve

    !
    ! FUNCTIONS
    !

    ! **********************************
    function gsh_i2l(i) result(l)

        ! Return the ith scale of the cascade
        ! in the cell i, gas is structured at scale (l) equal or lower than l(i)
        ! e.g cascade(i=1) stores gas structured at scale lower or equal than lstar = 0.1pc

        implicit none

        integer(kind=ikd), intent(in)  :: i

        real(kind=rkd)                 :: l

        l = stepFactor**(i-1)*lStar
        return

    end function gsh_i2l

    ! **********************************
    function gsh_l2i(l) result(i)

        ! Return the index (in the cascade array) associated to the scale l

        implicit none

        integer(kind=ikd)  :: i, s

        real(kind=rkd)     :: l    ! The scale

        do s = 1, nScales
            if (gsh_i2l(s) >= l) then
                exit
            end if
        end do
        i = min(s, nScales)

    end function gsh_l2i

        ! **********************************
    function gsh_escape_velocity(this) result(Vesc)

        ! Return the escape velovity of the current gas structuration
        ! The total mass (gsh%mass) of the gas is distributed into Nclouds formed
        ! at the largest occupied scale (gsh%l)

        implicit none

        real(kind=rkd)       :: l
        real(kind=rkd)       :: Vesc
        real(kind=rkd)       :: NClouds
        real(kind=rkd)       :: mass

        class(gsh)         :: this

        ! The mass is distributed in the NClouds formed 
        ! at the largest occupied scale
        NClouds = this%cascade(this%ih)%NClouds()
        l = this%cascade(this%ih)%l
        Vesc = 0.d0 ! Init
        if (NClouds > 0.d0) then
            mass = this%mass / NClouds
            !
            Vesc = sqrt(2.d0*GCst_CU*mass/l)  ! [CU]
        end if

    end function gsh_escape_velocity

    ! **********************************
    function gsh_dtoptim(this, st, dt, pStatus) result(adt)

        ! Return the optim time-step accoring to current status
        ! and mass

        implicit none

        integer(kind=ikd), intent(in) :: st
        integer(kind=ikd)             :: s

        real(kind=rkd), intent(in)    :: dt         ! time-step
        real(kind=rkd)                :: adt
        real(kind=rkd)                :: mass

        type(status), pointer         :: pStatus(:) ! Pointer to the current set of status
        type(status), pointer         :: pSclStatus ! Pointer to a current scale status

        type(scale), pointer          :: pScl       ! Pointer to a scale

        class(gsh), target            :: this

        ! Run through the different scales
        ! and deduce optimal time step
        ! from the largest one to the lowest one
        ! Init adaptative time-step, to current time-step
        adt = dt
        do s = nScales, 1, -1 
            !
            ! Get pointers
            pSclStatus => pStatus(s)
            pScl => this%cascade(s)
            mass = pScl%gas%mass
            if (mass > 0.d0) then
                adt = pSclStatus%dtMax(st, adt, pScl%gas%mass)
            end if
        end do

    end function gsh_dtoptim

    ! **********************************
    function gsh_instantaneous_SFR(this) result (sfr)

        ! Return the instantaneous SFR
        ! i.e the SFR computed according to the current state
        ! the value is different of the RK solver value
        ! that use a combination of instantaneous SFR

        implicit none

        type(gas)          :: sfr

        class(gsh)         :: this

        ! The instantaneous star formation rate is given by
        ! the "transfert" rate of the lowest scale of the gas structuration history
        sfr = this%cascade(1)%transfer()

    end function gsh_instantaneous_SFR

    ! **********************************
    function gsh_disrupt(this, Q) result(disruptRate)

        implicit none

        integer(kind=ikd)           :: s

        real(kind=rkd), intent(in)  :: Q              ! Instantaneous injected power
        real(kind=rkd)              :: Qm
        real(kind=rkd)              :: Qs
        real(kind=rkd)              :: mbin
        real(kind=rkd)              :: sV             ! Velocity dispersion of the s + 1 scale
        real(kind=rkd)              :: vRate

        type(gas), allocatable      :: disruptRate(:) ! Set of gas disrupt rates

        class(gsh)                  :: this           ! The current gas structuration

        ! Create and init outputs
        ! this array SHOULD be deallocate in stevolve procedure
        allocate(disruptRate(nScales))

        if (this%mass > 0.d0 .and. Q > 0.d0) then
            ! Power is injected at all scales
            Qm = Q / this%mass
            do s = 1, nScales
                ! Init
                call disruptRate(s)%create()
                mbin = this%cascade(s)%gas%mass
                if (mbin > 0.d0) then
                    Qs = Qm * mbin
                    if (s < this%ih) then  ! this%ih <= nScale
                        ! Lowest scales, lower than the 1-cloud scale: disruption with transfer to the higher scale
                        ! this%ih <= nScale
                        sV = abs(this%cascade(s + 1)%sV() - this%cascade(s)%sV())
                    else
                       sV = (Vwind - this%cascade(s)%sV())
                    end if
                    !
                    vRate = gasDisruptionEfficiency * real(2.d0, kind=rkd) * Qs / sV**2.
                    ! Gas conversion is performed with the signature
                    disruptRate(s) = vRate * this%cascade(s)%gas%signature()
                end if
            end do
        end if

    end function gsh_disrupt

end module gsh_mod