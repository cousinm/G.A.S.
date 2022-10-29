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
    use model_mod   ! Acces to model parameters and integration scheme configuration

    implicit none
    
    public

    ! SCALE TYPE DEFINITION

    type scale
        real(kind=rkd)  :: l        ! The current scale
        type(gas)     :: gas      ! The gas stored at this scale > l-1 and < l
    contains
        procedure   :: create => scale_create               ! Create a scale object
        procedure   :: delete => scale_delete               ! Delete a scale object
        procedure   :: copy => scale_copy                   ! Copy a scale object
        procedure   :: add => scale_add                     ! Add gas to the scale
        procedure   :: sub => scale_sub                     ! Substract gas to the scale
        procedure   :: isValid => scale_isValid             ! Test the validity of the scale
        procedure   :: nClouds => scale_get_N_clouds        ! Get the number of active clouds at this scale
        procedure   :: BE_mass => scale_get_BE_mass         ! Get the Bonnor Ebert Mass associated to the scale
        procedure   :: volume => scale_volume               ! Get the volume of a sphere at the scale
        procedure   :: sV => scale_get_velocity_dispersion  ! Get the scaled velocity dispersion
        procedure   :: mu => scale_get_mass_surface_density ! Get the scaled mass surface density
        procedure   :: status => scale_status               ! Get the output rate of the scale
        procedure   :: update => scale_update               ! Update the scale during dt
        procedure   :: solve => scale_solve                 ! Solve, for dt, the next step of the integration scheme
        procedure   :: stevolve => scale_stevolve           ! Evolve, for dt, the next step of the integration scheme
        procedure   :: evolve => scale_evolve               !         according to internal/external input/output
    end type scale

    ! TEMPORARY INTERMEDIATE STATUS
    ! Targets
    type(gas), target     :: sclStatus

    ! INTERFACE OPERATOR DECLARATIONS

    interface assignment (=)  ! allows to copy a gas component by using the symbol '='
        module procedure scale_copy_
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

        ! Init temporary intermediate status
        call sclStatus%create()

        ! Unit conversion
        ! Initialy lStar is given in pc
        lStar = lStar*pc2kpc
        ! Initialy muStar is given in Msun/pc^2
        muStar = muStar*MSun2mass/pc2kpc**2.
        ! Initialy sigmaStar is given in km/sec
        sigmaStar = sigmaStar*km2kpc/s2Gyr
        !
        ! The slope index of the velocity dispersion scalling relation is deduced from
        ! energy transfert conservation
        sigma_slope = 1./3.*(2.0 - mu_slope)
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

        this%l = l
        call this%gas%create()

    end subroutine

    ! **********************************
    subroutine scale_delete(this)

        ! Delete a scale object

        implicit none

        class(scale)   :: this

        this%l = 0.d0
        call this%gas%delete()

    end subroutine

    ! **********************************
    subroutine scale_copy(s1, s2)

        ! Copy the gas object g2 into the gas object g1
    
        implicit none

        type(scale), intent(in)          :: s2

        class(scale), intent(inout)      :: s1
        
        ! Copy fields
        s1%l = s2%l
        s1%gas = s2%gas

        return
    end subroutine scale_copy

    ! **********************************
    subroutine scale_copy_(s1, s2)

        ! Interface procedure to copy

        implicit none

        type(scale), intent(in)      :: s2
        class(scale), intent(inout)  :: s1

        call scale_copy(s1, s2)

    end subroutine scale_copy_

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
    subroutine scale_stevolve(this, dt, st, inRate, outRate, pStatus, pScl)

        ! Compute current status of the scale structure
        ! according to :
        ! - Internal evolution processes
        ! - A constant input rate (due to external process) "inRate" and
        ! - A constant output rate (due to external process) "outRate"
        ! Then
        ! Apply the next step "st" of the complete integration scheme

        implicit none

        integer(kind=ikd), intent(in)  :: st         ! Step index of evolution scheme
        real(kind=rkd), intent(in)     :: dt         ! The scale is evolve during dt

        type(gas), intent(in)          :: inRate     ! The (dt-)constant input rate
        type(gas), intent(in)          :: outRate    ! The (dt-)constant output rate
        type(gas)                      :: status     ! Current status
        type(gas), pointer             :: pStatus    ! Pointer to the corrected status

        type(scale), pointer           :: pScl        ! Pointer to the current intermediate evolved scale

        class(scale)                   :: this       ! The current scale

        ! Compute current global status due to
        ! Internal process + external input/output
        status = inRate - outRate - pScl%status()
        ! Apply one more solver step and update the current intermediate stage
        pScl = this%solve(st, dt, status, pStatus)

    end subroutine scale_stevolve

    ! **********************************
    subroutine scale_evolve(this, dt, inRate, outRate)

        ! Evolve, during dt, the current scale structure 
        ! according to :
        ! - Internal evolution processes
        ! - A constant input rate (due to external process) "inRate" and
        ! - A constant output rate (due to external process) "outRate"

        implicit none

        integer(kind=ikd)              :: st         ! Step index of evolution scheme
        real(kind=rkd), intent(in)     :: dt         ! The scale is evolve during dt

        type(gas), intent(in)          :: inRate     ! The (dt-)constant input rate
        type(gas), intent(in)          :: outRate    ! The (dt-)constant output rate
        type(gas), pointer             :: pStatus    ! Pointer to the corrected status

        type(scale), target            :: scl        ! The current intermediate evolved scale
        type(scale), pointer           :: pScl       ! Pointer to the curretn evolved scale

        class(scale)                   :: this       ! The current scale

        ! Init intermediate status with the current status
        scl = this
        pScl => scl
        ! Define complete corrected status
        pStatus => sclStatus
        do st = 1, nSolverStep
            call this%stevolve(dt, st, inRate, outRate, pStatus, pScl)
        end do
        ! Save complete evolved state
        this = scl

    end subroutine scale_evolve

    !
    ! FUNCTIONS
    !

    ! **********************************
    function scale_get_N_clouds(this) result(N)

        ! Return the number of active clouds at this scale

        implicit none

        real(kind=rkd) :: N

        class(scale) :: this

        N = floor(this%gas%mass/this%BE_mass())
    
    end function scale_get_N_clouds

    ! **********************************
    function scale_get_BE_mass(this) result(M_BE)

        ! Return the Bonnot Ebert mass associated to the scale

        implicit none

        real(kind=rkd)    :: M_BE

        class(scale)    :: this

        M_BE = 3.d0/2.d0 * sigmaStar**4./GCst_CU**2./muStar*(this%l/lStar)**(4.*sigma_slope - mu_slope)
    
    end function scale_get_BE_mass

    ! **********************************
    function scale_volume(this) result(V)

        ! Return the volume of a sphere of scale l

        real(kind=rkd)   :: V

        class(scale)   :: this

        V = 4./3.*pi*(this%l/2.d0)**3.

    end function scale_volume

    ! **********************************
    function scale_get_velocity_dispersion(this) result(sV)

        ! Return the mean 1D-velocity dispersion at this scale

        implicit none

        real(kind=rkd)    :: sV

        class(scale)    :: this

        sV = sigmaStar * (this%l / lStar)**sigma_slope

    end function scale_get_velocity_dispersion

    ! **********************************
    function scale_get_mass_surface_density(this) result(mu)

        ! Return the mean mass surface density at this scale

        real(kind=rkd)   :: mu

        class(scale)   :: this

        mu = muStar * (this%l / lStar)**mu_slope

    end function scale_get_mass_surface_density

    ! **********************************
    function scale_status(this) result(rate)

        ! Return the current status of the scale
        ! according to its internal evolution only

        real(kind=rkd)         :: vRate

        type(gas)              :: rate     ! Global evolution rate

        class(scale)           :: this

        ! Transfer mass to the lower scale is only possible
        ! if the current scale is unstable (this%nClouds() > 0)
        vRate = 3.d0/2.d0 * ETRV * this%volume() / this%sV()**2. * this%nClouds()
        !
        ! The output rate gas object is build according to the gas signature
        rate = vRate * this%gas%signature()

    end function scale_status

    ! **********************************
    function scale_update(this, dt, rate) result(scl)

        ! Update the scale structure during dt according to an
        ! global evolution rate "rate"

        implicit none

        character(MAXPATHSIZE)   :: calledBy

        real(kind=rkd)           :: dt       ! The scale is evolve during dt

        type(gas), intent(in)    :: rate     ! The global evolution rate (can be negative)

        type(scale)              :: scl      ! The evolved scale

        class(scale)             :: this     ! The current scale

        ! Copy of the current scale
        scl = this
        !
        ! Evolution
        call scl%add(dt * rate)
        !
        ! Test the validity of the current scale
        write(calledBy, '(a)') 'scale_update'
        call scl%isValid(calledBy)

    end function scale_update

    ! **********************************
    function scale_solve(this, it, dt, status, pStatus) result(scl)

        ! Apply one solver step

        implicit none

        integer(kind=ikd), intent(in)   :: it       ! Current solver iteration

        real(kind=rkd)                  :: w
        real(kind=rkd), intent(in)      :: dt       ! Full time step

        type(gas), intent(in)           :: status   ! Current global evolution rate
        type(gas), intent(in), pointer  :: pStatus  ! Pointer to the corrected status

        type(scale)                     :: scl      ! The new intermediate scale evolution

        class(scale)                    :: this     ! The current scale

        select case (trim(solver))
        case ('RK4')
            ! Range-Kutta 4th order
            select case (it)
            case (1)
                ! Reset corrected status with status 1
                w = real(1.d0, kind=rkd)/real(6.d0, kind=rkd)
                pStatus = w*status
                ! Update by dt/2. using status 1
                scl = this%update(dt/real(2.d0, kind=rkd), status)
            case (2)
                ! Update corrected status with status 2
                w = real(2.d0, kind=rkd)/real(6.d0, kind=rkd)
                pStatus = pStatus + w*status
                ! Update by dt/2. using status 2
                scl = this%update(dt/real(2.d0, kind=rkd), status)
            case (3)
                ! Update corrected status with status 3
                w = real(2.d0, kind=rkd)/real(6.d0, kind=rkd)
                pStatus = pStatus + w*status
                ! Update by dt using status 3
                scl = this%update(dt, status)
            case (4)
                ! Final step
                ! Update corrected status with status 4
                w = real(1.d0, kind=rkd)/real(6.d0, kind=rkd)
                pStatus = pStatus + w*status
                ! Update by dt using complete corrected sclStatus
                scl = this%update(dt, pStatus)
            end select
        case default
            ! Range-Kutta 2d order
            select case (it)
            case (1)
                ! Reset corrected status with status 1
                pStatus = status
                ! Update by dt/2. using status 1
                scl = this%update(dt/real(2.d0, kind=rkd), status)
            case (2)
                ! Reset corrected status with status 2
                pStatus = status
                ! Update by dt using complete corrected sclStatus
                scl = this%update(dt, pStatus)
            end select
        end select

    end function scale_solve

end module scale_mod

