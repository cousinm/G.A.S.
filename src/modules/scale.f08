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
        procedure   :: evolve => scale_evolve               ! Evolve the scale during dt
        procedure   :: solve => scale_solve                 ! Solve the input output system during dt
    end type scale

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
    subroutine scale_solve(this, dt, inRate)

        ! Solve the evolution of the current scale structure during dt
        !  according to a constant input rate "inRate"

        implicit none

        real(kind=rkd)        :: dt       ! The scale is evolve during dt

        type(gas), intent(in) :: inRate   ! The (dt-)constant input rate
        type(gas)             :: outRate  ! The corrected output rate
        type(gas)             :: outRate1 ! Intermediates output rates
        type(gas)             :: outRate2
        type(gas)             :: outRate3
        type(gas)             :: outRate4
        type(scale)           :: scl_tmp  ! Intermediate states of the scale

        class(scale)          :: this     ! The current scale

        ! Get curent status
        outRate1 = this%status()
        ! Performed evolution for dt/2
        scl_tmp = this%evolve(dt/2.d0, inRate, outRate1)
        !
        select case (trim(intScheme))
            !
            case ('RK4')
                ! Range-Kutta 4th order
                ! Get curent status
                outRate2 = scl_tmp%status()
                ! Get second intermediate state
                scl_tmp = this%evolve(dt/2.d0, inRate, outRate2)
                outRate3 = scl_tmp%status()
                ! Get last intermediate state
                scl_tmp = this%evolve(dt, inRate, outRate3)
                outRate4 = scl_tmp%status()
                ! Perform complete (corrected) evolution
                outRate = real(1.d0/6.d0,kind=rkd)*(outRate1 + &
                                     real(2.d0,kind=rkd)*outRate2 + &
                                     real(2.d0,kind=rkd)*outRate3 + outRate4)
                ! Final evolution with corrected output rate
                this = this%evolve(dt, inRate, outRate)
            case default
                ! Range-Kutta 2d order
                ! Get intermedite status
                outRate = scl_tmp%status()
                ! Perform complete (corrected) evolution
                this = this%evolve(dt, inRate, outRate)
        end select
        !
        ! Delete temporary gas object
        call outRate%delete()
        call outRate1%delete()
        call outRate2%delete()
        call outRate3%delete()
        call outRate4%delete()
        ! And scale temporary object
        call scl_tmp%delete()

    end subroutine scale_solve

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

        ! Return the current output rate of the scale

        real(kind=rkd)   :: vRate

        type(gas)        :: rate

        class(scale)     :: this

        ! Transfer mass to the lower scale is only possible
        ! if the current scale is unstable (this%nClouds() > 0)
        
        vRate = 3.d0/2.d0 * ETRV * this%volume() / this%sV()**2. * this%nClouds()
        !
        ! The output rate gas object is build according to the gas signature
        rate = vRate * this%gas%signature()

    end function scale_status

    ! **********************************
    function scale_evolve(this, dt, inRate, outRate) result(scl)

        ! Evolve the scale structure during dt according to an
        ! global evolution rate "rate"

        implicit none

        character(MAXPATHSIZE) :: calledBy !

        real(kind=rkd)           :: dt       ! The scale is evolve during dt

        type(gas), intent(in)  :: inRate   ! The input rate
        type(gas), intent(in)  :: outRate  ! The output rate

        class(scale)           :: this     ! The current scale
        type(scale)            :: scl      ! The evolved scale

        ! Copy of the current scale
        scl = this
        ! Evolution
        call scl%add(dt*inRate)  ! Add input
        call scl%sub(dt*outRate) ! Substract output
        !
        ! Test the validity of the current scale
        write(calledBy, '(a)') 'scale_evolve'
        call scl%isValid(calledBy)

    end function scale_evolve

end module scale_mod

