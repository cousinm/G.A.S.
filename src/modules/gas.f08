module gas_mod

    !*****************************************************************************************************************
    !
    !  OVERVIEW
    !
    !  The gas module defines all properties and procedures assossiated to a gas object.
    !  In G.A.S. a gas phase is defined by its total mass, metals mass
    !     and also the mass of six main ISM elements (H, He, C12, N14, O16 and Fe56)
    !  This module starts with the definition of this gas type
    !
    !*****************************************************************************************************************

    use parameters  ! Acces to global defintions and properties
    use log_mod     ! Acces to logging procedures
    use config_mod  ! Acces to configurations parameters (path)

    implicit none
    
    public

    ! GAS TYPE DEFINITION

    type gas
        real(kind=rkd)              :: mass       ! total mass  (in code unit)
        real(kind=rkd)              :: mZ         ! mass of metals (in code unit)
        real(kind=rkd)              :: Eint       ! Internal energy of the gas (in code unit)
        real(kind=rkd), allocatable :: elts(:)    ! mass of main ISM elements (e.g H1, He, C12, N14, O16, Fe56)
    contains
        procedure :: create => gas_create       ! Create and initialize a gas object
        procedure :: delete => gas_delete       ! Delete a gas object
        procedure :: isCreated => gas_isCreated ! Test if gas object is already created
        procedure :: isValid => gas_isValid     ! Test concistency of the gas object
        procedure :: copy => gas_copy           ! Copy a gas object from an other
        procedure :: add => gas_add             ! Add a gas object
        procedure :: sub => gas_sub             ! Substract a gas object
        procedure :: metalicity                 ! Return the metalicity (mass fraction) of the gas
        procedure :: signature                  ! Return a gas signature.
        procedure :: abundance                  ! Compute and return the abundance of a specific element
        procedure :: temperature                ! Return the mean temperature of the gas
        procedure :: molecular_mass             ! Return the mean molecular mass of the gas
        procedure :: setTemperature             ! Set the internal energy according to the temperature and the mass
    end type gas

    ! Define gas specific parameters
    integer(kind=ikd)               :: nElts         ! Number of Main ISM elements followed
    integer(kind=ikd)               :: nMetBins      ! Number of metallicity bins used
    !
    character(len=4), allocatable   :: eltNames(:)   ! List of main ISM elements followed (e.g H1, He, C12, N14, O16, Fe56)
    !
    real(kind=rkd), allocatable     :: metBins(:)    ! Mass fraction corresponding to metallicity bins
    type(gas), allocatable          :: initAbund(:)  ! Initial abundances (in mass fraction) associated to each metallicity bin
                                                     ! InitAbund is defined as a gas object, e.i InitAbund[Z]%mass = 1., InitAbund[Z]%mZ = metBins(Z)
                                                     ! and for each main ISM elements InitAbund[Z]%elts(e) = X0_elt

    ! INTERFACE OPERATOR DECLARATIONS

    interface assignment (=)  ! allows to copy a gas component by using the symbol '='
        module procedure gas_copy
    end interface assignment (=)

    interface operator (+)  ! allows to add two gas components by using the symbol '+'
        module procedure add
    end interface operator (+)

    interface operator (-)  ! allows to substract a gas component to an other by using the symbol '-'
        module procedure sub
    end interface operator (-)

    interface operator (*)  ! allows to ponderate a gas component (a* or *a) by using the symbol '*'
        module procedure gas_scalar_multiply_left, gas_scalar_multiply_right
    end interface operator (*)

contains
    
    !
    ! SUBROUTINES
    !

    ! **********************************
    subroutine gas_init()

        ! Initialize the gas module
        !   Load gas properties

        implicit none

        call gas_read_properties()
    
    end subroutine gas_init

    ! **********************************
    subroutine gas_read_properties()

        ! Load gas properties: abundance, tracking elements name ...

        implicit none
    
        integer(kind=ikd)          :: j,e ! loop indexes (metallicity and elements)
    
        character(MAXPATHSIZE)   :: filename
        character(MAXPATHSIZE)   :: line
        character(MAXPATHSIZE)   :: message
    
        write(filename,'(a,a,a,a)') trim(librariesPath), '/gas_properties.in'
        write(message,'(a,a,a,a)') 'Load data from: ', trim(filename)
        call log_message(message)
        !
        ! Open the "gas properties" file
        open(unit=gasPropertiesUnit, file=trim(filename), status='old')
        ! Read and load data
        do
            read(gasPropertiesUnit, '(a)', end=2) line
            if (trim(line) .eq. '----') then
                ! The number main ISM elements followed, the number of Metallicity bin used
                read(gasPropertiesUnit, *) nElts, nMetBins 
                call log_message('gas properties used: ', &
                                  paramNames=(/'nElts                    ','nMetBins                 '/), &
                                  intParams=(/nElts,nMetBins/))
                !
                ! Create the table of the Main ISM element names
                allocate(eltNames(nElts))
                ! Create the metallicity bin table
                allocate(metBins(nMetBins))
                ! Create the initial abundance table
                allocate(initAbund(nMetBins))
                !
                ! Read Main ISM elements list
                read(gasPropertiesUnit, *) (eltNames(e), e=1, nElts)
                write(message,'(a)') 'Chemical evolution is followed for: '
                
                do e = 1, nElts
                    write(message,'(a,a,a)') trim(message), trim(eltNames(e)), ','
                end do
                call log_message(message)
                ! Read metallicity bins and initial abundance table
                ! InitAbund is defined as a gas object,
                ! InitAbund[Z]%mass = 1., InitAbund[Z]%mZ = metBins(Z)
                ! and for each main ISM elements InitAbund[Z]%elts(e) = X0_elt
                do j = 1, nMetBins
                    call InitAbund(j)%create()
                    InitAbund(j)%mass = 1.d0
                    read(gasPropertiesUnit, *) MetBins(j), (InitAbund(j)%elts(e), e=1,nElts)
                    InitAbund(j)%mZ   = MetBins(j)
                end do
                exit  ! quit do loop
            end if
            if (line(1:1) .eq. '#') then
                cycle ! header or something like this (skip)
            else
                write(message, '(a, a)') 'Impossible to read a line in ', trim(filename)
                call log_message(message, &
                                 logLevel=LOG_ERROR, &
                                 calledBy='gas_read_properties')
            end if
        end do
2       close(gasPropertiesUnit)

        return
    end subroutine gas_read_properties

    ! **********************************
    subroutine gas_create(this)
    
        ! Create a new gas object, initialize all fields to 0.

        class(gas) :: this
        
        this%mass = 0.  ! Total mass
        this%mZ  = 0.   ! Metal mass
        this%Eint = 0.  ! Internal energy

        if (.not. allocated(this%elts)) then
            allocate(this%elts(nElts))  ! create
        end if
        this%elts = 0.  ! Specific element mass
    
    end subroutine gas_create

    ! **********************************
    subroutine gas_delete(this)
    
        ! Delete a gas object

        class(gas) :: this
        
        this%mass = 0.  ! Total mass
        this%mZ  = 0.   ! Metal mass
        this%Eint = 0.  ! Internal energy
        
        if (allocated(this%elts)) then
            deallocate(this%elts)  ! Delete
        end if
    
    end subroutine gas_delete

    ! **********************************
    subroutine gas_isValid(this, calledBy)

        ! Test a gas component
        !    gas%elts should be allocated
        !    Mass should be > 0. for all components
        !    Total internal energy should be > 0. too

        integer(kind=ikd)                     :: e

        character(MAXPATHSIZE)              :: message
        character(MAXPATHSIZE), intent(in)  :: calledBy

        class(gas)                          :: this

        ! Test mass
        if (this%mass < 0.) then
            write(message, '(a)') 'Current mass is not valid.'
            call log_message(message, &
                             logLevel=LOG_ERROR, &
                             paramNames=(/'this%mass                ', &
                                          'this%mZ                  ', &
                                          'this%Eint                '/), &
                             realParams=(/this%mass,this%mZ,this%Eint/), &
                             calledBy=calledBy)
        end if
        !
        ! Test metal mass
        if (this%mZ < 0.) then
            write(message, '(a)') 'Current metal mass is not valid.'
            call log_message(message, &
                             logLevel=LOG_ERROR, &
                             paramNames=(/'this%mass                ', &
                                          'this%mZ                  ', &
                                          'this%Eint                '/), &
                             realParams=(/this%mass,this%mZ,this%Eint/), &
                             calledBy=calledBy)
        end if
        !
        ! Test if elts is allocated
        if (.not. allocated(this%elts)) then
            write(message, *) 'Structure this%elts is not allocated'
            call log_message(message, &
                             logLevel=LOG_ERROR, &
                             paramNames=(/'this%mass                ', &
                                          'this%mZ                  ', &
                                          'this%Eint                '/), &
                             realParams=(/this%mass,this%mZ,this%Eint/), &
                             calledBy=calledBy)
        end if
        !
        ! Test specific element mass
        do e = 1, nElts
            if (this%elts(e) < 0.) then
                write(message, '(a,a,a)') 'Current ', eltNames(e), ' mass is not valid.'
                call log_message(message, &
                                 logLevel=LOG_ERROR, &
                                 paramNames=(/'this%mass                ', &
                                              'this%mZ                  ', &
                                              'this%Eint                ', &
                                              'this%elts(e)             '/), &
                                 realParams=(/this%mass,this%mZ,this%Eint,this%elts(e)/), &
                                 calledBy=calledBy)
                
            end if
        end do
        !
        ! Test internal Energy
        if (this%Eint < 0.) then
            write(message, '(a)') 'Current total internal ernergy is not valid.'
            call log_message(message, &
                             logLevel=LOG_ERROR, &
                             paramNames=(/'this%mass                ', &
                                          'this%mZ                  ', &
                                          'this%Eint                '/), &
                             realParams=(/this%mass,this%mZ,this%Eint/), &
                             calledBy=calledBy)
        end if

    end subroutine gas_isValid

    ! **********************************
    subroutine gas_copy(g1, g2)

        ! Copy the gas object g2 into the gas object g1
    
        implicit none

        type(gas), intent(in)            :: g2

        class(gas), intent(inout)        :: g1

        if (.not. g1%isCreated()) then
            ! Create g1
            call g1%create()
        end if

        ! Copy fields
        g1%mass = g2%mass  ! Total mass
        g1%mZ   = g2%mZ    ! Metal mass
        g1%Eint = g2%Eint  ! Internal energy
        g1%elts = g2%elts  ! Specific element mass

        return
    end subroutine gas_copy

    ! **********************************
    subroutine gas_add(this, g)

        ! Add the gas object g to the current gas object (this = this - g)
        ! For this procedure, this AND g have to be both created
    
        implicit none

        type(gas), intent(in)  :: g

        class(gas)             :: this

        ! Add
        ! Add fields from input gas objects
        this%mass = this%mass + g%mass  ! Total mass
        this%mZ   = this%mZ + g%mZ      ! Metal mass
        this%Eint = this%Eint + g%Eint  ! Internal energy
        this%elts = this%elts + g%elts  ! Specific element mass

    end subroutine gas_add

    ! **********************************
    subroutine gas_sub(this, g)

        ! Substract a gas object g to the current gas object (this = this - g)
        ! For this procedure, this AND g have to be both created

        implicit none

        type(gas), intent(in)            :: g

        class(gas)                       :: this

        ! @ this point all checks have been done, we can substract g to this
        ! Substract
        this%mass = this%mass - g%mass  ! Total mass
        this%mZ   = this%mZ - g%mZ      ! Metal mass
        this%Eint = this%Eint - g%Eint  ! Internal energy
        this%elts = this%elts - g%elts  ! Specific element mass

    end subroutine gas_sub

    ! **********************************
    subroutine setTemperature(this, T)

        ! Set the gas total internal energy
        ! according to the temperature and the total mass

        implicit none

        character(MAXPATHSIZE)           :: calledBy

        real(kind=rkd), intent(in)         :: T      ! The temperature [K]
        real(kind=rkd)                     :: eint   ! Internal energy per particle
        real(kind=rkd)                     :: mu     ! mean molecular mass

        class(gas)                       :: this

        write(calledBy, '(a)') 'setTemperature'

        ! The temperature allows to define the internal energy per particle
        eint = 3.d0 / 2.d0 * kb * T  ! [Joule]

        ! Set the total internal energy according to the mass and 
        ! the mean molecular mass

        mu = this%molecular_mass()
        this%Eint = eint * this%mass * Mass_kg / (mu * mp)  ! [Joule]
        this%Eint = this%Eint * Energy_CU                   ! [CU]
    end subroutine setTemperature

    !
    ! FUNCTIONS
    !

    ! **********************************
    function gas_isCreated(this) result(isValid)

        ! Test the validity of a gas object

        logical     :: isValid

        class(gas)  :: this

        isValid = allocated(this%elts)
    end function gas_isCreated

    ! **********************************
    function add(g1, g2) result(g)

        implicit none

        type(gas), intent(in)  :: g1
        type(gas), intent(in)  :: g2
        type(gas)              :: g

        ! Create the output gas object
        call g%create()
        call g%add(g1)
        call g%add(g2)

    end function add

    ! **********************************
    function sub(g1, g2) result(g)

        implicit none

        type(gas), intent(in)  :: g1
        type(gas), intent(in)  :: g2
        type(gas)              :: g

        ! Create the output gas object
        call g%create()
        call g%add(g1)
        call g%sub(g2)

    end function sub

    ! **********************************
    function metalicity(this) result(mZ)

        ! Return the gas metalicity (metal mass fraction)

        implicit none

        class(gas), intent(in)   :: this

        real(kind=rkd)             :: mZ

        mZ = this%mZ / this%mass

    end function metalicity

    ! **********************************
    function signature(this) result(g)

        ! Return a gas signature: 
        !   i.e a gas object normalized to gas%mass = 1.

        implicit none

        type(gas)                :: g

        class(gas), intent(in)   :: this

        if (this%mass > 0.d0) then
            ! Compute gas signature
            g = 1.d0/this%mass * this
        else
            g = this
        end if

    end function signature

    ! **********************************
    function abundance(this, component) result(X)

        ! Return the abundance (in mass ratio)
        ! of a given gas specific element

        implicit none

        integer(kind=ikd)          :: e

        character(len=*)         :: component

        real(kind=rkd)             :: X

        class(gas), intent(in)   :: this

        ! Loop oven specific element and compute abundance
        X = -1.d0
        do e = 1, nElts
            if (trim(component) .eq. trim(eltNames(e))) then
                X = this%elts(e) / this%mass
            exit
            end if
        end do

        ! Test ouput value
        if (X <= 0.d0) then
            call log_message('Abundance is not valid', &
                             logLevel = LOG_ERROR, &
                             paramNames=(/'X                        ', &
                                          'this%mass                ', &
                                          'this%mZ                  ', &
                                          'this%Eint                '/), &
                             realParams=(/X,this%mass,this%mZ,this%Eint/), &
                             calledBy='gas%abundance()')  
        end if

    end function abundance

    ! **********************************
    function molecular_mass(this) result(mu)

        ! Return the gas mean molecular mass
        ! according to the metalicity

        implicit none

        class(gas), intent(in)   :: this

        real(kind=rkd)             :: X, Y
        real(kind=rkd)             :: mu

        X = this%abundance('H1')
        Y = this%abundance('He4')
        mu = 4.d0 / (6.d0*X + Y + 2.d0)

        if (mu <= 0.) then
            call log_message('Mean molecular mass is not valid', &
                             logLevel=LOG_ERROR, &
                             paramNames=(/'mu                       ', &
                                          'this%mass                ', &
                                          'this%mZ                  ', &
                                          'this%Eint                '/), &
                             realParams=(/mu,this%mass,this%mZ,this%Eint/), &
                             calledBy='gas%molecular_mass')
        end if

    end function molecular_mass

    ! **********************************
    function temperature(this) result(T)

        ! Return the mean gas temperature

        implicit none

        real(kind=rkd)             :: T, N, mu

        class(gas), intent(in)   :: this

        ! Temperature is based on internal energy by particle
        mu = this%molecular_mass()
        N = this%mass * Mass_kg / (mu*mp)

        T = 2.d0*this%Eint*Energy_J / (3.d0*kb*N)  ! [K]

        if (T <= 0.) then
            call log_message('Temperature is not valid', &
                             logLevel=LOG_ERROR, &
                             paramNames=(/'T                        ', &
                                          'this%mass                ', &
                                          'this%mZ                  ', &
                                          'this%Eint                '/), &
                             realParams=(/T,this%mass,this%mZ,this%Eint/), &
                             calledBy='gas%temperature')
        end if

    end function temperature

    ! **********************************
    function gas_scalar_multiply_left(g, a) result(ga)

        ! Return gas * a
    
        implicit none

        real(kind=rkd), intent(in)  :: a
    
        type(gas), intent(in)       :: g
        type(gas)                   :: ga

        ! Create the return gas object
        call ga%create()
    
        ga%mass = g%mass * a
        ga%mZ   = g%mZ   * a
        ga%Eint = g%Eint * a
        ga%elts = g%elts * a

    end function gas_scalar_multiply_left

    ! **********************************
    function gas_scalar_multiply_right(a, g) result(ag)

        ! Return a * gas
    
        implicit none
    
        real(kind=rkd), intent(in)  :: a
    
        type(gas), intent(in)       :: g
        type(gas)                   :: ag

        ! Use comutativity property
        ag = gas_scalar_multiply_left(g, a)

    end function gas_scalar_multiply_right

end module gas_mod
