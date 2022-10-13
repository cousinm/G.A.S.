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
        real(kind=8)              :: mass     ! total mass  (in code unit)
        real(kind=8)              :: mZ       ! mass of metals (in code unit)
        real(kind=8)              :: Eint     ! Internal energy of the gas (in code unit)
        real(kind=8), allocatable :: elts(:)  ! mass of main ISM elements (e.g H1, He, C12, N14, O16, Fe56)
    contains
        procedure :: create         ! Create and initialize a gas object
        procedure :: metalicity     ! Return the metalicity (mass fraction) of the gas
        procedure :: signature      ! Return a gas signature.
        procedure :: abundance      ! Compute and return the abundance of a specific element
        procedure :: temperature    ! Return the mean temperature of the gas
        procedure :: molecular_mass ! Return the mean molecular mass of the gas
        procedure :: setTemperature ! Set the internal energy according to the temperature and the mass
    end type gas

    ! Define gas specific parameters
    integer(kind=4)                 :: nElts         ! Number of Main ISM elements followed
    integer(kind=4)                 :: nMetBins      ! Number of metallicity bins used
    !
    character(len=4), allocatable   :: eltNames(:)   ! List of main ISM elements followed (e.g H1, He, C12, N14, O16, Fe56)
    !
    real(kind=8), allocatable       :: metBins(:)    ! Mass fraction corresponding to metallicity bins
    type(gas), allocatable          :: initAbund(:)  ! Initial abundances (in mass fraction) associated to each metallicity bin
                                                     ! InitAbund is defined as a gas object, e.i InitAbund[Z]%mass = 1., InitAbund[Z]%mZ = metBins(Z)
                                                     ! and for each main ISM elements InitAbund[Z]%elts(e) = X0_elt
    !
    ! File unit for input parameter file "gas_properties.in"
    integer(kind=4), parameter     :: gasPropertiesUnit = 117

    ! INTERFACE OPERATOR DECLARATIONS

    interface assignment (=)  ! allows to copy a gas component by using the symbol '='
        module procedure gas_copy
    end interface assignment (=)

    interface operator (+)  ! allows to add two gas components by using the symbol '+'
        module procedure gas_add
    end interface operator (+)

    interface operator (-)  ! allows to substract a gas component to an other by using the symbol '-'
        module procedure gas_sub
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
    
        integer(kind=4)          :: j,e ! loop indexes (metallicity and elements)
    
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
            if (trim(line) .eq. 'START') then
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
                write(message, '(a)') 'Impossible to read a line in the gas_properties.in file'
                call log_message(message, &
                                 logLevel=LOG_ERROR, &
                                 calledBy='gas_read_properties')
                stop ! stop the program
            end if
        end do
2       close(gasPropertiesUnit)

        return
    end subroutine gas_read_properties

    ! **********************************
    subroutine create(this)
    
        ! Create a new gas object, initialize all fields to 0.

        class(gas) :: this
        
        this%mass = 0.  ! Total mass
        this%mZ  = 0.   ! Metal mass
        this%Eint = 0.  ! Internal energy
        
        allocate(this%elts(nElts))  ! create
        this%elts = 0.  ! Specific element mass
    
    end subroutine

    ! **********************************
    subroutine gas_copy(g1, g2)

        ! Copy the gas object g2 into the gas object g1
    
        implicit none
    
        class(gas), intent(inout) :: g1
        type(gas), intent(in)     :: g2
        
        ! Create g1
        call g1%create()

        ! Copy fields
        g1%mass = g2%mass  ! Total mass
        g1%mZ   = g2%mZ    ! Metal mass
        g1%Eint = g2%Eint  ! Internal energy
        
        ! Test if specific element array is allocated
        if (.not. allocated(g2%elts)) then
            call log_message('Try to copy an non allocated gas%elts component', &
                             logLevel=LOG_ERROR, &
                             paramNames=(/'g2%mass                  ', &
                                          'g2%mZ                    ', &
                                          'g2%Eint                  '/), &
                             realParams=(/g2%mass,g2%mZ,g2%Eint/), &
                             calledBy='gas_copy')
        end if
        g1%elts = g2%elts  ! copy

        return
    end subroutine gas_copy

    ! **********************************
    subroutine setTemperature(this, T)

        ! Set the gas total internal energy
        ! according to the temperature and the total mass

        implicit none

        class(gas)                :: this

        real(kind=8), intent(in)  :: T      ! The temperature [K]
        real(kind=8)              :: eint   ! Internal energy per particle
        real(kind=8)              :: mu     ! mean molecular mass

        ! The temperature allows to define the internal energy per particle
        eint = 3.d0 / 2.d0 * kb * T  ! [Joule]

        ! Set the total internal energy according to the mass and 
        ! the mean molecular mass

        mu = this%molecular_mass()
        if (mu > 0.) then
            this%Eint = eint * this%mass * Mass_kg / (mu * mp)  ! [Joule]
            this%Eint = this%Eint * Energy_CU                   ! [CU]
        else
            call log_message('Unvalid mean molecular mass', &
                             logLevel=LOG_ERROR, &
                             paramNames=(/'mu                       ', &
                                          'this%mass                ', &
                                          'this%mZ                  ', &
                                          'this%Eint                '/), &
                             realParams=(/mu,this%mass,this%mZ,this%Eint/), &
                             calledBy='setTemperature')
        end if

    end subroutine setTemperature

    !
    ! FUNCTIONS
    !

    ! **********************************
    function metalicity(this) result(mZ)

        ! Return the gas metalicity (metal mass fraction)

        implicit none

        class(gas), intent(in)   :: this

        real(kind=8)             :: mZ

        mZ = this%mZ / this%mass

    end function metalicity

    ! **********************************
    function signature(this) result(g)

        ! Return a gas signature: 
        !   i.e a gas object normalized to gas%mass = 1.

        implicit none

        class(gas), intent(in)   :: this
        type(gas)                :: g

        g = this
        g = 1.d0/this%mass * g

    end function signature

    ! **********************************
    function abundance(this, component) result(X)

        ! Return the abundance (in mass ratio)
        ! of a given gas specific element

        implicit none

        integer(kind=4)          :: e

        character(len=*)         :: component

        real(kind=8)             :: X

        class(gas), intent(in)   :: this

        ! Test if the total mass is valid
        if (this%mass <= 0.d0) then
            call log_message('Try to compute abundance of a null mass gas', &
                             logLevel=LOG_ERROR, &
                             paramNames=(/'this%mass                ', &
                                          'this%mZ                  ', &
                                          'this%Eint                '/), &
                             realParams=(/this%mass,this%mZ,this%Eint/), &
                             calledBy='gas%abundance()')  
        end if

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

        real(kind=8)             :: X, Y
        real(kind=8)             :: mu

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

        real(kind=8)             :: T, N, mu

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
    function gas_add(g1, g2) result(g)

        ! Add the gas object g2 to the gas object g1
        ! For this procedure, g1 AND g2 havr to be both created
    
        implicit none
    
        type(gas), intent(in)     :: g1
        type(gas), intent(in)     :: g2
        type(gas)                 :: g

        ! Test if specific element array is allocated
        ! for g1
        if (.not. allocated(g1%elts)) then
            call log_message('Try to use an non allocated gas%elts component', &
                             logLevel=LOG_ERROR, &
                             paramNames=(/'g1%mass                  ', &
                                          'g1%mZ                    ', &
                                          'g1%Eint                  '/), &
                             realParams=(/g1%mass,g1%mZ,g1%Eint/), &
                             calledBy='gas_add')
        end if
        ! and for g2
        if (.not. allocated(g2%elts)) then
            call log_message('Try to use an non allocated gas%elts component', &
                             logLevel=LOG_ERROR, &
                             paramNames=(/'g2%mass                  ', &
                                          'g2%mZ                    ', &
                                          'g2%Eint                  '/), &
                             realParams=(/g2%mass,g2%mZ,g2%Eint/), &
                             calledBy='gas_add')
        end if
        
        ! Create the return gas object
        call g%create()
        ! Add
        ! Add fields from input gas objects
        g%mass = g1%mass + g2%mass  ! Total mass
        g%mZ   = g1%mZ + g2%mZ      ! Metal mass
        g%Eint = g1%Eint + g2%Eint  ! Internal energy
        g%elts = g1%elts + g2%elts  ! Specific element mass

        return
    end function gas_add

    ! **********************************
    function gas_sub(g1, g2) result(g)

        ! Substract a gas object g2 to the gas object g1 (g = g1 - g2)

        implicit none

        integer(kind=4)           :: e

        character(MAXPATHSIZE)    :: message

        type(gas), intent(in)     :: g1
        type(gas), intent(in)     :: g2
        type(gas)                 :: g

        ! Test total mass of components
        if ((g2%mass .gt. g1%mass) .and. (abs(g2%mass - g1%mass) .gt. num_accuracy)) then
            call log_message('Substract too much gas', &
                             logLevel=LOG_ERROR, &
                             paramNames=(/'g1%mass                  ', &
                                          'g1%mZ                    ', &
                                          'g1%Eint                  ', &
                                          'g2%mass                  ', &
                                          'g2%mZ                    ', &
                                          'g2%Eint                  '/), &
                             realParams=(/g1%mass,g1%mZ,g1%Eint,g2%mass,g2%mZ,g2%Eint/), &
                             calledBy='gas_sub')
        end if

        ! Test metal mass of components
        if ((g2%mZ .gt. g1%mZ) .and. (abs(g2%mZ - g1%mZ) .gt. num_accuracy)) then
            call log_message('Substract too much gas', &
                             logLevel=LOG_ERROR, &
                             paramNames=(/'g1%mass                  ', &
                                          'g1%mZ                    ', &
                                          'g1%Eint                  ', &
                                          'g2%mass                  ', &
                                          'g2%mZ                    ', &
                                          'g2%Eint                  '/), &
                             realParams=(/g1%mass,g1%mZ,g1%Eint,g2%mass,g2%mZ,g2%Eint/), &
                             calledBy='gas_sub')
        end if

        ! Test total internal energy of components
        if ((g2%Eint .gt. g1%Eint) .and. (abs(g2%Eint - g1%Eint) .gt. num_accuracy)) then
            call log_message('Substract too much total internal energy', &
                             logLevel=LOG_ERROR, &
                             paramNames=(/'g1%mass                  ', &
                                          'g1%mZ                    ', &
                                          'g1%Eint                  ', &
                                          'g2%mass                  ', &
                                          'g2%mZ                    ', &
                                          'g2%Eint                  '/), &
                             realParams=(/g1%mass,g1%mZ,g1%Eint,g2%mass,g2%mZ,g2%Eint/), &
                             calledBy='gas_sub')
        end if

        ! Test if specific element array is allocated
        ! for g1
        if (.not. allocated(g1%elts)) then
            call log_message('Try to use an non allocated gas%elts component', &
                             logLevel=LOG_ERROR, &
                             paramNames=(/'g1%mass                  ', &
                                          'g1%mZ                    ', &
                                          'g1%Eint                  '/), &
                             realParams=(/g1%mass,g1%mZ,g1%Eint/), &
                             calledBy='gas_sub')
        end if
        ! and for g2
        if (.not. allocated(g2%elts)) then
            call log_message('Try to use an non allocated gas%elts component', &
                             logLevel=LOG_ERROR, &
                             paramNames=(/'g2%mass                  ', &
                                          'g2%mZ                    ', &
                                          'g2%Eint                  '/), &
                             realParams=(/g2%mass,g2%mZ,g2%Eint/), &
                             calledBy='gas_sub')
        end if

        ! @ this point specific element array are available
        do e = 1, nElts
            if ((g1%elts(e) < g2%elts(e)) .and. (abs(g2%elts(e) - g1%elts(e)) .gt. num_accuracy)) then
                write(message, '(a,a)') 'Substract too much gas for the specific element e = ', eltNames(e)
                call log_message(message, &
                                 logLevel=LOG_ERROR, &
                                 paramNames=(/'g1%mass                  ', &
                                              'g1%mZ                    ', &
                                              'g1%Eint                  ', &
                                              'g1%elts(e)               ', &
                                              'g2%mass                  ', &
                                              'g2%mZ                    ', &
                                              'g2%Eint                  ', &
                                              'g2%elts(e)               '/), &
                                 realParams=(/g1%mass,g1%mZ,g1%Eint,g1%elts(e), &
                                              g2%mass,g2%mZ,g2%Eint,g2%elts(e)/), &
                                 calledBy='gas_sub')
                
            end if
        end do

        ! @ this point all checks have been done, we can substract g2 to g1
        ! Create the return gas object
        call g%create()
        ! Substract
        g%mass = max(0., g1%mass - g2%mass)  ! Total mass
        g%mZ   = max(0., g1%mZ - g2%mZ)      ! Metal mass
        g%Eint = max(0., g1%Eint - g2%Eint)  ! Internal energy
        g%elts = max(0., g1%elts - g2%elts)  ! Specific element mass

    end function gas_sub

    ! **********************************
    function gas_scalar_multiply_left(g, a) result(ga)

        ! Return gas * a
    
        implicit none
    
        real(kind=8), intent(in)  :: a
    
        type(gas), intent(in)     :: g
        type(gas)                 :: ga

        ! Create the return gas object
        call ga%create()
    
        ga%mass  = g%mass * a
        ga%mZ    = g%mZ   * a
        ga%Eint  = g%Eint * a
    
        if (allocated(g%elts)) then
           ga%elts = g%elts * a
        else
            call log_message('Try to use an non allocated gas%elts component', &
                             logLevel=LOG_ERROR, &
                             paramNames=(/'g%mass                   ', &
                                          'g%mZ                     ', &
                                          'g%Eint                   '/), &
                            realParams=(/g%mass,g%mZ,g%Eint/), &
                            calledBy='gas_scalar_multiply_left')
        endif
    end function gas_scalar_multiply_left

    ! **********************************
    function gas_scalar_multiply_right(a, g) result(ag)

        ! Return a * gas
    
        implicit none
    
        real(kind=8), intent(in)  :: a
    
        type(gas), intent(in)     :: g
        type(gas)                 :: ag

        ! Use comutativity property
        ag = gas_scalar_multiply_left(g, a)

    end function gas_scalar_multiply_right

end module gas_mod
