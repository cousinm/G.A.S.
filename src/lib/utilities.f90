module utilities

  use global_variables  ! Acces to global definition and properties
  use IO                ! Input/Output definitions, display procedures and reading procedures for input files

  !*****************************************************************************************************************
  !
  ! OVERVIEW
  !
  ! General purpose (e.g. interpolation) routines used by (but not specifically related to) the G.A.S model
  !
  ! NO SUBROUTINES IN THIS MODULE
  !
  ! FUNCTION IN THIS MODULE
  !
  ! dt2dz                            :: compute redshift interval corresponding to a time interval at a given redshift
  !
  ! is_NaN                           :: Allows to test if a number is not defined correctly (Not a Number)
  !                                     Return .true. if x is not a number
  !
  ! BESSI0                           :: "I" Bessel function of rank 0
  !
  ! BESSI1                           :: "I" Bessel function of rank 1
  !
  ! BESSK0                           :: "K" Bessel function of rank 0
  !
  ! BESSK1                           :: "K" Bessel function of rank 1
  !
  ! slope_limiter(x,x0)              :: Reduce the evolution rate of a given properties (x)
  !                                        allows to slowly approach the resolution limit (x0)
  !
  ! Ronbint(f, a, b, err, param ...) :: Compute the integral of the function "f" between "a" and "b"
  !                                        The function "f" can be defined with some parameters given in the input array "param"
  !
  ! trap                             :: simple trapezium rule to integrate y over a pre-tabulated x array
  !
  ! Normal_distribution(mean,sigma)  :: Returns a random variable "x" following a normal distribution of mean "mean" and of standard deviation "sigma"
  !
  ! erf_func(x)                      :: Compute the standard error function
  !
  ! locate2D                         :: Returns the value stored in pre-tabulted 2D-table
  !
  !*****************************************************************************************************************

  contains

  !***********************************************************************************************************

  function dt2dz(dt,z)

    ! COMPUTE THE 'dz' KNOWING THE 'dt'
    ! at a given z

    implicit none

    real(kind=8),intent(in)      :: dt, z   ! in Gyr and without unit
    real(kind=8)                 :: dt2dz   ! without unit

    dt2dz = (Omega_m*z+1.d0)*(1.d0+z)**2. - Omega_L*z*(z+2.d0)
    dt2dz = -1.d2*h_0_code_unit*(1.d0+z)*sqrt(dt2dz)*dt

    return
  end function dt2dz

  !***********************************************************************************************************

  function is_NaN(x)

    ! Return .true. if x is not a number

    implicit none

    logical         :: is_NaN ! is true if

    real(kind=8)    :: x      ! is not a number

    if (x == x) then
      !
      is_NaN = .false.
    else
      !
      is_NaN = .true.
    endif

    return

  end function is_NaN

  !*************************************************************************************************************

  pure function BESSI0(X)

    ! "I" Bessel function of rank 0

    implicit none

    real(kind=8),intent(in) :: X
    real(kind=8)            :: BESSI0,Y,AX,BX
    real(kind=8),parameter  :: P1 = 1.d0
    real(kind=8),parameter  :: P2 = 3.5156229d0
    real(kind=8),parameter  :: P3 = 3.0899424d0
    real(kind=8),parameter  :: P4 = 1.2067429d0
    real(kind=8),parameter  :: P5 = 0.2659732d0
    real(kind=8),parameter  :: P6 = 0.360768d-1
    real(kind=8),parameter  :: P7 = 0.45813d-2
    real(kind=8),parameter  :: Q1 = 0.39894228d0
    real(kind=8),parameter  :: Q2 = 0.1328592d-1
    real(kind=8),parameter  :: Q3 = 0.225319d-2
    real(kind=8),parameter  :: Q4 = -0.157565d-2
    real(kind=8),parameter  :: Q5 = 0.916281d-2
    real(kind=8),parameter  :: Q6 = -0.2057706d-1
    real(kind=8),parameter  :: Q7 = 0.2635537d-1
    real(kind=8),parameter  :: Q8 = -0.1647633d-1
    real(kind=8),parameter  :: Q9 = 0.392377d-2

    if(abs(X) .lt. 3.75d0) then
      !
      Y=(X/3.75D0)**2
      BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
    else
      !
      AX=abs(X)
      Y=3.75D0/AX
      BX=exp(AX)/sqrt(AX)
      AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
      BESSI0=AX*BX
    endif

    return

  end function BESSI0

  !*****************************************************************************************************************

  pure function BESSI1(X)

    ! "I" Bessel function of rank 1

    implicit none

    real(kind=8), intent(in)  :: X
    real(kind=8)              :: BESSI1,Y,AX,BX
    real(kind=8),parameter    :: P1 = 0.5d0
    real(kind=8),parameter    :: P2 = 0.87890594d0
    real(kind=8),parameter    :: P3 = 0.51498869d0
    real(kind=8),parameter    :: P4 = 0.15084934d0
    real(kind=8),parameter    :: P5 = 0.2658733d-1
    real(kind=8),parameter    :: P6 = 0.301532d-2
    real(kind=8),parameter    :: P7 = 0.32411d-3
    real(kind=8),parameter    :: Q1 = 0.39894228d0
    real(kind=8),parameter    :: Q2 = -0.3988024d-1
    real(kind=8),parameter    :: Q3 = -0.362018d-2
    real(kind=8),parameter    :: Q4 = 0.163801d-2
    real(kind=8),parameter    :: Q5 = -0.1031555d-1
    real(kind=8),parameter    :: Q6 = 0.2282967d-1
    real(kind=8),parameter    :: Q7 = -0.2895312d-1
    real(kind=8),parameter    :: Q8 = 0.1787654d-1
    real(kind=8),parameter    :: Q9 = -0.420059d-2

    if(abs(X) .lt. 3.75D0) then
      !
      Y=(X/3.75D0)**2
      BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
    else
      !
      AX=abs(X)
      Y=3.75D0/AX
      BX=exp(AX)/sqrt(AX)
      AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
      BESSI1=AX*BX
    endif

    return

  end function BESSI1

  !*****************************************************************************************************************

  pure function BESSK0(X)

    ! "K" Bessel function of rank 0

    implicit none

    real(kind=8), intent(in)  :: X
    real(kind=8)              :: BESSK0,Y,AX
    real(kind=8),parameter    :: P1 = -0.57721566d0
    real(kind=8),parameter    :: P2 = 0.42278420d0
    real(kind=8),parameter    :: P3 = 0.23069756d0
    real(kind=8),parameter    :: P4 = 0.3488590d-1
    real(kind=8),parameter    :: P5 = 0.262698d-2
    real(kind=8),parameter    :: P6 = 0.10750d-3
    real(kind=8),parameter    :: P7 = 0.74d-5
    real(kind=8),parameter    :: Q1 = 1.25331414d0
    real(kind=8),parameter    :: Q2 = -0.7832358d-1
    real(kind=8),parameter    :: Q3 = 0.2189568d-1
    real(kind=8),parameter    :: Q4 = -0.1062446d-1
    real(kind=8),parameter    :: Q5 = 0.587872d-2
    real(kind=8),parameter    :: Q6 = -0.251540d-2
    real(kind=8),parameter    :: Q7 = 0.53208d-3

    if(X .eq. 0.D0) then
      !
      BESSK0=1.D30
      return
    endif

    if(X .le. 2.D0) then
      !
      Y=X*X/4.D0
      AX=-log(X/2.D0)*BESSI0(X)
      BESSK0=AX+(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
    else
      !
      Y=(2.D0/X)
      AX=exp(-X)/sqrt(X)
      BESSK0=AX*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
    endif

    return

  end function BESSK0

  !*****************************************************************************************************************

  pure function BESSK1(X)

    ! "K" Bessel function of rank 1

    implicit none

    real(kind=8),intent(in)  :: X
    real(kind=8)             :: BESSK1,Y,AX
    real(kind=8),parameter   :: P1 = 1.d0
    real(kind=8),parameter   :: P2 = 0.15443144d0
    real(kind=8),parameter   :: P3 = -0.67278579d0
    real(kind=8),parameter   :: P4 = -0.18156897d0
    real(kind=8),parameter   :: P5 = -0.1919402d-1
    real(kind=8),parameter   :: P6 = -0.110404d-2
    real(kind=8),parameter   :: P7 = -0.4686d-4
    real(kind=8),parameter   :: Q1 = 1.25331414d0
    real(kind=8),parameter   :: Q2 = 0.23498619d0
    real(kind=8),parameter   :: Q3 = -0.3655620d-1
    real(kind=8),parameter   :: Q4 = 0.1504268d-1
    real(kind=8),parameter   :: Q5 = -0.780353d-2
    real(kind=8),parameter   :: Q6 = 0.325614d-2
    real(kind=8),parameter   :: Q7 = -0.68245d-3

    if(X .eq. 0.D0) then
      !
      BESSK1=1.D32
      return
    endif

    if(X .le. 2.D0) then
      !
      Y=X*X/4.D0
      AX=log(X/2.D0)*BESSI1(X)
      BESSK1=AX+(1.D0/X)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
    else
      !
      Y=(2.D0/X)
      AX=exp(-X)/sqrt(X)
      BESSK1=AX*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
    endif

    return

  end function BESSK1

  !*****************************************************************************************************************

  function slope_limiter(x,x0)

    ! Allows to reduce the derivative function (slope) associated to a given properties
    ! Allows to slowly approach the resolution limit

    implicit none

    real(kind=8),intent(in)  :: x,x0  !
    real(kind=8)             :: slope_limiter

    slope_limiter = 0.d0

    if (x .le. x0) return

    slope_limiter = 1.d0 - exp(-abs(x-x0)/abs(2.d0*x0))

    return
  end function slope_limiter

  !*****************************************************************************************************************

  function Ronbint(f, a, b, param, error, final_order, called_by)

    ! Return the integral of the function f between a and b by using the
    ! Romberg method

    implicit none

    integer(kind=4)                               :: order
    integer(kind=4),parameter                     :: maxorder = 12
    integer(kind=4)                               :: loopk
    integer(kind=4)                               :: loopj
    integer(kind=4),intent(out),optional          :: final_order

    character(*),intent(in),optional              :: called_by    ! name of the calling subroutine

    real(kind=8),intent(in)                       :: a            ! min range integration
    real(kind=8),intent(in)                       :: b            ! max range integration
    real(kind=8)                                  :: Ronbint
    real(kind=8)                                  :: f
    real(kind=8),intent(in),optional,dimension(:) :: param        ! parameters array
    real(kind=8),intent(out),optional             :: error
    real(kind=8)                                  :: e
    real(kind=8)                                  :: dh
    real(kind=8)                                  :: R(0:maxorder, 0:maxorder)

    external                                      :: f

    if (a .gt. b) then
      !
      call IO_print_error_message('Ronbint : a > b',only_rank = rank,called_by = called_by)
      stop
    end if

    order  = 0        ! init
    dh     = b - a    ! init
    e      = 1.d0     ! init
    R(:,:) = 0.d0     ! init
    R(0,0) = dh/2.d0*(f(a,param) + f(b,param))  ! first term
    !
    do while (order .lt. maxorder)
      !
      order = order + 1
      dh = dh/2.d0
      R(order, 0) = 5.d-1*R(order-1, 0)
      !
      do loopk = 1, 2**(order-1)
        !
        R(order, 0) = R(order, 0) + dh*f(a + (2.d0*real(loopk,8)-1.d0)*dh, param)
      end do
      !
      do loopj = 1, order
        !
        R(order, loopj) = R(order, loopj-1) + (R(order, loopj -1) - R(order -1, loopj -1))/(4.d0**(real(loopj,8)) - 1.d0)
      end do
      !
      if (R(order, order) .gt. 0.d0) then
        !
        e = abs(R(order, order) - R(order -1, order -1))/(R(order, order))
        if ((e .lt. num_precision) .and. (order .gt. 6)) exit
      end if
      !
    end do
    !
    if (is_NaN(R(order, order))) then
        !
        call IO_print_warning_message('Ronbint : NAN',only_rank = rank,called_by = called_by)
    end if
    !
    if (R(order, order) .ne. 0.d0) then
      !
      e = e*R(order, order)
    end if
    !
    if (present(final_order)) final_order = order
    !
    if(present(error))              error = e
    !
    Ronbint = R(order, order)

    return
  end function Ronbint

  !*****************************************************************************************************************

  function trap(y,x)

    ! ARRAY TRAPEZIUM RULE TO INTEGRATE PRE-TABULATED Y OVER X

    implicit none

    integer(kind=4)         :: i           ! index loop
    integer(kind=4)         :: n

    real(kind=4),intent(in) :: x(:),y(:)   ! x and y tables
    real(kind=8)            :: trap

    if (size(x) .ne. size(y)) then
      !
      write(*,'(a)') '!!! ERROR in trap : nx /= ny'
      stop ! stop the program
    else
      !
      n = size(x)
    end if

    trap = 0.d0    ! init
    do i = 1,n-1
       !
       trap = trap + (y(i) + y(i+1))*(x(i+1) - x(i))
    end do
    trap = trap*0.5d0

    return

  end function trap

  !*****************************************************************************************************************

  function normal_distribution(mean,sigma)

    ! Allow to give a ramdom variable following a normal distribution of mean "mean" and of standard deviation "sigma"

    implicit none

    real(kind=8),intent(in)    :: mean               ! the mean of the normal distribution
    real(kind=8),intent(in)    :: sigma              ! the standard deviation of the normal distribution
    real(kind=8)               :: normal_distribution
    real(kind=8)               :: r1,r2

    call random_number(r1)   ! get a first uniformed distributed number [0:1]
    call random_number(r2)   ! get a second uniformed distributed number [0:1]

    ! convert to normal distribution

    normal_distribution = sigma*sqrt(-2.0d0*log(r1))*cos(2.0d0*pi*r2) + mean
    return

  end function normal_distribution

  !*****************************************************************************************************************

  function erf_func(x)

    ! erf_func_x -> -1 for x-> -infinity and erf_x -> 1 for x -> +infinity

    implicit none

    real(kind=8),intent(in)   :: x        ! the argument
    real(kind=8)              :: erf_func ! the result
    real(kind=8)              :: xx

    ! Computed with the fitting formula in Lu et al 2010
    xx = x * sqrt(2.0d0)
    if (x > 0.d0) then
       !
       erf_func = 1.0d0 - ( 2.0d0 * exp(-x**2)/( 1.64*xx+sqrt(0.76*xx**2+4.0d0)) )
    else
       !
       erf_func = - 1.0d0 + ( 2.0d0 * exp(-x**2)/(-1.64*xx+sqrt(0.76*xx**2+4.0d0)) )
    endif

    return

  end function erf_func

  !*****************************************************************************************************************

  function locate2D(x,y,tab,nx,ny,xmin,ymin,dx,dy)

    ! RETURN THE VALUE STORED IN A PRE-TABULATED 2D TABLE LINKED TO INPUT PARAMETERS x AND y
    !

    implicit none

    integer(kind=4),intent(in) ::  nx,ny  ! number of data point over x and y axis respectivelly
    integer(kind=4)            ::  ix,iy  ! loop index over x and y

    real(kind=8),intent(in)    :: x,y
    real(kind=8)               :: xv, yv
    real(kind=8),intent(in)    :: tab(nx,ny)
    real(kind=8),intent(in)    :: xmin,ymin   ! minimum values over x and y respectively
    real(kind=8)               :: xmax,ymax   ! maximum values over x and y respectively
    real(kind=8),intent(in)    :: dx,dy       ! stepes over x and y respectively
    real(kind=8)               :: locate2D

    ! x
    !
    xmax = xmin+real(nx-1,8)*dx
    xv   = min(max(x,xmin),xmax)
    !
    ! y
    !
    ymax = ymin+real(ny-1,8)*dy
    yv   = min(max(y,ymin),ymax)
    !
    ix = floor((xv-xmin)/dx) +1   ! Define the x index
    iy = floor((yv-ymin)/dy) +1   ! Define the y index
    locate2D = tab(ix,iy)

    return
  end function locate2D

  !*****************************************************************************************************************

end module utilities
