module tools_mod

    !*****************************************************************************************************************
    !
    !  OVERVIEW
    !
    ! The tools module defines some usefull transversal methods
    !   locate2D, intergation, etc ...
    !
    !*****************************************************************************************************************

    use parameters  ! Acces to global defintions and properties

    implicit none

    public
    
contains

    !
    ! SUBROUTINES
    !
    ! **********************************

    !
    ! FUNCTIONS
    !
    ! **********************************

    function locate2D(x, y, tab, nx, ny, xmin, ymin, dx, dy) result(val)

        !
        ! Return for parameters (x ,y) the specific value stored in a pre-tabultaed 2D table

        implicit none

        integer(kind=ikd), intent(in) ::  nx, ny      ! Dimensions (x, y)
        integer(kind=ikd)             ::  ix, iy      ! Loop index over x and y

        real(kind=rkd), intent(in)    :: x, y         ! Coordinates
        real(kind=rkd)                :: xv, yv
        real(kind=rkd),intent(in)     :: tab(nx, ny)  ! Input 2D array
        real(kind=rkd),intent(in)     :: xmin, ymin   ! Minimum values over x and y directions respectively
        real(kind=rkd)                :: xmax, ymax   ! Maximum values over x and y directions respectively
        real(kind=rkd),intent(in)     :: dx, dy       ! elementary step over x and y directions respectively
        real(kind=rkd)                :: val

        ! x
        !
        xmax = xmin + real(nx - 1, kind=rkd) * dx
        xv = min(max(x, xmin), xmax)
        !
        ! y
        !
        ymax = ymin + real(ny - 1, kind=rkd) * dy
        yv   = min(max(y, ymin), ymax)
        !
        ix = floor((xv - xmin) / dx) + 1   ! Define the x index
        iy = floor((yv - ymin) / dy) + 1   ! Define the y index
        !
        val = tab(ix,iy)

        return
    end function locate2D


end module tools_mod