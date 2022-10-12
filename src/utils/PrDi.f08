module PrDi_mod

    implicit none

    public

    !*****************************************************************************************************************
    !
    ! OVERVIEW
    !
    ! Here are defined global parameter used for process distribution.
    ! G.A.S. uses a set processes that allow to compute in parallel
    ! physical and luminous propreties of all modeled galaxies.
    ! The G.A.S model is fully compatible with MPI cluster
    ! the standard scheme uses N+1 processes
    !*****************************************************************************************************************

    ! DEFINITION OF GLOBAL VARIABLE LINKED TO PROCESS DISTRIBUTION
    !
    ! MPI process variable, rank, number of processes
    integer(kind=4)                      :: myRank             ! Rank of the process
    integer(kind=4)                      :: mainProcessRank    ! rank of the main process
    integer(kind=4)                      :: nbProc             ! Number of processus available for the computation

contains

end module PrDi_mod
