!!****f* source/Simulation/Simulation_sendOutputData
!!
!! NAME
!!  Simulation_sendOutputData
!!
!! SYNOPSIS
!! 
!!  Simulation_sendOutputData()
!!  
!! DESCRIPTION 
!!
!! This routine sends the scalar variables owned by the Simulation unit
!! to the IO unit, to be written to a checkpoint file.
!!
!!
!!***

subroutine Simulation_sendOutputData()
    use IO_interface, ONLY : IO_setScalar
    use Simulation_data, ONLY : sim_dynRefineMax
    implicit none

    call IO_setScalar("dynrefinemax", sim_dynRefineMax)
end subroutine Simulation_sendOutputData

