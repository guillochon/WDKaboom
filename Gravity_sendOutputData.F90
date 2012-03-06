!!****f* source/Gravity/Gravity_sendOutputData
!!
!! NAME
!!  Gravity_sendOutputData
!!
!! SYNOPSIS
!! 
!!  Gravity_sendOutputData()
!!  
!! DESCRIPTION 
!!
!! This routine sends the scalar variables owned by the Gravity unit
!! to the IO unit, to be written to a checkpoint file.
!!
!!
!!***

subroutine Gravity_sendOutputData()
    use IO_interface, ONLY : IO_setScalar
    use Gravity_data, ONLY: grv_dynRefineMax
    implicit none

    call IO_setScalar("dynrefinemax", grv_dynRefineMax)
end subroutine Gravity_sendOutputData

