!!****if* source/physics/Gravity/GravityMain/Poisson/Multipole/Gravity_init
!!
!! NAME
!!
!!  Gravity_init
!!
!! 
!! SYNOPSIS
!!
!!  Gravity_init(integer(IN) :: myPE)
!!
!!
!! DESCRIPTION
!!
!!  Initialize the multipole Poisson solver.  Read in any of the
!!  runtime parameters for this solver.  All solver common data
!!  is stored in the mpole_common module
!!
!!  ARGUMENTS
!!
!!  myPE - local processor number
!!
!!***

subroutine Gravity_init(myPE)

  use Gravity_data
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY: Grid_getNumProcs
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
    RuntimeParameters_mapStrToInt
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Driver_data, ONLY : dr_restart
  use tree, ONLY : lrefine_max
  use IO_interface, ONLY : IO_getScalar

  implicit none
  real,save :: newton

#include "constants.h"

  integer, intent(IN) :: myPE
  character(len=MAX_STRING_LENGTH) :: strGeometry

  ! Everybody should know these
  grv_myPE = myPE
  call Grid_getNumProcs(grv_numProcs)


  call RuntimeParameters_get("geometry", strGeometry)

  call RuntimeParameters_mapStrToInt(strGeometry, grav_geometry)
  !! DEV testing for invalid gravity geometries?  Perhaps it is done in Grid

  
  call RuntimeParameters_get("grav_boundary_type", grav_boundary_type)

  ! Can't use periodic b.c. with Multipole
  select case (grav_boundary_type)
     case ("periodic")
        call Driver_abortFlash('[Gravity_init] No periodic gravity boundary conditions with Multipole.')
     case ("isolated")
        !! Life is good here
     case default
        call Driver_abortFlash('[Gravity_init] Unsupported gravity boundary conditions, only isolated allowed.')
  end select

  call RuntimeParameters_get("useGravity", useGravity)
  call RuntimeParameters_get("updateGravity", updateGravity)
  call PhysicalConstants_get("Newton", newton)

  if (dr_restart) then
      call IO_getScalar("dynrefinemax", grv_dynRefineMax) 
  else
      grv_dynRefineMax = lrefine_max
  endif

  grav_poisfact = 4. * PI * Newton


  return
end subroutine Gravity_init
