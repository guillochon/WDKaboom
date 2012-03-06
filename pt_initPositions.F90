!!****if* source/Particles/localAPI/pt_initPositions
!!
!! NAME
!!    pt_initPositions
!!
!! SYNOPSIS
!!
!!    call pt_initPositions(integer(in)  :: blockID,
!!                          logical(out) :: success)
!!
!! DESCRIPTION
!!
!!    Initializes particle locations for one block in the grid.
!!
!! ARGUMENTS
!!
!!  blockID:        local block ID containing particles to create
!!
!!  success:        returns .TRUE. if positions for all particles
!!                  that should be assigned to this block have been
!!                  successfully initialized.
!!
!!***


subroutine pt_initPositions (blockID,success)

  use Particles_data, ONLY:  pt_numLocal, particles, pt_maxPerProc, &
       pt_posAttrib, pt_velNumAttrib, pt_velAttrib, pt_typeInfo, pt_numParticlesWanted
  use Simulation_data, ONLY : sim_xCenter, sim_yCenter, sim_zCenter
  use Grid_interface, ONLY : Grid_getBlkBoundBox, Grid_mapMeshToParticles

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Particles.h"

  integer, INTENT(in) :: blockID
  logical,intent(OUT) :: success

  integer :: part_props=NPART_PROPS, mapType
  integer :: p, i, j, k, minr, maxr, minth, maxth, minph, maxph, nr, nth, nph
  real    :: dr, dth, r, th, ph
  real    :: bxl, bxu, byl, byu, bzl, bzu, xpos, ypos, zpos
  real, dimension(8) :: corner_vals
  real, dimension(2,MDIM) :: boundBox

  nr = idnint(pt_numParticlesWanted**(1.d0/3.d0))
  nph = nr
  nth = nr

  dr = 1.1d10/nr
  dth = 2.0*PI/nth

  p = pt_numLocal

  ! Get grid geometry for this block 

  call Grid_getBlkBoundBox(blockID,boundBox)
  bxl = boundBox(LOW,IAXIS) - sim_xCenter
  bxu = boundBox(HIGH,IAXIS) - sim_xCenter
  
  if (NDIM >= 2) then
     byl = boundBox(LOW,JAXIS) - sim_yCenter
     byu = boundBox(HIGH,JAXIS) - sim_yCenter
  endif
  if (NDIM == 3) then
     bzl = boundBox(LOW,KAXIS) - sim_zCenter
     bzu = boundBox(HIGH,KAXIS) - sim_zCenter
  endif
  
  !Spherical coordinate convention taken from mathworld
  corner_vals(1) = sqrt(bxl**2.d0 + byl**2.d0 + bzl**2.d0)
  corner_vals(2) = sqrt(bxu**2.d0 + byl**2.d0 + bzl**2.d0)
  corner_vals(3) = sqrt(bxl**2.d0 + byu**2.d0 + bzl**2.d0)
  corner_vals(4) = sqrt(bxu**2.d0 + byu**2.d0 + bzl**2.d0)
  corner_vals(5) = sqrt(bxl**2.d0 + byl**2.d0 + bzu**2.d0)
  corner_vals(6) = sqrt(bxu**2.d0 + byl**2.d0 + bzu**2.d0)
  corner_vals(7) = sqrt(bxl**2.d0 + byu**2.d0 + bzu**2.d0)
  corner_vals(8) = sqrt(bxu**2.d0 + byu**2.d0 + bzu**2.d0)

  minr = int(ceiling(minval(corner_vals)/dr))
  maxr = int(floor(maxval(corner_vals)/dr))
  
  corner_vals(1) = acos(bzl/corner_vals(1))
  corner_vals(2) = acos(bzl/corner_vals(2))
  corner_vals(3) = acos(bzl/corner_vals(3))
  corner_vals(4) = acos(bzl/corner_vals(4))
  corner_vals(5) = acos(bzu/corner_vals(5))
  corner_vals(6) = acos(bzu/corner_vals(6))
  corner_vals(7) = acos(bzu/corner_vals(7))
  corner_vals(8) = acos(bzu/corner_vals(8))

  minph = int(ceiling(nph*0.5*minval(cos(corner_vals) + 1.0)))
  maxph = int(floor(nph*0.5*maxval(cos(corner_vals) + 1.0)))

  corner_vals(1) = atan2(byl, bxl)
  corner_vals(2) = atan2(byl, bxu)
  corner_vals(3) = atan2(byu, bxl)
  corner_vals(4) = atan2(byu, bxu)
  corner_vals(5) = atan2(byl, bxl)
  corner_vals(6) = atan2(byl, bxu)
  corner_vals(7) = atan2(byu, bxl)
  corner_vals(8) = atan2(byu, bxu)

  minth = int(ceiling((minval(corner_vals) + PI/2.0)/dth))
  maxth = int(floor((maxval(corner_vals) + PI/2.0)/dth))

  !! initialization in case of lower dimensionality
  zpos = 0.0
  ypos = 0.0
  xpos = 0.0

  do i = minr, maxr
     r = i*dr
     if (r .gt. dr*nr) cycle
     do j = minth, maxth
        th = j*dth
        do k = minph, maxph
           ph = acos(2.0*real(k)/nph - 1.0)

           xpos = r*cos(th)*sin(ph)
           ypos = r*sin(th)*sin(ph)
           zpos = r*cos(ph)
           
           if (xpos .lt. bxl .or. xpos .gt. bxu) cycle
           if (ypos .lt. byl .or. ypos .gt. byu) cycle
           if (zpos .lt. bzl .or. zpos .gt. bzu) cycle

           p = p + 1
           !! Check space allocation
           if (p > pt_maxPerProc) then
              print *,' '
              print *,'PARAMETER pt_maxPerProc is set to ',pt_maxPerProc
              print *,'  To avoid this crash, redimension bigger in your flash.par'
              call Driver_abortFlash &
                   ("pt_initPositionsLattice:  Exceeded max # of particles/processor!")
           endif
           !! particle is defined, set up data structure
           
           particles(BLK_PART_PROP,p) = real(blockID)
#ifdef MASS_PART_PROP
           particles(MASS_PART_PROP,p) = 1.
#endif
           particles(POSX_PART_PROP,p) = xpos
           particles(POSY_PART_PROP,p) = ypos
           particles(POSZ_PART_PROP,p) = zpos
        enddo
     enddo
  enddo

  pt_numLocal = p
  mapType=pt_typeInfo(PART_MAPMETHOD,1)
  call Grid_mapMeshToParticles(particles,&
       part_props,pt_numLocal,&
       pt_posAttrib,pt_velNumAttrib,pt_velAttrib,mapType)

  success=.true.

  return

end subroutine pt_initPositions


