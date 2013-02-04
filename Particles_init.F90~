!!****if* source/Particles/ParticlesInitialization/Particles_init
!!
!! NAME
!!    Particles_init
!!
!! SYNOPSIS
!!    Particles_init( integer(in) :: myPE,
!!                    integer(in) :: numProcs,
!!                    logical(in) :: restart )
!!
!! DESCRIPTION
!!
!!    General initialization routine for the particle module.
!!
!! ARGUMENTS
!!
!!    myPE     : current processor number
!!    numProcs : total number of processors
!!    restart  : indicates if run is starting from scratch or restarting
!!               from checkpoint file
!!
!! PARAMETERS
!!
!!    useParticles   BOOLEAN [TRUE]  Should particles be used in this simulation?
!!    pt_maxPerProc  INTEGER [100]   Maximum number of particles per processor. Allocates array space
!!                                   Particles are distributed per PROCESSOR rather than per BLOCK
!!    pt_dtFactor    REAL    [0.5]   Factor to make sure that time step is small enough that particles
!!    pt_dtChangeTolerance REAL [0.4] For uncorrected Estimated Midpoint propagation scheme:
!!                                    Do Euler step if change in time step is greater than this
!!                                    percentage.  Set to 0 to always do Euler, set to a huge
!!                                    number to always use estimated midpoint velocities
!!    pt_small       REAL    [1.0E-10] Used for general comparisons of real values 
!!                                   For example, IF (abs(real1 - real2) .lt. pt_small) THEN
!!                                   don't move farther than one block in each step
!!    
!!  NOTE 
!!
!!***

!!#define DEBUG_PARTICLES

subroutine Particles_init (myPE, numProcs, restart)
  
  use Particles_data 
  use Driver_interface, ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
    RuntimeParameters_mapStrToInt
  use Simulation_interface, ONLY : Simulation_mapStrToInt,Simulation_mapParticlesVar
  use pt_interface, ONLY : pt_mapStringParamToInt
  use Particles_interface, ONLY : Particles_specifyMethods
  implicit none
  
#include "constants.h"
#include "Flash.h"
#include "Particles.h"

  integer, INTENT(in) :: numProcs
  integer, INTENT(in) :: myPE
  logical, INTENT(in) :: restart

  integer :: ierr
  character (len=MAX_STRING_LENGTH) :: partAttrPrefix,partAttr
  integer :: i,j
  integer, save :: attributes
  logical :: isLattice, isWithDensity, isByPosition

!-------------------------------------------------------------------------------
  
  !! get out of here if we don't need Particles

  call RuntimeParameters_get ( "useParticles", useParticles)

  !DEV: Neccessary initializations for Particles_sendOutputData 
  !(called from IO_setScalar) so that log message is correct.
  pt_numLocal = 0   ! number of particles on this processor
  pt_myPe = myPE
  pt_numProcs = numProcs
  pt_restart = restart   ! Is this a restart?

     
  if (.NOT. useParticles) return

!We need a work around for pure active particle simulations.
!NONEXISTENT is known to be -1, which does not conflict with 
!real particle types.
!Note the undefine at the end of this subroutine.
#ifdef PASSIVE_PART_TYPE
#define PASSIVE_PART_TYPE_LOCAL PASSIVE_PART_TYPE
#else
#define PASSIVE_PART_TYPE_LOCAL NONEXISTENT
#endif

  
  ! Time stepping parameters, for Particles_computeDt
  call RuntimeParameters_get("pt_dtFactor", pt_dtFactor)
  
  ! Location bounds, needed for Particles_initPosition and Particles_mapParticles2Mesh
  call RuntimeParameters_get( "xmin", pt_xmin)
  call RuntimeParameters_get( "xmax", pt_xmax)
  call RuntimeParameters_get( "ymin", pt_ymin)
  call RuntimeParameters_get( "ymax", pt_ymax)
  call RuntimeParameters_get( "zmin", pt_zmin)
  call RuntimeParameters_get( "zmax", pt_zmax)

    ! space allocation, for Particles_init
  call RuntimeParameters_get ("pt_maxPerProc", pt_maxPerProc)

  ! geometry, for Particles_initPositions
  call RuntimeParameters_get("geometry",pt_str_geometry)
  call RuntimeParameters_mapStrToInt(pt_str_geometry, pt_geometry)

  ! defined variously for the advance methods, used in Particles_initAttributes

  !  Set a small comparison number, for timestep comparisons in some Particles_advance
  call RuntimeParameters_get ("pt_small", pt_small)
  call RuntimeParameters_get ("pt_dtChangeTolerance", pt_dtChangeTolerance)
  call RuntimeParameters_get ("pt_numAtOnce", pt_numAtOnce)
  call RuntimeParameters_get ("pt_logLevel", pt_logLevel)


  ! We create Particles_specifyMethods() at setup time, which is used 
  ! to populate PART_TYPE, PART_INITMETHOD, PART_MAPMETHOD in pt_typeInfo.
  pt_typeInfo=0
  call Particles_specifyMethods()
  pt_typeInfo(PART_LOCAL_MAX,:)=pt_maxPerProc


  !! If any of the particle types are initialized with the lattice method
  !! then the following parameters are needed.
  isLattice=.false.
  isWithDensity=.false.
  isByPosition=.false.

  do i=1,NPART_TYPES
     isLattice=isLattice.or.(pt_typeInfo(PART_INITMETHOD,i)==LATTICE)
     isWithDensity=isWithDensity.or.&
          (pt_typeInfo(PART_INITMETHOD,i)==CELLMASS).or.&
          (pt_typeInfo(PART_INITMETHOD,i)==REJECTION).or.&
          (pt_typeInfo(PART_INITMETHOD,i)==WITH_DENSITY)
     isByPosition=isByPosition.or.(pt_typeInfo(PART_INITMETHOD,i)==CUSTOM)
  end do

  if(isLattice) then
#ifndef PART_INITMETHOD_LATTICE
     call Driver_abortFlash("Particles are set up to use a Lattice initialization method,"// &
          " but a unit to implement that method has not been compiled in!")
#endif
     call RuntimeParameters_get ("pt_numX", pt_numX)
     call RuntimeParameters_get ("pt_numY", pt_numY)
     call RuntimeParameters_get ("pt_numZ", pt_numZ)
     call RuntimeParameters_get ("pt_initialXMin", pt_initialXMin)
     call RuntimeParameters_get ("pt_initialXMax", pt_initialXMax)
     call RuntimeParameters_get ("pt_initialYMin", pt_initialYMin)
     call RuntimeParameters_get ("pt_initialYMax", pt_initialYMax)
     call RuntimeParameters_get ("pt_initialZMin", pt_initialZMin)
     call RuntimeParameters_get ("pt_initialZMax", pt_initialZMax)
  end if

  if(isWithDensity) then
#ifndef PART_INITMETHOD_WITHDENSITY
     call Driver_abortFlash("Particles are set up to use a density-based initialization method,"// &
          " but a unit to implement that method has not been compiled in!")
#endif
     call RuntimeParameters_get ("pt_pRand", pt_pRand)
     call RuntimeParameters_get("pt_numParticlesWanted", pt_numParticlesWanted)
  end if

  if(isByPosition) then
     call RuntimeParameters_get("pt_numParticlesWanted", pt_numParticlesWanted)
  end if 



  ! Guard cell masks, with loads of ifdefs.  This first one is for
  ! use in advancing particles.

  pt_gcMaskForAdvance = .FALSE.
#ifdef DENS_VAR
  pt_gcMaskForAdvance(DENS_VAR) = .TRUE.
#endif

#ifdef VELX_VAR
  pt_gcMaskForAdvance(VELX_VAR) = .TRUE.
#endif

#if NDIM > 1
#ifdef VELY_VAR
  pt_gcMaskForAdvance(VELY_VAR) = .TRUE.
#endif
#endif
  
#if NDIM > 2
#ifdef VELZ_VAR
  pt_gcMaskForAdvance(VELZ_VAR) = .TRUE.
#endif
#endif

#ifdef GPOT_VAR
  pt_gcMaskForAdvance(GPOT_VAR) = .FALSE.
  do i = 1,NPART_TYPES
     pt_gcMaskForAdvance(GPOT_VAR) = pt_gcMaskForAdvance(GPOT_VAR).or.&
          (pt_typeInfo(PART_TYPE,i)/=PASSIVE_PART_TYPE_LOCAL) 
  end do
#endif

  pt_gcMaskForWrite = .FALSE.
  pt_numAttributes=0
  partAttrPrefix='particle_attribute_'

  do i = 1,PT_MAX_ATTRIBUTES

     call pt_mapStringParamToInt(attributes,partAttrPrefix,MAPBLOCK_PART,i)
     if(attributes>0) then
        pt_numAttributes=pt_numAttributes+1
        j=pt_numAttributes
        pt_attributes(j) = attributes
        call Simulation_mapParticlesVar(pt_attributes(j),&
             pt_meshVar(PT_VAR,j), pt_meshVar(PT_MAP,j))

#ifdef DEBUG_PARTICLES
        print*,'  attribute #',pt_numAttributes,' =',pt_attributes(j),'->meshVar',pt_meshVar(PT_VAR,j),pt_meshVar(PT_MAP,j)
#endif
        select case (pt_meshVar(PT_MAP,j))
        case(PARTICLEMAP_UNK)
           pt_gcMaskForWrite(pt_meshVar(PT_VAR,j))=.TRUE.
#if NSCRATCH_GRID_VARS > 0
        case(PARTICLEMAP_SCRATCH)
           ! DO NOTHING
#endif
#if NFACE_VARS > 0
        case(PARTICLEMAP_FACEX)
           pt_gcMaskForWrite(NUNK_VARS+pt_meshVar(PT_VAR,j))=.TRUE.
#if NDIM > 1
        case(PARTICLEMAP_FACEY)
           pt_gcMaskForWrite(NUNK_VARS+NFACE_VARS+pt_meshVar(PT_VAR,j))=.TRUE.
#if NDIM > 2
        case(PARTICLEMAP_FACEZ)
           pt_gcMaskForWrite(NUNK_VARS+NFACE_VARS+NFACE_VARS+pt_meshVar(PT_VAR,j))=.TRUE.
#endif
#endif
#endif
        case default
           if(pt_myPe==MASTER_PE) then
              call concatStringWithInt(partAttrPrefix,i,partAttr) !for diagonstic msg only
              print*,"WARNING: ",trim(partAttr)," does not map to any mesh variable, will be ignored: PARTICLEMAP "
99            format(1x,'Particles_init: particle property',i2,' from variable',i2,', map',i2,', attr #',i2,'.')
              print 99, pt_attributes(j), pt_meshVar(PT_VAR,j), pt_meshVar(PT_MAP,j), j
           end if
           pt_numAttributes = pt_numAttributes - 1 ! ignore this one.
        end select
     else
        if(pt_myPe==MASTER_PE)then
           if (attributes .NE. 0) then
              call concatStringWithInt(partAttrPrefix,i,partAttr) !for diagonstic msg only
              print*,"WARNING: ",trim(partAttr)," not recognized as a particle property, will be ignored"
           end if
        end if
     end if
  end do
#ifdef DEBUG_PARTICLES
  print*,'pt_gcMaskForAdvance:',pt_gcMaskForAdvance
  print*,'pt_gcMaskForWrite:  ',pt_gcMaskForWrite
#endif
  
  
  pt_posInitialized=.false.
  pt_velInitialized=.false.
  if (.not. restart) then !if we starting from scratch
     
     ! allocate data structures for particle storage
     if (.not. allocated(particles)) &
          allocate (particles(NPART_PROPS,pt_maxPerProc), stat=ierr)
     if (ierr /= 0) then
        call Driver_abortFlash("Particles_init:  could not allocate particle array")
     endif
          
     !initialize all particles to NONEXISTENT
     particles = NONEXISTENT
     
  end if  ! end of .not. restart
  
  pt_posAttrib(IAXIS)=POSX_PART_PROP
  pt_posAttrib(JAXIS)=POSY_PART_PROP
  pt_posAttrib(KAXIS)=POSZ_PART_PROP

  pt_posPredAttrib(:) = 1
#ifdef POSPREDX_PART_PROP
  pt_posPredAttrib(IAXIS)=POSPREDX_PART_PROP
  pt_posPredAttrib(JAXIS)=POSPREDY_PART_PROP
  pt_posPredAttrib(KAXIS)=POSPREDZ_PART_PROP
#endif

  pt_velNumAttrib=0
  pt_velAttrib=0
#ifdef VELX_PART_PROP
#ifdef VELX_VAR
  pt_velAttrib(PART_DS_IND,1)=VELX_PART_PROP
  pt_velAttrib(GRID_DS_IND,1)=VELX_VAR
  pt_velNumAttrib = 1

#ifdef VELY_PART_PROP
#ifdef VELY_VAR
  if(NDIM>1) then
     pt_velAttrib(PART_DS_IND,2)=VELY_PART_PROP
     pt_velAttrib(GRID_DS_IND,2)=VELY_VAR
     pt_velNumAttrib = 2
  else
     if (allocated(particles)) particles(VELY_PART_PROP,:)=0.0
  end if
#endif
#endif
  
#ifdef VELZ_PART_PROP
#ifdef VELZ_VAR
  if(NDIM>2) then
     pt_velAttrib(PART_DS_IND,3)=VELZ_PART_PROP
     pt_velAttrib(GRID_DS_IND,3)=VELZ_VAR
     pt_velNumAttrib = NDIM
  else
     if (allocated(particles)) particles(VELZ_PART_PROP,:)=0.0
  end if
#endif
#endif
#else
#ifdef DEBUG_PARTICLES
  if(pt_myPE==MASTER_PE)print*,"no velocity attributes defined for particles"
#endif
#endif
#endif

#ifdef DEBUG_PARTICLES
  if(pt_myPE==MASTER_PE)print*,"Particles_init: pt_velNumAttrib is", pt_velNumAttrib
  if(pt_myPE==MASTER_PE)print*,"Particles_init: pt_velAttrib is", pt_velAttrib
#endif

pt_velPredAttrib(:,:)=1
#ifdef VELPREDX_PART_PROP
#ifdef VELX_VAR
  pt_velPredAttrib(PART_DS_IND,1)=VELPREDX_PART_PROP
  pt_velPredAttrib(GRID_DS_IND,1)=VELX_VAR
  
#ifdef VELPREDY_PART_PROP
#ifdef VELY_VAR
  if(NDIM>1) then
     pt_velPredAttrib(PART_DS_IND,2)=VELPREDY_PART_PROP
     pt_velPredAttrib(GRID_DS_IND,2)=VELY_VAR
  else
     if (allocated(particles)) particles(VELPREDY_PART_PROP,:)=0.0
  end if
#endif
#endif
  
#ifdef VELPREDZ_PART_PROP
#ifdef VELZ_VAR
  if(NDIM>2) then
     pt_velPredAttrib(PART_DS_IND,3)=VELPREDZ_PART_PROP
     pt_velPredAttrib(GRID_DS_IND,3)=VELZ_VAR
  else
     if (allocated(particles)) particles(VELPREDZ_PART_PROP,:)=0.0
  end if
#endif
#endif
#ifdef DEBUG_PARTICLES
  if(pt_myPE==MASTER_PE)print*,"Particles_init: pt_velPredAttrib is", pt_velPredAttrib
#endif
#endif
#endif

#undef PASSIVE_PART_TYPE_LOCAL

  return

end subroutine Particles_init
