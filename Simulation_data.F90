!!****if* source/Simulation/SimulationMain/Sedov/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data 
!!
!!  DESCRIPTION
!!
!!  Stores the local data for Simulation setup: Sedov
!!  
!! PARAMETERS
!!
!!  sim_pAmbient       Initial ambient pressure
!!  sim_rhoAmbient     Initial ambient density
!!  sim_xctr           Explosion center coordinates
!!  sim_yctr           Explosion center coordinates
!!  sim_zctr           Explosion center coordinates
!!  sim_nsubzones      Number of `sub-zones' in cells for applying 1d profile
!!
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"
#include "Flash.h"

  !! *** Runtime Parameters *** !!

  double precision, save    :: sim_tAmbient, sim_rhoAmbient
  double precision, save    :: sim_gamma, sim_xCenter, sim_yCenter, sim_zCenter
  double precision, save    :: sim_smallX, sim_smallRho, sim_smallP, sim_pi
  double precision, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
  double precision, save    :: sim_tInitial, sim_massAcc, sim_profileSubdiv, sim_minProfDelta, &
                               sim_maxProfDelta, sim_rhoGuess, sim_axisRatio, sim_dbleDetTemp, &
                               sim_detDens, sim_detTemp, sim_detRadius, sim_critDens, sim_critKine
  integer, save             :: sim_nSubZones

  !! *** Variables pertaining to this Simulation *** !!

  integer, parameter        :: np = 200000, ns = 3
  double precision, save    :: sim_inSubZones, sim_inSubzm1
  double precision, save    :: sim_inszd, sim_tRelax, sim_relaxRate, sim_tExplode, sim_tSpinup
  double precision, save    :: sim_accMass, sim_accTemp, wd_radius, sim_torusMass, sim_virTemp, sim_rotFac

  double precision, save    :: radius(np), rhop(np, ns), eint(np, ns), m(np), &
                               temper(np, ns), omega(np, ns), theta(ns)
  double precision, save    :: cur_xn(SPECIES_BEGIN:SPECIES_END)

  double precision, save    :: core_xn(SPECIES_BEGIN:SPECIES_END)
  double precision, save    :: torus_xn(SPECIES_BEGIN:SPECIES_END)
  double precision, save    :: p_des, s_des
  integer, save             :: ipos(ns), core_pos, sim_maxBlocks, sim_detHeight
  logical, save             :: exploded, sim_explodeCore

  integer, save             :: sim_dynRefineMax

end module Simulation_data
