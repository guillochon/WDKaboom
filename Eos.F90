!!****if* source/physics/Eos/EosMain/Helmholtz/Eos
!!
!! NAME
!!
!! Eos
!!
!! SYNOPSIS
!!
!! subroutine Eos(integer(IN) :: mode,
!!                integer(IN) :: vecLen,
!!                real(INOUT) :: eosData(vecLen*EOS_NUM),
!!      optional, real(IN)    :: massFrac(vecLen*NSPECIES),
!!      optional, logical(IN) :: mask(EOS_VARS+1:EOS_NUM)    )
!!
!! DESCRIPTION
!!
!!   Driver for the Helmholtz and Nadyozhin equations of state.
!!   See the NOTES section for important information about this implementation.
!!
!!  This routine applies the equation of state to thermodynamic 
!!  quantities at one or more grid cells.  The number of cells is 
!!  determined by the argument veclen.  Data is packaged for this 
!!  routine in the 1d array, eosData.  The data in eosData is organized: 
!!  1:vecLen points contain the first variable, vecLen+1:2*vecLen points 
!!  contain the second variable, and so on. The number and order of
!!  variables in the array is determined by the constants defined in Eos.h.
!!  
!!  The routine takes different quantities as givens depending on the
!!  value of the mode variable: if mode=MODE_DENS_TEMP, density and
!!  temperature are taken as given, and pressure and energy are generated
!!  as output; if mode=MODE_DENS_EI, density and internal energy are taken as
!!  givens, and pressure and temperature are generated as output.  If
!!  mode=MODE_DENS_PRES, density and pressure are taken as givens, and
!!  energy and temperature are generated as output.
!!  
!!  In addition to pressure, temperature, and internal energy, which are
!!  always thermodynamically consistent after this call, other quantities
!!  such as the various thermodynamic partial derivatives can be
!!  calculated based on the values in the argument, mask.  mask is a
!!  logical array with one entry per quantity, with the order determined
!!  by constants defined in Eos.h (the same as those for the eosData
!!  argument); .true. means return the quantity, .false. means don't.
!!
!!  ARGUMENTS
!! 
!!  mode :    Selects the mode of operation of the Eos unit.
!!             The valid values are MODE_DENS_EI, MODE_DENS_PRES and  
!!             MODE_DENS_TEMP as decribed above.
!!
!!  vecLen   : number of points for each input variable
!!
!!  eosData  : This array is the data structure through which variable values are 
!!             passed in and out of the Eos routine. The arrays is sized as 
!!             EOS_NUM*vecLen. EOS_NUM, and individual input and output
!!             Eos variables are defined in Eos.h. The array is organizes such that
!!             the first 1:vecLen entries represent the first Eos variable, vecLen+1:
!!             2*vecLen represent the second Eos variable and so on. 
!!
!!  massFrac : Contains the mass fractions of the species included in
!!             the simulation. The array is sized as NSPECIES*vecLen.
!!
!!  mask     : Mask is a logical array the size of EOS_DERIVS (number
!!              of partial derivatives that can be computed, defined in
!!              Eos.h), where each index represents a specific partial derivative
!!              that can be calculated by the Eos unit. A .true. value in mask 
!!              results in the corresponding derivative being calculated and 
!!              returned. It should preferably be dimensioned as
!!              mask(EOS_VARS+1:EOS_NUM) in the calling routine 
!!              to exactly match the arguments declaration in Eos Unit.
!!             Note that the indexing of mask does not begin at 1, but rather at one past
!!             the number of variables.
!!
!!             The Helmholtz EOS kernel calculation ignores the mask setting and calculates
!!             all derivatives, whether needed or not.  This routine does not return
!!             derivatives if the mask is requested, but the calculation is not speeded up
!!             by setting the mask.
!!
!!
!! PARAMETERS
!!
!!  eos_tol    Controls the accuracy of the Newton Rhapson iterations for MODE_DENS_EI and 
!!             MODE_DENS_PRES.
!!
!! EXAMPLE
!!
!! --- A single-point at a time example, does not calculate derivatives (based on Cellular Simulation)---
!!
!!  #include "constants.h"   ! for MODE_DENS_TEMP
!!  #include "Flash.h"       ! for NSPECIES
!!  #include "Eos.h"         ! for EOS_VAR order
!!
!!  real  :: temp_zone, rho_zone, ptot, eint, entr, gamma
!!  real, dimension(EOS_NUM)  :: eosData
!!  real, dimension(SPECIES_BEGIN:SPECIES_END) ::  massFraction  
!!  integer, dimension(2,MDIM)                 :: blockRange,blockExtent
!!
!!
!!  massFraction(:) = 1.0e-12        
!!  massFraction(C12_SPEC) = 1.0
!!
!!  .... initiale temp_zone, rho_zone
!!
!!  call Grid_getBlkIndexLimits(blockId,blockRange,blockExtent)
!!  do k = blockRange(LOW,KAXIS), blockRange(HIGH,KAXIS)
!!     do j = blockRange(LOW,JAXIS),blockRange(HIGH,JAXIS)
!!        do i = blockRange(LOW,IAXIS),blockRange(HIGH,IAXIS)
!!
!!           eosData(EOS_TEMP) = temp_zone
!!           eosData(EOS_DENS) = rho_zone
!!
!!           call Eos(MODE_DENS_TEMP,1,eosData,massFraction)
!!
!!           ptot = eosData(EOS_PRES)
!!           eint = eosData(EOS_EINT)  
!!           entr = eosData(EOS_ENTR)
!!           gamma = eosData(EOS_GAMC)
!!           
!!           call Grid_putPointData(blockId,CENTER,TEMP_VAR,EXTERIOR,iPosition,temp_zone)
!!           call Grid_putPointData(blockId,CENTER,DENS_VAR,EXTERIOR,iPosition,rho_zone)
!!           call Grid_putPointData(blockId,CENTER,PRES_VAR,EXTERIOR,iPosition,ptot)
!!           call Grid_putPointData(blockId,CENTER,EINT_VAR,EXTERIOR,iPosition,eint)
!!               if you want ENER_VAR, calculate it from eint and kinetic energy
!!           call Grid_putPointData(blockId,CENTER,GAMC_VAR,EXTERIOR,iPosition,gamma)
!!           call Grid_putPointData(blockId,CENTER,GAME_VAR,EXTERIOR,iPosition,(ptot/(etot*sim_rhoAmbient) + 1.0))
!!
!!         enddo  ! end of k loop
!!     enddo     ! end of j loop
!!  enddo        ! end of i loop
!!
!! ------------------ Row at a time example, with derivates (based on Eos_unitTest) --------
!!
!!  #include "constants.h"   ! for MODE_DENS_TEMP
!!  #include "Flash.h"       ! for NSPECIES, EOS_NUM
!!  #include "Eos.h"         ! for EOS_VAR order
!!  integer veclen, isize, jsize, ksize, i,j,k, e
!!  real, dimension(:), allocatable :: eosData
!!  real, dimension(:), allocatable :: massFrac
!!  logical, dimension (EOS_VARS+1:EOS_NUM) :: mask
!!  real, allocatable, dimension(:,:,:,:) :: derivedVariables
!!  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
!!
!!   ! in the Eos_unitTest, this loops over all blocks.... here is a snippet from inside
!!     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
!!
!!    !  Allocate the necessary arrays for an entire block of data
!!    isize = (blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1)
!!    jsize = (blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1)
!!    ksize = (blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1)
!!    vecLen=isize
!!    allocate(derivedVariables(isize,jsize,ksize,EOS_NUM))
!!    allocate(eosData(vecLen*EOS_NUM))
!!    allocate(massFrac(vecLen*NSPECIES))
!!    mask = .true.
!!
!!    ! indices into the first location for these variables
!!    pres = (EOS_PRES-1)*vecLen
!!    dens = (EOS_DENS-1)*vecLen
!!    temp = (EOS_TEMP-1)*vecLen
!!
!!
!!    call Grid_getBlkPtr(blockID,solnData)
!!    do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
!!        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH, JAXIS)
!!           do i = 1,vecLen
!!              massFrac((i-1)*NSPECIES+1:i*NSPECIES) = &
!!                   solnData(SPECIES_BEGIN:SPECIES_END,ib+i-1,j,k)
!!           end do
!!
!!           eosData(pres+1:pres+vecLen) =  solnData(PRES_VAR,ib:ie,j,k)
!!           eosData(dens+1:dens+vecLen) =  solnData(DENS_VAR,ib:ie,j,k)
!!           ! Eos Helmholtz needs a good initial estimate of temperature no matter what the mode
!!           eosData(temp+1:temp+vecLen) =  solnData(TEMP_VAR,ib:ie,j,k)
!!
!!           call Eos(MODE_DENS_PRES,vecLen,eosData,massFrac,mask)
!!
!!           do e=EOS_VARS+1,EOS_NUM
!!              m = (e-1)*vecLen
!!              derivedVariables(1:vecLen,j-NGUARD,k-NGUARD,e) =  eosData(m+1:m+vecLen)
!!           end do
!!        end do
!!     end do
!!
!! NOTES
!!
!!  NSPECIES is defined in Flash.h.
!!
!!  EOS_VARS and EOS_NUM  are defined in Eos.h.
!!  Calling funtions should included Eos.h, in order to get the definitions of
!!  Eos-specific constants to be able to populate the eosData and mask arrays.
!!  
!!  MODE_DENS_TEMP, MODE_DENS_EI, and MODE_DENS_PRES are defined in constants.h.
!!
!!  All routines calling this routine should include a 
!!  use Eos_interface 
!!  statement, preferable with "ONLY" attribute e.g.
!!  use Eos_interface, ONLY:  Eos
!!
!!  The Helmholtz equation of state calculations are iterative for any mode other
!!  than MODE_DENS_TEMP.  Therefore, the intial estimates for temperature and density
!!  must be pretty good upon entering Eos with any other MODE_....or the calculations will
!!  not converge.
!!
!!  This algorithm uses a data table helm_table.dat which contains the coefficients for
!!  one of the interpolating algorithms.  Upon first entry to the Eos, a binary version of this
!!  table (helm_table.bdat) is created for speed of access.  This binary file should NOT be
!!  carried across machine architectures or even compilers.
!!
!!  When USE_EOS_YE is defined, this routine is replaced by the one in 
!!  physics/Eos/EosMain/Helmholtz/Ye
!!
!!  When operating in MODE_DENS_EI, the INPUT energy is updated.  This change of an input parameter
!!     can be overridden by setting the runtime parameter eos_forceConstantInput to true.
!!     Noted below, see comments prefaced with ConstantInput.
!!  Similarly, when operating in MODE_DENS_PRES, the INPUT pressure is updated.  Physicists need
!!     to be aware of this.  Similarly can be overridden with the runtime parameter/
!!
!!  The accuracy can be adjusted with the parameter eos_tol.
!!
!! SEE ALSO
!! 
!!  Eos.h    defines the variables used.
!!  Eos_wrapped  sets up the data structure.
!!
!!
!!*** 


subroutine Eos(mode,vecLen,eosData,massFrac,mask)

  use Driver_interface, ONLY : Driver_abortFlash
  use Multispecies_interface, ONLY : Multispecies_getSumInv, &
       Multispecies_getSumFrac
  use Grid_interface, ONLY:  Grid_getMyPE
  use Grid_data, ONLY: gr_smallrho
  use Logfile_interface, ONLY:  Logfile_stampMessage

  use eos_helmInterface, ONLY : eos_helm

  use Eos_data, ONLY: eos_smallt, eos_larget, eos_tol, eos_maxNewton, &
      eos_forceConstantInput, nr_eos_dens, nr_eos_eint, nr_eos_pres, nr_eos_entr, nr_eos_xn, &
      nr_eos_mode, nr_eos_k, eos_failed, nr_eos_scale, nr_eos_limited, nr_eos_lim_cnt, eos_smallgamc
  use eos_vecData, ONLY:  tempRow, denRow, etotRow, abarRow, zbarRow, &
      gamcRow, ptotRow, stotRow, dsdRow, dstRow, &
      detRow, dptRow, dpdRow, dedRow, cvRow, cpRow
  use newt_wrappers
  use eos_newt_functions
  use Hydro_data, ONLY: hy_smallp
  use Simulation_data, ONLY: core_xn, torus_xn, sim_tExplode
  use Driver_data, ONLY: dr_simTime

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#ifdef FLASH_MULTISPECIES
#include "Multispecies.h"
#endif

  !     Arguments
  integer, INTENT(in) :: mode, vecLen
  real, INTENT(inout), dimension(vecLen*EOS_NUM) :: eosData
  real, optional,INTENT(in), dimension(vecLen*NSPECIES) :: massFrac
  real, dimension(SPECIES_BEGIN:SPECIES_END) :: locFrac
  ! must correspond to dimensions of Eos_wrapped
  logical,optional, dimension(EOS_VARS+1:EOS_NUM),INTENT(in)::mask

  integer :: i, k, myPE
  integer :: vecBegin,vecEnd
  integer :: pres, temp, dens, gamc, eint
  integer :: abar, zbar
  integer :: entr, dst, dsd
  integer :: dpt, dpd, det, ded, pel, ne, eta, c_v, c_p
  real    :: abarInv, zbarFrac

  ! declare some local storage for the results of the Newton iteration
  real,dimension(1:vecLen) :: ewantRow, pwantRow, swantRow
  !  local storage for forcedConstantInput -- could be allocatable, but might be so slow
  !  that it's not worth the small storage save.
  real,dimension(1:vecLen) :: psaveRow, esaveRow, ssaveRow
  logical  :: trueMask = .false.
  logical  :: firstCall = .true.
  real :: maxch, tmin, tmax, errmax, error
  real, dimension(2) :: input, output

  !      Fill the pipe with the initial temperature, density, and composition.
  !      The composition is parametrized by abar and zbar, which are computed
  !      from the mass fractions xn().

  !  if you are working with electron abundance mass scalars, then you don't
  !  necessarily have to have mass fractions.
  if(.not.present(massFrac)) then
     call Driver_abortFlash("[Eos] Helmholtz needs mass fractions")
  end if

  eos_failed = .false.
  call Grid_getMyPE(myPE)

  ! warn user if mask is used. Helmholtz EOS doesn't use it
  !if (present(mask)) then
  !   if (firstCall) then
  !      call Grid_getMyPE(myPE)
  !      call Logfile_stampMessage(myPE,"[EOS Helmholtz] WARNING!  Mask setting does not speed up Eos Helmholtz calls")
  !      write(*,*)"[EOS Helmholtz] WARNING!  Mask setting does not speed up Eos Helmholtz calls"
  !      firstCall = .false.
  !   end if
  !end if

  vecBegin = 1
  vecEnd = vecLen

  ! These integers are indexes into the lowest location in UNK that contain the appropriate variable
  pres = (EOS_PRES-1)*vecLen
  dens = (EOS_DENS-1)*vecLen
  temp = (EOS_TEMP-1)*vecLen
  eint = (EOS_EINT-1)*vecLen   ! in flash2 eos_helm, this is etot
  gamc = (EOS_GAMC-1)*vecLen   ! in flash2 eos_helm, this is gamc
  abar = (EOS_ABAR-1)*vecLen   
  zbar = (EOS_ZBAR-1)*vecLen   
  entr = (EOS_ENTR-1)*vecLen

  !! For allocatable arrays, set them up now.
#ifndef FIXEDBLOCKSIZE
  call eos_vecAlloc(vecLen)
#endif

  do k = 1, vecLen

     tempRow(k)    = eosData(temp+k)
     denRow(k)     = eosData(dens+k)
     ewantRow(k)   = eosData(eint+k)   ! store desired internal energy for mode=2 case
     pwantRow(k)   = eosData(pres+k)   ! store desired pressure for mode=3 case
     swantRow(k)   = eosData(entr+k)   ! store desired pressure for mode=3 case
#ifdef FLASH_MULTISPECIES
     !Calculate the inverse in a way that allows for zero mass fractions
     locFrac = massFrac((k-1)*NSPECIES+1:k*NSPECIES)

     if (dr_simTime .lt. sim_tExplode) then
         if (locFrac(HE4_SPEC) .gt. locFrac(C12_SPEC) + locFrac(O16_SPEC)) then
             locFrac = torus_xn
         else
             locFrac = core_xn
         endif
     endif

     call Multispecies_getSumInv(A, abarInv, locFrac)
     abarRow(k) = 1.e0 / abarInv

     call Multispecies_getSumFrac(Z, zbarFrac, locFrac)
     zbarRow(k) = abarRow(k) * zbarFrac
#else
     ! No multispecies defined, use default values (same as Gamma formulation)
     abarRow(k) = eos_singleSpeciesA
     zbarRow(k) = eos_singleSpeciesZ
#endif

  enddo

  eosData(abar+1:abar+vecLen) = abarRow(1:vecLen) 
  eosData(zbar+1:zbar+vecLen) = zbarRow(1:vecLen)

  !! Save the input parameters if eos_forceConstantInput is defined
  if (eos_forceConstantInput) then
     esaveRow = ewantRow
     psaveRow = pwantRow
     ssaveRow = swantRow
  end if


  !==============================================================================

  !      MODE_DENS_TEMP  temperature and density given

  !      Crank the EOS on the pipes filled above, then fill the FLASH arrays
  !      with the thermodynamic quantities returned by the EOS.


  if (mode==MODE_DENS_TEMP) then

     call eos_helm(1,vecLen)

     eosData(pres+1:pres+vecLen)=ptotRow(1:vecLen)
     eosData(eint+1:eint+vecLen)=etotRow(1:vecLen)
     eosData(gamc+1:gamc+vecLen)=gamcRow(1:vecLen)
     eosData(entr+1:entr+vecLen)=stotRow(1:vecLen)

     !==============================================================================
     !      MODE_DENS_EI  internal energy and density given

  else if (mode==MODE_DENS_EI) then
     nr_eos_mode = mode

     ! Do the first eos call with all the zones in the pipe
     !  NOTE that eos_helm can ONLY operate in the equivalent of
     !  MODE_DENS_TEMP, as it returns pressure, energy and derivatives only
     !  So if you send in a crappy temperature here, you'll get a crappy starting
     !  position and the iteration won't converge.
     !  Initial temperature here is what is stored in the grid, even though we 
     !    SUSPECT this is not in equilibrium (or we wouldn't be calling Eos if it was fine)

     call eos_helm(vecBegin,vecEnd)
     !  Now eos_helm has returned ptotRow, etotRow, detRow, and gamcRow

     ! Loop over the zones individually now
     do k = vecBegin, vecEnd
        input = 1.0d0
        nr_eos_scale(1) = tempRow(k)
        nr_eos_eint = ewantRow(k)
        nr_eos_k = k
        nr_eos_lim_cnt = 0
        call run_annewt(newt_eos, 1, eos_tol, input(1), output(1), eos_failed)
        if (eos_failed .and. .not. nr_eos_limited) then
            input(1) = output(1)
            call run_newt(newt_eos, 1, eos_tol, input(1), output(1), error)
            if (error .gt. eos_tol) then
                print *, 'An EOS call failed both newton routines. Final error was:', error
            endif
        endif
        tempRow(k) = max(nr_eos_scale(1)*output(1), eos_smallt)
        if (tempRow(k) .ge. eos_larget) tempRow(k) = eos_smallt
     end do

     ! Crank through the entire eos one last time
     call eos_helm(vecBegin,vecEnd)

     if (ptotRow(k) .lt. hy_smallp) then
        ptotRow(k) = hy_smallp
        gamcRow(k) = eos_smallgamc
     endif

     ! Fill the FLASH arrays with the results.  

     !  In MODE_DENS_EI, we should be generating temperature and pressure (plus gamma and entropy)
     eosData(temp+1:temp+vecLen)=tempRow(1:vecLen)
     eosData(pres+1:pres+vecLen)=ptotRow(1:vecLen)
     eosData(gamc+1:gamc+vecLen)=gamcRow(1:vecLen)
     eosData(entr+1:entr+vecLen)=stotRow(1:vecLen)

     !  Update the energy to be the true energy, instead of the energy we were trying to meet
     !  ConstantInput LBR and KW believe this is WRONG -- the input arrays should not be changed
     if (eos_forceConstantInput)  then
        eosData(eint+1:eint+vecLen) = esaveRow(1:vecLen)
     else
        eosData(eint+1:eint+vecLen) = etotRow(1:vecLen)
     end if


     !==============================================================================

     !      MODE_DENS_PRES  pressure and density given

  else if (mode==MODE_DENS_PRES) then
     nr_eos_mode = mode

     call eos_helm(vecBegin,vecEnd)
     !  Now eos_helm has returned ptotRow, etotRow, detRow, and gamcRow

     ! Loop over the zones individually now
     do k = vecBegin, vecEnd
        input = 1.0d0
        nr_eos_scale(1) = tempRow(k)
        nr_eos_pres = pwantRow(k)
        nr_eos_k = k
        call run_annewt(newt_eos, 1, eos_tol, input(1), output(1), eos_failed)
        tempRow(k) = min(max(nr_eos_scale(1)*output(1), eos_smallt), eos_larget)
     end do

     ! Crank through the entire eos one last time
     call eos_helm(vecBegin,vecEnd)

     ! Fill the FLASH arrays with the results.  
     eosData(temp+1:temp+vecLen)=tempRow(1:vecLen)
     eosData(gamc+1:gamc+vecLen)=gamcRow(1:vecLen)
     eosData(eint+1:eint+vecLen)=etotRow(1:vecLen)
     eosData(entr+1:entr+vecLen)=stotRow(1:vecLen)

     ! Update the pressure to be the equilibrium pressure, instead of the pressure we were trying to meet
     !  ConstantInput LBR and KW believe this is wrong.  See notes at the top of the routine
     if (eos_forceConstantInput) then
        eosData(pres+1:pres+vecLen) = psaveRow(1:vecLen)
     else
        eosData(pres+1:pres+vecLen) = ptotRow(1:vecLen)
     end if
  else if (mode==MODE_PRES_TEMP) then
     nr_eos_mode = mode

     call eos_helm(vecBegin,vecEnd)
     !  Now eos_helm has returned ptotRow, etotRow, detRow, and gamcRow

     ! Loop over the zones individually now
     do k = vecBegin, vecEnd
        input = 1.0d0
        nr_eos_scale(1) = denRow(k)
        nr_eos_pres = pwantRow(k)
        nr_eos_k = k
        call run_annewt(newt_eos, 1, eos_tol, input(1), output(1), eos_failed)
        denRow(k) = max(nr_eos_scale(1)*output(1), gr_smallrho)
     end do

     ! Crank through the entire eos one last time
     call eos_helm(vecBegin,vecEnd)

     ! Fill the FLASH arrays with the results.  
     eosData(dens+1:dens+vecLen)=denRow(1:vecLen)
     eosData(gamc+1:gamc+vecLen)=gamcRow(1:vecLen)
     eosData(eint+1:eint+vecLen)=etotRow(1:vecLen)
     eosData(entr+1:entr+vecLen)=stotRow(1:vecLen)

     ! Update the pressure to be the equilibrium pressure, instead of the pressure we were trying to meet
     !  ConstantInput LBR and KW believe this is wrong.  See notes at the top of the routine
     if (eos_forceConstantInput) then
        eosData(pres+1:pres+vecLen) = psaveRow(1:vecLen)
     else
        eosData(pres+1:pres+vecLen) = ptotRow(1:vecLen)
     end if
  else if (mode==MODE_PRES_ENTR) then
     nr_eos_mode = mode

     call eos_helm(vecBegin,vecEnd)
     !  Now eos_helm has returned ptotRow, etotRow, detRow, and gamcRow

     ! Loop over the zones individually now
     do k = vecBegin, vecEnd
        input = 1.0d0
        nr_eos_scale(1) = denRow(k)
        nr_eos_scale(2) = tempRow(k)
        nr_eos_pres = pwantRow(k)
        nr_eos_entr = swantRow(k)
        nr_eos_k = k
        call run_annewt(newt_eos, 2, eos_tol, input, output, eos_failed)
        !if (output(1) .ne. output(1) .or. output(2) .ne. output(2)) then
        !    print *, output
        !    call Driver_abortFlash('outputs NaN!')
        !endif
        denRow(k) = max(nr_eos_scale(1)*output(1), gr_smallrho)
        tempRow(k) = max(nr_eos_scale(2)*output(2), eos_smallt)
        if (tempRow(k) .ge. eos_larget) tempRow(k) = eos_smallt
     end do

     ! Crank through the entire eos one last time
     call eos_helm(vecBegin,vecEnd)

     ! Fill the FLASH arrays with the results.  
     eosData(dens+1:dens+vecLen)=denRow(1:vecLen)
     eosData(temp+1:temp+vecLen)=tempRow(1:vecLen)
     eosData(gamc+1:gamc+vecLen)=gamcRow(1:vecLen)
     eosData(eint+1:eint+vecLen)=etotRow(1:vecLen)
     eosData(entr+1:entr+vecLen)=stotRow(1:vecLen)

     ! Update the pressure to be the equilibrium pressure, instead of the pressure we were trying to meet
     !  ConstantInput LBR and KW believe this is wrong.  See notes at the top of the routine
     if (eos_forceConstantInput) then
        eosData(pres+1:pres+vecLen) = psaveRow(1:vecLen)
        eosData(entr+1:entr+vecLen) = ssaveRow(1:vecLen)
     else
        eosData(pres+1:pres+vecLen) = ptotRow(1:vecLen)
        eosData(entr+1:entr+vecLen) = stotRow(1:vecLen)
     end if
  else
     call Driver_abortFlash('[Eos] Error: unknown input mode in routine Eos')
  end if


  ! Get the optional values
  if(present(mask)) then
     ! Entropy derivatives
     if(mask(EOS_DST)) then
        dst = (EOS_DST-1)*vecLen
        eosData(dst+1:dst+vecLen) = dstRow(1:vecLen)
     end if
     if(mask(EOS_DSD)) then
        dsd = (EOS_DSD-1)*vecLen
        eosData(dsd+1:dsd+vecLen) = dsdRow(1:vecLen)
     end if
     if(mask(EOS_DPT)) then
        dpt = (EOS_DPT-1)*vecLen
        eosData(dpt+1:dpt+vecLen) = dptRow(1:vecLen)
     end if
     if(mask(EOS_DPD)) then
        dpd = (EOS_DPD-1)*vecLen
        eosData(dpd+1:dpd+vecLen) = dpdRow(1:vecLen)
     end if
     if(mask(EOS_DET))then
        det = (EOS_DET-1)*vecLen
        eosData(det+1:det+vecLen) = detRow(1:vecLen)
     end if
     if(mask(EOS_DED))then 
        ded = (EOS_DED-1)*vecLen
        eosData(ded+1:ded+vecLen) = dedRow(1:vecLen)
     end if

     if(mask(EOS_CV))then
        if(mask(EOS_DET)) then
           c_v = (EOS_CV-1)*vecLen
           eosData(c_v+1:c_v+vecLen) = cvRow(1:vecLen)
        else
           call Driver_abortFlash("[Eos] cannot calculate C_V without DET.  Set mask appropriately.")
        end if
     end if

     if(mask(EOS_CP))then
        if(mask(EOS_CV).and.mask(EOS_DET)) then
           c_p = (EOS_CP-1)*vecLen
           eosData(c_p+1:c_p+vecLen) = cpRow(1:vecLen)
        else
           call Driver_abortFlash("[Eos] cannot calculate C_P without C_V and DET.  Set mask appropriately.")
        end if
     end if
  end if




  !! Close up arrays if previously allocated
#ifndef FIXEDBLOCKSIZE  
  call eos_vecDealloc()
#endif

  return

end subroutine Eos


