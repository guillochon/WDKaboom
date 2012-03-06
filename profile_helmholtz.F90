subroutine wd_intout(rho0_cgs,T0,wd_mass_tot,cnt,r,rho,m,eint,temper,omega,theta)
!***************************************************************
!                                                              *
!     integrate outward until pressure gets negative,          *
!     give back total mass to correct guess for central        *
!     density; 15.6.98                                         *
!     Version for Helmholtz-EOS; SKR 22.1.2003                 *
!                                                              *
!***************************************************************
      
      use Eos_interface, ONLY : Eos
      use Eos_data, ONLY : eos_smallt
      use Grid_interface, ONLY : Grid_getMinCellSize
      use PhysicalConstants_interface, ONLY : PhysicalConstants_get
      use RuntimeParameters_interface, ONLY : RuntimeParameters_mapStrToInt, RuntimeParameters_get
      use Simulation_data, ONLY : sim_accMass, p_des, s_des, cur_xn, core_xn, torus_xn, &
          sim_nSubZones, sim_virTemp, sim_rotFac, np, ns, sim_profileSubdiv, sim_minProfDelta, sim_maxProfDelta, &
          core_pos, sim_axisRatio
      implicit none
#include "Eos.h"
#include "constants.h"
#include "Flash.h"

      double precision, intent(in) :: T0
      double precision, intent(out) :: wd_mass_tot, r(np), rho(np, ns), m(np), &
          eint(np, ns), temper(np, ns), omega(np, ns), theta(ns)
      integer, intent(out) :: cnt(ns)

      integer i, j, torus_cnt
      double precision rho_cgs_cut
      double precision pr(np, ns),dr,dthick
      double precision rho0_cgs,jspe
      double precision p0,G,mp,kb,minpratio,minrratio
      double precision, dimension(EOS_NUM) :: eosData
      double precision :: s_coretor(ns)
      logical :: at_edge
      logical, dimension(ns) :: terminated
      double precision, dimension(ns) :: pratio, rratio

      call RuntimeParameters_get('sim_rhoAmbient',rho_cgs_cut)
      call PhysicalConstants_get("Newton", G)
      call PhysicalConstants_get("proton mass", mp)
      call PhysicalConstants_get("Boltzmann", kb)
      call Grid_getMinCellSize(dthick)

      eosData(EOS_DENS) = rho0_cgs
      eosData(EOS_TEMP) = T0
      call Eos(MODE_DENS_TEMP,1,eosData,core_xn)

      torus_cnt = 0
      p0 = eosData(EOS_PRES)
      r(1) = 0.d0
      pr = 0.d0
      pr(1,:) = p0
      rho = 0.d0
      rho(1,:) = rho0_cgs
      m(1) = 0.d0
      eint = 0.d0
      eint(1,:) = eosData(EOS_EINT)
      temper = 0.d0
      temper(1,:) = T0
      cnt = 0
      theta(1) = 0.d0
      do j = 2, ns
          theta(j) = theta(j-1) + (PI/2.d0 - 1.0d-4)/(ns - 1)
      enddo
      omega = 0.d0
      !write(*,*) 'central pr, rho, m', pr(0), rho(0), m(0)
!
!--give grid size
!     
      dr = sim_profileSubdiv

      terminated = .false.
      
!     calculate masses and densities at grid points
      do i=2,np
         r(i) = r(i-1) + dr
        
!     mass inside grid point i
         m(i) = m(i-1) + 4.d0*PI*(r(i-1)**2.d0 + r(i-1)*dr + dr**2.d0)*sum(rho(i-1,:))/ns*dr

         if (m(i-1) .le. sim_accMass) then
             cnt = cnt+1
             pr(i,:) = pr(i-1,1) - G*rho(i-1,1)*dr*m(i)/r(i)**2.d0
             if(pr(i,1).lt. 0.d0) exit

             rho(i,:) = rho(i-1,:)
             eosData(EOS_DENS) = rho(i,1)
             eosData(EOS_PRES) = pr(i,1)
             eosData(EOS_TEMP) = T0
             call Eos(MODE_PRES_TEMP,1,eosData,core_xn)
             !print *, i
             !print *, eosData
             rho(i,:) = eosData(EOS_DENS)
             if(rho(i,1).lt.rho_cgs_cut) exit
             temper(i,:) = T0
             eint(i,:) = eosData(EOS_EINT)
         else
             if (torus_cnt .eq. 0) then
                 core_pos = i
                 !jspe = sqrt(G*sim_accMass*r(i))*sim_rotFac
                 jspe = 0.0
                 !if (jspe .ne. 0.0) terminated(ns) = .true.
             endif
             at_edge = .true.
             do j = 1, ns
                 if(terminated(j)) cycle
                 if(pr(i-1,j) .le. 0.d0) then! .or. cos(theta(j)) .lt. sim_rotFac) then
                    terminated(j) = .true.
                    cycle
                 endif


                 !See Eriguchi and Muller 1985
                 !if(jspe .gt. 0.d0) omega(i,j) = (1.d0 + (1.d0/sim_axisRatio))*&
                 !   jspe/(1.d0 + (r(i)*cos(theta(j))/sim_axisRatio/r(core_pos))**2.d0)/r(core_pos)**2.d0
                 if(jspe .gt. 0.d0) omega(i,j) = jspe/(r(i)*cos(theta(j)))**2.d0
                 pr(i,j) = pr(i-1,j) - G*rho(i-1,j)*dr*m(i)/r(i)**2. &
                     - omega(i,j)**2.d0*r(i)*rho(i-1,j)*dr*cos(theta(j))**2.d0
                 if(pr(i,j) .le. 0.d0) then
                    terminated(j) = .true.
                    cycle
                 endif

                 if (torus_cnt .eq. 0) then
                     eosData(EOS_TEMP) = sim_virTemp
                     temper(i,j) = eosData(EOS_TEMP)
                     if (temper(i,j) .le. eos_smallt) then
                         terminated(j) = .true.
                         cycle
                     endif
                     rho(i,j) = rho(i-1,j)
                     eosData(EOS_DENS) = rho(i,j)
                     eosData(EOS_PRES) = pr(i,j)
                     call Eos(MODE_PRES_TEMP,1,eosData,torus_xn)
                     rho(i,j) = eosData(EOS_DENS)
                     if (rho(i,j) .le. rho_cgs_cut) then
                         terminated(j) = .true.
                         cycle
                     endif
                     s_coretor(j) = eosData(EOS_ENTR)
                     eint(i,j) = eosData(EOS_EINT)
                 else
                     if (torus_cnt .eq. 1) then
                         eosData(EOS_DENS) = rho(i-1,j)
                         eosData(EOS_TEMP) = temper(i-1,j)
                     else
                         eosData(EOS_DENS) = 2.0*rho(i-2,j)/(rho(i-2,j) + rho(i-1,j))*rho(i-1,j)
                         eosData(EOS_TEMP) = 2.0*temper(i-2,j)/(temper(i-2,j) + temper(i-1,j))*temper(i-1,j)
                     endif
                     eosData(EOS_ENTR) = s_coretor(j)
                     eosData(EOS_PRES) = pr(i,j)
                     call Eos(MODE_PRES_ENTR,1,eosData,torus_xn)
                     rho(i,j) = eosData(EOS_DENS)
                     temper(i,j) = eosData(EOS_TEMP)
                     if (rho(i,j) .le. rho_cgs_cut) then ! .or. temper(i,j) .le. eos_smallt) then
                         terminated(j) = .true.
                         cycle
                     endif
                     eint(i,j) = eosData(EOS_EINT)
                 endif
                 at_edge = .false.
                 cnt(j) = cnt(j) + 1
                 !print *, eosData
             enddo
             torus_cnt = torus_cnt + 1
             if (at_edge) exit
         endif
         where (terminated .ne. .true.)
             pratio = pr(i,:)/pr(i-1,:)
             rratio = rho(i,:)/rho(i-1,:)
         elsewhere
             pratio = 1.0d0
             rratio = 1.0d0
         endwhere
         minpratio = min(minval(pratio), minval(rratio))
         if (minpratio .lt. sim_maxProfDelta) dr = dr/1.1
         if (minpratio .gt. sim_minProfDelta) dr = 1.1*dr
         !print *, i, dr
      enddo
   
      wd_mass_tot= m(maxval(cnt))

      return
end

subroutine wd_guess_rho_0(mdes,rho0_cgs)

!******************************************************************
!                                                                 *
!     provide guess value for central density of WD of mass mdes  *
!     SKR 22.7.2005                                               *
!                                                                 *
!******************************************************************

      implicit none

      integer ilo,ihi,i
      double precision, intent(in) :: mdes
      double precision, intent(out) :: rho0_cgs
      double precision mgrid(5),lrhogrid(5)     
      double precision dlrhodm,lrho0_cgs

!     mass grid
      mgrid(1)=0.02*2d33
      mgrid(2)=0.2*2d33
      mgrid(3)=0.6*2d33
      mgrid(4)=1.2*2d33
      mgrid(5)=1.35*2d33

!     log-density grid
      lrhogrid(1)=3.36
      lrhogrid(2)=5.38
      lrhogrid(3)=6.61
      lrhogrid(4)=8.17
      lrhogrid(5)=8.98

!     get lower index
      
      if (mdes .lt. mgrid(1)) then
          lrho0_cgs = lrhogrid(1)
      elseif (mdes .gt. mgrid(size(mgrid))) then
          lrho0_cgs = lrhogrid(5)
      else
          do i = 1, size(mgrid)-1
              if (mdes .gt. mgrid(i)) then
                  dlrhodm= (lrhogrid(i+1)-lrhogrid(i))/(mgrid(i+1)-mgrid(i))
                  lrho0_cgs= lrhogrid(i)+ dlrhodm*(mdes-mgrid(i))
              endif
          enddo
      endif
      rho0_cgs= 10**lrho0_cgs

      return
end

