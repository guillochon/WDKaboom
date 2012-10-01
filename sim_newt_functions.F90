MODULE sim_newt_functions

implicit none

CONTAINS

    FUNCTION newt_wd_mass_1d(x)
        use nrtype
        use Simulation_data
        use Logfile_interface
        implicit none
        real, dimension(:), intent(in) :: x
        real, dimension(size(x)) :: newt_wd_mass_1d
        double precision    :: wd_mass_tot
        character(len=32) :: int_to_str
    
        !loop until mass is correct
        wd_mass_tot = 0.0
    
        radius = 0.0e0
        ipos = 0
        rhop = 0.0e0
        call wd_intout(x(1),sim_accTemp,wd_mass_tot,ipos,radius,rhop,m,eint,temper,omega,theta)
    
        call Logfile_stampMessage('WD Mass 1D')
        write (int_to_str, '(e15.9)') radius(ipos(1))
        call Logfile_stampMessage('Radius (equator): ' // trim(adjustl(int_to_str)))
        write (int_to_str, '(e15.9)') wd_mass_tot
        call Logfile_stampMessage('Mass: ' // trim(adjustl(int_to_str)))
        write (int_to_str, '(i10)') ipos(1)
        call Logfile_stampMessage('Profile length (equator): ' // trim(adjustl(int_to_str)))
        write (int_to_str, '(e15.9)') abs(wd_mass_tot-sim_accMass-sim_torusMass)/(sim_accMass + sim_torusMass)
        call Logfile_stampMessage('Mass Error: ' // trim(adjustl(int_to_str)))
        write (int_to_str, '(e15.9)') rhop(1,1)
        call Logfile_stampMessage('Central Density: ' // trim(adjustl(int_to_str)))
        call Logfile_stampMessage('')
    
        newt_wd_mass_1d(1) = (wd_mass_tot - sim_accMass - sim_torusMass)/(sim_accMass + sim_torusMass)
    END FUNCTION newt_wd_mass_1d
    
END MODULE sim_newt_functions
