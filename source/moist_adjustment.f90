module moist_adjustment_module
    use dimensions, only : nLay, nLev, nGas
    use composition, only: mu_H2O, mu_He, mu_H2, iGas_H2O
    use fund_consts, only: Rstar
    use thermodynamics, only  : get_psat_H2O, get_Tsat_CO2
    use thermodynamics, only  : therm, tdpsdt, cp_CO2, cp_H2O

    implicit none

    ! Generalize this
    real(8), parameter :: R_specific = 461  ! Specific gas constant for dry air in J/(kg·K)
    real(8) grav
    real(8), parameter :: cp_dry = 5192.6 ! Helium

    namelist/atmos_structure_nml/grav

contains

    ! determine latent heat for given T
    subroutine latent_heat_H2O(T,L)
        real(8), intent(in) :: T
        real(8), intent(out) :: L
        real(8) :: E0, Rvapor, cv, cl, T0

        if(T < 250)then
            L = 0.0
        else
            E0 = 6008/(mu_H2O*1e-3)           !J/kg
            Rvapor = Rstar/(mu_H2O*1e-3)        !J/kg/K
            cv = 4184                          !J/kg/K
            cl = 1900                       !J/kg/K
            T0 = 273.16                        !K

            L = E0 + Rvapor + (cv - cl)*(T - T0)
        endif

    end subroutine latent_heat_H2O

    function potential_temperature(cp_heat, T, P, ps, rv, Lv) result(theta)
        real(8), intent(in) :: cp_heat, T, P, ps, rv, Lv
        real(8) :: theta, Rspec, Rvapor, Rmoist

        Rspec = Rstar/(mu_He*1e-3)
        Rvapor = Rstar/(mu_H2O*1e-3)            !J/kg/K
        Rmoist = Rspec*(1-rv)+Rvapor*rv 

        theta = T*(ps/P)**(Rmoist/cp_heat)*exp((Lv*rv)/(cp_heat*T))     ! moist

    end function potential_temperature

    subroutine adia_adjustment(T1, T2, ps, P1, P2, dp1, dp2, T1_adj, T2_adj)
        real(8), intent(in) :: T1, T2, P1, P2, dp1, dp2, ps
        real(8), intent(out) :: T1_adj, T2_adj

        real(8) :: theta1, theta2
        real(8) :: theta, theta_mw_init, theta_mw_final
        real(8) :: cp_heat1, cp_heat2, rv1, rv2, psat1, psat2, Lv1, Lv2
        real(8) :: Rvapor, Rspec, eps, Rmoist1, Rmoist2

        cp_heat1 = cp_H2O(T1)*1e3
        cp_heat2 = cp_H2O(T2)*1e3
        call get_psat_H2O(T1,psat1)
        call get_psat_H2O(T2,psat2)
        call latent_heat_H2O(T1,Lv1)
        call latent_heat_H2O(T2,Lv2)

        Rvapor = Rstar/(mu_H2O*1e-3)                !J/kg/K
        Rspec = Rstar/(mu_He*1e-3)                     !J/kg/K for He atmosphere
        eps = Rspec/Rvapor

        rv1 = eps*psat1/(P1-psat1) 
        rv2 = eps*psat2/(P2-psat2) 

        Rmoist1 = Rspec*(1-rv1)+Rvapor*rv1 
        Rmoist2 = Rspec*(1-rv2)+Rvapor*rv2 

        theta1 = potential_temperature(cp_heat1, T1, P1, ps, rv1, Lv1) 
        theta2 = potential_temperature(cp_heat1, T2, P2, ps, rv2, Lv2)

        theta_mw_init = (theta1*dp1 + theta2*dp2)
        if ((theta2 - theta1) < 1e-6) then
            ! if (theta2 < theta1) then ! potential temp decreasing
            ! Calculate potential temperatures for each layer
            theta = theta_mw_init/(dp1+dp2)
            ! Adjust temperatures to conserve potential temperature
            T1_adj = theta*(P1/ps)**(Rmoist1/cp_heat1)/exp((Lv1*rv1)/(cp_heat1*T1))
            T2_adj = theta*(P2/ps)**(Rmoist2/cp_heat2)/exp((Lv2*rv2)/(cp_heat2*T2))
        else
            T1_adj = T1
            T2_adj = T2
        end if
        theta_mw_final = potential_temperature(cp_heat1,T1_adj,P1,ps,rv1,Lv1)*dp1+potential_temperature(cp_heat2,T2_adj,P2,ps,rv2,Lv2)*dp2

        if (abs((theta_mw_final - theta_mw_init)/theta_mw_final) > 1e-4) then
            print *, "Potential Temp of column not conserved"
        end if

    end subroutine adia_adjustment

    subroutine process_layers(T, plev, play, adjusted_T)
        real(8), intent(in) :: T(:), plev(:), play(:)
        real(8), intent(out) :: adjusted_T(:)
        integer :: i
        real(8) :: T1, T2, dp1, dp2, ps
    
        ps = plev(1)
        adjusted_T = T
        do i = 1, size(play) - 1
            T1 = adjusted_T(i)
            T2 = adjusted_T(i+1)
            dp1 = plev(i+1) - plev(i)
            dp2 = plev(i+2) - plev(i+1)

            call adia_adjustment(T1, T2, ps, play(i), play(i+1), dp1, dp2, adjusted_T(i), adjusted_T(i+1))
        end do
    end subroutine process_layers

    subroutine interpolate(Tlay_new, play, plev, Tlev, Tlev_new, flay)
        real(8), intent(in) :: Tlay_new(:), plev(:), play(:), Tlev(:)
        real(8), intent(in) :: flay(nLay,nGas)
        real(8), intent(out) :: Tlev_new(:)
        real(8) :: gamma_moist, dry_gradient, gradient
        real(8) :: psat, cp_heat, Lv, rv, psurf_v, ps, cp_dry
        real(8) :: Rvapor, Rspec, eps
        integer :: i

        
        ! Calculate the gradient between the first two known data points (Tlay(1) and Tlay(2))
        dry_gradient = (Tlay_new(2) - Tlay_new(1)) / (play(2) - play(1))

        call get_psat_H2O(Tlay_new(1), psat)
        cp_heat = cp_H2O(Tlay_new(1))*1e3
        call latent_heat_H2O(Tlay_new(1), Lv)

        rv = psat/(plev(1)-psat)
        psurf_v = flay(1,iGas_H2O)*ps

        Rvapor = Rstar/(mu_H2O*1e-3) !J/kg/K
        Rspec = Rstar/(mu_He*1e-3) !J/kg/K for He as background
        eps = Rspec/Rvapor
        
        ! approximate moist lapse rate
        gamma_moist = grav*(1+(Lv*rv)/(Rspec*Tlay_new(1)))/(cp_dry+Lv**2*rv*eps/(Rspec*Tlay_new(1)**2))

        ! dry or moist
        if(psurf_v < psat)then
            gradient = dry_gradient
        else
            gradient = gamma_moist
        end if
       
    
        ! Linear interpolation to get values for Tlev_new corresponding to plev pressure array
        do i = 1, nLev
            if (i == 1) then
                ! Estimate the temperature at the first pressure level assuming a fixed gradient
                Tlev_new(i) = Tlay_new(1) + gradient * (plev(i) - play(1))
            else if (i == nLev) then
                ! Interpolate between Tlay_new(nLev-1) and Tlev(nLev)
                Tlev_new(i) = Tlev_new(i-1) + ((Tlev_new(i-1) - Tlev_new(i-2)) / (plev(i-1) - plev(i-2))) * &
                              (plev(i) - plev(i-1))
            else
                ! Interpolate between Tlay_new(i-1) and Tlay_new(i)
                Tlev_new(i) = Tlay_new(i-1) + (Tlay_new(i) - Tlay_new(i-1)) * &
                              (plev(i) - play(i-1)) / (play(i) - play(i-1))
            end if
        end do
    end subroutine interpolate

    subroutine convective_adjustment(initial_T, plev, play, flay, adjusted_T, Tlev_new)
        real(8), intent(in) :: initial_T(:), plev(:), play(:)
        real(8), intent(in) :: flay(nLay,nGas)
        real(8), intent(out) :: adjusted_T(:)
        real(8), intent(out) :: Tlev_new(:)  
        real(8) :: old_adjusted_T(size(initial_T))
        real(8), allocatable :: adjusted_curves(:,:) 
        integer :: num_iterations, max_iterations
    
        adjusted_T = initial_T
        num_iterations = 0
        max_iterations = 500  ! Set a maximum number of iterations
    
        allocate(adjusted_curves(size(initial_T), max_iterations))
    
        do while (.true.)
            num_iterations = num_iterations + 1
            old_adjusted_T = adjusted_T
            call process_layers(old_adjusted_T, plev, play, adjusted_T)
            if (all(abs(old_adjusted_T - adjusted_T) < 1e-3)) then
                exit
            end if
    
            ! Store the adjusted curve from this iteration
            if (num_iterations <= max_iterations) then
                adjusted_curves(:, num_iterations) = adjusted_T(:)
            else
                print *, "Warning: Exceeded maximum number of iterations"
                exit
            end if
        end do

        ! Print or process adjusted_curves as needed here
    
        call interpolate(adjusted_T, play, plev, initial_T, Tlev_new,flay)

    end subroutine convective_adjustment

    subroutine write_data_to_files(play, Tlay, adjusted)
        real(8), intent(in) :: play(:), Tlay(:), adjusted(:)
        integer :: i

        ! Write play data to file
        open(unit=20, file='play_data.txt', status='replace')
        do i = 1, size(play)
            write(20, *) play(i)
        end do
        close(unit=20)

        ! Write Tlay data to file
        open(unit=40, file='Tlay_data.txt', status='replace')
        do i = 1, size(Tlay)
            write(40, *) Tlay(i)
        end do
        close(unit=40)

        ! Write adjusted data to file
        open(unit=50, file='adjusted_data.txt', status='replace')
        do i = 1, size(adjusted)
            write(50, *) adjusted(i)
        end do
        close(unit=50)

    end subroutine write_data_to_files

end module moist_adjustment_module