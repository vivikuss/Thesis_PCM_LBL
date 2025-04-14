module adiabatic_adjustment_module
    use dimensions, only : nLay, nLev

    implicit none

    ! Generalize this
    real(8), parameter :: R_specific = 461  ! Specific gas constant for dry air in J/(kg·K)
    ! real(8), parameter :: Cp = 4000.0       ! Specific heat capacity of dry air at constant pressure in J/(kg·K)
    ! real(8), parameter :: ps = 1e5          ! Reference pressure in Pa
    ! real(8), parameter :: ptop = 1e3        ! Ending pressure level
    ! integer, parameter :: nLev = 101       ! Number of pressure levels

contains

    function potential_temperature(cp_heat, T, P, ps) result(theta)
        real(8), intent(in) :: cp_heat, T, P, ps
        real(8) :: theta

        theta = T * (ps / P) ** (R_specific / cp_heat)
    end function potential_temperature

    subroutine adia_adjustment(cp_heat, T1, T2, ps, P1, P2, dp1, dp2, T1_adj, T2_adj)
        real(8), intent(in) :: cp_heat, T1, T2, P1, P2, dp1, dp2, ps
        real(8), intent(out) :: T1_adj, T2_adj

        real(8) :: theta1, theta2
        real(8) :: theta, theta_mw_init, theta_mw_final

        theta1 = potential_temperature(cp_heat, T1, P1, ps)
        theta2 = potential_temperature(cp_heat, T2, P2, ps)

        theta_mw_init = (theta1*dp1 + theta2*dp2)
        if ((theta2 - theta1) < 1e-6) then
            ! if (theta2 < theta1) then ! potential temp decreasing
            ! Calculate potential temperatures for each layer
            theta = theta_mw_init/(dp1+dp2)
            ! Adjust temperatures to conserve potential temperature
            T1_adj = theta*(P1/ps)**(R_specific/cp_heat)
            T2_adj = theta*(P2/ps)**(R_specific/cp_heat)
        else
            T1_adj = T1
            T2_adj = T2
        end if
        theta_mw_final = potential_temperature(cp_heat, T1_adj, P1, ps)*dp1 + potential_temperature(cp_heat, T2_adj, P2, ps)*dp2

        if (abs((theta_mw_final - theta_mw_init)/theta_mw_final) > 1e-6) then
            print *, "Potential Temp of column not conserved"
        end if
    end subroutine adia_adjustment

    subroutine process_layers(cp_heat, T, plev, play, adjusted_T)
        real(8), intent(in) :: cp_heat, T(:), plev(:), play(:)
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

            call adia_adjustment(cp_heat, T1, T2, ps, play(i), play(i+1), dp1, dp2, adjusted_T(i), adjusted_T(i+1))
        end do
    end subroutine process_layers

    subroutine interpolate(Tlay_new, play, plev, Tlev, Tlev_new)
        real(8), intent(in) :: Tlay_new(:), plev(:), play(:), Tlev(:)
        real(8), intent(out) :: Tlev_new(:)
        real(8) :: gradient
        integer :: i
        ! Calculate the gradient between the first two known data points (Tlay(1) and Tlay(2))
        gradient = (Tlay_new(2) - Tlay_new(1)) / (play(2) - play(1))
    
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

    subroutine convective_adjustment(cp_heat, initial_T, plev, play, adjusted_T, Tlev_new)
        real(8), intent(in) :: cp_heat, initial_T(:), plev(:), play(:)
        real(8), intent(out) :: adjusted_T(:)
        real(8), intent(out) :: Tlev_new(:)  ! Declare Tlev_new with intent(out) here
        real(8) :: old_adjusted_T(size(initial_T))
        real(8), allocatable :: adjusted_curves(:,:) ! Declare as allocatable
        integer :: num_iterations, max_iterations
    
        adjusted_T = initial_T
        num_iterations = 0
        max_iterations = 500  ! Set a maximum number of iterations
    
        allocate(adjusted_curves(size(initial_T), max_iterations))
    
        do while (.true.)
            num_iterations = num_iterations + 1
            old_adjusted_T = adjusted_T
            call process_layers(cp_heat, old_adjusted_T, plev, play, adjusted_T)
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
    
        call interpolate(adjusted_T, play, plev, initial_T, Tlev_new)

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

end module adiabatic_adjustment_module