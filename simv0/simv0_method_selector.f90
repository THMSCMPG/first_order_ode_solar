!=============================================================================
!  simv0_method_selector.f90  --  Intelligent ODE Integration Method Selector
!
!  Automatically chooses the most appropriate numerical integration method
!  for the solar-panel ODE based on:
!    1. Simulation mode (decoupled / coupled BTE-NS / spatial 1D PDE)
!    2. Thermal stiffness ratio  tau_thermal / h_step
!    3. User-specified accuracy tolerance
!
!  Selection rules (evaluated in priority order, highest score wins):
!    Rule 1: Spatial 1D PDE (MODE_SPATIAL)       → RK45  (adaptive error control)
!    Rule 2: Very high accuracy (tol < 0.01 K)   → RK45  (minimum local error)
!    Rule 3: Tight stiffness  (tau/h < 5)        → RK4   (stability + accuracy)
!    Rule 4: Coupled BTE-NS system (MODE_COUPLED) → RK4  (nonlinear coupling)
!    Rule 5: Moderate accuracy (tol < 0.1 K)     → RK4
!    Rule 6: Any coupling present                 → Heun  (minimum for stability)
!    Rule 7: Default (decoupled, low accuracy)    → Euler
!
!  Thermal time constant:
!    tau = m*cp / (h_conv_ref * A)   where h_conv_ref uses McAdams at WS_geo
!    tau ≈ 12*900 / ((5.7 + 3.8*3.0)*1.6) ≈ 419 s  (characteristic panel time)
!
!  Compile (after simv0_config):  gfortran -O2 -c simv0_method_selector.f90
!=============================================================================
module simv0_method_selector
    use simv0_config
    implicit none

contains

    !--------------------------------------------------------------------------
    !  select_method
    !
    !  Inputs:
    !    sim_mode    integer  MODE_DECOUPLED | MODE_COUPLED | MODE_SPATIAL
    !    h_step      real(8)  proposed (or default) step size in seconds
    !    tol_in      real(8)  desired absolute accuracy tolerance in Kelvin
    !
    !  Outputs:
    !    method_code integer  METHOD_EULER | METHOD_HEUN | METHOD_RK4 | METHOD_RK45
    !    method_name char(32) human-readable method name
    !    reason      char(80) brief rationale string (for diagnostic output)
    !--------------------------------------------------------------------------
    subroutine select_method(sim_mode, h_step, tol_in, method_code, method_name, reason)
        integer,       intent(in)  :: sim_mode
        real(8),       intent(in)  :: h_step
        real(8),       intent(in)  :: tol_in
        integer,       intent(out) :: method_code
        character(32), intent(out) :: method_name
        character(80), intent(out) :: reason

        real(8) :: h_conv_ref, tau_thermal, tau_ratio
        integer :: score

        !  Reference convection coefficient: McAdams at geostrophic wind speed
        h_conv_ref = h_forced_0 + h_forced_1 * WS_geo   ! ≈ 5.7 + 3.8*3.0 = 17.1 W/m^2/K

        !  Thermal time constant: tau = m*cp / (h_conv * A)
        tau_thermal = m_pan * cp_pan / (h_conv_ref * A_panel)
        !  tau = 12*900 / (17.1*1.6) ≈ 395 s

        !  Stiffness ratio: tau / h_step  (large = well-resolved = less stiff numerically)
        tau_ratio = tau_thermal / h_step

        !  -------  scoring accumulator  (take max across all rules)  ----------
        score = METHOD_EULER    ! start at lowest method

        !  Rule 1 – Spatial 1D PDE: always use RK45 (error control over many states)
        if (sim_mode == MODE_SPATIAL) then
            score = max(score, METHOD_RK45)
        end if

        !  Rule 2 – Very high accuracy (tol < 0.01 K): adaptive RK45
        if (tol_in < 0.01d0) then
            score = max(score, METHOD_RK45)
        end if

        !  Rule 3 – Tight stiffness (tau/h < 5): RK4 for accuracy
        if (tau_ratio < 5.0d0) then
            score = max(score, METHOD_RK4)
        end if

        !  Rule 4 – Coupled BTE-NS (2-state [T,WS]): RK4 for nonlinear coupling
        if (sim_mode == MODE_COUPLED) then
            score = max(score, METHOD_RK4)
        end if

        !  Rule 5 – Moderate accuracy (tol < 0.1 K): at least RK4
        if (tol_in < 0.1d0) then
            score = max(score, METHOD_RK4)
        end if

        !  Rule 6 – Any coupling present (not pure decoupled): at least Heun
        if (sim_mode /= MODE_DECOUPLED) then
            score = max(score, METHOD_HEUN)
        end if

        !  Rule 7 – Default: Euler (decoupled, low accuracy, tau/h >> 1)
        !           (score already initialised to METHOD_EULER)

        method_code = score

        !  -------  human-readable output  ------------------------------------
        select case (method_code)
            case (METHOD_EULER)
                method_name = 'Forward Euler'
                write(reason, '(A,F6.1,A,F7.3,A)') &
                    'Decoupled baseline, tau/h=', tau_ratio, &
                    ', tol=', tol_in, ' K: speed > accuracy'

            case (METHOD_HEUN)
                method_name = 'Improved Euler (Heun)'
                write(reason, '(A,F6.1,A,F7.3,A)') &
                    'Light coupling, tau/h=', tau_ratio, &
                    ', tol=', tol_in, ' K: O(h^2) sufficient'

            case (METHOD_RK4)
                method_name = 'Classical RK4'
                if (sim_mode == MODE_COUPLED) then
                    write(reason, '(A,F6.1,A,F7.3,A)') &
                        'BTE-NS coupling, tau/h=', tau_ratio, &
                        ', tol=', tol_in, ' K: O(h^4) accuracy'
                else if (tau_ratio < 5.0d0) then
                    write(reason, '(A,F6.1,A,F7.3,A)') &
                        'Stiffness: tau/h=', tau_ratio, &
                        ' < 5, tol=', tol_in, ' K: O(h^4) accuracy'
                else
                    write(reason, '(A,F7.3,A)') &
                        'High accuracy tol=', tol_in, ' K: O(h^4) needed'
                end if

            case (METHOD_RK45)
                method_name = 'RK45 Dormand-Prince (adaptive)'
                if (sim_mode == MODE_SPATIAL) then
                    write(reason, '(A,F6.1,A)') &
                        'Spatial 1D PDE, tau/h=', tau_ratio, &
                        ': adaptive error control essential'
                else
                    write(reason, '(A,F7.3,A)') &
                        'High-accuracy tol=', tol_in, &
                        ' K: adaptive step-size control'
                end if
        end select

    end subroutine select_method

end module simv0_method_selector
