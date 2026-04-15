!=============================================================================
!  simv0_spatial_1d.f90  --  1D Panel Thermal PDE  (Method of Lines)
!
!  Solves the transient 1D heat equation along the panel depth (z-direction)
!  using Method of Lines (MOL) discretisation into N_layers cell-centred slabs.
!
!  ---- PDE formulation -------------------------------------------------------
!
!  For z in [0, L_panel]  (z=0: top surface, z=L_panel: bottom/back):
!
!    rho * cp * dT/dt  =  d/dz( kappa_cond * dT/dz )  +  S(z,t)
!
!  Source term S(z,t) = kappa_abs * G_eff(t) * exp(-kappa_abs * z)  [W/m^3]
!    This is the volumetric solar heat generation from Beer-Lambert absorption.
!
!  Boundary conditions:
!    z = 0 (top face):
!      -kappa_cond * dT/dz = h_eff*(T_1-T_amb) + eps*sigma*(T_1^4 - T_sky^4)
!                           [convective + radiative loss at the top surface]
!    z = L_panel (bottom face):
!       kappa_cond * dT/dz = 0    [insulated back-sheet]
!
!  ---- MOL discretisation (finite-volume, cell-centred) ----------------------
!
!  N_layers cells; cell i centred at z_i = (i - 0.5)*dz, i = 1..N_layers
!  dz = L_panel / N_layers = 0.002 m per layer
!
!  Cell 1 (top, adjacent to atmosphere):
!    rho*cp*dz * dT_1/dt =
!        kappa_cond/dz * (T_2 - T_1)         <- conduction from cell 2
!      - h_eff*(T_1 - T_amb)                  <- convective BC at top face
!      - eps*sigma*(T_1^4 - T_sky^4)          <- radiative BC at top face
!      - G_up_corr                            <- two-way BTE upwelling loss
!      + Q_solar(1)                            <- solar absorbed in cell 1 [W/m^2]
!
!  Cell i, 1 < i < N_layers:
!    rho*cp*dz * dT_i/dt =
!        kappa_cond/dz * (T_{i-1} - 2*T_i + T_{i+1})
!      + Q_solar(i)
!
!  Cell N_layers (bottom, insulated):
!    rho*cp*dz * dT_N/dt =
!        kappa_cond/dz * (T_{N-1} - T_N)     <- only one conduction flux
!      + Q_solar(N)
!
!  Wind speed state (appended as element N_layers+1):
!    dWS/dt = -lambda_BL*(WS - WS_geo) + gamma_BL*(T_1 - T_amb)  [NS BL ODE]
!
!  ---- State vector layout ---------------------------------------------------
!
!    y(1)           : T_1   (top layer temperature, K)
!    y(2..N_layers) : T_2..T_N
!    y(N_layers+1)  : WS    (wind speed, m/s)
!
!  ---- Subroutines provided --------------------------------------------------
!
!    spatial_rhs(t, y, n, dydt)  –  matches rhs_interface in simv0_ode_suite
!    init_spatial_y(y_out)       –  initialise state vector at t=t_start
!
!  Compile (after simv0_config, simv0_bte_coupling):
!    gfortran -O2 -c simv0_spatial_1d.f90
!=============================================================================
module simv0_spatial_1d
    use simv0_config
    use simv0_bte_coupling
    implicit none

contains

    !--------------------------------------------------------------------------
    !  spatial_rhs  –  RHS of the spatial 1D MOL system
    !
    !  Matches the abstract interface rhs_interface in simv0_ode_suite so it
    !  can be passed directly to euler_step / heun_step / rk4_step / rk45_step.
    !
    !  Expected n = N_layers + 1
    !--------------------------------------------------------------------------
    subroutine spatial_rhs(t, y, n, dydt)
        real(8), intent(in)  :: t
        integer, intent(in)  :: n
        real(8), intent(in)  :: y(n)
        real(8), intent(out) :: dydt(n)

        !  Extract state variables
        real(8) :: T_lay(N_layers)   ! layer temperatures (K)
        real(8) :: WS                ! wind speed (m/s)

        !  Physics quantities
        real(8) :: G_surf            ! solar irradiance at panel top (W/m^2)
        real(8) :: Q_solar(N_layers) ! solar absorbed per layer (W/m^2)
        real(8) :: G_up_corr         ! two-way upwelling correction (W/m^2)
        real(8) :: h_eff             ! mixed convection coefficient (W/m^2/K)

        !  Layer thermal constants (per unit area, per unit depth)
        !    rcp_dz = rho*cp*dz [J/m^2/K]; divide source [W/m^2] to get [K/s]
        !    Numerically: 750 * 900 * 0.002 = 1350 J/m^2/K
        real(8) :: rcp_dz

        integer :: i

        !  Compute rcp_dz (cannot be a parameter; depends on module variables)
        rcp_dz = rho_pan * cp_pan * dz_layer

        !  Unpack state vector
        T_lay(1:N_layers) = y(1:N_layers)
        WS                = y(N_layers+1)   ! unclamped; h_eff_mixed handles this

        !  ----  Solar source term  (Beer-Lambert through the panel)  ----------
        G_surf = geff(t)
        call bte_solar_per_layer(G_surf, Q_solar)

        !  ----  Two-way BTE upwelling correction  ----------------------------
        !  G_up_corr is the upwelling thermal flux (W/m^2) escaping the top
        !  surface from interior emission.  It reduces the effective energy
        !  available to the top layer.
        G_up_corr = upwelling_top(T_lay)

        !  ----  Mixed convection at top surface  -----------------------------
        h_eff = h_eff_mixed(T_lay(1), WS)

        !  ====  Layer energy balances (all terms already in K/s) =============

        !  --- Cell 1: top layer (exposed to atmosphere) ----------------------
        !
        !  dT_1/dt = [kappa_cond/dz * (T_2 - T_1)             ] (conduction inward)
        !            / (rho*cp*dz)
        !          - h_eff * (T_1 - T_amb)  / (rho*cp*dz)    (convective loss)
        !          - eps*sigma*(T_1^4 - T_sky^4) / (rho*cp*dz)  (radiative loss)
        !          + Q_solar(1) / (rho*cp*dz)                 (solar gain)
        !
        !  Note: G_up_corr is the two-way BTE diagnostic (upwelling from interior
        !  layers escaping through the top surface). It is not added here as a
        !  separate loss term because each layer's thermal emission is already
        !  captured by the Stefan-Boltzmann top-surface BC (which uses T_1) plus
        !  heat conduction from interior layers to T_1. Adding G_up_corr here
        !  would double-count the interior emission already routed to T_1 via
        !  conduction. G_up_corr is written to output for diagnostic purposes.
        dydt(1) = (  (kappa_cond / dz_layer) * (T_lay(2) - T_lay(1))                &
                   - h_eff      * (T_lay(1) - T_amb)                                 &
                   - eps_panel  * SIGMA * (T_lay(1)**4 - T_sky**4)                   &
                   + Q_solar(1) * (alpha_panel - eta_T_func(T_lay(1)))               &
                 ) / rcp_dz

        !  --- Cells 2 to N_layers-1: interior layers -------------------------
        do i = 2, N_layers - 1
            dydt(i) = (  (kappa_cond / dz_layer) * (T_lay(i-1) - 2.0d0*T_lay(i) + T_lay(i+1)) &
                       + Q_solar(i) * alpha_panel                                               &
                     ) / rcp_dz
        end do

        !  --- Cell N_layers: bottom layer (insulated back-sheet) --------------
        !
        !  dT_N/dt = kappa_cond/dz * (T_{N-1} - T_N) / (rho*cp*dz) + Q_solar(N)/(rho*cp*dz)
        !  (no flux through the insulated bottom face)
        dydt(N_layers) = (  (kappa_cond / dz_layer) * (T_lay(N_layers-1) - T_lay(N_layers)) &
                          + Q_solar(N_layers) * alpha_panel                                   &
                        ) / rcp_dz

        !  --- Wind speed ODE: NS boundary-layer dynamics ---------------------
        !  dWS/dt = -lambda_BL*(WS - WS_geo) + gamma_BL*(T_1 - T_amb)
        dydt(N_layers+1) = -lambda_BL * (WS - WS_geo) + gamma_BL * (T_lay(1) - T_amb)

    end subroutine spatial_rhs

    !--------------------------------------------------------------------------
    !  init_spatial_y  –  Initialise the spatial state vector at t = t_start
    !
    !  All layers start at T_0 (= T_amb); wind speed starts at WS_0.
    !--------------------------------------------------------------------------
    subroutine init_spatial_y(y_out)
        real(8), intent(out) :: y_out(N_layers+1)

        integer :: i

        do i = 1, N_layers
            y_out(i) = T_0
        end do
        y_out(N_layers+1) = WS_0

    end subroutine init_spatial_y

end module simv0_spatial_1d
