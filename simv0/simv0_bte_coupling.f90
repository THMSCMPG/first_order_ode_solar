!=============================================================================
!  simv0_bte_coupling.f90  --  Two-Way BTE / NS Physics Functions
!
!  Implements the full physics toolkit used by every RHS in the driver:
!
!  BTE (Beer-Lambert / Boltzmann photon transport):
!    geff          – diurnal effective surface irradiance  (W/m^2)
!    downwelling   – Beer-Lambert downwelling flux profile through panel (W/m^2)
!    upwelling_top – two-way upwelling flux escaping the top surface (W/m^2)
!    bte_solar_per_layer – net solar energy deposited in each layer (W/m^2)
!
!  NS (Navier-Stokes / convection):
!    h_eff_mixed   – Churchill-Usagi mixed convection coefficient  (W/m^2/K)
!
!  Thermal-electrical coupling:
!    eta_T_func    – Skoplaki-Palyvos temperature-dependent efficiency
!
!  ---- Two-way BTE coupling model ----------------------------------------
!
!  One-way (standard Beer-Lambert):
!    G_down(z) = G_eff * exp(-kappa_abs * z)
!    Solar absorbed in layer i: dG_i = G_down(z_{i-1}) - G_down(z_i)  [W/m^2]
!
!  Two-way correction (backward BTE, new in simv0):
!    Each interior layer emits upward as a grey body.
!    Upwelling flux at the top surface from interior layers 1..N:
!
!      G_up_top = eps * sigma * sum_{j=1}^{N} T_j^4
!                 * (1 - exp(-kappa_abs*dz))    <- emittance of layer j
!                 * exp(-kappa_abs*(j-1)*dz)    <- transmittance to top face
!
!    This upwelling is subtracted from the net energy balance of the top
!    layer, closing the irradiance-temperature feedback loop:
!      Higher interior T → larger G_up_top → reduced net top-layer heating.
!
!  Compile (after simv0_config):  gfortran -O2 -c simv0_bte_coupling.f90
!=============================================================================
module simv0_bte_coupling
    use simv0_config
    implicit none

contains

    !--------------------------------------------------------------------------
    !  geff  –  Beer-Lambert diurnal effective surface irradiance (W/m^2)
    !
    !  G_eff(t) = G_toa * cos_z(t) * exp(-tau_tot / cos_z(t))
    !    cos_z(t) = sin(pi * t / T_day)  [solar elevation angle proxy]
    !  Returns 0 outside the photoperiod [0, T_day].
    !--------------------------------------------------------------------------
    real(8) function geff(t)
        real(8), intent(in) :: t

        real(8) :: cos_z

        if (t <= 0.0d0 .or. t >= T_day) then
            geff = 0.0d0
            return
        end if

        cos_z = sin(PI * t / T_day)           ! solar elevation proxy

        if (cos_z < 1.0d-6) then
            geff = 0.0d0
        else
            geff = G_toa * cos_z * exp(-tau_tot / cos_z)
        end if

    end function geff

    !--------------------------------------------------------------------------
    !  h_eff_mixed  –  Churchill-Usagi mixed convection coefficient (W/m^2/K)
    !
    !  McAdams forced:   h_f = h_forced_0 + h_forced_1 * max(WS, 0)
    !  Churchill natural: h_n = C_nc * |T - T_amb|^(1/3)
    !  Mixed (n=3):      h   = (h_f^3 + h_n^3)^(1/3)
    !
    !  WS is clamped to >= 0 before computing h_forced to prevent the
    !  cube-root of a negative argument (NaN) when the NS boundary-layer
    !  ODE drives WS below zero at night.  The WS dynamics themselves are
    !  NOT clamped so the BL ODE converges to its physical night-time
    !  equilibrium WS_eq = WS_geo + (gamma_BL/lambda_BL)*(T-T_amb) < 0.
    !--------------------------------------------------------------------------
    real(8) function h_eff_mixed(T, WS)
        real(8), intent(in) :: T, WS

        real(8) :: h_f, h_n, WS_safe

        WS_safe = max(0.0d0, WS)    ! clamp: h_forced requires non-negative WS
        h_f = h_forced_0 + h_forced_1 * WS_safe
        h_n = C_nc * abs(T - T_amb)**(1.0d0/3.0d0)
        h_eff_mixed = (h_f**3 + h_n**3)**(1.0d0/3.0d0)

    end function h_eff_mixed

    !--------------------------------------------------------------------------
    !  eta_T_func  –  Skoplaki-Palyvos temperature-dependent efficiency
    !
    !  eta(T) = eta_stc * [1 - beta_T * (T - T_stc)]   clamped to [0, 1]
    !--------------------------------------------------------------------------
    real(8) function eta_T_func(T)
        real(8), intent(in) :: T

        eta_T_func = eta_stc * (1.0d0 - beta_T * (T - T_stc))
        eta_T_func = max(0.0d0, min(1.0d0, eta_T_func))

    end function eta_T_func

    !--------------------------------------------------------------------------
    !  downwelling  –  Beer-Lambert flux at each layer BOUNDARY  (W/m^2)
    !
    !  G_down(0) = G_surf    (incident at top surface)
    !  G_down(i) = G_surf * exp(-kappa_abs * i * dz_layer)   i = 1..N_layers
    !
    !  Output array G_down has indices 0:N_layers (N_layers+1 elements).
    !--------------------------------------------------------------------------
    subroutine downwelling(G_surf, G_down)
        real(8), intent(in)  :: G_surf
        real(8), intent(out) :: G_down(0:N_layers)

        integer :: i

        G_down(0) = G_surf
        do i = 1, N_layers
            G_down(i) = G_surf * exp(-kappa_abs * real(i, 8) * dz_layer)
        end do

    end subroutine downwelling

    !--------------------------------------------------------------------------
    !  upwelling_top  –  Two-way backward BTE: upwelling flux leaving the top
    !                    surface due to thermal emission from each layer  (W/m^2)
    !
    !  Each layer j (j=1 at top) acts as a grey body with:
    !    emittance     = (1 - exp(-kappa_abs * dz_layer))
    !    transmittance = exp(-kappa_abs * (j-1) * dz_layer)  [from layer j to top]
    !
    !  G_up_top = eps_panel * sigma * sum_{j=1}^{N_layers} T_j^4
    !             * (1 - exp(-kappa_abs*dz_layer))
    !             * exp(-kappa_abs*(j-1)*dz_layer)
    !--------------------------------------------------------------------------
    real(8) function upwelling_top(T_layers)
        real(8), intent(in) :: T_layers(N_layers)

        real(8) :: emittance, contrib
        integer :: j

        emittance = 1.0d0 - exp(-kappa_abs * dz_layer)

        upwelling_top = 0.0d0
        do j = 1, N_layers
            contrib = eps_panel * SIGMA * T_layers(j)**4 &
                    * emittance                           &
                    * exp(-kappa_abs * real(j-1, 8) * dz_layer)
            upwelling_top = upwelling_top + contrib
        end do

    end function upwelling_top

    !--------------------------------------------------------------------------
    !  bte_solar_per_layer  –  Solar energy absorbed in each layer  (W/m^2)
    !
    !  Uses Beer-Lambert downwelling.  Layer i (1-indexed from top) absorbs:
    !    Q_solar(i) = G_down(i-1) - G_down(i)
    !               = G_surf * exp(-kappa_abs*(i-1)*dz_layer)
    !                         * (1 - exp(-kappa_abs*dz_layer))
    !
    !  This is a one-way solar absorption profile.  The two-way upwelling
    !  correction is applied to the top layer energy balance in the driver.
    !--------------------------------------------------------------------------
    subroutine bte_solar_per_layer(G_surf, Q_solar)
        real(8), intent(in)  :: G_surf
        real(8), intent(out) :: Q_solar(N_layers)

        real(8) :: G_down(0:N_layers)
        integer :: i

        call downwelling(G_surf, G_down)

        do i = 1, N_layers
            Q_solar(i) = G_down(i-1) - G_down(i)
        end do

    end subroutine bte_solar_per_layer

end module simv0_bte_coupling
