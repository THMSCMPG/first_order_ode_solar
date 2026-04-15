!=============================================================================
!  simv0_config.f90  --  Configuration Module for simv0 Integration Platform
!
!  Central repository of all physical constants, panel parameters, BTE/NS
!  coefficients, spatial model dimensions, and solver control knobs used by
!  every other module in simv0.  Import with:
!
!    use simv0_config
!
!  Compile: gfortran -O2 -c simv0_config.f90
!=============================================================================
module simv0_config
    implicit none

    !--------------------------------------------------------------------------
    !  Mathematical / physical constants
    !--------------------------------------------------------------------------
    real(8), parameter :: PI    = 3.141592653589793d0
    real(8), parameter :: SIGMA = 5.670374419d-8      ! W/m^2/K^4  (CODATA 2018)

    !--------------------------------------------------------------------------
    !  Panel thermal / optical parameters  (from rk4_ode / bte_ns_ode baseline)
    !--------------------------------------------------------------------------
    real(8), parameter :: alpha_panel = 0.9d0         ! shortwave absorptivity
    real(8), parameter :: A_panel     = 1.6d0         ! panel area (m^2)
    real(8), parameter :: m_pan       = 12.0d0        ! panel mass (kg)
    real(8), parameter :: cp_pan      = 900.0d0       ! specific heat (J/kg/K)
    real(8), parameter :: eps_panel   = 0.85d0        ! longwave emissivity
    real(8), parameter :: T_amb       = 298.15d0      ! ambient temperature (K)
    real(8), parameter :: T_sky       = 278.15d0      ! effective sky temperature (K)

    !--------------------------------------------------------------------------
    !  Thermal-electrical coupling  (Skoplaki-Palyvos, from coupled_thermal_elec)
    !--------------------------------------------------------------------------
    real(8), parameter :: eta_stc  = 0.18d0           ! STC conversion efficiency
    real(8), parameter :: beta_T   = 0.0045d0         ! temperature coefficient /K  (c-Si)
    real(8), parameter :: T_stc    = 298.15d0         ! STC reference temperature (K)

    !--------------------------------------------------------------------------
    !  BTE atmospheric parameters  (from bte_ns_ode.f90)
    !--------------------------------------------------------------------------
    real(8), parameter :: G_toa    = 1361.0d0         ! top-of-atm solar irradiance (W/m^2)
    real(8), parameter :: tau_ray  = 0.09d0            ! Rayleigh scattering optical depth
    real(8), parameter :: tau_wv   = 0.03d0            ! water vapour absorption
    real(8), parameter :: AOD      = 0.20d0            ! aerosol optical depth (moderate)
    real(8), parameter :: tau_tot  = tau_ray + tau_wv + AOD  ! total one-way optical depth = 0.32
    real(8), parameter :: T_day    = 43200.0d0         ! photoperiod = 12 hours (s)

    !--------------------------------------------------------------------------
    !  NS boundary-layer parameters  (from bte_ns_ode.f90)
    !--------------------------------------------------------------------------
    real(8), parameter :: h_forced_0 = 5.7d0           ! McAdams intercept (W/m^2/K)
    real(8), parameter :: h_forced_1 = 3.8d0           ! McAdams slope (W/m^2/K per m/s)
    real(8), parameter :: C_nc       = 1.31d0          ! Churchill natural-conv constant
    real(8), parameter :: WS_geo     = 3.0d0           ! geostrophic wind speed (m/s)
    real(8), parameter :: lambda_BL  = 1.0d0/600.0d0  ! BL relaxation rate (1/s)
    real(8), parameter :: gamma_BL   = 3.0d-3          ! thermal buoyancy coupling (m/s/K)

    !--------------------------------------------------------------------------
    !  Spatial 1D panel model  (new in simv0)
    !
    !  Panel cross-section divided into N_layers cell-centred slabs along z:
    !    Layer 1 (top)    : anti-reflective coating + glass (exposed to atmosphere)
    !    Layer 2          : EVA encapsulant (front)
    !    Layer 3 (centre) : silicon cell  (primary absorber)
    !    Layer 4          : EVA encapsulant (rear)
    !    Layer 5 (bottom) : back sheet (insulated outer face)
    !
    !  Thermal conductivity uses a volume-averaged glass/Si/polymer composite.
    !  Optical absorption coefficient: 100 /m → 63% solar absorbed over 1 cm
    !    (exp(-100 * 0.01) = exp(-1) ≈ 0.37 transmitted; ~63% absorbed).
    !--------------------------------------------------------------------------
    integer, parameter :: N_layers   = 5              ! number of panel layers
    real(8), parameter :: L_panel    = 0.01d0         ! total panel thickness (m)  = 1 cm
    real(8), parameter :: kappa_cond = 1.0d0          ! thermal conductivity (W/m/K)
    real(8), parameter :: kappa_abs  = 100.0d0        ! optical absorption coeff (1/m)

    !  Derived spatial quantities (computed from above; kept as parameters)
    real(8), parameter :: dz_layer  = L_panel / real(N_layers, 8)  ! layer thickness (m) = 0.002 m
    real(8), parameter :: rho_pan   = m_pan / (A_panel * L_panel)  ! effective density kg/m^3
    !   rho_pan = 12.0 / (1.6 * 0.01) = 750 kg/m^3
    real(8), parameter :: alpha_diff = kappa_cond / (rho_pan * cp_pan)  ! thermal diffusivity m^2/s
    !   alpha_diff = 1.0 / (750 * 900) = 1.481e-6 m^2/s

    !--------------------------------------------------------------------------
    !  Integration / simulation parameters
    !--------------------------------------------------------------------------
    real(8), parameter :: t_start   = 0.0d0           ! sunrise (s)
    real(8), parameter :: t_end     = 86400.0d0        ! 24-hour diurnal cycle (s)
    real(8), parameter :: T_0       = 298.15d0         ! initial panel temperature (K)
    real(8), parameter :: WS_0      = 1.5d0            ! initial wind speed (m/s)
    real(8), parameter :: h_default = 100.0d0          ! default fixed step size (s)

    !--------------------------------------------------------------------------
    !  RK45 Dormand-Prince adaptive control
    !--------------------------------------------------------------------------
    real(8), parameter :: atol        = 0.1d0           ! absolute tolerance (K)
    real(8), parameter :: rtol        = 1.0d-6          ! relative tolerance
    real(8), parameter :: rk45_safety = 0.9d0           ! conservative step-size factor
    real(8), parameter :: h_min       = 1.0d0           ! minimum step (s)
    real(8), parameter :: h_max       = 600.0d0         ! maximum step (s)
    integer, parameter :: max_steps   = 50000           ! safety limit

    !--------------------------------------------------------------------------
    !  Simulation mode codes (used by driver and method selector)
    !--------------------------------------------------------------------------
    integer, parameter :: MODE_DECOUPLED = 1   ! 1-state T; const G, h, eta  (baseline)
    integer, parameter :: MODE_COUPLED   = 2   ! 2-state [T,WS]; BTE+NS+eta(T)
    integer, parameter :: MODE_SPATIAL   = 3   ! (N_layers+1)-state; 1D PDE + two-way BTE

    !--------------------------------------------------------------------------
    !  Integration method codes (returned by method selector)
    !--------------------------------------------------------------------------
    integer, parameter :: METHOD_EULER = 1
    integer, parameter :: METHOD_HEUN  = 2
    integer, parameter :: METHOD_RK4   = 3
    integer, parameter :: METHOD_RK45  = 4

end module simv0_config
