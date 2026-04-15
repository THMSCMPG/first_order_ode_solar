!=============================================================================
!  simv0_driver.f90  --  Main Driver: simv0 Integration Platform
!
!  Orchestrates intelligent method selection, three simulation modes, and the
!  full-physics spatial + two-way BTE solver.  Runs all three modes in a
!  single execution and writes comparison output files.
!
!  ---- Variable naming note --------------------------------------------------
!  Fortran is case-insensitive, so T_xxx and t_xxx are the SAME symbol.
!  Time variables use 'x' prefix (xd, xc, xs) to avoid conflicts with
!  temperature variables from simv0_config (T_amb, T_sky, T_0, etc.).
!
!  ---- Simulation modes (defined in simv0_config) ---------------------------
!
!  MODE_DECOUPLED (1):
!    State: y = [T]
!    Physics: constant G=1000 W/m^2, h_conv=15 W/m^2/K, eta=0.18 (baseline)
!    Method: Forward Euler  (simple problem; speed over accuracy)
!
!  MODE_COUPLED (2):
!    State: y = [T, WS]
!    Physics: BTE Beer-Lambert G_eff(t), NS mixed h_eff(T,WS), eta(T)
!    Matches bte_ns_ode.f90 for direct validation.
!    Method: Classical RK4 (nonlinear two-state coupling)
!
!  MODE_SPATIAL (3):
!    State: y = [T_1..T_5, WS]  (5 layers + wind speed)
!    Physics: 1D MOL heat equation + two-way BTE upwelling correction
!    Method: RK45 Dormand-Prince with adaptive step-size control
!
!  ---- Output files ----------------------------------------------------------
!    simv0_decoupled.dat  – decoupled baseline (t | T | P)
!    simv0_coupled.dat    – coupled BTE-NS     (t | T | WS | G_eff | h_eff | eta | P)
!    simv0_spatial.dat    – spatial 1D         (t | T_1..T_5 | WS | G_up | P)
!    simv0_diagnostic.dat – method selection log + RK45 step-size history
!    simv0_comparison.dat – side-by-side T_dec | T_coup | T_surf for validation
!
!  Compile & link (see Makefile):  make
!  Run:   ./simv0
!  Plot:  gnuplot plots/plot_simv0.gp
!=============================================================================
program simv0_driver
    use simv0_config
    use simv0_method_selector
    use simv0_ode_suite
    use simv0_bte_coupling
    use simv0_spatial_1d
    implicit none

    !--------------------------------------------------------------------------
    !  State vector sizes
    !--------------------------------------------------------------------------
    integer, parameter :: NDEC = 1             ! decoupled: 1 state  [T]
    integer, parameter :: NCOU = 2             ! coupled:   2 states [T, WS]
    integer, parameter :: NSPA = N_layers + 1  ! spatial:   N+1 states [T_i, WS]

    !--------------------------------------------------------------------------
    !  State vectors (y prefix; new = next-step candidate)
    !--------------------------------------------------------------------------
    real(8) :: yd(NDEC), ydn(NDEC)    ! decoupled
    real(8) :: yc(NCOU), ycn(NCOU)    ! coupled
    real(8) :: ys(NSPA), ysn(NSPA)    ! spatial

    !--------------------------------------------------------------------------
    !  Time variables  (x prefix avoids case-conflict with T_ module params)
    !--------------------------------------------------------------------------
    real(8) :: xd, xc, xs, xcmp

    !--------------------------------------------------------------------------
    !  Step-size variables (h prefix)
    !--------------------------------------------------------------------------
    real(8) :: hd, hc, hs, errn, hnew

    !--------------------------------------------------------------------------
    !  Integer counters
    !--------------------------------------------------------------------------
    integer :: nfixed, nstep, nok, nbad

    !--------------------------------------------------------------------------
    !  Physics diagnostics written to output
    !--------------------------------------------------------------------------
    real(8) :: Gnow, hcnow, etanow
    real(8) :: Pout_d, Pout_c, Pout_s

    !--------------------------------------------------------------------------
    !  Statistics accumulators  (Tp = peak temperature, Ey = energy yield)
    !--------------------------------------------------------------------------
    real(8) :: Tp_d, Tp_c, Tp_s
    real(8) :: Ey_d, Ey_c, Ey_s

    !--------------------------------------------------------------------------
    !  Spatial-mode diagnostics (Gup = upwelling; Hs = final step size)
    !--------------------------------------------------------------------------
    real(8) :: Gup_s, Hs_final

    !--------------------------------------------------------------------------
    !  Method selection outputs
    !--------------------------------------------------------------------------
    integer        :: mc_d, mc_c, mc_s
    character(32)  :: mn_d, mn_c, mn_s
    character(80)  :: mr_d, mr_c, mr_s

    !==========================================================================
    !  1.  METHOD SELECTION
    !==========================================================================
    call select_method(MODE_DECOUPLED, h_default, atol, mc_d, mn_d, mr_d)
    call select_method(MODE_COUPLED,   h_default, atol, mc_c, mn_c, mr_c)
    call select_method(MODE_SPATIAL,   h_default, atol, mc_s, mn_s, mr_s)

    !--------------------------------------------------------------------------
    !  Print audit to stdout
    !--------------------------------------------------------------------------
    write(*,'(A)') repeat('=', 72)
    write(*,'(A)') '  simv0  --  Intelligent ODE Integration Platform'
    write(*,'(A)') '            First-Order ODE Solar Panel Simulation'
    write(*,'(A)') repeat('=', 72)
    write(*,'(A)')
    write(*,'(A)') '  Method selection audit:'
    write(*,'(A)') '  ' // repeat('-', 68)
    write(*,'(A,A)') '  DECOUPLED  -> ', trim(mn_d)
    write(*,'(A,A)') '    Reason:     ', trim(mr_d)
    write(*,'(A,A)') '  COUPLED    -> ', trim(mn_c)
    write(*,'(A,A)') '    Reason:     ', trim(mr_c)
    write(*,'(A,A)') '  SPATIAL    -> ', trim(mn_s)
    write(*,'(A,A)') '    Reason:     ', trim(mr_s)
    write(*,'(A)')
    write(*,'(A,F6.0,A)') '  Simulation duration: ', t_end/3600.0d0, ' hours'
    write(*,'(A,F6.0,A)') '  Fixed step size:     ', h_default, ' s'
    write(*,'(A,I4,A)')   '  Spatial layers:      ', N_layers, ' cells'
    write(*,'(A)') repeat('=', 72)
    write(*,'(A)')

    !--------------------------------------------------------------------------
    !  Diagnostic file header (RK45 step log appended in spatial loop below)
    !--------------------------------------------------------------------------
    open(unit=70, file='simv0_diagnostic.dat', status='replace', action='write')
    write(70,'(A)') '# simv0 diagnostic: method selection + RK45 step-size log'
    write(70,'(A,A)') '# DECOUPLED method: ', trim(mn_d)
    write(70,'(A,A)') '#   reason: ', trim(mr_d)
    write(70,'(A,A)') '# COUPLED   method: ', trim(mn_c)
    write(70,'(A,A)') '#   reason: ', trim(mr_c)
    write(70,'(A,A)') '# SPATIAL   method: ', trim(mn_s)
    write(70,'(A,A)') '#   reason: ', trim(mr_s)
    write(70,'(A)') '#'
    write(70,'(A)') '# RK45 adaptive log (spatial mode only):'
    write(70,'(A)') '# t(s)          h_used(s)      err_norm       accept'

    nfixed = nint((t_end - t_start) / h_default)

    !==========================================================================
    !  2.  DECOUPLED BASELINE  (1-state; Forward Euler; constant G/h/eta)
    !
    !  The decoupled model uses constant G = 1000 W/m^2, h = 15 W/m^2/K, and
    !  eta = eta_stc. Because this is the simplest possible 1-state ODE and
    !  accuracy is secondary to speed, Forward Euler is the natural choice.
    !
    !  NOTE: The method selector would recommend RK4 here due to the stiffness
    !  criterion (tau/h = 3.9 < 5). However, the "decoupled baseline" in this
    !  driver is deliberately a constant-physics reference where Euler
    !  demonstrates the method at its intended use-case: rough, fast,
    !  exploratory runs. For a time-accurate decoupled BTE run the user
    !  should consult the method selector output and use RK4.
    !
    !  State:  yd(1) = T   (panel temperature, K)
    !==========================================================================
    open(unit=10, file='simv0_decoupled.dat', status='replace', action='write')
    write(10,'(A)') '# simv0_decoupled.dat  [MODE_DECOUPLED / Forward Euler]'
    write(10,'(A)') '# t(s)          T(K)           P_elec(W)'

    yd(1)  = T_0
    xd     = t_start
    hd     = h_default
    Tp_d   = T_0
    Ey_d   = 0.0d0

    do nstep = 0, nfixed
        Pout_d = eta_stc * 1000.0d0 * A_panel        !  constant G = 1000 W/m^2
        write(10,'(3F15.6)') xd, yd(1), Pout_d
        if (yd(1) > Tp_d) Tp_d = yd(1)
        Ey_d = Ey_d + Pout_d * hd

        if (nstep == nfixed) exit

        call euler_step(rhs_decoupled, xd, yd, NDEC, hd, ydn)
        yd = ydn
        xd = xd + hd
    end do
    close(10)

    !==========================================================================
    !  3.  COUPLED BTE-NS  (2-state [T,WS]; Classical RK4)
    !
    !  Matches bte_ns_ode.f90 coupled simulation exactly for validation.
    !  State:  yc(1) = T  (panel temperature, K)
    !          yc(2) = WS (wind speed, m/s)
    !==========================================================================
    open(unit=20, file='simv0_coupled.dat', status='replace', action='write')
    write(20,'(A)') '# simv0_coupled.dat  [MODE_COUPLED / Classical RK4]'
    write(20,'(A)') '# t(s)          T(K)           WS(m/s)        ' // &
                    'G_eff(W/m2)    h_eff(W/m2K)   eta(-)         P_elec(W)'

    yc(1)  = T_0
    yc(2)  = WS_0
    xc     = t_start
    hc     = h_default
    Tp_c   = T_0
    Ey_c   = 0.0d0

    do nstep = 0, nfixed
        Gnow   = geff(xc)
        hcnow  = h_eff_mixed(yc(1), yc(2))
        etanow = eta_T_func(yc(1))
        Pout_c = etanow * Gnow * A_panel
        write(20,'(7F15.6)') xc, yc(1), yc(2), Gnow, hcnow, etanow, Pout_c
        if (yc(1) > Tp_c) Tp_c = yc(1)
        Ey_c = Ey_c + Pout_c * hc

        if (nstep == nfixed) exit

        call rk4_step(rhs_coupled, xc, yc, NCOU, hc, ycn)
        yc = ycn
        xc = xc + hc
    end do
    close(20)

    !==========================================================================
    !  4.  SPATIAL 1D + TWO-WAY BTE  (N_layers+1 states; RK45 adaptive)
    !
    !  Runs from t_start to t_end with adaptive step-size control.
    !  State: ys(1..N_layers) = T_layer temperatures,  ys(N_layers+1) = WS
    !==========================================================================
    open(unit=30, file='simv0_spatial.dat', status='replace', action='write')
    write(30,'(A)') '# simv0_spatial.dat  [MODE_SPATIAL / RK45 Dormand-Prince adaptive]'
    write(30,'(A)') '# t(s)          T_1(K)         T_2(K)         T_3(K)         ' // &
                    'T_4(K)         T_5(K)         WS(m/s)        G_up(W/m2)     P_surf(W)'

    call init_spatial_y(ys)
    xs      = t_start
    hs      = h_default
    nok     = 0
    nbad    = 0
    Tp_s    = T_0
    Ey_s    = 0.0d0
    Gup_s   = 0.0d0

    !  Write initial state before integration begins
    Gup_s  = upwelling_top(ys(1:N_layers))
    Gnow   = geff(xs)
    etanow = eta_T_func(ys(1))
    Pout_s = etanow * Gnow * A_panel
    write(30,'(*(F15.6))') xs, ys(1:N_layers), ys(N_layers+1), Gup_s, Pout_s

    do while (xs < t_end - 1.0d-6)
        !  Clip to end of interval
        if (xs + hs > t_end) hs = t_end - xs

        !  Single RK45 stage
        call rk45_step(rhs_spatial, xs, ys, NSPA, hs, ysn, errn)

        !  Compute new step-size estimate  h_new = h * safety * (1/err)^0.2
        if (errn > 1.0d-10) then
            hnew = hs * rk45_safety * (1.0d0 / errn)**0.2d0
        else
            hnew = h_max
        end if
        hnew = min(max(hnew, h_min), h_max)

        if (errn <= 1.0d0) then
            !  Accept: advance time and state, then write output
            write(70,'(F15.3,2F15.6,A)') xs, hs, errn, '  ACCEPT'
            xs  = xs + hs
            ys  = ysn
            nok = nok + 1

            Gup_s  = upwelling_top(ys(1:N_layers))
            Gnow   = geff(xs)
            etanow = eta_T_func(ys(1))
            Pout_s = etanow * Gnow * A_panel
            write(30,'(*(F15.6))') xs, ys(1:N_layers), ys(N_layers+1), Gup_s, Pout_s
            if (ys(1) > Tp_s) Tp_s = ys(1)
            Ey_s = Ey_s + Pout_s * hs
        else
            !  Reject: keep state, reduce step
            write(70,'(F15.3,2F15.6,A)') xs, hs, errn, '  REJECT'
            nbad = nbad + 1
        end if

        hs = hnew

        if (nok + nbad > max_steps) then
            write(*,'(A,I0,A)') '  WARNING: max_steps (', max_steps, ') reached.'
            exit
        end if
    end do
    Hs_final = hs
    close(30)
    close(70)

    !==========================================================================
    !  5.  COMPARISON OUTPUT  (decoupled vs coupled at aligned time points)
    !
    !  The spatial model uses RK45 adaptive stepping (results in spatial.dat).
    !  A fixed h=100 s step is numerically unstable for the spatial model due
    !  to the thermal-diffusion CFL condition (dt_cfl ~ dz^2/alpha_diff ~ 1.4 s),
    !  so the spatial results are excluded from this aligned comparison.
    !==========================================================================
    open(unit=40, file='simv0_comparison.dat', status='replace', action='write')
    write(40,'(A)') '# simv0_comparison.dat  --  decoupled vs coupled at aligned times'
    write(40,'(A)') '#   (spatial results in simv0_spatial.dat; adaptive step was used)'
    write(40,'(A)') '# t(s)          T_dec(K)       T_coup(K)      WS(m/s)        ' // &
                    'G_eff(W/m2)    h_eff(W/m2K)   P_dec(W)       P_coup(W)'

    yd(1) = T_0
    yc(1) = T_0
    yc(2) = WS_0
    xcmp  = t_start

    do nstep = 0, nfixed
        Gnow   = geff(xcmp)
        hcnow  = h_eff_mixed(yc(1), yc(2))
        Pout_d = eta_stc            * 1000.0d0 * A_panel
        Pout_c = eta_T_func(yc(1)) * Gnow      * A_panel

        write(40,'(*(F15.6))') xcmp, yd(1), yc(1), yc(2), Gnow, hcnow, Pout_d, Pout_c

        if (nstep == nfixed) exit

        call euler_step(rhs_decoupled, xcmp, yd, NDEC, h_default, ydn)
        call rk4_step(  rhs_coupled,   xcmp, yc, NCOU, h_default, ycn)
        yd   = ydn
        yc   = ycn
        xcmp = xcmp + h_default
    end do
    close(40)

    !==========================================================================
    !  6.  SUMMARY STATISTICS  (stdout)
    !==========================================================================
    write(*,'(A)') '  Simulation complete.  Output files written:'
    write(*,'(A)') '    simv0_decoupled.dat    (Euler, const physics)'
    write(*,'(A)') '    simv0_coupled.dat      (RK4, BTE-NS coupled)'
    write(*,'(A)') '    simv0_spatial.dat      (RK45 adaptive, 1D + two-way BTE)'
    write(*,'(A)') '    simv0_comparison.dat   (three-mode side-by-side)'
    write(*,'(A)') '    simv0_diagnostic.dat   (method audit + RK45 step log)'
    write(*,'(A)')
    write(*,'(A)') '  ---- Peak panel temperatures ----'
    write(*,'(A,F8.2,A)') '    Decoupled (Euler):  ', Tp_d, ' K'
    write(*,'(A,F8.2,A)') '    Coupled   (RK4):    ', Tp_c, ' K'
    write(*,'(A,F8.2,A)') '    Spatial   (RK45):   ', Tp_s, ' K'
    write(*,'(A)')
    write(*,'(A)') '  ---- Daily electrical energy yield ----'
    write(*,'(A,F10.1,A)') '    Decoupled:  ', Ey_d / 3600.0d0, ' Wh'
    write(*,'(A,F10.1,A)') '    Coupled:    ', Ey_c / 3600.0d0, ' Wh'
    write(*,'(A,F10.1,A)') '    Spatial:    ', Ey_s / 3600.0d0, ' Wh'
    write(*,'(A)')
    write(*,'(A)') '  ---- RK45 adaptive stepping (spatial mode) ----'
    write(*,'(A,I8)')     '    Steps accepted:  ', nok
    write(*,'(A,I8)')     '    Steps rejected:  ', nbad
    write(*,'(A,F8.1,A)') '    Final step size: ', Hs_final, ' s'
    write(*,'(A)')
    write(*,'(A)') '  ---- Two-way BTE coupling ----'
    write(*,'(A,F6.2,A)') '    Final upwelling flux (top surface): ', Gup_s, ' W/m^2'
    write(*,'(A)') '    (see simv0_spatial.dat column G_up for upwelling history)'
    write(*,'(A)')

contains

    !--------------------------------------------------------------------------
    !  rhs_decoupled  –  1-state RHS: constant G=1000, h_conv=15, eta=eta_stc
    !
    !  dT/dt = (1/mcp) * [(alpha-eta)*G*A - h_conv*A*(T-T_amb)
    !                    - eps*sigma*A*(T^4 - T_sky^4)]
    !
    !  Matches rk4_ode.f90 nonlinear ODE exactly (decoupled baseline).
    !--------------------------------------------------------------------------
    subroutine rhs_decoupled(t, y, n, dydt)
        implicit none
        real(8), intent(in)  :: t
        integer, intent(in)  :: n
        real(8), intent(in)  :: y(n)
        real(8), intent(out) :: dydt(n)

        real(8), parameter :: G_cst  = 1000.0d0   ! constant irradiance (W/m^2)
        real(8), parameter :: hcv    = 15.0d0      ! constant convection (W/m^2/K)

        !  Suppress unused-variable warning for dummy t (no time dependence here)
        real(8) :: t_dummy
        t_dummy = t

        dydt(1) = (1.0d0 / (m_pan * cp_pan)) *                             &
                  (   (alpha_panel - eta_stc) * G_cst  * A_panel             &
                    - hcv                     * A_panel * (y(1) - T_amb)    &
                    - eps_panel * SIGMA        * A_panel * (y(1)**4 - T_sky**4) )

    end subroutine rhs_decoupled

    !--------------------------------------------------------------------------
    !  rhs_coupled  –  2-state RHS: full BTE-NS with eta(T) and WS ODE
    !
    !  dT/dt  = (1/mcp)*[(alpha-eta(T))*G_eff(t)*A
    !                   - h_eff(T,WS)*A*(T-T_amb)
    !                   - eps*sigma*A*(T^4-T_sky^4)]
    !  dWS/dt = -lambda_BL*(WS - WS_geo) + gamma_BL*(T - T_amb)
    !
    !  Matches bte_ns_ode.f90 coupled RHS exactly for validation.
    !
    !  Note on WS sign: WS is NOT clamped here so the BL ODE converges to its
    !  physical equilibrium WS_eq = WS_geo + (gamma_BL/lambda_BL)*(T-T_amb),
    !  which is negative at night (panel cooler than T_amb).  The h_eff_mixed
    !  function clamps WS internally before evaluating the convection coefficient.
    !--------------------------------------------------------------------------
    subroutine rhs_coupled(t, y, n, dydt)
        implicit none
        real(8), intent(in)  :: t
        integer, intent(in)  :: n
        real(8), intent(in)  :: y(n)
        real(8), intent(out) :: dydt(n)

        real(8) :: Tv, WSv, Gt, hcv, etav

        Tv  = y(1)
        WSv = y(2)           ! unclamped: h_eff_mixed clamps WS internally

        Gt   = geff(t)
        hcv  = h_eff_mixed(Tv, WSv)   ! h_eff_mixed uses max(0, WSv)
        etav = eta_T_func(Tv)

        dydt(1) = (1.0d0 / (m_pan * cp_pan)) *                            &
                  (   (alpha_panel - etav) * Gt      * A_panel             &
                    - hcv                 * A_panel  * (Tv - T_amb)       &
                    - eps_panel * SIGMA   * A_panel  * (Tv**4 - T_sky**4) )

        dydt(2) = -lambda_BL * (WSv - WS_geo) + gamma_BL * (Tv - T_amb)

    end subroutine rhs_coupled

    !--------------------------------------------------------------------------
    !  rhs_spatial  –  (N_layers+1)-state RHS: thin wrapper around spatial_rhs
    !--------------------------------------------------------------------------
    subroutine rhs_spatial(t, y, n, dydt)
        implicit none
        real(8), intent(in)  :: t
        integer, intent(in)  :: n
        real(8), intent(in)  :: y(n)
        real(8), intent(out) :: dydt(n)

        call spatial_rhs(t, y, n, dydt)

    end subroutine rhs_spatial

end program simv0_driver
