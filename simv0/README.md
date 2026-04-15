# simv0 — Intelligent ODE Integration Platform

`simv0` is a comprehensive Fortran integration platform for the solar-panel
thermal model.  It unifies all ODE methods developed across the `Code_Dev`
subdirectories, adds intelligent method selection, implements Dormand-Prince
RK45 with adaptive step-size control, and introduces two new physical models:

- **Spatial 1D thermal PDE** — temperature stratification through the panel depth  
- **Two-way BTE coupling** — downwelling + upwelling radiation fluxes linked to
  interior layer temperatures

---

## Quick start

```bash
cd simv0/
make        # compile (requires gfortran)
./simv0     # run the default 24-hour simulation
```

Output files appear in the working directory; sample results are in `data/`.

---

## Module structure

| File | Description |
|------|-------------|
| `simv0_config.f90` | Physical constants and tunable parameters |
| `simv0_method_selector.f90` | Automatic integration-method selection |
| `simv0_ode_suite.f90` | Euler, Heun, RK4, RK45 (Dormand-Prince) |
| `simv0_bte_coupling.f90` | BTE physics: G_eff, h_eff, η(T), downwelling/upwelling |
| `simv0_spatial_1d.f90` | 1D panel heat PDE (Method of Lines) |
| `simv0_driver.f90` | Main program: orchestrates all modules |

Compilation order: `simv0_config` → `simv0_method_selector` →
`simv0_ode_suite` → `simv0_bte_coupling` → `simv0_spatial_1d` →
`simv0_driver`.

---

## Intelligent method selection

The `select_method` subroutine in `simv0_method_selector.f90` scores the
problem against seven rules and assigns the **minimum sufficient** method:

| Priority | Rule | Selected method |
|----------|------|-----------------|
| 1 | Spatial 1D PDE (MODE_SPATIAL) | **RK45** |
| 2 | Very high accuracy (tol < 0.01 K) | **RK45** |
| 3 | Stiffness: τ/h < 5 | **RK4** |
| 4 | Coupled BTE-NS (MODE_COUPLED) | **RK4** |
| 5 | Moderate accuracy (tol < 0.1 K) | **RK4** |
| 6 | Any coupling (MODE_COUPLED or higher) | **Heun** |
| 7 | Default | **Euler** |

The thermal time constant τ = m·c_p / (h_ref·A) ≈ 395 s.  With the default
step h = 100 s, τ/h ≈ 3.9 < 5 and the stiffness rule promotes the method to
at least RK4 for all modes.  To obtain Euler, use h ≤ 50 s and tol ≥ 0.5 K
with MODE_DECOUPLED (see `configs/fast.cfg`).

---

## Simulation modes

### MODE_DECOUPLED (constant-physics baseline)
- 1-state ODE: dT/dt = f(T)  with constant G = 1000 W/m², h = 15 W/m²/K, η = η_STC
- Demonstrates Forward Euler as the method of choice for the simplest physics
- Intentionally does *not* use RK4 here to contrast Euler vs RK4 outputs

### MODE_COUPLED (BTE-NS tight coupling, 2-state)
- 2-state ODE: [dT/dt, dWS/dt]
- BTE: Beer-Lambert diurnal G_eff(t), Skoplaki-Palyvos η(T)
- NS: Churchill-Usagi mixed convection h_eff(T,WS), BL relaxation ODE for WS
- Uses Classical RK4 (h = 100 s)
- **Validates against `bte_ns_ode.f90`**: peak T = 306.2 K, yield = 1857.1 Wh ✓

### MODE_SPATIAL (1D spatial PDE + two-way BTE, adaptive)
- (N_layers + 1)-state system: [T₁, T₂, T₃, T₄, T₅, WS]
- Method of Lines discretisation of the 1D heat equation:
  ρ·c_p·(∂T/∂t) = κ·(∂²T/∂z²) + Q_solar(z,t)
- Beer-Lambert solar absorption per layer (down-going flux)
- Top surface: convection + Stefan-Boltzmann radiation
- Bottom surface: insulated (zero-flux BC)
- Two-way BTE: upwelling flux G_up from interior thermal emission (diagnostic)
- Uses **RK45 Dormand-Prince with adaptive step size** (h_min = 1 s, h_max = 600 s)

---

## RK45 Dormand-Prince (adaptive)

Implemented in `simv0_ode_suite.f90` using the standard DOPRI5 Butcher tableau:

- **7 stages** per step (6 for the 5th-order solution + 1 for error estimate)
- Error estimated as the difference between the embedded 4th- and 5th-order
  solutions using the coefficients `e₁..e₇`
- **WRMS norm**: `err = sqrt(Σ(eᵢ / scaleᵢ)² / n)`, `scaleᵢ = atol + rtol·|yᵢ|`
- **Step accepted** if `err ≤ 1`; **step rejected** (and retried) if `err > 1`
- **New step size**: `h_new = h · 0.9 · (1/err)^0.2`, clamped to [h_min, h_max]

Default tolerances: `atol = 0.1 K`, `rtol = 1×10⁻⁶`.

For the spatial system (N_layers = 5 cells, dz = 2 mm) the thermal-diffusion
CFL stability limit is approximately:

    dt_CFL = dz² / (2·α) = (0.002)² / (2 × 1.48×10⁻⁶) ≈ 1.35 s

The adaptive RK45 automatically honours this limit (typical accepted step
≈ 2–3 s during peak irradiance).

---

## Spatial 1D model — physics

The finite-volume Method of Lines gives one ODE per cell:

**Layer 1 (top, z = 0)** — exposed to atmosphere:
```
ρ·cp·dz·(dT₁/dt) = κ/dz·(T₂−T₁) − h_eff·(T₁−T_amb) − ε·σ·(T₁⁴−T_sky⁴) + Q_sol,1
```

**Interior layers i = 2 … N−1**:
```
ρ·cp·dz·(dT_i/dt) = κ/dz·(T_{i-1} − 2·T_i + T_{i+1}) + Q_sol,i
```

**Layer N (bottom, insulated)**:
```
ρ·cp·dz·(dT_N/dt) = κ/dz·(T_{N-1} − T_N) + Q_sol,N
```

Solar absorbed per layer (Beer-Lambert, down-going):
```
Q_sol,i = G_eff(t) · [exp(−κ_abs·(i−1)·dz) − exp(−κ_abs·i·dz)]   [W/m²]
```

Wind speed ODE (NS boundary-layer dynamics):
```
dWS/dt = −λ·(WS − WS_geo) + γ·(T₁ − T_amb)
```

---

## Two-way BTE coupling

The upwelling flux leaving the top surface is computed from the
Stefan-Boltzmann emission of every layer attenuated by the optical
transmittance above it:

```
G_up = ε·σ · Σᵢ Tᵢ⁴ · (1 − exp(−κ_abs·dz)) · exp(−κ_abs·(i−1)·dz)
```

This represents the thermal radiation from the panel interior that escapes
through the top glass.  Its magnitude (≈ 215–240 W/m²) quantifies the
**upwelling feedback**: higher interior temperatures → more upwelling
→ increased net radiative loss from the system.

The `G_up` column in `simv0_spatial.dat` provides a time-series of this
two-way BTE diagnostic.

---

## Output files

| File | Contents |
|------|----------|
| `simv0_decoupled.dat` | Euler baseline: t, T, P_elec |
| `simv0_coupled.dat` | RK4 BTE-NS: t, T, WS, G_eff, h_eff, η, P_elec |
| `simv0_spatial.dat` | RK45 spatial: t, T₁…T₅, WS, G_up, P_surf |
| `simv0_comparison.dat` | Aligned decoupled vs coupled: t, T_dec, T_coup, WS, G_eff, h_eff, P_dec, P_coup |
| `simv0_diagnostic.dat` | Method audit + RK45 step-size history |

---

## Visualisation

```bash
gnuplot plots/plot_simv0.gp
```

Produces five PNG images in `plots/`:
1. `01_temperature_comparison.png` — three-mode temperature comparison
2. `02_wind_and_irradiance.png` — dynamic WS and G_eff
3. `03_spatial_layers.png` — all five layer temperatures
4. `04_upwelling_bte.png` — two-way BTE upwelling flux
5. `05_power_output.png` — electrical power (all modes)

---

## Example configurations

| Config | h (s) | tol (K) | Mode | Method selected |
|--------|-------|---------|------|-----------------|
| `fast.cfg` | 50 | 0.5 | DECOUPLED | **Euler** |
| `balanced.cfg` | 100 | 0.1 | COUPLED | **RK4** (stiffness) |
| `high_accuracy.cfg` | 10 | 0.01 | COUPLED | **RK4** → RK45 if tol < 0.01 |
| `full_physics.cfg` | adaptive | 0.1 | SPATIAL | **RK45** always |

---

## Validation

The coupled simulation (`bte_ns_ode.f90` reference → `simv0` RK4 coupled):

| Quantity | `bte_ns_ode.f90` | `simv0` (RK4) | Match |
|----------|-----------------|---------------|-------|
| Peak panel temperature | 306.22 K | 306.22 K | ✓ |
| Daily energy yield | 1857.07 Wh | 1857.1 Wh | ✓ |
| Final WS (night) | −12.416 m/s | −12.416 m/s | ✓ |

Spatial model (RK45, 5 layers):
- Peak surface temperature: 304.75 K (slightly lower due to heat redistribution)
- Daily energy yield: 1868.7 Wh (marginally higher; lower surface T → higher η)
- Upwelling G_up at end of day: 214.3 W/m²

---

## Open research problems addressed

### 1. Spatial effects (1D PDE along panel depth)
Implemented in `simv0_spatial_1d.f90`.  The Method of Lines discretisation
produces a coupled system of ODEs—one per layer—capturing:
- Temperature stratification within the panel
- Radial (depth-wise) heat conduction between layers
- Layer-specific solar absorption via Beer-Lambert
- Improved top-surface boundary condition (uses live T₁ not a fixed T)

### 2. Two-way irradiance–temperature coupling
Implemented in `simv0_bte_coupling.f90` (`upwelling_top` function).
The forward Beer-Lambert pass distributes solar to each layer; the backward
pass computes G_up from interior temperatures.  The resulting upwelling flux
is a direct measure of the panel-microclimate radiative feedback loop.

---

## Building

```
gfortran -O2 -Wall -Wno-unused-variable \
  simv0_config.f90 simv0_method_selector.f90 simv0_ode_suite.f90 \
  simv0_bte_coupling.f90 simv0_spatial_1d.f90 simv0_driver.f90 \
  -o simv0
```

Or simply:
```
make
```

Tested with `gfortran` 9–14.
