# ===========================================================================
#  plot_simv0.gp  --  Gnuplot visualisation scripts for simv0 output
#
#  Usage (run from the plots/ directory):
#    cd plots && gnuplot plot_simv0.gp
#
#  Produces PNG images in plots/:
#    01_temperature_comparison.png  : T_dec vs T_coup vs T_surf(spatial)
#    02_wind_and_irradiance.png     : WS and G_eff from coupled run
#    03_spatial_layers.png          : Layer temperature profiles (day & night)
#    04_upwelling_bte.png           : Two-way BTE upwelling flux
#    05_power_output.png            : Electrical power (all modes)
# ===========================================================================

set terminal png size 1200,700
set style data lines
set grid
set xrange [0:86400]
set format x "%g"

t_hours(t) = t / 3600.0

# ---------------------------------------------------------------------------
# 1.  Temperature comparison  (decoupled vs coupled vs spatial surface)
# ---------------------------------------------------------------------------
set output "01_temperature_comparison.png"
set title "Panel Surface Temperature – Three-Mode Comparison\n(simv0: Euler / RK4 / RK45-Adaptive)"
set xlabel "Time (hours)"
set ylabel "Temperature (K)"
set xtics ("0h" 0, "3h" 10800, "6h" 21600, "9h" 32400, "12h" 43200, \
           "15h" 54000, "18h" 64800, "21h" 75600, "24h" 86400)

# Comparison file columns: t, T_dec, T_coup, WS, G_eff, h_eff, P_dec, P_coup
plot "../data/simv0_comparison.dat" using 1:2 title "Decoupled / Euler (const G=1000)" \
         lw 2 lc rgb "#CC0000" dt 2, \
     "../data/simv0_comparison.dat" using 1:3 title "Coupled BTE-NS / RK4"  \
         lw 2 lc rgb "#0055CC", \
     "../data/simv0_spatial.dat"    using 1:2 title "Spatial layer 1 / RK45" \
         lw 2 lc rgb "#00AA44"

# ---------------------------------------------------------------------------
# 2.  Wind speed and G_eff  (from coupled run)
# ---------------------------------------------------------------------------
set output "02_wind_and_irradiance.png"
set title "Dynamic Wind Speed and Effective Irradiance\n(simv0 coupled BTE-NS / RK4)"
set ylabel "G_{eff}  (W/m^2) / WS (m/s)"
set yrange [-20:1100]

plot "../data/simv0_comparison.dat" using 1:5 title "G_{eff}  (W/m^2)" \
         lw 2 lc rgb "#FF8800" axes x1y1, \
     "../data/simv0_comparison.dat" using 1:4 title "WS  (m/s) [×40 scale]" \
         lw 2 lc rgb "#0055CC" axes x1y1

# ---------------------------------------------------------------------------
# 3.  Spatial layer temperature profiles  (selected snapshots)
# ---------------------------------------------------------------------------
set output "03_spatial_layers.png"
set title "1D Spatial Panel Temperature (all 5 layers)\n(simv0 MODE_SPATIAL / RK45 adaptive)"
set ylabel "Temperature (K)"
unset yrange

# Spatial file columns: t, T1, T2, T3, T4, T5, WS, G_up, P
plot "../data/simv0_spatial.dat" using 1:2 title "Layer 1 (glass, top)"  lw 2 lc rgb "#CC0000", \
     "../data/simv0_spatial.dat" using 1:3 title "Layer 2 (EVA top)"     lw 1 lc rgb "#FF6600" dt 2, \
     "../data/simv0_spatial.dat" using 1:4 title "Layer 3 (silicon cell)" lw 2 lc rgb "#228B22", \
     "../data/simv0_spatial.dat" using 1:5 title "Layer 4 (EVA bot)"     lw 1 lc rgb "#0055CC" dt 2, \
     "../data/simv0_spatial.dat" using 1:6 title "Layer 5 (back sheet)"   lw 2 lc rgb "#880088"

# ---------------------------------------------------------------------------
# 4.  Two-way BTE upwelling flux at top surface
# ---------------------------------------------------------------------------
set output "04_upwelling_bte.png"
set title "Two-Way BTE: Upwelling Thermal Flux from Panel Interior\n(escaping through top surface; simv0 spatial mode)"
set ylabel "Upwelling flux  G_{up}  (W/m^2)"
set yrange [*:*]

# Spatial file column 8 = G_up
plot "../data/simv0_spatial.dat" using 1:8 title "G_{up} (upwelling, two-way BTE)" \
         lw 2 lc rgb "#CC0000"

# ---------------------------------------------------------------------------
# 5.  Electrical power output  (all modes)
# ---------------------------------------------------------------------------
set output "05_power_output.png"
set title "Electrical Power Output\n(simv0: Euler / RK4 / RK45-Adaptive spatial)"
set ylabel "Power  (W)"
set yrange [0:*]

# Comparison: col 7=P_dec, col 8=P_coup; Spatial: col 9=P_surf
plot "../data/simv0_comparison.dat" using 1:7 title "Decoupled / Euler" \
         lw 2 lc rgb "#CC0000" dt 2, \
     "../data/simv0_comparison.dat" using 1:8 title "Coupled BTE-NS / RK4" \
         lw 2 lc rgb "#0055CC", \
     "../data/simv0_spatial.dat"    using 1:9 title "Spatial layer 1 / RK45" \
         lw 2 lc rgb "#00AA44"

print "Plots written to plots/ directory."
