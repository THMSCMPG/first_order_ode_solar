#!/usr/bin/env python3
"""
run_simv0.py  --  Entry-point for the simv0 ODE Integration Platform

Usage
-----
    python run_simv0.py              # build (if needed) + run all modes
    python run_simv0.py --build-only # only compile Fortran sources

Options
-------
    --build-only        Only compile the Fortran sources (make).
    --simv0-dir DIR     Path to the simv0 source directory.
                        Defaults to the directory containing this script.
    -h, --help          Show this help message and exit.

Examples
--------
    python run_simv0.py
    python run_simv0.py --build-only

To visualise results:
    gnuplot plots/plot_simv0.gp
"""

import argparse
import datetime
import json
import os
import platform
import shutil
import subprocess
import sys

# ---------------------------------------------------------------------------
#  Helpers
# ---------------------------------------------------------------------------

BANNER = """\
==============================================================
  simv0 — Intelligent ODE Integration Platform
  Solar Panel Thermal Model  |  Euler · Heun · RK4 · RK45
=============================================================="""

OUTPUT_FILES = [
    "simv0_decoupled.dat",
    "simv0_coupled.dat",
    "simv0_spatial.dat",
    "simv0_comparison.dat",
    "simv0_diagnostic.dat",
]


def _binary_path(simv0_dir):
    return os.path.join(simv0_dir, "bin", "simv0")


def _run(cmd, cwd, label=""):
    """Run *cmd* in *cwd*, inheriting stdout/stderr so output streams live
    to the terminal, and return the exit code."""
    if label:
        print(f"\n  {label}")
        print("  " + "-" * 56)
    result = subprocess.run(cmd, cwd=cwd)
    return result.returncode


def _build(simv0_dir):
    """Compile Fortran sources via `make`."""
    print("  Building simv0 …")
    rc = _run(["make"], cwd=simv0_dir, label="make")
    if rc != 0:
        print(f"\n  ERROR: make returned exit code {rc}", file=sys.stderr)
        return False
    print("  ✓ Build successful.")
    return True


def _run_simulation(simv0_dir):
    """Execute the simulation via `make run`."""
    print("  Running simulation …")
    rc = _run(["make", "run"], cwd=simv0_dir, label="make run")
    if rc != 0:
        print(f"\n  ERROR: simulation returned exit code {rc}", file=sys.stderr)
        return False
    print("  ✓ Simulation completed.")
    return True


def _parse_dat_file(path):
    """
    Parse key metrics from a simv0 .dat output file.

    Returns a dict mapping column name → last-row value (float), or an empty
    dict if the file is missing or contains no numeric rows.
    """
    if not os.path.isfile(path):
        return {}

    header = []
    last_row = []

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                cols = line.lstrip("#").split()
                if cols:
                    header = cols
            else:
                try:
                    last_row = [float(x) for x in line.split()]
                except ValueError:
                    pass

    if not last_row:
        return {}

    if header and len(header) == len(last_row):
        return dict(zip(header, last_row))
    return {f"col_{i}": val for i, val in enumerate(last_row)}


def _peak_value(path, col_index):
    """Return the maximum value of column *col_index* across all data rows."""
    if not os.path.isfile(path):
        return None
    peak = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                vals = [float(x) for x in line.split()]
                if col_index < len(vals):
                    v = vals[col_index]
                    if peak is None or v > peak:
                        peak = v
            except ValueError:
                pass
    return peak


def _parse_diagnostic(path):
    """
    Parse `simv0_diagnostic.dat` for RK45 step statistics.

    The diagnostic file records one row per RK45 step attempt (ACCEPT/REJECT):
      t(s)  h_used(s)  err_norm  status
    Returns a dict with rk45_steps_accepted, rk45_steps_rejected, and
    h_final_s.  Missing file → empty dict.
    """
    if not os.path.isfile(path):
        return {}

    accepted = 0
    rejected = 0
    h_final = None

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            # The last token is ACCEPT or REJECT; h_used is the second token
            # (may have '***' appended when err_norm overflows the format field)
            if len(parts) < 2:
                continue
            status = parts[-1].upper()
            if status not in ("ACCEPT", "REJECT"):
                continue
            try:
                # Strip any trailing non-numeric characters from h field
                h_str = parts[1].rstrip("*")
                h = float(h_str)
            except (ValueError, IndexError):
                # Count rejected even if we can't parse h (first overflow row)
                if status == "REJECT":
                    rejected += 1
                continue
            if status == "ACCEPT":
                accepted += 1
                h_final = h
            elif status == "REJECT":
                rejected += 1

    result = {
        "rk45_steps_accepted": accepted,
        "rk45_steps_rejected": rejected,
    }
    if h_final is not None:
        result["h_final_s"] = round(h_final, 4)
    return result


def _write_metadata(simv0_dir):
    """
    Collect run context and key metrics, then write `data/simv0_metadata.json`.
    """
    data_dir = os.path.join(simv0_dir, "data")
    if not os.path.isdir(data_dir):
        print("  (data/ not found; skipping metadata write)")
        return

    # ── Parameters (from simv0_config.f90 defaults) ─────────────────────────
    parameters = {
        "h_default_s": 100,
        "atol_K": 0.1,
        "rtol": 1e-6,
        "N_layers": 5,
        "t_start_s": 0,
        "t_end_s": 86400,
    }

    # ── Summary statistics from output files ─────────────────────────────────
    dec_path  = os.path.join(data_dir, "simv0_decoupled.dat")
    coup_path = os.path.join(data_dir, "simv0_coupled.dat")
    spat_path = os.path.join(data_dir, "simv0_spatial.dat")
    diag_path = os.path.join(data_dir, "simv0_diagnostic.dat")

    # Peak temperatures (column index 1 = T)
    tp_dec  = _peak_value(dec_path,  1)
    tp_coup = _peak_value(coup_path, 1)
    tp_spat = _peak_value(spat_path, 1)

    # Energy yield: last row P_elec column, then integrate via trapezoidal
    # approximation using the last-row entry as a proxy for the final power
    # (full integration is done in the Fortran driver; we read the last value
    # of energy columns if present, otherwise we skip).
    last_dec  = _parse_dat_file(dec_path)
    last_coup = _parse_dat_file(coup_path)
    last_spat = _parse_dat_file(spat_path)

    diag_stats = _parse_diagnostic(diag_path)

    summary = {}
    if tp_dec  is not None:
        summary["peak_T_decoupled_K"]  = round(tp_dec,  2)
    if tp_coup is not None:
        summary["peak_T_coupled_K"]    = round(tp_coup, 2)
    if tp_spat is not None:
        summary["peak_T_spatial_K"]    = round(tp_spat, 2)

    # Last-row P_elec as the final instantaneous power (W)
    if last_dec and "P_elec(W)" in last_dec:
        summary["final_P_decoupled_W"] = round(last_dec["P_elec(W)"], 2)
    if last_coup and "P_elec(W)" in last_coup:
        summary["final_P_coupled_W"]   = round(last_coup["P_elec(W)"], 2)

    summary.update(diag_stats)

    # ── Compiler info ─────────────────────────────────────────────────────────
    gfortran_version = None
    gfc = shutil.which("gfortran")
    if gfc:
        try:
            out = subprocess.check_output(
                ["gfortran", "--version"], stderr=subprocess.STDOUT, text=True
            )
            gfortran_version = out.splitlines()[0].strip()
        except Exception:
            pass

    metadata = {
        "timestamp": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "hostname": platform.node(),
        "fortran_compiler": gfortran_version,
        "simulation_mode": "all (decoupled + coupled + spatial)",
        "parameters": parameters,
        "summary": summary,
        "output_files": OUTPUT_FILES,
    }

    out_path = os.path.join(data_dir, "simv0_metadata.json")
    with open(out_path, "w") as fh:
        json.dump(metadata, fh, indent=2)
        fh.write("\n")

    print(f"  ✓ Metadata written → {out_path}")


def _print_results(simv0_dir):
    """Print a brief metrics summary from the .dat output files."""
    data_dir = os.path.join(simv0_dir, "data")
    if not os.path.isdir(data_dir):
        print("  (data/ directory not found; skipping results summary)")
        return

    dat_files = [
        "simv0_decoupled.dat",
        "simv0_coupled.dat",
        "simv0_spatial.dat",
        "simv0_comparison.dat",
    ]

    print("\n  ---- Key metrics from output files ----")
    for fname in dat_files:
        path = os.path.join(data_dir, fname)
        metrics = _parse_dat_file(path)
        if metrics:
            print(f"\n  {fname}")
            for col, val in metrics.items():
                print(f"    {col:<28s} {val:.4g}")


# ---------------------------------------------------------------------------
#  CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        prog="run_simv0.py",
        description="Run the simv0 ODE integration platform (headless).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--build-only",
        action="store_true",
        help="Only compile — do not run the simulation.",
    )
    parser.add_argument(
        "--simv0-dir",
        default=None,
        help="Path to the simv0 source directory (default: directory of this script).",
    )

    args = parser.parse_args()

    simv0_dir = args.simv0_dir or os.path.dirname(os.path.abspath(__file__))

    print(BANNER)
    print(f"  simv0 directory : {simv0_dir}")
    binary = _binary_path(simv0_dir)
    binary_exists = os.path.isfile(binary)
    print(f"  Binary          : {'found' if binary_exists else 'not found (will build)'}")
    print()

    # Build if needed
    if not binary_exists or args.build_only:
        ok = _build(simv0_dir)
        if not ok:
            sys.exit(1)
        if args.build_only:
            print("\n  ✓ Build-only run complete.")
            return

    # Run simulation
    ok = _run_simulation(simv0_dir)
    if not ok:
        sys.exit(1)

    # Write metadata JSON
    _write_metadata(simv0_dir)

    # Print metrics summary
    _print_results(simv0_dir)

    print("\n" + "=" * 62)
    print("  simv0 run complete.")
    print("  To visualise: gnuplot plots/plot_simv0.gp")
    print("=" * 62)


if __name__ == "__main__":
    main()
