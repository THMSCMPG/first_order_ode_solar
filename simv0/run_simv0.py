#!/usr/bin/env python3
"""
run_simv0.py  --  Entry-point for the simv0 ODE Integration Platform

Usage
-----
    python run_simv0.py              # launch Tkinter GUI (default)
    python run_simv0.py --no-gui     # headless, build + run via make

Options
-------
    --no-gui            Run without a graphical interface.
    --build-only        Only compile the Fortran sources (make).
    --scenarios LIST    Comma-separated list of scenarios to run in
                        headless mode (euler,rk4,atmospheric,diurnal,
                        coupled,mc,bte_ns).  Default: all.
    --simv0-dir DIR     Path to the simv0 source directory.
                        Defaults to the directory containing this script.
    -h, --help          Show this help message and exit.

Examples
--------
    python run_simv0.py
    python run_simv0.py --no-gui
    python run_simv0.py --no-gui --scenarios coupled,bte_ns
    python run_simv0.py --no-gui --build-only
"""

import argparse
import os
import subprocess
import sys

# ---------------------------------------------------------------------------
#  Helpers
# ---------------------------------------------------------------------------

SCENARIOS_ALL = ["euler", "rk4", "atmospheric", "diurnal", "coupled", "mc", "bte_ns"]

BANNER = """\
==============================================================
  simv0 — Intelligent ODE Integration Platform
  Solar Panel Thermal Model  |  Euler · Heun · RK4 · RK45
=============================================================="""


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

    Returns a dict mapping column name → last-row value string, or an empty
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
        return {name: f"{val:.4g}" for name, val in zip(header, last_row)}
    return {f"col_{i}": f"{val:.4g}" for i, val in enumerate(last_row)}


def _print_results(simv0_dir, scenarios):
    """Print a brief metrics summary from the .dat output files."""
    data_dir = os.path.join(simv0_dir, "data")
    if not os.path.isdir(data_dir):
        print("  (data/ directory not found; skipping results summary)")
        return

    file_map = {
        "euler":       "simv0_decoupled.dat",
        "rk4":         "simv0_decoupled.dat",
        "atmospheric": "simv0_coupled.dat",
        "diurnal":     "simv0_coupled.dat",
        "coupled":     "simv0_coupled.dat",
        "mc":          "simv0_comparison.dat",
        "bte_ns":      "simv0_spatial.dat",
    }

    seen = set()
    print("\n  ---- Key metrics from output files ----")
    for key in scenarios:
        fname = file_map.get(key)
        if fname and fname not in seen:
            seen.add(fname)
            metrics = _parse_dat_file(os.path.join(data_dir, fname))
            if metrics:
                print(f"\n  {fname}")
                for col, val in metrics.items():
                    print(f"    {col:<28s} {val}")


# ---------------------------------------------------------------------------
#  CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        prog="run_simv0.py",
        description="Launch the simv0 ODE integration platform.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--no-gui",
        action="store_true",
        help="Run headless (no Tkinter window).",
    )
    parser.add_argument(
        "--build-only",
        action="store_true",
        help="Only compile — do not run the simulation.",
    )
    parser.add_argument(
        "--scenarios",
        default=",".join(SCENARIOS_ALL),
        help=(
            "Comma-separated scenarios for headless mode "
            "(default: all).  Ignored when launching the GUI."
        ),
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

    # ── GUI mode ─────────────────────────────────────────────────────────────
    if not args.no_gui:
        try:
            from simv0_launcher_gui import launch
        except ImportError:
            # Allow running from a different working directory
            sys.path.insert(0, simv0_dir)
            from simv0_launcher_gui import launch

        print("  Launching GUI …\n")
        launch(simv0_dir=simv0_dir)
        return

    # ── Headless mode ────────────────────────────────────────────────────────
    scenarios = [s.strip() for s in args.scenarios.split(",") if s.strip()]
    unknown = [s for s in scenarios if s not in SCENARIOS_ALL]
    if unknown:
        print(
            f"  WARNING: unknown scenario(s) ignored: {', '.join(unknown)}",
            file=sys.stderr,
        )
        scenarios = [s for s in scenarios if s in SCENARIOS_ALL]

    if not scenarios:
        print("  ERROR: no valid scenarios selected.", file=sys.stderr)
        sys.exit(1)

    print(f"  Scenarios       : {', '.join(scenarios)}")
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

    # Print metrics summary
    _print_results(simv0_dir, scenarios)

    print("\n" + "=" * 62)
    print("  simv0 run complete.")
    print("=" * 62)


if __name__ == "__main__":
    main()
