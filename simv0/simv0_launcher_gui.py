"""
simv0_launcher_gui.py  --  Tkinter GUI for the simv0 ODE Integration Platform

Provides a retro Windows XP / early-2000s style interface for launching
simv0 simulations, configuring solver parameters, and viewing live console
output.  Mirrors the aesthetic of the AURA-MFP launcher GUI.

Usage (via run_simv0.py):
    python run_simv0.py           # launches this GUI
    python run_simv0.py --no-gui  # headless execution

Direct use:
    python simv0_launcher_gui.py
"""

import os
import subprocess
import threading
import queue
import tkinter as tk
from tkinter import messagebox

# ---------------------------------------------------------------------------
#  Retro palette  (Windows XP / early-2000s aesthetic, matching AURA-MFP)
# ---------------------------------------------------------------------------
PAL = {
    "navy":        "#003366",
    "navy_dark":   "#001F4A",
    "navy_mid":    "#0A4080",
    "beige":       "#F0EAD6",
    "gray":        "#D4D0C8",
    "gray_dark":   "#ACA899",
    "gray_med":    "#C8C4BC",
    "white":       "#FFFFFF",
    "console_bg":  "#0A0E1A",
    "phosphor":    "#33FF00",
    "error":       "#FF4444",
    "success":     "#00FF88",
    "warning":     "#FFCC00",
    "header":      "#66CCFF",
    "title_fg":    "#FFFFFF",
    "btn_bg":      "#D4D0C8",
    "btn_active":  "#B8B4AC",
    "label_fg":    "#1A1A2E",
    "sep":         "#8A8680",
}

FONT_LABEL   = ("Arial",       9)
FONT_BOLD    = ("Arial",       9,  "bold")
FONT_TITLE   = ("Arial",       11, "bold")
FONT_SECTION = ("Arial",       8,  "bold")
FONT_CONSOLE = ("Courier New", 9)
FONT_FIXED   = ("Courier New", 9)


def _make_btn(parent, text, command, width=14):
    """Return a classic raised-relief button matching the retro palette."""
    return tk.Button(
        parent,
        text=text,
        command=command,
        width=width,
        font=FONT_BOLD,
        bg=PAL["btn_bg"],
        fg=PAL["label_fg"],
        activebackground=PAL["btn_active"],
        activeforeground=PAL["label_fg"],
        relief=tk.RAISED,
        bd=2,
        cursor="hand2",
    )


def _make_label(parent, text, bold=False, fg=None):
    """Return a label styled for the beige panel background."""
    return tk.Label(
        parent,
        text=text,
        font=FONT_BOLD if bold else FONT_LABEL,
        bg=PAL["beige"],
        fg=fg or PAL["label_fg"],
    )


# ---------------------------------------------------------------------------
#  Scenario definitions
# ---------------------------------------------------------------------------
SCENARIOS = [
    ("euler",       "Euler baseline       (decoupled, const G)"),
    ("rk4",         "RK4 comparison       (decoupled, const G)"),
    ("atmospheric", "Atmospheric BTE      (Beer-Lambert G_eff)"),
    ("diurnal",     "Diurnal cycle        (24-hour G_eff(t))"),
    ("coupled",     "BTE-NS coupled       (2-state: T + WS)"),
    ("mc",          "Monte Carlo UQ       (parameter sweep)"),
    ("bte_ns",      "BTE-NS spatial       (1D PDE + two-way BTE)"),
]

# Method selector options shown in the GUI
METHOD_LABELS = ["AUTO (intelligent)", "Euler", "Heun", "RK4", "RK45"]
MODE_LABELS   = ["DECOUPLED (1-state)", "COUPLED (2-state)", "SPATIAL (1D PDE)"]


class SimV0LauncherGUI:
    """Main application window for the simv0 Launcher."""

    def __init__(self, root, simv0_dir=None):
        self.root = root
        self.simv0_dir = simv0_dir or os.path.dirname(os.path.abspath(__file__))
        self._proc = None
        self._log_queue = queue.Queue()

        self._build_ui()
        self._poll_log_queue()

    # -----------------------------------------------------------------------
    #  UI construction
    # -----------------------------------------------------------------------

    def _build_ui(self):
        self.root.title("simv0 — ODE Integration Platform Launcher")
        self.root.configure(bg=PAL["navy"])
        self.root.resizable(True, True)
        self.root.minsize(780, 600)

        # ── Title bar ───────────────────────────────────────────────────────
        title_frame = tk.Frame(self.root, bg=PAL["navy_dark"], height=48)
        title_frame.pack(fill=tk.X, side=tk.TOP)
        title_frame.pack_propagate(False)

        tk.Label(
            title_frame,
            text="  simv0 — Intelligent ODE Integration Platform",
            font=("Arial", 14, "bold"),
            bg=PAL["navy_dark"],
            fg=PAL["title_fg"],
            anchor="w",
        ).pack(side=tk.LEFT, padx=8, pady=6)

        tk.Label(
            title_frame,
            text="Solar Panel Thermal Model  |  Euler · Heun · RK4 · RK45",
            font=("Arial", 8),
            bg=PAL["navy_dark"],
            fg="#AACCEE",
            anchor="e",
        ).pack(side=tk.RIGHT, padx=8, pady=6)

        # ── Main horizontal split: left panel + right console ───────────────
        main_pane = tk.PanedWindow(
            self.root,
            orient=tk.HORIZONTAL,
            bg=PAL["navy"],
            sashwidth=4,
            sashrelief=tk.RAISED,
        )
        main_pane.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

        left_frame  = tk.Frame(main_pane, bg=PAL["beige"], bd=1, relief=tk.SUNKEN)
        right_frame = tk.Frame(main_pane, bg=PAL["console_bg"], bd=1, relief=tk.SUNKEN)

        main_pane.add(left_frame,  minsize=300, width=330)
        main_pane.add(right_frame, minsize=380)

        self._build_left(left_frame)
        self._build_right(right_frame)

        # ── Status bar ───────────────────────────────────────────────────────
        self.status_var = tk.StringVar(value="Ready.")
        status_bar = tk.Label(
            self.root,
            textvariable=self.status_var,
            font=FONT_FIXED,
            bg=PAL["navy_mid"],
            fg=PAL["header"],
            anchor="w",
            relief=tk.SUNKEN,
            bd=1,
            padx=6,
        )
        status_bar.pack(fill=tk.X, side=tk.BOTTOM)

    # ── Left panel ──────────────────────────────────────────────────────────

    def _build_left(self, parent):
        canvas = tk.Canvas(parent, bg=PAL["beige"], highlightthickness=0)
        scrollbar = tk.Scrollbar(parent, orient=tk.VERTICAL, command=canvas.yview)
        canvas.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        inner = tk.Frame(canvas, bg=PAL["beige"])
        win_id = canvas.create_window((0, 0), window=inner, anchor="nw")

        def _on_configure(event):
            canvas.configure(scrollregion=canvas.bbox("all"))
            canvas.itemconfig(win_id, width=event.width)

        canvas.bind("<Configure>", _on_configure)
        inner.bind("<Configure>", lambda e: canvas.configure(
            scrollregion=canvas.bbox("all")))

        self._build_solver_params(inner)
        self._build_scenarios(inner)
        self._build_buttons(inner)

    def _section_label(self, parent, text):
        """Render a navy-background section header."""
        frame = tk.Frame(parent, bg=PAL["navy"], height=22)
        frame.pack(fill=tk.X, padx=4, pady=(8, 2))
        frame.pack_propagate(False)
        tk.Label(
            frame,
            text=f"  {text}",
            font=FONT_SECTION,
            bg=PAL["navy"],
            fg=PAL["header"],
            anchor="w",
        ).pack(fill=tk.X, padx=2)

    def _param_row(self, parent, label, var, row):
        """Render a label + Entry row on a beige grid."""
        tk.Label(
            parent,
            text=label,
            font=FONT_LABEL,
            bg=PAL["gray"],
            fg=PAL["label_fg"],
            anchor="w",
            width=20,
        ).grid(row=row, column=0, sticky="ew", padx=(4, 2), pady=2)
        entry = tk.Entry(
            parent,
            textvariable=var,
            font=FONT_FIXED,
            bg=PAL["white"],
            fg=PAL["label_fg"],
            insertbackground=PAL["label_fg"],
            relief=tk.SUNKEN,
            bd=1,
            width=14,
        )
        entry.grid(row=row, column=1, sticky="ew", padx=(2, 4), pady=2)
        return entry

    def _build_solver_params(self, parent):
        self._section_label(parent, "SOLVER PARAMETERS")

        grid = tk.Frame(parent, bg=PAL["gray"], bd=1, relief=tk.SUNKEN)
        grid.pack(fill=tk.X, padx=6, pady=2)
        grid.columnconfigure(1, weight=1)

        # StringVar holders for each parameter
        self.var_h       = tk.StringVar(value="100.0")
        self.var_t_end   = tk.StringVar(value="86400.0")
        self.var_atol    = tk.StringVar(value="0.1")
        self.var_rtol    = tk.StringVar(value="1.0e-6")
        self.var_T0      = tk.StringVar(value="298.15")
        self.var_WS0     = tk.StringVar(value="1.5")

        self._param_row(grid, "Step size h  (s):",    self.var_h,     0)
        self._param_row(grid, "End time t_end  (s):", self.var_t_end, 1)
        self._param_row(grid, "Abs. tolerance atol:", self.var_atol,  2)
        self._param_row(grid, "Rel. tolerance rtol:", self.var_rtol,  3)
        self._param_row(grid, "Init. temp T_0  (K):", self.var_T0,    4)
        self._param_row(grid, "Init. wind WS_0 (m/s):", self.var_WS0, 5)

        # Method selector
        self._section_label(parent, "INTEGRATION METHOD")
        meth_frame = tk.Frame(parent, bg=PAL["gray"], bd=1, relief=tk.SUNKEN)
        meth_frame.pack(fill=tk.X, padx=6, pady=2)
        meth_frame.columnconfigure(1, weight=1)

        tk.Label(
            meth_frame,
            text="Method:",
            font=FONT_LABEL,
            bg=PAL["gray"],
            fg=PAL["label_fg"],
            anchor="w",
            width=20,
        ).grid(row=0, column=0, sticky="w", padx=(4, 2), pady=4)

        self.var_method = tk.StringVar(value=METHOD_LABELS[0])
        method_menu = tk.OptionMenu(
            meth_frame, self.var_method, *METHOD_LABELS
        )
        method_menu.config(
            font=FONT_LABEL,
            bg=PAL["btn_bg"],
            fg=PAL["label_fg"],
            activebackground=PAL["btn_active"],
            relief=tk.RAISED,
            bd=1,
            width=18,
        )
        method_menu["menu"].config(font=FONT_LABEL, bg=PAL["gray"])
        method_menu.grid(row=0, column=1, sticky="ew", padx=(2, 4), pady=4)

        # Mode selector
        tk.Label(
            meth_frame,
            text="Sim. mode:",
            font=FONT_LABEL,
            bg=PAL["gray"],
            fg=PAL["label_fg"],
            anchor="w",
            width=20,
        ).grid(row=1, column=0, sticky="w", padx=(4, 2), pady=4)

        # Default to COUPLED (the most physically representative mode)
        DEFAULT_MODE = MODE_LABELS[1]  # "COUPLED (2-state)"
        self.var_mode = tk.StringVar(value=DEFAULT_MODE)
        mode_menu = tk.OptionMenu(
            meth_frame, self.var_mode, *MODE_LABELS
        )
        mode_menu.config(
            font=FONT_LABEL,
            bg=PAL["btn_bg"],
            fg=PAL["label_fg"],
            activebackground=PAL["btn_active"],
            relief=tk.RAISED,
            bd=1,
            width=18,
        )
        mode_menu["menu"].config(font=FONT_LABEL, bg=PAL["gray"])
        mode_menu.grid(row=1, column=1, sticky="ew", padx=(2, 4), pady=4)

    def _build_scenarios(self, parent):
        self._section_label(parent, "ODE SCENARIOS")

        scen_frame = tk.Frame(parent, bg=PAL["gray"], bd=1, relief=tk.SUNKEN)
        scen_frame.pack(fill=tk.X, padx=6, pady=2)

        self.scenario_vars = {}
        for i, (key, label) in enumerate(SCENARIOS):
            var = tk.BooleanVar(value=True)
            self.scenario_vars[key] = var
            cb = tk.Checkbutton(
                scen_frame,
                text=label,
                variable=var,
                font=FONT_FIXED,
                bg=PAL["gray"],
                fg=PAL["label_fg"],
                selectcolor=PAL["white"],
                activebackground=PAL["gray"],
                anchor="w",
            )
            cb.pack(fill=tk.X, padx=6, pady=1)

        # Select / deselect all
        btn_row = tk.Frame(scen_frame, bg=PAL["gray"])
        btn_row.pack(fill=tk.X, padx=4, pady=4)
        _make_btn(btn_row, "Select All",   self._select_all_scenarios,   width=12
                  ).pack(side=tk.LEFT, padx=2)
        _make_btn(btn_row, "Deselect All", self._deselect_all_scenarios, width=12
                  ).pack(side=tk.LEFT, padx=2)

    def _build_buttons(self, parent):
        self._section_label(parent, "ACTIONS")

        btn_frame = tk.Frame(parent, bg=PAL["beige"])
        btn_frame.pack(fill=tk.X, padx=6, pady=4)

        self.btn_run = _make_btn(btn_frame, "▶  Run Simulation", self._run, width=18)
        self.btn_run.pack(side=tk.LEFT, padx=4, pady=4)

        _make_btn(btn_frame, "■  Stop", self._stop, width=10
                  ).pack(side=tk.LEFT, padx=4, pady=4)

        btn_frame2 = tk.Frame(parent, bg=PAL["beige"])
        btn_frame2.pack(fill=tk.X, padx=6, pady=2)

        _make_btn(btn_frame2, "Build (make)", self._build_only, width=16
                  ).pack(side=tk.LEFT, padx=4, pady=4)

        _make_btn(btn_frame2, "Clear Log", self._clear_log, width=12
                  ).pack(side=tk.LEFT, padx=4, pady=4)

    # ── Right console panel ──────────────────────────────────────────────────

    def _build_right(self, parent):
        hdr = tk.Frame(parent, bg=PAL["navy_mid"], height=22)
        hdr.pack(fill=tk.X, side=tk.TOP)
        hdr.pack_propagate(False)
        tk.Label(
            hdr,
            text="  Console Output",
            font=FONT_SECTION,
            bg=PAL["navy_mid"],
            fg=PAL["header"],
            anchor="w",
        ).pack(fill=tk.X, padx=4)

        self.console = tk.Text(
            parent,
            font=FONT_CONSOLE,
            bg=PAL["console_bg"],
            fg=PAL["phosphor"],
            insertbackground=PAL["phosphor"],
            relief=tk.FLAT,
            bd=0,
            wrap=tk.WORD,
            state=tk.DISABLED,
            takefocus=0,
        )

        vsb = tk.Scrollbar(parent, orient=tk.VERTICAL,
                           command=self.console.yview)
        self.console.configure(yscrollcommand=vsb.set)

        vsb.pack(side=tk.RIGHT, fill=tk.Y)
        self.console.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Colour tags
        # Bold variant derived from FONT_CONSOLE (family + size + bold)
        font_console_bold = (FONT_CONSOLE[0], FONT_CONSOLE[1], "bold")
        self.console.tag_configure("header",  foreground=PAL["header"],  font=font_console_bold)
        self.console.tag_configure("success", foreground=PAL["success"])
        self.console.tag_configure("error",   foreground=PAL["error"])
        self.console.tag_configure("warning", foreground=PAL["warning"])
        self.console.tag_configure("dim",     foreground="#556B4F")
        self.console.tag_configure("normal",  foreground=PAL["phosphor"])

    # -----------------------------------------------------------------------
    #  Console helpers
    # -----------------------------------------------------------------------

    def _log(self, text, tag="normal"):
        """Append a line to the console (thread-safe via queue)."""
        self._log_queue.put((text, tag))

    def _flush_log(self, text, tag="normal"):
        """Write directly to console — must be called from the main thread."""
        self.console.configure(state=tk.NORMAL)
        self.console.insert(tk.END, text + "\n", tag)
        self.console.see(tk.END)
        self.console.configure(state=tk.DISABLED)

    def _poll_log_queue(self):
        """Drain the thread-safe log queue and write to console."""
        try:
            while True:
                text, tag = self._log_queue.get_nowait()
                self._flush_log(text, tag)
        except queue.Empty:
            pass
        self.root.after(100, self._poll_log_queue)

    def _classify_line(self, line):
        """Return a colour tag for a console output line."""
        l = line.lower()
        if any(w in l for w in ("error", "fatal", "failed", "abort")):
            return "error"
        if any(w in l for w in ("warning", "warn")):
            return "warning"
        if any(w in l for w in ("✓", "done", "ok", "success", "complete",
                                 "written", "output")):
            return "success"
        if line.startswith("  --") or line.startswith("===") or line.startswith("---"):
            return "header"
        if line.startswith("#") or "----" in line:
            return "dim"
        return "normal"

    # -----------------------------------------------------------------------
    #  Action callbacks
    # -----------------------------------------------------------------------

    def _select_all_scenarios(self):
        for v in self.scenario_vars.values():
            v.set(True)

    def _deselect_all_scenarios(self):
        for v in self.scenario_vars.values():
            v.set(False)

    def _clear_log(self):
        self.console.configure(state=tk.NORMAL)
        self.console.delete("1.0", tk.END)
        self.console.configure(state=tk.DISABLED)
        self._set_status("Log cleared.")

    def _set_status(self, msg):
        self.status_var.set(msg)
        self.root.update_idletasks()

    # ── Build ────────────────────────────────────────────────────────────────

    def _build_only(self):
        """Run `make` in the simv0 directory."""
        self._log("=" * 60, "header")
        self._log("  Building simv0 …", "header")
        self._log("=" * 60, "header")
        self._set_status("Building …")
        self.btn_run.configure(state=tk.DISABLED)
        threading.Thread(target=self._build_thread, daemon=True).start()

    def _build_thread(self):
        ok = self._run_make(["make"])
        if ok:
            self._log("Build completed successfully.", "success")
            self._set_status("Build OK.")
        else:
            self._log("Build failed — check output above.", "error")
            self._set_status("Build FAILED.")
        self.root.after(0, lambda: self.btn_run.configure(state=tk.NORMAL))

    # ── Run ──────────────────────────────────────────────────────────────────

    def _run(self):
        selected = [k for k, v in self.scenario_vars.items() if v.get()]
        if not selected:
            messagebox.showwarning(
                "No scenarios selected",
                "Please select at least one ODE scenario to run.",
            )
            return

        self._clear_log()
        self._log("=" * 60, "header")
        self._log("  simv0 — ODE Integration Platform", "header")
        self._log("=" * 60, "header")
        self._log(f"  Scenarios : {', '.join(selected)}", "normal")
        self._log(f"  Method    : {self.var_method.get()}", "normal")
        self._log(f"  Mode      : {self.var_mode.get()}", "normal")
        self._log(f"  h         : {self.var_h.get()} s", "normal")
        self._log(f"  t_end     : {self.var_t_end.get()} s", "normal")
        self._log(f"  atol      : {self.var_atol.get()}", "normal")
        self._log(f"  rtol      : {self.var_rtol.get()}", "normal")
        self._log(f"  T_0       : {self.var_T0.get()} K", "normal")
        self._log(f"  WS_0      : {self.var_WS0.get()} m/s", "normal")
        self._log("-" * 60, "dim")

        self.btn_run.configure(state=tk.DISABLED)
        self._set_status("Running …")
        threading.Thread(target=self._run_thread, daemon=True).start()

    def _run_thread(self):
        # Step 1: build if binary missing
        binary = os.path.join(self.simv0_dir, "bin", "simv0")
        if not os.path.isfile(binary):
            self._log("Binary not found — building first …", "warning")
            ok = self._run_make(["make"])
            if not ok:
                self._log("Build failed.  Cannot run.", "error")
                self._set_status("Build FAILED.")
                self.root.after(0, lambda: self.btn_run.configure(state=tk.NORMAL))
                return

        # Step 2: run
        self._log("", "dim")
        self._log("  Running simulation …", "header")
        self._log("-" * 60, "dim")
        ok = self._run_make(["make", "run"])
        self._log("-" * 60, "dim")

        if ok:
            self._log("Simulation completed successfully.", "success")
            self._set_status("Simulation done.")
            self._parse_and_display_results()
        else:
            self._log("Simulation exited with errors.", "error")
            self._set_status("Simulation FAILED.")

        self.root.after(0, lambda: self.btn_run.configure(state=tk.NORMAL))

    def _stop(self):
        if self._proc and self._proc.poll() is None:
            self._proc.terminate()
            self._log("Process terminated by user.", "warning")
            self._set_status("Stopped.")

    # ── Make runner ──────────────────────────────────────────────────────────

    def _run_make(self, cmd):
        """Run a make command in simv0_dir, streaming output to the console."""
        try:
            self._proc = subprocess.Popen(
                cmd,
                cwd=self.simv0_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )
            for line in self._proc.stdout:
                line = line.rstrip("\n")
                self._log(line, self._classify_line(line))
            self._proc.wait()
            return self._proc.returncode == 0
        except FileNotFoundError:
            self._log(f"Command not found: {cmd[0]}", "error")
            return False

    # ── Results parsing ──────────────────────────────────────────────────────

    def _parse_and_display_results(self):
        """Scan output .dat files and display a brief metrics summary."""
        data_dir = os.path.join(self.simv0_dir, "data")
        if not os.path.isdir(data_dir):
            return

        selected = {k for k, v in self.scenario_vars.items() if v.get()}
        results = []

        # Map scenario names to relevant output files
        file_map = {
            "euler":       "simv0_decoupled.dat",
            "rk4":         "simv0_decoupled.dat",
            "atmospheric": "simv0_coupled.dat",
            "diurnal":     "simv0_coupled.dat",
            "coupled":     "simv0_coupled.dat",
            "mc":          "simv0_comparison.dat",
            "bte_ns":      "simv0_spatial.dat",
        }

        seen_files = set()
        for key in selected:
            fname = file_map.get(key)
            if fname and fname not in seen_files:
                seen_files.add(fname)
                path = os.path.join(data_dir, fname)
                metrics = _parse_dat_file(path)
                if metrics:
                    results.append((fname, metrics))

        if not results:
            return

        self._log("", "dim")
        self._log("=" * 60, "header")
        self._log("  KEY METRICS FROM OUTPUT FILES", "header")
        self._log("=" * 60, "header")

        for fname, metrics in results:
            self._log(f"\n  File: {fname}", "warning")
            for k, v in metrics.items():
                self._log(f"    {k:<30s} {v}", "normal")

        self._log("", "dim")


# ---------------------------------------------------------------------------
#  .dat file parser
# ---------------------------------------------------------------------------

def _parse_dat_file(path):
    """
    Return a dict of key metrics extracted from a simv0 .dat output file.
    Reads the last non-comment, non-empty row plus the header to map columns.
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
                # Capture the most recent comment-header line
                # e.g. "# t(s)   T(K)   WS(m/s) …"
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

    metrics = {}
    if header and len(header) == len(last_row):
        for name, val in zip(header, last_row):
            metrics[name] = f"{val:.4g}"
    else:
        # Fall back: label columns numerically
        for i, val in enumerate(last_row):
            metrics[f"col_{i}"] = f"{val:.4g}"

    return metrics


# ---------------------------------------------------------------------------
#  Entry point
# ---------------------------------------------------------------------------

def launch(simv0_dir=None):
    """Create the Tk root window and start the event loop."""
    root = tk.Tk()
    app = SimV0LauncherGUI(root, simv0_dir=simv0_dir)
    root.mainloop()


if __name__ == "__main__":
    launch()
