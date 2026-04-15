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
from tkinter import ttk, messagebox, filedialog

import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib import image as mpimg

from plotting_service import (
    PRESET_TEMPLATES,
    PlotRequest,
    PlottingService,
    SeriesSpec,
    apply_preset_to_request,
    discover_dat_files,
    format_stats_table,
)

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
        self._plot_service = PlottingService(
            allowed_data_dir=os.path.join(self.simv0_dir, "data")
        )
        self._pb_series = []
        self._pb_last_result = None
        self._pb_current_request = None

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

        def _on_mousewheel(event):
            if getattr(event, "num", None) == 4:
                delta = -1
            elif getattr(event, "num", None) == 5:
                delta = 1
            elif getattr(event, "delta", 0):
                delta = -1 if event.delta > 0 else 1
            else:
                return
            canvas.yview_scroll(delta, "units")

        def _bind_mousewheel(_event=None):
            canvas.bind_all("<MouseWheel>", _on_mousewheel)
            canvas.bind_all("<Button-4>", _on_mousewheel)
            canvas.bind_all("<Button-5>", _on_mousewheel)

        def _unbind_mousewheel(_event=None):
            canvas.unbind_all("<MouseWheel>")
            canvas.unbind_all("<Button-4>")
            canvas.unbind_all("<Button-5>")

        canvas.bind("<Configure>", _on_configure)
        inner.bind("<Configure>", lambda e: canvas.configure(
            scrollregion=canvas.bbox("all")))
        canvas.bind("<Enter>", _bind_mousewheel)
        canvas.bind("<Leave>", _unbind_mousewheel)
        inner.bind("<Enter>", _bind_mousewheel)
        inner.bind("<Leave>", _unbind_mousewheel)

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

    # ── Right panel: tabbed notebook ─────────────────────────────────────────

    def _build_right(self, parent):
        # Apply retro styling to the ttk.Notebook
        style = ttk.Style(parent)
        style.theme_use("default")
        style.configure(
            "Retro.TNotebook",
            background=PAL["navy"],
            borderwidth=0,
            tabmargins=[2, 2, 0, 0],
        )
        style.configure(
            "Retro.TNotebook.Tab",
            background=PAL["navy_mid"],
            foreground=PAL["header"],
            font=FONT_SECTION,
            padding=[8, 3],
            borderwidth=1,
        )
        style.map(
            "Retro.TNotebook.Tab",
            background=[("selected", PAL["navy_dark"]), ("active", PAL["navy"])],
            foreground=[("selected", PAL["title_fg"]), ("active", PAL["header"])],
        )

        self._notebook = ttk.Notebook(parent, style="Retro.TNotebook")
        self._notebook.pack(fill=tk.BOTH, expand=True)

        # Tab 1 – Data Viewer
        data_frame = tk.Frame(self._notebook, bg=PAL["console_bg"])
        self._notebook.add(data_frame, text="  Data Viewer  ")
        self._build_data_viewer_tab(data_frame)

        # Tab 2 – Plot Builder
        plot_frame = tk.Frame(self._notebook, bg=PAL["console_bg"])
        self._notebook.add(plot_frame, text="  Plot Builder  ")
        self._build_plot_builder_tab(plot_frame)

        # Tab 3 – Console Output  (existing behaviour)
        console_frame = tk.Frame(self._notebook, bg=PAL["console_bg"])
        self._notebook.add(console_frame, text="  Console  ")
        self._build_console_tab(console_frame)

    # ── Console tab ──────────────────────────────────────────────────────────

    def _build_console_tab(self, parent):
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
        vsb = tk.Scrollbar(parent, orient=tk.VERTICAL, command=self.console.yview)
        self.console.configure(yscrollcommand=vsb.set)
        vsb.pack(side=tk.RIGHT, fill=tk.Y)
        self.console.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        font_console_bold = (FONT_CONSOLE[0], FONT_CONSOLE[1], "bold")
        self.console.tag_configure("header",  foreground=PAL["header"],  font=font_console_bold)
        self.console.tag_configure("success", foreground=PAL["success"])
        self.console.tag_configure("error",   foreground=PAL["error"])
        self.console.tag_configure("warning", foreground=PAL["warning"])
        self.console.tag_configure("dim",     foreground="#556B4F")
        self.console.tag_configure("normal",  foreground=PAL["phosphor"])

    # ── Data Viewer tab ──────────────────────────────────────────────────────

    # Fallback output files available for viewing
    _DAT_FILES = [
        "simv0_decoupled.dat",
        "simv0_coupled.dat",
        "simv0_comparison.dat",
        "simv0_spatial.dat",
        "simv0_diagnostic.dat",
    ]

    def _list_dat_files(self):
        data_dir = os.path.join(self.simv0_dir, "data")
        files = discover_dat_files(data_dir)
        return files or list(self._DAT_FILES)

    def _build_data_viewer_tab(self, parent):
        # ── Toolbar ──────────────────────────────────────────────────────────
        toolbar = tk.Frame(parent, bg=PAL["navy_mid"], height=28)
        toolbar.pack(fill=tk.X, side=tk.TOP)
        toolbar.pack_propagate(False)

        tk.Label(
            toolbar,
            text="  File:",
            font=FONT_SECTION,
            bg=PAL["navy_mid"],
            fg=PAL["header"],
        ).pack(side=tk.LEFT, padx=(4, 0), pady=4)

        self._dv_files = self._list_dat_files()
        self._dv_file_var = tk.StringVar(value=self._dv_files[0])
        dv_menu = tk.OptionMenu(
            toolbar, self._dv_file_var, *self._dv_files,
            command=lambda _: self._dv_load(),
        )
        dv_menu.config(
            font=FONT_LABEL,
            bg=PAL["btn_bg"],
            fg=PAL["label_fg"],
            activebackground=PAL["btn_active"],
            relief=tk.RAISED,
            bd=1,
        )
        dv_menu["menu"].config(font=FONT_LABEL, bg=PAL["gray"])
        dv_menu.pack(side=tk.LEFT, padx=4, pady=3)
        self._dv_file_menu = dv_menu

        _make_btn(toolbar, "↺  Refresh", self._dv_refresh_and_load, width=10
                  ).pack(side=tk.LEFT, padx=4, pady=3)

        # Search / filter
        tk.Label(
            toolbar,
            text="Filter:",
            font=FONT_SECTION,
            bg=PAL["navy_mid"],
            fg=PAL["header"],
        ).pack(side=tk.LEFT, padx=(12, 0), pady=4)
        self._dv_filter_var = tk.StringVar()
        self._dv_filter_var.trace_add("write", lambda *_: self._dv_apply_filter())
        filter_entry = tk.Entry(
            toolbar,
            textvariable=self._dv_filter_var,
            font=FONT_FIXED,
            bg=PAL["white"],
            fg=PAL["label_fg"],
            insertbackground=PAL["label_fg"],
            relief=tk.SUNKEN,
            bd=1,
            width=16,
        )
        filter_entry.pack(side=tk.LEFT, padx=4, pady=4)

        # ── Treeview ─────────────────────────────────────────────────────────
        tv_frame = tk.Frame(parent, bg=PAL["console_bg"])
        tv_frame.pack(fill=tk.BOTH, expand=True)

        style = ttk.Style(parent)
        style.configure(
            "Retro.Treeview",
            background=PAL["console_bg"],
            foreground=PAL["phosphor"],
            fieldbackground=PAL["console_bg"],
            font=FONT_CONSOLE,
            rowheight=16,
        )
        style.configure(
            "Retro.Treeview.Heading",
            background=PAL["navy_mid"],
            foreground=PAL["header"],
            font=FONT_SECTION,
            relief="flat",
        )
        style.map(
            "Retro.Treeview",
            background=[("selected", PAL["navy"])],
            foreground=[("selected", PAL["title_fg"])],
        )

        self._dv_tree = ttk.Treeview(
            tv_frame,
            style="Retro.Treeview",
            show="headings",
            selectmode="browse",
        )

        vsb = tk.Scrollbar(tv_frame, orient=tk.VERTICAL,
                           command=self._dv_tree.yview)
        hsb = tk.Scrollbar(tv_frame, orient=tk.HORIZONTAL,
                           command=self._dv_tree.xview)
        self._dv_tree.configure(
            yscrollcommand=vsb.set,
            xscrollcommand=hsb.set,
        )

        hsb.pack(side=tk.BOTTOM, fill=tk.X)
        vsb.pack(side=tk.RIGHT,  fill=tk.Y)
        self._dv_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Status label for data viewer
        self._dv_status_var = tk.StringVar(value="No data loaded.")
        tk.Label(
            parent,
            textvariable=self._dv_status_var,
            font=FONT_FIXED,
            bg=PAL["navy_dark"],
            fg=PAL["header"],
            anchor="w",
            padx=6,
        ).pack(fill=tk.X, side=tk.BOTTOM)

        # Internal cache of all rows for filtering
        self._dv_all_rows = []

    def _dv_load(self):
        """Parse the selected .dat file and populate the Treeview."""
        fname = self._dv_file_var.get()
        data_dir = os.path.join(self.simv0_dir, "data")
        path = os.path.join(data_dir, fname)

        headers, rows = _parse_dat_full(path)

        # Clear existing columns/rows
        self._dv_tree.delete(*self._dv_tree.get_children())
        self._dv_tree["columns"] = headers if headers else []

        if not headers:
            self._dv_status_var.set(
                f"File not found or empty: {fname}" if not os.path.isfile(path)
                else f"No data in {fname}"
            )
            self._dv_all_rows = []
            return

        # Configure columns
        for col in headers:
            self._dv_tree.heading(col, text=col, anchor="e")
            # Width heuristic: longer headers get more space
            width = max(90, len(col) * 9 + 16)
            self._dv_tree.column(col, width=width, anchor="e", minwidth=60, stretch=False)

        self._dv_all_rows = rows
        self._dv_filter_var.set("")  # reset filter
        self._dv_populate_rows(rows)
        self._dv_status_var.set(
            f"Loaded {len(rows)} rows × {len(headers)} columns  ·  {fname}"
        )

    def _dv_refresh_files(self):
        files = self._list_dat_files()
        # Rebuild OptionMenu entries in-place
        if hasattr(self, "_dv_file_var") and hasattr(self, "_dv_file_menu"):
            menu = self._dv_file_menu["menu"]
            menu.delete(0, tk.END)
            def _select_file(v):
                self._dv_file_var.set(v)
                self._dv_load()
            for f in files:
                menu.add_command(label=f, command=lambda v=f: _select_file(v))
            if self._dv_file_var.get() not in files:
                self._dv_file_var.set(files[0])
        self._dv_files = files

    def _dv_refresh_and_load(self):
        self._dv_refresh_files()
        self._dv_load()

    def _dv_populate_rows(self, rows):
        """Insert a list of string-tuples into the Treeview."""
        self._dv_tree.delete(*self._dv_tree.get_children())
        for row in rows:
            self._dv_tree.insert("", tk.END, values=row)

    def _dv_apply_filter(self):
        """Filter rows to those containing the search string in any column."""
        term = self._dv_filter_var.get().strip().lower()
        if not term:
            self._dv_populate_rows(self._dv_all_rows)
            return
        filtered = [r for r in self._dv_all_rows
                    if any(term in str(v).lower() for v in r)]
        self._dv_populate_rows(filtered)
        self._dv_status_var.set(
            f"Showing {len(filtered)} / {len(self._dv_all_rows)} rows  (filter: '{term}')"
        )

    # ── Plot Builder tab ─────────────────────────────────────────────────────

    _PLOT_TYPES  = ["line", "scatter", "bar"]
    _PLOT_MODES  = ["2d", "3d"]
    _BACKENDS    = ["auto", "gnuplot", "matplotlib"]
    _COLORS      = ["#33FF00", "#66CCFF", "#FF4444", "#FFCC00", "#FF88FF",
                    "#00FFCC", "#FFFFFF"]
    _COLOR_NAMES = ["Green", "Cyan", "Red", "Yellow", "Magenta", "Teal", "White"]
    _MARKERS     = ["None", "o", "s", "^", "D", "x", "+"]
    _LINE_WIDTHS = ["0.5", "1.0", "1.5", "2.0", "2.5"]

    def _build_plot_builder_tab(self, parent):
        # Split: controls on left, canvas on right
        pane = tk.PanedWindow(
            parent,
            orient=tk.HORIZONTAL,
            bg=PAL["navy"],
            sashwidth=4,
            sashrelief=tk.RAISED,
        )
        pane.pack(fill=tk.BOTH, expand=True)

        ctrl_frame   = tk.Frame(pane, bg=PAL["gray"], bd=0)
        canvas_frame = tk.Frame(pane, bg=PAL["console_bg"], bd=0)
        pane.add(ctrl_frame,   minsize=200, width=220)
        pane.add(canvas_frame, minsize=300)

        # ── Controls ─────────────────────────────────────────────────────────
        self._pb_build_controls(ctrl_frame)

        # ── Embedded canvas ───────────────────────────────────────────────────
        self._pb_fig = Figure(figsize=(5, 4), dpi=96,
                              facecolor=PAL["console_bg"])
        self._pb_ax  = self._pb_fig.add_subplot(111)
        _style_axes(self._pb_ax)

        self._pb_canvas = FigureCanvasTkAgg(self._pb_fig, master=canvas_frame)
        self._pb_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self._pb_canvas.draw()

        # Save / clear buttons below the canvas
        btn_bar = tk.Frame(canvas_frame, bg=PAL["navy_mid"], height=30)
        btn_bar.pack(fill=tk.X, side=tk.BOTTOM)
        btn_bar.pack_propagate(False)
        _make_btn(btn_bar, "Save Plot", self._pb_save, width=12
                  ).pack(side=tk.LEFT, padx=6, pady=3)
        _make_btn(btn_bar, "Save Stats", self._pb_save_stats, width=12
                  ).pack(side=tk.LEFT, padx=4, pady=3)
        _make_btn(btn_bar, "Clear",     self._pb_clear, width=8
                  ).pack(side=tk.LEFT, padx=4, pady=3)
        self._pb_save_status = tk.StringVar(value="")
        tk.Label(
            btn_bar,
            textvariable=self._pb_save_status,
            font=FONT_FIXED,
            bg=PAL["navy_mid"],
            fg=PAL["success"],
            anchor="w",
        ).pack(side=tk.LEFT, padx=6)

        stats_frame = tk.Frame(canvas_frame, bg=PAL["console_bg"])
        stats_frame.pack(fill=tk.X, side=tk.BOTTOM)
        tk.Label(
            stats_frame,
            text="Statistics:",
            font=FONT_SECTION,
            bg=PAL["console_bg"],
            fg=PAL["header"],
            anchor="w",
        ).pack(fill=tk.X, padx=6, pady=(4, 2))
        self._pb_stats_text = tk.Text(
            stats_frame,
            height=8,
            font=FONT_FIXED,
            bg=PAL["console_bg"],
            fg=PAL["phosphor"],
            insertbackground=PAL["phosphor"],
            relief=tk.SUNKEN,
            bd=1,
            wrap=tk.NONE,
        )
        self._pb_stats_text.pack(fill=tk.X, padx=6, pady=(0, 4))
        self._pb_set_stats_text("No statistics yet.")

    def _pb_build_controls(self, parent):
        """Build the left-side control panel for the Plot Builder."""

        def _sec(text):
            fr = tk.Frame(parent, bg=PAL["navy"], height=20)
            fr.pack(fill=tk.X, padx=4, pady=(6, 1))
            fr.pack_propagate(False)
            tk.Label(fr, text=f"  {text}", font=FONT_SECTION,
                     bg=PAL["navy"], fg=PAL["header"], anchor="w"
                     ).pack(fill=tk.X, padx=2)

        def _row(parent_fr, label, widget_factory, pady=2):
            tk.Label(parent_fr, text=label, font=FONT_LABEL,
                     bg=PAL["gray"], fg=PAL["label_fg"], anchor="w"
                     ).pack(fill=tk.X, padx=6, pady=(pady, 0))
            w = widget_factory(parent_fr)
            w.pack(fill=tk.X, padx=6, pady=(0, pady))
            return w

        def _omenu(parent_fr, var, choices, cmd=None):
            m = tk.OptionMenu(parent_fr, var, *choices, command=cmd)
            m.config(font=FONT_LABEL, bg=PAL["btn_bg"], fg=PAL["label_fg"],
                     activebackground=PAL["btn_active"], relief=tk.RAISED, bd=1)
            m["menu"].config(font=FONT_LABEL, bg=PAL["gray"])
            return m

        def _plain_omenu(parent_fr, var, choices, cmd=None):
            kw = {"command": cmd} if cmd else {}
            m = tk.OptionMenu(parent_fr, var, *choices, **kw)
            m.config(font=FONT_LABEL, bg=PAL["btn_bg"], fg=PAL["label_fg"],
                     activebackground=PAL["btn_active"], relief=tk.RAISED, bd=1)
            m["menu"].config(font=FONT_LABEL, bg=PAL["gray"])
            return m

        # Data files
        _sec("DATA SOURCE")
        self._pb_files = self._list_dat_files()
        files_frame = tk.Frame(parent, bg=PAL["gray"])
        files_frame.pack(fill=tk.BOTH, padx=6, pady=2)
        self._pb_file_list = tk.Listbox(
            files_frame,
            selectmode=tk.EXTENDED,
            height=4,
            exportselection=False,
            font=FONT_FIXED,
            bg=PAL["white"],
            fg=PAL["label_fg"],
        )
        for f in self._pb_files:
            self._pb_file_list.insert(tk.END, f)
        self._pb_file_list.pack(fill=tk.BOTH, expand=True)
        self._pb_file_list.bind("<<ListboxSelect>>", lambda *_: self._pb_update_col_menus())
        if self._pb_files:
            self._pb_file_list.selection_set(0)
        _make_btn(parent, "↺ Refresh files", self._pb_refresh_files, width=18
                  ).pack(padx=8, pady=(0, 4), fill=tk.X)

        # Axes
        _sec("AXES")
        _placeholder = ["(load file first)"]
        self._pb_xcol_var = tk.StringVar(value=_placeholder[0])
        self._pb_ycol_var = tk.StringVar(value=_placeholder[0])
        self._pb_yerr_var = tk.StringVar(value="None")
        self._pb_zcol_var = tk.StringVar(value="None")

        self._pb_xcol_menu = _row(parent, "X axis:", lambda p: _plain_omenu(p, self._pb_xcol_var, _placeholder))
        self._pb_ycol_menu = _row(parent, "Y axis (series):", lambda p: _plain_omenu(p, self._pb_ycol_var, _placeholder))
        self._pb_yerr_menu = _row(parent, "Y error (opt.):", lambda p: _plain_omenu(p, self._pb_yerr_var, ["None"] + _placeholder))
        self._pb_zcol_menu = _row(parent, "Z axis (3D):", lambda p: _plain_omenu(p, self._pb_zcol_var, ["None"] + _placeholder))

        self._pb_series_label_var = tk.StringVar(value="")
        _row(parent, "Series label (opt.):", lambda p: tk.Entry(
            p, textvariable=self._pb_series_label_var, font=FONT_FIXED,
            bg=PAL["white"], fg=PAL["label_fg"], relief=tk.SUNKEN, bd=1
        ))

        # Plot type
        _sec("STYLE")
        self._pb_mode_var  = tk.StringVar(value=self._PLOT_MODES[0])
        _row(parent, "Plot mode:", lambda p: _plain_omenu(p, self._pb_mode_var, self._PLOT_MODES))

        self._pb_type_var  = tk.StringVar(value=self._PLOT_TYPES[0])
        _row(parent, "Plot type:", lambda p: _plain_omenu(p, self._pb_type_var, self._PLOT_TYPES))

        self._pb_backend_var = tk.StringVar(value=self._BACKENDS[0])
        _row(parent, "Backend:", lambda p: _plain_omenu(p, self._pb_backend_var, self._BACKENDS))

        self._pb_color_var = tk.StringVar(value=self._COLOR_NAMES[0])
        _row(parent, "Color:", lambda p: _plain_omenu(p, self._pb_color_var, self._COLOR_NAMES))

        self._pb_marker_var = tk.StringVar(value=self._MARKERS[0])
        _row(parent, "Marker:", lambda p: _plain_omenu(p, self._pb_marker_var, self._MARKERS))

        self._pb_lw_var = tk.StringVar(value="1.5")
        _row(parent, "Line width:", lambda p: _plain_omenu(p, self._pb_lw_var, self._LINE_WIDTHS))

        self._pb_title_var = tk.StringVar(value="")
        _row(parent, "Title:", lambda p: tk.Entry(
            p, textvariable=self._pb_title_var, font=FONT_FIXED,
            bg=PAL["white"], fg=PAL["label_fg"], relief=tk.SUNKEN, bd=1
        ))

        self._pb_include_stats_var = tk.BooleanVar(value=True)
        tk.Checkbutton(
            parent,
            text="Compute statistics",
            variable=self._pb_include_stats_var,
            font=FONT_LABEL,
            bg=PAL["gray"],
            fg=PAL["label_fg"],
            selectcolor=PAL["white"],
            activebackground=PAL["gray"],
            anchor="w",
        ).pack(fill=tk.X, padx=6, pady=(2, 1))

        # Presets
        _sec("PRESETS")
        preset_names = list(PRESET_TEMPLATES.keys())
        self._pb_preset_var = tk.StringVar(value=preset_names[0] if preset_names else "None")
        _row(parent, "Preset:", lambda p: _omenu(p, self._pb_preset_var, preset_names or ["None"]))
        _make_btn(parent, "Apply preset", self._pb_apply_preset, width=18
                  ).pack(padx=8, pady=2, fill=tk.X)

        # Series list
        _sec("SERIES")
        self._pb_series_list = tk.Listbox(
            parent,
            selectmode=tk.BROWSE,
            height=5,
            exportselection=False,
            font=FONT_FIXED,
            bg=PAL["white"],
            fg=PAL["label_fg"],
        )
        self._pb_series_list.pack(fill=tk.BOTH, padx=6, pady=2)
        sbtn = tk.Frame(parent, bg=PAL["gray"])
        sbtn.pack(fill=tk.X, padx=6, pady=(0, 2))
        _make_btn(sbtn, "+ Add series", self._pb_add_series, width=12).pack(side=tk.LEFT, padx=2)
        _make_btn(sbtn, "- Remove", self._pb_remove_series, width=10).pack(side=tk.LEFT, padx=2)

        # Build button
        tk.Frame(parent, bg=PAL["gray"], height=8).pack()
        _make_btn(parent, "▶  Build Plot", self._pb_build, width=18
                  ).pack(padx=8, pady=6, fill=tk.X)

        # Seed first series for convenience
        self._pb_update_col_menus()

    def _pb_update_col_menus(self):
        """Reload column names from the first selected file into axis menus."""
        selected_files = self._pb_selected_files()
        if not selected_files:
            cols = ["(no data)"]
            self._pb_set_menu_choices(self._pb_xcol_var, self._pb_xcol_menu, cols)
            self._pb_set_menu_choices(self._pb_ycol_var, self._pb_ycol_menu, cols)
            self._pb_set_menu_choices(self._pb_yerr_var, self._pb_yerr_menu, ["None"] + cols)
            self._pb_set_menu_choices(self._pb_zcol_var, self._pb_zcol_menu, ["None"] + cols)
            return
        path = os.path.join(self.simv0_dir, "data", selected_files[0])
        pdata = self._plot_service.parse_file(path, max_rows=0)
        headers = pdata.headers

        cols = headers if headers else ["(no data)"]
        cols_with_none = ["None"] + cols

        self._pb_set_menu_choices(self._pb_xcol_var, self._pb_xcol_menu, cols)
        self._pb_set_menu_choices(self._pb_ycol_var, self._pb_ycol_menu, cols)
        self._pb_set_menu_choices(self._pb_yerr_var, self._pb_yerr_menu, cols_with_none)
        self._pb_set_menu_choices(self._pb_zcol_var, self._pb_zcol_menu, cols_with_none)

    def _pb_set_menu_choices(self, var, menu_widget, choices):
        menu = menu_widget["menu"]
        menu.delete(0, tk.END)
        for c in choices:
            menu.add_command(label=c, command=lambda v=c, sv=var: sv.set(v))
        if choices and var.get() not in choices:
            var.set(choices[0])

    def _pb_selected_files(self):
        return [self._pb_file_list.get(i) for i in self._pb_file_list.curselection()]

    def _pb_refresh_files(self):
        files = self._list_dat_files()
        self._pb_file_list.delete(0, tk.END)
        for f in files:
            self._pb_file_list.insert(tk.END, f)
        if files:
            self._pb_file_list.selection_set(0)
        self._plot_service.clear_cache()
        self._pb_update_col_menus()

    def _pb_add_series(self):
        files = self._pb_selected_files()
        if not files:
            messagebox.showwarning("Plot Builder", "Select at least one data file.")
            return
        ycol = self._pb_ycol_var.get()
        if not ycol or ycol in ("(no data)", "(load file first)"):
            messagebox.showwarning("Plot Builder", "Select a valid Y column.")
            return

        marker = self._pb_marker_var.get()
        marker = None if marker == "None" else marker
        y_error = self._pb_yerr_var.get()
        y_error = None if y_error in ("None", "(no data)", "(load file first)") else y_error
        zcol = self._pb_zcol_var.get()
        zcol = None if zcol in ("None", "(no data)", "(load file first)") else zcol

        try:
            line_width = float(self._pb_lw_var.get())
        except ValueError:
            messagebox.showerror("Plot Builder", "Line width must be numeric.")
            return

        for fname in files:
            label = self._pb_series_label_var.get().strip() or f"{fname}:{ycol}"
            if len(files) > 1 and self._pb_series_label_var.get().strip():
                label = f"{self._pb_series_label_var.get().strip()} [{fname}]"
            spec = SeriesSpec(
                file_path=os.path.join(self.simv0_dir, "data", fname),
                y_column=ycol,
                label=label,
                color=self._COLORS[self._COLOR_NAMES.index(self._pb_color_var.get())],
                marker=marker,
                line_width=line_width,
                yerr_column=y_error,
                z_column=zcol,
            )
            self._pb_series.append(spec)
            self._pb_series_list.insert(tk.END, self._pb_series_desc(spec))

    def _pb_remove_series(self):
        sel = self._pb_series_list.curselection()
        if not sel:
            return
        idx = sel[0]
        del self._pb_series[idx]
        self._pb_series_list.delete(idx)

    def _pb_series_desc(self, spec):
        eb = f" ±{spec.yerr_column}" if spec.yerr_column else ""
        z = f" z={spec.z_column}" if spec.z_column else ""
        return f"{os.path.basename(spec.file_path)}:{spec.y_column}{eb}{z}"

    def _pb_apply_preset(self):
        if not self._pb_series:
            messagebox.showinfo("Plot Builder", "Add at least one series before applying presets.")
            return
        preset_name = self._pb_preset_var.get()
        if preset_name not in PRESET_TEMPLATES:
            messagebox.showwarning("Plot Builder", f"Unknown preset: {preset_name}")
            return
        req = self._pb_make_request(output_path=None)
        req = apply_preset_to_request(req, preset_name)
        self._pb_mode_var.set(req.plot_mode)
        self._pb_type_var.set(req.plot_type)
        self._pb_title_var.set(req.title)
        if req.x_column:
            self._pb_xcol_var.set(req.x_column)
        self._set_status(f"Preset applied: {preset_name}")

    def _pb_make_request(self, output_path=None):
        title = self._pb_title_var.get().strip() or "simv0 plot"
        return PlotRequest(
            plot_mode=self._pb_mode_var.get(),
            plot_type=self._pb_type_var.get(),
            backend=self._pb_backend_var.get(),
            x_column=self._pb_xcol_var.get(),
            series=list(self._pb_series),
            title=title,
            x_label=self._pb_xcol_var.get(),
            y_label=self._pb_series[0].y_column if self._pb_series else "",
            z_label=self._pb_zcol_var.get() if self._pb_mode_var.get() == "3d" else None,
            output_path=output_path,
            output_format="png",
            include_stats=self._pb_include_stats_var.get(),
        )

    def _pb_build(self):
        """Build the plot from selected options via plotting service."""
        if not self._pb_series:
            messagebox.showwarning("Plot Builder", "Add at least one series before building.")
            return
        req = self._pb_make_request(output_path=None)
        if req.plot_mode == "3d" and any(not s.z_column for s in req.series):
            messagebox.showerror("Plot Builder", "3D mode requires a Z column for all series.")
            return

        result = self._plot_service.render(req)
        self._pb_current_request = req
        self._pb_last_result = result
        if not result.ok:
            messagebox.showerror("Plot Builder", result.error)
            self._pb_set_stats_text(format_stats_table(result.stats))
            self._set_status("Plot build failed.")
            return

        # Display rendered image
        self._pb_fig.clear()
        ax = self._pb_fig.add_subplot(111)
        _style_axes(ax)
        try:
            img = mpimg.imread(result.output_path)
            ax.imshow(img)
            ax.set_axis_off()
        except Exception:
            ax.text(0.5, 0.5, "Plot rendered.\n(Preview unavailable)", ha="center", va="center", color=PAL["phosphor"])
        self._pb_fig.tight_layout(pad=0.2)
        self._pb_canvas.draw()
        self._pb_set_stats_text(format_stats_table(result.stats))
        self._pb_save_status.set("")
        self._set_status(f"Plot built via {result.backend_used}.")
        self._log(f"Plot backend: {result.backend_used}", "header")

    def _pb_save(self):
        """Export the current plot via plotting service to PNG/PDF."""
        if not self._pb_current_request:
            messagebox.showwarning("Save Plot", "Build a plot first.")
            return
        plots_dir = os.path.join(self.simv0_dir, "data", "plots")
        os.makedirs(plots_dir, exist_ok=True)

        path = filedialog.asksaveasfilename(
            initialdir=plots_dir,
            defaultextension=".png",
            filetypes=[("PNG image", "*.png"), ("PDF document", "*.pdf")],
            title="Save Plot",
        )
        if not path:
            return
        req = self._pb_make_request(output_path=path)
        req.output_format = "pdf" if path.lower().endswith(".pdf") else "png"
        result = self._plot_service.render(req)
        if not result.ok:
            messagebox.showerror("Save Plot", f"Could not save file:\n{result.error}")
            return
        short = os.path.basename(path)
        self._pb_save_status.set(f"Saved: {short}")
        self._set_status(f"Plot saved to {path}")

    def _pb_save_stats(self):
        if not self._pb_last_result:
            messagebox.showwarning("Save Stats", "Build a plot first.")
            return
        stats_dir = os.path.join(self.simv0_dir, "data", "plots")
        os.makedirs(stats_dir, exist_ok=True)
        path = filedialog.asksaveasfilename(
            initialdir=stats_dir,
            defaultextension=".txt",
            filetypes=[("Text file", "*.txt")],
            title="Save Statistics",
        )
        if not path:
            return
        text = self._pb_stats_text.get("1.0", tk.END).strip()
        try:
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(text + "\n")
            self._set_status(f"Stats saved to {path}")
            self._pb_save_status.set(f"Saved: {os.path.basename(path)}")
        except OSError as exc:
            messagebox.showerror("Save Stats", f"Could not save file:\n{exc}")

    def _pb_set_stats_text(self, text):
        self._pb_stats_text.configure(state=tk.NORMAL)
        self._pb_stats_text.delete("1.0", tk.END)
        self._pb_stats_text.insert("1.0", text)
        self._pb_stats_text.configure(state=tk.DISABLED)

    def _pb_clear(self):
        """Clear the embedded plot."""
        self._pb_fig.clear()
        ax = self._pb_fig.add_subplot(111)
        _style_axes(ax)
        self._pb_canvas.draw()
        self._pb_series = []
        self._pb_series_list.delete(0, tk.END)
        self._pb_last_result = None
        self._pb_current_request = None
        self._pb_set_stats_text("No statistics yet.")
        self._pb_save_status.set("")
        self._set_status("Plot cleared.")

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
            self.root.after(0, self._dv_refresh_and_load)   # auto-refresh data viewer
            if hasattr(self, "_pb_refresh_files"):
                self.root.after(0, self._pb_refresh_files)
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
#  Full .dat parser — returns (headers, rows) for Treeview / plot use
# ---------------------------------------------------------------------------

def _parse_dat_full(path, max_rows=None):
    """
    Parse a simv0 .dat file and return (headers, rows).

    headers : list[str]  — column names from the last #-comment header line
    rows    : list[tuple[str]]  — formatted numeric strings, right-aligned

    Comment lines (# …) and blank lines are skipped from data rows.
    If max_rows==0, only headers are returned (no rows read).
    """
    if not os.path.isfile(path):
        return [], []

    headers = []
    rows = []

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                parts = line.lstrip("#").split()
                if parts:
                    headers = parts
                continue
            if max_rows == 0:
                continue
            try:
                vals = [float(x) for x in line.split()]
            except ValueError:
                continue
            # Format for display: up to 8 sig-figs, avoid clutter
            rows.append(tuple(f"{v:.6g}" for v in vals))
            if max_rows and len(rows) >= max_rows:
                break

    return headers, rows


# ---------------------------------------------------------------------------
#  Matplotlib helpers
# ---------------------------------------------------------------------------

def _style_axes(ax, spine_color=None):
    """Apply the retro dark theme to a matplotlib Axes."""
    spine_color = spine_color or PAL["navy_mid"]
    ax.set_facecolor(PAL["console_bg"])
    ax.tick_params(colors=PAL["phosphor"], labelsize=7)
    ax.xaxis.label.set_color(PAL["header"])
    ax.yaxis.label.set_color(PAL["phosphor"])
    for spine in ax.spines.values():
        spine.set_edgecolor(spine_color)
    ax.grid(color="#1A2E1A", linestyle="--", linewidth=0.5, alpha=0.7)


def _plot_series(ax, xdata, ydata, ptype, color, marker, lw, label):
    """Dispatch a plot call based on ptype."""
    kw_shared = {"label": label, "color": color}
    if ptype == "scatter":
        ax.scatter(xdata, ydata, marker=marker or "o", s=36, **kw_shared)
    elif ptype == "bar":
        ax.bar(xdata, ydata, color=color, alpha=0.85, label=label)
    else:  # line (default)
        ax.plot(xdata, ydata, linewidth=lw, marker=marker,
                markersize=4 if marker else 0, **kw_shared)




def launch(simv0_dir=None):
    """Create the Tk root window and start the event loop."""
    root = tk.Tk()
    app = SimV0LauncherGUI(root, simv0_dir=simv0_dir)
    root.mainloop()


if __name__ == "__main__":
    launch()
