"""Plotting service for simv0 GUI.

Provides backend-agnostic plot request models, dynamic .dat ingestion,
unit parsing/validation, statistical summaries, and rendering via gnuplot
(primary) or matplotlib (fallback).
"""

from __future__ import annotations

import math
import os
import re
import shutil
import statistics
import subprocess
import tempfile
from dataclasses import dataclass, field, replace
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
from matplotlib.figure import Figure


_HEADER_RE = re.compile(r"^(?P<name>[^()]+?)(?:\((?P<unit>[^()]+)\))?$")


@dataclass(frozen=True)
class ColumnMeta:
    raw: str
    name: str
    unit: str


@dataclass
class ParsedData:
    path: str
    headers: List[str]
    columns: List[ColumnMeta]
    rows: List[List[float]]

    def header_index(self, header: str) -> int:
        return self.headers.index(header)


@dataclass
class SeriesSpec:
    file_path: str
    y_column: str
    label: str
    color: str = "#33FF00"
    marker: Optional[str] = None
    line_width: float = 1.5
    yerr_column: Optional[str] = None
    z_column: Optional[str] = None


@dataclass
class PlotRequest:
    plot_mode: str  # "2d" | "3d"
    plot_type: str  # "line" | "scatter" | "bar"
    backend: str  # "auto" | "gnuplot" | "matplotlib"
    x_column: str
    series: List[SeriesSpec]
    title: str = ""
    x_label: Optional[str] = None
    y_label: Optional[str] = None
    z_label: Optional[str] = None
    output_path: Optional[str] = None
    output_format: str = "png"
    include_stats: bool = True


@dataclass
class SeriesStats:
    label: str
    n: int
    mean: float
    std: float
    minimum: float
    maximum: float
    ci95_low: float
    ci95_high: float
    correlation_r: Optional[float] = None
    slope: Optional[float] = None
    intercept: Optional[float] = None


@dataclass
class PlotResult:
    ok: bool
    backend_used: str
    output_path: Optional[str]
    error: str = ""
    stats: List[SeriesStats] = field(default_factory=list)


PRESET_TEMPLATES: Dict[str, Dict[str, str]] = {
    "Temperature comparison (classic)": {
        "plot_mode": "2d",
        "plot_type": "line",
        "x_column": "t(s)",
        "y_label": "Temperature (K)",
        "title": "Temperature comparison",
    },
    "Power comparison": {
        "plot_mode": "2d",
        "plot_type": "line",
        "x_column": "t(s)",
        "y_label": "Power (W)",
        "title": "Power comparison",
    },
    "Wind and irradiance": {
        "plot_mode": "2d",
        "plot_type": "line",
        "x_column": "t(s)",
        "y_label": "Value",
        "title": "Wind speed and irradiance",
    },
    "Error bars (measurement style)": {
        "plot_mode": "2d",
        "plot_type": "line",
        "x_column": "t(s)",
        "title": "Series with error bars",
    },
    "3D surface/scatter": {
        "plot_mode": "3d",
        "plot_type": "scatter",
        "x_column": "t(s)",
        "title": "3D series view",
    },
}


def discover_dat_files(data_dir: str) -> List[str]:
    if not os.path.isdir(data_dir):
        return []
    files = [f for f in os.listdir(data_dir) if f.lower().endswith(".dat")]
    return sorted(files)


def parse_header_token(token: str) -> ColumnMeta:
    tok = token.strip()
    m = _HEADER_RE.match(tok)
    if not m:
        return ColumnMeta(raw=tok, name=tok, unit="")
    name = (m.group("name") or tok).strip()
    unit = (m.group("unit") or "").strip()
    return ColumnMeta(raw=tok, name=name, unit=unit)


def parse_dat_full_numeric(path: str, max_rows: Optional[int] = None) -> ParsedData:
    if not os.path.isfile(path):
        return ParsedData(path=path, headers=[], columns=[], rows=[])

    headers: List[str] = []
    rows: List[List[float]] = []

    with open(path, encoding="utf-8") as fh:
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
            rows.append(vals)
            if max_rows and len(rows) >= max_rows:
                break

    cols = [parse_header_token(h) for h in headers]
    return ParsedData(path=path, headers=headers, columns=cols, rows=rows)


def validate_series_units(parsed_by_file: Dict[str, ParsedData], series: Sequence[SeriesSpec]) -> Tuple[bool, str]:
    if not series:
        return False, "No series selected."

    units = set()
    for s in series:
        pdata = parsed_by_file.get(s.file_path)
        if not pdata or s.y_column not in pdata.headers:
            return False, f"Column '{s.y_column}' not found in {os.path.basename(s.file_path)}"
        col = pdata.columns[pdata.header_index(s.y_column)]
        units.add(col.unit)

    if len(units) > 1:
        shown = ", ".join(sorted(u or "(unitless)" for u in units))
        return False, f"Incompatible Y units across series: {shown}"
    return True, ""


def _safe_std(vals: Sequence[float]) -> float:
    return statistics.stdev(vals) if len(vals) > 1 else 0.0


def _extract_xy(rows: Sequence[Sequence[float]], xi: int, yi: int) -> Tuple[List[float], List[float]]:
    x: List[float] = []
    y: List[float] = []
    need = max(xi, yi)
    for r in rows:
        if len(r) <= need:
            continue
        x.append(r[xi])
        y.append(r[yi])
    return x, y


def _extract_xye(rows: Sequence[Sequence[float]], xi: int, yi: int, ei: int) -> Tuple[List[float], List[float], List[float]]:
    x: List[float] = []
    y: List[float] = []
    e: List[float] = []
    need = max(xi, yi, ei)
    for r in rows:
        if len(r) <= need:
            continue
        x.append(r[xi])
        y.append(r[yi])
        e.append(r[ei])
    return x, y, e


def _gp_quote(text: str) -> str:
    safe = text.replace("\\", "\\\\").replace("\"", "\\\"").replace("\n", " ").replace("\r", " ")
    return f"\"{safe}\""


def _linreg(x: Sequence[float], y: Sequence[float]) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    if len(x) < 2 or len(x) != len(y):
        return None, None, None
    xmean = statistics.mean(x)
    ymean = statistics.mean(y)
    sxx = sum((xi - xmean) ** 2 for xi in x)
    syy = sum((yi - ymean) ** 2 for yi in y)
    sxy = sum((xi - xmean) * (yi - ymean) for xi, yi in zip(x, y))
    if sxx == 0:
        return None, None, None
    slope = sxy / sxx
    intercept = ymean - slope * xmean
    if sxx > 0 and syy > 0:
        r = sxy / math.sqrt(sxx * syy)
    else:
        r = None
    return r, slope, intercept


def compute_stats(label: str, x: Sequence[float], y: Sequence[float]) -> SeriesStats:
    n = len(y)
    if n == 0:
        return SeriesStats(label=label, n=0, mean=float("nan"), std=float("nan"),
                           minimum=float("nan"), maximum=float("nan"),
                           ci95_low=float("nan"), ci95_high=float("nan"))

    mean = statistics.mean(y)
    std = _safe_std(y)
    sem = std / math.sqrt(n) if n > 1 else 0.0
    ci = 1.96 * sem  # Approximate 95% CI using standard normal z-score.
    r, slope, intercept = _linreg(x, y)

    return SeriesStats(
        label=label,
        n=n,
        mean=mean,
        std=std,
        minimum=min(y),
        maximum=max(y),
        ci95_low=mean - ci,
        ci95_high=mean + ci,
        correlation_r=r,
        slope=slope,
        intercept=intercept,
    )


def format_stats_table(stats: Sequence[SeriesStats]) -> str:
    if not stats:
        return "No statistics available."
    lines = [
        "label | n | mean | std | min | max | 95% CI | r | slope | intercept",
        "-" * 96,
    ]
    for s in stats:
        ci_txt = f"[{s.ci95_low:.5g}, {s.ci95_high:.5g}]"
        r_txt = "-" if s.correlation_r is None else f"{s.correlation_r:.5g}"
        m_txt = "-" if s.slope is None else f"{s.slope:.5g}"
        b_txt = "-" if s.intercept is None else f"{s.intercept:.5g}"
        lines.append(
            f"{s.label} | {s.n} | {s.mean:.5g} | {s.std:.5g} | {s.minimum:.5g} | "
            f"{s.maximum:.5g} | {ci_txt} | {r_txt} | {m_txt} | {b_txt}"
        )
    return "\n".join(lines)


class PlottingService:
    def __init__(self):
        self._cache: Dict[str, ParsedData] = {}

    def parse_file(self, path: str, max_rows: Optional[int] = None) -> ParsedData:
        key = f"{path}:{max_rows}"
        if max_rows is None and path in self._cache:
            return self._cache[path]
        parsed = parse_dat_full_numeric(path, max_rows=max_rows)
        if max_rows is None:
            self._cache[path] = parsed
        return parsed

    def clear_cache(self):
        self._cache.clear()

    def render(self, req: PlotRequest) -> PlotResult:
        parsed_by_file: Dict[str, ParsedData] = {}
        for s in req.series:
            parsed_by_file[s.file_path] = self.parse_file(s.file_path)

        ok, msg = validate_series_units(parsed_by_file, req.series)
        if not ok:
            return PlotResult(ok=False, backend_used="none", output_path=None, error=msg)

        for p in parsed_by_file.values():
            if req.x_column not in p.headers:
                return PlotResult(ok=False, backend_used="none", output_path=None,
                                  error=f"X column '{req.x_column}' not found in {os.path.basename(p.path)}")

        out_path = req.output_path
        if not out_path:
            ext = req.output_format.lower()
            if ext not in ("png", "pdf"):
                ext = "png"
            fd, out_path = tempfile.mkstemp(prefix="simv0_plot_", suffix=f".{ext}")
            os.close(fd)

        stats: List[SeriesStats] = []
        if req.include_stats:
            for s in req.series:
                pdata = parsed_by_file[s.file_path]
                xi = pdata.header_index(req.x_column)
                yi = pdata.header_index(s.y_column)
                x, y = _extract_xy(pdata.rows, xi, yi)
                stats.append(compute_stats(s.label, x, y))

        backend_order = [req.backend]
        if req.backend == "auto":
            backend_order = ["gnuplot", "matplotlib"]

        errors = []
        for backend in backend_order:
            if backend == "gnuplot":
                res = self._render_gnuplot(req, parsed_by_file, out_path)
            else:
                res = self._render_matplotlib(req, parsed_by_file, out_path)
            if res.ok:
                res.stats = stats
                return res
            errors.append(f"{backend}: {res.error}")

        return PlotResult(ok=False, backend_used="none", output_path=None,
                          error="; ".join(errors), stats=stats)

    def _render_gnuplot(self, req: PlotRequest, parsed_by_file: Dict[str, ParsedData], out_path: str) -> PlotResult:
        if shutil.which("gnuplot") is None:
            return PlotResult(ok=False, backend_used="gnuplot", output_path=None, error="gnuplot not found in PATH")

        ext = os.path.splitext(out_path)[1].lower()
        if ext == ".pdf":
            term = "set terminal pdfcairo enhanced color\n"
        else:
            term = "set terminal pngcairo size 1400,900 enhanced\n"

        script_lines = [
            term,
            f"set output {_gp_quote(out_path.replace('\\', '/'))}",
            "set grid",
            "set key left top",
            f"set title {_gp_quote(req.title or 'simv0 plot')}",
            f"set xlabel {_gp_quote(req.x_label or req.x_column)}",
        ]

        ylabel = req.y_label or req.series[0].y_column
        script_lines.append(f"set ylabel {_gp_quote(ylabel)}")

        if req.plot_mode == "3d":
            zlabel = req.z_label or "Z"
            script_lines.append(f"set zlabel {_gp_quote(zlabel)}")

        plots: List[str] = []
        for s in req.series:
            pdata = parsed_by_file[s.file_path]
            xi = pdata.header_index(req.x_column) + 1
            yi = pdata.header_index(s.y_column) + 1
            real_path = os.path.realpath(s.file_path)
            if (not os.path.isfile(real_path)) or (not real_path.lower().endswith(".dat")):
                return PlotResult(ok=False, backend_used="gnuplot", output_path=None,
                                  error=f"Invalid data source for series '{s.label}'")
            path = real_path.replace("\\", "/")
            if req.plot_mode == "3d":
                if not s.z_column or s.z_column not in pdata.headers:
                    return PlotResult(ok=False, backend_used="gnuplot", output_path=None,
                                      error=f"Z column missing for series '{s.label}'")
                zi = pdata.header_index(s.z_column) + 1
                style = "with points pt 7" if req.plot_type == "scatter" else "with lines"
                plots.append(f"{_gp_quote(path)} using {xi}:{yi}:{zi} {style} lw {s.line_width} title {_gp_quote(s.label)}")
            else:
                if s.yerr_column and s.yerr_column in pdata.headers:
                    ei = pdata.header_index(s.yerr_column) + 1
                    style = "with yerrorbars"
                    plots.append(f"{_gp_quote(path)} using {xi}:{yi}:{ei} {style} lw {s.line_width} title {_gp_quote(s.label)}")
                else:
                    if req.plot_type == "scatter":
                        style = "with points pt 7"
                    elif req.plot_type == "bar":
                        style = "with boxes"
                    else:
                        style = "with lines"
                    plots.append(f"{_gp_quote(path)} using {xi}:{yi} {style} lw {s.line_width} title {_gp_quote(s.label)}")

        cmd = "splot" if req.plot_mode == "3d" else "plot"
        script_lines.append(f"{cmd} " + ", \\\n  ".join(plots))
        script = "\n".join(script_lines) + "\n"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".gp", delete=False, encoding="utf-8") as tf:
            tf.write(script)
            gp_path = tf.name

        try:
            run = subprocess.run(["gnuplot", gp_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if run.returncode != 0:
                err = (run.stderr or run.stdout or "gnuplot failed").strip()
                return PlotResult(ok=False, backend_used="gnuplot", output_path=None, error=err)
            return PlotResult(ok=True, backend_used="gnuplot", output_path=out_path)
        finally:
            try:
                os.remove(gp_path)
            except OSError:
                pass

    def _render_matplotlib(self, req: PlotRequest, parsed_by_file: Dict[str, ParsedData], out_path: str) -> PlotResult:
        fig = Figure(figsize=(10, 6), dpi=120)
        if req.plot_mode == "3d":
            ax = fig.add_subplot(111, projection="3d")
        else:
            ax = fig.add_subplot(111)

        for s in req.series:
            pdata = parsed_by_file[s.file_path]
            xi = pdata.header_index(req.x_column)
            yi = pdata.header_index(s.y_column)
            x, y = _extract_xy(pdata.rows, xi, yi)

            if req.plot_mode == "3d":
                if not s.z_column or s.z_column not in pdata.headers:
                    return PlotResult(ok=False, backend_used="matplotlib", output_path=None,
                                      error=f"Z column missing for series '{s.label}'")
                zi = pdata.header_index(s.z_column)
                z = [r[zi] for r in pdata.rows if len(r) > max(xi, yi, zi)]
                if req.plot_type == "scatter":
                    ax.scatter(x, y, z, label=s.label, c=s.color, marker=s.marker or "o", s=14)
                else:
                    ax.plot(x, y, z, label=s.label, color=s.color, linewidth=s.line_width,
                            marker=s.marker, markersize=3 if s.marker else 0)
            else:
                if s.yerr_column and s.yerr_column in pdata.headers:
                    ei = pdata.header_index(s.yerr_column)
                    x_err, y_err, y_error = _extract_xye(pdata.rows, xi, yi, ei)
                    ax.errorbar(x_err, y_err, yerr=y_error, label=s.label, color=s.color,
                                linewidth=s.line_width, marker=s.marker or "o", markersize=3, capsize=2)
                elif req.plot_type == "scatter":
                    ax.scatter(x, y, label=s.label, c=s.color, marker=s.marker or "o", s=14)
                elif req.plot_type == "bar":
                    ax.bar(x, y, label=s.label, color=s.color, alpha=0.75)
                else:
                    ax.plot(x, y, label=s.label, color=s.color, linewidth=s.line_width,
                            marker=s.marker, markersize=3 if s.marker else 0)

        ax.set_title(req.title or "simv0 plot")
        ax.set_xlabel(req.x_label or req.x_column)
        if req.plot_mode == "3d":
            ax.set_ylabel(req.y_label or req.series[0].y_column)
            ax.set_zlabel(req.z_label or "Z")
        else:
            ax.set_ylabel(req.y_label or req.series[0].y_column)
        ax.grid(True, linestyle="--", alpha=0.4)
        ax.legend(loc="best", fontsize=8)
        fig.tight_layout()

        try:
            fig.savefig(out_path, bbox_inches="tight")
        except OSError as exc:
            return PlotResult(ok=False, backend_used="matplotlib", output_path=None, error=str(exc))
        return PlotResult(ok=True, backend_used="matplotlib", output_path=out_path)


def apply_preset_to_request(req: PlotRequest, preset_name: str) -> PlotRequest:
    preset = PRESET_TEMPLATES.get(preset_name)
    if not preset:
        return req
    merged = req
    for key, val in preset.items():
        if hasattr(merged, key):
            merged = replace(merged, **{key: val})
    return merged
