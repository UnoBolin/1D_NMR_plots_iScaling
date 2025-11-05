"""
TopSpin‑processed 1D plotter — multi‑dataset with **per‑trace scaling**

Adds:
- SCALE_FACTORS: list of per‑trace multipliers matching INPUT_PATHS length
  (e.g., [1.0, 1.0, 1.25, 0.9]). Traces without an entry fall back to SCALE_FACTOR.
- SCALE_FACTORS_MODE: "override" (default) or "multiply"
  • override → uses SCALE_FACTORS[i] if present, otherwise SCALE_FACTOR
  • multiply → uses SCALE_FACTOR * SCALE_FACTORS[i] (if present)

Everything else from your multi script is preserved: overlay/stack modes, ppm
limits, baseline offset, colormap, hidden Y‑axis, custom filenames, etc.

Install (venv):
    pip install nmrglue matplotlib numpy
"""
from __future__ import annotations

import os
import sys
from datetime import datetime
from typing import Optional, List, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize

try:
    import nmrglue as ng
except Exception as e:
    raise SystemExit("nmrglue is required. Install with `pip install nmrglue`. Error: " + str(e))

# =====================
# CONFIG — EDIT ME (drop in your values; keys are superset of your previous script)
# =====================
CONFIG = {
    # Input (single or multiple)
    "INPUT_PATH": r"/Users/katja/Desktop/Uno/NMR-data/20250630_CUC25nt_1mM_283K_Copy/2/pdata/1",  # kept for backward compatibility
    "INPUT_PATHS": [
        r"/Users/katja/Desktop/Uno/NMR-data/20250630_CUC25nt_1mM_1mM_Mg_283K/2/pdata/1",
        r"/Users/katja/Desktop/Uno/NMR-data/20250702_CUC25nt_1mM_2mM_Mg_283K/2/pdata/1",
        r"/Users/katja/Desktop/Uno/NMR-data/20250704_CUC25nt_1mM_5mM_Mg_283K/2/pdata/1",
        # r"/path/to/sampleA/…/pdata/1",
        # r"/path/to/sampleB/…/pdata/1",
        # r"/path/to/sampleC/…/pdata/1",
    ],
    "AUTO_FIND_PDATA": True,
    "USE_GUI_IF_INPUT_EMPTY": False,

    # Output
    "OUTPUT_DIR": r"/Users/katja/Desktop/Uno/Uno_Test_1Ds",                   # If empty, uses <first_pdata>/plots
    "SAVE_FORMATS": ["png"],
    "DPI": 300,

    # Filename controls (individual files still use these; combined has its own template below)
    "OUTPUT_BASENAME": r"Testinggg",              # e.g. "CUC25_Test"; empty → default per dataset
    "FILENAME_TEMPLATE": "{basename}_S{scale:g}_{ppmmin}-{ppmmax}_{date}_{time}",
    "AVOID_OVERWRITE": True,

    # Combined figure naming (for overlay/stack)
    "COMBINED_BASENAME": r"",            # if empty, derived from first dataset basename + mode
    "COMBINED_TEMPLATE": "{basename}_{mode}_{n}traces_{ppmmin}-{ppmmax}_{date}_{time}",

    # Appearance
    "FIGSIZE": (10, 6),
    "LINEWIDTH": 2.5,
    "TITLE_TEMPLATE": "",               # Empty → no title
    "SHOW_GRID": False,

    # Axis / scaling
    "REVERSE_PPM": True,
    "NORMALIZE": True,
    "SCALE_FACTOR": 15.0,                 # global default
    # NEW: per‑trace scaling controls
    "SCALE_FACTORS": [1, 1, 1, 4],                 # e.g., [1.0, 1.0, 1.25, 0.9]
    "SCALE_FACTORS_MODE": "multiply",   # "override" or "multiply"

    # Baseline offset (applied to every trace)
    "BASELINE_OFFSET": -1.07,
    "APPLY_OFFSET_AFTER_SCALING": True,

    # Y‑axis visibility
    "HIDE_Y_AXIS": True,                 # hide label, ticks, tick labels, spines

    # Limits — fixed y‑limits for reproducible visuals
    "AUTO_YLIM": False,
    "YLIM": None,                        # if None and AUTO_YLIM=False → (-1.1, 1.1)

    # PPM window
    "PPM_MIN": 14,
    "PPM_MAX": 9.5,

    # Colormap
    "USE_COLORMAP": True,
    "COLORMAP": "viridis",
    "COLORMAP_FIXED_NORM": None,

    # Multi‑dataset controls
    "MULTI_MODE": "Stack",            # "overlay" or "stack"
    "LABELS": [],                        # optional labels per path; fallback to expno/procno
    "COLORS": [],                        # optional colors per path; fallback to matplotlib cycle
    "LEGEND": True,                      # show legend for overlay

    # Stacking parameters (used when MULTI_MODE="stack")
    "STACK_STEP": 2.15,                   # vertical separation between traces
    "STACK_DIRECTION": "up",          # "down" or "up"

    # Debug prints
    "VERBOSE": False,
}
# =====================


def _is_pdata_dir(path: str) -> bool:
    return os.path.isfile(os.path.join(path, "procs")) and (
        os.path.isfile(os.path.join(path, "1r")) or os.path.isfile(os.path.join(path, "1rr"))
    )


def _discover_pdata_dir(start: str) -> Optional[str]:
    if _is_pdata_dir(start):
        return start

    candidates: List[str] = []

    def collect_pdata(root: str):
        pdata_base = os.path.join(root, "pdata")
        if os.path.isdir(pdata_base):
            for name in sorted(os.listdir(pdata_base)):
                p = os.path.join(pdata_base, name)
                if os.path.isdir(p) and _is_pdata_dir(p):
                    candidates.append(p)

    collect_pdata(start)
    parent = os.path.dirname(start.rstrip(os.sep))
    if parent and parent != start:
        collect_pdata(parent)
    try:
        for child in sorted(os.listdir(start)):
            ch_path = os.path.join(start, child)
            if os.path.isdir(ch_path) and child.isdigit():
                collect_pdata(ch_path)
    except Exception:
        pass

    def procno_key(p: str) -> int:
        try:
            return int(os.path.basename(p))
        except ValueError:
            return -1

    candidates.sort(key=procno_key, reverse=True)
    return candidates[0] if candidates else None


def _maybe_pick_dir() -> Optional[str]:
    if not CONFIG.get("USE_GUI_IF_INPUT_EMPTY", False):
        return None
    try:
        import tkinter as tk
        from tkinter import filedialog
    except Exception:
        return None
    root = tk.Tk(); root.withdraw()
    return filedialog.askdirectory(title="Select Bruker dataset folder (pdata/<procno> or parent)") or None


def _ppm_axis_from_procs(dic: dict, npoints: int) -> np.ndarray:
    procs = dic.get("procs", {})
    try:
        offset = float(procs.get("OFFSET"))
    except Exception:
        offset = 0.0
    try:
        sw_hz = float(procs.get("SW_p", procs.get("SW", 0.0)))
    except Exception:
        sw_hz = 0.0
    try:
        sf_mhz = float(procs.get("SF", procs.get("SFO1", 1.0)))
    except Exception:
        sf_mhz = 1.0
    try:
        si = int(procs.get("SI", npoints))
    except Exception:
        si = npoints

    width_ppm = (sw_hz / max(sf_mhz, 1e-12)) if sw_hz else 0.0
    step = width_ppm / max(si, 1)
    idx = np.arange(npoints)
    if npoints != si and si > 1:
        step *= si / npoints
    return offset - idx * step


def _load_processed_bruker(pdata_dir: str):
    dic, data = ng.bruker.read_pdata(pdata_dir)
    y = np.asarray(data).squeeze().astype(float)

    x_ppm: np.ndarray
    make_uc = getattr(ng.bruker, "make_uc", None)
    if callable(make_uc):
        try:
            uc = make_uc(dic, y)
            x_ppm = uc.ppm_scale()
        except Exception:
            x_ppm = _ppm_axis_from_procs(dic, y.size)
    else:
        x_ppm = _ppm_axis_from_procs(dic, y.size)

    procno = os.path.basename(pdata_dir)
    expno = os.path.basename(os.path.dirname(os.path.dirname(pdata_dir)))
    label = f"{expno}/pdata/{procno}"
    return x_ppm, y, dic, expno, procno, label


def _scale_for_index(i: int) -> float:
    base = float(CONFIG.get("SCALE_FACTOR", 1.0))
    lst = CONFIG.get("SCALE_FACTORS") or []
    mode = str(CONFIG.get("SCALE_FACTORS_MODE", "override")).lower()
    if i < len(lst) and lst[i] not in (None, ""):
        v = float(lst[i])
        return base * v if mode.startswith("mult") else v
    return base


def _prepare_y(y: np.ndarray, i: int) -> np.ndarray:
    cfg = CONFIG
    y_out = y.astype(float)

    if cfg.get("NORMALIZE", True):
        denom = np.max(np.abs(y_out)) or 1.0
        y_out = y_out / denom

    y_out *= _scale_for_index(i)

    if cfg.get("APPLY_OFFSET_AFTER_SCALING", True):
        y_out = y_out + float(cfg.get("BASELINE_OFFSET", 0.0))
    else:
        y_out = (y_out + float(cfg.get("BASELINE_OFFSET", 0.0)))

    return y_out


def _hide_y_axis(ax):
    ax.set_ylabel("")
    ax.yaxis.set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _determine_ylims(y: np.ndarray) -> Tuple[float, float]:
    cfg = CONFIG
    if cfg.get("AUTO_YLIM", False):
        ymin, ymax = float(np.min(y)), float(np.max(y))
        pad = 0.05 * (ymax - ymin if ymax > ymin else 1.0)
        return ymin - pad, ymax + pad
    if isinstance(cfg.get("YLIM"), (list, tuple)) and len(cfg["YLIM"]) == 2:
        return float(cfg["YLIM"][0]), float(cfg["YLIM"][1])
    return -1.1, 1.1


def _safe_filename(name: str) -> str:
    bad = "/\\?%*:|\"<>"
    for ch in bad:
        name = name.replace(ch, "-")
    return name.strip().strip(".")


def _build_basename(expno: str, procno: str) -> str:
    base_cfg = CONFIG.get("OUTPUT_BASENAME", "")
    if base_cfg:
        return base_cfg
    return f"{expno}_pdata_{procno}"


def _filename_from_template(expno: str, procno: str, scale_used: float) -> str:
    tpl = CONFIG.get("FILENAME_TEMPLATE", "{basename}")
    basename = _build_basename(expno, procno)
    ppm_min = CONFIG.get("PPM_MIN")
    ppm_max = CONFIG.get("PPM_MAX")
    if ppm_min is None or ppm_max is None:
        ppmmin_s, ppmmax_s = "auto", "auto"
    else:
        ppmmin_s, ppmmax_s = str(ppm_min), str(ppm_max)
    now = datetime.now()
    res = tpl.format(
        basename=basename,
        expno=expno,
        procno=procno,
        scale=scale_used,
        ppmmin=ppmmin_s,
        ppmmax=ppmmax_s,
        date=now.strftime("%Y%m%d"),
        time=now.strftime("%H%M%S"),
    )
    return _safe_filename(res)


def _combined_filename(n: int, mode: str, first_expno: str, first_procno: str) -> str:
    tpl = CONFIG.get("COMBINED_TEMPLATE", "{basename}_{mode}_{n}traces_{ppmmin}-{ppmmax}_{date}_{time}")
    base = CONFIG.get("COMBINED_BASENAME") or _build_basename(first_expno, first_procno)
    ppm_min = CONFIG.get("PPM_MIN")
    ppm_max = CONFIG.get("PPM_MAX")
    if ppm_min is None or ppm_max is None:
        ppmmin_s, ppmmax_s = "auto", "auto"
    else:
        ppmmin_s, ppmmax_s = str(ppm_min), str(ppm_max)
    now = datetime.now()
    res = tpl.format(
        basename=_safe_filename(base),
        mode=mode,
        n=n,
        ppmmin=ppmmin_s,
        ppmmax=ppmmax_s,
        date=now.strftime("%Y%m%d"),
        time=now.strftime("%H%M%S"),
    )
    return _safe_filename(res)


def _next_nonconflicting_path(base_noext: str, fmt: str, out_dir: str) -> str:
    path = os.path.join(out_dir, f"{base_noext}.{fmt}")
    if not CONFIG.get("AVOID_OVERWRITE", True):
        return path
    if not os.path.exists(path):
        return path
    for k in range(1, 1000):
        candidate = os.path.join(out_dir, f"{base_noext}_{k:03d}.{fmt}")
        if not os.path.exists(candidate):
            return candidate
    return path


def _resolve_paths() -> List[str]:
    paths = list(CONFIG.get("INPUT_PATHS", []))
    if len(sys.argv) > 1 and os.path.exists(sys.argv[1]):
        paths.insert(0, sys.argv[1])
    elif CONFIG.get("INPUT_PATH"):
        paths.insert(0, CONFIG["INPUT_PATH"])

    # Fallback to GUI if still empty and enabled
    if not paths:
        p = _maybe_pick_dir()
        if p:
            paths = [p]

    if not paths:
        raise SystemExit("No input path(s) provided. Fill CONFIG['INPUT_PATHS'] or 'INPUT_PATH', or pass argv[1].")

    # Auto‑discover pdata for each
    final_paths: List[str] = []
    for p in paths:
        if CONFIG.get("AUTO_FIND_PDATA", True):
            pdata = _discover_pdata_dir(p) or p
        else:
            pdata = p
        if not _is_pdata_dir(pdata):
            raise SystemExit(f"Not a processed TopSpin directory: {p}")
        final_paths.append(pdata)
    return final_paths


def _load_all(paths: List[str]):
    loaded = []
    for p in paths:
        x_ppm, y, dic, expno, procno, label = _load_processed_bruker(p)
        loaded.append((x_ppm, y, dic, expno, procno, label))
    return loaded


def _prep_all(loaded):
    prepped = []
    for i, (x_ppm, y, dic, expno, procno, label) in enumerate(loaded):
        y_p = _prepare_y(y, i)
        prepped.append((x_ppm, y_p, dic, expno, procno, label, _scale_for_index(i)))
    return prepped


def _apply_ppm_limits(ax):
    if CONFIG.get("REVERSE_PPM", True):
        ax.invert_xaxis()
    if CONFIG.get("PPM_MIN") is not None and CONFIG.get("PPM_MAX") is not None:
        ax.set_xlim(CONFIG["PPM_MIN"], CONFIG["PPM_MAX"])


def _common_ylim_from_all(prepped) -> Tuple[float, float]:
    # Compute limits across all prepped y arrays
    ymins = [float(np.min(y_p)) for (_, y_p, *_) in prepped]
    ymaxs = [float(np.max(y_p)) for (_, y_p, *_) in prepped]
    ymin, ymax = min(ymins), max(ymaxs)
    if CONFIG.get("AUTO_YLIM", False):
        pad = 0.05 * (ymax - ymin if ymax > ymin else 1.0)
        return ymin - pad, ymax + pad
    # If fixed, respect CONFIG['YLIM']
    if isinstance(CONFIG.get("YLIM"), (list, tuple)) and len(CONFIG["YLIM"]) == 2:
        return float(CONFIG["YLIM"][0]), float(CONFIG["YLIM"][1])
    # Default fixed box that fits all traces reasonably
    rng = max(abs(ymin), abs(ymax))
    rng = 1.1 if rng == 0 else rng * 1.05
    return -rng, rng


def _plot_overlay(prepped, out_dir: str):
    n = len(prepped)
    fig = plt.figure(figsize=CONFIG.get("FIGSIZE", (10, 6)))
    ax = fig.add_subplot(111)

    ymin, ymax = _common_ylim_from_all(prepped)

    colors = CONFIG.get("COLORS", [])
    labels = CONFIG.get("LABELS", [])

    for i, (x_ppm, y_p, dic, expno, procno, default_label, scale_used) in enumerate(prepped):
        color = colors[i] if i < len(colors) and colors[i] else None
        label = labels[i] if i < len(labels) and labels[i] else default_label

        if CONFIG.get("USE_COLORMAP", False):
            norm_cfg = CONFIG.get("COLORMAP_FIXED_NORM")
            if isinstance(norm_cfg, (list, tuple)) and len(norm_cfg) == 2:
                norm = Normalize(vmin=float(norm_cfg[0]), vmax=float(norm_cfg[1]))
            else:
                norm = Normalize(vmin=ymin, vmax=ymax)
            pts = np.column_stack([x_ppm, y_p])
            segs = np.stack([pts[:-1], pts[1:]], axis=1)
            lc = LineCollection(segs, cmap=plt.get_cmap(CONFIG.get("COLORMAP", "viridis")), norm=norm, linewidths=CONFIG.get("LINEWIDTH", 2.0))
            lc.set_array((y_p[:-1] + y_p[1:]) / 2.0)
            ax.add_collection(lc)
        else:
            ax.plot(x_ppm, y_p, lw=CONFIG.get("LINEWIDTH", 2.0), color=color, label=f"{label} (×{scale_used:g})")

    _apply_ppm_limits(ax)
    ax.set_ylim(ymin, ymax)
    if CONFIG.get("HIDE_Y_AXIS", True):
        _hide_y_axis(ax)

    if CONFIG.get("TITLE_TEMPLATE", ""):
        t = CONFIG["TITLE_TEMPLATE"].format(expno=prepped[0][3], procno=prepped[0][4])
        ax.set_title(t)

    if CONFIG.get("SHOW_GRID", False):
        ax.grid(True, ls=":", alpha=0.6)

    if CONFIG.get("LEGEND", True) and not CONFIG.get("USE_COLORMAP", False):
        ax.legend(frameon=False, loc="upper right")

    plt.tight_layout()

    # Save combined
    first_expno, first_procno = prepped[0][3], prepped[0][4]
    base_noext = _combined_filename(n, "overlay", first_expno, first_procno)
    os.makedirs(out_dir, exist_ok=True)
    for fmt in CONFIG.get("SAVE_FORMATS", ["png"]):
        path = _next_nonconflicting_path(base_noext, fmt, out_dir)
        plt.savefig(path, dpi=CONFIG.get("DPI", 300), bbox_inches="tight")
        print(f"Saved: {path}")

    plt.show()


def _plot_stack(prepped, out_dir: str):
    n = len(prepped)
    fig = plt.figure(figsize=CONFIG.get("FIGSIZE", (10, 6)))
    ax = fig.add_subplot(111)

    # Base y limits from the **first** (for scale), then compute stack envelope
    y0 = prepped[0][1]
    ymin0, ymax0 = _determine_ylims(y0)
    step = float(CONFIG.get("STACK_STEP", 0.6))
    direction = -1.0 if str(CONFIG.get("STACK_DIRECTION", "down")).lower().startswith("d") else 1.0

    colors = CONFIG.get("COLORS", [])
    labels = CONFIG.get("LABELS", [])

    for i, (x_ppm, y_p, dic, expno, procno, default_label, scale_used) in enumerate(prepped):
        y_off = y_p + direction * i * step
        color = colors[i] if i < len(colors) and colors[i] else None
        label = labels[i] if i < len(labels) and labels[i] else default_label

        if CONFIG.get("USE_COLORMAP", False):
            norm_cfg = CONFIG.get("COLORMAP_FIXED_NORM")
            if isinstance(norm_cfg, (list, tuple)) and len(norm_cfg) == 2:
                norm = Normalize(vmin=float(norm_cfg[0]), vmax=float(norm_cfg[1]))
            else:
                norm = Normalize(vmin=ymin0, vmax=ymax0)
            pts = np.column_stack([x_ppm, y_off])
            segs = np.stack([pts[:-1], pts[1:]], axis=1)
            lc = LineCollection(segs, cmap=plt.get_cmap(CONFIG.get("COLORMAP", "viridis")), norm=norm, linewidths=CONFIG.get("LINEWIDTH", 2.0))
            lc.set_array((y_p[:-1] + y_p[1:]) / 2.0)
            ax.add_collection(lc)
        else:
            ax.plot(x_ppm, y_off, lw=CONFIG.get("LINEWIDTH", 2.0), color=color, label=f"{label} (×{scale_used:g})")

    _apply_ppm_limits(ax)

    # Compute stacked y‑limits to include all offsets
    ymin_all = ymin0 + direction * 0 * step
    ymax_all = ymax0 + direction * (n - 1) * step
    ylo, yhi = (min(ymin_all, ymax_all) - 0.05, max(ymin_all, ymax_all) + 0.05)
    ax.set_ylim(ylo, yhi)

    if CONFIG.get("HIDE_Y_AXIS", True):
        _hide_y_axis(ax)

    if CONFIG.get("TITLE_TEMPLATE", ""):
        t = CONFIG["TITLE_TEMPLATE"].format(expno=prepped[0][3], procno=prepped[0][4])
        ax.set_title(t)

    if CONFIG.get("SHOW_GRID", False):
        ax.grid(True, ls=":", alpha=0.6)

    # Right-edge labels for stacked plot
    try:
        for i, (x_ppm, y_p, dic, expno, procno, default_label, scale_used) in enumerate(prepped):
            label = labels[i] if i < len(labels) and labels[i] else default_label
            x_end = x_ppm[-1]
            y_end = (y_p + direction * i * step)[-1]
            ax.text(x_end, y_end, f"  {label} (×{scale_used:g})", va="center", fontsize=9)
    except Exception:
        pass

    plt.tight_layout()

    # Save combined
    first_expno, first_procno = prepped[0][3], prepped[0][4]
    base_noext = _combined_filename(n, "stack", first_expno, first_procno)
    os.makedirs(out_dir, exist_ok=True)
    for fmt in CONFIG.get("SAVE_FORMATS", ["png"]):
        path = _next_nonconflicting_path(base_noext, fmt, out_dir)
        plt.savefig(path, dpi=CONFIG.get("DPI", 300), bbox_inches="tight")
        print(f"Saved: {path}")

    plt.show()


def main():
    paths = _resolve_paths()
    loaded = _load_all(paths)
    prepped = _prep_all(loaded)

    out_dir = CONFIG.get("OUTPUT_DIR") or os.path.join(paths[0], "plots")

    mode = str(CONFIG.get("MULTI_MODE", "overlay")).lower()
    if mode.startswith("over"):
        _plot_overlay(prepped, out_dir)
    elif mode.startswith("stack"):
        _plot_stack(prepped, out_dir)
    else:
        raise SystemExit("Unknown MULTI_MODE. Use 'overlay' or 'stack'.")


if __name__ == "__main__":
    main()
