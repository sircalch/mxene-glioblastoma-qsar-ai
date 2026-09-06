"""
_pubstyle.py - one house style for every matplotlib figure in this manuscript.

Call `apply()` once at import time in each figure generator. Gives a consistent
Q1-journal look: Arial, tight sizing, colour-blind-safe palette, clean spines,
300 dpi raster + matching vector PDF, and helpers for panel labels / saving.
"""
from __future__ import annotations
import os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# ---- fonts -----------------------------------------------------------------
_AVAIL = {f.name for f in fm.fontManager.ttflist}
_SANS = next((n for n in ("Arial", "Helvetica", "TeX Gyre Heros",
                          "Liberation Sans", "DejaVu Sans") if n in _AVAIL),
             "DejaVu Sans")

# ---- palette (Okabe-Ito, colour-blind safe) ------------------------------
OKABE_ITO = ["#0072B2", "#D55E00", "#009E73", "#CC79A7",
             "#E69F00", "#56B4E9", "#F0E442", "#000000"]
# sequential + diverging for heatmaps
SEQ = "mako"
DIVERGING = "vlag"
# semantic
INK = "#1a1a1a"
MUTED = "#5f6b7a"
GRID = "#d9dde3"
ACCENT = "#0072B2"
GOOD = "#009E73"
WARN = "#D55E00"


def apply():
    matplotlib.rcParams.update({
        "figure.dpi": 120,
        "savefig.dpi": 400,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,

        "font.family": "sans-serif",
        "font.sans-serif": [_SANS, "DejaVu Sans"],
        "font.size": 8.5,
        "axes.titlesize": 9.5,
        "axes.titleweight": "bold",
        "axes.labelsize": 8.5,
        "axes.labelweight": "regular",
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7.5,
        "figure.titlesize": 11,
        "figure.titleweight": "bold",

        "axes.edgecolor": "#3a3f47",
        "axes.linewidth": 0.8,
        "axes.labelcolor": INK,
        "axes.titlecolor": INK,
        "text.color": INK,
        "xtick.color": "#3a3f47",
        "ytick.color": "#3a3f47",
        "xtick.labelcolor": INK,
        "ytick.labelcolor": INK,

        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "axes.grid.axis": "both",
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "grid.alpha": 0.9,
        "axes.axisbelow": True,

        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,

        "lines.linewidth": 1.6,
        "lines.markersize": 5.0,
        "lines.markeredgewidth": 0.6,
        "scatter.edgecolors": "white",

        "legend.frameon": True,
        "legend.framealpha": 0.92,
        "legend.edgecolor": "#c7ccd3",
        "legend.borderpad": 0.5,
        "legend.handlelength": 1.4,

        "axes.prop_cycle": plt.cycler(color=OKABE_ITO),
        "image.cmap": "mako",
        "errorbar.capsize": 2.5,
        "pdf.fonttype": 42,   # embed as TrueType, editable text
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })
    try:
        import seaborn as sns
        sns.set_theme(context="paper", style="ticks", font=_SANS,
                      rc=matplotlib.rcParams)
        sns.set_palette(OKABE_ITO)
    except Exception:
        pass


def panel_label(ax, letter, dx=-0.02, dy=1.04, size=12):
    ax.text(dx, dy, f"({letter})", transform=ax.transAxes, ha="right",
            va="bottom", fontsize=size, fontweight="bold", color=INK)


def finish(ax):
    """Small consistent tidy-ups after plotting on an Axes."""
    ax.tick_params(length=3, width=0.8)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#3a3f47")
    return ax


def save(fig, path, also_pdf=True, dpi=400):
    """Save PNG (raster) and, when feasible, a vector PDF alongside it."""
    d = os.path.dirname(os.path.abspath(path))
    if d:
        os.makedirs(d, exist_ok=True)
    fig.savefig(path, dpi=dpi)
    if also_pdf:
        root, _ = os.path.splitext(path)
        try:
            fig.savefig(root + ".pdf")
        except Exception:
            pass
    plt.close(fig)
    return path
