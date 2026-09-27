"""
style.py - house style for every figure of the manuscript.

Sized for Springer (Journal of Molecular Modeling) artwork rules: figures are
drawn at final print size (single column 84 mm, double column 174 mm, max
height 234 mm) so 7-8 pt Arial lettering is exactly 7-8 pt on the page. No
figure titles are baked into the images - the caption lives in the manuscript.

Outputs per figure: vector PDF + 600 dpi PNG (+ 600 dpi LZW TIFF for upload).
"""
from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

MM = 1 / 25.4
SINGLE = 84 * MM
ONEHALF = 129 * MM
DOUBLE = 174 * MM - 0.04                # Springer double column, minus the 2 x 0.02 in save padding

_AVAIL = {f.name for f in fm.fontManager.ttflist}
SANS = next((n for n in ("Arial", "Helvetica", "Liberation Sans", "DejaVu Sans")
             if n in _AVAIL), "DejaVu Sans")

# ---- colours --------------------------------------------------------------
INK = "#1d2330"
MUTED = "#6b7280"
FAINT = "#e5e7eb"
GRID = "#eceef1"
PANEL_BG = "#f7f8fa"

# endpoints (kept identical in every figure)
DOCK = "#2166ac"        # docking endpoint
ADS = "#1b7837"         # adsorption endpoint
# adsorption regimes
PHYS = "#4d9de0"
CHEM = "#d1495b"
# carriers
PRISTINE = "#5b6770"
DOPED = "#b2182b"
# categorical groups (Paul Tol "muted", colour-blind safe), assigned in order
GROUPS = ["#882255", "#ddaa33", "#44aa99", "#332288", "#117733", "#cc6677"]


def apply():
    matplotlib.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 600,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,

        "font.family": "sans-serif",
        "font.sans-serif": [SANS, "DejaVu Sans"],
        "mathtext.fontset": "custom",
        "mathtext.rm": SANS,
        "mathtext.it": f"{SANS}:italic",
        "mathtext.bf": f"{SANS}:bold",
        "mathtext.default": "regular",
        "font.size": 7.5,
        "axes.titlesize": 7.5,
        "axes.titleweight": "bold",
        "axes.titlepad": 4,
        "axes.titlelocation": "left",
        "axes.labelsize": 7.5,
        "axes.labelpad": 2.5,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "legend.fontsize": 6.5,
        "legend.title_fontsize": 6.5,

        "axes.edgecolor": INK,
        "axes.linewidth": 0.6,
        "axes.labelcolor": INK,
        "axes.titlecolor": INK,
        "text.color": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "grid.color": GRID,
        "grid.linewidth": 0.5,
        "axes.axisbelow": True,

        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "xtick.minor.size": 1.5,
        "ytick.minor.size": 1.5,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.pad": 2,
        "ytick.major.pad": 2,

        "lines.linewidth": 1.0,
        "lines.markersize": 4,
        "patch.linewidth": 0.5,

        "legend.frameon": False,
        "legend.borderaxespad": 0.3,
        "legend.handlelength": 1.2,
        "legend.handletextpad": 0.4,
        "legend.labelspacing": 0.3,
        "legend.columnspacing": 1.0,

        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def panel(ax, letter, x=-0.02, y=1.0, fig=None):
    """Bold lower-case panel letter at the top-left corner of an axes."""
    ax.text(x, y, letter, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=9, fontweight="bold", color=INK)


def fig_panel(fig, letter, x, y):
    fig.text(x, y, letter, ha="left", va="top", fontsize=9,
             fontweight="bold", color=INK)


def scatter_kw(color, size=18, lw=0.4, edge="white", alpha=1.0, zorder=3):
    return dict(s=size, color=color, edgecolor=edge, linewidth=lw,
                alpha=alpha, zorder=zorder)


MIN_PT = 7.0          # Springer: lettering about 2-3 mm at final size; 7 pt = 2.5 mm


def finish(fig):
    """Last pass over every text of a figure: true minus signs (U+2212) for negative
    numbers, and no lettering below MIN_PT."""
    import re
    from matplotlib.text import Text
    for t in fig.findobj(Text):
        s = t.get_text()
        if s:
            s2 = re.sub(r"(?<![A-Za-z0-9_$\\{])-(?=\d)", "\u2212", s)
            if s2 != s:
                t.set_text(s2)
        if t.get_fontsize() < MIN_PT and s:
            t.set_fontsize(MIN_PT)


def save(fig, out_dir, stem, tiff=True):
    """Write <stem>.pdf (vector), <stem>.png and <stem>.tif (600 dpi)."""
    finish(fig)
    os.makedirs(out_dir, exist_ok=True)
    base = os.path.join(out_dir, stem)
    fig.savefig(base + ".pdf")
    fig.savefig(base + ".png", dpi=600)
    if tiff:
        fig.savefig(base + ".tif", dpi=600,
                    pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    return base + ".png"


def imshow_clean(ax, img):
    ax.imshow(img, interpolation="lanczos")
    ax.set_axis_off()


def autocrop(path, pad=12, out=None):
    """Trim uniform white/transparent borders of a render in place."""
    from PIL import Image, ImageChops
    im = Image.open(path)
    if im.mode == "RGBA":
        bg = Image.new("RGBA", im.size, (255, 255, 255, 0))
        diff = ImageChops.difference(im, bg).split()[-1]
        bbox = diff.getbbox()
    else:
        rgb = im.convert("RGB")
        bbox = ImageChops.difference(rgb, Image.new("RGB", im.size, "white")).getbbox()
    if bbox:
        l, t, r, b = bbox
        im = im.crop((max(0, l - pad), max(0, t - pad),
                      min(im.width, r + pad), min(im.height, b + pad)))
    im.save(out or path)
    return out or path
