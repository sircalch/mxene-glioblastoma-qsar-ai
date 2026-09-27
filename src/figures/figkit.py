"""
figkit.py - figure building blocks shared by the manuscript figures.

Every function draws into axes it is given (or builds a complete figure) from
plain pandas/numpy data, in the house style of style.py. Nothing here reads
project files; make_figures.py wires data to panels.
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch

import style as S


# ------------------------------------------------------------------ helpers
def label_points(ax, x, y, labels, fontsize=5.8, color=None, offsets=None, arrows=True):
    """Place short text labels next to points without overlapping each other."""
    try:
        from adjustText import adjust_text
    except ImportError:
        adjust_text = None
    texts = []
    for i, (xi, yi, t) in enumerate(zip(x, y, labels)):
        dx, dy = (offsets or {}).get(t, (0, 0))
        texts.append(ax.text(xi + dx, yi + dy, t, fontsize=fontsize, color=color or S.INK,
                             ha="left", va="bottom", zorder=6,
                             path_effects=[pe.withStroke(linewidth=1.8, foreground="white")]))
    if adjust_text and texts:
        adjust_text(texts, x=list(x), y=list(y), ax=ax, expand=(1.15, 1.3))
    if arrows:
        # leader lines as annotations: they are anchored to the data point and are
        # redrawn with the final layout, so they always end at their point
        out = []
        for t, xi, yi in zip(texts, x, y):
            tx, ty = t.get_position()
            lab = ax.annotate(t.get_text(), xy=(xi, yi), xytext=(tx, ty), fontsize=t.get_fontsize(),
                              color=t.get_color(), ha=t.get_ha(), va=t.get_va(), zorder=6,
                              path_effects=[pe.withStroke(linewidth=1.8, foreground="white")],
                              arrowprops=dict(arrowstyle="-", color=S.MUTED, lw=0.5, shrinkA=1.5, shrinkB=3))
            t.remove()
            out.append(lab)
        texts = out
    return texts


def stat_box(ax, lines, loc="upper left"):
    x, ha = (0.04, "left") if "left" in loc else (0.96, "right")
    y, va = (0.96, "top") if "upper" in loc else (0.04, "bottom")
    ax.text(x, y, "\n".join(lines), transform=ax.transAxes, ha=ha, va=va, fontsize=6.3,
            linespacing=1.35, color=S.INK,
            bbox=dict(boxstyle="round,pad=0.35,rounding_size=0.15", fc="white", ec=S.FAINT, lw=0.6))


def light_grid(ax, axis="both"):
    ax.grid(True, axis=axis, color=S.GRID, lw=0.5, zorder=0)


# ------------------------------------------------------------------ parity / QSPR
def log_ticks(ax, axis="y"):
    """Readable ticks (5, 10, 20, 50 ...) on a log axis instead of 10^n."""
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
    a = ax.yaxis if axis == "y" else ax.xaxis
    lo, hi = (ax.get_ylim() if axis == "y" else ax.get_xlim())
    a.set_major_locator(FixedLocator([t for t in (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000) if lo <= t <= hi]))
    a.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    a.set_minor_formatter(NullFormatter())


def parity(ax, y, yhat, color, n_label=None, stats_lines=None, unit="kcal mol$^{-1}$", names=None,
           label_extremes=3, xlabel="GFN2-xTB value"):
    lo = min(y.min(), yhat.min())
    hi = max(y.max(), yhat.max())
    pad = 0.08 * (hi - lo)
    lo, hi = lo - pad, hi + pad
    ax.fill_between([lo, hi], [lo - 0.1 * (hi - lo), hi - 0.1 * (hi - lo)],
                    [lo + 0.1 * (hi - lo), hi + 0.1 * (hi - lo)], color=S.PANEL_BG, zorder=0, lw=0)
    ax.plot([lo, hi], [lo, hi], color=S.MUTED, lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.scatter(y, yhat, **S.scatter_kw(color, size=22, edge="white", lw=0.5))
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.set_xlabel(f"{xlabel} ({unit})")
    ax.set_ylabel(f"Out-of-fold prediction ({unit})")
    if stats_lines:
        stat_box(ax, stats_lines)
    if names is not None and label_extremes:
        err = np.abs(np.asarray(yhat) - np.asarray(y))
        idx = np.argsort(err)[-label_extremes:]
        label_points(ax, np.asarray(y)[idx], np.asarray(yhat)[idx], np.asarray(names)[idx])


def williams(ax, h, r, h_star, color, names=None):
    inside = (h <= h_star) & (np.abs(r) <= 3)
    ax.axhspan(-3, 3, xmin=0, xmax=1, color=S.PANEL_BG, zorder=0, lw=0)
    ax.axhline(3, color=S.CHEM, lw=0.7, ls=(0, (4, 3)))
    ax.axhline(-3, color=S.CHEM, lw=0.7, ls=(0, (4, 3)))
    ax.axhline(0, color=S.MUTED, lw=0.5)
    ax.axvline(h_star, color=S.INK, lw=0.7, ls=(0, (1, 2)))
    ax.scatter(h[inside], r[inside], **S.scatter_kw(color, size=20))
    ax.scatter(h[~inside], r[~inside], **S.scatter_kw("white", size=22, edge=S.CHEM, lw=0.9))
    top = max(4, np.abs(r).max() * 1.15)
    ax.set_ylim(-top, top)
    ax.set_xlim(0, max(h.max(), h_star) * 1.12)
    ax.text(h_star, -top * 0.97, f" $h^*$ = {h_star:.2f}", fontsize=6.2, va="bottom", ha="left", color=S.INK)
    ax.set_xlabel("Leverage $h_i$")
    ax.set_ylabel("Standardised residual")
    if names is not None and (~inside).any():
        label_points(ax, h[~inside], r[~inside], np.asarray(names)[~inside], color=S.CHEM)


def scrambling(ax, perm, q2, color):
    ax.hist(perm, bins=30, color=S.FAINT, edgecolor="white", lw=0.4)
    ax.axvline(q2, color=color, lw=1.4)
    lo, hi = ax.get_xlim()
    below = q2 < np.median(perm)
    if below:
        ax.set_xlim(min(lo, q2 - 0.35), hi)
        txt, ha = f"model  \n$Q^2_{{CV}}$ = {q2:.2f}  ", "right"
    else:
        ax.set_xlim(lo, max(hi, q2 + 0.45))
        txt, ha = f"  model\n  $Q^2_{{CV}}$ = {q2:.2f}", "left"
    ax.text(q2, ax.get_ylim()[1] * 0.95, txt, color=color, fontsize=6.3, va="top", ha=ha)
    ax.set_xlabel("$Q^2_{CV}$ with permuted target")
    ax.set_ylabel("Permutations")


# ------------------------------------------------------------------ docking
def ranked_dots(ax, names, values, colors, xlabel, highlight=None):
    order = np.argsort(values)[::-1]          # most negative at the top
    y = np.arange(len(values))
    v = np.asarray(values)[order]
    ax.hlines(y, v, v.max() + 0.3, color=S.FAINT, lw=0.6, zorder=1)
    ax.scatter(v, y, s=16, c=np.asarray(colors)[order], edgecolor="white", lw=0.4, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(np.asarray(names)[order], fontsize=5.6)
    if highlight:
        for t in ax.get_yticklabels():
            if t.get_text() in highlight:
                t.set_fontweight("bold")
    ax.set_xlabel(xlabel)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    light_grid(ax, "x")
    ax.invert_xaxis()


def group_strip(ax, groups, values, colors, ylabel, order=None, seed=0):
    rng = np.random.default_rng(seed)
    order = order or list(dict.fromkeys(groups))
    for i, g in enumerate(order):
        v = np.asarray(values)[np.asarray(groups) == g]
        box = ax.boxplot(v, positions=[i], widths=0.55, showfliers=False, patch_artist=True,
                         medianprops=dict(color=S.INK, lw=1.0), whiskerprops=dict(color=S.MUTED, lw=0.6),
                         capprops=dict(color=S.MUTED, lw=0.6), boxprops=dict(lw=0.6, ec=S.MUTED))
        box["boxes"][0].set_facecolor(colors[i] + "26")
        ax.scatter(i + rng.uniform(-0.16, 0.16, len(v)), v, **S.scatter_kw(colors[i], size=15))
    ax.set_xticks(range(len(order)))
    ax.set_xlim(-0.6, len(order) - 0.4)
    ax.set_ylabel(ylabel)
    light_grid(ax, "y")


def contact_map(ax, matrix, drugs, residues, cmap_color, vmax=None):
    """Binary/count drug x residue contact matrix."""
    from matplotlib.colors import LinearSegmentedColormap
    cm = LinearSegmentedColormap.from_list("c", ["#ffffff", cmap_color])
    ax.imshow(matrix, aspect="auto", cmap=cm, vmin=0, vmax=vmax or matrix.max(), interpolation="nearest")
    ax.set_xticks(range(len(residues)))
    ax.set_xticklabels(residues, rotation=90, fontsize=5.4)
    ax.set_yticks(range(len(drugs)))
    ax.set_yticklabels(drugs, fontsize=5.4)
    ax.set_xticks(np.arange(-0.5, len(residues)), minor=True)
    ax.set_yticks(np.arange(-0.5, len(drugs)), minor=True)
    ax.grid(which="minor", color="white", lw=0.6)
    ax.tick_params(which="both", length=0)
    for s in ax.spines.values():
        s.set_visible(False)


# ------------------------------------------------------------------ adsorption
def landscape(ax, contact, energy, chem_mask, names=None, label="chem", max_labels=8, ylog=True):
    """Closest drug-carrier contact vs -dE_int, covalent-contact band shaded.
    label: 'chem' labels chemisorbers, 'all' every point, None nothing (at most
    max_labels, strongest binders first)."""
    from matplotlib.ticker import FixedLocator, NullFormatter, FuncFormatter
    x = np.asarray(contact)
    e = -np.asarray(energy)
    ax.axvspan(1.2, 1.9, color=S.CHEM, alpha=0.07, lw=0, zorder=0)
    ax.text(1.55, 0.02, "covalent\ncontact", transform=ax.get_xaxis_transform(), ha="center",
            va="bottom", fontsize=5.8, color=S.CHEM)
    ax.scatter(x[~chem_mask], e[~chem_mask], **S.scatter_kw(S.PHYS, size=20))
    ax.scatter(x[chem_mask], e[chem_mask], **S.scatter_kw(S.CHEM, size=20))
    if ylog:
        ax.set_yscale("log")
        lo, hi = e.min() * 0.8, e.max() * 1.25
        ticks = [t for t in (1, 2, 3, 5, 10, 20, 30, 50, 75, 100, 150, 200, 300, 500) if lo <= t <= hi]
        ax.set_ylim(lo, hi)
        ax.yaxis.set_major_locator(FixedLocator(ticks))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("Closest drug–carrier contact (Å)")
    ax.set_ylabel("−Δ$E_{int}$ (kcal mol$^{-1}$)")
    if names is not None and label:
        mask = chem_mask if label == "chem" else np.ones_like(chem_mask)
        idx = np.where(mask)[0]
        idx = idx[np.argsort(-e[idx])][:max_labels]
        label_points(ax, x[idx], e[idx], np.asarray(names)[idx], fontsize=5.5,
                     color=S.CHEM if label == "chem" else S.INK)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", mfc=S.PHYS, mec="white", ms=4.5, label="physisorption"),
                       Line2D([], [], marker="o", ls="", mfc=S.CHEM, mec="white", ms=4.5, label="chemisorption")],
              loc="upper right", frameon=False, borderaxespad=0.1)


# ------------------------------------------------------------------ renders
def render_panel(ax, png, title=None, sub=None):
    import matplotlib.image as mpimg
    ax.imshow(mpimg.imread(png), interpolation="lanczos")
    ax.set_axis_off()
    ax.set_box_aspect(0.9)
    if title:
        ax.set_title(title, fontsize=7, fontweight="bold", loc="center", pad=2)
    if sub:
        ax.text(0.5, -0.02, sub, transform=ax.transAxes, ha="center", va="top", fontsize=6.2, color=S.MUTED)


def residue_labels(ax, positions, img_shape, fontsize=6.0, color=S.INK, pad_frac=0.035):
    """Typeset residue labels at projected pixel positions (from render3d)."""
    H, W = img_shape[:2]
    for name, (x, y) in positions.items():
        if not (0 <= x <= W and 0 <= y <= H):
            continue
        ax.text(x, y - pad_frac * H, name, fontsize=fontsize, ha="center", va="bottom", color=color,
                path_effects=[pe.withStroke(linewidth=2.2, foreground="white")], zorder=5)


def legend_row(fig, items, y=0.01, fontsize=6.3):
    """items: [(kind, colour, text)], kind in {'dot', 'line', 'dash', 'patch'}."""
    handles = []
    for kind, col, txt in items:
        if kind == "dot":
            handles.append(Line2D([], [], marker="o", ls="", mfc=col, mec="white", ms=5, label=txt))
        elif kind == "line":
            handles.append(Line2D([], [], color=col, lw=2.2, label=txt))
        elif kind == "dash":
            handles.append(Line2D([], [], color=col, lw=1.0, ls=(0, (2, 1.5)), label=txt))
        else:
            from matplotlib.patches import Patch
            handles.append(Patch(fc=col, ec="none", label=txt))
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), frameon=False,
               bbox_to_anchor=(0.5, y), fontsize=fontsize, handletextpad=0.4, columnspacing=1.4)


# ------------------------------------------------------------------ workflow
def workflow(fig, stages, outcome, images=None, accent="#2166ac"):
    """Horizontal pipeline of cards. stages: [(title, [lines])]; images: {i: png}."""
    import matplotlib.image as mpimg
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    n = len(stages)
    gap, left, right = 0.022, 0.012, 0.80
    w = (right - left - gap * (n - 1)) / n
    for i, (title, lines) in enumerate(stages):
        x = left + i * (w + gap)
        ax.add_patch(FancyBboxPatch((x, 0.06), w, 0.88, boxstyle="round,pad=0,rounding_size=0.018",
                                    fc="white", ec="#d5dae1", lw=0.7, transform=ax.transAxes))
        ax.add_patch(FancyBboxPatch((x, 0.90), w, 0.04, boxstyle="square,pad=0", fc=accent, ec="none",
                                    transform=ax.transAxes, alpha=0.9))
        ax.text(x + 0.012, 0.855, f"{i + 1}", fontsize=10, fontweight="bold", color=accent, va="top",
                transform=ax.transAxes)
        ax.text(x + 0.034, 0.853, title, fontsize=7, fontweight="bold", color=S.INK, va="top",
                transform=ax.transAxes, wrap=True)
        if images and i in images:
            img = mpimg.imread(images[i])
            ih = 0.36
            iw = ih * img.shape[1] / img.shape[0] * fig.get_figheight() / fig.get_figwidth()
            iw = min(iw, w - 0.02)
            ih = iw / (img.shape[1] / img.shape[0] * fig.get_figheight() / fig.get_figwidth())
            iax = fig.add_axes([x + (w - iw) / 2, 0.40, iw, ih])
            iax.imshow(img, interpolation="lanczos")
            iax.set_axis_off()
        ax.text(x + 0.012, 0.35, "\n".join(lines), fontsize=6.0, color=S.INK, va="top", linespacing=1.45,
                transform=ax.transAxes)
        if i < n - 1:
            ax.annotate("", xy=(x + w + gap - 0.003, 0.5), xytext=(x + w + 0.003, 0.5),
                        xycoords=ax.transAxes, arrowprops=dict(arrowstyle="-|>", color=S.MUTED, lw=0.9,
                                                               mutation_scale=7))
    xo = right + gap
    ax.annotate("", xy=(xo - 0.003, 0.5), xytext=(right + 0.003, 0.5), xycoords=ax.transAxes,
                arrowprops=dict(arrowstyle="-|>", color=S.MUTED, lw=0.9, mutation_scale=7))
    ax.add_patch(FancyBboxPatch((xo, 0.06), 1 - xo - 0.01, 0.88, boxstyle="round,pad=0,rounding_size=0.018",
                                fc=S.PANEL_BG, ec="#d5dae1", lw=0.7, transform=ax.transAxes))
    ax.text(xo + 0.012, 0.86, outcome[0], fontsize=7, fontweight="bold", va="top", color=accent,
            transform=ax.transAxes)
    ax.text(xo + 0.012, 0.76, "\n".join(outcome[1]), fontsize=6.0, va="top", linespacing=1.5,
            color=S.INK, transform=ax.transAxes)
    return ax
