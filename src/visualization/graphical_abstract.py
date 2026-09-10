"""
graphical_abstract.py
Composed graphical abstract -> figures/fig1_graphical_abstract.png
(title bar + real 3D render + inline results panel + KEY FINDING band).
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.ticker

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
KEY = 'gbm'
CFG = {'color': '#1565C0', 'accent': '#00838F', 'title': 'Explainable AI + Quantum Screening of 2D Ti$_3$C$_2$T$_x$ MXene for Glioblastoma', 'hero': 'figures/fig9_gbm_3d_spatial_binding_modes.png', 'hero_crop': (0.34, 0.66), 'hero_cap': 'Kinase inhibitor physisorbed on the pristine Ti$_3$C$_2$O$_2$ MXene surface', 'data': 'data/processed/dataset_drug_mxene_pristine.csv', 'mode': 'xy', 'x': 'vina_4ZAU_kcal_mol', 'y': 'delta_Eint_SP_kcal_mol', 'panel_title': 'EGFR docking vs. MXene adsorption (n = 33)', 'xlab': 'AutoDock Vina on EGFR 4ZAU (kcal/mol)', 'ylab': 'ΔE$_{int,SP}$ on Ti$_3$C$_2$O$_2$ (kcal/mol)', 'takeaway': 'All 35 drugs physisorb on pristine Ti$_3$C$_2$O$_2$ (−0.9 to −15.5 kcal/mol); EGFR affinity and surface adsorption are only weakly (negatively) correlated. The descriptor QSPR is at best weakly predictive (Q$^2_{CV}$ = 0.31 docking / 0.10 adsorption).'}

def build():
    key, c = KEY, CFG
    fig = plt.figure(figsize=(12.0, 5.6), dpi=300)
    fig.patch.set_facecolor("white")

    # ---- title bar
    tb = fig.add_axes([0, 0.90, 1, 0.10]); tb.axis("off")
    tb.add_patch(plt.Rectangle((0, 0), 1, 1, transform=tb.transAxes, color=c["color"]))
    tb.text(0.5, 0.5, c["title"], transform=tb.transAxes, ha="center", va="center",
            color="white", fontsize=13.5, fontweight="bold")

    # ---- hero 3D image (left) -- crop out fig9's own panel title/caption
    ax_h = fig.add_axes([0.02, 0.16, 0.45, 0.70]); ax_h.axis("off")
    img = mpimg.imread(os.path.join(BASE, c["hero"]))
    h, w = img.shape[:2]
    lo, hi = c["hero_crop"]
    ax_h.imshow(img[int(h*0.21):int(h*0.84), int(w*lo):int(w*hi)])
    ax_h.text(0.02, 1.02, "a", transform=ax_h.transAxes, fontsize=13, fontweight="bold",
              va="bottom", color=c["color"])
    # caption sits inside the image frame, bottom-left, on a translucent plate
    ax_h.text(0.5, 0.015, c["hero_cap"], transform=ax_h.transAxes, ha="center", va="bottom",
              fontsize=8.3, style="italic", color="#222222",
              bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.75", alpha=0.9))

    # ---- results panel (right)
    ax_p = fig.add_axes([0.58, 0.27, 0.38, 0.53])
    df = pd.read_csv(os.path.join(BASE, c["data"]))
    df = df.dropna(subset=[c["x"], c["y"]])
    if c["mode"] == "regime":
        for lab, col in [("physisorption", "#1565C0"), ("chemisorption", c["accent"])]:
            s = df[df["adsorption_mode"] == lab]
            ax_p.scatter(s[c["x"]], s[c["y"]], s=42, c=col, edgecolor="k", lw=0.4,
                         label=f"{lab} (n={len(s)})", zorder=3)
        ax_p.axvspan(1.2, 1.9, color="0.90", zorder=0)
        ax_p.set_yscale("symlog")
        ax_p.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
            lambda v, p: "" if v == 0 else f"{int(v)}"))
        ax_p.set_yticks([-10, -30, -100, -200])
        ax_p.legend(fontsize=7.5, loc="upper left", framealpha=0.95)
    else:
        ax_p.scatter(df[c["x"]], df[c["y"]], s=42, c=c["color"], edgecolor="k", lw=0.4, zorder=3)
        r = np.corrcoef(df[c["x"]], df[c["y"]])[0, 1]
        ax_p.text(0.05, 0.93, f"Pearson r = {r:.2f}", transform=ax_p.transAxes, fontsize=8.5,
                  va="top", bbox=dict(boxstyle="round", fc="white", ec="0.7"))
    ax_p.set_title(c["panel_title"], fontsize=9.5, fontweight="bold")
    ax_p.set_xlabel(c["xlab"], fontsize=8.5)
    ax_p.set_ylabel(c["ylab"], fontsize=8.5)
    ax_p.tick_params(labelsize=7.5)
    ax_p.grid(alpha=0.25)
    ax_p.text(-0.28, 1.10, "b", transform=ax_p.transAxes, fontsize=13, fontweight="bold",
              va="top", color=c["color"])

    # ---- takeaway strip (distinct band: solid fill + coloured top rule + KEY label)
    st = fig.add_axes([0, 0, 1, 0.135]); st.axis("off")
    st.add_patch(plt.Rectangle((0, 0), 1, 1, transform=st.transAxes, fc=c["color"], alpha=0.12, ec="none"))
    st.add_patch(plt.Rectangle((0, 0.92), 1, 0.08, transform=st.transAxes, fc=c["color"], ec="none"))
    st.text(0.012, 0.5, "KEY\nFINDING", transform=st.transAxes, ha="left", va="center",
            fontsize=8.5, fontweight="bold", color=c["color"])
    st.text(0.53, 0.46, c["takeaway"], transform=st.transAxes, ha="center", va="center",
            fontsize=8.7, color="#111111", wrap=True)

    out = os.path.join(BASE, "figures", "fig1_graphical_abstract.png")
    fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)



if __name__ == "__main__":
    build()
