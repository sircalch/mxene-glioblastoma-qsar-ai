"""
make_figures.py - every figure of the GBM / Ti3C2O2 MXene manuscript, drawn from
the pipeline outputs at Springer print size (see style.py).

usage: python src/figures/make_figures.py [fig ...]      (default: all)
writes figures/FigN.{pdf,png,tif}; 3D renders are cached in figures/_renders
Fig. 5 (DFT adsorption) is drawn only when results/quantum/dft_adsorption.csv
holds completed calculations.
"""
import json
import sys
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import figkit as K  # noqa: E402
import render3d as R  # noqa: E402
import style as S  # noqa: E402

BASE = HERE.parents[1]
FIG = BASE / "figures"
REN = FIG / "_renders"
POSES = BASE / "results" / "docking" / "real_poses"
DFT = BASE / "calculations" / "gbm_dft"

FAMILY = {"Alkylating agents": S.GROUPS[0], "EGFR inhibitors": S.GROUPS[1],
          "Other kinase inhibitors": S.GROUPS[2], "Other targeted agents": S.GROUPS[3]}
SHORT = {"Alkylating agents": "Alkylating", "EGFR inhibitors": "EGFR", "Other kinase inhibitors": "Other\nkinase",
         "Other targeted agents": "Other\ntargeted"}
OTHER = {"mTOR Inhibitor", "Proteasome Inhibitor", "HDAC Inhibitor"}


def family(cls):
    return ("Alkylating agents" if cls == "Alkylating Agent" else "EGFR inhibitors" if cls.startswith("EGFR")
            else "Other targeted agents" if cls in OTHER else "Other kinase inhibitors")


def data():
    lib = pd.read_csv(BASE / "data" / "processed" / "compound_library_pubchem.csv")
    dock = pd.read_csv(BASE / "results" / "docking" / "real_vina_docking_summary.csv")
    m = lib[["name", "class", "smiles", "n_heavy_atoms"]].merge(
        dock[["name", "vina_1M17_kcal_mol", "ligand_efficiency"]], on="name")
    m["family"] = m["class"].map(family)
    m["color"] = m.family.map(FAMILY)
    d = {"m": m, "docked": m.dropna(subset=["vina_1M17_kcal_mol"]).reset_index(drop=True)}
    d["redock"] = pd.read_csv(BASE / "results" / "docking" / "redocking_validation.csv")
    d["modes"] = pd.read_csv(BASE / "results" / "docking" / "redocking_all_modes.csv")
    d["contacts"] = pd.read_csv(BASE / "results" / "docking" / "residue_contacts.csv")
    q = BASE / "results" / "qspr"
    d["q"] = json.loads((q / "vina_summary.json").read_text())
    d["oof"] = pd.read_csv(q / "vina_oof.csv")
    d["perm"] = pd.read_csv(q / "vina_y_scrambling.csv")
    f = BASE / "results" / "quantum" / "xtb_lattice_scan.csv"
    d["scan"] = pd.read_csv(f) if f.exists() else None
    f = BASE / "results" / "quantum" / "dft_adsorption.csv"
    d["dft"] = pd.read_csv(f) if f.exists() else None
    return d


def render(name, fn, *a, **kw):
    REN.mkdir(parents=True, exist_ok=True)
    png, meta = REN / f"{name}.png", REN / f"{name}.json"
    if not png.exists():
        out = fn(*a, out_png=str(png), **kw)
        meta.write_text(json.dumps(out if isinstance(out, dict) else {}))
        S.autocrop(png)
    return png, json.loads(meta.read_text()) if meta.exists() else {}


def family_legend(fig, y=0.0):
    K.legend_row(fig, [("dot", c, f) for f, c in FAMILY.items()], y=y, fontsize=6.2)


def model_pdb(pdbqt, out_pdb, k):
    """MODEL k (1-based) of a Vina .pdbqt -> plain PDB."""
    blocks = Path(pdbqt).read_text().split("ENDMDL")
    tmp = Path(out_pdb).with_suffix(".tmp.pdbqt")
    tmp.write_text(blocks[k - 1].split("MODEL")[-1] + "ENDMDL\n")
    R.first_model_pdb(str(tmp), str(out_pdb))
    tmp.unlink()
    return out_pdb


def best_mode(d, control):
    x = d["modes"][d["modes"].control == control]
    return x.loc[x.rmsd_A.idxmin()]


def pocket_render(d):
    REN.mkdir(parents=True, exist_ok=True)
    ctl = "self-redock, crystal conformation"
    bm = best_mode(d, ctl)
    top, good = REN / "redock_top.pdb", REN / "redock_best.pdb"
    model_pdb(POSES / "AQ4_redock_out.pdbqt", top, 1)
    model_pdb(POSES / "AQ4_redock_out.pdbqt", good, int(bm["mode"]))
    return render("egfr_pocket", R.pocket_closeup, str(BASE / "data" / "raw" / "1M17.pdb"),
                  [{"sel": "resn AQ4", "color": (0.70, 0.72, 0.75), "radius": 0.30, "name": "xtal"},
                   {"path": str(top), "color": (0.91, 0.55, 0.16), "radius": 0.16, "name": "top"},
                   {"path": str(good), "color": (0.20, 0.62, 0.60), "radius": 0.16, "name": "best"}],
                  size=(1600, 1300), slab=24, zoom=2.4)


# ------------------------------------------------------------------ Fig. 1
def fig1(d):
    """Workflow."""
    m, q = d["m"], d["q"]
    pocket, _ = pocket_render(d)
    slab, _ = render("slab_top", R.molecule, str(slab_for_render()), tilt=0, size=(1300, 1000))
    mol2d = REN / "tmz_2d.png"
    if not mol2d.exists():
        from rdkit import Chem
        from rdkit.Chem import Draw
        Draw.MolToFile(Chem.MolFromSmiles(m.set_index("name").loc["Temozolomide", "smiles"]), str(mol2d),
                       size=(600, 460))
        S.autocrop(mol2d)
    mini = REN / "qspr_mini.png"
    o = d["oof"]
    fm, am = plt.subplots(figsize=(1.2, 1.2))
    lo, hi = min(o.vina.min(), o.oof_pred.min()), max(o.vina.max(), o.oof_pred.max())
    am.plot([lo, hi], [lo, hi], color=S.MUTED, lw=0.8, ls=(0, (3, 2)))
    am.scatter(o.vina, o.oof_pred, s=9, color=S.DOCK, edgecolor="white", lw=0.3)
    am.set_xticks([]); am.set_yticks([]); am.set_xlabel("observed", fontsize=6); am.set_ylabel("predicted", fontsize=6)
    fm.savefig(mini, dpi=400, bbox_inches="tight"); plt.close(fm)
    ndft = 0 if d["dft"] is None else int(d["dft"].get("delta_Eads_kcal_mol", pd.Series(dtype=float)).notna().sum())
    stages = [
        ("Drug set", [f"{len(m)} GBM drugs", "4 families", "from PubChem"]),
        ("EGFR docking", ["PDB 1M17", "Vina 1.2.7", "all-mode RMSD"]),
        ("Ti$_3$C$_2$O$_2$ MXene", ["periodic 4×4 slab", "a = 3.03 Å", "xTB tested: fails"]),
        ("DFT adsorption", ["PBE-D3 (QE 7.5)", "4 alkylating agents", "relaxed complexes"]),
        ("QSPR", ["Vina score, ridge", "nested 5×5 CV", "Y-scrambling, AD"]),
    ]
    imgs = {0: mol2d, 1: pocket, 2: slab, 4: mini}
    tmz = DFT / "cplx_Temozolomide" / "final.xyz"
    if tmz.exists():                                    # DFT card: the relaxed temozolomide complex
        imgs[3], _ = render("dft_Temozolomide", R.molecule, str(tmz), tilt=75, size=(1200, 1000))
    fig = plt.figure(figsize=(S.DOUBLE, 58 * S.MM))
    K.workflow(fig, stages, ("Outcome", [f"{len(d['docked'])} drugs docked", f"{ndft} DFT complexes",
                                         f"Vina $Q^2_{{CV}}$ = {q['Q2_CV']:.2f}"]),
               images=imgs, accent=S.DOCK)
    S.save(fig, FIG, "Fig1")


# ------------------------------------------------------------------ Fig. 2
def fig2(d):
    """Docking controls: pocket, RMSD of every mode, contact residues."""
    png, _ = pocket_render(d)
    md = d["modes"]
    ct = d["contacts"]
    n = ct.name.nunique()
    freq = ct.groupby("residue").name.nunique().sort_values(ascending=False).head(12)
    polar = ct[ct.polar].groupby("residue").name.nunique().reindex(freq.index).fillna(0)
    fig = plt.figure(figsize=(S.DOUBLE, 70 * S.MM))
    gs = GridSpec(1, 3, width_ratios=[1.25, 0.9, 0.85], wspace=0.45, left=0.01, right=0.99, top=0.93, bottom=0.22)
    ax = fig.add_subplot(gs[0])
    ax.imshow(mpimg.imread(png))
    ax.set_axis_off()
    S.panel(ax, "a", x=0.02, y=0.97)
    ctl = "self-redock, crystal conformation"
    bm = best_mode(d, ctl)
    top = md[(md.control == ctl) & (md["mode"] == 1)].iloc[0]
    K.legend_row(fig, [("line", "#b3b8bf", "crystal erlotinib (AQ4)"),
                       ("line", "#e88c29", f"top-ranked pose ({top.rmsd_A:.1f} Å)"),
                       ("line", "#339e99", f"mode {int(bm['mode'])} ({bm.rmsd_A:.2f} Å)")], y=0.0, fontsize=6)
    ax = fig.add_subplot(gs[1])
    for c, col, mk in (("self-redock, crystal conformation", S.DOCK, "o"),
                       ("production protocol, from SMILES", S.ADS, "s")):
        x = md[md.control == c]
        ax.scatter(x.vina_kcal_mol, x.rmsd_A, s=18, color=col, marker=mk, edgecolor="white", lw=0.4, zorder=3,
                   label="crystal conf." if "crystal" in c else "from SMILES")
    ax.axhline(2.0, color=S.CHEM, lw=0.7, ls=(0, (4, 3)))
    ax.text(ax.get_xlim()[0], 2.1, " 2 Å", fontsize=6, color=S.CHEM, va="bottom")
    ax.set_xlabel("Vina score of mode (kcal mol$^{-1}$)")
    ax.set_ylabel("RMSD to crystal pose (Å)")
    ax.invert_xaxis()
    ax.legend(loc="upper right", frameon=False, fontsize=6)
    K.light_grid(ax)
    S.panel(ax, "b", x=-0.3)
    ax2 = fig.add_subplot(gs[2])
    y = np.arange(len(freq))[::-1]
    ax2.barh(y, freq.values / n * 100, color=S.DOCK, alpha=0.25, height=0.7, label="any")
    ax2.barh(y, polar.values / n * 100, color=S.DOCK, height=0.7, label="polar")
    ax2.set_yticks(y)
    ax2.set_yticklabels(freq.index, fontsize=6.2)
    ax2.set_xlabel(f"Drugs in contact (% of {n})")
    ax2.set_xlim(0, 100)
    ax2.tick_params(axis="y", length=0)
    K.light_grid(ax2, "x")
    ax2.legend(loc="lower right", frameon=False, fontsize=6)
    S.panel(ax2, "c", x=-0.42)
    S.save(fig, FIG, "Fig2")


# ------------------------------------------------------------------ Fig. 3
def fig3(d):
    """Docking scores by drug and ligand efficiency by family."""
    m = d["docked"]
    fig, axs = plt.subplots(1, 2, figsize=(S.DOUBLE, 100 * S.MM),
                            gridspec_kw=dict(width_ratios=[1, 1], wspace=0.35))
    fig.subplots_adjust(bottom=0.14)
    K.ranked_dots(axs[0], m.name.values, m.vina_1M17_kcal_mol.values, m.color.values,
                  "Vina score (kcal mol$^{-1}$)", highlight={"Erlotinib", "Temozolomide"})
    S.panel(axs[0], "a", x=-0.42)
    fams = list(FAMILY)
    K.group_strip(axs[1], m.family.values, m.ligand_efficiency.values, [FAMILY[f] for f in fams],
                  "Ligand efficiency (kcal mol$^{-1}$ per heavy atom)", order=fams)
    axs[1].set_xticklabels([SHORT[f] for f in fams], fontsize=6)
    S.panel(axs[1], "b", x=-0.22)
    family_legend(fig, y=0.0)
    S.save(fig, FIG, "Fig3")


# ------------------------------------------------------------------ Fig. 4
def slab_for_render():
    """Slab start geometry without the edge atoms that have fewer than two
    neighbours inside the (unwrapped) cell - for the picture only."""
    L = (DFT / "slab_start.xyz").read_text().splitlines()
    n = int(L[0])
    el = [l.split()[0] for l in L[2:2 + n]]
    x = np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]])
    nb = ((np.linalg.norm(x[:, None] - x[None], axis=2) < 2.4).sum(1) - 1)
    keep = nb >= 2
    out = REN / "slab_render.xyz"
    REN.mkdir(parents=True, exist_ok=True)
    lines = [f"{e} {p[0]:.4f} {p[1]:.4f} {p[2]:.4f}" for e, p, k_ in zip(el, x, keep) if k_]
    out.write_text("\n".join([str(len(lines)), "slab (render)"] + lines) + "\n")
    return out


def fig4(d):
    """MXene carrier: slab model and why tight binding was not used."""
    slab, _ = render("slab_top", R.molecule, str(slab_for_render()), tilt=0, size=(1300, 1000))
    fig = plt.figure(figsize=(S.DOUBLE, 64 * S.MM))
    gs = GridSpec(1, 2, width_ratios=[1, 1.1], wspace=0.3, left=0.01, right=0.98, top=0.9, bottom=0.18)
    ax = fig.add_subplot(gs[0])
    K.render_panel(ax, slab, sub="Ti$_3$C$_2$O$_2$, 4×4 cell, a = 3.03 Å")
    S.panel(ax, "a", x=0.02, y=1.02)
    ax = fig.add_subplot(gs[1])
    sc = d["scan"]
    if sc is not None:
        for meth, col, mk in (("GFN2-xTB", S.DOCK, "o"), ("GFN1-xTB", S.ADS, "s")):
            x = sc[(sc.method == meth)].copy()
            ok = x.E_eV_per_fu.notna() & (x.E_eV_per_fu > -600)
            e = x.E_eV_per_fu[ok]
            ax.plot(x.a_A[ok], e - e.min(), marker=mk, ms=3.5, lw=1, color=col, label=meth)
            bad = x[~ok]
            if len(bad):
                ax.scatter(bad.a_A, np.zeros(len(bad)) - 2, marker="x", color=col, s=16, lw=0.8)
        ax.axvline(3.03, color=S.INK, lw=0.8, ls=(0, (1, 2)))
        ax.text(3.03, ax.get_ylim()[1] * 0.95, " DFT/exp.\n a = 3.03 Å", fontsize=6, va="top")
        ax.set_xlabel("In-plane lattice constant a (Å)")
        ax.set_ylabel("E − E$_{min}$ (eV per formula unit)")
        ax.legend(handles=[Line2D([], [], color=S.DOCK, marker="o", ms=3.5, label="GFN2-xTB"),
                           Line2D([], [], color=S.ADS, marker="s", ms=3.5, label="GFN1-xTB"),
                           Line2D([], [], color=S.INK, marker="x", ls="", label="SCF failed / unphysical")],
                  loc="upper center", bbox_to_anchor=(0.36, 1.0), frameon=False, fontsize=6)
        K.light_grid(ax)
    S.panel(ax, "b", x=-0.18)
    S.save(fig, FIG, "Fig4")


# ------------------------------------------------------------------ Fig. 5
def fig5(d):
    """DFT adsorption of the alkylating agents on Ti3C2O2."""
    t = d["dft"]
    if t is None or "delta_Eads_kcal_mol" not in t or t.delta_Eads_kcal_mol.notna().sum() == 0:
        print("Fig5 skipped: no completed DFT results")
        return
    t = t.dropna(subset=["delta_Eads_kcal_mol"]).sort_values("delta_Eads_kcal_mol")
    n = len(t)
    fig = plt.figure(figsize=(S.DOUBLE, 70 * S.MM))
    gs = GridSpec(1, 1 + min(n, 3), width_ratios=[1.1] + [1] * min(n, 3), wspace=0.1, left=0.07, right=0.99,
                  top=0.88, bottom=0.2)
    ax = fig.add_subplot(gs[0])
    cols = [S.CHEM if mo == "chemisorption" else S.PHYS for mo in t.adsorption_mode]
    ax.barh(np.arange(n)[::-1], -t.delta_Eads_kcal_mol, color=cols, height=0.65)
    ax.set_yticks(np.arange(n)[::-1])
    ax.set_yticklabels(t.name, fontsize=6.4)
    ax.set_xlabel("−Δ$E_{ads}$, PBE-D3 (kcal mol$^{-1}$)")
    K.light_grid(ax, "x")
    S.panel(ax, "a", x=-0.45)
    for i, (r, l) in enumerate(zip(t.head(3).itertuples(), "bcd")):
        f = DFT / f"cplx_{r.name}" / "final.xyz"
        if not f.exists():
            continue
        png, _ = render(f"dft_{r.name}", R.molecule, str(f), tilt=75, size=(1200, 1000))
        a2 = fig.add_subplot(gs[i + 1])
        K.render_panel(a2, png, title=r.name, sub=f"Δ$E_{{ads}}$ = {r.delta_Eads_kcal_mol:.1f} kcal mol$^{{-1}}$")
        S.panel(a2, l, x=0.02, y=1.02)
    S.save(fig, FIG, "Fig5")


# ------------------------------------------------------------------ Fig. 6
def fig6(d):
    """QSPR of the Vina score."""
    q, oof, perm = d["q"], d["oof"], d["perm"]
    fig, axs = plt.subplots(1, 3, figsize=(S.DOUBLE, 62 * S.MM),
                            gridspec_kw=dict(wspace=0.55, width_ratios=[1, 1, 0.9]))
    K.parity(axs[0], oof.vina.values, oof.oof_pred.values, S.DOCK, stats_lines=None, names=oof.name.values,
             label_extremes=2, xlabel="Vina score")
    K.stat_box(axs[0], [f"$n$ = {q['n']}, $p$ = {q['p']}", f"$Q^2_{{CV}}$ = {q['Q2_CV']:.2f}",
                        f"RMSE = {q['RMSE']:.2f}"], loc="upper left")
    K.williams(axs[1], oof.leverage.values, oof.std_residual.values, q["AD"]["h_star"], S.DOCK,
               names=oof.name.values)
    K.scrambling(axs[2], perm.Q2_perm.values, q["Q2_CV"], S.DOCK)
    fig.canvas.draw()
    top = max(ax.get_position().y1 for ax in axs)
    for ax, l in zip(axs, "abc"):
        fig.text(ax.get_position().x0 - 0.055, top + 0.02, l, fontsize=9, fontweight="bold", va="bottom")
    S.save(fig, FIG, "Fig6")


FIGS = {"1": fig1, "2": fig2, "3": fig3, "4": fig4, "5": fig5, "6": fig6}

if __name__ == "__main__":
    S.apply()
    d = data()
    for key in (sys.argv[1:] or FIGS):
        FIGS[key](d)
        print(f"Fig{key} done")
