"""
generate_gbm_master_figures.py
Master 9-Figure Q1 Scientific Visualization Engine at 300+ DPI for Article 2:
Explainable AI and Quantum Chemical Exploration of 2D Titanium Carbide MXene (Ti3C2Tx) 
Nanosheets as Targeted Nanovehicles for Glioblastoma Therapeutics.
"""

import os
import sys
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.ensemble import ExtraTreesRegressor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _pubstyle
_pubstyle.apply()
try:
    import _mol3d
except Exception:
    _mol3d = None

def get_dirs():
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    fig_dir = os.path.join(base_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)
    return base_dir, fig_dir

def make_graphical_abstract(base_dir, fig_dir):
    fig, ax = plt.subplots(figsize=(12, 6.5), dpi=300)
    ax.axis('off')
    
    # Header Banner
    ax.fill_between([0, 1], [0.88, 0.88], [1.0, 1.0], color='#0D47A1', transform=ax.transAxes)
    ax.text(0.5, 0.94, "GRAPHICAL ABSTRACT: 2D MXENE NANOVEHICLES IN GLIOBLASTOMA", 
            ha='center', va='center', fontsize=13, fontweight='bold', color='white', transform=ax.transAxes)
    
    # 3 Main Pillars
    panels = [
        ("A. 2D Titanium Carbide MXene\n\n"
         "Pristine Ti3C2O2 cluster\n"
         "Real GFN2-xTB adsorption\n"
         "Angiopep-2 / LRP-1 route:\nfuture work\n"
         "(no real data for the\nfunctionalized carrier)", 0.03, 0.12, 0.29, 0.70, "#E3F2FD", "#1565C0"),
        ("B. Physical Docking\n(AutoDock Vina v1.2.7)\n\n"
         "Human EGFR kinase\n(PDB 4ZAU, 2.80 A;\n2J6M control)\n"
         "35 CNS/GBM therapeutics\n"
         "Real Vina -4.0 to -8.9\nkcal/mol (exploratory)", 0.355, 0.12, 0.29, 0.70, "#E8F5E9", "#2E7D32"),
        ("C. Explainable AI & OECD QSAR\n\n"
         "Leak-free nested 5x5\nRidge CV\n"
         "Q2_CV = 0.65 (isolated),\n0.10 (pristine MXene)\n"
         "Top feature: MolWt / MolMR\n"
         "All compounds inside the\nWilliams domain", 0.68, 0.12, 0.29, 0.70, "#FFF3E0", "#E65100"),
    ]

    for text, x, y, w, h, bg_c, border_c in panels:
        rect = patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02",
                                      facecolor=bg_c, edgecolor=border_c, lw=2.0, transform=ax.transAxes)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=8.5, fontweight='bold', color='#1A237E', transform=ax.transAxes)

    # Flow Arrows
    arrow_props = dict(facecolor='#0D47A1', edgecolor='#0D47A1', width=3.0, headwidth=10, shrink=0.05)
    ax.annotate('', xy=(0.352, 0.47), xytext=(0.322, 0.47), xycoords='axes fraction', arrowprops=arrow_props)
    ax.annotate('', xy=(0.678, 0.47), xytext=(0.648, 0.47), xycoords='axes fraction', arrowprops=arrow_props)
    
    out_p = os.path.join(fig_dir, "fig1_graphical_abstract.png")
    plt.savefig(out_p, bbox_inches='tight')
    plt.close()
    print(f"Generated Graphical Abstract: {out_p}")

def make_fig1_workflow(base_dir, fig_dir):
    fig, ax = plt.subplots(figsize=(14, 7), dpi=300)
    ax.axis('off')
    
    boxes = [
        ("1. 2D Titanium Carbide MXene\n(pristine Ti3C2O2 cluster)", 0.05, 0.55, 0.25, 0.35, "#E3F2FD", "#1565C0"),
        ("2. Blood-Brain Barrier (BBB)\nLRP-1 transcytosis route\n(proposed; not modelled here)", 0.38, 0.55, 0.25, 0.35, "#E8F5E9", "#2E7D32"),
        ("3. Glioblastoma molecular target\nHuman EGFR kinase\n(PDB 4ZAU, 2.80 A; 2J6M control)", 0.70, 0.55, 0.25, 0.35, "#FCE4EC", "#AD1457"),
        ("4. Quantum tight-binding (GFN2-xTB)\nReal interaction energies (pristine)\n+ CDFT indices\n(Delta_E_int,SP -0.9 to -15.5 kcal/mol)", 0.04, 0.10, 0.27, 0.35, "#FFF8E1", "#F57F17"),
        ("5. Real physical docking\nAutoDock Vina v1.2.7\n(catalytic pocket; exploratory)\n(35 GBM clinical drugs)", 0.375, 0.10, 0.25, 0.35, "#EDE7F6", "#4A148C"),
        ("6. Explainable AI & OECD QSAR\nLeak-free nested Ridge CV\n(Q2_CV up to 0.65; Williams domain)", 0.685, 0.10, 0.27, 0.35, "#E0F2F1", "#00695C"),
    ]
    
    for title, x, y, w, h, bg_c, border_c in boxes:
        rect = patches.Rectangle((x, y), w, h, facecolor=bg_c, edgecolor=border_c, lw=2.0, transform=ax.transAxes, zorder=2)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, title, ha='center', va='center', fontsize=9.0, fontweight='bold', color='#1A237E', transform=ax.transAxes, zorder=3)
        
    arrow_props = dict(facecolor='#37474F', edgecolor='#37474F', width=2.5, headwidth=8, shrink=0.05)
    ax.annotate('', xy=(0.37, 0.72), xytext=(0.31, 0.72), xycoords='axes fraction', arrowprops=arrow_props)
    ax.annotate('', xy=(0.69, 0.72), xytext=(0.64, 0.72), xycoords='axes fraction', arrowprops=arrow_props)
    ax.annotate('', xy=(0.37, 0.27), xytext=(0.31, 0.27), xycoords='axes fraction', arrowprops=arrow_props)
    ax.annotate('', xy=(0.69, 0.27), xytext=(0.64, 0.27), xycoords='axes fraction', arrowprops=arrow_props)
    ax.annotate('', xy=(0.50, 0.48), xytext=(0.50, 0.54), xycoords='axes fraction', arrowprops=dict(facecolor='#1565C0', width=2.0, headwidth=7))
    
    plt.title("Figure 1: Multi-Scale Computational Workflow for 2D Ti3C2Tx MXenes in Glioblastoma Oncology", fontsize=13, fontweight='bold', pad=15)
    out_p = os.path.join(fig_dir, "fig1_gbm_workflow_methodology.png")
    plt.savefig(out_p, bbox_inches='tight')
    plt.close()
    print(f"Generated Figure 1: {out_p}")

def make_fig2_quantum(base_dir, fig_dir):
    # Was entirely hardcoded arrays (homo/lumo/eta/omega for 3 "systems"),
    # never computed from any real xTB output -- and no real complex-level
    # FMO calculation exists at all for either MXene variant (the pristine
    # dataset's E_HOMO_eV/E_LUMO_eV are the isolated-drug orbitals, reused,
    # not a distinct complex electronic structure; the Angiopep-2
    # functionalized carrier has no real data whatsoever). Now shows the
    # real per-compound distribution of GFN2-xTB frontier orbitals and CDFT
    # indices for the 35-compound isolated cohort (real *_drug_sp.out
    # single points), with no fabricated "complex" comparison implied.
    homo_lumo_csv = os.path.join(base_dir, "data", "processed", "gbm_isolated_real_homo_lumo.csv")
    df = pd.read_csv(homo_lumo_csv)
    homo = df["E_HOMO_real_eV"].values
    lumo = df["E_LUMO_real_eV"].values
    gap = lumo - homo
    mu = -(homo + lumo) / 2.0
    eta = gap / 2.0
    omega = mu ** 2 / (2.0 * eta)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=300)
    plt.subplots_adjust(top=0.86, wspace=0.28)

    ax0 = axes[0]
    ax0.hist(homo, bins=12, color='#1565C0', alpha=0.75, edgecolor='k', label=f'E_HOMO (mean={homo.mean():.2f} eV)')
    ax0.hist(lumo, bins=12, color='#D84315', alpha=0.75, edgecolor='k', label=f'E_LUMO (mean={lumo.mean():.2f} eV)')
    ax0.set_xlabel("Electronic Energy (eV)", fontsize=11)
    ax0.set_ylabel("Compound Count", fontsize=11)
    ax0.set_title(f"(a) Real GFN2-xTB Frontier Molecular Orbitals (n={len(df)})", fontsize=11.5, fontweight='bold', pad=10)
    ax0.grid(True, linestyle=':', alpha=0.6)
    ax0.legend(loc='upper left', frameon=True, fontsize=9)

    ax1 = axes[1]
    ax1.scatter(eta, omega, color='#2E7D32', s=70, edgecolor='k', alpha=0.85)
    ax1.set_xlabel(r"Chemical Hardness $\eta$ (eV)", fontsize=11)
    ax1.set_ylabel(r"Electrophilicity Index $\omega$ (eV)", fontsize=11)
    ax1.set_title("(b) Real Conceptual DFT Global Reactivity Indices", fontsize=11.5, fontweight='bold', pad=10)
    ax1.grid(True, linestyle=':', alpha=0.6)

    plt.suptitle("Figure 2: Real Quantum CDFT Electronic Reactivity of the Isolated GBM Therapeutics Cohort", fontsize=12.5, fontweight='bold', y=0.98)
    out_p = os.path.join(fig_dir, "fig2_gbm_quantum_cdft_architecture.png")
    plt.savefig(out_p, bbox_inches='tight')
    plt.close()
    print(f"Generated Figure 2: {out_p}")

def make_fig3_docking_profiles(base_dir, fig_dir):
    # Uses the validated primary target (PDB 4ZAU, 2.80 A) from the QSAR master
    # dataset -- NOT the deprecated results/docking/real_vina_docking_summary.csv,
    # which was docked against the superseded 4UV7 receptor and is inconsistent
    # with redocking_validation.csv / the QSAR features used everywhere else in
    # this project (vina_4ZAU_kcal_mol, vina_2J6M_kcal_mol).
    data_csv = os.path.join(base_dir, "data", "processed", "dataset_drug_mxene_pristine.csv")
    if not os.path.exists(data_csv):
        return
    df = pd.read_csv(data_csv).rename(columns={"vina_4ZAU_kcal_mol": "Real_Vina_Docking_Score_kcal_mol"})

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.8), dpi=300)
    plt.subplots_adjust(top=0.86, wspace=0.30, bottom=0.15)

    # Panel A: Distribution of Real Vina Scores
    ax0 = axes[0]
    sns.histplot(df['Real_Vina_Docking_Score_kcal_mol'], kde=True, color='#1565C0', bins=12, ax=ax0, edgecolor='k')
    ax0.axvline(df['Real_Vina_Docking_Score_kcal_mol'].mean(), color='r', linestyle='--', lw=2.0,
                label=f"Mean Delta_G = {df['Real_Vina_Docking_Score_kcal_mol'].mean():.2f} kcal/mol")
    ax0.set_xlabel("AutoDock Vina Real Binding Energy (kcal/mol)", fontsize=10.5, fontweight='bold')
    ax0.set_ylabel("Therapeutic Compound Count", fontsize=10.5, fontweight='bold')
    ax0.set_title("(a) Binding Affinity Distribution on EGFR (PDB: 4ZAU, 2.80 Å)", fontsize=11.5, fontweight='bold', pad=10)
    ax0.legend(loc='upper left', frameon=True)
    ax0.grid(True, linestyle=':', alpha=0.6)

    # Panel B: Top 10 Best Docked GBM Therapeutics
    ax1 = axes[1]
    df_sorted = df.sort_values(by='Real_Vina_Docking_Score_kcal_mol', ascending=True).head(10)
    colors = sns.color_palette("viridis_r", n_colors=10)
    bars = ax1.barh(df_sorted['name'], df_sorted['Real_Vina_Docking_Score_kcal_mol'], color=colors, edgecolor='k')
    ax1.set_xlim(right=0)
    ax1.tick_params(axis='y', pad=8)
    ax1.set_xlabel("Real AutoDock Vina Score (kcal/mol)", fontsize=10.5, fontweight='bold')
    ax1.set_ylabel("GBM / CNS Therapeutic", fontsize=10.5, fontweight='bold')
    ax1.set_title("(b) Top 10 High-Affinity EGFR Inhibitors (PDB 4ZAU)", fontsize=11.5, fontweight='bold', pad=10)
    ax1.invert_yaxis()
    ax1.grid(True, linestyle=':', alpha=0.6)

    for bar in bars:
        w = bar.get_width()
        ax1.text(w - 0.1, bar.get_y() + bar.get_height()/2, f"{w:.2f}",
                 va='center', ha='right', fontsize=9, fontweight='bold', color='white')

    plt.suptitle("Figure 3: Molecular Docking Statistical Profiles on Human EGFR Kinase (PDB 4ZAU, exploratory: redocking RMSD 5.32 Å)", fontsize=12.5, fontweight='bold', y=0.98)
    out_p = os.path.join(fig_dir, "fig3_gbm_docking_vina_statistical_profiles.png")
    plt.savefig(out_p, bbox_inches='tight')
    plt.close()
    print(f"Generated Figure 3: {out_p}")

def make_fig4_residues(base_dir, fig_dir):
    freq_csv = os.path.join(base_dir, "results", "docking", "residue_frequency_ranking.csv")
    if not os.path.exists(freq_csv):
        return
    df = pd.read_csv(freq_csv).head(12)
    
    fig, ax = plt.subplots(figsize=(10, 5.5), dpi=300)
    colors = sns.color_palette("rocket", n_colors=len(df))
    bars = ax.bar(df['Residue'], df['Contact_Frequency'], color=colors, edgecolor='k', lw=1.2)
    
    ax.set_xlabel("Human EGFR Catalytic Pocket Residue (PDB 4ZAU)", fontsize=11, fontweight='bold')
    ax.set_ylabel("Atomic Contact Frequency (d <= 3.8 Å)", fontsize=11, fontweight='bold')
    ax.set_title("Figure 4: Residue-Level Interaction Fingerprints & Engagement Frequencies on EGFR Kinase", fontsize=12.5, fontweight='bold', pad=12)
    ax.grid(True, linestyle=':', alpha=0.6)
    
    for bar in bars:
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, h + 0.5, str(int(h)), 
                ha='center', va='bottom', fontsize=9.5, fontweight='bold')
                
    ax.set_ylim(0, max(df['Contact_Frequency']) + 4)
    out_p = os.path.join(fig_dir, "fig4_gbm_residue_contact_frequency.png")
    plt.savefig(out_p, bbox_inches='tight')
    plt.close()
    print(f"Generated Figure 4: {out_p}")

def make_fig5_parity(base_dir, fig_dir):
    # Panels (b)/(c) were fit on `Target_DeltaG_bind` from
    # dataset_drug_Ti3C2O2_pristine.csv / dataset_drug_Ti3C2_functionalized.csv,
    # whose Delta_E_ads_kcal_mol was FABRICATED by train_gbm_qsar_models.py from
    # an empirical formula over RDKit descriptors, never a real xTB calculation
    # (despite being printed as "100% REAL"). Real GFN2-xTB single-point
    # interaction energies for all 35 compounds on the pristine MXene already
    # exist (dataset_drug_mxene_pristine.csv, delta_Eint_SP_kcal_mol -- the same
    # data used by the leak-free QSPR in scripts/run_nested_cv_leakfree.py), so
    # panel (b) now uses that real data with a leak-free nested CV. No real
    # structural/quantum data exists at all for the Ti3C2-Angiopep-2
    # functionalized carrier (no complex geometries were ever built for it) --
    # that panel is intentionally omitted rather than left on fabricated data;
    # building it requires new structural modeling of the functionalized
    # carrier, not a rewiring fix.
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import RidgeCV
    from sklearn.model_selection import KFold, cross_val_predict
    from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

    alpha_grid = np.array([0.001, 0.01, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0])

    systems = [
        ("Isolated GBM Drugs", os.path.join(base_dir, "data", "processed", "dataset_isolated_gbm_drugs.csv"),
         ["MW", "LogP", "Polarizability_alpha", "Electrophilicity_omega"], "Real_Vina_Docking_Score_kcal_mol"),
        ("Ti3C2O2 Pristine (real xTB)", os.path.join(base_dir, "data", "processed", "dataset_drug_mxene_pristine.csv"),
         ["MolWt", "MolMR", "E_HOMO_eV", "Omega_eV"], "delta_Eint_SP_kcal_mol"),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.5), dpi=300)
    plt.subplots_adjust(top=0.80, wspace=0.28, bottom=0.15)
    colors = ["#1565C0", "#2E7D32"]

    for ax_idx, (sys_name, f_path, desc_cols, target_col) in enumerate(systems):
        if not os.path.exists(f_path):
            continue
        df = pd.read_csv(f_path).dropna(subset=desc_cols + [target_col])
        X = df[desc_cols].values
        y = df[target_col].values
        n, p = X.shape

        outer_cv = KFold(n_splits=5, shuffle=True, random_state=42)
        inner_cv = KFold(n_splits=5, shuffle=True, random_state=42)
        pipe = Pipeline([("scaler", StandardScaler()), ("ridge", RidgeCV(alphas=alpha_grid, cv=inner_cv))])
        y_pred = cross_val_predict(pipe, X, y, cv=outer_cv)
        rmse = mean_squared_error(y, y_pred) ** 0.5
        mae = mean_absolute_error(y, y_pred)
        r2 = r2_score(y, y_pred)

        ax = axes[ax_idx]
        ax.scatter(y, y_pred, color=colors[ax_idx], alpha=0.85, s=70, edgecolor='k', label=f'Out-of-Fold (n={n})')
        min_v = min(y.min(), y_pred.min()) - 0.5
        max_v = max(y.max(), y_pred.max()) + 0.5
        ax.plot([min_v, max_v], [min_v, max_v], 'r--', lw=2.0, label='Ideal 1:1 Parity')

        stats_txt = f"Leak-free nested 5x5 CV (n={n}, p={p})\nRMSE = {rmse:.2f} kcal/mol\nMAE = {mae:.2f} kcal/mol\n$Q^2_{{CV}}$ = {r2:.3f}"
        ax.text(0.05, 0.95, stats_txt, transform=ax.transAxes, fontsize=8.5, va='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.85, edgecolor='#B0BEC5'))

        ax.set_title(f"({chr(97+ax_idx)}) {sys_name}", fontsize=11.5, fontweight='bold', pad=10)
        ax.set_xlabel("Real Observed (kcal/mol)", fontsize=10.5)
        if ax_idx == 0:
            ax.set_ylabel("Out-of-Fold Predicted (kcal/mol)", fontsize=10.5)
        ax.grid(True, linestyle=':', alpha=0.6)
        ax.legend(loc='lower right', fontsize=8.5, frameon=True)

    plt.suptitle("Figure 5: Leak-Free Nested CV Parity for Nano-QSAR Models (real data only)", fontsize=13, fontweight='bold', y=0.98)
    out_p = os.path.join(fig_dir, "fig5_gbm_parity_models_evaluation.png")
    plt.savefig(out_p, bbox_inches='tight')
    plt.close()
    print(f"Generated Figure 5: {out_p}")

def make_fig6_shap(base_dir, fig_dir):
    # Was fit on dataset_drug_Ti3C2_functionalized.csv's FABRICATED
    # Target_DeltaG_bind (no real structural/quantum data exists at all for
    # the Ti3C2-Angiopep-2 functionalized carrier -- no complex geometries
    # were ever built for it). Refit on the real GFN2-xTB
    # delta_Eint_SP_kcal_mol for the pristine MXene (all 35 compounds,
    # dataset_drug_mxene_pristine.csv), the same real data used in Figure 5.
    f_path = os.path.join(base_dir, "data", "processed", "dataset_drug_mxene_pristine.csv")
    if not os.path.exists(f_path):
        return
    df = pd.read_csv(f_path)
    feature_cols = ["MolWt", "MolMR", "E_HOMO_eV", "E_LUMO_eV", "Gap_eV", "Eta_eV", "Mu_eV", "Omega_eV"]
    df = df.dropna(subset=feature_cols + ["delta_Eint_SP_kcal_mol"])
    X = df[feature_cols]
    y = df['delta_Eint_SP_kcal_mol']

    model = ExtraTreesRegressor(n_estimators=100, random_state=42)
    model.fit(X, y)
    
    importances = model.feature_importances_
    indices = np.argsort(importances)[::-1][:10]
    
    top_features = [feature_cols[i] for i in indices]
    top_importances = importances[indices]
    
    fig, ax = plt.subplots(figsize=(10, 5.5), dpi=300)
    colors = sns.color_palette("Blues_r", n_colors=len(top_features))
    bars = ax.barh(top_features[::-1], top_importances[::-1], color=colors, edgecolor='k')
    
    ax.set_xlabel("Mean Absolute SHAP Value / Gini Feature Importance", fontsize=11, fontweight='bold')
    ax.set_ylabel("Molecular / Quantum CDFT Descriptor", fontsize=11, fontweight='bold')
    ax.set_title("Figure 6: Exploratory Feature Importance Rankings for 2D MXene Delivery (real ΔE_int, pristine)", fontsize=11, fontweight='bold', pad=12)
    ax.grid(True, linestyle=':', alpha=0.6)
    
    for bar in bars:
        w = bar.get_width()
        ax.text(w + 0.005, bar.get_y() + bar.get_height()/2, f"{w:.3f}", 
                va='center', ha='left', fontsize=9, fontweight='bold')
                
    ax.set_xlim(0, max(top_importances) + 0.06)
    out_p = os.path.join(fig_dir, "fig6_gbm_shap_xai_importance_rankings.png")
    plt.savefig(out_p, bbox_inches='tight')
    plt.close()
    print(f"Generated Figure 6: {out_p}")

try:
    import _pymol
except Exception:
    _pymol = None


def _pm_panel(ax, png, title=None, subtitle=None):
    import matplotlib.image as mpimg
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    if png and os.path.exists(png):
        ax.imshow(mpimg.imread(png))
    else:
        ax.text(0.5, 0.5, "render unavailable", ha="center", va="center", transform=ax.transAxes)
    if title:
        ax.set_title(title, fontsize=9.5, fontweight="bold", pad=5)
    if subtitle:
        ax.text(0.5, -0.03, subtitle, ha="center", va="top", fontsize=8.0,
                color=_pubstyle.MUTED, transform=ax.transAxes)


def make_fig9_3d_spatial(base_dir, fig_dir):
    """Figure 9 - PyMOL ray-traced renders of the MXene adsorption modes."""
    df = pd.read_csv(os.path.join(base_dir, "data", "processed",
                     "dataset_drug_mxene_pristine.csv")).set_index("name")
    ads = df["delta_Eint_SP_kcal_mol"]
    v4z = df["vina_4ZAU_kcal_mol"]
    top_vina = v4z.idxmin()
    strong_ads = ads.idxmin()
    calc = os.path.join(base_dir, "calculations", "gbm")
    C = os.path.join(fig_dir, "_pm_cache"); os.makedirs(C, exist_ok=True)

    jobs = [
        (os.path.join(calc, "Ti3C2O2_pristine.xyz"), os.path.join(C, "g9_a.png"),
         "(a)  Pristine Ti$_3$C$_2$O$_2$ MXene cluster", "GFN2-xTB optimised carrier model"),
        (os.path.join(calc, top_vina, f"{top_vina}_Ti3C2O2_complex.xyz"), os.path.join(C, "g9_b.png"),
         f"(b)  {top_vina} on Ti$_3$C$_2$O$_2$",
         f"top Vina binder (4ZAU {v4z[top_vina]:.2f} kcal/mol) · $\\Delta E_{{int,SP}}$ = {ads[top_vina]:.2f} kcal/mol"),
        (os.path.join(calc, strong_ads, f"{strong_ads}_Ti3C2O2_complex.xyz"), os.path.join(C, "g9_c.png"),
         f"(c)  {strong_ads} on Ti$_3$C$_2$O$_2$",
         f"strongest GFN2-xTB interaction · $\\Delta E_{{int,SP}}$ = {ads[strong_ads]:.2f} kcal/mol"),
    ]
    if _pymol and _pymol.AVAILABLE:
        for i, (src, png, _, _) in enumerate(jobs):
            try:
                _pymol.complex_figure(src, png, size=(1400, 1150), carbon="grey55",
                                      tilt=22, mode=("cpk" if i == 0 else "ball_stick"))
            except Exception as exc:
                print(f"[fig9 PyMOL {os.path.basename(src)}] {exc}")

    fig, axes = plt.subplots(1, 3, figsize=(11.4, 4.1))
    fig.subplots_adjust(wspace=0.05, top=0.85, bottom=0.15, left=0.02, right=0.98)
    for ax, (_, png, title, sub) in zip(axes, jobs):
        _pm_panel(ax, png, title, sub)
    fig.suptitle("Figure 9. Representative drug-carrier adsorption modes on 2D Ti$_3$C$_2$O$_2$ MXene (real GFN2-xTB geometries)",
                 fontsize=10.5, fontweight="bold", y=0.99)
    out_p = os.path.join(fig_dir, "fig9_gbm_3d_spatial_binding_modes.png")
    _pubstyle.save(fig, out_p, also_pdf=False)
    print(f"Generated Figure 9 (PyMOL ray-traced): {out_p}")

def make_fig7_correlation(base_dir, fig_dir):
    """Real Pearson inter-descriptor correlation heat-map."""
    csv_p = os.path.join(base_dir, "data", "processed", "gbm_isolated_descriptors.csv")
    if not os.path.exists(csv_p):
        return
    df = pd.read_csv(csv_p)
    cols = [c for c in ["MW", "LogP", "LogS", "WS_mg_mL", "HBA", "HBD", "PSA",
                        "RBC", "NOR", "AromRings", "Polarizability_alpha",
                        "Fraction_Csp3", "E_HOMO", "E_LUMO", "Gap_eV",
                        "Hardness_eta", "Softness_S", "Electronegativity_chi",
                        "Chemical_Potential_mu", "Electrophilicity_omega"]
            if c in df.columns]
    corr = df[cols].corr()
    fig, ax = plt.subplots(figsize=(9.6, 8.0))
    sns.heatmap(corr, annot=True, fmt=".2f", cmap="vlag", center=0, vmin=-1, vmax=1,
                cbar_kws={"label": "Pearson correlation $r$", "shrink": 0.8},
                ax=ax, annot_kws={"size": 6.0}, linewidths=0.4, linecolor="white",
                square=True)
    ax.set_title(f"Pearson inter-descriptor correlation ({len(cols)} descriptors, "
                 f"{len(df)} GBM therapeutics)")
    ax.tick_params(labelsize=6.5)
    out_p = os.path.join(fig_dir, "fig7_gbm_descriptor_correlation_matrix.png")
    _pubstyle.save(fig, out_p, also_pdf=False)
    print(f"Generated Figure 7: {out_p}")


def make_fig10_deltarho(base_dir, fig_dir):
    """Figure 10 - charge-density difference for the strongest-adsorbing MXene
    complex (Larotrectinib / Ti3C2O2), real GFN2-xTB densities. Cube + build
    script in results/quantum/drho/."""
    try:
        import _drho_fig
    except Exception as exc:
        print(f"[fig10 drho] helper unavailable: {exc}")
        return
    drho_dir = os.path.join(base_dir, "results", "quantum", "drho")
    dEint = None
    try:
        df = pd.read_csv(os.path.join(base_dir, "data", "processed",
                         "dataset_drug_mxene_pristine.csv")).set_index("name")
        dEint = float(df.loc["Larotrectinib", "delta_Eint_SP_kcal_mol"])
    except Exception:
        pass
    render = os.path.join(drho_dir, "gbm_deltarho_render.png")
    render = _drho_fig.render_isosurface(drho_dir, "gbm", render, level=0.006,
                                         turn=(-55, -12, 0))
    out_p = os.path.join(fig_dir, "fig10_gbm_charge_density_difference.png")
    _drho_fig.compose(out_p, render, 10,
                      "Interfacial charge redistribution on the Ti$_3$C$_2$O$_2$ MXene carrier",
                      "Larotrectinib", "Ti$_3$C$_2$O$_2$", 0.006, dEint_kcal=dEint)
    print(f"Generated Figure 10 (charge-density difference): {out_p}")


def generate_master_suite():
    base_dir, fig_dir = get_dirs()
    make_graphical_abstract(base_dir, fig_dir)
    make_fig1_workflow(base_dir, fig_dir)
    make_fig2_quantum(base_dir, fig_dir)
    make_fig3_docking_profiles(base_dir, fig_dir)
    make_fig4_residues(base_dir, fig_dir)
    make_fig5_parity(base_dir, fig_dir)
    make_fig6_shap(base_dir, fig_dir)
    make_fig7_correlation(base_dir, fig_dir)
    make_fig9_3d_spatial(base_dir, fig_dir)
    make_fig10_deltarho(base_dir, fig_dir)
    print("Master figure suite for Article 2 (GBM) generated successfully.")

if __name__ == "__main__":
    generate_master_suite()
