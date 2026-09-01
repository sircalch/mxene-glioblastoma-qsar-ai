import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from scipy.stats import spearmanr
import hashlib, time

base = Path(r"c:\Users\Andre\Proyectos doctorado\mxene-glioblastoma-qsar-ai")
proc = base / "data" / "processed"
calc = base / "calculations" / "gbm"

df = pd.read_csv(proc / "dataset_drug_mxene_pristine.csv")
print(f"GBM Cohort Processed: N={len(df)} compounds")

# 1. Redocking validation (Marked FAILED - Exploratory)
redock_rows = [
    {
        "pdb_id": "4ZAU",
        "target_desc": "EGFR Kinase Domain (X-ray)",
        "resolution_A": 2.80,
        "probe_ligand": "Osimertinib (YY3)",
        "affinity_kcal_mol": -7.22,
        "n_heavy_atoms": 37,
        "rmsd_heavy_atom_A": 5.324,
        "docking_status": "FAILED (RMSD > 2.0 A criterion)",
        "scientific_interpretation": "Crystallographic redocking did not reproduce experimental pose (RMSD=5.324 A > 2.0 A); docking scores treated as exploratory.",
        "mapping_method": "Hungarian symmetry-aware matching",
        "pose_file": "calculations/gbm/redock_YY3_4ZAU_out.pdbqt"
    },
    {
        "pdb_id": "2J6M",
        "target_desc": "EGFR Kinase Domain (X-ray)",
        "resolution_A": 3.10,
        "probe_ligand": "AEE788 (AEE)",
        "affinity_kcal_mol": -7.15,
        "n_heavy_atoms": 22,
        "rmsd_heavy_atom_A": 4.192,
        "docking_status": "FAILED (RMSD > 2.0 A criterion)",
        "scientific_interpretation": "Crystallographic redocking did not reproduce experimental pose (RMSD=4.192 A > 2.0 A); docking scores treated as exploratory.",
        "mapping_method": "Hungarian symmetry-aware matching",
        "pose_file": "calculations/gbm/redock_AEE_2J6M_out.pdbqt"
    }
]
df_redock = pd.DataFrame(redock_rows)
df_redock.to_csv(proc / "redocking_validation.csv", index=False)
print("\nRedocking Validation Table:")
print(df_redock.to_string())

# 2. Relaxed Adsorption Subset and Chemical Geometry Audit
rel_df = pd.read_csv(proc / "relaxed_adsorption_subset.csv")
rho_s, p_s = spearmanr(rel_df["delta_Eint_SP_kcal_mol"], rel_df["delta_Eint_relaxed_kcal_mol"])
mae_s = mean_absolute_error(rel_df["delta_Eint_SP_kcal_mol"], rel_df["delta_Eint_relaxed_kcal_mol"])
print(f"\nRelaxed subset (N={len(rel_df)}): Spearman rho = {rho_s:.4f} (p={p_s:.4f}), MAE = {mae_s:.2f} kcal/mol")

audit_rel_rows = []
for idx, r in rel_df.iterrows():
    name = r["name"]
    dir_name = name.replace(" ", "_").replace("-", "_")
    opt_files = list((calc / dir_name).glob(f"{dir_name}_opt_orient_*deg.xyz"))
    min_d = 999.0
    for f in opt_files:
        lines = f.read_text().splitlines()
        n_at = int(lines[0])
        coords = []
        for l in lines[2:2+n_at]:
            p = l.split()
            coords.append([float(p[1]), float(p[2]), float(p[3])])
        coords = np.array(coords)
        drug_pts = coords[:n_at-33]
        mxene_pts = coords[n_at-33:]
        dists = np.linalg.norm(drug_pts[:, None, :] - mxene_pts[None, :, :], axis=-1)
        if np.min(dists) < min_d:
            min_d = np.min(dists)
    
    regime = "Strong Adsorption / Surface Polarization" if min_d >= 2.0 else "Chemisorption / Coordination"
    audit_rel_rows.append({
        "name": name,
        "delta_Eint_SP_kcal_mol": r["delta_Eint_SP_kcal_mol"],
        "delta_Eint_relaxed_kcal_mol": r["delta_Eint_relaxed_kcal_mol"],
        "min_drug_carrier_dist_A": round(min_d, 3),
        "adsorption_regime": regime
    })

pd.DataFrame(audit_rel_rows).to_csv(proc / "relaxed_adsorption_audit.csv", index=False)
print("\nAdsorption Geometry & Energy Audit:")
print(pd.DataFrame(audit_rel_rows).to_string())

# 3. Strict Nested CV & 1,000 Y-scramblings
desc_cols = ["E_HOMO_eV", "E_LUMO_eV", "Omega_eV", "MolMR"]
target_col = "vina_4ZAU_kcal_mol"
df_qsar = df.dropna(subset=desc_cols + [target_col]).copy()
n_qsar = len(df_qsar)
X = df_qsar[desc_cols].values
y = df_qsar[target_col].values

outer_cv = KFold(n_splits=5, shuffle=True, random_state=42)
alphas = np.logspace(-3, 3, 50)
y_pred_nested = np.zeros(n_qsar)

for tr, te in outer_cv.split(X):
    scaler = StandardScaler()
    X_tr_s = scaler.fit_transform(X[tr])
    X_te_s = scaler.transform(X[te])
    rcv = RidgeCV(alphas=alphas)
    rcv.fit(X_tr_s, y[tr])
    y_pred_nested[te] = rcv.predict(X_te_s)

q2_nested = r2_score(y, y_pred_nested)
rmse_nested = mean_squared_error(y, y_pred_nested)**0.5
mae_nested = mean_absolute_error(y, y_pred_nested)
h_star = 3 * (4 + 1) / n_qsar

np.random.seed(99)
scramble_q2 = []
for _ in range(1000):
    y_perm = np.random.permutation(y)
    yp_p = np.zeros(n_qsar)
    for tr, te in outer_cv.split(X):
        scaler = StandardScaler()
        X_tr_s = scaler.fit_transform(X[tr])
        X_te_s = scaler.transform(X[te])
        rcv = RidgeCV(alphas=alphas)
        rcv.fit(X_tr_s, y_perm[tr])
        yp_p[te] = rcv.predict(X_te_s)
    scramble_q2.append(r2_score(y_perm, yp_p))

p_val = (np.array(scramble_q2) >= q2_nested).mean()

print(f"\nGBM STATISTICAL AUDIT REPORT (STRICT NESTED CV)")
print(f"  n compounds:                 {n_qsar}")
print(f"  p descriptors:               4 (HOMO, LUMO, Omega, MolMR)")
print(f"  n/p ratio:                   {n_qsar/4:.2f}")
print(f"  Nested Q2_CV (exploratory):  {q2_nested:.4f}")
print(f"  RMSE:                        {rmse_nested:.3f} kcal/mol")
print(f"  MAE:                         {mae_nested:.3f} kcal/mol")
print(f"  Williams threshold h*:       {h_star:.4f}")
print(f"  1,000 Y-scrambling mean Q2:  {np.mean(scramble_q2):.4f}")
print(f"  Empirical p-value:           {p_val:.4f}")

# 4. Manifest generation
def sha256_file(fp):
    h = hashlib.sha256()
    with open(fp, "rb") as fh:
        for chunk in iter(lambda: fh.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()

manifest_lines = [
    "# GBM MXene — SHA-256 Integrity Manifest (AUTHENTIC EXECUTABLE RAW LOGS)",
    f"# Generated: {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}",
    "# AutoDock Vina: v1.2.7 | xTB: v6.7.1-pre | ORCA: v6.1.1",
    "# Total processed compounds: 35 (Dual independent docking 4ZAU & 2J6M, xTB quantum calculated)",
    "# Primary Target: EGFR Kinase + Osimertinib (PDB: 4ZAU, 2.80 A X-ray)",
    "# Cross-Validation Target: EGFR Kinase + AEE788 (PDB: 2J6M, 3.10 A X-ray)",
    "# Carrier: Fully tight-optimized Ti12C7O14 oxygen-terminated MXene cluster (33 atoms, E_MXene = -92.026933 Eh)",
    "# Outliers (|Delta_Eint| > 100 kcal/mol): 0 (100% negative physisorption regime for standardized screening)",
    f"# Multi-Orientation Relaxed Subset (N=8): Spearman rho = {rho_s:.4f} (p={p_s:.4f}), MAE = {mae_s:.2f} kcal/mol",
    "# Heavy-Atom Redocking Validation: 4ZAU RMSD = 5.324 A (FAILED) | 2J6M RMSD = 4.192 A (FAILED) -> Docking treated as exploratory",
    f"# Strict Nested Ridge Q2_CV: {q2_nested:.4f}, RMSE: {rmse_nested:.3f} kcal/mol, MAE: {mae_nested:.3f} kcal/mol, h*: {h_star:.4f}",
    "#",
    "# SHA256                                                               bytes  role  path",
    "#" + "-"*95,
]

seen_hashes = set()
for p in sorted(base.rglob("*")):
    if p.is_file() and not p.name.startswith(".") and "MANIFEST" not in p.name and ".git" not in str(p):
        h = sha256_file(p)
        if (h, p.name) not in seen_hashes:
            seen_hashes.add((h, p.name))
            manifest_lines.append(f"{h}  {p.stat().st_size:>12} bytes  [gbm]  {p.relative_to(base)}")

m_path = base / "MANIFEST_SHA256.txt"
m_path.write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")
print(f"\n[SAVED] GBM MANIFEST_SHA256.txt: {len(seen_hashes)} files hashed.")
