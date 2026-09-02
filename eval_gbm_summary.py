"""
eval_gbm_summary.py
===================
Executes fully leak-free nested CV (per-inner-fold standard scaling) for GBM and updates MANIFEST_SHA256.txt.
"""

import time, hashlib
import numpy as np, pandas as pd
from pathlib import Path
from sklearn.pipeline import Pipeline
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold, GridSearchCV
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error

base = Path(__file__).resolve().parent if "__file__" in globals() else Path(r"c:\Users\Andre\Proyectos doctorado\nano-qsar-ai-papers\mxene-glioblastoma-qsar-ai")
proc = base / "data" / "processed"

df_main = pd.read_csv(proc / "dataset_drug_mxene_pristine.csv")
desc_cols = ["E_HOMO_eV", "E_LUMO_eV", "Omega_eV", "MolMR"]
target_col = "vina_4ZAU_kcal_mol"
df_qsar = df_main.dropna(subset=desc_cols + [target_col]).copy()

n_qsar = len(df_qsar)
X = df_qsar[desc_cols].values
y = df_qsar[target_col].values
h_star = 3 * (4 + 1) / n_qsar

outer_cv = KFold(n_splits=5, shuffle=True, random_state=42)
inner_cv = KFold(n_splits=5, shuffle=True, random_state=42)
alphas = np.logspace(-3, 3, 50)
param_grid = {"ridge__alpha": alphas}

y_pred_nested = np.zeros(n_qsar)
alphas_selected = []

for tr_idx, te_idx in outer_cv.split(X):
    pipe = Pipeline([("scaler", StandardScaler()), ("ridge", Ridge())])
    gscv = GridSearchCV(pipe, param_grid=param_grid, cv=inner_cv, scoring="neg_mean_squared_error")
    gscv.fit(X[tr_idx], y[tr_idx])
    y_pred_nested[te_idx] = gscv.predict(X[te_idx])
    alphas_selected.append(gscv.best_params_["ridge__alpha"])

q2_nested = r2_score(y, y_pred_nested)
rmse = mean_squared_error(y, y_pred_nested)**0.5
mae = mean_absolute_error(y, y_pred_nested)

# Save OOF predictions
oof_df = df_qsar[["name", target_col]].copy()
oof_df["y_pred_nested_cv"] = np.round(y_pred_nested, 3)
oof_df["residual"] = np.round(oof_df[target_col] - y_pred_nested, 3)
oof_df.to_csv(proc / "nested_cv_oof_predictions.csv", index=False)

# 1,000 Y-scramblings with identical inner fold standard scaling
outer_splits_list = list(outer_cv.split(X))
fold_structures = []
for tr_idx, te_idx in outer_splits_list:
    X_tr = X[tr_idx]
    inner_splits_list = list(inner_cv.split(X_tr))
    inner_mats = []
    for in_tr, in_te in inner_splits_list:
        scaler_in = StandardScaler()
        X_in_tr_s = np.hstack([np.ones((len(in_tr), 1)), scaler_in.fit_transform(X_tr[in_tr])])
        X_in_te_s = np.hstack([np.ones((len(in_te), 1)), scaler_in.transform(X_tr[in_te])])
        p_dim = X_in_tr_s.shape[1]
        XtX = X_in_tr_s.T @ X_in_tr_s
        reg_invs = []
        for a in alphas:
            reg_mat = XtX.copy()
            for i in range(1, p_dim): reg_mat[i, i] += a
            inv_m = np.linalg.pinv(reg_mat) @ X_in_tr_s.T
            reg_invs.append(inv_m)
        inner_mats.append((in_tr, in_te, X_in_te_s, reg_invs))
        
    scaler_out = StandardScaler()
    X_tr_s = np.hstack([np.ones((len(tr_idx), 1)), scaler_out.fit_transform(X_tr)])
    X_te_s = np.hstack([np.ones((len(te_idx), 1)), scaler_out.transform(X[te_idx])])
    p_dim = X_tr_s.shape[1]
    XtX_out = X_tr_s.T @ X_tr_s
    out_reg_invs = []
    for a in alphas:
        reg_mat = XtX_out.copy()
        for i in range(1, p_dim): reg_mat[i, i] += a
        inv_m = np.linalg.pinv(reg_mat) @ X_tr_s.T
        out_reg_invs.append(inv_m)
    fold_structures.append((tr_idx, te_idx, inner_mats, X_te_s, out_reg_invs))

np.random.seed(99)
scramble_q2 = []
for _ in range(1000):
    y_perm = np.random.permutation(y)
    yp_p = np.zeros(n_qsar)
    for tr_idx, te_idx, inner_mats, X_te_s, out_reg_invs in fold_structures:
        y_tr_perm = y_perm[tr_idx]
        best_alpha_idx = 0
        best_mse = 1e9
        for a_idx, a in enumerate(alphas):
            mse_tot = 0.0
            for in_tr, in_te, X_in_te_s, reg_invs in inner_mats:
                beta_in = reg_invs[a_idx] @ y_tr_perm[in_tr]
                y_in_pred = X_in_te_s @ beta_in
                mse_tot += np.mean((y_tr_perm[in_te] - y_in_pred)**2)
            if mse_tot < best_mse:
                best_mse = mse_tot
                best_alpha_idx = a_idx
        beta_out = out_reg_invs[best_alpha_idx] @ y_tr_perm
        yp_p[te_idx] = X_te_s @ beta_out
    scramble_q2.append(r2_score(y_perm, yp_p))

p_val = (np.array(scramble_q2) >= q2_nested).mean()

print("\n" + "="*70)
print(f"GBM FULLY LEAK-FREE NESTED CV REPORT (N={n_qsar}, h*={h_star:.4f})")
print(f"  Strict Nested Q2_CV:        {q2_nested:.4f}")
print(f"  RMSE:                       {rmse:.3f} kcal/mol")
print(f"  MAE:                        {mae:.3f} kcal/mol")
print(f"  Alphas per outer fold:      {[round(a, 4) for a in alphas_selected]}")
print(f"  1,000 Y-scrambling mean Q2: {np.mean(scramble_q2):.4f}")
print(f"  Empirical p-value:          {p_val:.4f}")
print("="*70)

# Manifest Generation
def sha256_file(fp):
    h = hashlib.sha256()
    with open(fp, "rb") as fh:
        for chunk in iter(lambda: fh.read(8192), b""): h.update(chunk)
    return h.hexdigest()

seen_hashes = set()
m_lines = [
    "# Glioblastoma MXene — SHA-256 Integrity Manifest (AUTHENTIC EXECUTABLE RAW LOGS)",
    f"# Generated: {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}",
    "# AutoDock Vina: v1.2.7 | xTB: v6.7.1-pre | ORCA: v6.1.1",
    f"# Total processed compounds: {len(df_main)} (34 organic + cisplatin control, Vina docking on 4ZAU & 2J6M, xTB calculated)",
    "# Primary Target: EGFR Kinase Domain (PDB: 4ZAU, 2.80 A X-ray) | Secondary: (PDB: 2J6M, 3.10 A X-ray)",
    "# Heavy-Atom Redocking Validation: 4ZAU (5.324 A, FAILED) | 2J6M (4.192 A, FAILED) -> Docking treated as exploratory",
    "# Carrier: Fully tight-optimized 33-atom metallic Ti12C7O14 cluster (E_carrier = -92.026933 Eh)",
    "# Adsorption Energetics: Subsystem deformation & frozen interaction decomposed on true xtbopt.xyz geometries",
    f"# Strict Leak-Free Nested Ridge Q2_CV: {q2_nested:.4f}, RMSE: {rmse:.3f} kcal/mol, MAE: {mae:.3f} kcal/mol, h*: {h_star:.4f}",
    "#",
    "# SHA256                                                               bytes  role  path",
    "#" + "-"*95,
]
for p in sorted(base.rglob("*")):
    if p.is_file() and not p.name.startswith(".") and "MANIFEST" not in p.name and ".git" not in str(p):
        h = sha256_file(p)
        if (h, p.name) not in seen_hashes:
            seen_hashes.add((h, p.name))
            m_lines.append(f"{h}  {p.stat().st_size:>12} bytes  [gbm]  {p.relative_to(base)}")
(base / "MANIFEST_SHA256.txt").write_text("\n".join(m_lines) + "\n", encoding="utf-8")
print(f"[SAVED] GBM MANIFEST_SHA256.txt: {len(seen_hashes)} files hashed.")
