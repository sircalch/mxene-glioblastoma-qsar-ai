import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.linear_model import Ridge, RidgeCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from scipy.stats import spearmanr

proc = Path(r"c:\Users\Andre\Proyectos doctorado\mxene-glioblastoma-qsar-ai\data\processed")
df = pd.read_csv(proc / "dataset_drug_mxene_pristine.csv")
print(f"GBM Cohort Processed: N={len(df)} compounds")
tmz = df[df["name"] == "Temozolomide"]["delta_Eint_SP_kcal_mol"].values[0]
print(f"Temozolomide: SP Delta_Eint = {tmz:.2f} kcal/mol")

# Outlier check
outliers = df[df["delta_Eint_SP_kcal_mol"].abs() > 100]
print(f"Outliers (|Delta_Eint| > 100): {len(outliers)}")

# Relaxed subset
df_rel = pd.read_csv(proc / "relaxed_adsorption_subset.csv")
print(f"\nRelaxed subset: N={len(df_rel)} compounds")
rho_s, p_s = spearmanr(df_rel["delta_Eint_SP_kcal_mol"], df_rel["delta_Eint_relaxed_kcal_mol"])
mae_s = mean_absolute_error(df_rel["delta_Eint_SP_kcal_mol"], df_rel["delta_Eint_relaxed_kcal_mol"])
print(f"Spearman rho = {rho_s:.4f} (p={p_s:.4f}), MAE = {mae_s:.2f} kcal/mol")

# Redocking validation
df_redock = pd.read_csv(proc / "redocking_validation.csv")
print("\nRedocking Validation:")
print(df_redock.to_string())

# Nested CV
desc_cols = ["E_HOMO_eV", "E_LUMO_eV", "Omega_eV", "MolMR"]
target_col = "vina_4ZAU_kcal_mol"
df_qsar = df.dropna(subset=desc_cols + [target_col]).copy()
n = len(df_qsar)
X = df_qsar[desc_cols].values
y = df_qsar[target_col].values

outer_cv = KFold(n_splits=5, shuffle=True, random_state=42)
alphas = np.logspace(-3, 3, 50)
y_pred_nested = np.zeros(n)

for tr, te in outer_cv.split(X):
    scaler = StandardScaler()
    X_tr_s = scaler.fit_transform(X[tr])
    X_te_s = scaler.transform(X[te])
    rcv = RidgeCV(alphas=alphas, cv=5)
    rcv.fit(X_tr_s, y[tr])
    y_pred_nested[te] = rcv.predict(X_te_s)

q2_nested = r2_score(y, y_pred_nested)
rmse_nested = mean_squared_error(y, y_pred_nested)**0.5
mae_nested = mean_absolute_error(y, y_pred_nested)
h_star = 3 * (4 + 1) / n

# 1,000 Y-scramblings
best_alpha = 1.0
scaler_p = StandardScaler()
X_s_full = scaler_p.fit_transform(X)

np.random.seed(99)
scramble_q2 = []
for _ in range(1000):
    y_perm = np.random.permutation(y)
    yp_p = np.zeros(n)
    for tr, te in outer_cv.split(X_s_full):
        r_mod = Ridge(alpha=best_alpha)
        r_mod.fit(X_s_full[tr], y_perm[tr])
        yp_p[te] = r_mod.predict(X_s_full[te])
    scramble_q2.append(r2_score(y_perm, yp_p))

p_val = (np.array(scramble_q2) >= q2_nested).mean()

print(f"\nNested Q2_CV (exploratory): {q2_nested:.4f}")
print(f"RMSE: {rmse_nested:.3f} kcal/mol, MAE: {mae_nested:.3f} kcal/mol")
print(f"Williams threshold h*: {h_star:.4f}")
print(f"1,000 Y-scrambling mean Q2: {np.mean(scramble_q2):.4f}")
print(f"Empirical p-value: {p_val:.4f}")
