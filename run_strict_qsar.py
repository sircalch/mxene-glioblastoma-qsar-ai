import numpy as np, pandas as pd
from pathlib import Path


def _project_root(marker="MANIFEST_SHA256.txt"):
    from pathlib import Path as _P
    here = _P(__file__).resolve()
    for anc in [here.parent, *here.parents]:
        if (anc / marker).exists() or ((anc / "data").is_dir() and (anc / "README.md").exists()):
            return anc
    return here.parent


def _find_xtb():
    import shutil
    from pathlib import Path as _P
    w = shutil.which("xtb") or shutil.which("xtb.exe")
    if w:
        return _P(w)
    for anc in [_P(__file__).resolve().parent, *_P(__file__).resolve().parents]:
        hits = list(anc.glob("**/xtb-*/bin/xtb.exe")) or list(anc.glob("**/xtb-*/bin/xtb"))
        if hits:
            return hits[0]
    return _P("xtb")


from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
import time

proc = _project_root() / "data" / "processed"
df = pd.read_csv(proc / "dataset_drug_mxene_pristine.csv")
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
    rcv = RidgeCV(alphas=alphas)
    rcv.fit(X_tr_s, y[tr])
    y_pred_nested[te] = rcv.predict(X_te_s)

q2_nested = r2_score(y, y_pred_nested)
rmse_nested = mean_squared_error(y, y_pred_nested)**0.5
mae_nested = mean_absolute_error(y, y_pred_nested)
h_star = 3 * (4 + 1) / n

print(f"GBM Strict Nested Q2_CV = {q2_nested:.4f}")
print(f"RMSE = {rmse_nested:.3f} kcal/mol, MAE = {mae_nested:.3f} kcal/mol, h* = {h_star:.4f}")

t0 = time.time()
np.random.seed(99)
scramble_q2 = []
for _ in range(1000):
    y_perm = np.random.permutation(y)
    yp_p = np.zeros(n)
    for tr, te in outer_cv.split(X):
        scaler = StandardScaler()
        X_tr_s = scaler.fit_transform(X[tr])
        X_te_s = scaler.transform(X[te])
        rcv = RidgeCV(alphas=alphas)
        rcv.fit(X_tr_s, y_perm[tr])
        yp_p[te] = rcv.predict(X_te_s)
    scramble_q2.append(r2_score(y_perm, yp_p))
t1 = time.time()

p_val = (np.array(scramble_q2) >= q2_nested).mean()
print(f"1,000 Strict Y-scramblings in {t1-t0:.2f}s: Mean Q2 = {np.mean(scramble_q2):.4f}, Empirical p = {p_val:.4f}")
