"""
qspr_core.py - one leak-free QSPR protocol for every endpoint of the study.

  model       : StandardScaler + RidgeCV in one sklearn Pipeline
  validation  : nested CV - outer 5-fold (reported), inner 5-fold RidgeCV;
                pooled out-of-fold predictions -> Q2_CV, RMSE, MAE (+ per fold)
  chance      : Y-permutations through the identical nested procedure,
                p = (1 + #{Q2_perm >= Q2}) / (n_perm + 1)
  domain      : leverage from the standardised descriptor matrix with an
                intercept column, h* = 3(p+1)/n; standardised OOF residuals
  external    : model refit on all training rows, applied unchanged
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

ALPHAS = np.logspace(-3, 3, 25)
SEED = 42


def model():
    return Pipeline([("scale", StandardScaler()),
                     ("ridge", RidgeCV(alphas=ALPHAS, cv=KFold(5, shuffle=True, random_state=SEED)))])


def nested_oof(X, y):
    return cross_val_predict(model(), X, y, cv=KFold(5, shuffle=True, random_state=SEED))


def leverage(X_train, X_query):
    sc = StandardScaler().fit(X_train)
    A = np.column_stack([np.ones(len(X_train)), sc.transform(X_train)])
    Q = np.column_stack([np.ones(len(X_query)), sc.transform(X_query)])
    return np.einsum("ij,jk,ik->i", Q, np.linalg.pinv(A.T @ A), Q)


def run(df, features, target, out_dir, tag, external=None, n_perm=1000):
    """df: training rows (name + features + target). external: optional frame
    with name, features and the reference value in column `target`."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    X, y = df[features].to_numpy(float), df[target].to_numpy(float)
    n, p = X.shape
    oof = nested_oof(X, y)
    q2 = r2_score(y, oof)
    folds = [float(r2_score(y[te], oof[te])) for _, te in KFold(5, shuffle=True, random_state=SEED).split(X)]
    rng = np.random.default_rng(SEED)
    perm = np.array([r2_score(yp, nested_oof(X, yp)) for yp in (rng.permutation(y) for _ in range(n_perm))])
    p_perm = (1 + np.sum(perm >= q2)) / (n_perm + 1)
    h = leverage(X, X)
    h_star = 3 * (p + 1) / n
    res = y - oof
    std_res = res / res.std(ddof=1)
    inside = (h <= h_star) & (np.abs(std_res) <= 3)
    base = df[["name"] + features + [target]].copy()
    base.assign(oof_pred=oof, residual=res, leverage=h, std_residual=std_res, inside_AD=inside).to_csv(
        out / f"{tag}_oof.csv", index=False)
    pd.DataFrame({"Q2_perm": perm}).to_csv(out / f"{tag}_y_scrambling.csv", index=False)
    final = model().fit(X, y)
    summary = {
        "target": target, "features": features, "n": n, "p": p,
        "Q2_CV": float(q2), "RMSE": float(np.sqrt(mean_squared_error(y, oof))),
        "MAE": float(mean_absolute_error(y, oof)), "Q2_folds": folds,
        "alpha_final": float(final.named_steps["ridge"].alpha_),
        "coef_std": dict(zip(features, np.round(final.named_steps["ridge"].coef_, 4).tolist())),
        "Y_scrambling": {"n": n_perm, "mean_Q2": float(perm.mean()), "p": float(p_perm)},
        "AD": {"h_star": float(h_star), "n_inside": int(inside.sum()),
               "outside": df.loc[~inside, "name"].tolist()},
    }
    if external is not None and len(external):
        Xe = external[features].to_numpy(float)
        pred = final.predict(Xe)
        he = leverage(X, Xe)
        e = pd.DataFrame({"name": external["name"], "reference": external[target], "predicted": pred,
                          "error": pred - external[target].to_numpy(float), "leverage": he,
                          "inside_AD": he <= h_star})
        e.to_csv(out / f"{tag}_external.csv", index=False)
        summary["external"] = {"n": len(e), "MAE": float(e.error.abs().mean()),
                               "RMSE": float(np.sqrt((e.error ** 2).mean())),
                               "n_inside_AD": int(e.inside_AD.sum())}
    (out / f"{tag}_summary.json").write_text(json.dumps(summary, indent=2))
    return summary
