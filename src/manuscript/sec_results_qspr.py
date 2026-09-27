"""sec_results_qspr.py - QSPR results of the GBM manuscript (docking score model)."""
import docx_kit as k


def f1(x):
    return f"{x:.1f}".replace("-", "−")


def f2(x):
    return "0.00" if abs(x) < 0.005 else f"{x:.2f}".replace("-", "−")


def results_qspr(doc, d, c, fig_path, table_no=1):
    q, oof = d["q"], d["oof"]
    coef = q["coef_std"]
    ev = oof.set_index("name").loc["Everolimus"]
    k.heading(doc, "QSPR model of the docking score", 2)
    k.para(doc,
           f"For the {q['n']} docked drugs, the ridge model on four RDKit descriptors reaches an out-of-fold "
           f"*Q*^{{2}}_{{CV}} of {f2(q['Q2_CV'])} (RMSE {f2(q['RMSE'])} kcal mol^{{−1}}, against a standard "
           f"deviation of the score of {f2(oof.vina.std())}; Table {table_no}, Fig. 6). Only "
           f"{int((d['perm'].Q2_perm >= q['Q2_CV']).sum())} of 1,000 permuted targets gave a *Q*^{{2}}_{{CV}} "
           f"as high (*p* = {q['Y_scrambling']['p']:.3f}), so the relation is not chance, but it is weak and fold-dependent: the outer-fold *Q*^{{2}} values range from "
           f"{f2(min(q['Q2_folds']))} to {f2(max(q['Q2_folds']))}. The largest standardised coefficients are "
           f"those of log*P* ({f2(coef['logP'])}) and molecular weight ({f2(coef['MW'])}), which reproduce the "
           "size dependence of the score. Everolimus, a 68-heavy-atom macrolide, lies outside the "
           f"applicability domain (leverage {f2(ev.leverage)} > *h*^{{*}} = {f2(q['AD']['h_star'])}). Extrapolating "
           f"from its size, the model predicts {f1(ev.oof_pred)} kcal mol^{{−1}}, whereas Vina gives "
           f"{f1(ev.vina)}, consistent with a molecule too large to be accommodated fully in the pocket. The model summarises the "
           "size bias of the score and should not be used to rank new compounds.", indent=True)
    k.figure(doc, fig_path, 6,
             "Validation of the QSPR model of the Vina score. **a** Out-of-fold predictions from nested 5×5 "
             "cross-validation; the shaded band is ±10% of the range. **b** Williams plot (dotted line, "
             "*h*^{*}; open symbol, outside the applicability domain). **c** Distribution of *Q*^{2}_{CV} over "
             "1,000 permuted targets, with the value of the model.")
    names = {"MW": "MW", "TPSA": "TPSA", "logP": "log*P*", "RotB": "rotatable bonds"}
    k.table(doc, (table_no, "Ridge QSPR model of the Vina score (MW, TPSA, log*P*, rotatable bonds) under nested 5×5 "
                     "cross-validation."),
            ["Quantity", "Value"],
            [["*n*, *p*", f"{q['n']}, {q['p']}"],
             ["*Q*^{2}_{CV} (out-of-fold)", f2(q["Q2_CV"])],
             ["*Q*^{2} of the outer folds", ", ".join(f2(x) for x in q["Q2_folds"])],
             ["RMSE / MAE (kcal mol^{−1})", f"{f2(q['RMSE'])} / {f2(q['MAE'])}"],
             ["Ridge penalty λ (inner CV)", f"{q['alpha_final']:g}"],
             ["Standardised coefficients", ", ".join(f"{names[a]} {f2(b)}" for a, b in coef.items())],
             ["Y-scrambling, 1,000 runs: mean *Q*^{2}; *p*",
              f"{f2(q['Y_scrambling']['mean_Q2'])}; {q['Y_scrambling']['p']:.3f}"],
             ["Applicability domain: *h*^{*}; outside", f"{f2(q['AD']['h_star'])}; {', '.join(q['AD']['outside'])}"]],
            align="ll", font=8.5)
