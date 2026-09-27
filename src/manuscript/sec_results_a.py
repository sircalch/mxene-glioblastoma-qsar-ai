"""sec_results_a.py - Results of the GBM / Ti3C2O2 manuscript that do not depend
on the DFT campaign: docking controls and scores, and the tight-binding tests."""
import numpy as np
import pandas as pd
from scipy import stats

import docx_kit as k

FIG = None  # set by the caller


def f1(x):
    return f"{x:.1f}".replace("-", "−")


def f2(x):
    return "0.00" if abs(x) < 0.005 else f"{x:.2f}".replace("-", "−")


def dn(n):
    return n[0].lower() + n[1:] if n[0].isupper() and n[1:2].islower() else n


def docking_stats(d):
    md, dk, ct = d["modes"], d["docked"], d["contacts"]
    s = {}
    for key, ctl in (("x", "self-redock, crystal conformation"), ("s", "production protocol, from SMILES")):
        x = md[md.control == ctl].reset_index(drop=True)
        s[f"top_{key}"], s[f"best_{key}"] = x.loc[0], x.loc[x.rmsd_A.idxmin()]
        s[f"span_{key}"] = x.vina_kcal_mol.max() - x.vina_kcal_mol.min()
        s[f"n_ok_{key}"] = int((x.rmsd_A <= 2.0).sum())
    n = ct.name.nunique()
    s["freq"] = (ct.groupby("residue").name.nunique() / n * 100).sort_values(ascending=False)
    s["polar"] = (ct[ct.polar].groupby("residue").name.nunique() / n * 100).sort_values(ascending=False)
    s["vina"] = dk.set_index("name").vina_1M17_kcal_mol.sort_values()
    s["fam_v"] = dk.groupby("family").vina_1M17_kcal_mol.median()
    s["fam_le"] = dk.groupby("family").ligand_efficiency.median()
    s["rho_size"] = stats.spearmanr(dk.n_heavy_atoms, dk.vina_1M17_kcal_mol)
    s["kw"] = stats.kruskal(*[g.vina_1M17_kcal_mol for _, g in dk.groupby("family")])
    return s


def scan_stats(d):
    sc = d["scan"]
    g2 = sc[sc.method == "GFN2-xTB"].dropna(subset=["E_eV_per_fu"]).sort_values("a_A").reset_index(drop=True)
    g1 = sc[sc.method == "GFN1-xTB"].sort_values("a_A").reset_index(drop=True)
    jumps = np.abs(np.diff(g2.E_eV_per_fu))
    i = int(jumps.argmax())
    ok1 = g1[g1.E_eV_per_fu.notna() & (g1.E_eV_per_fu > -600)]
    wild = g1[g1.E_eV_per_fu < -600]
    ft = pd.read_csv(d["flake_tests"])
    flake = ft[~ft.structure.str.contains("discarded")]
    neutral = flake[flake.charge == 0]
    return {"g2_min": float(g2.loc[g2.E_eV_per_fu.idxmin(), "a_A"]), "g2_jump": float(jumps[i]),
            "g2_jump_a": (float(g2.a_A[i]), float(g2.a_A[i + 1])),
            "g2_at_303": float(np.interp(3.03, g2.a_A, g2.E_eV_per_fu) - g2.E_eV_per_fu.min()),
            "g1_fail": int(g1.E_eV_per_fu.isna().sum()), "g1_n": len(g1),
            "g1_fail_max": float(g1[g1.E_eV_per_fu.isna()].a_A.max()),
            "g1_wild": wild.a_A.tolist(), "g1_wild_E": wild.E_eV_per_fu.tolist(),
            "g1_ok_min": float(ok1.a_A.min()), "g1_ok_ref": float(ok1.E_eV_per_fu.median()),
            "n_flake": len(flake), "n_neutral": len(neutral),
            "n_neutral_ok": int(neutral.scf_converged.sum()),
            "neutral_ok_note": neutral[neutral.scf_converged].extra.tolist(),
            "flake_opt_ok": bool(flake[flake.run == "opt"].scf_converged.any())}


def results_docking(doc, d, c):
    s = docking_stats(d)
    tx, bx, ts, bs = s["top_x"], s["best_x"], s["top_s"], s["best_s"]
    fr, po, v, fv = s["freq"], s["polar"], s["vina"], s["fam_v"]
    k.heading(doc, "Docking into the erlotinib site of EGFR", 2)
    k.para(doc,
           "The redocking controls show a specific weakness of the score rather than of the search. The "
           f"top-ranked pose of erlotinib redocked from its crystal conformation lies {f1(tx.rmsd_A)} Å from the "
           f"crystal pose, and that of erlotinib rebuilt from SMILES {f1(ts.rmsd_A)} Å. The crystal-like pose is "
           f"nevertheless found in both runs: mode {int(bx['mode'])} of the first control ({f2(bx.rmsd_A)} Å) and "
           f"mode {int(bs['mode'])} of the second ({f2(bs.rmsd_A)} Å), scored only "
           f"{f2(bx.vina_kcal_mol - tx.vina_kcal_mol)} and {f2(bs.vina_kcal_mol - ts.vina_kcal_mol)} kcal "
           f"mol^{{−1}} above the top poses (Fig. 2a, b; Table S4). All nine modes of each control fall within "
           f"{f2(max(s['span_x'], s['span_s']))} kcal mol^{{−1}}, so the Vina score does not discriminate the "
           "crystal-like pose from the alternatives in this pocket. The ranking of the cohort is therefore an "
           "exploratory measure of fit in the ATP site, not a set of predicted binding modes.", indent=True)
    k.para(doc,
           f"The top poses of the {len(d['docked'])} docked drugs occupy the ATP site. {fr.index[0]} and "
           f"{fr.index[1]} are contacted by {'every drug' if fr.iloc[1] == 100 else f'{fr.iloc[0]:.0f}% and {fr.iloc[1]:.0f}% of drugs'}, the hinge residue "
           f"Met769 by {fr['Met769']:.0f}% and the gatekeeper Thr766 by {fr['Thr766']:.0f}%, with polar contacts "
           f"to Thr766 in {po['Thr766']:.0f}% of drugs (Fig. 2c). Scores range from {f1(v.iloc[-1])} to "
           f"{f1(v.iloc[0])} kcal mol^{{−1}} (Fig. 3a). The best-scoring drugs are the kinase inhibitors "
           f"{dn(v.index[0])}, {dn(v.index[1])} and {dn(v.index[2])}, and the worst are the small alkylating "
           f"agents (family median {f1(fv['Alkylating agents'])} kcal mol^{{−1}}, against "
           f"{f1(fv['Other kinase inhibitors'])} for the other kinase inhibitors). Erlotinib itself scores "
           f"{f1(v['Erlotinib'])} kcal mol^{{−1}}, ranking {list(v.index).index('Erlotinib') + 1}th of {len(v)}. The score becomes more "
           "favourable with molecular size (Spearman "
           f"ρ = {f2(s['rho_size'].statistic)} between heavy-atom count and score). Ligand efficiency, which "
           "divides the score by the number of heavy atoms, is highest for the small alkylating agents "
           "(Fig. 3b), a known tendency of this metric for small ligands " + c("hopkins2014") + ", even "
           "though these drugs do not act on EGFR. These trends describe how the drugs fit a kinase pocket, which is relevant to off-target "
           "binding, but they carry no information about the alkylating mechanism.", indent=True)
    return s


def results_tightbinding(doc, d, c):
    t = scan_stats(d)
    k.heading(doc, "Tight-binding models of Ti_{3}C_{2}O_{2}", 2)
    note = ("restarted from the charges of the dianion" if any("dianion" in x for x in t["neutral_ok_note"])
            else "")
    k.para(doc,
           f"Neither tight-binding model described the carrier. Of {t['n_neutral']} attempts on neutral finite "
           f"flakes (Ti_{{55}}C_{{24}}O_{{20}} and smaller cuts), {t['n_neutral_ok']} converged, and only as a single point "
           f"{note}; the geometry optimisation never converged (Table S2). Charged flakes converged, but they "
           "are not the physical system. In the periodic slab (Fig. 4b, Table S5), GFN2-xTB places the energy minimum "
           f"at a = {f2(t['g2_min'])} Å, {100 * (3.03 - t['g2_min']) / 3.03:.0f}% below the lattice constant of Ti_{{3}}C_{{2}}O_{{2}} "
           f"(3.03 Å), which lies {f1(t['g2_at_303'])} eV per formula unit above that minimum, and the energy "
           f"jumps by {f1(t['g2_jump'])} eV per formula unit between {f2(t['g2_jump_a'][0])} and "
           f"{f2(t['g2_jump_a'][1])} Å. GFN1-xTB did not converge at {t['g1_fail']} of the "
           f"{t['g1_n']} lattice constants (every value up to {f2(t['g1_fail_max'])} Å), returned energies about "
           f"{min(t['g1_wild_E']) / t['g1_ok_ref']:.0f} times too negative at {' and '.join(f2(a) for a in t['g1_wild'])} Å, and in its converged range "
           f"(a ≥ {f2(t['g1_ok_min'])} Å) decreases monotonically toward the smallest converged value, so it "
           "has no defined minimum. A description of adsorption on this carrier therefore needs density "
           "functional theory.", indent=True)
    return t
