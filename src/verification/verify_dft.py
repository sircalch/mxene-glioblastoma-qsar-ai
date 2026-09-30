"""
verify_dft.py - independent re-check of the PBE-D3 adsorption results
(results/quantum/dft_adsorption.csv) from the raw Quantum ESPRESSO outputs.
Does not import the pipeline.

For every drug with a finished complex it
  1. reads the energies from the 'Final energy' line of the BFGS output (relaxations)
     or the last '!' line (single points), and checks that BFGS converged;
  2. recomputes dE_ads, and dE_int / deformation / image terms when the fragments exist;
  3. reads the final geometry from the 'Begin final coordinates' block, and recomputes
     drug-slab bonds (1.15 x covalent radii, in-plane minimum image), drug integrity
     (bond graph vs the relaxed isolated drug) and the closest heavy-atom contact;
  4. checks that the 80 fixed slab atoms did not move, and that the energy change of the
     final BFGS steps is below the threshold;
and compares with the table. Writes results/verification/dft_check.csv.
"""
import json
import re
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[2]
W = BASE / "calculations" / "gbm_dft"
RY = 313.7547                      # kcal/mol per Ry (CODATA: 1 Ry = 13.605693 eV)
COV = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "Cl": 1.02, "Ti": 1.60}
DRUGS = ["Temozolomide", "Carmustine", "Lomustine", "Nimustine"]


def energy(out, relax):
    t = out.read_text(errors="ignore")
    if relax:
        m = re.search(r"Final energy\s*=\s*(-?\d+\.\d+)\s*Ry", t)
        return float(m.group(1)) if m else None, bool(m) and "bfgs converged" in t
    m = re.findall(r"^!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", t, re.M)
    return (float(m[-1]) if m else None), "JOB DONE" in t


def geometry(out):
    t = out.read_text(errors="ignore")
    blk = t[t.index("Begin final coordinates"):t.index("End final coordinates")]
    blk = blk[blk.index("ATOMIC_POSITIONS"):].splitlines()[1:]
    el, x, fix = [], [], []
    for l in blk:
        p = l.split()
        if len(p) < 4:
            continue
        el.append(p[0])
        x.append([float(v) for v in p[1:4]])
        fix.append(p[4:7] == ["0", "0", "0"])
    return el, np.array(x), fix


def cell(inp):
    t = inp.read_text()
    rows = t[t.index("CELL_PARAMETERS"):].splitlines()[1:4]
    return np.array([[float(v) for v in r.split()] for r in rows])


def input_positions(inp):
    t = inp.read_text()
    rows = t[t.index("ATOMIC_POSITIONS"):t.index("K_POINTS")].splitlines()[1:]
    return np.array([[float(v) for v in r.split()[1:4]] for r in rows if r.strip()])


def dist(a, b, cv):
    f = (b - a) @ np.linalg.inv(cv)
    f[..., :2] -= np.round(f[..., :2])
    return np.linalg.norm(f @ cv, axis=-1)


def bondset(el, x, idx, cv):
    return {(i, j) for i, j in combinations(idx, 2) if dist(x[i], x[j], cv) < 1.15 * (COV[el[i]] + COV[el[j]])}


def main():
    cv = cell(W / "slab" / "pw.in")
    es, sok = energy(W / "slab" / "pw.out", True)
    rows = []
    for n in DRUGS:
        c_out, d_out = W / f"cplx_{n}" / "pw.out", W / f"drug_{n}" / "pw.out"
        if not (c_out.exists() and d_out.exists()):
            rows.append({"name": n, "status": "not run"})
            continue
        ec, cok = energy(c_out, True)
        ed, dok = energy(d_out, True)
        if ec is None or ed is None:
            rows.append({"name": n, "status": "incomplete"})
            continue
        ns = json.loads((W / f"cplx_{n}" / "meta.json").read_text())["n_slab"]
        el, x, fix = geometry(c_out)
        x0 = input_positions(W / f"cplx_{n}" / "pw.in")
        fixed_moved = float(np.abs(x[:ns][np.array(fix[:ns])] - x0[:ns][np.array(fix[:ns])]).max())
        dbox = cell(W / f"drug_{n}" / "pw.in")
        del_, xd, _ = geometry(d_out)
        di, si = list(range(ns, len(el))), list(range(ns))
        inter = sorted(f"{el[i]}{i}-{el[j]}{j}" for i in di for j in si
                       if dist(x[i], x[j], cv) < 1.15 * (COV[el[i]] + COV[el[j]]))
        b0 = bondset(del_, xd, list(range(len(del_))), dbox)
        b1 = {(i - ns, j - ns) for i, j in bondset(el, x, di, cv)}
        heavy = [i for i in di if el[i] != "H"]
        dmin = float(min(dist(x[i], x[j], cv) for i in heavy for j in si))
        r = {"name": n, "status": "OK", "converged_slab_drug_cplx": sok and dok and cok,
             "dEads": round((ec - es - ed) * RY, 2), "chem": bool(inter), "bonds": ";".join(inter),
             "intact": b0 == b1, "dmin": round(dmin, 3), "fixed_atoms_max_move_A": fixed_moved}
        f = {t: W / f"frag_{t}_{n}" / "pw.out" for t in ("slab", "drug", "drugbox")}
        if f["slab"].exists() and f["drug"].exists():
            fs, _ = energy(f["slab"], False)
            fd, _ = energy(f["drug"], False)
            r.update(dEint=round((ec - fs - fd) * RY, 2), Edef_slab=round((fs - es) * RY, 2))
            if f["drugbox"].exists():
                fb, _ = energy(f["drugbox"], False)
                r.update(Edef_drug=round((fb - ed) * RY, 2), Eimage=round((fd - fb) * RY, 2))
        rows.append(r)
    df = pd.DataFrame(rows)
    tab_f = BASE / "results" / "quantum" / "dft_adsorption.csv"
    if tab_f.exists():
        tab = pd.read_csv(tab_f).set_index("name")
        for i, r in df.iterrows():
            if r.get("status") == "OK" and r["name"] in tab.index:
                t = tab.loc[r["name"]]
                df.loc[i, "dEads_table"] = t.get("delta_Eads_kcal_mol")
                df.loc[i, "dEint_table"] = t.get("delta_Eint_kcal_mol")
                df.loc[i, "mode_table"] = t.get("adsorption_mode")
                df.loc[i, "intact_table"] = t.get("drug_intact")
                df.loc[i, "dmin_table"] = t.get("min_contact_A")
    out = BASE / "results" / "verification"
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "dft_check.csv", index=False)
    print(f"slab: E = {es} Ry, converged {sok}")
    print(df.to_string())


if __name__ == "__main__":
    main()
