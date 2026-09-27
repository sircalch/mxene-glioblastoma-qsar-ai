"""
run_adsorption.py
=================
GFN2-xTB adsorption of the cohort (PubChem structures) on the validated carrier
(protocol: src/quantum/adsorption_engine.py), plus the isolated-drug
electronic descriptors from the same relaxed geometries.

usage: python src/quantum/run_adsorption.py [--workers 4]
writes results/quantum/adsorption_results.csv and
       results/quantum/isolated_drugs_qm_results.csv
"""
import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import adsorption_engine as ae  # noqa: E402

BASE = Path(__file__).resolve().parents[2]
LIB = BASE / "data" / "processed" / "compound_library_pubchem.csv"
WORK = BASE / "calculations" / "adsorption"
QM = BASE / "results" / "quantum"
CARRIER = json.loads((BASE / "data" / "processed" / "carrier.json").read_text())
EXCLUDE = set(CARRIER.get("exclude_drugs", []))


def job(name, smiles):
    d = ae.relaxed_drug(name, smiles, WORK)
    c = {k: CARRIER[k] for k in ("name", "E_Eh", "uhf", "xtb_args")}
    c["xyz"] = str(BASE / CARRIER["xyz"])
    return d, ae.adsorb(d, c, WORK)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    lib = pd.read_csv(LIB)
    lib = lib[~lib.name.isin(EXCLUDE)]
    drugs, ads = [], []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = [ex.submit(job, r.name, r.smiles) for r in lib.itertuples()]
        for k, f in enumerate(as_completed(futs), 1):
            d, r = f.result()
            drugs.append(d)
            ads.append({x: y for x, y in r.items() if x != "poses"})
            print(f"[{k:02d}/{len(futs)}] {d['name']:<22} {r.get('delta_Eint_kcal_mol', r['status'])} "
                  f"{r.get('adsorption_mode', '')} intact={r.get('drug_intact')}", flush=True)
    QM.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(ads).to_csv(QM / "adsorption_results.csv", index=False)
    iso = pd.DataFrame(drugs).drop(columns="xyz")
    iso["Gap_eV"] = iso.E_LUMO_eV - iso.E_HOMO_eV
    iso["Eta_eV"] = iso.Gap_eV / 2
    iso["Mu_eV"] = (iso.E_HOMO_eV + iso.E_LUMO_eV) / 2
    iso["Omega_eV"] = iso.Mu_eV ** 2 / (2 * iso.Eta_eV)
    iso.to_csv(QM / "isolated_drugs_qm_results.csv", index=False)


if __name__ == "__main__":
    main()
