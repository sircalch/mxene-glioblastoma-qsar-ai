"""
qspr_nested_cv.py - QSPR model of the AutoDock Vina score in EGFR (PDB 1M17)
(protocol: qspr_core.py - ridge, nested 5x5 CV, 1,000 Y-permutations,
leverage applicability domain).

  target      : Vina score of the 32 docked drugs (results/docking/real_vina_docking_summary.csv)
  descriptors : fixed before any fit, RDKit on the PubChem structure -
                MW, TPSA, Crippen logP, rotatable bonds

Outputs in results/qspr/: vina_{summary.json, oof.csv, y_scrambling.csv}
"""
import json
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors

sys.path.insert(0, str(Path(__file__).resolve().parent))
import qspr_core  # noqa: E402

BASE = Path(__file__).resolve().parents[2]
FEATURES = ["MW", "TPSA", "logP", "RotB"]


def desc(smiles):
    m = Chem.MolFromSmiles(smiles)
    return {"MW": Descriptors.MolWt(m), "TPSA": Descriptors.TPSA(m), "logP": Crippen.MolLogP(m),
            "RotB": rdMolDescriptors.CalcNumRotatableBonds(m)}


def main():
    lib = pd.read_csv(BASE / "data" / "processed" / "compound_library_pubchem.csv")
    dock = pd.read_csv(BASE / "results" / "docking" / "real_vina_docking_summary.csv")
    df = lib[["name", "smiles"]].merge(dock[["name", "vina_1M17_kcal_mol"]], on="name")
    df = pd.concat([df, pd.DataFrame([desc(s) for s in df.smiles])], axis=1)
    df = df.rename(columns={"vina_1M17_kcal_mol": "vina"}).dropna(subset=["vina"]).reset_index(drop=True)
    s = qspr_core.run(df, FEATURES, "vina", BASE / "results" / "qspr", "vina")
    (BASE / "results" / "qspr" / "vina_summary.json").write_text(json.dumps(s, indent=2))
    print(json.dumps(s, indent=2))


if __name__ == "__main__":
    main()
