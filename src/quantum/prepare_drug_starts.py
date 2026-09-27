"""Lowest-energy RDKit ETKDG + MMFF94 conformer (of 8) of every GBM drug, as xyz
start geometries for periodic_adsorption.py (which runs in an env without RDKit)."""
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

BASE = Path(__file__).resolve().parents[2]
OUT = BASE / "calculations" / "gbm_periodic" / "drug_starts"
OUT.mkdir(parents=True, exist_ok=True)
for r in pd.read_csv(BASE / "data" / "processed" / "compound_library_pubchem.csv").itertuples():
    m = Chem.AddHs(Chem.MolFromSmiles(r.smiles))
    cids = list(AllChem.EmbedMultipleConfs(m, numConfs=8, randomSeed=0xC0FFEE))
    res = AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=2000)
    c = m.GetConformer(cids[min(range(len(cids)), key=lambda k: res[k][1])])
    with open(OUT / f"{r.name.replace(' ', '_')}.xyz", "w") as f:
        f.write(f"{m.GetNumAtoms()}\n{r.name} charge={Chem.GetFormalCharge(m)}\n")
        for a in m.GetAtoms():
            p = c.GetAtomPosition(a.GetIdx())
            f.write(f"{a.GetSymbol()} {p.x:.6f} {p.y:.6f} {p.z:.6f}\n")
print("done", len(list(OUT.glob('*.xyz'))))
