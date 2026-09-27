"""
control_modes.py - heavy-atom RMSD of every Vina mode of the two erlotinib (AQ4)
redocking controls to the crystal pose (PDB 1M17), and the mode scores.

writes results/docking/redocking_all_modes.csv
"""
import re
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolAlign

sys.path.insert(0, str(Path(__file__).resolve().parent))
import vina_protocol as vp  # noqa: E402

BASE = Path(__file__).resolve().parents[2]
POSES = BASE / "results" / "docking" / "real_poses"
PDB = BASE / "data" / "raw" / "1M17.pdb"
SMILES = "COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC"


def main():
    from meeko import PDBQTMolecule, RDKitMolCreate
    ref = vp.crystal_mol(PDB, "AQ4", "A", SMILES)
    rows = []
    for control, f in (("self-redock, crystal conformation", "AQ4_redock_out.pdbqt"),
                       ("production protocol, from SMILES", "AQ4_smiles_out.pdbqt")):
        txt = (POSES / f).read_text()
        scores = [float(x) for x in re.findall(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)", txt)]
        mol = RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(txt, skip_typing=True))[0]
        mol = Chem.RemoveHs(mol)
        for i, conf in enumerate(mol.GetConformers()):
            m1 = Chem.Mol(mol, confId=conf.GetId())
            one = Chem.Mol(mol)
            one.RemoveAllConformers()
            one.AddConformer(Chem.Conformer(conf), assignId=True)
            rows.append({"control": control, "mode": i + 1, "vina_kcal_mol": scores[i],
                         "rmsd_A": round(rdMolAlign.CalcRMS(one, ref), 3)})
    df = pd.DataFrame(rows)
    df.to_csv(BASE / "results" / "docking" / "redocking_all_modes.csv", index=False)
    print(df.to_string())


if __name__ == "__main__":
    main()
