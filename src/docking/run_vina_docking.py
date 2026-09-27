"""
run_vina_docking.py
===================
AutoDock Vina docking of the GBM cohort into the ATP site of the human EGFR
kinase domain, PDB 1M17 (bound to erlotinib = AQ4), chain A, following
src/docking/vina_protocol.py.

The earlier receptor, 4ZAU, holds osimertinib, a covalent (Cys797) inhibitor:
redocking it non-covalently failed even from its crystal conformation (8.6 A),
so the classic non-covalent erlotinib complex is used to validate the protocol.

Receptor: chain A, waters and erlotinib removed, missing side-chain atoms
and hydrogens (pH 7.4) added with PDBFixer and relaxed by a short restrained
Amber14 minimisation (deposited heavy atoms held in place), Meeko typing. Box: 22 A cube centred on the crystallographic
erlotinib. Controls: redocking of erlotinib from its crystal conformation
and from SMILES through the production protocol (every mode is compared with
the crystal pose by control_modes.py). Meeko could not type two C-terminal
residues (His964, Leu977, 22-27 A from erlotinib), which were left out.

Usage:  python src/docking/run_vina_docking.py [--controls-only]
Outputs
  data/raw/1M17_A_receptor.pdbqt (+ _H.pdb, _prep.log)
  results/docking/real_poses/<drug>_out.pdbqt, <drug>_vina.log
  results/docking/real_vina_docking_summary.csv, redocking_validation.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parent))
import vina_protocol as vp  # noqa: E402

BASE = Path(__file__).resolve().parents[2]
RAW = BASE / "data" / "raw"
LIB = BASE / "data" / "processed" / "compound_library_pubchem.csv"
POSES = BASE / "results" / "docking" / "real_poses"
LIGS = RAW / "ligands_pdbqt"
VINA = (BASE / "src" / "docking" / "vina.exe").resolve()

PDB = RAW / "1M17.pdb"
CHAIN, REF = "A", "AQ4"
BOX = 22.0


def main(controls_only=False):
    POSES.mkdir(parents=True, exist_ok=True)
    LIGS.mkdir(parents=True, exist_ok=True)
    lib = pd.read_csv(LIB)
    center = vp.box_center(PDB, REF, CHAIN)
    _, receptor = vp.prepare_receptor(PDB, CHAIN, RAW / f"1M17_{CHAIN}_receptor", relax_added=True,
                                      pocket_center=center)
    print(f"receptor {receptor.name}; box centre {np.round(center, 3)} ({REF})", flush=True)
    ref_smiles = lib.set_index("name").loc["Erlotinib", "smiles"]
    if controls_only or not (BASE / "results" / "docking" / "redocking_validation.csv").exists():
        controls, nha = vp.redock(VINA, PDB, receptor, center, BOX, REF, CHAIN, ref_smiles, POSES)
        pd.DataFrame([{
            "pdb_id": "1M17", "chain": CHAIN, "target_desc": "EGFR kinase domain (X-ray, 2.6 A)",
            "probe_ligand": "erlotinib (AQ4)", "control": c, "affinity_kcal_mol": s, "n_heavy_atoms": nha,
            "rmsd_heavy_atom_A": round(r, 3),
            "docking_status": "PASSED (RMSD <= 2.0 A)" if r <= 2.0 else "FAILED (RMSD > 2.0 A)",
            "mapping_method": "RDKit CalcRMS (symmetry-aware, no re-alignment)",
            "pose_file": f"results/docking/real_poses/{p}"} for c, s, r, p in controls]).to_csv(
            BASE / "results" / "docking" / "redocking_validation.csv", index=False)
        for c, s, r, _ in controls:
            print(f"redocking {REF} [{c}]: {s:.2f} kcal/mol, RMSD {r:.2f} A", flush=True)
    if controls_only:
        return
    rows = []
    for i, r in enumerate(lib.itertuples(), 1):
        if "B" in {a.GetSymbol() for a in Chem.MolFromSmiles(r.smiles).GetAtoms()}:
            # boronic acids (bortezomib): Vina has no boron atom type
            rows.append({"name": r.name, "class": r[2], "vina_1M17_kcal_mol": float("nan"),
                         "ligand_efficiency": float("nan"), "n_ring_conformers": 0,
                         "pose_file": "not docked: boron has no AutoDock Vina atom type"})
            print(f"[{i:02d}/{len(lib)}] {r.name:<16} not docked (boron)", flush=True)
            continue
        s, n_conf = vp.dock_ensemble(VINA, receptor, vp.ligand_pdbqts(r.name, r.smiles, LIGS),
                                     center, BOX, POSES, vp.slug(r.name))
        rows.append({"name": r.name, "class": r[2], "vina_1M17_kcal_mol": s,
                     "ligand_efficiency": round(-s / r.n_heavy_atoms, 4), "n_ring_conformers": n_conf,
                     "pose_file": f"results/docking/real_poses/{vp.slug(r.name)}_out.pdbqt"})
        print(f"[{i:02d}/{len(lib)}] {r.name:<16} {s:7.3f} kcal/mol  ({n_conf} conformer(s))", flush=True)
    pd.DataFrame(rows).to_csv(BASE / "results" / "docking" / "real_vina_docking_summary.csv", index=False)


if __name__ == "__main__":
    main(controls_only="--controls-only" in sys.argv)
