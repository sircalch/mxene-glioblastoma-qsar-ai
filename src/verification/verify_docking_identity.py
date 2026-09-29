"""
verify_docking_identity.py - independent re-check of the GBM docking numbers and of
the identity of every structure. Does not import the pipeline.

  1. Vina score of every drug: best 'REMARK VINA RESULT' in its pose file vs
     results/docking/real_vina_docking_summary.csv (what the manuscript reads).
  2. Redocking: RMSD of every mode of the two erlotinib controls recomputed with RDKit
     against the crystal ligand AQ4 of 1M17 chain A (symmetry-aware, no re-alignment),
     vs redocking_all_modes.csv and redocking_validation.csv.
  3. Identity: InChIKey from the stored SMILES (RDKit, largest fragment) vs the table
     and vs the InChIKey that PubChem returns for the stored CID (one batch request).
writes results/verification/docking_identity_check.csv and prints a summary.
"""
import re
import time
from pathlib import Path

import pandas as pd
import requests
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.Chem.MolStandardize import rdMolStandardize

BASE = Path(__file__).resolve().parents[2]
POSES = BASE / "results" / "docking" / "real_poses"
ERLOTINIB = "COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC"
ERLOTINIB_CID = 176870


def scores(pdbqt):
    return [float(x) for x in re.findall(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)", Path(pdbqt).read_text())]


def crystal_ligand():
    lines = [l for l in open(BASE / "data" / "raw" / "1M17.pdb")
             if l.startswith("HETATM") and l[17:20] == "AQ4" and l[21] == "A"]
    ref = Chem.MolFromPDBBlock("".join(lines) + "END\n", removeHs=True)
    return AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(ERLOTINIB), ref)


def mode_rmsds(pose, ref):
    from meeko import PDBQTMolecule, RDKitMolCreate
    mol = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule.from_file(str(pose), skip_typing=True))[0])
    out = []
    for conf in mol.GetConformers():
        one = Chem.Mol(mol)
        one.RemoveAllConformers()
        one.AddConformer(Chem.Conformer(conf), assignId=True)
        out.append(rdMolAlign.CalcRMS(one, ref))
    return out


def pubchem_inchikeys(cids):
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + ",".join(map(str, cids)) +
           "/property/InChIKey/JSON")
    for _ in range(5):
        try:
            props = requests.get(url, timeout=60).json()["PropertyTable"]["Properties"]
            return {int(p["CID"]): p["InChIKey"] for p in props}
        except Exception:
            time.sleep(5)
    return {}


def ik(smi):
    return Chem.MolToInchiKey(rdMolStandardize.FragmentParent(Chem.MolFromSmiles(smi)))


def main():
    lib = pd.read_csv(BASE / "data" / "processed" / "compound_library_pubchem.csv")
    dock = pd.read_csv(BASE / "results" / "docking" / "real_vina_docking_summary.csv")
    pc = pubchem_inchikeys(list(lib.pubchem_cid.astype(int)) + [ERLOTINIB_CID])
    m = lib.merge(dock[["name", "vina_1M17_kcal_mol", "pose_file"]], on="name", how="left")
    rows = []
    for r in m.itertuples():
        pose = BASE / r.pose_file if isinstance(r.pose_file, str) else None
        rows.append({"name": r.name, "cid": r.pubchem_cid, "vina_table": r.vina_1M17_kcal_mol,
                     "vina_pose_file": scores(pose)[0] if pose is not None and pose.exists() else float("nan"),
                     "inchikey_smiles": ik(r.smiles), "inchikey_table": r.inchikey,
                     "inchikey_pubchem_cid": pc.get(int(r.pubchem_cid))})
    df = pd.DataFrame(rows)
    df["vina_ok"] = (df.vina_table - df.vina_pose_file).abs() < 1e-3
    df["identity_ok"] = (df.inchikey_smiles == df.inchikey_table) & (df.inchikey_smiles == df.inchikey_pubchem_cid)
    df["skeleton_ok"] = df.inchikey_smiles.str[:14] == df.inchikey_pubchem_cid.fillna("").str[:14]
    out = BASE / "results" / "verification"
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "docking_identity_check.csv", index=False)

    ref = crystal_ligand()
    print(f"erlotinib SMILES = PubChem CID {ERLOTINIB_CID}: {ik(ERLOTINIB) == pc.get(ERLOTINIB_CID)}; "
          f"crystal ligand heavy atoms: {ref.GetNumAtoms()}")
    allm = pd.read_csv(BASE / "results" / "docking" / "redocking_all_modes.csv")
    val = pd.read_csv(BASE / "results" / "docking" / "redocking_validation.csv")
    worst = 0.0
    for i, (control, f) in enumerate((("self-redock, crystal conformation", "AQ4_redock_out.pdbqt"),
                                      ("production protocol, from SMILES", "AQ4_smiles_out.pdbqt"))):
        rm, sc = mode_rmsds(POSES / f, ref), scores(POSES / f)
        tab = allm[allm.control == control].sort_values("mode")
        worst = max(worst, max(abs(a - b) for a, b in zip(rm, tab.rmsd_A)),
                    max(abs(a - b) for a, b in zip(sc, tab.vina_kcal_mol)))
        print(f"{f}: {len(rm)} modes; mode-1 RMSD {rm[0]:.3f} (table {val.rmsd_heavy_atom_A[i]:.3f}), "
              f"score {sc[0]:.3f} (table {val.affinity_kcal_mol[i]:.3f}); "
              f"RMSDs {', '.join(f'{x:.2f}' for x in rm)}; best mode RMSD {min(rm):.2f} A (mode {rm.index(min(rm)) + 1})")
    print(f"max |recomputed - redocking_all_modes.csv| (RMSD A or score): {worst:.4f}")
    print(f"Vina scores matching the pose files: {df.vina_ok.sum()}/{len(df)}")
    print(f"identity (SMILES = table = PubChem CID InChIKey): {df.identity_ok.sum()}/{len(df)}; "
          f"same skeleton: {df.skeleton_ok.sum()}/{len(df)}; PubChem unreachable: {df.inchikey_pubchem_cid.isna().sum()}")
    bad = df[~(df.vina_ok & df.identity_ok)]
    if len(bad):
        print(bad[["name", "vina_table", "vina_pose_file", "inchikey_smiles", "inchikey_pubchem_cid"]].to_string())


if __name__ == "__main__":
    main()
