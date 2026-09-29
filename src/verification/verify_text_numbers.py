"""
verify_text_numbers.py - recomputes, from the pose files and the receptor, the numbers
the Results and Methods state about docking, without the pipeline's code:
families, score range and medians, rank of erlotinib, Spearman rho(heavy atoms, score),
ligand efficiency by family, residue contacts of the top poses (heavy atoms within 4.0 A;
N/O pairs within 3.5 A), the mode-score gaps of the redocking controls, the distance of
His964 / Leu977 from erlotinib, and the length of procarbazine.
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.stats import spearmanr

BASE = Path(__file__).resolve().parents[2]
POSES = BASE / "results" / "docking" / "real_poses"


def pdbqt_heavy(f):
    """first MODEL of a Vina output: heavy atoms (element, xyz)."""
    el, xyz = [], []
    for l in Path(f).read_text().splitlines():
        if l.startswith("ENDMDL"):
            break
        if l.startswith(("ATOM", "HETATM")):
            t = l[77:79].strip()
            if t in ("H", "HD"):
                continue
            e = {"A": "C", "NA": "N", "OA": "O", "SA": "S", "N": "N", "C": "C", "O": "O", "S": "S", "F": "F",
                 "Cl": "Cl", "CL": "Cl", "Br": "Br", "I": "I", "P": "P"}.get(t, t)
            el.append(e)
            xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return np.array(el), np.array(xyz)


def receptor():
    rows = []
    for l in open(BASE / "data" / "raw" / "1M17_A_receptor.pdbqt"):
        if l.startswith(("ATOM", "HETATM")) and l[77:79].strip() not in ("H", "HD"):
            rows.append((f"{l[17:20].strip().title()}{int(l[22:26])}", l[77:79].strip()[0],
                         float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return rows


def main():
    lib = pd.read_csv(BASE / "data" / "processed" / "compound_library_pubchem.csv")
    d = pd.read_csv(BASE / "results" / "docking" / "real_vina_docking_summary.csv").dropna(subset=["vina_1M17_kcal_mol"])
    fam = lib["class"].value_counts()
    print("families (all 33):", fam.to_dict())
    d = d.merge(lib[["name", "n_heavy_atoms"]], on="name")
    d["score"] = d.vina_1M17_kcal_mol
    d = d.sort_values("score")
    print(f"docked: {len(d)}; range {d.score.max():.2f} to {d.score.min():.2f}; best three: {list(d.name[:3])}; "
          f"worst three: {list(d.name[-3:])}")
    print("family medians:", d.groupby("class").score.median().round(2).to_dict())
    d["rank"] = np.arange(1, len(d) + 1)
    print("erlotinib:", d.loc[d.name == "Erlotinib", ["score", "rank"]].values.tolist())
    print(f"Spearman rho(heavy atoms, score) = {spearmanr(d.n_heavy_atoms, d.score)[0]:.3f}")
    d["LE"] = -d.score / d.n_heavy_atoms
    print("LE median by family:", d.groupby("class").LE.median().round(3).to_dict(),
          "| top 5 LE:", list(d.sort_values("LE", ascending=False).name[:5]))

    rec = receptor()
    rxyz = np.array([r[2:] for r in rec])
    cont, pol = {}, {}
    for r in d.itertuples():
        el, xyz = pdbqt_heavy(BASE / r.pose_file)
        D = np.linalg.norm(xyz[:, None] - rxyz[None], axis=2)
        res = {rec[j][0] for j in np.where(D.min(axis=0) <= 4.0)[0]}
        lig_no = np.isin(el, ["N", "O"])
        rec_no = np.array([x[1] in "NO" for x in rec])
        pres = {rec[j][0] for j in np.where((D[lig_no][:, rec_no].min(axis=0) if lig_no.any() else np.array([9.0])) <= 3.5)[0]} \
            if lig_no.any() else set()
        pres = {[rec[k] for k in np.where(rec_no)[0]][j][0] for j in np.where(D[lig_no][:, rec_no].min(axis=0) <= 3.5)[0]} \
            if lig_no.any() else set()
        for x in res:
            cont[x] = cont.get(x, 0) + 1
        for x in pres:
            pol[x] = pol.get(x, 0) + 1
    n = len(d)
    for x in ("Leu820", "Val702", "Met769", "Thr766"):
        print(f"{x}: contact {100 * cont.get(x, 0) / n:.0f}%, polar {100 * pol.get(x, 0) / n:.0f}%")
    top = sorted(cont.items(), key=lambda t: -t[1])[:12]
    print("12 most frequent:", [(k, v) for k, v in top])
    tab = pd.read_csv(BASE / "results" / "docking" / "residue_contacts.csv")
    tc = tab.groupby("residue").name.nunique()
    diff = {k: (cont.get(k, 0), int(tc.get(k, 0))) for k in set(cont) | set(tc.index) if cont.get(k, 0) != tc.get(k, 0)}
    print("contact counts recomputed vs residue_contacts.csv, differing residues:", diff)

    for f in ("AQ4_redock_out.pdbqt", "AQ4_smiles_out.pdbqt"):
        s = [float(x) for x in re.findall(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)", (POSES / f).read_text())]
        print(f"{f}: spread of the 9 modes {max(s) - min(s):.3f}; mode6-mode1 {s[5] - s[0]:.3f}; mode7-mode1 {s[6] - s[0]:.3f}")

    aq4 = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in open(BASE / "data" / "raw" / "1M17.pdb")
                    if l.startswith("HETATM") and l[17:20] == "AQ4" and l[21] == "A"])
    for resn, num in (("HIS", 964), ("LEU", 977)):
        at = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in open(BASE / "data" / "raw" / "1M17.pdb")
                       if l.startswith("ATOM") and l[17:20] == resn and l[21] == "A" and int(l[22:26]) == num])
        print(f"{resn}{num}: min distance to erlotinib {np.linalg.norm(at[:, None] - aq4[None], axis=2).min():.1f} A")

    smi = lib.set_index("name").loc["Procarbazine", "smiles"]
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    cids = AllChem.EmbedMultipleConfs(m, 8, randomSeed=42)
    e = AllChem.MMFFOptimizeMoleculeConfs(m)
    best = int(np.argmin([x[1] for x in e]))
    X = m.GetConformer(cids[best]).GetPositions()
    print(f"procarbazine longest atom-atom distance (lowest of 8 MMFF conformers, with H): "
          f"{np.linalg.norm(X[:, None] - X[None], axis=2).max():.1f} A; 4x4 cell edge {4 * 3.03:.2f} A")


if __name__ == "__main__":
    main()
