"""
contacts.py - residue contacts of the top-ranked Vina pose of every drug.

For each pose: receptor residues with any heavy atom within CUTOFF (4.0 A) of
a ligand heavy atom, and polar contacts (ligand N/O to receptor N/O within
3.5 A). Receptor = the prepared, protonated chain used for docking.

usage: python contacts.py <receptor_H.pdb> <poses_dir> <names.csv> <out.csv>
   names.csv needs a 'name' column; pose file = <slug(name)>_out.pdbqt
writes a long table: name, residue, chain, min_dist_A, polar
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from vina_protocol import slug  # noqa: E402

CUTOFF = 4.0
POLAR = 3.5
THREE = {"ALA": "Ala", "ARG": "Arg", "ASN": "Asn", "ASP": "Asp", "CYS": "Cys", "GLN": "Gln",
         "GLU": "Glu", "GLY": "Gly", "HIS": "His", "HID": "His", "HIE": "His", "HIP": "His",
         "ILE": "Ile", "LEU": "Leu", "LYS": "Lys", "MET": "Met", "PHE": "Phe", "PRO": "Pro",
         "SER": "Ser", "THR": "Thr", "TRP": "Trp", "TYR": "Tyr", "VAL": "Val"}


def receptor_atoms(pdb):
    rows = []
    for l in open(pdb):
        if l.startswith(("ATOM", "HETATM")):
            el = (l[76:78].strip() or l[12:14].strip()).upper()
            if el == "H":
                continue
            rows.append((THREE.get(l[17:20], l[17:20].strip().capitalize()) + l[22:26].strip(),
                         l[21], el, float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return rows


def pose_atoms(pdbqt):
    out = []
    for l in open(pdbqt):
        if l.startswith("ENDMDL"):
            break
        if l.startswith(("ATOM", "HETATM")):
            t = l[77:79].strip()
            if t in ("H", "HD"):
                continue
            el = {"A": "C", "OA": "O", "NA": "N", "SA": "S"}.get(t, t).upper()
            out.append((el, float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return out


def contacts(rec, lig):
    R = np.array([r[3:] for r in rec])
    L = np.array([a[1:] for a in lig])
    D = np.linalg.norm(R[:, None] - L[None], axis=2)
    rel = np.array([r[2] in ("N", "O") for r in rec])
    lpol = np.array([a[0] in ("N", "O") for a in lig])
    res = {}
    for i in np.where(D.min(1) <= CUTOFF)[0]:
        key = (rec[i][0], rec[i][1])
        d = D[i].min()
        pol = bool(rel[i] and (D[i][lpol] <= POLAR).any()) if lpol.any() else False
        m = res.setdefault(key, [9.9, False])
        m[0] = min(m[0], d)
        m[1] = m[1] or pol
    return res


def main(rec_pdb, poses, names_csv, out_csv):
    rec = receptor_atoms(rec_pdb)
    rows = []
    for name in pd.read_csv(names_csv)["name"]:
        f = Path(poses) / f"{slug(name)}_out.pdbqt"
        if not f.exists():
            continue
        for (resid, chain), (d, pol) in contacts(rec, pose_atoms(f)).items():
            rows.append({"name": name, "residue": resid, "chain": chain, "min_dist_A": round(d, 2),
                         "polar": pol})
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    print(f"wrote {out_csv}: {len(rows)} contacts")


if __name__ == "__main__":
    main(*sys.argv[1:5])
