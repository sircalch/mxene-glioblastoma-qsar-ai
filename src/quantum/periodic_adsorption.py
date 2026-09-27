"""
periodic_adsorption.py
======================
Drug adsorption on a periodic, O-terminated Ti3C2O2 MXene slab with GFN2-xTB
(tblite + ASE; run in the 'tb' conda environment:
    micromamba run -r C:/Users/Andre/mm -n tb python src/quantum/periodic_adsorption.py <cmd>)

MXenes are 2D metals: a finite flake has metallic edge states and its SCF
does not converge at the GFN2-xTB level, so the carrier is modelled as an
infinite slab (Gamma point, 25 A vacuum, Fermi smearing 1500 K).

  lattice   : in-plane constant fixed at a = 3.03 A (DFT/experimental value for
              Ti3C2O2); a GFN2-xTB scan of a on a 4x4 cell (2.95-3.15 A) has no
              minimum in that range, so the tight-binding lattice is not used
  slab      : 6x6 supercell, Ti3C2O2 (252 atoms), O over the middle-layer Ti;
              atoms of the lower half (z below the middle Ti layer) fixed
  drug      : GFN2-xTB relaxed molecule (same method, same cell, periodic)
  complex   : drug plane-parallel over the slab centre, lowest atom 3.0 A above
              the top O layer; 0 and 90 deg rotations relaxed (FIRE, fmax
              0.08 eV/A, <= 400 steps); lowest-energy bound pose kept
  energies  : dE_int = E(complex) - E(slab@complex) - E(drug@complex),
              dE_ads = E(complex) - E(slab, relaxed) - E(drug, relaxed)
              (kcal/mol); chemisorption = any drug-slab bond (< 1.15 x sum of
              covalent radii); the drug's own bond graph is compared before/after

commands:  lattice | slab | drug <name> | adsorb <name> | all
Results: calculations/gbm_periodic/...  and results/quantum/periodic_adsorption.csv
"""
import json
import sys
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.constraints import FixAtoms
from ase.io import read, write
from ase.optimize import FIRE
from tblite.ase import TBLite

BASE = Path(__file__).resolve().parents[2]
WORK = BASE / "calculations" / "gbm_periodic"
EV2KCAL = 23.0605
VAC = 25.0
N = 6
ETEMP = 1500
A_FIXED = 3.03
Z = {"C": 1.15, "Ti_out": 2.33, "O": 3.30}
COV = {"H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07, "S": 1.05,
       "Cl": 1.02, "Br": 1.20, "I": 1.39, "Ti": 1.60}


def calc():
    return TBLite(method="GFN2-xTB", electronic_temperature=ETEMP, max_iterations=500, verbosity=0)


def slab(a, n=N):
    a1, a2 = np.array([a, 0, 0]), np.array([a / 2, a * np.sqrt(3) / 2, 0])
    B, C = (a1 + a2) / 3, 2 * (a1 + a2) / 3
    sym, pos = [], []
    for i in range(n):
        for j in range(n):
            o = i * a1 + j * a2
            for e, x in (("Ti", o), ("C", o + B + [0, 0, Z["C"]]), ("C", o + C - [0, 0, Z["C"]]),
                         ("Ti", o + C + [0, 0, Z["Ti_out"]]), ("Ti", o + B - [0, 0, Z["Ti_out"]]),
                         ("O", o + [0, 0, Z["O"]]), ("O", o - [0, 0, Z["O"]])):
                sym.append(e)
                pos.append(x)
    at = Atoms(sym, positions=np.array(pos) + [0, 0, VAC / 2], cell=[n * a1, n * a2, [0, 0, VAC]], pbc=True)
    return at


def fix_lower_half(at, n_slab):
    zmid = np.median(at.positions[:n_slab][np.array(at.get_chemical_symbols()[:n_slab]) == "Ti", 2])
    at.set_constraint(FixAtoms(indices=[i for i in range(n_slab) if at.positions[i, 2] < zmid - 0.1]))


def relax(at, traj, fmax=0.08, steps=400):
    at.calc = calc()
    converged = FIRE(at, logfile=str(traj) + ".log").run(fmax=fmax, steps=steps)
    return at.get_potential_energy(), bool(converged)


def sp(at):
    at = at.copy()
    at.set_constraint()
    at.calc = calc()
    return at.get_potential_energy()


def cmd_lattice():
    WORK.mkdir(parents=True, exist_ok=True)
    rows = []
    for a in np.arange(2.95, 3.151, 0.025):
        at = slab(a, 4)
        at.calc = calc()
        rows.append((float(a), float(at.get_potential_energy()) / 16))
        print(f"a = {a:.3f}  E/f.u. = {rows[-1][1]:.4f} eV", flush=True)
    E_ = np.array(rows)[:, 1]
    json.dump({"scan": rows, "minimum_in_range": bool(0 < int(E_.argmin()) < len(E_) - 1),
               "a_used": A_FIXED}, open(WORK / "lattice.json", "w"), indent=2)
    print("energy minimum inside scan range:", 0 < int(E_.argmin()) < len(E_) - 1, "; a used =", A_FIXED)


def cmd_slab():
    a0 = A_FIXED
    at = slab(a0)
    fix_lower_half(at, len(at))
    e, ok = relax(at, WORK / "slab_opt")
    write(WORK / "slab_opt.xyz", at, format="extxyz")
    json.dump({"a": a0, "n": N, "E_eV": e, "converged": ok, "natoms": len(at),
               "formula": at.get_chemical_formula()}, open(WORK / "slab.json", "w"), indent=2)
    print(f"slab {at.get_chemical_formula()} E = {e:.4f} eV converged={ok}")


def drug_xyz(name):
    """RDKit ETKDG + MMFF start geometry (written by prepare_drug_starts.py)."""
    return WORK / "drug_starts" / f"{name.replace(' ', '_')}.xyz"


def cmd_drug(name):
    wd = WORK / "drugs" / name.replace(" ", "_")
    wd.mkdir(parents=True, exist_ok=True)
    if (wd / "drug.json").exists():
        return json.load(open(wd / "drug.json"))
    s = read(WORK / "slab_opt.xyz")
    mol = read(drug_xyz(name))
    mol.set_cell(s.cell)
    mol.pbc = True
    mol.center()
    e, ok = relax(mol, wd / "drug_opt")
    write(wd / "drug_opt.xyz", mol, format="extxyz")
    rec = {"name": name, "E_eV": e, "converged": ok}
    json.dump(rec, open(wd / "drug.json", "w"), indent=2)
    return rec


def plane_frame(x):
    return np.linalg.svd(x - x.mean(0))[2]


def bonds(sym, pos, idx_a, idx_b=None, cell=None, scale=1.15):
    from ase.geometry import get_distances
    idx_b = idx_a if idx_b is None else idx_b
    _, D = get_distances(pos[idx_a], pos[idx_b], cell=cell, pbc=True)
    out = set()
    for ii, i in enumerate(idx_a):
        for jj, j in enumerate(idx_b):
            if i < j or (idx_b is not idx_a and i != j):
                if D[ii, jj] < scale * (COV[sym[i]] + COV[sym[j]]):
                    out.add((min(i, j), max(i, j)))
    return out


def cmd_adsorb(name):
    wd = WORK / "complexes" / name.replace(" ", "_")
    wd.mkdir(parents=True, exist_ok=True)
    res = wd / "result.json"
    if res.exists():
        return json.load(open(res))
    s = read(WORK / "slab_opt.xyz")
    sl = json.load(open(WORK / "slab.json"))
    d = cmd_drug(name)
    mol = read(WORK / "drugs" / name.replace(" ", "_") / "drug_opt.xyz")
    ns = len(s)
    dp = mol.positions - mol.positions.mean(0)
    dp = dp @ plane_frame(dp).T
    top = s.positions[:, 2].max()
    centre = s.cell[0] / 2 + s.cell[1] / 2
    best, poses = None, []
    for ang in (0, 90):
        th = np.radians(ang)
        R = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
        p = dp @ R.T
        p[:, 2] += top + 3.0 - p[:, 2].min()
        p[:, :2] += centre[:2]
        cx = s.copy() + Atoms(mol.get_chemical_symbols(), positions=p)
        cx.set_cell(s.cell)
        cx.pbc = True
        fix_lower_half(cx, ns)
        e, ok = relax(cx, wd / f"o{ang}")
        write(wd / f"o{ang}.xyz", cx, format="extxyz")
        sym = cx.get_chemical_symbols()
        from ase.geometry import get_distances
        heavy_d = [i for i in range(ns, len(cx)) if sym[i] != "H"]
        _, D = get_distances(cx.positions[heavy_d], cx.positions[:ns], cell=cx.cell, pbc=True)
        mc = float(D.min())
        poses.append({"angle": ang, "E_eV": e, "converged": ok, "min_contact_A": round(mc, 3)})
        if 1.25 <= mc <= 4.5 and (best is None or e < best[0]):
            best = (e, ang, cx.copy(), mc)
    rec = {"name": name, "poses": poses}
    if best is None:
        rec["status"] = "NO_BOUND_POSE"
        json.dump(rec, open(res, "w"), indent=2)
        return rec
    e_c, ang, cx, mc = best
    write(wd / "complex_opt.xyz", cx, format="extxyz")
    sym = cx.get_chemical_symbols()
    e_sf = sp(cx[:ns])
    e_df = sp(cx[ns:])
    dmol = read(WORK / "drugs" / name.replace(" ", "_") / "drug_opt.xyz")
    b0 = bonds(dmol.get_chemical_symbols(), dmol.positions, list(range(len(dmol))), cell=dmol.cell)
    didx = list(range(ns, len(cx)))
    b1 = {(i - ns, j - ns) for i, j in bonds(sym, cx.positions, didx, cell=cx.cell)}
    inter = bonds(sym, cx.positions, didx, list(range(ns)), cell=cx.cell)
    rec.update(status="OK", best_orientation_deg=ang, min_contact_A=round(mc, 3),
               delta_Eint_kcal_mol=round((e_c - e_sf - e_df) * EV2KCAL, 3),
               delta_Eads_kcal_mol=round((e_c - sl["E_eV"] - d["E_eV"]) * EV2KCAL, 3),
               drug_intact=b0 == b1, drug_slab_bonds=len(inter),
               adsorption_mode="chemisorption" if inter else "physisorption")
    json.dump(rec, open(res, "w"), indent=2)
    return rec


if __name__ == "__main__":
    c = sys.argv[1]
    if c == "lattice":
        cmd_lattice()
    elif c == "slab":
        cmd_slab()
    elif c == "drug":
        print(cmd_drug(sys.argv[2]))
    elif c == "adsorb":
        print(json.dumps(cmd_adsorb(sys.argv[2]), indent=2))
