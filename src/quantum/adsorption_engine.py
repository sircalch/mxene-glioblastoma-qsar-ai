"""
adsorption_engine.py
====================
GFN2-xTB adsorption of drugs on a finite 2D-material flake (one protocol for
every carrier of this study).

For each drug (PubChem SMILES) and carrier:
  1. drug: RDKit ETKDG + MMFF (8 conformers, lowest MMFF kept), GFN2-xTB opt;
     single point on the relaxed drug -> E_drug, HOMO/LUMO.
  2. carrier: pre-built, validated and relaxed structure (xyz + energy + spin).
  3. placement: carrier rotated into the xy plane (its best-fit plane normal
     along z), drug rotated so its own best-fit plane is parallel to it,
     centred over the carrier centroid with its lowest atom GAP A above the
     highest carrier atom; four rotations about z (0/90/180/270 deg).
  4. each pose: GFN2-xTB opt (same charge; unpaired electrons of the carrier
     kept). Converged poses whose closest drug-carrier heavy-atom contact is
     1.25-4.0 A are bound; the lowest-energy bound pose is kept.
  5. on that pose: SP of the complex and of both fragments frozen at the
     complex geometry:
        dE_int = E_complex - E_carrier@complex - E_drug@complex
        dE_ads = E_complex - E_carrier(relaxed) - E_drug(relaxed)
     chemisorption = at least one drug-carrier bond (distance < 1.15 x sum of
     covalent radii); the drug's own bond graph is compared before/after.

Results are cached per (carrier, drug) in <work>/<carrier>/<drug>/result.json,
so the run is resumable.
"""
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np

HARTREE = 627.509474
GAP = 3.2
ANGLES = (0, 90, 180, 270)


def xtb_setup():
    xtb = os.environ.get("XTB_EXE") or shutil.which("xtb") or "C:/Users/Andre/mm/xtb/Library/bin/xtb.exe"
    env = dict(os.environ)
    share = Path(xtb).parent.parent / "share" / "xtb"
    if share.is_dir():
        env["XTBPATH"] = str(share)
    env.setdefault("OMP_NUM_THREADS", "3")
    return xtb, env


def read_xyz(p):
    L = Path(p).read_text().splitlines()
    n = int(L[0].split()[0])
    return ([x.split()[0] for x in L[2:2 + n]],
            np.array([[float(v) for v in x.split()[1:4]] for x in L[2:2 + n]]))


def write_xyz(p, el, xyz, comment=""):
    with open(p, "w") as f:
        f.write(f"{len(el)}\n{comment}\n")
        for e, (x, y, z) in zip(el, xyz):
            f.write(f"{e:2s} {x:15.8f} {y:15.8f} {z:15.8f}\n")


def run_xtb(args, cwd, timeout=3600):
    xtb, env = xtb_setup()
    t0 = time.monotonic()
    try:
        p = subprocess.run([xtb, *args], cwd=str(cwd), env=env, capture_output=True, text=True,
                           errors="replace", timeout=timeout)
    except subprocess.TimeoutExpired:   # treated as a non-converged run
        return "", time.monotonic() - t0
    return p.stdout + "\n" + p.stderr, time.monotonic() - t0


def parse(out):
    d = {"E": None, "HOMO": None, "LUMO": None, "alpha": None,
         "converged": "GEOMETRY OPTIMIZATION CONVERGED" in out and "FAILED TO CONVERGE" not in out}
    for l in out.splitlines():
        if "TOTAL ENERGY" in l:
            d["E"] = float(l.split()[3])
        elif "(HOMO)" in l:
            d["HOMO"] = float(l.split()[-2])
        elif "(LUMO)" in l:
            d["LUMO"] = float(l.split()[-2])
        elif "Mol. " in l and "(0) /au" in l:        # static polarizability alpha(0), bohr^3
            d["alpha"] = float(l.split()[-1])
    return d


def plane_frame(xyz):
    """Rotation matrix whose rows are the principal axes (last = plane normal)."""
    c = xyz - xyz.mean(0)
    return np.linalg.svd(c)[2]


def drug_from_smiles(smiles, out_xyz, seed=0xC0FFEE):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    cids = list(AllChem.EmbedMultipleConfs(m, numConfs=8, randomSeed=seed))
    if not cids:
        AllChem.EmbedMolecule(m, randomSeed=1, useRandomCoords=True)
        cids = [0]
    res = AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=2000)
    best = min(range(len(cids)), key=lambda k: res[k][1])
    xyz = m.GetConformer(cids[best]).GetPositions()
    el = [a.GetSymbol() for a in m.GetAtoms()]
    write_xyz(out_xyz, el, xyz, "RDKit ETKDG+MMFF")
    return Chem.GetFormalCharge(m)


def relaxed_drug(name, smiles, work):
    """GFN2-xTB relaxed isolated drug (cached, shared by all carriers)."""
    wd = Path(work) / "_drugs" / name.replace(" ", "_")
    wd.mkdir(parents=True, exist_ok=True)
    res = wd / "drug.json"
    if res.exists():
        return json.loads(res.read_text())
    q = drug_from_smiles(smiles, wd / "drug_raw.xyz")
    run_xtb(["drug_raw.xyz", "--opt", "--gfn", "2", "--chrg", str(q), "--uhf", "0",
             "--iterations", "500", "--namespace", "dopt"], wd)
    out, _ = run_xtb(["dopt.xtbopt.xyz", "--sp", "--gfn", "2", "--chrg", str(q), "--uhf", "0",
                      "--namespace", "dsp"], wd, 600)
    d = parse(out)
    shutil.copy(wd / "dopt.xtbopt.xyz", wd / "drug_opt.xyz")
    rec = {"name": name, "charge": q, "E_drug_Eh": d["E"], "E_HOMO_eV": d["HOMO"],
           "E_LUMO_eV": d["LUMO"], "alpha_au": d["alpha"], "xyz": str(wd / "drug_opt.xyz")}
    res.write_text(json.dumps(rec, indent=2))
    return rec


def place(d_el, d_xyz, c_el, c_xyz, angle):
    c = c_xyz - c_xyz.mean(0)
    Rc = plane_frame(c)
    c = c @ Rc.T                                  # carrier in xy plane
    d = d_xyz - d_xyz.mean(0)
    d = d @ plane_frame(d).T                      # drug plane parallel to xy
    th = np.radians(angle)
    Rz = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
    d = d @ Rz.T
    d[:, 2] += c[:, 2].max() + GAP - d[:, 2].min()
    return list(d_el) + list(c_el), np.vstack([d, c]), len(d_el)


COV = {"H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07, "S": 1.05,
       "Cl": 1.02, "Br": 1.20, "I": 1.39, "Ti": 1.60, "Pt": 1.36}


def bond_set(el, xyz, idx, scale=1.15):
    idx = list(idx)
    out = set()
    for a in range(len(idx)):
        for b in range(a + 1, len(idx)):
            i, j = idx[a], idx[b]
            if np.linalg.norm(xyz[i] - xyz[j]) < scale * (COV[el[i]] + COV[el[j]]):
                out.add((i, j))
    return out


def integrity(d_el, d_xyz, fel, fxyz, nd):
    """Internal bonds of the drug and of the carrier that changed on adsorption
    (tautomerisation, bond breaking, carrier reconstruction)."""
    drug_before = bond_set(d_el, d_xyz, range(nd))
    drug_after = bond_set(fel, fxyz, range(nd))
    inter = {(i, j) for i, j in bond_set(fel, fxyz, range(len(fel))) if (i < nd) != (j < nd)}
    fmt = lambda s, el: sorted(f"{el[i]}{i}-{el[j]}{j}" for i, j in s)
    return {"drug_bonds_broken": fmt(drug_before - drug_after, d_el),
            "drug_bonds_formed": fmt(drug_after - drug_before, fel),
            "drug_intact": drug_before == drug_after,
            "drug_carrier_bonds": fmt(inter, fel)}


def min_contact(el, xyz, nd):
    dh = [i for i in range(nd) if el[i] != "H"]
    ch = [i for i in range(nd, len(el)) if el[i] != "H"]
    return float(np.linalg.norm(xyz[dh][:, None] - xyz[ch][None], axis=2).min())


def adsorb(drug, carrier, work):
    """drug: dict from relaxed_drug(); carrier: {'name','xyz','E_Eh','uhf'[, 'xtb_args']}
    xtb_args (e.g. ['--etemp', '1500'] for metallic carriers) are applied to
    every calculation that contains the carrier."""
    wd = Path(work) / carrier["name"] / drug["name"].replace(" ", "_")
    wd.mkdir(parents=True, exist_ok=True)
    res = wd / "result.json"
    if res.exists():
        r = json.loads(res.read_text())
        if r.get("status") == "OK":
            return r
    t0 = time.monotonic()
    q, uhf = drug["charge"], carrier["uhf"]
    extra = list(carrier.get("xtb_args", []))
    seed = carrier.get("seed_charge")      # e.g. -2 for metallic MXene flakes

    def seeded(xyz, ns, chrg):
        """Converge the charged system first; its <ns>.xtbrestart seeds the next xtb call."""
        if seed:
            run_xtb([xyz, "--sp", "--gfn", "2", "--chrg", str(chrg + seed), "--uhf", str(uhf),
                     "--iterations", "1000", "--namespace", ns, *extra], wd, 3600)
    d_el, d_xyz = read_xyz(drug["xyz"])
    c_el, c_xyz = read_xyz(carrier["xyz"])
    rec = {"name": drug["name"], "carrier": carrier["name"]}
    best, poses = None, []
    for ang in ANGLES:
        el, xyz, nd = place(d_el, d_xyz, c_el, c_xyz, ang)
        write_xyz(wd / f"in_{ang}.xyz", el, xyz, f"{drug['name']} on {carrier['name']} {ang} deg")
        seeded(f"in_{ang}.xyz", f"o{ang}", q)
        out, _ = run_xtb([f"in_{ang}.xyz", "--opt", "--gfn", "2", "--chrg", str(q), "--uhf", str(uhf),
                          "--iterations", "500", "--cycles", "800", "--namespace", f"o{ang}", *extra], wd)
        d = parse(out)
        f = wd / f"o{ang}.xtbopt.xyz"
        if not (d["converged"] and f.exists()):
            poses.append({"angle": ang, "status": "not converged"})
            continue
        fel, fxyz = read_xyz(f)
        mc = min_contact(fel, fxyz, nd)
        poses.append({"angle": ang, "E_Eh": d["E"], "min_contact_A": round(mc, 3)})
        if 1.25 <= mc <= 4.0 and (best is None or d["E"] < best["E"]):
            best = {"angle": ang, "E": d["E"], "contact": mc, "path": f}
    rec["poses"] = poses
    if best is None:
        rec.update(status="NO_BOUND_POSE", seconds=round(time.monotonic() - t0, 1))
        res.write_text(json.dumps(rec, indent=2))
        return rec
    fel, fxyz = read_xyz(best["path"])
    nd = len(d_el)
    write_xyz(wd / "complex_opt.xyz", fel, fxyz, f"{drug['name']}/{carrier['name']} GFN2 opt")
    write_xyz(wd / "frag_carrier.xyz", fel[nd:], fxyz[nd:], "carrier @complex")
    write_xyz(wd / "frag_drug.xyz", fel[:nd], fxyz[:nd], "drug @complex")
    seeded("complex_opt.xyz", "csp", q)
    seeded("frag_carrier.xyz", "cf", 0)
    e_c = parse(run_xtb(["complex_opt.xyz", "--sp", "--gfn", "2", "--chrg", str(q), "--uhf", str(uhf),
                         "--namespace", "csp", *extra], wd, 900)[0])["E"]
    e_cf = parse(run_xtb(["frag_carrier.xyz", "--sp", "--gfn", "2", "--chrg", "0", "--uhf", str(uhf),
                          "--namespace", "cf", *extra], wd, 900)[0])["E"]
    e_df = parse(run_xtb(["frag_drug.xyz", "--sp", "--gfn", "2", "--chrg", str(q), "--uhf", "0",
                          "--namespace", "df"], wd, 900)[0])["E"]
    rec.update(integrity(d_el, d_xyz, fel, fxyz, nd))
    rec.update(status="OK", best_orientation_deg=best["angle"], min_contact_A=round(best["contact"], 3),
               E_complex_Eh=e_c, E_carrier_frozen_Eh=e_cf, E_drug_frozen_Eh=e_df,
               delta_Eint_kcal_mol=round((e_c - e_cf - e_df) * HARTREE, 3),
               delta_Eads_kcal_mol=round((e_c - carrier["E_Eh"] - drug["E_drug_Eh"]) * HARTREE, 3),
               adsorption_mode="chemisorption" if rec["drug_carrier_bonds"] else "physisorption",
               seconds=round(time.monotonic() - t0, 1))
    for f in wd.glob("*"):          # keep inputs/outputs that matter, drop xtb scratch
        if f.suffix in (".wbo", ".charges", ".xtbrestart") or f.name.endswith("xtbtopo.mol"):
            f.unlink()
    res.write_text(json.dumps(rec, indent=2))
    return rec
