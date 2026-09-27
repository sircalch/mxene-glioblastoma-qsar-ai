"""
qe_mxene_dft.py
===============
Periodic DFT (Quantum ESPRESSO, PBE-D3) for drug adsorption on O-terminated
Ti3C2O2 MXene.

GFN2-xTB and GFN1-xTB were tested first and are not usable for this carrier:
finite flakes do not reach SCF convergence, and in periodic form GFN2-xTB puts
the lattice minimum near a = 2.70 A (experiment/DFT 3.03 A), with energy jumps
of tens of eV between neighbouring lattice constants. DFT is therefore used,
for a subset of drugs that fits the 4x4 cell.

  cell      : 4x4 Ti3C2O2 (112 atoms), a = 3.03 A, 20 A vacuum
  method    : PBE + DFT-D3, SSSP 1.3 efficiency pseudopotentials,
              ecutwfc/ecutrho = 50/400 Ry, Marzari-Vanderbilt smearing 0.01 Ry,
              Gamma point for slab and complex (a 2x2x1 grid needs ~23 GB, more
              than the machine offers; slab and complex share the cell, so the
              k-point error largely cancels in dE_ads); the isolated drug is
              relaxed at Gamma in a 22 A cubic box (no image interactions)
  scf       : local-TF mixing (beta 0.2, 10 vectors) - plain mixing sloshed for
              the metallic slab; a relaxation interrupted for this reason
              continues from its last geometry (cplx_<drug>/restart_from.out)
  relax     : BFGS, forces < 2e-3 Ry/bohr; only the surface that meets the drug
              relaxes (top O layer and outer Ti layer); the Ti3C2 core and the
              lower surface stay at the ideal lattice, identical in slab and
              complexes (relaxing the whole upper half took ~45 min per BFGS
              step for the metallic slab); symmetry off (nosym), since fixed atoms and the adsorbate break
              it and QE otherwise stops in checkallsym
  complex   : drug plane-parallel over the slab centre, lowest atom 2.6 A above
              the top O layer
  energies  : dE_ads = E(complex) - E(slab) - E(drug), all relaxed (kcal/mol);
              chemisorption and drug integrity as in the xTB engine

usage (Windows side writes inputs, WSL side runs pw.x):
  python src/quantum/qe_mxene_dft.py inputs            # relaxations (run_all.sh)
  python src/quantum/qe_mxene_dft.py fragments         # frozen-fragment SPs (run_fragments.sh)
  python src/quantum/qe_mxene_dft.py checks            # cutoff / k-point checks (run_checks.sh)
  python src/quantum/qe_mxene_dft.py collect           # -> results/quantum/dft_adsorption.csv
  python src/quantum/qe_mxene_dft.py collect-checks    # -> results/quantum/dft_convergence_checks.csv
  wsl -d Ubuntu-24.04 -- bash -l calculations/gbm_dft/run_all.sh
"""
import json
import re
import sys
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parents[2]
WORK = BASE / "calculations" / "gbm_dft"
STARTS = BASE / "calculations" / "gbm_periodic" / "drug_starts"
OUT = BASE / "results" / "quantum" / "dft_adsorption.csv"

A = 3.03
N = 4
VAC = 20.0
Z = {"C": 1.15, "Ti_out": 2.33, "O": 3.30}
GAP = 2.6
DRUGS = ["Temozolomide", "Carmustine", "Lomustine", "Nimustine"]
# procarbazine (13.8 A long) overlaps its own periodic image in the 12.1 A cell
# (0.6 A) and was dropped; MIN_IMAGE is the smallest acceptable drug-image distance
MIN_IMAGE = 3.0
PSEUDO_DIR = "/home/andres/pseudo/SSSP_1.3.0_PBE_efficiency"
PP = {"Ti": "ti_pbe_v1.4.uspp.F.UPF", "C": "C.pbe-n-kjpaw_psl.1.0.0.UPF", "O": "O.pbe-n-kjpaw_psl.0.1.UPF",
      "N": "N.pbe-n-radius_5.UPF", "H": "H.pbe-rrkjus_psl.1.0.0.UPF", "Cl": "cl_pbe_v1.4.uspp.F.UPF"}
MASS = {"Ti": 47.867, "C": 12.011, "O": 15.999, "N": 14.007, "H": 1.008, "Cl": 35.45}
RY2KCAL = 313.7547
DBOX = np.diag([22.0, 22.0, 22.0])
COV = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "Cl": 1.02, "Ti": 1.60}


def cell():
    a1, a2 = np.array([A, 0, 0]), np.array([A / 2, A * np.sqrt(3) / 2, 0])
    thick = 2 * Z["O"]
    return a1, a2, np.array([N * a1, N * a2, [0, 0, thick + VAC]])


def slab():
    a1, a2, cv = cell()
    B, C = (a1 + a2) / 3, 2 * (a1 + a2) / 3
    z0 = cv[2, 2] / 2
    el, xyz = [], []
    for i in range(N):
        for j in range(N):
            o = i * a1 + j * a2
            for e, x in (("Ti", o), ("C", o + B + [0, 0, Z["C"]]), ("C", o + C - [0, 0, Z["C"]]),
                         ("Ti", o + C + [0, 0, Z["Ti_out"]]), ("Ti", o + B - [0, 0, Z["Ti_out"]]),
                         ("O", o + [0, 0, Z["O"]]), ("O", o - [0, 0, Z["O"]])):
                el.append(e)
                xyz.append(x + [0, 0, z0])
    return el, np.array(xyz), cv


def read_xyz(f):
    L = Path(f).read_text().split("\n")
    n = int(L[0])
    el = [l.split()[0] for l in L[2:2 + n]]
    xyz = np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]])
    return el, xyz


def image_distance(x, cv):
    """Shortest distance between the drug and its in-plane periodic images."""
    best = np.inf
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            if i or j:
                sh = i * cv[0] + j * cv[1]
                best = min(best, np.linalg.norm(x[:, None] - (x[None] + sh), axis=2).min())
    return best


def place(del_, dxyz, sel, sxyz, cv):
    """Drug mean plane parallel to the slab, lowest atom GAP above the top O
    layer. If this leaves less than MIN_IMAGE A to the drug's periodic images,
    the drug is rotated about the surface normal (5 deg steps) to the angle
    that maximises that distance."""
    d = dxyz - dxyz.mean(0)
    _, _, vt = np.linalg.svd(d)
    d = d @ vt.T                       # mean plane -> xy
    if image_distance(d, cv) < MIN_IMAGE:
        def rot(t):
            c_, s_ = np.cos(np.radians(t)), np.sin(np.radians(t))
            return d @ np.array([[c_, -s_, 0], [s_, c_, 0], [0, 0, 1]]).T
        d = rot(max(range(0, 180, 5), key=lambda t: image_distance(rot(t), cv)))
    d[:, 2] += sxyz[:, 2].max() + GAP - d[:, 2].min()
    d[:, :2] += (cv[0, :2] + cv[1, :2]) / 2
    return sel + del_, np.vstack([sxyz, d])


def pw_input(prefix, el, xyz, cv, kpts, fix_below=None, calc="relax", ecut=(50, 400)):
    species = sorted(set(el), key=lambda e: ["Ti", "C", "O", "N", "H", "Cl"].index(e))
    lines = [
        "&CONTROL", f"  calculation = '{calc}'", f"  prefix = '{prefix}'", "  outdir = './tmp'",
        f"  pseudo_dir = '{PSEUDO_DIR}'", "  forc_conv_thr = 2.0d-3", "  etot_conv_thr = 1.0d-4",
        "  nstep = 200", "/",
        "&SYSTEM", "  ibrav = 0", f"  nat = {len(el)}", f"  ntyp = {len(species)}",
        f"  ecutwfc = {ecut[0]}", f"  ecutrho = {ecut[1]}", "  occupations = 'smearing'", "  smearing = 'mv'",
        "  degauss = 0.01", "  vdw_corr = 'dft-d3'", "  nosym = .true.", "/",
        "&ELECTRONS", "  conv_thr = 1.0d-7", "  mixing_beta = 0.2", "  mixing_mode = 'local-TF'",
        "  mixing_ndim = 10", "  electron_maxstep = 300", "/",
        "&IONS", "  ion_dynamics = 'bfgs'", "/",
        "ATOMIC_SPECIES"] + [f"  {e} {MASS[e]} {PP[e]}" for e in species] + [
        "CELL_PARAMETERS angstrom"] + [f"  {v[0]:.8f} {v[1]:.8f} {v[2]:.8f}" for v in cv] + [
        "ATOMIC_POSITIONS angstrom"]
    for e, x in zip(el, xyz):
        fix = fix_below is not None and e != "H" and x[2] < fix_below
        lines.append(f"  {e} {x[0]:.6f} {x[1]:.6f} {x[2]:.6f}" + ("  0 0 0" if fix else ""))
    lines.append("K_POINTS gamma" if kpts == "gamma" else f"K_POINTS automatic\n  {kpts} {kpts} 1 0 0 0")
    return "\n".join(lines) + "\n"


def cmd_inputs(write_runner=True):
    sel, sxyz, cv = slab()
    zmid = np.median(sxyz[np.array(sel) == "Ti", 2])
    fix = zmid + Z["Ti_out"] - 0.1            # atoms below the upper outer-Ti layer are fixed
    WORK.mkdir(parents=True, exist_ok=True)
    jobs = []
    (WORK / "slab").mkdir(exist_ok=True)
    (WORK / "slab" / "pw.in").write_text(pw_input("slab", sel, sxyz, cv, "gamma", fix))
    jobs.append("slab")
    for name in DRUGS:
        del_, dxyz = read_xyz(STARTS / f"{name}.xyz")
        d = WORK / f"drug_{name}"
        d.mkdir(exist_ok=True)
        dx = dxyz - dxyz.mean(0) + DBOX.sum(0) / 2
        (d / "pw.in").write_text(pw_input(f"d_{name}", del_, dx, DBOX, "gamma"))
        el, xyz = place(del_, dxyz, sel, sxyz, cv)
        c = WORK / f"cplx_{name}"
        c.mkdir(exist_ok=True)
        if (c / "restart_from.out").exists():   # continue a relaxation from its last geometry
            last = parse(c / "restart_from.out")
            assert last["el"] == el, "restart geometry does not match the complex"
            xyz = last["xyz"]
        (c / "pw.in").write_text(pw_input(f"c_{name}", el, xyz, cv, "gamma", fix))
        (c / "meta.json").write_text(json.dumps({"n_slab": len(sel), "drug_elements": del_}))
        jobs += [f"drug_{name}", f"cplx_{name}"]
    if not write_runner:                # never rewrite run_all.sh while bash is executing it
        print(f"wrote {len(jobs)} inputs to {WORK} (run_all.sh untouched)")
        return
    wsl_dir = "/mnt/c" + str(WORK).replace("\\", "/").split(":", 1)[1]
    sh = ["#!/bin/bash", "export PATH=$HOME/miniforge3/envs/qe/bin:$PATH", "export OMP_NUM_THREADS=1",
          f"NP=${{NP:-12}}", f"cd \"{wsl_dir}\""]
    for j in jobs:
        sh.append(f"if ! grep -q 'JOB DONE' {j}/pw.out 2>/dev/null; then (cd {j} && mpirun --oversubscribe -np $NP "
                  f"pw.x -nk 1 -in pw.in > pw.out 2>&1; rm -rf tmp); fi")
    (WORK / "run_all.sh").write_text("\n".join(sh) + "\n", newline="\n")
    print(f"wrote {len(jobs)} inputs to {WORK}")


def parse(out):
    t = Path(out).read_text(errors="ignore")
    E = [float(x) for x in re.findall(r"!\s+total energy\s+=\s+(-?\d+\.\d+)", t)]
    done = "JOB DONE" in t
    conv = "bfgs converged" in t or "End of BFGS" in t
    blocks = t.split("ATOMIC_POSITIONS (angstrom)")
    last = blocks[-1].split("End final coordinates")[0] if len(blocks) > 1 else None
    el, xyz = [], []
    if last:
        for l in last.strip().split("\n"):
            p = l.split()
            try:
                row = [float(v) for v in p[1:4]]
            except ValueError:
                row = None
            if len(p) >= 4 and p[0].isalpha() and row is not None:
                el.append(p[0])
                xyz.append(row)
            elif el:                       # first non-atom line after the block
                break
    return {"E_Ry": E[-1] if E else None, "done": done, "converged": conv, "el": el, "xyz": np.array(xyz)}


def bonds(el, xyz, idx_a, idx_b, cv):
    inv = np.linalg.inv(cv)
    out = set()
    for i in idx_a:
        for j in idx_b:
            if i >= j and idx_a is idx_b:
                continue
            d = xyz[j] - xyz[i]
            f = d @ inv
            f[:2] -= np.round(f[:2])
            r = np.linalg.norm(f @ cv)
            if r < 1.15 * (COV[el[i]] + COV[el[j]]):
                out.add((min(i, j), max(i, j)))
    return out


def min_contact(xyz, di, si, cv):
    """Closest drug heavy atom - slab atom distance with in-plane minimum image."""
    inv = np.linalg.inv(cv)
    dd = xyz[di][:, None] - xyz[si][None]
    f = dd @ inv
    f[..., :2] -= np.round(f[..., :2])
    return float(np.linalg.norm(f @ cv, axis=2).min())


def write_xyz(path, el, xyz, comment=""):
    Path(path).write_text("\n".join([str(len(el)), comment] +
                                    [f"{e} {x[0]:.6f} {x[1]:.6f} {x[2]:.6f}" for e, x in zip(el, xyz)]) + "\n")


def runner(name, jobs):
    wsl_dir = "/mnt/c" + str(WORK).replace("\\", "/").split(":", 1)[1]
    sh = ["#!/bin/bash", "export PATH=$HOME/miniforge3/envs/qe/bin:$PATH", "export OMP_NUM_THREADS=1",
          "NP=${NP:-12}", f'cd "{wsl_dir}"']
    for j in jobs:
        sh.append(f"if ! grep -q 'JOB DONE' {j}/pw.out 2>/dev/null; then (cd {j} && mpirun --oversubscribe -np $NP "
                  f"pw.x -nk 1 -in pw.in > pw.out 2>&1; rm -rf tmp); fi")
    (WORK / name).write_text("\n".join(sh) + "\n", newline="\n")


def finished(job):
    f = WORK / job / "pw.out"
    return f.exists() and parse(f)["done"]


def cmd_fragments():
    """Single points of the slab and of the drug frozen at their geometries in each
    relaxed complex, in the same cell (for dE_int: drug-image interactions cancel), and
    of the deformed drug alone in the isolated-drug box (for its deformation energy;
    the difference to the in-cell drug is the drug-image interaction)."""
    _, _, cv = cell()
    jobs = []
    for name in DRUGS:
        if not finished(f"cplx_{name}"):
            continue
        c = parse(WORK / f"cplx_{name}" / "pw.out")
        ns = json.loads((WORK / f"cplx_{name}" / "meta.json").read_text())["n_slab"]
        write_xyz(WORK / f"cplx_{name}" / "final.xyz", c["el"], c["xyz"], f"{name} on Ti3C2O2, PBE-D3 relaxed")
        for tag, idx in (("slab", slice(0, ns)), ("drug", slice(ns, None))):
            j = f"frag_{tag}_{name}"
            (WORK / j).mkdir(exist_ok=True)
            (WORK / j / "pw.in").write_text(pw_input(f"f{tag}_{name}", c["el"][idx], c["xyz"][idx], cv, "gamma",
                                                     calc="scf"))
            jobs.append(j)
        j = f"frag_drugbox_{name}"            # deformed drug alone: separates deformation from image interaction
        (WORK / j).mkdir(exist_ok=True)
        (WORK / j / "pw.in").write_text(pw_input(f"fdbox_{name}", c["el"][ns:], c["xyz"][ns:], DBOX, "gamma",
                                                 calc="scf"))
        jobs.append(j)
    runner("run_fragments.sh", jobs)
    print(f"{len(jobs)} fragment single points")


def cmd_checks(name="Temozolomide"):
    """Convergence checks of dE_ads on one complex: cutoff (60/480 vs 50/400 Ry,
    Gamma) and k-points (Gamma vs 2x2x1, both at 35/280 Ry to fit in memory)."""
    _, _, cv = cell()
    s_, d_, c_ = (parse(WORK / j / "pw.out") for j in ("slab", f"drug_{name}", f"cplx_{name}"))
    geo = {"slab": (s_["el"], s_["xyz"], cv), "drug": (d_["el"], d_["xyz"], DBOX),
           "cplx": (c_["el"], c_["xyz"], cv)}
    jobs = []
    for tag, ecut, kp in (("e60", (60, 480), "gamma"), ("e35g", (35, 280), "gamma"), ("e35k", (35, 280), 2)):
        for part in ("slab", "drug", "cplx"):
            if part == "drug" and tag == "e35k":
                continue                          # isolated molecule: k-independent, reuse e35g
            el, xyz, box = geo[part]
            j = f"chk_{tag}_{part}"
            (WORK / j).mkdir(exist_ok=True)
            (WORK / j / "pw.in").write_text(pw_input(f"{tag}{part}", el, xyz, box,
                                                     kp if part != "drug" else "gamma", calc="scf", ecut=ecut))
            jobs.append(j)
    runner("run_checks.sh", jobs)
    print(f"{len(jobs)} check single points for {name}")


def cmd_collect_checks(name="Temozolomide"):
    import pandas as pd
    jobs = ["chk_e60_slab", "chk_e60_drug", "chk_e60_cplx", "chk_e35g_slab", "chk_e35g_drug", "chk_e35g_cplx",
            "chk_e35k_slab", "chk_e35k_cplx"]
    E = {j: parse(WORK / j / "pw.out")["E_Ry"] for j in jobs}
    base = (parse(WORK / f"cplx_{name}" / "pw.out")["E_Ry"] - parse(WORK / "slab" / "pw.out")["E_Ry"]
            - parse(WORK / f"drug_{name}" / "pw.out")["E_Ry"])
    rows = [("50/400 Ry, Gamma (production)", base),
            ("60/480 Ry, Gamma", E["chk_e60_cplx"] - E["chk_e60_slab"] - E["chk_e60_drug"]),
            ("35/280 Ry, Gamma", E["chk_e35g_cplx"] - E["chk_e35g_slab"] - E["chk_e35g_drug"]),
            ("35/280 Ry, 2x2x1 k-points", E["chk_e35k_cplx"] - E["chk_e35k_slab"] - E["chk_e35g_drug"])]
    df = pd.DataFrame([{"setting": a, "dEads_kcal_mol": round(b * RY2KCAL, 2)} for a, b in rows])
    df.to_csv(BASE / "results" / "quantum" / "dft_convergence_checks.csv", index=False)
    print(df.to_string())


def cmd_collect():
    import pandas as pd
    _, _, cv = cell()
    s = parse(WORK / "slab" / "pw.out")
    rows = []
    for name in DRUGS:
        d = parse(WORK / f"drug_{name}" / "pw.out") if (WORK / f"drug_{name}" / "pw.out").exists() else None
        c = parse(WORK / f"cplx_{name}" / "pw.out") if (WORK / f"cplx_{name}" / "pw.out").exists() else None
        r = {"name": name}
        if d and c and s["E_Ry"] and d["E_Ry"] and c["E_Ry"] and c["done"]:
            ns = json.loads((WORK / f"cplx_{name}" / "meta.json").read_text())["n_slab"]
            el, xyz = c["el"], c["xyz"]
            di, si = list(range(ns, len(el))), list(range(ns))
            b0 = bonds(d["el"], d["xyz"], list(range(len(d["el"]))), list(range(len(d["el"]))), DBOX)
            b1 = {(i - ns, j - ns) for i, j in bonds(el, xyz, di, di, cv)}
            inter = sorted(f"{el[i]}{i}-{el[j]}{j}" for i, j in bonds(el, xyz, di, si, cv))
            heavy = [i for i in di if el[i] != "H"]
            D = min_contact(xyz, heavy, si, cv)
            write_xyz(WORK / f"cplx_{name}" / "final.xyz", el, xyz, f"{name} on Ti3C2O2, PBE-D3 relaxed")
            if finished(f"frag_slab_{name}") and finished(f"frag_drug_{name}"):
                fs, fd = (parse(WORK / f"frag_{t}_{name}" / "pw.out")["E_Ry"] for t in ("slab", "drug"))
                r["delta_Eint_kcal_mol"] = round((c["E_Ry"] - fs - fd) * RY2KCAL, 2)
                if finished(f"frag_drugbox_{name}"):
                    fb = parse(WORK / f"frag_drugbox_{name}" / "pw.out")["E_Ry"]
                    r["E_def_drug_kcal_mol"] = round((fb - d["E_Ry"]) * RY2KCAL, 2)
                    r["E_image_drug_kcal_mol"] = round((fd - fb) * RY2KCAL, 2)
                r["E_def_slab_kcal_mol"] = round((fs - s["E_Ry"]) * RY2KCAL, 2)
            r.update(delta_Eads_kcal_mol=round((c["E_Ry"] - s["E_Ry"] - d["E_Ry"]) * RY2KCAL, 2),
                     min_contact_A=round(float(D), 3), drug_intact=b0 == b1, drug_slab_bonds=inter,
                     adsorption_mode="chemisorption" if inter else "physisorption",
                     relax_converged=c["converged"] and d["converged"] and s["converged"])
        else:
            r["status"] = "incomplete"
        rows.append(r)
    pd.DataFrame(rows).to_csv(OUT, index=False)
    print(pd.DataFrame(rows).to_string())


if __name__ == "__main__":
    if sys.argv[1] == "inputs":
        cmd_inputs(write_runner="--keep-runner" not in sys.argv)
    else:
        {"collect": cmd_collect, "fragments": cmd_fragments, "checks": cmd_checks,
         "collect-checks": cmd_collect_checks}[sys.argv[1]]()
