"""
build_mxene_flake.py
====================
Builds a stoichiometric, O-terminated Ti3C2O2 MXene flake cut from the
crystalline lattice and relaxes it with GFN2-xTB.

Lattice (P-3m1, a = 3.03 A): middle Ti at A (0,0,0); C at B (+z) and C (-z)
hollows, outer Ti at C (+z) and B (-z); O terminations above/below the middle
Ti (the most stable "hollow over middle Ti" site). Layer heights: C +-1.15 A,
outer Ti +-2.33 A, O +-3.30 A.

The flake is cut around carbon: every C within a 5.3 A in-plane radius (both
C layers), the Ti atoms of their Ti6C octahedra, and O over every hollow with
all three outer-Ti neighbours (Ti55C24O20, 99 atoms). Every C is therefore
six-coordinate; the edges are Ti-rich, as in exfoliated flakes. MXenes are
metallic, so GFN2-xTB is run with Fermi smearing (electronic temperature
ETEMP K). The neutral flake did not reach SCF convergence at 300-3000 K, in
several spin states, with damping or other initial guesses, nor with GFN1-xTB;
a single point converged only when restarted from the charges of the dianion,
and the geometry optimisation never converged (results/quantum/
mxene_flake_scf_tests.csv). The finite flake was therefore abandoned: see
periodic_adsorption.py, xtb_lattice_scan.py and qe_mxene_dft.py.

Replaces the earlier 33-atom "Ti12C7O14" cluster, which was not
stoichiometric and contained C-C bonds.

Validation: no C-C, O-O or C-O bonds; each C bonded to six Ti.
Output: calculations/gbm/Ti3C2O2_flake_optimized.xyz, Ti3C2O2_flake.json
"""
import itertools
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parents[2]
OUT = BASE / "calculations" / "gbm"
A = 3.03
Z = {"C": 1.15, "Ti_out": 2.33, "O": 3.30}
RADIUS = 5.3                  # A, in-plane radius of the carbon cut -> Ti55C24O20
ETEMP = 1500
COV = {"Ti": 1.60, "C": 0.76, "O": 0.66}


def flake(radius=None):
    """Carbon-centred cut: every C of both C layers within `radius` (in-plane)
    of the axis, the Ti atoms of their Ti6C octahedra, and O over every hollow
    whose three outer-Ti neighbours are present."""
    radius = radius or RADIUS
    a1 = np.array([A, 0, 0])
    a2 = np.array([A / 2, A * np.sqrt(3) / 2, 0])
    B = (a1 + a2) / 3
    C = 2 * (a1 + a2) / 3
    cells = [n1 * a1 + n2 * a2 for n1, n2 in itertools.product(range(-8, 9), repeat=2)]
    all_ti = np.array(cells + [o + C + [0, 0, Z["Ti_out"]] for o in cells]
                      + [o + B - [0, 0, Z["Ti_out"]] for o in cells])
    all_c = np.array([o + B + [0, 0, Z["C"]] for o in cells] + [o + C - [0, 0, Z["C"]] for o in cells])
    cs = all_c[np.linalg.norm(all_c[:, :2], axis=1) <= radius]
    ti = np.array([t for t in all_ti if (np.linalg.norm(cs - t, axis=1) < 2.4).any()])
    ti_out = ti[np.abs(ti[:, 2]) > 1]
    ox = [o + [0, 0, s * Z["O"]] for o in cells for s in (1, -1)
          if (np.linalg.norm(ti_out - (o + [0, 0, s * Z["O"]]), axis=1) < 2.3).sum() == 3]
    el = ["Ti"] * len(ti) + ["C"] * len(cs) + ["O"] * len(ox)
    return el, np.vstack([ti, cs, np.array(ox)])


def census(el, X):
    c = {}
    for i, j in itertools.combinations(range(len(el)), 2):
        if np.linalg.norm(X[i] - X[j]) < 1.2 * (COV[el[i]] + COV[el[j]]):
            k = "-".join(sorted((el[i], el[j])))
            c[k] = c.get(k, 0) + 1
    return c


def write_xyz(p, el, X, comment=""):
    with open(p, "w") as f:
        f.write(f"{len(el)}\n{comment}\n")
        for e, x in zip(el, X):
            f.write(f"{e:2s} {x[0]:15.8f} {x[1]:15.8f} {x[2]:15.8f}\n")


def read_xyz(p):
    L = Path(p).read_text().splitlines()
    n = int(L[0])
    return [l.split()[0] for l in L[2:2 + n]], np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]])


def validate(el, X, label):
    c = census(el, X)
    bad = {"C-C", "O-O", "C-O"} & set(c)
    if bad:
        sys.exit(f"{label}: forbidden bonds {bad}")
    for i, e in enumerate(el):
        if e == "C":
            nti = sum(1 for j, f in enumerate(el) if f == "Ti" and np.linalg.norm(X[i] - X[j]) < 2.6)
            if nti != 6:
                sys.exit(f"{label}: C{i} has {nti} Ti neighbours")
    return c


def main():
    el, X = flake()
    print(f"{len(el)} atoms: " + ", ".join(f"{e}{el.count(e)}" for e in ("Ti", "C", "O")))
    print("initial census", validate(el, X, "initial"))
    xtb = os.environ.get("XTB_EXE") or shutil.which("xtb") or "C:/Users/Andre/mm/xtb/Library/bin/xtb.exe"
    env = dict(os.environ)
    env["XTBPATH"] = str(Path(xtb).parent.parent / "share" / "xtb")
    env.setdefault("OMP_NUM_THREADS", "4")
    wd = OUT / "mxene_build"
    wd.mkdir(parents=True, exist_ok=True)
    write_xyz(wd / "start.xyz", el, X, "Ti3C2O2 flake from lattice")
    # The neutral metallic flake does not reach SCF convergence from the default guess; the
    # dianion does. Its converged charges (xtbrestart) seed the neutral calculation.
    subprocess.run([xtb, "start.xyz", "--sp", "--gfn", "2", "--chrg", "-2", "--etemp", str(ETEMP),
                    "--iterations", "1000"], cwd=wd, env=env, capture_output=True, text=True, errors="replace")
    out = subprocess.run([xtb, "start.xyz", "--opt", "--gfn", "2", "--etemp", str(ETEMP), "--iterations", "1500",
                          "--cycles", "1500"], cwd=wd, env=env, capture_output=True, text=True,
                         errors="replace").stdout
    (wd / "opt.out").write_text(out, encoding="utf-8")
    if "GEOMETRY OPTIMIZATION CONVERGED" not in out:
        sys.exit("optimisation did not converge")
    el2, X2 = read_xyz(wd / "xtbopt.xyz")
    c = validate(el2, X2, "relaxed")
    formula = "".join(f"{e}{el2.count(e)}" for e in ("Ti", "C", "O"))
    e = [float(l.split()[3]) for l in out.splitlines() if "TOTAL ENERGY" in l][-1]
    gap = [float(l.split()[3]) for l in out.splitlines() if "HOMO-LUMO GAP" in l][-1]
    write_xyz(OUT / "Ti3C2O2_flake_optimized.xyz", el2, X2,
              f"{formula} flake, GFN2-xTB opt, etemp {ETEMP} K, E = {e:.8f} Eh")
    (OUT / "Ti3C2O2_flake.json").write_text(json.dumps(
        {"E_Eh": e, "gap_eV": gap, "etemp_K": ETEMP, "bond_census": c, "formula": formula}, indent=2))
    print(f"relaxed: E = {e:.6f} Eh, gap {gap:.3f} eV, census {c}")


if __name__ == "__main__":
    main()
