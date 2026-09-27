"""
xtb_lattice_scan.py - energy of a periodic 4x4 Ti3C2O2 slab versus the in-plane
lattice constant a, with GFN2-xTB and GFN1-xTB (tblite, Gamma point, 1500 K
Fermi smearing). Documents why the tight-binding methods were not used for the
MXene carrier: the DFT/experimental value is a = 3.03 A.

run in the 'tb' environment:
  micromamba run -r C:/Users/Andre/mm -n tb python src/quantum/xtb_lattice_scan.py
writes results/quantum/xtb_lattice_scan.csv
"""
import sys
from pathlib import Path

import numpy as np
import csv
from tblite.ase import TBLite

sys.path.insert(0, str(Path(__file__).resolve().parent))
import periodic_adsorption as P  # noqa: E402

BASE = Path(__file__).resolve().parents[2]


def main():
    rows = []
    for method in ("GFN2-xTB", "GFN1-xTB"):
        for a in np.round(np.arange(2.60, 3.251, 0.05), 3):
            at = P.slab(a, 4)
            at.calc = TBLite(method=method, electronic_temperature=1500, max_iterations=500, verbosity=0)
            try:
                e = at.get_potential_energy() / 16
                st = "ok"
            except Exception as exc:      # SCF not converged
                e, st = np.nan, str(exc)[:60]
            rows.append({"method": method, "a_A": a, "E_eV_per_fu": e, "status": st})
            print(method, a, e, st, flush=True)
    with open(BASE / "results" / "quantum" / "xtb_lattice_scan.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
