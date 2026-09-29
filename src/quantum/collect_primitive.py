"""
collect_primitive.py - lattice constant and metallic character of Ti3C2O2 from the
PBE-D3 calculations on the primitive 1x1 cell (same settings as the slab campaign:
SSSP 1.3, 50/400 Ry, Marzari-Vanderbilt 0.01 Ry, D3; vc-relax 2Dxy at 12x12x1 k,
then SCF + tetrahedron NSCF at 24x24x1 k and dos.x).

inputs : calculations/gbm_dft/lattice_primitive/{pw.out, dos.dat}
         (run with pw.in / run.sh and scf.in, nscf.in, dos.in / run_dos.sh there)
writes : results/quantum/ti3c2o2_primitive_pbe_d3.json
"""
import json
import re
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parents[2]
D = BASE / "calculations" / "gbm_dft" / "lattice_primitive"


def main():
    out = (D / "pw.out").read_text(encoding="utf-8", errors="replace")
    assert "bfgs converged" in out and "JOB DONE" in out
    blk = out[out.rfind("Begin final coordinates"):out.rfind("End final coordinates")]
    a1 = [float(x) for x in re.search(r"CELL_PARAMETERS \(angstrom\)\n\s*(\S+)\s+(\S+)\s+(\S+)", blk).groups()]
    a = float(np.linalg.norm(a1))
    p = float(re.findall(r"P=\s*(-?\d+\.\d+)", out)[-1])
    first = (D / "dos.dat").open().readline()
    ef = float(re.search(r"EFermi\s*=\s*(-?\d+\.\d+)", first).group(1))
    dos = np.loadtxt(D / "dos.dat")
    at_ef = float(np.interp(ef, dos[:, 0], dos[:, 1]))
    rec = {"a_A": round(a, 4), "residual_pressure_kbar": p, "E_Fermi_eV": ef,
           "DOS_at_EF_states_per_eV_per_cell": round(at_ef, 3), "metallic": at_ef > 0.1,
           "method": "PBE-D3, SSSP 1.3 efficiency, 50/400 Ry, MV 0.01 Ry; vc-relax 2Dxy 12x12x1; DOS tetrahedra 24x24x1"}
    out_f = BASE / "results" / "quantum" / "ti3c2o2_primitive_pbe_d3.json"
    out_f.write_text(json.dumps(rec, indent=2))
    print(rec)


if __name__ == "__main__":
    main()
