"""
summarize_flake_tests.py - tabulates the xtb single points tried on the finite
Ti3C2O2 flakes (calculations/gbm/mxene_tests, mxene_build): method, charge,
spin, electronic temperature and whether the SCF converged. These tests are the
record of why a finite flake could not be used as the carrier model.

writes results/quantum/mxene_flake_scf_tests.csv
"""
import re
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parents[2]
ROOT = BASE / "calculations" / "gbm"


def arg(call, flag, default):
    m = re.search(rf"{flag}\s+(-?\S+)", call)
    return m.group(1) if m else default


def main():
    rows = []
    for f in sorted(list((ROOT / "mxene_tests").rglob("*.out")) + [ROOT / "mxene_build" / "opt.out"]):
        t = f.read_text(errors="ignore")
        m = re.search(r"program call\s+:\s+(.*)", t)
        if not m:
            continue
        call = m.group(1)
        gfn = "GFN-FF pre-opt + GFN2" if "--gfnff" in call else f"GFN{arg(call, '--gfn', '2')}"
        conv = "convergence criteria satisfied" in t and "abnormal termination" not in t
        rows.append({"test": str(f.relative_to(ROOT)).replace("\\", "/"),
                     "run": "opt" if "--opt" in call else "sp", "method": gfn,
                     "charge": int(arg(call, "--chrg", "0")), "uhf": int(arg(call, "--uhf", "0")),
                     "etemp_K": int(arg(call, "--etemp", "300")),
                     "extra": " ".join(x for x in ("--acc 0.1" if "--acc" in call else "",
                                                   "damped" if "damp.inp" in call else "",
                                                   "gasteiger guess" if "gasteiger" in call else "",
                                                   "restarted from dianion charges" if "restarted?                       true" in t
                                                   else "") if x),
                     "structure": ("earlier non-stoichiometric Ti12C7O14 cluster (discarded)"
                                   if "Ti12C7O14" in call else "Ti55C24O20 flake or smaller cut"),
                     "scf_converged": conv})
    df = pd.DataFrame(rows)
    df.to_csv(BASE / "results" / "quantum" / "mxene_flake_scf_tests.csv", index=False)
    print(df.to_string())


if __name__ == "__main__":
    main()
