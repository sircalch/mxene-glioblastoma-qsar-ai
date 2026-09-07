"""Reproduce the charge-density-difference (Delta-rho) cube for this system.

    Delta-rho(r) = rho_complex(r) - rho_carrier(r) - rho_drug(r)

all three densities evaluated with GFN2-xTB at the *complex* geometry, on one
shared grid. Requires xtb (GFN2-xTB) and Multiwfn 3.8 on PATH.

Inputs (this directory):
    {KEY}_complex.xyz   - GFN2-xTB optimised drug/2D-material complex
    {KEY}_carrier.xyz   - the 2D-material fragment, complex coordinates
    {KEY}_drug.xyz      - the drug fragment, complex coordinates
Output:
    {KEY}_deltarho.cub(.gz)

Fragments were separated by chemical identity (2D-material marker element vs the
rest); H atoms assigned to the nearer heavy atom. No geometry is changed.
"""
import os, sys, subprocess, shutil, gzip
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
KEY = os.environ.get("DRHO_KEY", "system")
XTB = shutil.which("xtb") or os.environ.get("XTB_EXE", "xtb")
MW = shutil.which("Multiwfn") or os.environ.get("MULTIWFN_EXE", "Multiwfn")
ETEMP = os.environ.get("DRHO_ETEMP", "400")


def xtb_molden(xyz, ns):
    for et in (ETEMP, "1500", "5000"):
        p = subprocess.run([XTB, os.path.basename(xyz), "--sp", "--molden",
                            "--etemp", et, "--iterations", "500", "--namespace", ns],
                           cwd=HERE, capture_output=True, text=True)
        m = os.path.join(HERE, ns + ".molden.input")
        if p.returncode == 0 and os.path.exists(m):
            return m
    raise SystemExit(f"xtb failed for {xyz}")


def mw_cube(molden, out_name, grid_from=None):
    stdin = (f"\n5\n1\n8\n{grid_from}\n2\n0\nq\n" if grid_from
             else "\n5\n1\n2\n2\n0\nq\n")
    subprocess.run([MW, os.path.basename(molden)], cwd=HERE, input=stdin,
                   capture_output=True, text=True)
    src = os.path.join(HERE, "density.cub")
    if not os.path.exists(src):
        raise SystemExit(f"Multiwfn produced no cube for {molden}")
    shutil.move(src, os.path.join(HERE, out_name))


def read_cube(fn):
    L = open(fn).read().splitlines()
    na = int(L[2].split()[0])
    n = np.array([int(L[3 + i].split()[0]) for i in range(3)])
    return L[2:6 + na], n, np.array(" ".join(L[6 + na:]).split(), float).reshape(n)


def main():
    key = sys.argv[1] if len(sys.argv) > 1 else KEY
    mcx = xtb_molden(os.path.join(HERE, f"{key}_complex.xyz"), "cx")
    mca = xtb_molden(os.path.join(HERE, f"{key}_carrier.xyz"), "ca")
    mdr = xtb_molden(os.path.join(HERE, f"{key}_drug.xyz"), "dr")
    mw_cube(mcx, "rho_complex.cub")
    mw_cube(mca, "rho_carrier.cub", grid_from="rho_complex.cub")
    mw_cube(mdr, "rho_drug.cub", grid_from="rho_complex.cub")
    h, n, dc = read_cube(os.path.join(HERE, "rho_complex.cub"))
    _, _, da = read_cube(os.path.join(HERE, "rho_carrier.cub"))
    _, _, db = read_cube(os.path.join(HERE, "rho_drug.cub"))
    diff = dc - da - db
    nx, ny, nz = n
    body = [" ".join("%12.5E" % v for v in diff[i, j, k:k + 6])
            for i in range(nx) for j in range(ny) for k in range(0, nz, 6)]
    out = os.path.join(HERE, f"{key}_deltarho.cub")
    with open(out, "w") as f:
        f.write(" Delta-rho = rho(complex) - rho(carrier) - rho(drug)\n GFN2-xTB\n")
        f.write("\n".join(h) + "\n" + "\n".join(body) + "\n")
    with open(out, "rb") as fi, gzip.open(out + ".gz", "wb") as fo:
        shutil.copyfileobj(fi, fo)
    dv = abs(np.linalg.det(np.array([[float(x) for x in h[1 + i].split()[1:4]]
                                     for i in range(3)])))
    print(f"wrote {out}(.gz)  transferred +{diff[diff>0].sum()*dv:.3f} e")


if __name__ == "__main__":
    main()
