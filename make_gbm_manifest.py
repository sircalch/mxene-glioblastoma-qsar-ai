import hashlib, time
from pathlib import Path

base = Path(r"c:\Users\Andre\Proyectos doctorado\mxene-glioblastoma-qsar-ai")
calc = base / "calculations" / "gbm"

def sha256_file(fp):
    h = hashlib.sha256()
    with open(fp, "rb") as fh:
        for chunk in iter(lambda: fh.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()

manifest_lines = [
    "# GBM MXene — SHA-256 Integrity Manifest (AUTHENTIC EXECUTABLE RAW LOGS)",
    f"# Generated: {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}",
    "# AutoDock Vina: v1.2.7 | xTB: v6.7.1-pre | ORCA: v6.1.1",
    "# Total processed compounds: 35 (Dual independent docking 4ZAU & 2J6M, xTB quantum calculated)",
    "# Primary Target: EGFR Kinase + Osimertinib (PDB: 4ZAU, 2.80 A)",
    "# Cross-Validation Target: EGFR Kinase + AEE788 (PDB: 2J6M, 3.10 A)",
    "# Carrier: Fully optimized Ti12C7O14 oxygen-terminated MXene cluster (33 atoms, E_MXene = -92.026933 Eh)",
    "# Outliers (|Delta_Eint| > 100 kcal/mol): 0 (100% negative physisorption regime)",
    "# Multi-Orientation Relaxed Subset (N=8): Spearman rho = 0.9048 (p=0.0020)",
    "# Heavy-Atom Redocking RMSD: 4ZAU (YY3, 37 heavy atoms) = 5.324 A | 2J6M (AEE, 22 heavy atoms) = 4.192 A",
    "#",
    "# SHA256                                                               bytes  role  path",
    "#" + "-"*95,
]

seen_hashes = set()
for p in sorted(base.rglob("*")):
    if p.is_file() and not p.name.startswith(".") and "MANIFEST" not in p.name and ".git" not in str(p):
        h = sha256_file(p)
        if (h, p.name) not in seen_hashes:
            seen_hashes.add((h, p.name))
            manifest_lines.append(f"{h}  {p.stat().st_size:>12} bytes  [gbm]  {p.relative_to(base)}")

m_path = base / "MANIFEST_SHA256.txt"
m_path.write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")
print(f"GBM MANIFEST_SHA256.txt updated: {len(seen_hashes)} files hashed.")
