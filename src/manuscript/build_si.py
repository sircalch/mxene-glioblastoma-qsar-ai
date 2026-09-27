"""
build_si.py - Supporting Information (Online Resource 1, Word; converted to PDF for
submission) of the GBM / Ti3C2O2 study, generated from the pipeline outputs.

writes manuscript/submission/Supporting_Information_GBM_MXene.docx
Tables are numbered in the order of their first citation in the manuscript:
  S1 compounds, S2 finite-flake xtb tests, S3 DFT convergence checks,
  S4 every redocking mode, S5 periodic lattice scan, S6 per-drug DFT details.
S3 and S6 depend on the DFT campaign and are added when their result files exist.
"""
import ast
import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import docx_kit as k  # noqa: E402
from build_manuscript import AFFIL, AUTHOR, EMAIL, TITLE, load  # noqa: E402
from sec_results_dft import complete  # noqa: E402

BASE = HERE.parents[1]
OUT = BASE / "manuscript" / "submission"
Q = BASE / "results" / "quantum"


def num(x, nd=2):
    return "–" if pd.isna(x) else f"{x:.{nd}f}".replace("-", "−")


def formula(f):
    return re.sub(r"(\d+)", r"_{\1}", str(f))


def bond_types(x):
    """"['O120-Ti40', 'N118-Ti9']" -> 'O–Ti, N–Ti'."""
    out = []
    for b in ast.literal_eval(x) if isinstance(x, str) else []:
        a, c = (re.sub(r"\d", "", u) for u in b.split("-"))
        if f"{a}–{c}" not in out:
            out.append(f"{a}–{c}")
    return ", ".join(out) or "–"


def main():
    d = load()
    m = d["m"]
    lib = pd.read_csv(BASE / "data" / "processed" / "compound_library_pubchem.csv")
    doc = k.new_document()
    k.si_header(doc, TITLE, "Journal of Molecular Modeling", AUTHOR, AFFIL, EMAIL)

    # S1 compounds
    fam = m.set_index("name")
    rows = [[r.name, fam.loc[r.name, "family"], str(r.pubchem_cid), formula(r.formula),
             num(fam.loc[r.name, "vina_1M17_kcal_mol"]), num(fam.loc[r.name, "ligand_efficiency"], 3)]
            for r in lib.sort_values("name").itertuples()]
    k.table(doc, ("S1", "Compounds (structures from PubChem, identity checked by InChIKey) with their Vina score "
                        "in EGFR (PDB 1M17, kcal mol^{−1}) and ligand efficiency (kcal mol^{−1} per heavy atom). "
                        "Bortezomib was not docked: AutoDock Vina has no boron atom type."),
            ["Compound", "Family", "PubChem CID", "Formula", "Vina", "LE"], rows, align="llllcc", font=7.5)

    # S2 finite-flake tests
    ft = pd.read_csv(Q / "mxene_flake_scf_tests.csv")
    rows = [[r.structure, r.run, r.method, str(r.charge), str(r.uhf), str(r.etemp_K),
             r.extra if isinstance(r.extra, str) else "", "yes" if r.scf_converged else "no"]
            for r in ft.itertuples()]
    k.table(doc, ("S2", "xtb calculations attempted on finite Ti_{3}C_{2}O_{2} flakes. uhf, number of unpaired "
                        "electrons; *T*_{el}, electronic temperature (K)."),
            ["Structure", "Run", "Method", "Charge", "uhf", "*T*_{el}", "Options", "SCF converged"],
            rows, align="lllccccc", font=7)

    # S3 DFT convergence checks
    f = Q / "dft_convergence_checks.csv"
    if f.exists():
        ch = pd.read_csv(f)
        rows = [[r.setting.replace("x", "×").replace("Gamma", "Γ point"), num(r.dEads_kcal_mol)] for r in ch.itertuples()]
        k.table(doc, ("S3", "Convergence of the PBE-D3 adsorption energy of temozolomide with the plane-wave cutoffs "
                            "(wavefunctions/density) and the k-point sampling, as single points on the relaxed "
                            "geometries. The k-point test is run at 35/280 Ry to fit in memory and is compared with "
                            "the Γ point at the same cutoffs."),
                ["Setting", "Δ*E*_{ads} (kcal mol^{−1})"], rows, align="lc")
    else:
        k.placeholder(doc, "[Table S3 (DFT convergence checks) is generated when "
                           "results/quantum/dft_convergence_checks.csv exists.]")

    # S4 every redocking mode
    md = d["modes"]
    rows = [[r.control, str(r.mode), num(r.vina_kcal_mol, 3), num(r.rmsd_A)] for r in md.itertuples()]
    k.table(doc, ("S4", "Every output mode of the two erlotinib redocking controls in PDB 1M17: Vina score and "
                        "heavy-atom RMSD to the crystal pose (symmetry-aware, no re-alignment)."),
            ["Control", "Mode", "Vina (kcal mol^{−1})", "RMSD (Å)"], rows, align="lccc", font=7.5)

    # S5 periodic lattice scan
    sc = pd.read_csv(Q / "xtb_lattice_scan.csv")
    g1 = sc[sc.method == "GFN1-xTB"]
    wild = g1[g1.E_eV_per_fu < -600]
    ref = g1[(g1.E_eV_per_fu > -600)].E_eV_per_fu.median()
    rows = [[r.method, f"{r.a_A:.2f}", num(r.E_eV_per_fu, 3), r.status] for r in sc.itertuples()]
    note = ""
    if len(wild):
        note = (f" The GFN1-xTB values at {' and '.join(f'{a:.2f}' for a in wild.a_A)} Å "
                f"({num(wild.E_eV_per_fu.min(), 0)} to {num(wild.E_eV_per_fu.max(), 0)} eV per f.u.) are unphysical, "
                f"about {wild.E_eV_per_fu.mean() / ref:.1f} times the energy at the other lattice constants.")
    k.table(doc, ("S5", "Energy of the periodic 4×4 Ti_{3}C_{2}O_{2} slab per formula unit against the in-plane "
                        "lattice constant, GFN2-xTB and GFN1-xTB (tblite, Γ point, 1500 K)." + note),
            ["Method", "*a* (Å)", "*E* (eV per f.u.)", "Status"], rows, align="lccl", font=7.5)

    # S6 per-drug DFT details
    t = d["dft"]
    if complete(t):
        rows = []
        for r in t.itertuples():
            rows.append([r.name, num(r.delta_Eads_kcal_mol, 1), num(r.delta_Eint_kcal_mol, 1),
                         num(r.E_def_drug_kcal_mol, 1), num(r.E_def_slab_kcal_mol, 1),
                         num(getattr(r, "E_image_drug_kcal_mol", float("nan")), 1), num(r.min_contact_A),
                         bond_types(r.drug_slab_bonds), "yes" if bool(r.drug_intact) else "no"])
        k.table(doc, ("S6", "PBE-D3 adsorption of each alkylating agent on the 4×4 Ti_{3}C_{2}O_{2} slab. "
                            "Δ*E*_{ads} = Δ*E*_{int} + *E*_{def,drug} + *E*_{def,slab} + *E*_{image}, where "
                            "*E*_{image} is the interaction of the deformed drug with its own periodic images "
                            "(in-cell drug minus the same drug in a 22 Å box). Energies in kcal mol^{−1}; "
                            "*d*_{min}, closest drug–surface heavy-atom contact (Å); bonds, drug–surface pairs "
                            "closer than 1.15 times the sum of the covalent radii."),
                ["Drug", "Δ*E*_{ads}", "Δ*E*_{int}", "*E*_{def,drug}", "*E*_{def,slab}", "*E*_{image}",
                 "*d*_{min}", "Bonds", "Intact"], rows, align="lcccccclc", font=7)
    else:
        k.placeholder(doc, "[Table S6 (per-drug DFT details) is generated when the DFT campaign is complete.]")

    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / "Supporting_Information_GBM_MXene.docx"
    doc.save(out)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
