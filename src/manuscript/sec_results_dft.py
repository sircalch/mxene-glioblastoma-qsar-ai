"""sec_results_dft.py - periodic DFT adsorption results of the GBM manuscript.

Every number comes from results/quantum/dft_adsorption.csv and
dft_convergence_checks.csv (written by src/quantum/qe_mxene_dft.py collect and
collect-checks). While the campaign is running, a marked placeholder is written
instead, so that the manuscript can be assembled at any time.
"""
import ast

import pandas as pd

import docx_kit as k


def f1(x):
    return f"{x:.1f}".replace("-", "−")


def f2(x):
    return "0.00" if abs(x) < 0.005 else f"{x:.2f}".replace("-", "−")


def dn(n):
    return n[0].lower() + n[1:] if n[0].isupper() and n[1:2].islower() else n


def complete(t):
    """True when every drug has a converged relaxation and its fragment energies."""
    need = {"delta_Eads_kcal_mol", "delta_Eint_kcal_mol", "E_def_drug_kcal_mol", "E_def_slab_kcal_mol"}
    return (t is not None and need <= set(t.columns) and t[list(need)].notna().all().all()
            and t.relax_converged.astype(bool).all())


def bond_types(x):
    """"['O120-Ti40']" -> 'O–Ti'."""
    import re
    out = []
    for b in ast.literal_eval(x) if isinstance(x, str) else []:
        a, c = (re.sub(r"\d", "", u) for u in b.split("-"))
        if f"{a}–{c}" not in out:
            out.append(f"{a}–{c}")
    return ", ".join(out)


def results_dft(doc, d, c, fig_path, conv):
    t = d["dft"]
    k.heading(doc, "Periodic DFT adsorption of the alkylating agents", 2)
    if not complete(t):
        k.placeholder(doc, "[DFT RESULTS PENDING: the PBE-D3 campaign (slab, drugs, complexes, fragments and "
                           "convergence checks) is running. This subsection, Fig. 5 and the DFT sentences of the "
                           "abstract and conclusions are generated from results/quantum/dft_adsorption.csv once "
                           "all four complexes are complete.]")
        return None
    t = t.set_index("name").sort_values("delta_Eads_kcal_mol")
    chem = t[t.adsorption_mode == "chemisorption"]
    react = t[~t.drug_intact.astype(bool)]
    n = len(t)
    s = {"t": t, "chem": chem, "react": react}
    mode = ("All four drugs physisorbed: no drug–surface bond formed" if len(chem) == 0 else
            f"{len(chem)} of the {n} drugs formed a bond to the surface "
            f"({'; '.join(f'{dn(i)}: {bond_types(r.drug_slab_bonds)}' for i, r in chem.iterrows())})")
    intact = ("every drug kept its own bonding" if len(react) == 0 else
              f"{', '.join(dn(i) for i in react.index)} changed {'its' if len(react) == 1 else 'their'} own bonding")
    k.para(doc,
           f"{mode}, and {intact} (Fig. 5). The adsorption energies range from "
           f"{f1(t.delta_Eads_kcal_mol.max())} to {f1(t.delta_Eads_kcal_mol.min())} kcal mol^{{−1}}, in the order "
           + ", ".join(f"{dn(i)} ({f1(v)})" for i, v in t.delta_Eads_kcal_mol.items()) +
           ". The closest drug–surface heavy-atom contacts are "
           f"{f2(t.min_contact_A.min())}–{f2(t.min_contact_A.max())} Å (Table 1; per-drug details in Table S6).",
           indent=True)
    k.para(doc,
           "The interaction energies, which refer to both partners frozen at their geometries in the complex, "
           f"range from {f1(t.delta_Eint_kcal_mol.max())} to {f1(t.delta_Eint_kcal_mol.min())} kcal mol^{{−1}}. "
           "The difference to Δ*E*_{ads} is the deformation of the two partners and the lateral interaction "
           "of the drug with its periodic images. The deformation energy of the drugs is "
           f"{f1(t.E_def_drug_kcal_mol.min())}–{f1(t.E_def_drug_kcal_mol.max())} kcal mol^{{−1}}, that of the "
           f"surface {f1(t.E_def_slab_kcal_mol.min())}–{f1(t.E_def_slab_kcal_mol.max())} kcal mol^{{−1}}"
           + (f", and the drug–image interaction {f1(t.E_image_drug_kcal_mol.min())} to "
              f"{f1(t.E_image_drug_kcal_mol.max())} kcal mol^{{−1}}" if "E_image_drug_kcal_mol" in t else "") +
           ".", indent=True)
    if conv is not None:
        cv = conv.set_index("setting").dEads_kcal_mol
        prod = cv["50/400 Ry, Gamma (production)"]
        k.para(doc,
               "For temozolomide, raising the cutoffs to 60/480 Ry changes Δ*E*_{ads} by "
               f"{f2(cv['60/480 Ry, Gamma'] - prod)} kcal mol^{{−1}}, and replacing the Γ point by a 2×2×1 grid "
               f"(at 35/280 Ry) by {f2(cv['35/280 Ry, 2x2x1 k-points'] - cv['35/280 Ry, Gamma'])} kcal "
               "mol^{−1} (Table S3).", indent=True)
    k.placeholder(doc, "[AUTHOR/ANALYSIS TO COMPLETE AFTER INSPECTING THE RELAXED COMPLEXES: binding geometry "
                       "(which groups face the surface), comparison between the drugs, and what the values "
                       "imply for loading and release.]")
    if fig_path.exists():
        top = ", ".join(t.index[:3])
        k.figure(doc, fig_path, 5,
                 "PBE-D3 adsorption of the alkylating agents on Ti_{3}C_{2}O_{2}. **a** −Δ*E*_{ads} of the four "
                 "drugs (relaxed complex against the separately relaxed slab and drug), coloured by regime. "
                 f"**b**–**d** Relaxed complexes of {top} (side views). Colours: Ti silver, C grey, O red, "
                 "N blue, Cl green, H white")
    rows = [[i, f1(r.delta_Eads_kcal_mol), f1(r.delta_Eint_kcal_mol), f1(r.E_def_drug_kcal_mol),
             f1(r.E_def_slab_kcal_mol), f2(r.min_contact_A), r.adsorption_mode,
             "yes" if bool(r.drug_intact) else "no"] for i, r in t.iterrows()]
    k.table(doc, (1, "PBE-D3 adsorption of the alkylating agents on the 4×4 Ti_{3}C_{2}O_{2} slab."),
            ["Drug", "Δ*E*_{ads}", "Δ*E*_{int}", "*E*_{def,drug}", "*E*_{def,slab}", "*d*_{min}", "Regime",
             "Drug intact"], rows, align="lccccclc", font=8,
            note="Energies in kcal mol^{−1}. *d*_{min}, closest drug–surface heavy-atom contact (Å). Γ point, "
                 "50/400 Ry; only the upper O and outer Ti layers and the drug relaxed.")
    return s
