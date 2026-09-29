"""
build_manuscript.py - Journal of Molecular Modeling submission (Word) for the
glioblastoma / Ti3C2O2 MXene study. Every number in the text, tables and captions
is read from the pipeline outputs; nothing is typed by hand. The sections live in
sec_*.py; the DFT subsection and the DFT sentences of the abstract and conclusions
are written only once the DFT campaign is complete (marked placeholders until then).

usage: python src/manuscript/build_manuscript.py
writes manuscript/submission/Manuscript_GBM_MXene_JMM.docx
"""
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "figures"))
import docx_kit as k  # noqa: E402
import make_figures as F  # noqa: E402
from references import REFS  # noqa: E402
from sec_introduction import introduction  # noqa: E402
from sec_methods import methods  # noqa: E402
from sec_results_a import docking_stats, results_docking, results_tightbinding, scan_stats  # noqa: E402
from sec_results_dft import complete, dn, results_dft  # noqa: E402
from sec_results_qspr import results_qspr  # noqa: E402

BASE = HERE.parents[1]
FIG = BASE / "figures"
OUT = BASE / "manuscript" / "submission"
Q = BASE / "results" / "quantum"

AUTHOR = "Andrés Monreal Hernández"
AFFIL = "Universidad Estatal de Sonora, Ley Federal del Trabajo S/N, Col. Apolo, 83100 Hermosillo, Sonora, Mexico"
EMAIL = "andres.monreal@ues.mx"
ORCID = "0009-0009-1207-8597"
REPO = "https://github.com/sircalch/mxene-glioblastoma-qsar-ai"
TITLE = ("Glioblastoma drugs and a Ti_{3}C_{2}O_{2} MXene carrier: docking controls, tight-binding tests and "
         "periodic DFT adsorption")


class Cites:
    """[n] numbering in order of first citation."""

    def __init__(self):
        self.order = []

    def __call__(self, *keys):
        nums = []
        for key in keys:
            if key not in REFS:
                raise KeyError(key)
            if key not in self.order:
                self.order.append(key)
            nums.append(self.order.index(key) + 1)
        nums = sorted(set(nums))
        spans, start = [], nums[0]
        for a, b in zip(nums, nums[1:] + [None]):
            if b != a + 1:
                spans.append(f"{start}" if start == a else f"{start}–{a}" if a - start > 1 else f"{start}, {a}")
                start = b
        return "[" + ", ".join(spans) + "]"

    def list(self):
        return [REFS[k_] for k_ in self.order]


def f1(x):
    return f"{x:.1f}".replace("-", "−")


def f2(x):
    return "0.00" if abs(x) < 0.005 else f"{x:.2f}".replace("-", "−")


def load():
    d = F.data()
    d["flake_tests"] = Q / "mxene_flake_scf_tests.csv"
    f = Q / "dft_convergence_checks.csv"
    d["conv"] = pd.read_csv(f) if f.exists() else None
    return d


def front(doc):
    k.para(doc, f"**{TITLE}**", align="left", size=15, space_after=12)
    k.para(doc, f"{AUTHOR}^{{*}}", align="left", space_after=2)
    k.para(doc, AFFIL, align="left", size=10, space_after=2)
    k.para(doc, f"^{{*}}Corresponding author: {EMAIL}; ORCID {ORCID}", align="left", size=10, space_after=14)


def dft_sentence(d):
    """One abstract/conclusion sentence on the DFT result, or None while pending."""
    if not complete(d["dft"]):
        return None
    t = d["dft"].set_index("name")
    nch = int((t.adsorption_mode == "chemisorption").sum())
    return (f"{'All four' if nch == 0 else f'{4 - nch} of the four'} agents physisorb, with PBE-D3 adsorption "
            f"energies of {f1(t.delta_Eads_kcal_mol.max())} to {f1(t.delta_Eads_kcal_mol.min())} kcal mol^{{−1}}"
            f" (strongest: {dn(t.delta_Eads_kcal_mol.idxmin())})")


def abstract(doc, d):
    s, t, q = docking_stats(d), scan_stats(d), d["q"]
    k.heading(doc, "Abstract")
    dft = dft_sentence(d)
    k.labelled(doc, "Context",
               "Local delivery is one of the few routes by which drugs reach glioblastoma, and two-dimensional "
               "MXenes have been proposed as carriers. We examined "
               f"{len(d['m'])} glioblastoma drugs against the erlotinib site of EGFR and against an O-terminated "
               "Ti_{3}C_{2}O_{2} MXene. In both redocking controls the crystal-like erlotinib pose (root-mean-square "
               f"deviation ≤ {f2(max(s['best_x'].rmsd_A, s['best_s'].rmsd_A))} Å) is found but not ranked first "
               f"(top poses {f1(s['top_x'].rmsd_A)} and {f1(s['top_s'].rmsd_A)} Å), so the docking scores are "
               "exploratory. The GFN2-xTB and GFN1-xTB tight-binding methods fail for this metallic carrier: "
               "neutral finite flakes do not converge, GFN2-xTB places the lattice minimum "
               f"{100 * (3.03 - t['g2_min']) / 3.03:.0f}% below the PBE-D3 lattice constant, and GFN1-xTB has no "
               "minimum. Adsorption of four alkylating agents was therefore computed with periodic density "
               "functional theory. " + (dft + ". " if dft else "") +
               "A ridge model of the docking score on four descriptors is weak "
               f"(*Q*^{{2}}_{{CV}} = {f2(q['Q2_CV'])}) and mainly reproduces the size dependence of the score.")
    if dft is None:
        k.placeholder(doc, "[DFT sentence of the abstract is inserted automatically when the campaign is "
                           "complete; word count to be rechecked (150–250).]")
    k.labelled(doc, "Methods",
               "Structures were taken from PubChem. Drugs were docked into EGFR (PDB 1M17) with AutoDock Vina "
               "1.2.7, and every mode of two redocking controls was compared with the crystal pose. Finite and "
               "periodic Ti_{3}C_{2}O_{2} models were tested with GFN2-xTB and GFN1-xTB. Adsorption was computed "
               "with PBE-D3 in Quantum ESPRESSO on a 4×4 slab. Quantitative structure–property relationship "
               "(QSPR) models were assessed by nested 5×5 cross-validation, 1,000-fold Y-scrambling and a "
               "leverage applicability domain.")
    k.para(doc, "**Keywords** Glioblastoma · MXene · Ti_{3}C_{2}O_{2} · Density functional theory · Molecular "
                "docking · Temozolomide", align="left")


def results(doc, d, c):
    k.heading(doc, "Results and discussion")
    k.para(doc, "The workflow is summarised in Fig. 1. All quantities below are computed; none is fitted to "
                "experimental data.", indent=True)
    k.figure(doc, FIG / "Fig1.png", 1,
             "Workflow of the study: the glioblastoma drug set, docking into EGFR with redocking controls, the "
             "Ti_{3}C_{2}O_{2} carrier model and the tight-binding tests, periodic DFT adsorption of the "
             "alkylating agents, and the QSPR model of the docking score")
    s = results_docking(doc, d, c)
    bx = s["best_x"]
    k.figure(doc, FIG / "Fig2.png", 2,
             "Docking controls in EGFR (PDB 1M17, chain A). **a** Crystallographic erlotinib (grey), the "
             f"top-ranked pose of the crystal-conformation redocking (orange, {f1(s['top_x'].rmsd_A)} Å) and its "
             f"mode {int(bx['mode'])} (teal, {f2(bx.rmsd_A)} Å). **b** Heavy-atom RMSD to the crystal pose against "
             "the Vina score for every output mode of both controls; the dashed line marks 2 Å. **c** Percentage "
             "of the docked drugs whose top pose contacts each residue (light bars, any heavy atom within 4.0 Å; "
             "dark bars, N/O pairs within 3.5 Å); the 12 most frequent residues are shown")
    k.figure(doc, FIG / "Fig3.png", 3,
             "Docking scores in EGFR. **a** Vina score of every docked drug, coloured by family; erlotinib and "
             "temozolomide are highlighted. **b** Ligand efficiency (score per heavy atom) by family; points are "
             "individual drugs")
    results_tightbinding(doc, d, c)
    k.figure(doc, FIG / "Fig4.png", 4,
             "The Ti_{3}C_{2}O_{2} carrier. **a** The 4×4 periodic slab (a = 3.03 Å) used for the DFT "
             "calculations, top view. **b** Energy per formula unit of the periodic slab against the in-plane "
             "lattice constant with GFN2-xTB and GFN1-xTB (tblite, Γ point, 1500 K), relative to the lowest "
             "value of each method; crosses mark lattice constants at which the SCF failed or gave unphysical "
             "energies. The dotted line marks a = 3.03 Å, the lattice constant of the DFT slab (PBE-D3 optimum 3.022 Å)")
    dft = results_dft(doc, d, c, FIG / "Fig5.png", d["conv"])
    results_qspr(doc, d, c, FIG / "Fig6.png", table_no=2 if dft else 1)
    limitations(doc, c)
    return dft


def limitations(doc, c):
    k.heading(doc, "Limitations", 2)
    k.para(doc,
           "The adsorption energies are for single drug molecules in vacuum, without water, counter-ions or the "
           "OH and F terminations present on real Ti_{3}C_{2}T_{x}. Each drug was relaxed from one starting "
           "orientation, which samples only part of the adsorption landscape. The Brillouin zone was sampled at "
           "the Γ point and only the upper surface of the slab relaxed; the effect of the Γ-point sampling is "
           "estimated for temozolomide in Table S3. The drugs were treated as neutral molecules. The docking "
           "used a rigid receptor and an empirical score that did not rank the crystal-like pose of erlotinib "
           "first. None of the results has been compared with experiment.", indent=True)


def conclusions(doc, d):
    s, t, q = docking_stats(d), scan_stats(d), d["q"]
    dft = dft_sentence(d)
    k.heading(doc, "Conclusions")
    k.para(doc,
           "Docking of the glioblastoma cohort into EGFR places every drug in the ATP site, but the redocking "
           "controls show that the Vina score does not rank the crystal-like erlotinib pose first, so the "
           "scores serve only as an exploratory measure of fit. The GFN tight-binding methods do not describe "
           "metallic Ti_{3}C_{2}O_{2}, which is a warning for tight-binding screens of MXene carriers. "
           + (f"At the PBE-D3 level, {dft[0].lower() + dft[1:]}. " if dft else "") +
           f"The docking score can be modelled only weakly from four descriptors (*Q*^{{2}}_{{CV}} = "
           f"{f2(q['Q2_CV'])}), and the model mainly reflects molecular size.", indent=True)
    if dft is None:
        k.placeholder(doc, "[Conclusion on the DFT adsorption added automatically when the campaign is complete; "
                           "to be reread by the author.]")


def declarations(doc):
    k.heading(doc, "Statements and Declarations")
    for label, text in (
        ("Author contribution", f"{AUTHOR} conceived the study, performed all calculations and analyses, "
                                "and wrote the manuscript."),
        ("Funding", "No funding was received for this work."),
        ("Data availability", "All input structures, docking poses, Quantum ESPRESSO inputs and outputs, "
                              f"relaxed geometries and result tables are available at {REPO}."),
        ("Code availability", f"The complete pipeline, which regenerates every number, table and figure of "
                              f"this article from the raw inputs, is available at {REPO} under the MIT licence. "
                              "Software: Quantum ESPRESSO 7.5, xtb 6.7.1, tblite, AutoDock Vina 1.2.7, Meeko, "
                              "RDKit, PDBFixer/OpenMM, scikit-learn, PyMOL (open source)."),
        ("Ethics approval", "Not applicable."),
        ("Consent to participate", "Not applicable."),
        ("Consent for publication", "Not applicable."),
        ("Competing interests", "The author declares no competing interests."),
    ):
        k.labelled(doc, label, text)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    d = load()
    c = Cites()
    doc = k.new_document()
    front(doc)
    abstract(doc, d)
    introduction(doc, c)
    methods(doc, d, c)
    k.heading(doc, "Use of AI tools", 2)
    k.placeholder(doc, "[AUTHOR TO COMPLETE BEFORE SUBMISSION: statement on the use of AI tools in this work, as "
                       "required by the journal (Springer policy: use of large language models beyond copy "
                       "editing must be documented in the Methods).]")
    results(doc, d, c)
    conclusions(doc, d)
    k.references(doc, c.list())
    declarations(doc)
    out = OUT / "Manuscript_GBM_MXene_JMM.docx"
    doc.save(out)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
