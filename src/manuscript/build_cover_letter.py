"""build_cover_letter.py - cover letter for Journal of Molecular Modeling (Word). Numbers
from the pipeline outputs; the DFT sentence is written only when the campaign is complete."""
import sys
from datetime import date
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import docx_kit as k  # noqa: E402
from build_manuscript import AFFIL, AUTHOR, EMAIL, ORCID, TITLE, dft_sentence, f1, f2, load  # noqa: E402
from sec_results_a import docking_stats, scan_stats  # noqa: E402

OUT = HERE.parents[1] / "manuscript" / "submission"
ZENODO_V2 = "10.5281/zenodo.22700227"        # supporting files of the pre-rebuild version (2026-09-11)


def main():
    d = load()
    s, t, q = docking_stats(d), scan_stats(d), d["q"]
    dft = dft_sentence(d)
    doc = k.new_document()
    for x in (AUTHOR, AFFIL, EMAIL, "", date.today().strftime("%d %B %Y"), "",
              "The Editor-in-Chief", "Journal of Molecular Modeling", ""):
        k.para(doc, x, align="left", space_after=0)
    k.para(doc, "Dear Editor,", align="left")
    k.para(doc, f"I submit the manuscript \"{TITLE}\" for consideration as an Original Paper in the Journal of "
                "Molecular Modeling.")
    k.para(doc,
           f"The study examines {len(d['m'])} glioblastoma drugs against the erlotinib site of EGFR and against "
           "an O-terminated Ti_{3}C_{2}O_{2} MXene proposed as a carrier. Two results concern methods that are "
           "widely used in screens of this kind. First, every docking mode of two redocking controls was "
           "compared with the crystal pose: the crystal-like erlotinib pose is found "
           f"(≤ {f2(max(s['best_x'].rmsd_A, s['best_s'].rmsd_A))} Å) but is not ranked first "
           f"({f1(s['top_x'].rmsd_A)} and {f1(s['top_s'].rmsd_A)} Å), so the docking scores are reported as "
           "exploratory. Second, the GFN2-xTB and GFN1-xTB tight-binding methods fail for this metallic carrier: "
           f"GFN2-xTB places the lattice minimum at {f2(t['g2_min'])} Å instead of 3.03 Å, GFN1-xTB has no "
           f"minimum, and only {t['n_neutral_ok']} of {t['n_neutral']} attempts on neutral finite flakes "
           "converged. Adsorption was therefore computed with periodic density functional theory (PBE-D3). "
           + (f"{dft}. " if dft else "[DFT RESULT SENTENCE ADDED AUTOMATICALLY WHEN THE CAMPAIGN IS COMPLETE.] ") +
           "A descriptor model of the docking score, validated by nested cross-validation and Y-scrambling, is "
           f"weak (*Q*^{{2}}_{{CV}} = {f2(q['Q2_CV'])}) and mainly reflects molecular size, and it is reported as such.")
    k.para(doc,
           "Related work. Supporting files of an earlier version of this study were deposited on Zenodo "
           f"(https://doi.org/{ZENODO_V2}); that version was superseded when the study was rebuilt from its raw "
           "inputs, and the results reported here replace it. A methods paper by the author, in preparation for "
           "the Journal of Chemical Information and Modeling, uses this study as one of four case studies of "
           "errors found and corrected during such rebuilds, and cites some of its summary numbers; the study, "
           "including its docking and DFT adsorption results, is reported in full only in this manuscript.")
    k.para(doc,
           "All structures, docking poses, Quantum ESPRESSO inputs and outputs and the complete pipeline that "
           "regenerates every number and figure are openly available. The manuscript has not been published and "
           "is not under consideration elsewhere. The author declares no competing interests.")
    k.para(doc, "Sincerely,", align="left", space_after=0)
    k.para(doc, f"{AUTHOR} (ORCID {ORCID})", align="left")
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / "Cover_Letter_JMM.docx"
    doc.save(out)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
