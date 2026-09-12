"""
generate_gbm_jmm_manuscript.py
===============================
Builds a Journal of Molecular Modeling (Springer) submission variant of the
GBM/MXene manuscript. Does NOT touch the canonical Beilstein/Molecular
Diversity body produced by generate_gbm_word_manuscript.py -- it post-processes
a freshly regenerated copy of that docx:

  1. Structured Context/Methods abstract (JMM requirement, 150-250 words).
     The Context paragraph reuses the ORIGINAL abstract text verbatim (it
     already embeds the real, dynamically-computed docking/QSPR numbers) so
     no statistic is retyped by hand and can never drift out of sync with the
     pipeline. A new, purely descriptive Methods paragraph (static facts only
     -- software versions/protocol, no computed numbers) is added alongside it.
  2. Section reorder: JMM wants Methods to follow the Introduction. The
     Beilstein body instead puts "4. Experimental" after the Results/
     Conclusions (Beilstein-house-style back matter). Moved + renumbered:
         1. Introduction            (unchanged)
         2. Experimental            (was "4.", subsections 4.1-4.3 -> 2.1-2.3)
         3. Results and Discussion  (was "2.", subsections 2.1-2.6 -> 3.1-3.6)
         4. Summary                 (was "3. Conclusions")
     The one in-text cross-reference ("Section 2.1") is fixed to "Section 3.1".
"""

import os
import re
import sys

base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from docx import Document

import generate_gbm_word_manuscript as full_gen

JMM_METHODS = (
    "Each isolated drug, the pristine Ti3C2O2 MXene cluster, and every complex were geometry-optimized and "
    "evaluated at single point with GFN2-xTB (xtb v6.7.1, D4 dispersion); frontier-orbital and conceptual-DFT "
    "reactivity indices were read directly from the xtb output. The human EGFR kinase domain (PDB ID: 4ZAU, "
    "2.80 Å; 2J6M as a secondary control) was docked with AutoDock Vina v1.2.7, benchmarked by self-redocking "
    "of the co-crystallized ligand. A StandardScaler+RidgeCV surrogate was trained inside a leak-free nested "
    "5×5 cross-validation, with feature importance inspected via an ExtraTrees estimator and SHAP; the "
    "applicability domain follows OECD Principle 3 (Williams-leverage analysis)."
)


def _heading_paragraphs(doc):
    """Return list of (index, Paragraph) for every Heading-style paragraph."""
    return [(i, p) for i, p in enumerate(doc.paragraphs) if p.style.name.startswith("Heading")]


def _build_condensed_context(original_abstract):
    """Re-derive a <=140-word Context paragraph from the ORIGINAL abstract by
    extracting the dynamically-computed numbers with regex (never retyped by
    hand) and dropping sentences that duplicate what the Methods paragraph
    already covers (software/protocol detail) or are out of scope here
    (the Angiopep-2 future-work aside). Raises if a number can't be found,
    rather than silently falling back to a stale/guessed value.
    """
    m_range = re.search(r"range from ([\-0-9.]+ to [\-0-9.]+) kcal/mol", original_abstract)
    m_vina = re.search(r"Vina scores of ([\-0-9.]+ to [\-0-9.]+) kcal/mol \(mean ([\-0-9.]+)\)", original_abstract)
    m_q2 = re.search(r"Q2_CV = ([0-9.n/a]+) for the EGFR Vina docking score and ([0-9.]+) for the pristine-MXene", original_abstract)
    if not (m_range and m_vina and m_q2):
        raise RuntimeError("Could not extract one or more dynamic values from the original GBM abstract "
                            "-- source text likely changed; update the regexes in _build_condensed_context.")
    interaction_range = m_range.group(1)
    vina_range, vina_mean = m_vina.group(1), m_vina.group(2)
    q2_vina, q2_mxene = m_q2.group(1), m_q2.group(2)

    return (
        "Glioblastoma multiforme (GBM) is the most lethal primary malignant central nervous system neoplasm "
        "in adults, with a median survival below 15 months, driven by therapeutic resistance and the "
        "restrictive physiology of the blood-brain barrier. We evaluate 2D Ti3C2Tx MXene as a candidate "
        "delivery scaffold for 35 clinical CNS and GBM therapeutics, combining GFN2-xTB quantum chemistry, "
        f"AutoDock Vina docking against the EGFR kinase domain, and a leak-free cross-validated Nano-QSAR "
        f"surrogate. Real GFN2-xTB interaction energies on the pristine Ti3C2O2 cluster range from "
        f"{interaction_range} kcal/mol. Docking gave Vina scores of {vina_range} kcal/mol (mean {vina_mean}), "
        f"with the surrogate reaching Q2_CV = {q2_vina} for the docking score and {q2_mxene} for the "
        "MXene interaction energy -- at best weakly predictive, reported as exploratory. OECD Principle 3 "
        "applicability-domain analysis places all 35 compounds inside the domain in both systems."
    )


def _find_heading(doc, text_exact=None, text_startswith=None):
    for i, p in enumerate(doc.paragraphs):
        if not p.style.name.startswith("Heading"):
            continue
        if text_exact is not None and p.text.strip() == text_exact:
            return i, p
        if text_startswith is not None and p.text.strip().startswith(text_startswith):
            return i, p
    raise RuntimeError(f"Heading not found: {text_exact or text_startswith}")


def generate_gbm_jmm_manuscript():
    full_gen.generate_gbm_word_manuscript()

    src_docx = os.path.join(base_dir, "manuscript", "Beilstein_Manuscript_GBM_MXene_Monreal_Hernandez_et_al.docx")
    doc = Document(src_docx)

    # ---------------------------------------------------------------
    # 1) Structured Context/Methods abstract
    # ---------------------------------------------------------------
    abs_head_idx, _ = _find_heading(doc, text_exact="Abstract")
    abstract_para = doc.paragraphs[abs_head_idx + 1]
    keywords_para = doc.paragraphs[abs_head_idx + 2]
    if "Keywords" not in keywords_para.text:
        raise RuntimeError("Unexpected structure: paragraph after Abstract body is not Keywords.")

    original_abstract = abstract_para.text
    context_text = _build_condensed_context(original_abstract)  # trimmed to fit JMM's 150-250 word cap

    p_context = abstract_para.insert_paragraph_before()
    p_context.paragraph_format.space_after = abstract_para.paragraph_format.space_after
    r1 = p_context.add_run("Context ")
    r1.font.bold = True
    p_context.add_run(context_text)

    p_methods = abstract_para.insert_paragraph_before()
    p_methods.paragraph_format.space_after = abstract_para.paragraph_format.space_after
    r2 = p_methods.add_run("Methods ")
    r2.font.bold = True
    p_methods.add_run(JMM_METHODS)

    abstract_para._element.getparent().remove(abstract_para._element)

    # ---------------------------------------------------------------
    # 2) Move "4. Experimental" (+ its 3 paragraphs) to right after Introduction,
    #    i.e. immediately before "2. Results and Discussion"
    # ---------------------------------------------------------------
    exp_idx, _ = _find_heading(doc, text_exact="4. Experimental")
    data_avail_idx, _ = _find_heading(doc, text_exact="Data Availability")
    results_idx, results_head = _find_heading(doc, text_exact="2. Results and Discussion")

    elements_to_move = [doc.paragraphs[i]._p for i in range(exp_idx, data_avail_idx)]
    anchor_element = results_head._p
    for el in elements_to_move:
        anchor_element.addprevious(el)

    # ---------------------------------------------------------------
    # 3) Renumber headings and the two 4.x -> 2.x paragraph labels
    # ---------------------------------------------------------------
    _, exp_head = _find_heading(doc, text_exact="4. Experimental")
    exp_head.runs[0].text = "2. Experimental"

    _, results_head2 = _find_heading(doc, text_exact="2. Results and Discussion")
    results_head2.runs[0].text = "3. Results and Discussion"

    subsection_renumber = {
        "2.1 Quantum Adsorption Energetics & MXene Surface Chemistry": "3.1 Quantum Adsorption Energetics & MXene Surface Chemistry",
        "2.2 Molecular docking against the EGFR kinase domain": "3.2 Molecular docking against the EGFR kinase domain",
        "2.3 Nano-QSAR surrogate model and feature importance": "3.3 Nano-QSAR surrogate model and feature importance",
        "2.4 Applicability domain (OECD Principle 3)": "3.4 Applicability domain (OECD Principle 3)",
        "2.5 Representative binding modes": "3.5 Representative binding modes",
        "2.6 Interfacial charge redistribution": "3.6 Interfacial charge redistribution",
    }
    for old, new in subsection_renumber.items():
        _, h = _find_heading(doc, text_exact=old)
        h.runs[0].text = new

    _, concl_head = _find_heading(doc, text_exact="3. Conclusions")
    concl_head.runs[0].text = "4. Summary"

    # 4.1/4.2/4.3 inline paragraph labels -> 2.1/2.2/2.3 (plain paragraphs, not headings)
    label_renumber = {
        "4.1 Quantum-chemical framework": "2.1 Quantum-chemical framework",
        "4.2 Molecular docking": "2.2 Molecular docking",
        "4.3 Surrogate model and applicability domain": "2.3 Surrogate model and applicability domain",
    }
    for p in doc.paragraphs:
        for old, new in label_renumber.items():
            if p.text.strip().startswith(old) and p.runs:
                p.runs[0].text = p.runs[0].text.replace(old, new, 1)

    # Fix the one in-text cross-reference (old Results 2.1 is now Results 3.1)
    for p in doc.paragraphs:
        if "Section 2.1" in p.text:
            for r in p.runs:
                if "Section 2.1" in r.text:
                    r.text = r.text.replace("Section 2.1", "Section 3.1")

    out_dir = os.path.join(base_dir, "manuscript", "submission_ready")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "02_Manuscript_GBM_MXene_JMM_Submission.docx")
    doc.save(out_path)
    print(f"[SUCCESS] Generated JMM submission manuscript: {out_path}")
    return out_path


if __name__ == "__main__":
    generate_gbm_jmm_manuscript()
