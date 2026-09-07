"""Shared Beilstein back-matter block for the four nano-QSAR manuscripts.

append(doc, add_heading, si_file, si_contents, funding=None)
    inserts, in Beilstein order, the sections that go between the main text and
    References: Supporting Information, Funding, Acknowledgements, Author
    Contributions, Conflict of Interest, ORCID iDs.
"""
try:
    from docx.shared import Inches, Pt
except Exception:
    Inches = Pt = lambda x: x

_AUTHORS = [
    ("Andrés Monreal Hernández", "0009-0009-1207-8597"),
    ("Sara Lizbeth Franco Amaya", None),
    ("Carlos Ivanhoe Martínez Osorio", None),
]
_DEFAULT_FUNDING = ("This work was supported by Universidad Estatal de Sonora and the "
                    "Doctorado en Nanotecnología, Universidad de Sonora. No external grant "
                    "funding was received.")


def append(doc, add_heading, si_file, si_contents, funding=None,
           already_has_conflict=False, already_has_ack=False):
    add_heading(doc, "Supporting Information", level=1)
    p = doc.add_paragraph()
    p.add_run("Supporting Information File 1: ").bold = True
    p.add_run(f"{si_file} — {si_contents}")
    doc.add_paragraph(
        "The complete input/output archive (curated dataset, GFN2-xTB and docking "
        "outputs, cross-validation predictions, figure-generation and manuscript "
        "scripts) is in the public repository and Zenodo deposit cited under Data "
        "Availability; every value in this article is reproduced by "
        "run_entire_*_study.py.")

    if not already_has_ack:
        add_heading(doc, "Funding", level=1)
        doc.add_paragraph(funding or _DEFAULT_FUNDING)
        add_heading(doc, "Acknowledgements", level=1)
        doc.add_paragraph(
            "The authors thank the computational resources of Universidad Estatal "
            "de Sonora and Universidad de Sonora.")
    else:
        add_heading(doc, "Funding", level=1)
        doc.add_paragraph(funding or _DEFAULT_FUNDING)

    add_heading(doc, "Author Contributions", level=1)
    doc.add_paragraph(
        "Andrés Monreal Hernández: conceptualization, methodology, software, "
        "formal analysis, investigation, data curation, visualization, "
        "writing – original draft. Sara Lizbeth Franco Amaya: validation, "
        "writing – review and editing. Carlos Ivanhoe Martínez Osorio: "
        "supervision, writing – review and editing. All authors read and "
        "approved the final manuscript.")

    if not already_has_conflict:
        add_heading(doc, "Conflict of Interest", level=1)
        doc.add_paragraph("The authors declare no competing financial or non-financial interest.")

    add_heading(doc, "ORCID iDs", level=1)
    for name, orcid in _AUTHORS:
        if orcid:
            doc.add_paragraph(f"{name} – https://orcid.org/{orcid}")
