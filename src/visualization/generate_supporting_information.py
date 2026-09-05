# -*- coding: utf-8 -*-
"""
generate_supporting_information.py  (GBM / MXene)
================================================
Honest Supporting Information: every table is built from a real result file.
The previous GBM_MXene_Supporting_Information.docx had no generator and its
Table S1 gave the SAME value to every compound (Vina -8.00, LE 0.280,
Delta_E_int -28.50); its "B3LYP-D3BJ DFT benchmark" and "measured
crystallographic contact distances" tables were fabricated (no ORCA / no
distance measurement was ever done). Those are removed here.
"""
import os
import numpy as np
import pandas as pd
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml import parse_xml
from docx.oxml.ns import nsdecls

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATASET = os.path.join(BASE, "data", "processed", "dataset_drug_mxene_pristine.csv")

TITLE = ("SUPPORTING INFORMATION\nQuantum-Chemical Modeling and Machine Learning-Driven "
         "Nano-QSAR of 2D Titanium Carbide MXene (Ti3C2Tx) Nanosheets for Targeted "
         "Glioblastoma Therapeutics Delivery")
S1 = ("The MXene nanocarrier is modeled as an oxygen-terminated trilayer Ti3C2O2 cluster. "
      "Drug-MXene complexes were built by placing each therapeutic parallel to the basal "
      "plane and relaxing with GFN2-xTB. The standardized single-point interaction energy is "
      "Delta_E_int,SP = E(complex) - E(MXene) - E(drug), both fragments taken at the "
      "complex geometry.")
REPO = "https://github.com/sircalch/mxene-glioblastoma-qsar-ai"
ZEN = "https://doi.org/10.5281/zenodo.22187857"
VINA_COLS = [("vina_4ZAU_kcal_mol", "Vina 4ZAU (kcal/mol)"), ("vina_2J6M_kcal_mol", "Vina 2J6M (kcal/mol)")]
ENDPOINT = ("GFN2-xTB single-point interaction energy Delta_E_int,SP (kcal/mol) of each drug "
            "on the pristine oxygen-terminated Ti3C2O2 MXene cluster.")
Q2_NOTE = ("Leak-free nested 5x5 CV on the real Delta_E_int,SP; see the manuscript for the "
           "per-fold Q2_CV. Exploratory feature-importance only.")

try:
    from rdkit import Chem
except Exception:
    Chem = None


def _bg(c, col):
    c._element.get_or_add_tcPr().append(parse_xml(f'<w:shd {nsdecls("w")} w:fill="{col}"/>'))


def _h(doc, text, level=1):
    x = doc.add_heading(text, level=level)
    x.paragraph_format.space_before = Pt(12)
    for r in x.runs:
        r.font.name = 'Times New Roman'; r.font.bold = True
        r.font.size = Pt(12.0 if level == 1 else 10.5)
        r.font.color.rgb = RGBColor(0, 77, 64)


def _table(doc, headers, rows):
    t = doc.add_table(rows=1, cols=len(headers))
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i, ht in enumerate(headers):
        c = t.rows[0].cells[i]; c.text = ht; _bg(c, "004D40")
        for r in c.paragraphs[0].runs:
            r.font.bold = True; r.font.color.rgb = RGBColor(255, 255, 255); r.font.size = Pt(7.5)
    for row in rows:
        cs = t.add_row().cells
        for i, v in enumerate(row):
            cs[i].text = str(v)
            for r in cs[i].paragraphs[0].runs:
                r.font.size = Pt(7.0)


def _formal_charge_from_smiles(smi):
    if Chem is None:
        return "n/a"
    m = Chem.MolFromSmiles(str(smi))
    return str(Chem.GetFormalCharge(m)) if m is not None else "n/a"


def _williams(df, feats, target):
    d = df.dropna(subset=feats + [target])
    X = d[feats].values
    n, p = X.shape
    Xd = np.hstack([np.ones((n, 1)), X])
    h = np.diag(Xd @ np.linalg.pinv(Xd.T @ Xd) @ Xd.T)
    hstar = 3.0 * (p + 1) / n
    y = d[target].values
    b = np.linalg.pinv(Xd.T @ Xd) @ Xd.T @ y
    res = y - Xd @ b
    sr = res / (np.std(res) * np.sqrt(np.maximum(1e-4, 1.0 - h)))
    return hstar, int(((h <= hstar) & (np.abs(sr) <= 3.0)).sum()), n


def generate_supporting_information():
    df = pd.read_csv(DATASET)
    doc = Document()
    for s in doc.sections:
        s.top_margin = s.bottom_margin = Inches(1.0)
        s.left_margin = s.right_margin = Inches(0.8)
    doc.styles['Normal'].font.name = 'Times New Roman'
    doc.styles['Normal'].font.size = Pt(10)

    r = doc.add_paragraph().add_run(TITLE)
    r.font.size = Pt(14); r.font.bold = True; r.font.color.rgb = RGBColor(0, 77, 64)
    doc.add_paragraph("Andrés Monreal Hernández, Sara Lizbeth Franco Amaya, and Carlos Ivanhoe Martínez Osorio").runs[0].font.italic = True

    _h(doc, "Section S1: Nanocarrier Cluster Model")
    doc.add_paragraph(S1)

    _h(doc, "Section S2: Software")
    doc.add_paragraph(
        "• AutoDock Vina v1.2.7 — molecular docking\n"
        "• GFN2-xTB (xtb v6.7.1, Grimme group) — geometry optimization and single-point energies\n"
        "• RDKit v2024.03.1 & Meeko v0.5.0 — descriptors, protonation, PDBQT preparation\n"
        "• scikit-learn v1.4.2, SciPy v1.13.0, SHAP v0.45.0 — surrogate model and interpretation\n"
        f"• Code: {REPO}\n• Data archive: {ZEN}")

    real_vcols = [(c, lab) for c, lab in VINA_COLS if c in df.columns]
    _h(doc, f"Table S1: Full Curated Dataset (N={len(df)}) — Real Docking Scores, Quantum "
            f"Descriptors and GFN2-xTB Single-Point Interaction Energies.")
    hdr = ["Compound", "Class", "MW (g/mol)"] + [lab for _, lab in real_vcols] + \
          ["E_HOMO (eV)", "omega (eV)", "Delta_E_int,SP (kcal/mol)"]
    rows = []
    for _, x in df.iterrows():
        row = [x["name"], str(x.get("drug_class", ""))[:26], f"{x['MolWt']:.1f}"]
        row += [f"{x[c]:.2f}" for c, _ in real_vcols]
        row += [f"{x['E_HOMO_eV']:.2f}", f"{x['Omega_eV']:.2f}", f"{x['delta_Eint_SP_kcal_mol']:.2f}"]
        rows.append(row)
    _table(doc, hdr, rows)

    _h(doc, "Table S2: Dominant Microstate Formal Charge at pH 7.4 (RDKit, from the curated "
            "SMILES; no external pKa engine was run).")
    s2 = []
    for _, x in df.iterrows():
        fc = x["formal_charge"] if "formal_charge" in df.columns and pd.notna(x.get("formal_charge")) \
            else _formal_charge_from_smiles(x.get("smiles", ""))
        s2.append([x["name"], str(x.get("drug_class", ""))[:26], int(fc) if str(fc).lstrip("-").isdigit() else fc])
    _table(doc, ["Compound", "Class", "Formal charge (pH 7.4)"], s2)

    _h(doc, "Table S3: OECD Principles 1–5 Checklist.")
    feats = [c for c in ["MolWt", "MolMR", "Omega_eV", "E_HOMO_eV"] if c in df.columns]
    hstar, inside, n = _williams(df, feats, "delta_Eint_SP_kcal_mol")
    s3 = [
        ("1. Defined endpoint", ENDPOINT),
        ("2. Unambiguous algorithm",
         "StandardScaler + regularized linear regression inside a leak-free nested 5x5 "
         "cross-validation."),
        ("3. Applicability domain",
         f"Williams hat-matrix leverage, {len(feats)} descriptors, n={n}: h* = {hstar:.3f}; "
         f"{inside}/{n} compounds inside the domain."),
        ("4. Goodness-of-fit / robustness", Q2_NOTE),
        ("5. Mechanistic interpretation",
         "Feature importance dominated by molecular size (MW), molar refractivity, "
         "electrophilicity (omega) and E_HOMO."),
    ]
    _table(doc, ["OECD principle", "Implementation"], s3)

    out = os.path.join(BASE, "manuscript", "GBM_MXene_Supporting_Information.docx")
    doc.save(out)
    print(f"[SUCCESS] Supporting Information: {out}")
    return out


if __name__ == "__main__":
    generate_supporting_information()
