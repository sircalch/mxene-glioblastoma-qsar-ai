"""
generate_gbm_word_manuscript.py
Builds the complete, publication-grade Microsoft Word (.docx) manuscript
with all 9 figures embedded, formatted tables, and 45 verified citations for Article 2.
"""

import os
import json
import pandas as pd
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT, WD_ALIGN_VERTICAL
from docx.oxml import OxmlElement, parse_xml
from docx.oxml.ns import nsdecls, qn

def set_cell_background(cell, fill_color):
    tcPr = cell._element.get_or_add_tcPr()
    shd = parse_xml(f'<w:shd {nsdecls("w")} w:fill="{fill_color}"/>')
    tcPr.append(shd)

def set_cell_margins(cell, top=100, bottom=100, left=150, right=150):
    tcPr = cell._element.get_or_add_tcPr()
    tcMar = OxmlElement('w:tcMar')
    for m, val in [('top', top), ('bottom', bottom), ('left', left), ('right', right)]:
        node = OxmlElement(f'w:{m}')
        node.set(qn('w:w'), str(val))
        node.set(qn('w:type'), 'dxa')
        tcMar.append(node)
    tcPr.append(tcMar)

def add_heading_styled(doc, text, level=1):
    h = doc.add_heading(text, level=level)
    h.paragraph_format.space_before = Pt(14)
    h.paragraph_format.space_after = Pt(6)
    h.paragraph_format.keep_with_next = True
    for r in h.runs:
        r.font.name = 'Times New Roman'
        r.font.bold = True
        if level == 1:
            r.font.size = Pt(14)
            r.font.color.rgb = RGBColor(21, 101, 192)
        elif level == 2:
            r.font.size = Pt(12)
            r.font.color.rgb = RGBColor(13, 71, 161)
        else:
            r.font.size = Pt(11)
            r.font.color.rgb = RGBColor(33, 33, 33)
    return h

def add_image_if_exists(doc, img_path, caption_text, width=Inches(6.2)):
    if os.path.exists(img_path):
        p_img = doc.add_paragraph()
        p_img.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p_img.paragraph_format.space_before = Pt(10)
        p_img.paragraph_format.space_after = Pt(4)
        run = p_img.add_run()
        run.add_picture(img_path, width=width)
        
        p_cap = doc.add_paragraph()
        p_cap.paragraph_format.space_after = Pt(12)
        p_cap.paragraph_format.line_spacing = 1.15
        r_num = p_cap.add_run(caption_text.split(':')[0] + ": ")
        r_num.font.bold = True
        r_num.font.size = Pt(9.5)
        r_num.font.color.rgb = RGBColor(21, 101, 192)
        
        r_desc = p_cap.add_run(':'.join(caption_text.split(':')[1:]))
        r_desc.font.size = Pt(9.5)
        r_desc.font.italic = True
    else:
        print(f"Warning: image {img_path} not found.")

def generate_gbm_word_manuscript():
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    fig_dir = os.path.join(base_dir, "figures")
    doc = Document()
    
    for s in doc.sections:
        s.top_margin = Inches(1.0)
        s.bottom_margin = Inches(1.0)
        s.left_margin = Inches(1.0)
        s.right_margin = Inches(1.0)
        
    style = doc.styles['Normal']
    font = style.font
    font.name = 'Times New Roman'
    font.size = Pt(11)
    font.color.rgb = RGBColor(33, 33, 33)
    
    # 1. Manuscript Header & Title
    p_title = doc.add_paragraph()
    p_title.paragraph_format.space_after = Pt(12)
    p_title.paragraph_format.line_spacing = 1.15
    r_title = p_title.add_run("Explainable AI and Quantum Chemical Exploration of 2D Titanium Carbide MXene (Ti3C2Tx) Nanosheets as Targeted Nanovehicles for Glioblastoma Therapeutics Across the Blood-Brain Barrier")
    r_title.font.size = Pt(16)
    r_title.font.bold = True
    r_title.font.color.rgb = RGBColor(21, 101, 192)
    
    p_auth = doc.add_paragraph()
    p_auth.paragraph_format.space_after = Pt(4)
    r_a1 = p_auth.add_run("Andrés Monreal Hernández")
    r_a1.font.bold = True
    p_auth.add_run("1,*, ")
    r_a2 = p_auth.add_run("Sara Lizbeth Franco Amaya")
    r_a2.font.bold = True
    p_auth.add_run("2, and ")
    r_a3 = p_auth.add_run("Carlos Ivanhoe Martínez Osorio")
    r_a3.font.bold = True
    p_auth.add_run("3")
    
    p_aff = doc.add_paragraph()
    p_aff.paragraph_format.space_after = Pt(14)
    p_aff.add_run(
        "1 Universidad Estatal de Sonora, Hermosillo, Sonora, Mexico. ORCID: 0009-0009-1207-8597\n"
        "2 Doctorado en Nanotecnología, Universidad de Sonora, Hermosillo, Sonora, Mexico. ORCID: 0009-0005-0272-0241\n"
        "3 Doctorado en Ciencia de Materiales, Universidad de Sonora, Hermosillo, Sonora, Mexico. ORCID: 0009-0003-7872-4965\n"
        "* Corresponding author: andres.monreal@ues.mx"
    )
    p_aff.runs[0].font.size = Pt(9.5)
    p_aff.runs[0].font.italic = True
    
    # Embedded Graphical Abstract
    add_image_if_exists(doc, os.path.join(fig_dir, "fig1_graphical_abstract.png"),
                        "Graphical Abstract: Atomistic, Quantum, and Machine Learning Framework for 2D Ti3C2Tx MXene Targeted Delivery Across the BBB in Glioblastoma.")
    
    # 2. Abstract & Keywords
    add_heading_styled(doc, "Abstract", level=1)
    p_abs = doc.add_paragraph()
    p_abs.paragraph_format.space_after = Pt(8)
    p_abs.paragraph_format.line_spacing = 1.15
    p_abs.add_run(
        "Glioblastoma multiforme (GBM) is the most lethal primary malignant central nervous system neoplasm in adults, with a median survival below 15 months, "
        "driven by therapeutic resistance and the restrictive physiology of the blood-brain barrier (BBB) [1,3]. Here we present an integrated computational "
        "framework combining GFN2-xTB tight-binding quantum chemistry (with D4 dispersion) [26,27], physical molecular docking (AutoDock Vina v1.2.7 [28,29] "
        "against the human EGFR kinase domain, PDB ID: 4ZAU, 2.80 Å, with 2J6M as a secondary control), and a leak-free cross-validated explainable Nano-QSAR "
        "surrogate, for a curated set of 35 clinical CNS and GBM therapeutics. Real GFN2-xTB single-point interaction energies of all 35 drugs on the pristine "
        "oxygen-terminated Ti3C2O2 MXene cluster range from -0.9 to -15.5 kcal/mol. An Angiopep-2-functionalized MXene for LRP-1-mediated transcytosis is "
        "discussed only as future work, since no real structural or quantum data for it exist in this study. Docking against EGFR (exploratory; redocking "
        "heavy-atom RMSD 5.32 Å) gave Vina scores of -3.96 to -8.94 kcal/mol (mean -7.16), with recurrent contacts at Asp392, His394, Arg427 and Thr391. "
        "A leak-free nested 5x5 cross-validated RidgeCV surrogate on the real data reached Q2_CV = 0.65 (isolated descriptors) and 0.10 (pristine-MXene "
        "interaction energy); the feature-importance analysis is reported as exploratory. OECD Principle 3 applicability-domain analysis (Williams leverage) "
        "places all 35 compounds inside the domain in both real-data systems. Every value reported is computed from the deposited pipeline; no descriptor "
        "or energy is estimated from an empirical formula."
    )
    
    p_kw = doc.add_paragraph()
    p_kw.paragraph_format.space_after = Pt(14)
    r_kwt = p_kw.add_run("Keywords: ")
    r_kwt.font.bold = True
    p_kw.add_run("2D MXene (Ti3C2Tx); Glioblastoma; Blood-Brain Barrier; EGFR kinase; AutoDock Vina; Explainable AI (SHAP); OECD Validation.")
    
    # 3. Section 1: Introduction
    add_heading_styled(doc, "1. Introduction", level=1)
    doc.add_paragraph(
        "Glioblastoma (GBM, WHO grade IV astrocytoma) is the most aggressive primary brain tumor in adults. "
        "Standard-of-care - maximal surgical resection, temozolomide chemotherapy and radiotherapy [1,2] - is almost invariably followed by recurrence, and "
        "MGMT-promoter methylation status modulates the temozolomide benefit [8]. "
        "A foundational barrier to clinical efficacy is the blood-brain barrier (BBB), which excludes the large majority of small-molecule therapeutics and "
        "nearly all large biologics from the brain parenchyma [9,10]."
    )
    doc.add_paragraph(
        "Overexpression and genomic amplification of the epidermal growth factor receptor (EGFR) and its constitutively active deletion mutant EGFRvIII "
        "occur in a large fraction of GBM patients and drive proliferation, invasion and neo-angiogenesis [4,6,7]. Systemic administration of EGFR tyrosine "
        "kinase inhibitors (TKIs) such as osimertinib, gefitinib and erlotinib is compromised by poor BBB penetration and rapid systemic clearance [10,13]."
    )
    doc.add_paragraph(
        "Two-dimensional transition-metal carbides and nitrides (MXenes), represented by titanium carbide (Ti3C2Tx) [21,22], combine metallic conductivity, "
        "hydrophilic surface terminations (-O, -OH, -F) and rich transition-metal coordination chemistry, and have been explored for nanomedicine [23-25]. "
        "Surface modification with Angiopep-2 peptides targeting low-density lipoprotein receptor-related protein 1 (LRP-1) has been proposed to promote "
        "transcytosis across the BBB; in the present study this functionalized carrier is treated only as a prospective extension (see Conclusions), and all "
        "quantum and docking results refer to the pristine Ti3C2O2 surface."
    )
    
    # Workflow Figure 1
    add_image_if_exists(doc, os.path.join(fig_dir, "fig1_gbm_workflow_methodology.png"),
                        "Figure 1: Multi-scale computational workflow: GFN2-xTB quantum-chemical adsorption on pristine Ti3C2O2, real AutoDock Vina docking (PDB 4ZAU / 2J6M), and a leak-free cross-validated explainable Nano-QSAR surrogate for 2D Ti3C2Tx MXene glioblastoma delivery.")
    
    # 4. Section 2: Computational and Experimental Methodology
    add_heading_styled(doc, "2. Results and Discussion", level=1)
    
    add_heading_styled(doc, "2.1 Quantum Adsorption Energetics & MXene Surface Chemistry", level=2)
    doc.add_paragraph(
        "Real GFN2-xTB single-point interaction energies (Delta_E_int,SP) across all 35 therapeutics on the pristine Ti3C2O2 cluster range from -0.9 kcal/mol "
        "for the most weakly interacting compounds to -15.5 kcal/mol for the most strongly stabilized (larotrectinib -15.5; temozolomide -13.4; lomustine "
        "-12.8 kcal/mol), consistent with dispersion-dominated physisorption of the drug pi-systems on the oxygen-terminated titanium-carbide surface."
    )

    add_image_if_exists(doc, os.path.join(fig_dir, "fig2_gbm_quantum_cdft_architecture.png"),
                        "Figure 2: Real Quantum CDFT Electronic Reactivity of the Isolated GBM Therapeutics (real GFN2-xTB single points, n=35): (a) Frontier Molecular Orbital (HOMO/LUMO) distribution; (b) Chemical hardness vs. electrophilicity index. No real complex-level frontier-orbital calculation exists for either MXene variant.")

    add_heading_styled(doc, "2.2 Molecular docking against the EGFR kinase domain", level=2)
    doc.add_paragraph(
        "AutoDock Vina v1.2.7 screening of the 35 therapeutics against EGFR (PDB 4ZAU) gave scores from -3.96 to -8.94 kcal/mol (mean -7.16 kcal/mol). "
        "The highest-ranked compounds were entrectinib (-8.94), sorafenib (-8.66), trametinib (-8.50), cobimetinib (-8.49) and cabozantinib (-8.29 kcal/mol). "
        "Because self-redocking reproduced the native pose only within 5.32 Å RMSD, this ranking is treated as exploratory and is not used as a QSAR endpoint."
    )

    # Docking Figures 3 and 4
    add_image_if_exists(doc, os.path.join(fig_dir, "fig3_gbm_docking_vina_statistical_profiles.png"),
                        "Figure 3: Molecular docking statistical profiles on the EGFR kinase domain (PDB 4ZAU; exploratory, redocking RMSD 5.32 Å): (a) distribution of real Vina scores; (b) ranking of the top-10 compounds.")

    add_image_if_exists(doc, os.path.join(fig_dir, "fig4_gbm_residue_contact_frequency.png"),
                        "Figure 4: Residue-level contact frequencies on the EGFR kinase domain (real Vina poses, contact distance <= 3.8 Å): most frequent contacts are Asp392, His394, Arg427, Thr391 and Arg390.")
    
    # Embed Table 1: Descriptors Summary
    # MW/LogP/PSA are real RDKit descriptors (always computed from SMILES).
    # E_HOMO/omega previously came from gbm_isolated_descriptors.csv, whose
    # E_HOMO was an empirical-formula placeholder ("-5.10 - 0.22*LogP - ...")
    # never overwritten with real data; merged here with real GFN2-xTB
    # frontier orbitals parsed from calculations/gbm/*/*_drug_sp.out.
    desc_csv = os.path.join(base_dir, "data", "processed", "gbm_isolated_descriptors.csv")
    homo_lumo_csv = os.path.join(base_dir, "data", "processed", "gbm_isolated_real_homo_lumo.csv")
    if os.path.exists(desc_csv):
        df_desc = pd.read_csv(desc_csv)
        if os.path.exists(homo_lumo_csv):
            df_real = pd.read_csv(homo_lumo_csv)
            df_desc = df_desc.merge(df_real, on="name", how="inner")
            df_desc["E_HOMO"] = df_desc["E_HOMO_real_eV"]
            gap = df_desc["E_LUMO_real_eV"] - df_desc["E_HOMO_real_eV"]
            mu = -(df_desc["E_HOMO_real_eV"] + df_desc["E_LUMO_real_eV"]) / 2.0
            df_desc["Electrophilicity_omega"] = mu ** 2 / (2.0 * (gap / 2.0))
        doc.add_paragraph()
        p_t1 = doc.add_paragraph()
        r_t1 = p_t1.add_run("Table 1: Physicochemical, Topological, and Quantum CDFT Descriptors for Representative GBM Therapeutics.")
        r_t1.font.bold = True
        r_t1.font.size = Pt(10)
        
        table1 = doc.add_table(rows=1, cols=7)
        table1.alignment = WD_TABLE_ALIGNMENT.CENTER
        hdr_cells = table1.rows[0].cells
        hdr_titles = ["Compound", "Class", "MW (g/mol)", "LogP", "PSA (Å²)", "E_HOMO (eV)", "omega (eV)"]
        for idx, title in enumerate(hdr_titles):
            hdr_cells[idx].text = title
            set_cell_background(hdr_cells[idx], "1565C0")
            set_cell_margins(hdr_cells[idx], 80, 80, 100, 100)
            for r in hdr_cells[idx].paragraphs[0].runs:
                r.font.bold = True
                r.font.color.rgb = RGBColor(255, 255, 255)
                r.font.size = Pt(9)
                
        for _, row in df_desc.head(10).iterrows():
            row_cells = table1.add_row().cells
            row_vals = [
                str(row['name']), str(row['drug_class'])[:22], f"{row['MW']:.1f}",
                f"{row['LogP']:.2f}", f"{row['PSA']:.1f}", f"{row['E_HOMO']:.2f}", f"{row['Electrophilicity_omega']:.2f}"
            ]
            for c_idx, val in enumerate(row_vals):
                row_cells[c_idx].text = val
                set_cell_margins(row_cells[c_idx], 60, 60, 80, 80)
                for r in row_cells[c_idx].paragraphs[0].runs:
                    r.font.size = Pt(8.5)
                    
    add_heading_styled(doc, "2.3 Nano-QSAR surrogate model and feature importance", level=2)
    doc.add_paragraph(
        "A StandardScaler + RidgeCV surrogate evaluated by leak-free nested 5x5 cross-validation on the real observed data (isolated Vina scores; and the real "
        "GFN2-xTB Delta_E_int,SP on the pristine Ti3C2O2 MXene) reached Q2_CV = 0.65 for the isolated-descriptor model and 0.10 for the pristine-MXene "
        "interaction-energy model, with RMSE of 0.54 and 3.95 kcal/mol respectively (n = 35, four descriptors each: MolWt, MolMR, E_HOMO, omega). The "
        "pristine-MXene model is therefore only weakly predictive, and the exploratory ExtraTrees feature-importance ranking (Figure 6) - led by molecular "
        "weight and molar refractivity - is reported as a qualitative indication rather than a validated structure-property relationship [39,40]."
    )

    # ML Parity and SHAP Figures 5 and 6
    add_image_if_exists(doc, os.path.join(fig_dir, "fig5_gbm_parity_models_evaluation.png"),
                        "Figure 5: Leak-free nested 5x5 CV parity plots (real observed vs out-of-fold predicted) for Isolated and Pristine-MXene systems. No real structural/quantum data exists for the functionalized Ti3C2-Angiopep-2 system, so it is not shown.")

    add_image_if_exists(doc, os.path.join(fig_dir, "fig6_gbm_shap_xai_importance_rankings.png"),
                        "Figure 6: Exploratory Feature Importance Rankings on the real GFN2-xTB pristine-MXene interaction energy, identifying molecular weight and molar refractivity as the leading descriptors.")
    
    # Inter-descriptor Correlation Figure 7
    add_image_if_exists(doc, os.path.join(fig_dir, "fig7_gbm_descriptor_correlation_matrix.png"),
                        "Figure 7: Pearson inter-descriptor correlation heatmap (real descriptor matrix, 35 GBM therapeutics).")
    
    add_heading_styled(doc, "2.4 Applicability domain (OECD Principle 3)", level=2)
    doc.add_paragraph(
        "The applicability domain was assessed by hat-matrix leverage on the real descriptor matrix (Williams plot) [31-33]. With the full descriptor set the "
        "warning leverage is h* = 1.80 for the isolated-drug system (20 descriptors) and h* = 0.77 for the pristine-MXene system (8 descriptors); all 35 "
        "compounds fall inside the domain (leverage below h* and standardized residual within +/-3sigma) in both cases."
    )

    # Williams Domain Figure 8
    add_image_if_exists(doc, os.path.join(fig_dir, "fig8_gbm_williams_applicability_domain.png"),
                        "Figure 8: OECD Principle 3: Williams Plots Defining the Applicability Domain for GBM Therapeutics on Ti3C2Tx MXene Nanosheets (real data only; Isolated and Pristine-MXene systems).")
    
    add_heading_styled(doc, "2.5 Representative binding modes", level=2)
    doc.add_paragraph(
        "Inspection of the top-ranked docked poses on EGFR (PDB 4ZAU) shows the inhibitors occupying the ATP cleft with contacts to the residues in Figure 4 "
        "(Asp392, His394, Arg427, Thr391). On the pristine Ti3C2O2 surface the strongest-binding drugs lie flat against the oxygen termination, consistent with "
        "the dispersion-dominated interaction energies of Section 2.1."
    )

    # 3D Spatial Figure 9
    add_image_if_exists(doc, os.path.join(fig_dir, "fig9_gbm_3d_spatial_binding_modes.png"),
                        "Figure 9: Representative binding modes (schematic): (a-b) top-ranked inhibitors in the EGFR ATP cleft (PDB 4ZAU); (c) a drug on the pristine Ti3C2O2 MXene surface with its real GFN2-xTB Delta_E_int,SP.")

    add_heading_styled(doc, "2.6 Interfacial charge redistribution", level=2)
    doc.add_paragraph(
        "The charge-density difference Delta_rho = rho(complex) - rho(carrier) - rho(drug) was computed from the real GFN2-xTB "
        "densities of the strongest-adsorbing complex (Larotrectinib / Ti3C2O2), all fragments at the bound geometry on a common grid "
        "(Figure 10). Accumulation (yellow) and depletion (blue) lobes concentrate at the drug oxygen / ester functionality facing "
        "the oxygen-terminated MXene, where the drug-surface separation is shortest, while the rest of the drug is essentially "
        "unperturbed. This localised reorganisation is consistent with the moderate interaction energy of this complex "
        "(Delta_E_int,SP = -15.5 kcal/mol), the largest in the set. The Delta_rho cube and its build script are in results/quantum/drho/."
    )
    add_image_if_exists(doc, os.path.join(fig_dir, "fig10_gbm_charge_density_difference.png"),
                        "Figure 10: Charge-density difference (real GFN2-xTB densities) for the Larotrectinib / Ti3C2O2 MXene complex. "
                        "Isovalue +/-0.006 e bohr^-3; yellow = electron accumulation, blue = electron depletion. "
                        "Delta_E_int,SP = -15.5 kcal/mol.")

    # 6. Section 3: Conclusions
    add_heading_styled(doc, "3. Conclusions", level=1)
    doc.add_paragraph(
        "We report a quantum-informed, explainable Nano-QSAR analysis of 2D Ti3C2Tx MXene as a candidate delivery scaffold for glioblastoma therapeutics. "
        "Combining real GFN2-xTB single-point interaction energies, exploratory AutoDock Vina docking against the EGFR kinase domain (PDB 4ZAU), and a "
        "leak-free cross-validated RidgeCV surrogate assessed under OECD guidelines, we find that the pristine Ti3C2O2 surface binds all 35 screened "
        "therapeutics by dispersion-dominated physisorption (Delta_E_int,SP = -0.9 to -15.5 kcal/mol). The surrogate is only weakly predictive for the "
        "MXene interaction energy, so the descriptor rankings are presented as exploratory. An Angiopep-2-functionalized MXene for LRP-1-mediated BBB "
        "transcytosis is a natural extension but is outside the present scope: it has no real structural or quantum data here and would require dedicated "
        "complex-geometry modeling."
    )
    
    # 7. Statements & References
    add_heading_styled(doc, "4. Experimental", level=1)
    doc.add_paragraph(
        "4.1 Quantum-chemical framework: "
        "Each isolated drug, the pristine oxygen-terminated Ti3C2O2 MXene cluster, and every drug-MXene complex were geometry-optimized and evaluated at single "
        "point with GFN2-xTB (xtb v6.7.1) including the D4 dispersion correction [26,27]. The standardized single-point interaction energy is "
        "Delta_E_int,SP = E(complex) - E(MXene) - E(drug), with both fragments taken at the complex geometry. Frontier-orbital energies (E_HOMO, E_LUMO) and "
        "conceptual-DFT global reactivity indices - chemical hardness (eta = gap/2), softness (S), electronegativity (chi) and electrophilicity index "
        "(omega = mu^2/2eta) [36,37] - were read directly from the xtb output; no descriptor is estimated from an empirical formula."
    )
    doc.add_paragraph(
        "4.2 Molecular docking: "
        "The X-ray crystal structure of the human EGFR kinase domain (PDB ID: 4ZAU, 2.80 Å) was used as the primary receptor, with 2J6M as a secondary control. "
        "Receptor and ligands were protonated at pH 7.4 and converted to PDBQT with Meeko; rigid-receptor flexible-ligand docking used AutoDock Vina v1.2.7 [28,29] "
        "over a 22 x 22 x 22 Å grid centred on the ATP pocket. Because self-redocking of the co-crystallized ligand reproduced the native pose only within "
        "5.32 Å heavy-atom RMSD, the docking scores are reported as an exploratory ranking rather than a quantitative endpoint."
    )
    doc.add_paragraph(
        "4.3 Surrogate model and applicability domain: "
        "A regularized linear model (StandardScaler + RidgeCV) was trained inside a leak-free nested 5x5 cross-validation on the real observed data. Feature "
        "importance was inspected with an ExtraTrees estimator and Shapley Additive Explanations (SHAP) [35] and is reported as exploratory only. The "
        "applicability domain follows OECD Principle 3 [31-33] via hat-matrix leverage (Williams plot) with warning leverage h* = 3(p+1)/n and +/-3sigma "
        "standardized-residual limits."
    )
    
    # Quantum Figure 2
    
    # 5. Section 3: Results and Discussion

    add_heading_styled(doc, "Data Availability", level=1)
    doc.add_paragraph(
        "All code, the curated dataset, the real GFN2-xTB and AutoDock Vina outputs, the leak-free cross-validation predictions and the figure/manuscript "
        "generators are in the public repository https://github.com/sircalch/mxene-glioblastoma-qsar-ai. run_entire_gbm_study.py reproduces every value and figure."
    )

    import _backmatter
    add_heading_styled(doc, "Conflict of Interest", level=1)
    doc.add_paragraph("The authors declare no competing financial or non-financial interest.")
    _backmatter.append(doc, add_heading_styled,
                       "GBM_MXene_Supporting_Information.docx",
                       "curated dataset (N = 35), formal charges, OECD checklist and per-residue contact frequencies",
                       already_has_conflict=True)

    add_heading_styled(doc, "References", level=1)
    from build_gbm_verified_references import GBM_VERIFIED_REFERENCES as VERIFIED_REFERENCES
    for idx, ref in enumerate(VERIFIED_REFERENCES, 1):
        p_ref = doc.add_paragraph()
        p_ref.paragraph_format.left_indent = Inches(0.4)
        p_ref.paragraph_format.space_after = Pt(3)
        r_num = p_ref.add_run(f"{idx}. ")
        r_num.font.bold = True
        p_ref.add_run(ref['citation'] + " ")
        if ref.get('doi'):
            r_doi = p_ref.add_run(f"doi:{ref['doi']}")
            r_doi.font.italic = True
            r_doi.font.size = Pt(9.0)
            r_doi.font.color.rgb = RGBColor(21, 101, 192)
        
    out_docx = os.path.join(base_dir, "manuscript", "Beilstein_Manuscript_GBM_MXene_Monreal_Hernandez_et_al.docx")
    doc.save(out_docx)
    print(f"Generated Comprehensive GBM Word Manuscript: {out_docx}")
    return out_docx

if __name__ == "__main__":
    generate_gbm_word_manuscript()
