# Explainable AI and Quantum Chemical Exploration of 2D Titanium Carbide MXene (Ti3C2Tx) Nanocarriers for Glioblastoma Therapeutics

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22187857.svg)](https://doi.org/10.5281/zenodo.22187857)
[![Version](https://img.shields.io/badge/version-2.0.0-blue.svg)](https://github.com/sircalch/mxene-glioblastoma-qsar-ai/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-teal.svg)](https://opensource.org/licenses/MIT)
[![Python 3.12](https://img.shields.io/badge/Python-3.12-blue.svg)](https://www.python.org/)
[![AutoDock Vina](https://img.shields.io/badge/Docking-AutoDock%20Vina-orange.svg)](https://github.com/ccsb-scripps/AutoDock-Vina)
[![XAI: SHAP](https://img.shields.io/badge/Explainability-SHAP-purple.svg)](https://github.com/shap/shap)

**Authors**: Andrés Monreal Hernández, Sara Lizbeth Franco Amaya, Carlos Ivanhoe Martínez Osorio  
**Affiliation**: Universidad Estatal de Sonora, Hermosillo, Sonora, México  

---

## 📌 Abstract

Glioblastoma multiforme (GBM) remains the most lethal primary central nervous system malignancy, characterized by diffuse cerebral infiltration, aggressive therapeutic resistance, and severe restriction by the Blood-Brain Barrier (BBB). In this work, we present an integrated **atomistic quantum chemical, molecular docking, and explainable machine learning (Nano-QSAR / XAI)** investigation of **two-dimensional (2D) titanium carbide MXene ($\text{Ti}_3\text{C}_2\text{T}_x$) nanosheets** as multifunctional nanocarriers for targeted delivery of 37 clinical-stage and FDA-approved glioblastoma therapeutics.

### Methodological Highlights:
- **Oncogenic Kinase Target**: Molecular docking against human Epidermal Growth Factor Receptor (EGFR kinase domain, PDB ID: 4UV7).
- **Curated GBM Therapeutic Library**: 37 blood-brain barrier penetrating alkylating agents, receptor tyrosine kinase inhibitors (TKIs), and HDAC inhibitors.
- **Quantum Conceptual DFT (CDFT)**: Electronic hardness ($\eta$), chemical potential ($\mu$), electrophilicity ($\omega$), and adsorption energetics across pristine ($\text{Ti}_3\text{C}_2\text{O}_2$) and functionalized MXene monolayers.
- **Explainable Machine Learning (XAI / SHAP)**: Benchmarking across Random Forest, Extra Trees, Gradient Boosting, SVR, and Ridge regression with strict nested cross-validation and OECD Principle 3 applicability domain (Williams plot).
- **Surface Functionalization Physics**: Unraveling how surface termination chemistry ($\text{-O}$, $\text{-OH}$, $\text{-F}$) modulates charge transfer and therapeutic retention.

---

## 🔬 Repository Architecture

```
├── data/
│   ├── processed/                             # Processed datasets and descriptor matrices
│   └── raw/                                   # PDB 4UV7 receptor and 37 ligand PDBQT coordinates
├── figures/                                   # High-resolution publication figures (300 DPI)
├── manuscript/
│   └── Beilstein_Manuscript_GBM_MXene_Monreal_Hernandez_et_al.docx
├── results/
│   ├── docking/                               # Real Vina binding scores and contact residues
│   └── models/                                # QSAR benchmark summaries and SHAP rankings
├── src/
│   ├── descriptors/                           # CDFT & molecular descriptor computation
│   ├── docking/                               # Docking execution & residue contact extraction
│   ├── ml_models/                             # QSAR regression & applicability domain scripts
│   └── visualization/                         # Manuscript & figure compilation pipelines
├── run_entire_gbm_study.py                    # Master execution workflow
└── README.md
```

---

## ⚙️ Quickstart & Execution

```bash
git clone https://github.com/sircalch/mxene-glioblastoma-qsar-ai.git
cd mxene-glioblastoma-qsar-ai

# Create virtual environment & install requirements
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install rdkit scikit-learn pandas numpy scipy matplotlib seaborn shap python-docx

# Execute end-to-end reproducible pipeline
python run_entire_gbm_study.py
```

---

## 🗓️ v2.0.0 (2026-09-10)

- **Docking source unified**: `run_gbm_real_docking.py` now docks the human
  EGFR kinase domain **PDB 4ZAU** (the primary receptor cited throughout the
  manuscript; previously the script targeted 4UV7) with a fixed seed, and writes
  `vina_4ZAU_kcal_mol` straight into `dataset_drug_mxene_pristine.csv`.
  `analyze_gbm_interactions.py` reads 4ZAU to match the poses.
- Seeded run: n = 33 dockable compounds, -3.8 to -9.2 kcal/mol, mean -7.1;
  descriptor QSPR Q2_CV = 0.31 (EGFR docking, at best weakly predictive) / 0.10
  (MXene interaction energy).
- **Quantum step hardened**: `execute_gbm_energetics_rigorous.py` now catches
  xtb timeouts, resumes per orientation, and honours a `GBM_ADS_REBUILD` guard.
- `run_entire_gbm_study.py` verified end-to-end (9/9 steps); manuscript docking
  statistics computed from the CSV.

---

## 📜 Citation

```bibtex
@article{MonrealHernandez2026_GBM_MXene,
  title={Explainable AI and Quantum Chemical Exploration of 2D Titanium Carbide MXene (Ti3C2Tx) Nanocarriers for Glioblastoma Therapeutics},
  author={Monreal Hern{\'a}ndez, Andr{\'e}s and Franco Amaya, Sara Lizbeth and Mart{\'i}nez Osorio, Carlos Ivanhoe},
  journal={Beilstein Journal of Nanotechnology / Submitted},
  year={2026},
  version={2.0.0},
  doi={10.5281/zenodo.22187857},
  url={https://github.com/sircalch/mxene-glioblastoma-qsar-ai}
}
```

## 📄 License
Released under the [MIT License](LICENSE).
