"""
run_entire_gbm_study.py
Master End-to-End Pipeline Runner for 100% Reproducibility of Article 2:
Glioblastoma Therapeutics & 2D Ti3C2Tx MXene Nanocarriers.
"""

import os
import sys
import time

def run_step(step_num, title, script_rel_path):
    base_dir = os.path.dirname(os.path.abspath(__file__))
    script_path = os.path.join(base_dir, script_rel_path)
    print(f"\n=======================================================")
    print(f"  [Step {step_num}/8] {title}")
    print(f"=======================================================")
    t0 = time.time()
    ret = os.system(f'python "{script_path}"')
    t_elapsed = time.time() - t0
    if ret != 0:
        print(f"[ERROR] Step {step_num}: {title} (Exit Code: {ret})")
        return False
    print(f"[OK] Step {step_num} completed in {t_elapsed:.2f} seconds.")
    return True

def main():
    print("=" * 65)
    print("  MXENE-GLIOBLASTOMA-QSAR-AI: MASTER REPRODUCIBILITY PIPELINE")
    print("  Authors: Andrés Monreal Hernández et al.")
    print("=" * 65)
    
    steps = [
        (1, "GBM Drug Library Curation", "src/descriptors/curate_gbm_dataset.py"),
        (2, "20-Descriptor RDKit & Quantum Calculation", "src/descriptors/compute_gbm_descriptors.py"),
        (3, "Parallel Real AutoDock Vina Docking (PDB 4UV7)", "src/docking/run_gbm_real_docking.py"),
        (4, "Residue-Level Contact Analysis", "src/docking/analyze_gbm_interactions.py"),
        (5, "Machine Learning Training & SHAP XAI", "src/ml_models/train_gbm_qsar_models.py"),
        (6, "OECD Applicability Domain (Williams Plot)", "src/ml_models/compute_gbm_oecd_applicability_domain.py"),
        (7, "Publication-Grade Figures Suite (300+ DPI)", "src/visualization/generate_gbm_master_figures.py"),
        (8, "Word Manuscript Compilation & Submission Packaging", "src/visualization/generate_gbm_word_manuscript.py")
    ]
    
    for s_num, title, path in steps:
        success = run_step(s_num, title, path)
        if not success:
            sys.exit(1)
            
    print("\n" + "=" * 65)
    print(">>> FULL REPRODUCIBILITY PIPELINE EXECUTED SUCCESSFULLY! <<<")
    print("  Manuscript Word File: manuscript/Beilstein_Manuscript_GBM_MXene_Monreal_Hernandez_et_al.docx")
    print("  Submission ZIP File:  mxene-glioblastoma-qsar-ai-FINAL-SUBMISSION-READY.zip")
    print("=" * 65)

if __name__ == "__main__":
    main()
