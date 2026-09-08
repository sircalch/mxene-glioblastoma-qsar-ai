"""
run_entire_gbm_study.py
Master end-to-end pipeline for Article 2 (Glioblastoma / 2D Ti3C2Tx MXene).
Reproduces every real number and figure in the manuscript from raw inputs.
"""
import os
import sys
import time

BASE = os.path.dirname(os.path.abspath(__file__))


def run_step(n, total, title, rel_path, args=""):
    script = os.path.join(BASE, rel_path)
    print(f"\n{'='*70}\n  [Step {n}/{total}] {title}\n{'='*70}")
    t0 = time.time()
    ret = os.system(f'python "{script}" {args}')
    if ret != 0:
        print(f"[ERROR] Step {n}: {title} (exit {ret})")
        return False
    print(f"[OK] Step {n} in {time.time()-t0:.1f}s")
    return True


def main():
    print("=" * 70)
    print("  GLIOBLASTOMA / Ti3C2Tx MXene : MASTER REPRODUCIBILITY PIPELINE")
    print("=" * 70)
    steps = [
        ("Drug-library curation", "src/descriptors/curate_gbm_dataset.py"),
        ("RDKit + GFN2-xTB descriptors", "src/descriptors/compute_gbm_descriptors.py"),
        ("Real AutoDock Vina docking (EGFR, PDB 4ZAU)", "src/docking/run_gbm_real_docking.py"),
        ("Residue-level contact analysis", "src/docking/analyze_gbm_interactions.py"),
        ("GFN2-xTB adsorption on Ti3C2O2 (relaxed complexes)", "execute_gbm_energetics_rigorous.py"),
        ("OECD applicability domain (Williams)", "src/ml_models/compute_gbm_oecd_applicability_domain.py"),
        ("Master figure suite", "src/visualization/generate_gbm_master_figures.py"),
        ("Word manuscript", "src/visualization/generate_gbm_word_manuscript.py"),
        ("Supporting information", "src/visualization/generate_supporting_information.py"),
    ]
    for i, (title, path) in enumerate(steps, 1):
        if not run_step(i, len(steps), title, path):
            sys.exit(1)
    print("\n" + "=" * 70)
    print(">>> PIPELINE COMPLETE <<<")
    print("  manuscript/Beilstein_Manuscript_GBM_MXene_Monreal_Hernandez_et_al.docx")
    print("=" * 70)


if __name__ == "__main__":
    main()
