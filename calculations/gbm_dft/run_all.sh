#!/bin/bash
export PATH=$HOME/miniforge3/envs/qe/bin:$PATH
export OMP_NUM_THREADS=1
NP=${NP:-12}
cd "/mnt/c/Users/Andre/Proyectos doctorado/nano-qsar-ai-papers/mxene-glioblastoma-qsar-ai/calculations/gbm_dft"
if ! grep -q 'JOB DONE' slab/pw.out 2>/dev/null; then (cd slab && mpirun --oversubscribe -np $NP pw.x -nk 1 -in pw.in > pw.out 2>&1; rm -rf tmp); fi
if ! grep -q 'JOB DONE' drug_Temozolomide/pw.out 2>/dev/null; then (cd drug_Temozolomide && mpirun --oversubscribe -np $NP pw.x -nk 1 -in pw.in > pw.out 2>&1; rm -rf tmp); fi
if ! grep -q 'JOB DONE' cplx_Temozolomide/pw.out 2>/dev/null; then (cd cplx_Temozolomide && mpirun --oversubscribe -np $NP pw.x -nk 1 -in pw.in > pw.out 2>&1; rm -rf tmp); fi
if ! grep -q 'JOB DONE' drug_Carmustine/pw.out 2>/dev/null; then (cd drug_Carmustine && mpirun --oversubscribe -np $NP pw.x -nk 1 -in pw.in > pw.out 2>&1; rm -rf tmp); fi
if ! grep -q 'JOB DONE' cplx_Carmustine/pw.out 2>/dev/null; then (cd cplx_Carmustine && mpirun --oversubscribe -np $NP pw.x -nk 1 -in pw.in > pw.out 2>&1; rm -rf tmp); fi
if ! grep -q 'JOB DONE' drug_Lomustine/pw.out 2>/dev/null; then (cd drug_Lomustine && mpirun --oversubscribe -np $NP pw.x -nk 1 -in pw.in > pw.out 2>&1; rm -rf tmp); fi
if ! grep -q 'JOB DONE' cplx_Lomustine/pw.out 2>/dev/null; then (cd cplx_Lomustine && mpirun --oversubscribe -np $NP pw.x -nk 1 -in pw.in > pw.out 2>&1; rm -rf tmp); fi
if ! grep -q 'JOB DONE' drug_Nimustine/pw.out 2>/dev/null; then (cd drug_Nimustine && mpirun --oversubscribe -np $NP pw.x -nk 1 -in pw.in > pw.out 2>&1; rm -rf tmp); fi
if ! grep -q 'JOB DONE' cplx_Nimustine/pw.out 2>/dev/null; then (cd cplx_Nimustine && mpirun --oversubscribe -np $NP pw.x -nk 1 -in pw.in > pw.out 2>&1; rm -rf tmp); fi
