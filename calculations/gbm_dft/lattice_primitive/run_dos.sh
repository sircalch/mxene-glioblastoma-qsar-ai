#!/bin/bash
export PATH=$HOME/miniforge3/envs/qe/bin:$PATH
export OMP_NUM_THREADS=1
cd "/mnt/c/Users/Andre/Proyectos doctorado/nano-qsar-ai-papers/mxene-glioblastoma-qsar-ai/calculations/gbm_dft/lattice_primitive"
mpirun --oversubscribe -np 8 pw.x -nk 2 -in scf.in > scf.out 2>&1
mpirun --oversubscribe -np 8 pw.x -nk 2 -in nscf.in > nscf.out 2>&1
mpirun --oversubscribe -np 4 dos.x -in dos.in > dos.out 2>&1
rm -rf tmp
echo done
