#!/bin/bash
# waits for run_all.sh (complexes) to finish, then relaxes the slab (nosym input)
cd "/mnt/c/Users/Andre/Proyectos doctorado/nano-qsar-ai-papers/mxene-glioblastoma-qsar-ai/calculations/gbm_dft"
export PATH=$HOME/miniforge3/envs/qe/bin:$PATH OMP_NUM_THREADS=1
while pgrep -f "bash run_all.sh" > /dev/null; do sleep 300; done
if ! grep -q 'JOB DONE' slab/pw.out 2>/dev/null; then
  (cd slab && mpirun --oversubscribe -np 12 pw.x -nk 1 -in pw.in > pw.out 2>&1; rm -rf tmp)
fi
echo "followup done $(date)" >> followup.log
