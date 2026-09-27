#!/bin/bash
cd "/mnt/c/Users/Andre/Proyectos doctorado/nano-qsar-ai-papers/mxene-glioblastoma-qsar-ai/calculations/gbm_dft"
NP=12 setsid nohup bash run_all.sh > run_all.log 2>&1 < /dev/null &
sleep 90
echo "pw.x procs: $(pgrep -c pw.x)"; free -g | head -2
grep -E "Estimated|RAM" slab/pw.out | head -4
