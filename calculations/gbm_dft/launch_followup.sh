#!/bin/bash
cd "/mnt/c/Users/Andre/Proyectos doctorado/nano-qsar-ai-papers/mxene-glioblastoma-qsar-ai/calculations/gbm_dft"
setsid nohup bash followup.sh > followup.log 2>&1 < /dev/null &
sleep 3; pgrep -fa followup.sh | head -2
