#!/bin/bash
cd "/mnt/c/Users/Andre/Proyectos doctorado/nano-qsar-ai-papers/mxene-glioblastoma-qsar-ai/calculations/gbm_dft"
pgrep -f "bash campaign.sh" > /dev/null && { echo "campaign already running"; exit 0; }
setsid nohup bash campaign.sh > campaign.out 2>&1 < /dev/null &
sleep 60; echo "pw.x procs: $(pgrep -c pw.x)"; free -g | sed -n 2p
