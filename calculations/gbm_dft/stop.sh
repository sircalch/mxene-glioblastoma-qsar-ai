#!/bin/bash
pkill -f campaign.sh; pkill -f followup.sh; pkill -f run_fragments.sh; pkill -f run_checks.sh; pkill -f "bash run_all.sh"; pkill -f pw.x; pkill -f prterun; sleep 3
echo "pw.x left: $(pgrep -c pw.x)"
