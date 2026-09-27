#!/bin/bash
pkill -f followup.sh; pkill -f "bash run_all.sh"; pkill -f pw.x; pkill -f prterun; sleep 3
echo "pw.x left: $(pgrep -c pw.x)"
