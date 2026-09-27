#!/bin/bash
echo "uptime: $(uptime)"
echo "pw.x procs: $(pgrep -c pw.x)"
ps -eo pid,etime,cmd | grep -E "run_all|mpirun|sleep" | grep -v grep | head
