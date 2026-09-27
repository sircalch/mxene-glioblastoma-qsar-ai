#!/bin/bash
free -m; vmstat 1 3 | tail -2; top -bn1 | head -15
