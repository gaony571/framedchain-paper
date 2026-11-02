#!/bin/bash
# round41 CPU 대기열 2 (2026-10-21): g2 cpu 조각 2 (낮은 우선순위)
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs
cd /content
nice -n 10 python3 $R/code/r41_dcomp.py $R cpu 2 0 >> $L/r41_dcomp_cpu0.log 2>&1 &
nice -n 10 python3 $R/code/r41_dcomp.py $R cpu 2 1 >> $L/r41_dcomp_cpu1.log 2>&1 &
wait
echo R41_CPU2_DONE >> $L/r41_dcomp_cpu0.log
