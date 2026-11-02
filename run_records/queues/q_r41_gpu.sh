#!/bin/bash
# round41 GPU 대기열 (2026-10-21): 줄기 336 D 끝(R39_GPU_DONE) 뒤 g2 gpu → g2 sum → g4 → E6 run
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs
while ! grep -q R39_GPU_DONE $L/r39_D_stem.log 2>/dev/null; do sleep 60; done
cd /content
python3 $R/code/r41_dcomp.py $R gpu >> $L/r41_dcomp_gpu.log 2>&1
while ! grep -q G2_CPU_DONE_0 $L/r41_dcomp_cpu0.log 2>/dev/null || ! grep -q G2_CPU_DONE_1 $L/r41_dcomp_cpu1.log 2>/dev/null; do sleep 60; done
CUDA_VISIBLE_DEVICES= python3 $R/code/r41_dcomp.py $R sum >> $L/r41_dcomp_gpu.log 2>&1
python3 $R/code/r41_dff50.py $R >> $L/r41_dff50.log 2>&1
python3 $R/code/r39_D_unfold.py $R run >> $L/r39_unfold_D.log 2>&1
echo R41_GPU_DONE >> $L/r41_dff50.log
