#!/bin/bash
# round39 E6 위치 판정 (2026-10-21): E6 run 끝(R41_GPU_DONE) 뒤 CPU
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs
while ! grep -q R41_GPU_DONE $L/r41_dff50.log 2>/dev/null; do sleep 120; done
cd /content && CUDA_VISIBLE_DEVICES= python3 $R/code/r39_e6.py $R >> $L/r39_e6.log 2>&1
