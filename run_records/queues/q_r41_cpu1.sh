#!/bin/bash
# round41 CPU 대기열 1 (2026-10-21): g1 static → tm → rdyn
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs
while [ ! -f /content/BASE_READY ]; do sleep 30; done
cd /content && CUDA_VISIBLE_DEVICES= nice -n 5 python3 $R/code/r41_gapinfo.py $R static,tm,rdyn 4 >> $L/r41_gapinfo.log 2>&1
echo R41_CPU1_DONE >> $L/r41_gapinfo.log
