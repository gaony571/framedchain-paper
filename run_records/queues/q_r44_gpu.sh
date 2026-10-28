#!/bin/bash
# (2026-10-28 UTC) round44 GPU 대기열: r44_inv frm (60 도메인, 프레임 50) → CPU 의 rep 이 끝나면 sum
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs; C=$R/code
cd /content
echo "$(date -u '+%F %T') GPU 대기열 시작 $(nvidia-smi --query-gpu=name --format=csv,noheader)" >> $L/r44.log
python3 $C/r44_inv.py $R frm >> $L/r44_inv_frm.log 2>&1
until grep -q R44_CPU_INV_DONE $L/r44_inv.log; do sleep 60; done
python3 $C/r44_inv.py $R sum >> $L/r44_inv.log 2>&1
echo "$(date -u '+%F %T') GPU 대기열 끝" >> $L/r44.log
