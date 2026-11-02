#!/bin/bash
# round41 g5 연산량 측정 대기열 (2026-10-22 UTC): 준비 끝 → 단계 차례로. 측정 중 다른 계산 없음.
R=/content/drive/MyDrive/atlas_framed_curves; T=/content/ct
while ! grep -q PREP_DONE $T/prep.log 2>/dev/null; do sleep 15; done
cd /content
taskset -c 1 python3 $R/code/r41_compute.py $R cpu >> $T/run.log 2>&1
python3 $R/code/r41_compute.py $R foldseek >> $T/run.log 2>&1
python3 $R/code/r41_compute.py $R usava >> $T/run.log 2>&1
python3 $R/code/r41_compute.py $R gpu >> $T/run.log 2>&1
python3 $R/code/r41_compute.py $R ens >> $T/run.log 2>&1
python3 $R/code/r41_compute.py $R sum >> $T/run.log 2>&1
cp $T/run.log $T/prep.log $R/results/logs/ 2>/dev/null; mv $R/results/logs/run.log $R/results/logs/r41_compute.log; mv $R/results/logs/prep.log $R/results/logs/r41_compute_prep.log
echo CT_ALL_DONE >> $R/results/logs/r41_compute.log
