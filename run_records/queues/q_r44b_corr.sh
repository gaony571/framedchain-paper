#!/bin/bash
# (2026-10-28 UTC) round44 2 판: r44_corr_v2 corr 8 조각 동시 → 모두 끝나면 stats (등록 = results/prereg_2026-10-28_addendum2.json 05:56:13 UTC)
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs; C=$R/code
cd /content; export CUDA_VISIBLE_DEVICES=
echo "$(date -u '+%F %T') corr_v2 8 조각 시작" >> $L/r44.log
P=""
for k in 0 1 2 3 4 5 6 7; do nice -n 5 python3 $C/r44_corr_v2.py $R corr 8 $k >> $L/r44_corr_v2_s$k.log 2>&1 & P="$P $!"; done
wait $P
echo "$(date -u '+%F %T') corr_v2 조각 모두 끝" >> $L/r44.log
python3 $C/r44_corr_v2.py $R stats >> $L/r44_corr_v2.log 2>&1
echo "$(date -u '+%F %T') corr_v2 stats 끝" >> $L/r44.log
