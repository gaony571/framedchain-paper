#!/bin/bash
# (2026-10-29 UTC) round45: r45_f1 corr 8 조각 동시 → 모두 끝나면 sanity → stats (등록 = results/prereg_2026-10-29_addendum7.json)
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs; C=$R/code
cd /content; export CUDA_VISIBLE_DEVICES=
echo "$(date -u '+%F %T') r45_f1 corr 8 조각 시작 (r45_f1.py md5 $(md5sum $C/r45_f1.py | cut -c1-8))" >> $L/r45.log
P=""
for k in 0 1 2 3 4 5 6 7; do nice -n 5 python3 $C/r45_f1.py $R corr 8 $k >> $L/r45_f1_s$k.log 2>&1 & P="$P $!"; done
wait $P
echo "$(date -u '+%F %T') r45_f1 corr 조각 모두 끝" >> $L/r45.log
python3 $C/r45_f1.py $R sanity >> $L/r45_f1_sanity.log 2>&1
echo "$(date -u '+%F %T') r45_f1 sanity 끝" >> $L/r45.log
python3 $C/r45_f1.py $R stats >> $L/r45_f1.log 2>&1
echo "$(date -u '+%F %T') r45_f1 stats 끝" >> $L/r45.log
