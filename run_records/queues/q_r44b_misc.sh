#!/bin/bash
# (2026-10-28 UTC) round44 2 판: r44_pairw → r44_p1win spread,stats (작업자 3)
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs; C=$R/code
cd /content; export CUDA_VISIBLE_DEVICES=
echo "$(date -u '+%F %T') pairw·p1win 시작" >> $L/r44.log
nice -n 5 python3 $C/r44_pairw.py $R >> $L/r44_pairw.log 2>&1
nice -n 5 python3 $C/r44_p1win.py $R spread,stats 3 >> $L/r44_p1win.log 2>&1
echo "$(date -u '+%F %T') pairw·p1win 끝" >> $L/r44.log
