#!/bin/bash
# (2026-10-28 UTC) round44 6 판: r44_redsame run → sum (등록 = results/prereg_2026-10-28_addendum6.json 11:50:14 UTC)
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs; C=$R/code
cd /content
echo "$(date -u '+%F %T') redsame 시작 (r44_redsame.py md5 $(md5sum $C/r44_redsame.py | cut -c1-8))" >> $L/r44.log
python3 $C/r44_redsame.py $R run,sum >> $L/r44_redsame.log 2>&1
echo "$(date -u '+%F %T') redsame 끝" >> $L/r44.log
