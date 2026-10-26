# (2026-10-20 UTC) round39 CPU 2: 줄기 336 D (R39_GPU_DONE) 와 CPU 1 (pr) 이 끝나면 E0′·E1·E2·E3·E4 → 6 절 stats.
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs; C=$R/code
cd /content
until grep -q R39_GPU_DONE $L/r39_D_stem.log 2>/dev/null && grep -q R39_CPU1_DONE $L/r39_core.log 2>/dev/null; do sleep 120; done
/usr/bin/python3 $C/r39_core.py $R e0,e1,e2,e3,e4 >> $L/r39_core.log 2>&1
/usr/bin/python3 $C/r39_sec6.py $R stats >> $L/r39_sec6.log 2>&1
echo R39_CPU2_DONE >> $L/r39_core.log
