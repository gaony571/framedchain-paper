# (2026-10-20 UTC) round39 CPU 1: E5 재료(R3′·R3″) → E5 방향 검정 → E1 유효 차원(pr). 줄기 336 D 를 쓰는 e0–e4 는 R39_GPU_DONE 뒤.
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs; C=$R/code
cd /content
/usr/bin/python3 $C/r39_e5prep.py $R r3p,r3pp >> $L/r39_e5prep.log 2>&1
/usr/bin/python3 $C/r39_core.py $R e5 >> $L/r39_core.log 2>&1
/usr/bin/python3 $C/r39_core.py $R pr >> $L/r39_core.log 2>&1
echo R39_CPU1_DONE >> $L/r39_core.log
