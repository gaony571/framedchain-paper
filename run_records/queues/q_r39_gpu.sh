# (2026-10-20 UTC) round39 D 판 (사용자 결정 "(가), 맨 마지막에"): q_full_gpu9 가 다 끝난 뒤(FULL_GPU_DONE) → 줄기 336 D → E6 접힘 대 풀림 D
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs; C=$R/code
cd /content
until grep -q FULL_GPU_DONE $L/r38_alphaflow_gen.log 2>/dev/null; do sleep 120; done
/usr/bin/python3 $C/r39_D_stem.py $R run,sum >> $L/r39_D_stem.log 2>&1
/usr/bin/python3 $C/r39_D_unfold.py $R run,sum >> $L/r39_D_unfold.log 2>&1
echo R39_GPU_DONE >> $L/r39_D_unfold.log
