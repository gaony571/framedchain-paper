# (2026-10-20 UTC) round39 D 판(맨 끝): FULL_GPU_DONE → 줄기 336 D (E6 접힘 대 풀림 D 는 뺌, 사용자 "있으면 좋음 까지만")
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs; C=$R/code
cd /content
until grep -q FULL_GPU_DONE $L/r38_alphaflow_gen.log 2>/dev/null; do sleep 120; done
/usr/bin/python3 $C/r39_D_stem.py $R run,sum >> $L/r39_D_stem.log 2>&1
echo R39_GPU_DONE >> $L/r39_D_stem.log
