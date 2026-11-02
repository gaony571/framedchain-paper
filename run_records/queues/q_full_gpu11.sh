# (2026-10-20 UTC) q_full_gpu10 대체(사용자 "있으면 좋음 까지만"): AlphaFlow D 칸(PID 182677) 끝 → 집계 D2 → 정적 축 → FULL_GPU_DONE (조건부 생성 뺌)
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs; C=$R/code
cd /content
while kill -0 182677 2>/dev/null; do sleep 30; done
for t in "q3_p2_67_ordered p2_67_D2" "p2_sf23 p2_sf23_D2" "p2_no_big4 p2_no_big4_D2" "p2_p80 p2_p80_D2"; do set -- $t; /usr/bin/python3 $C/r38_p2_aggregate_D2.py $R $R/derived/r38/$1.txt $2 full 0.75 >> $L/r38_p2_D_af.log 2>&1; done
echo P2_D2_DONE >> $L/r38_p2_D_af.log
until grep -q STATIC_DOM_DONE $L/r38_static.log 2>/dev/null; do sleep 60; done
/usr/bin/python3 $C/r38_static.py $R align,bnb,stats 1 >> $L/r38_static.log 2>&1
echo G7_DONE >> $L/r38_static.log
echo FULL_GPU_DONE >> $L/r38_alphaflow_gen.log
