# (2026-10-20) 트랙 B5 1차 P1 CPU 쪽: mean → align → spread → sim → pre (판정 stats 는 점검표를 채운 뒤 따로)
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs; C=$R/code
cd /content && CUDA_VISIBLE_DEVICES= nice -n 5 /usr/bin/python3 $C/r38_p1.py $R mean,align,spread,sim,pre 3 >> $L/r38_p1.log 2>&1
echo P1_PRE_DONE_Q >> $L/r38_p1.log
