# (2026-10-20 되살림) P1 mean → align → spread → sim → pre, 작업자 2
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs; C=$R/code
cd /content
until [ -f /content/BASE_READY ]; do sleep 30; done
CUDA_VISIBLE_DEVICES= nice -n 5 /usr/bin/python3 $C/r38_p1.py $R mean,align,spread,sim,pre 2 >> $L/r38_p1.log 2>&1
echo P1_PRE_DONE_Q >> $L/r38_p1.log
