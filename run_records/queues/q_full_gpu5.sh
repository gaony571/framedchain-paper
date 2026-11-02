# (2026-10-19 18:4x) 첫 세션 → E7 30 → P2 D 판(67) → P3 π̄ (트랙 B4; CPU dom 끝을 기다림) → 4 단계 AlphaFlow → 5 → 조건부. 첫 세션 python 은 끊지 않음.
while kill -0 44493 2>/dev/null; do sleep 30; done
echo S0_GPU_DONE >> /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_s0.log
cd /content
/usr/bin/python3 /content/drive/MyDrive/atlas_framed_curves/code/r38_e7_30.py /content/drive/MyDrive/atlas_framed_curves >> /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_e7_30.log 2>&1
echo E7_30_DONE_Q >> /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_e7_30.log
/usr/bin/python3 /content/drive/MyDrive/atlas_framed_curves/code/r38_p2_D.py /content/drive/MyDrive/atlas_framed_curves /content/drive/MyDrive/atlas_framed_curves/derived/r38/q3_p2_67_ordered.txt >> /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_p2_D.log 2>&1
/usr/bin/python3 /content/drive/MyDrive/atlas_framed_curves/code/r38_p2_aggregate_D.py /content/drive/MyDrive/atlas_framed_curves /content/drive/MyDrive/atlas_framed_curves/derived/r38/q3_p2_67_ordered.txt p2_67_D full 0.75 >> /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_p2_D.log 2>&1
echo P2_D_DONE_Q >> /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_p2_D.log
until grep -q DOM_DONE /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_p3.log 2>/dev/null; do sleep 60; done
/usr/bin/python3 /content/drive/MyDrive/atlas_framed_curves/code/r38_p3.py /content/drive/MyDrive/atlas_framed_curves pibar >> /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_p3.log 2>&1
echo PIBAR_DONE_Q >> /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_p3.log
until [ -f /content/AFENV_READY ]; do sleep 60; done
/usr/bin/python3 /content/drive/MyDrive/atlas_framed_curves/code/r38_alphaflow_gen.py /content/drive/MyDrive/atlas_framed_curves /content/drive/MyDrive/atlas_framed_curves/derived/r38/q3_p2_67_ordered.txt >> /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_alphaflow_gen.log 2>&1
echo STEP4_DONE >> /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_alphaflow_gen.log
/usr/bin/python3 /content/drive/MyDrive/atlas_framed_curves/code/r38_asamt_gen.py /content/drive/MyDrive/atlas_framed_curves /content/drive/MyDrive/atlas_framed_curves/derived/r38/q5a_450_17.txt 450 50,25,10 >> /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_asamt_gen.log 2>&1
/usr/bin/python3 /content/drive/MyDrive/atlas_framed_curves/code/r38_asamt_gen.py /content/drive/MyDrive/atlas_framed_curves /content/drive/MyDrive/atlas_framed_curves/derived/r38/q5b_p90_4.txt 320 50,25,10 >> /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_asamt_gen.log 2>&1
/usr/bin/python3 /content/drive/MyDrive/atlas_framed_curves/code/r38_asamt_gen.py /content/drive/MyDrive/atlas_framed_curves /content/drive/MyDrive/atlas_framed_curves/derived/r38/rest18.txt 320 50,25,10 >> /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_asamt_gen.log 2>&1
echo STEP5_DONE >> /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_asamt_gen.log
# 조건부(시간 남으면; 트랙 B 뒤)
/usr/bin/python3 /content/drive/MyDrive/atlas_framed_curves/code/r38_asamt_gen.py /content/drive/MyDrive/atlas_framed_curves /content/drive/MyDrive/atlas_framed_curves/derived/r38/q_cond_450_50.txt 450 50,25,10 >> /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_asamt_gen.log 2>&1
/usr/bin/python3 /content/drive/MyDrive/atlas_framed_curves/code/r38_asamt_gen.py /content/drive/MyDrive/atlas_framed_curves /content/drive/MyDrive/atlas_framed_curves/derived/r38/q3_p2_67_ordered.txt 348,379,413 50,25,10 >> /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_asamt_gen.log 2>&1
/usr/bin/python3 /content/drive/MyDrive/atlas_framed_curves/code/r38_alphaflow_gen.py /content/drive/MyDrive/atlas_framed_curves /content/drive/MyDrive/atlas_framed_curves/derived/r38/rest18.txt >> /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_alphaflow_gen.log 2>&1
echo FULL_GPU_DONE >> /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_alphaflow_gen.log
