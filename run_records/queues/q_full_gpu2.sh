# (2026-10-19) 트랙 B 첫 세션(code/r38_s0.py)을 7c 3 단계 뒤·4 단계 앞에 끼움(연구 계획 3.1: 트랙 B 가 준비되면 GPU 를 넘김). 3 단계는 끊지 않고 이어 돎.
while kill -0 8579 2>/dev/null; do sleep 30; done
echo STEP3_DONE >> /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_asamt_gen.log
cd /content
/usr/bin/python3 /content/drive/MyDrive/atlas_framed_curves/code/r38_s0.py /content/drive/MyDrive/atlas_framed_curves b0,e7,trunc,ens,bias,ustat,unc >> /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_s0.log 2>&1
echo S0_GPU_DONE >> /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_s0.log
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
