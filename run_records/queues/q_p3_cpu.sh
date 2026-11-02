# (2026-10-19) 트랙 B4 P3 CPU 쪽: π̄ 끝 → 쌍 계산 → 통계
until grep -q PIBAR_DONE_Q /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_p3.log 2>/dev/null; do sleep 120; done
cd /content && CUDA_VISIBLE_DEVICES= nice -n 5 /usr/bin/python3 /content/drive/MyDrive/atlas_framed_curves/code/r38_p3.py /content/drive/MyDrive/atlas_framed_curves pair,stats >> /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_p3.log 2>&1
echo P3_CPU_DONE >> /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_p3.log
