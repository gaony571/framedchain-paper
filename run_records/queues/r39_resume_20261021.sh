#!/bin/bash
# (2026-10-21 UTC) Colab 다시 켠 뒤 재개: base 패키지 → 줄기 336 D 이어 돌기(끝난 도메인 건너뜀) → round39 CPU 2 (e0–e4 + 6 절 stats)
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs; C=$R/code
cd /content
echo "$(date '+%F %T') 재개 시작 GPU=$(nvidia-smi --query-gpu=name --format=csv,noheader) CPU=$(nproc)" >> $L/r39_resume.log
pip install -q mdtraj pot tmtools >> $L/r39_resume.log 2>&1
python3 -c "import mdtraj, h5py, numpy, tmtools, ot, torch; print('BASECHECK', mdtraj.__version__, numpy.__version__, torch.__version__, torch.cuda.is_available())" >> $L/r39_resume.log 2>&1 && touch /content/BASE_READY
rm -f $R/derived/r39/D/*.tmp
cp $C/q_r39_gpu2.restore_20261020.sh /content/q_r39_gpu2.sh; cp $C/q_r39_cpu2.restore_20261020.sh /content/q_r39_cpu2.sh
setsid nohup bash /content/q_r39_gpu2.sh > /content/q_r39_gpu2.out 2>&1 < /dev/null &
setsid nohup bash /content/q_r39_cpu2.sh > /content/q_r39_cpu2.out 2>&1 < /dev/null &
echo "$(date '+%F %T') 대기열 띄움 (줄기 D 끝난 도메인 $(ls $R/derived/r39/D/*.json | wc -l))" >> $L/r39_resume.log
