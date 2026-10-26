#!/bin/bash
# 고침 3: 현재 sam2 코드(ce076ed)는 PEP 701 f-string(같은 따옴표 겹침)을 써서 Python ≥ 3.12 가 필요 → 3.12 환경으로 다시. torch 2.6.0(sam2 시험 상한)
set -x
export MAMBA_ROOT_PREFIX=/content/mamba
E=/content/samenv; rm -rf $E
/content/bin/micromamba create -y -p $E -c conda-forge python=3.12
P=$E/bin/pip
$P install torch==2.6.0 --index-url https://download.pytorch.org/whl/cu124
$P install -e /content/sam2_src
$E/bin/python -c "import torch, sam, diffusers, mdtraj; print('SAMCHECK', torch.__version__, torch.version.cuda, torch.cuda.is_available(), diffusers.__version__, mdtraj.__version__)"
$P freeze > /content/drive/MyDrive/atlas_framed_curves/results/logs/samenv_freeze.txt
echo SAMENV3_DONE
