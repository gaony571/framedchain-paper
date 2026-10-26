#!/bin/bash
# 고침 2: nvcc 11.8 은 gcc ≤ 11 만 받음(이 런타임은 gcc 13.3) → gcc-11 로 빌드
set -x
apt-get -qq update && apt-get -qq install -y gcc-11 g++-11
E=/content/afenv; export PATH=$E/bin:$PATH
export CC=/usr/bin/gcc-11 CXX=/usr/bin/g++-11 CUDA_HOME=$E TORCH_CUDA_ARCH_LIST="8.0"
cd /content/openfold_src && git log -1 --oneline && rm -rf build dist
$E/bin/pip install --no-deps --no-build-isolation . 2>&1 | tail -15
cd /content && $E/bin/python -c "import torch, openfold; print('TORCH', torch.__version__, torch.cuda.is_available(), torch.cuda.get_device_name(0)); import attn_core_inplace_cuda; print('OPENFOLD_CUDA_OK')"
$E/bin/pip freeze > /content/drive/MyDrive/atlas_framed_curves/results/logs/afenv_freeze.txt
echo AFENV3_DONE
