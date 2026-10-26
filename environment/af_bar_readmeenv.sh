#!/bin/bash
cd /content/drive/MyDrive/atlas_framed_curves/code/ext/alphaflow && OMP_NUM_THREADS=1 /content/afenv/bin/python -m scripts.analyze_ensembles --atlas_dir /content/r38_bar/atlas --pdbdir /content/r38_bar_env/full/alphaflow_md_base_202402 --num_workers 8
cp /content/r38_bar_env/full/alphaflow_md_base_202402/out.pkl /content/drive/MyDrive/atlas_framed_curves/derived/r38/af_bar/alphaflow_md_base_202402_readmeenv_out.pkl
cd /content/drive/MyDrive/atlas_framed_curves/code/ext/alphaflow && /content/afenv/bin/python -m scripts.print_analysis /content/drive/MyDrive/atlas_framed_curves/derived/r38/af_bar/alphaflow_md_base_202402_readmeenv_out.pkl > /content/drive/MyDrive/atlas_framed_curves/results/r38_af_bar_readmeenv.txt 2>&1
mkdir -p /content/drive/MyDrive/atlas_framed_curves/results/_backup/$(date +%F) && cp /content/drive/MyDrive/atlas_framed_curves/results/r38_af_bar_readmeenv.txt /content/drive/MyDrive/atlas_framed_curves/results/_backup/$(date +%F)/$(date +%H%M%S)_r38_af_bar_readmeenv.txt
echo AFBAR_READMEENV_DONE
