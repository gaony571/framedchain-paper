cd /content
/usr/bin/python3 /content/drive/MyDrive/atlas_framed_curves/code/r38_alphaflow_gen.py /content/drive/MyDrive/atlas_framed_curves /content/drive/MyDrive/atlas_framed_curves/derived/r38/q4_af_par.txt > /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_alphaflow_par.log 2>&1
echo STEP4_DONE >> /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_alphaflow_par.log
# 조건부
/usr/bin/python3 /content/drive/MyDrive/atlas_framed_curves/code/r38_alphaflow_gen.py /content/drive/MyDrive/atlas_framed_curves /content/drive/MyDrive/atlas_framed_curves/derived/r38/rest18.txt >> /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_alphaflow_par.log 2>&1
echo AF_ALL_DONE >> /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_alphaflow_par.log
