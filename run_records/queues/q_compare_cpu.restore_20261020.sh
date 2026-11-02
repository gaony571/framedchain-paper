until [ -f /content/BASE_READY ]; do sleep 30; done
until grep -q STEP4_DONE /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_alphaflow_gen.log 2>/dev/null; do sleep 120; done
cd /content && /usr/bin/python3 /content/drive/MyDrive/atlas_framed_curves/code/r38_compare.py /content/drive/MyDrive/atlas_framed_curves /content/drive/MyDrive/atlas_framed_curves/derived/r38/q3_p2_67_ordered.txt > /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_compare_67.log 2>&1
echo COMPARE67_DONE >> /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_compare_67.log
until grep -q STEP5_DONE /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_asamt_gen.log 2>/dev/null; do sleep 120; done
/usr/bin/python3 /content/drive/MyDrive/atlas_framed_curves/code/r38_compare.py /content/drive/MyDrive/atlas_framed_curves /content/drive/MyDrive/atlas_framed_curves/derived/r38/rest18.txt > /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_compare_rest18.log 2>&1
echo COMPARE_REST18_DONE >> /content/drive/MyDrive/atlas_framed_curves/results/logs/r38_compare_rest18.log
