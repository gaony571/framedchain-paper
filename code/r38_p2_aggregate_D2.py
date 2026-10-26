"""동결된 r38_p2_aggregate.py (md5 bd8657fd…) 를 고치지 않고, 입력만 derived/r38/compare_D2 (P2 D 판 + AlphaFlow D 칸, r38_p2_D_af.py) 와
results/r38_dff_D.json (R3′ d_ff 의 D 판, 7g (b)) 으로 바꿔 돌린다 — 7c 2차 2·3 의 D 판. 1차 P2·2차 1·4 는 compare_D 와 같은 칸이라 r38_p2_aggregate_D.py 와 같은 값이다.
사용: python r38_p2_aggregate_D2.py <ROOT> <목록> <TAG> [창] [Q] → results/r38_p2_<TAG>.json"""
import sys, hashlib
ROOT = sys.argv[1]; src = open(ROOT + '/code/r38_p2_aggregate.py').read()
assert hashlib.md5(src.encode()).hexdigest() == 'bd8657fdba6eabe318ac2f428c5a9222', '집계 스크립트가 동결 판과 다름'
a = "CMP=ROOT+'/derived/r38/compare'"; b = "'/results/r38_dff.json'"; assert src.count(a) == 1 and src.count(b) == 2
exec(compile(src.replace(a, "CMP=ROOT+'/derived/r38/compare_D2'").replace(b, "'/results/r38_dff_D.json'"), 'r38_p2_aggregate.py (입력 compare_D2, d_ff D)', 'exec'))
