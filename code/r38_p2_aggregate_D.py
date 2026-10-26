"""동결된 r38_p2_aggregate.py (md5 bd8657fd…) 를 고치지 않고, 입력만 derived/r38/compare_D (P2 D 판, r38_p2_D.py) 로 바꿔 돌린다.
고침 round36_amend_p2_D_2026-10-19. 사용: python r38_p2_aggregate_D.py <ROOT> <목록> <TAG> [창] [Q]  → results/r38_p2_<TAG>.json"""
import sys, hashlib
ROOT = sys.argv[1]; src = open(ROOT + '/code/r38_p2_aggregate.py').read()
assert hashlib.md5(src.encode()).hexdigest() == 'bd8657fdba6eabe318ac2f428c5a9222', '집계 스크립트가 동결 판과 다름'
a = "CMP=ROOT+'/derived/r38/compare'"; assert src.count(a) == 1
exec(compile(src.replace(a, "CMP=ROOT+'/derived/r38/compare_D'"), 'r38_p2_aggregate.py (입력 compare_D)', 'exec'))
