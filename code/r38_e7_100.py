"""(prereg round36 7g (a), 2026-10-20 UTC) 온도 눈금을 D 로: E7 30 과 첫 세션 E7 5 에 더해 무작위 65 도메인, 합 100. 계산은 r38_e7_30.py 의 one() 그대로(그 파일은 바꾸지 않음).
65 = 줄기 derived/mdcath_5T 336 − (첫 세션의 모든 도메인: E7 5, 잘라 낸 사슬 10, 28 쌍의 도메인 ∪ 7c 바름 점검 5) − E7 30 − E7 30 에서 건너뛴 도메인을 정렬해
default_rng(20261060).permutation 순서로, avail(320 K 접힌 프레임 ≥ 100 인 반복 ≥ 2)인 도메인 65. E7 5 는 같은 계산으로 다시 낸다(첫 세션 e7 부분과 출력 꼴이 달라서).
판정 없음(7g: 그림 5 의 공통 축을 D 로 그린다). 요약은 r38_e7_30.summarize 와 같은 비교를 100 도메인에서 + 도메인·칸별 값 표(그림 5 용).
사용: python r38_e7_100.py <ROOT> [pick,run,sum] 출력: derived/r38/s0/e7_100/{pick.json, <dom>.json}, results/r38_e7_100.json"""
import sys, os, json
ROOT = sys.argv[1]; PART = sys.argv[2].split(',') if len(sys.argv) > 2 else ['pick', 'run', 'sum']
_e = open(ROOT + '/code/r38_e7_30.py').read(); _e = _e[:_e.index('\nfor d in SEL:')]
sys.argv = ['r38_e7_30', ROOT]; exec(_e)
OUT3 = OUT + '/e7_100'; os.makedirs(OUT3, exist_ok=True)
_s = _e[_e.index('def summarize():'):]
assert _s.count("sorted(glob.glob(OUT2 + '/*.json'))") == 1 and _s.count("RES + '/r38_e7_30.json'") == 1
exec(_s.replace('def summarize():', 'def summarize100():').replace("sorted(glob.glob(OUT2 + '/*.json'))", 'FILES').replace("RES + '/r38_e7_30.json'", 'FO').replace("log('E7_30 요약'", "log('E7_100 요약'"))
FILES, FO = [], RES + '/r38_e7_100.json'


def part_pick():
    fo = OUT3 + '/pick.json'
    if os.path.exists(fo): log('pick 있음'); return
    pk = json.load(open(OUT + '/pick.json')); e30 = json.load(open(OUT + '/e7_30.json'))
    s5 = open(ROOT + '/derived/r38/asamt_sanity5.txt').read().split()
    first = set(pk['e7']) | set(pk['trunc']) | {x for p in pk['pairs'] for x in (p['a'], p['b'])} | set(s5)
    stem = sorted(os.path.basename(f)[:-4] for f in glob.glob(ROOT + '/derived/mdcath_5T/*.npz'))
    assert len(stem) == 336 and set(stem) & first == set(e30['excluded_in_trunk']), (len(stem), len(set(stem) & first)) # E7 30 을 고른 규칙과 같은 '첫 세션 도메인'
    pool = sorted(d for d in stem if d not in first and d not in set(e30['e7_30']) and d not in set(e30['skipped']))
    new, skipped = [], []
    for o in np.random.default_rng(20261060).permutation(len(pool)):
        d = pool[int(o)]; ok = avail(d); _runs.clear(); _md.clear()
        (new if ok else skipped).append(d)
        if len(new) == 65: break
    jdump(dict(written=time.strftime('%F %T'), seed=20261060, n_pool=len(pool), new65=sorted(new), skipped=skipped, e7_5=pk['e7'], e7_30=e30['e7_30'],
               rule='줄기 336 − 첫 세션 도메인(E7 30 의 excluded_in_trunk 20 과 같음) − E7 30 − E7 30 건너뜀을 정렬해 default_rng(20261060).permutation 순서로 avail 인 65'), fo, bk=True)
    log('E7_100 pick', len(new), '건너뜀', len(skipped), '풀', len(pool))


def part_run():
    pk = json.load(open(OUT3 + '/pick.json'))
    for d in pk['e7_5'] + pk['new65']:
        fo = f'{OUT3}/{d}.json'
        if os.path.exists(fo): continue
        try:
            r = one(d); jdump(r, fo); log('e7_100', d, len(r['cells']), '칸', round(r['sec']), 's')
        except Exception as e:
            traceback.print_exc(); log('ERR e7_100', d, repr(e)[:200])
    log('E7_100_RUN_DONE')


def part_sum():
    global FILES
    pk = json.load(open(OUT3 + '/pick.json')); src = {**{d: 'e7_30' for d in pk['e7_30']}, **{d: 'e7_5' for d in pk['e7_5']}, **{d: 'new65' for d in pk['new65']}}
    FILES = sorted([f'{OUT2}/{d}.json' for d in pk['e7_30'] if os.path.exists(f'{OUT2}/{d}.json')] + [f'{OUT3}/{d}.json' for d in pk['e7_5'] + pk['new65'] if os.path.exists(f'{OUT3}/{d}.json')])
    summarize100(); out = json.load(open(FO)); tab = []
    for f in FILES:
        j = json.load(open(f))
        for k, c in j.get('cells', {}).items():
            if k.startswith('MD'):
                tab.append(dict(dom=j['dom'], src=src[j['dom']], cell=k, L=j.get('L'), D2_D=c['D']['D2'], D2_A=c['Aall']['D2'], wbar_D=c['D']['wbar'], wbar_A=c['Aall']['wbar'],
                                ratio_D=c['D']['ratio'], ratio_A=c['Aall']['ratio'], dmu2_D=c['dmu2_D'], dmu_A=c['dmu_A']))
    doms = {r['dom'] for r in tab}
    out.update(table=tab, n_by_src={s: sum(1 for d in doms if src[d] == s) for s in ('e7_30', 'e7_5', 'new65')},
               verdict_note='7g: 판정 없음. verdict 칸은 3 절 규칙을 100 도메인에 대 본 서술이다(정식 판정은 E7 30: 실패).', pick_seed=pk['seed'])
    jdump(out, FO, bk=True); log('E7_100_DONE', len(doms), '도메인')


if __name__ == '__main__':
    for part in PART:
        log('부분 시작', part)
        try: globals()['part_' + part](); log('부분 끝', part)
        except Exception as e: traceback.print_exc(); log('ERR 부분', part, repr(e)[:300])
