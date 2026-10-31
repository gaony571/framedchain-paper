"""(고침 = Drive results/prereg_2026-10-30_addendum10_fix2.json; 원 등록 addendum10 "r47_ties", 고침 fix1) 비교 단계(combine)만 고쳐 다시 낸다.
원 combine 은 부분 파일의 글자 칸(written, part, mode)까지 추정치로 읽으려다 TypeError 로 멈췄다. 고친 것: 추정치 칸(dict)만 읽는다.
fix2: fix1 (b18f1d03) 이 r47_ties.py 를 처음 나오는 "if __name__" 에서 잘라(part_anm 안의 문자열) 문법 오류로 멈췄으므로 마지막 자리(rindex)에서 자른다. 그 밖은 fix1 과 같다.
사용: python r47_combine_fix2.py <ROOT>   출력: results/r47_ties.json"""
import sys, os, json, time
ROOT = sys.argv[1]
_src = open(ROOT + '/code/r47_ties.py').read(); _src = _src[:_src.rindex("if __name__ == '__main__':")]
sys.argv = ['r47_ties', ROOT, 'combine']; exec(_src)   # r47_ties 의 triple, zside, dmax, log, OUT47, RES


def part_combine_fix1():
    J = lambda f: json.load(open(f'{OUT47}/{f}.json'))
    est_only = lambda d: {k: v for k, v in d.items() if isinstance(v, dict)}   # 고침: 글자 칸(written, part, mode) 제외
    av = {**est_only(J('p3_average')), **est_only(J('nl_average')), **est_only(J('anm_average'))}
    ag = {**est_only(J('p3_argsort')), **est_only(J('nl_argsort')), **est_only(J('anm_argsort'))}
    st = {}
    try: st['first_1b'] = triple(json.load(open(RES + '/r38_p3.json'))['main'])
    except Exception: st['first_1b'] = None
    try:
        e = json.load(open(RES + '/r38_p3_explore2.json'))['real_sf_only']['main']; st['first_sf_only'] = dict(est=e['est'], lo=e['ci_sf'][0], hi=e['ci_sf'][1])
    except Exception: st['first_sf_only'] = None
    try:
        e = json.load(open(RES + '/r38_c_nl.json'))['H_main']
        for k in ('lin2', 'linall'): st[f'H_{k}'] = dict(est=e['est'][k], lo=e['ci95'][k][0], hi=e['ci95'][k][1])
    except Exception: st['H_lin2'] = st['H_linall'] = None
    try:
        e = json.load(open(RES + '/r43_sens.json'))['s3']['H']['judge']
        for v in ('anm10', 'anm13', 'anm15', 'gnm'): st[f'H_C2_{v}'] = triple(e[v]['main'])
    except Exception:
        for v in ('anm10', 'anm13', 'anm15', 'gnm'): st.setdefault(f'H_C2_{v}', None)
    rows = {}
    for k in sorted(av):
        a, q, s = triple(av[k]), triple(ag[k]), st.get(k)
        sides = dict(stored=zside(s), argsort=zside(q), average=zside(a))
        rows[k] = dict(n=av[k].get('n'), stored=s, argsort=q, average=a, d_avg_argsort=dmax(a, q), d_avg_stored=dmax(a, s), d_argsort_stored=dmax(q, s), sides=sides,
                       crossing=bool(len({x for x in sides.values() if x is not None}) > 1))
    fx = json.load(open(RES + '/r46_align_fix1.json'))['check']; c2a, c2b = fx['check_2a'], fx['check_2b']
    for c in ('C1', 'C2'):
        reg, q, a = triple(c2a[f'registered_{c}']), triple(c2a[f'argsort_{c}']), triple(c2b[f'avg_{c}'])
        rows[f'C_{c}'] = dict(n=c2a[f'argsort_{c}']['n'], stored=reg, argsort=q, average=a, d_avg_argsort=dmax(a, q), d_avg_stored=dmax(a, reg), d_argsort_stored=dmax(q, reg),
                              sides=dict(stored=zside(reg), argsort=zside(q), average=zside(a)), crossing=bool(len({zside(reg), zside(q), zside(a)}) > 1), source='r46_align_fix1 check')
    mx = lambda key: max((r[key] for r in rows.values() if r[key] is not None), default=None)
    out = dict(written=time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime()), estimates=rows, ties=J('ties_none'),
               max_d_avg_stored=mx('d_avg_stored'), max_d_avg_argsort=mx('d_avg_argsort'), max_d_argsort_stored=mx('d_argsort_stored'),
               crossings=[k for k, r in rows.items() if r['crossing']], fix='addendum10_fix1·fix2: combine 이 추정치 칸만 읽도록 고침, 등록 스크립트를 마지막 __main__ 자리에서 자름',
               rule='등록 r47_ties: 원고 본문은 등록된 argsort 값(저장값)을 그대로 둔다. 방법 절에 평균 순위와의 최대 |차|(max_d_avg_stored)와 동점 순서에 따른 최대 |차|를 "at most …" 로 적는다. '
                    '0 에 대한 구간 위치가 저장값·argsort·평균 순위 사이에서 바뀌는 추정치(crossings)가 있으면 원고의 그 자리에 두 판(저장값과 평균 순위)을 함께 보고한다.')
    json.dump(out, open(RES + '/r47_ties.json', 'x'), ensure_ascii=False, indent=1); log('combine_fix2 끝', out['max_d_avg_stored'], out['crossings'])


if __name__ == '__main__':
    part_combine_fix1(); log('R47F_DONE')
