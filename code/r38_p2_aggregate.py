"""(prereg round36 7c) 1차 P2 집계와 2차 1–4. 규칙은 prereg 7c(2026-10-19 판) 그대로 옮긴다. 바름 점검 5 로 먼저 돌리고 md5 를 동결 점검표에 적은 뒤 67 에 돌린다.
입력: derived/r38/compare/<dom>.json (r38_compare.py), results/r38_dff.json (r38_dff.py), CATH 상과(raw/cath v4.4.0, 없으면 v4.2.0).
P2 (순서 판): 도메인마다 aSAMt 320 K 대 MD 320 K 의 Δμ²/변(dmu) 를 같은 창·같은 Q 의 MD 온도 반응 dmu(T) 와 견준다.
  - T 가 "채워짐" = 그 창·Q 에서 N=100 을 채운 MD 반복이 둘 이상. 320 K 가 안 채워지거나 온도 점이 둘 미만이면 뺀다.
  - 구간 경계: 379 K(+59), 413 K(+93), 450 K(+130). 모델 dmu 보다 반응이 처음 커지는(채워진 가장 낮은) 온도로 구간을 정한다(비단조 규칙).
  - 채워진 온도 어디서도 넘지 못하면: 맨 위가 450 K 면 "> +130", 아니면 맨 위 너머로 오른쪽 검열(가능한 가장 낮은 구간 = 맨 위 다음 경계의 구간).
  - 4 구간 [<+59, [+59,+93), [+93,+130), >+130] 과 3 구간 [<+59, [+59,+93), ≥+93] 을 함께 낸다(450 K 도메인이 절반 미만이면 3 구간이 판정).
  - 검열 합산: 검열 도메인을 가능한 가장 낮은 구간에 둔 판과 가장 높은 구간에 둔 판. 상과 안 아래 중앙값 → 상과 사이 아래 중앙값.
    두 판이 같으면 그 구간, 다르면 합친 구간. 상과 부트스트랩 2000(seed 7)의 95% 구간도 두 판의 합집합.
  - 판정: 구간이 한 구간(또는 합친 구간) 안이면 그 구간, 걸치면 "경계에 걸침".
  - p_T: T ∈ {379, 413, 450} 마다 그 T 가 채워진 도메인 가운데 모델 dmu > 반응 dmu(T) 인 비율(상과마다 비율 → 상과 평균, 상과 부트스트랩 95%).
2차 1: aSAMt 대 MD(320 K) d2n 이 그 도메인 귀무(null320d2n 10 개)의 최댓값 위인 도메인 비율 ≥ 0.70.
2차 2: AlphaFlow 대 MD 320 K dmu 의 상과 중앙값 > R3′ d_ff dmu 의 중앙값(고친 코드, r38_dff.json).
2차 3: 도메인마다 d = √max(dmu,0). d_ff 는 길이가 가장 가까운 R3′ 10 사슬에서 고르게 뽑는다("크기 맞춤", 2026-10-19 이 스크립트로 고정), 1000 회.
  괄호 [|d_AF − d_ff|, d_AF + d_ff], aSAMt d_err 의 상대 위치 p 와 폭 w 의 도메인 중앙값 → 상과 중앙값. 0 ≤ p ≤ 1 이면 "모순되지 않는다",
  거기에 w < d_err 이면 "같은 곳을 가리킨다", p 가 밖이면 "어림 밖".
2차 4: 450 K 기준이 있는 도메인에서 aSAMt 450 대 MD 450 d2n 과 aSAMt 320 대 MD 320 d2n 의 상과 중앙값 비 > 1 (서술).
사용: python r38_p2_aggregate.py <ROOT> <도메인 목록 파일> <출력 꼬리표> [창=full] [Q=0.75]"""
import os, sys, json, glob, time
import numpy as np
ROOT=sys.argv[1]; doms=open(sys.argv[2]).read().split(); TAG=sys.argv[3]; W=sys.argv[4] if len(sys.argv)>4 else 'full'; Q=sys.argv[5] if len(sys.argv)>5 else '0.75'
CMP=ROOT+'/derived/r38/compare'; TS=['348','379','413','450']; BOUND={'379':0,'413':1,'450':2}
sf={}
for f in (ROOT+'/raw/cath/cath-domain-list-v4_2_0.txt', ROOT+'/raw/cath/cath-domain-list.txt'):   # v4.4.0 이 뒤에 와서 우선
    for l in open(f):
        if l.startswith('#') or not l.strip(): continue
        p=l.split(); sf[p[0]]='.'.join(p[1:5])
rng=np.random.default_rng(7)
def lower_median(v):
    v=sorted(v); return v[(len(v)-1)//2]
def sf_median(rows, key):
    g={}
    for r in rows: g.setdefault(r['sf'],[]).append(r[key])
    return lower_median([lower_median(v) for v in g.values()]), g
def boot_ord(rows, key, B=2000):
    g={}
    for r in rows: g.setdefault(r['sf'],[]).append(r[key])
    S=list(g); out=[]
    for _ in range(B):
        pick=[S[i] for i in rng.integers(0,len(S),len(S))]; out.append(lower_median([lower_median(g[s]) for s in pick]))
    return int(np.percentile(out,2.5,method='lower')), int(np.percentile(out,97.5,method='higher'))
def sf_mean_frac(rows,key):
    g={}
    for r in rows: g.setdefault(r['sf'],[]).append(r[key])
    S=list(g); est=float(np.mean([np.mean(g[s]) for s in S])); bs=[]
    for _ in range(2000):
        pick=[S[i] for i in rng.integers(0,len(S),len(S))]; bs.append(np.mean([np.mean(g[s]) for s in pick]))
    return dict(n_dom=len(rows),n_sf=len(S),frac=round(est,3),ci=[round(float(np.percentile(bs,2.5)),3),round(float(np.percentile(bs,97.5)),3)])
rows=[]; skipped={}
for d in doms:
    f=f'{CMP}/{d}.json'
    if not os.path.exists(f): skipped[d]='비교 없음'; continue
    j=json.load(open(f)); c=j['cells']; nr=j['nruns']
    filled=lambda T: nr.get(f'MD{T}|{W}|{Q}',0)>=2
    mk=f'aSAMt320-MD320|{W}|{Q}'
    if not filled('320') or mk not in c: skipped[d]='320 K 기준 없음 또는 aSAMt 칸 없음'; continue
    ft=[T for T in TS if filled(T) and f'MD{T}-MD320|{W}|{Q}' in c]
    if len(ft)<1: skipped[d]='온도 점 둘 미만'; continue
    m=c[mk]['dmu']; resp={T:c[f'MD{T}-MD320|{W}|{Q}']['dmu'] for T in ft}; top=max(ft,key=int)
    lo=hi=None
    for T in ('379','413','450'):
        if T in resp and resp[T]>m: lo=hi=BOUND[T]; break
    cens=False
    if lo is None:
        if top=='450': lo=hi=3
        else:
            cens=True; lo={'348':0,'379':1,'413':2}[top]; hi=3
            if top=='348' and resp['348']>m: lo=hi=0; cens=False
    r=dict(dom=d,sf=sf.get(d,'?'),L=j['L'],top=top,model_dmu=m,resp=resp,lo4=lo,hi4=hi,lo3=min(lo,2),hi3=min(hi,2),cens=cens)
    for T in ('379','413','450'):
        if T in resp: r[f'exceed{T}']=int(m>resp[T])
    nl=c.get(f'null320d2n|{W}|{Q}'); r['above_null']=int(c[mk]['d2n']>max(nl)) if nl else None
    afk=f'AlphaFlow-MD320|{W}|{Q}'; r['af_dmu']=c[afk]['dmu'] if afk in c else None
    k450=f'aSAMt450-MD450|{W}|{Q}'; r['d2n320']=c[mk]['d2n']; r['d2n450']=c[k450]['d2n'] if k450 in c else None
    rows.append(r)
res=dict(written=time.strftime('%F %T'),tag=TAG,window=W,Q=Q,n_input=len(doms),n_used=len(rows),skipped=skipped,
         top_counts={T:sum(r['top']==T for r in rows) for T in TS},n_450=sum(r['top']=='450' for r in rows))
names4=['<+59','[+59,+93)','[+93,+130)','>+130']; names3=['<+59','[+59,+93)','≥+93']
for nb,names in ((4,names4),(3,names3)):
    if not rows: break
    for v in ('lo','hi'):
        for r in rows: r[f'_{v}']=r[f'{v}{nb}']
    mlo,_=sf_median(rows,'_lo'); mhi,_=sf_median(rows,'_hi'); blo=boot_ord(rows,'_lo'); bhi=boot_ord(rows,'_hi')
    a,b=min(blo[0],bhi[0]),max(blo[1],bhi[1])
    point=names[mlo] if mlo==mhi else f'{names[mlo]}–{names[mhi]} (검열 구역)'
    verdict=(names[a] if a==b else (point if (mlo!=mhi and mlo<=a and b<=mhi) else '경계에 걸침'))   # 구간이 한 구간 또는 검열 합산의 합친 구간 안이면 그것
    res[f'P2_{nb}bins']=dict(median_lo=names[mlo],median_hi=names[mhi],point=point,ci95=[names[a],names[b]],verdict=verdict,
                            n_censored=sum(r['cens'] for r in rows),bins_lo=[sum(r['_lo']==i for r in rows) for i in range(nb)],bins_hi=[sum(r['_hi']==i for r in rows) for i in range(nb)])
res['P2_verdict_rule']='3 구간' if res['n_450']*2<len(rows) else '4 구간'
res['p_T']={T:sf_mean_frac([r for r in rows if f'exceed{T}' in r],f'exceed{T}') for T in ('379','413','450') if any(f'exceed{T}' in r for r in rows)}
an=[r for r in rows if r['above_null'] is not None]
res['sec1_above_null']=dict(n=len(an),frac=round(float(np.mean([r['above_null'] for r in an])),3) if an else None,pass_=bool(an and np.mean([r['above_null'] for r in an])>=0.70))
dff=json.load(open(ROOT+'/results/r38_dff.json')) if os.path.exists(ROOT+'/results/r38_dff.json') else None
af=[r for r in rows if r['af_dmu'] is not None]
if af and dff and dff.get('rows'):
    ffm=float(np.median([x['dmu'] for x in dff['rows']])); afm,_=sf_median([dict(sf=r['sf'],v=r['af_dmu']) for r in af],'v')
    res['sec2']=dict(af_dmu_sf_median=round(afm,5),dff_dmu_median=round(ffm,5),n_ff=len(dff['rows']),pass_=bool(afm>ffm))
    FF=sorted(dff['rows'],key=lambda x:x['L']); pw=[]
    for r in af:
        near=sorted(FF,key=lambda x:abs(x['L']-r['L']))[:10]; dd=np.sqrt(np.maximum([x['dmu'] for x in near],0))
        dff_s=dd[rng.integers(0,len(dd),1000)]; dobs=np.sqrt(max(r['af_dmu'],0)); derr=np.sqrt(max(r['model_dmu'],0))
        lo_=np.abs(dobs-dff_s); hi_=dobs+dff_s; w=hi_-lo_; p=(derr-lo_)/np.where(w>0,w,np.nan)
        pw.append(dict(sf=r['sf'],p=float(np.nanmedian(p)),w=float(np.median(w)),derr=float(derr)))
    pm,_=sf_median(pw,'p'); wm,_=sf_median(pw,'w'); dm,_=sf_median(pw,'derr')
    res['sec3']=dict(p_sf_median=round(pm,3),w_sf_median=round(wm,4),derr_sf_median=round(dm,4),
                     verdict=('같은 곳을 가리킨다' if (0<=pm<=1 and wm<dm) else ('모순되지 않는다' if 0<=pm<=1 else '어림 밖(모델에 따라 다르다)')))
h=[r for r in rows if r['d2n450'] is not None]
if h:
    a450,_=sf_median([dict(sf=r['sf'],v=r['d2n450']) for r in h],'v'); a320,_=sf_median([dict(sf=r['sf'],v=r['d2n320']) for r in h],'v')
    res['sec4']=dict(n=len(h),ratio_450_320=round(a450/a320,3) if a320 else None)
res['rows']=[{k:v for k,v in r.items() if not k.startswith('_')} for r in rows]
P=f'{ROOT}/results/r38_p2_{TAG}.json'; json.dump(res,open(P,'w'),indent=1,ensure_ascii=False)
import shutil; b=ROOT+'/results/_backup/'+time.strftime('%Y-%m-%d'); os.makedirs(b,exist_ok=True); shutil.copy(P,b+'/'+time.strftime('%H%M%S_')+os.path.basename(P))
print(json.dumps({k:v for k,v in res.items() if k!='rows'},ensure_ascii=False,indent=1))
