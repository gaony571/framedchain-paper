# Interactive checks: code and printed output

Seven results were computed in interactive Colab cells, or by a script whose printed output was not saved to a file:
the comparison of 12 test pairs with TM-align (S2 Text 8.1; section 1), the reversed chains (S2 Text 8.5; section 2),
the equipment check before the computational-cost measurement (S2 Text 7; section 3), the decomposition of the
unshared part between mean points (section 4), where correspondence shifts occur (section 5), the chain-break check
(section 6) and the reproduction of the as-run code by framedchain 0.4.1 (section 7). The code and printed output
survive in the client-side log of the working session through which the cells were
run, which records every cell run with the UTC time, by the local clock, at which it was requested and at which
its output returned.
This file copies them from that log verbatim; only the parts noted were left out. Comments and printed labels in
the copied code and output are in Korean, the working language of the project; the labels needed to read the results
are glossed where they occur.

Names used by the cells but defined earlier in the same runtime: `ROOT = '/content/drive/MyDrive/atlas_framed_curves'`;
`z = np.load(ROOT + '/derived/xprot_structs.npz')` (N, CA, C and O of the 379 non-redundant ATLAS chains);
`names` = its chain names, sorted; `names_o` = its chain names in file order; `s0` = `framespace_wip/stage0`;
`fs` = the framespace package (`framespace_pkg`); `tma` and `tma_orig` = `tmtools.tm_align` (at 15:36 `tma` was
redefined to raise an error so that later cells could not call it by mistake, and `tma_orig` kept the original).
In sections 4–6, `R` is the same project root and `J(name)` reads `results/<name>.json`.

## 1. Twelve test pairs against TM-align (16 October 2026)

**Module version.** The runs used here (15:32 and 16:18 UTC) imported `framespace_wip/stage0/rulers.py` with md5
`9156347202b9e91ca2227c654b5fc7b4` (written to Drive at 15:30:45 UTC). The deposited `rulers.py`
(`eee47ef47ac7f783201202d081c77990`) differs from it in two later edits only: the relative tolerance of
`triangle_rate` (16:19 UTC) and the removal of zero terms in `wl_align` (17:28 UTC). `tm_align`, `tm_cost`,
`d0_tm`, `start_poses` and `sub_fragment_poses` are byte-identical in the two versions.

**Pair selection.** The 12 pairs were drawn with `np.random.default_rng(0)` from the 1,379 rows of
`results/v_residue_unit_graded.tsv` whose two chains were in `xprot_structs.npz` (1,371 non-redundant pairs and
8 same-superfamily pairs); all 12 drawn pairs are among the 1,371.

**First run (superseded).** The first run used the first, torch-based version of the module (md5
`15f125a50824afa155152f1ef3300aa5`, 0.4–84 s per pair). It was rewritten with numba, the sub-fragment search and
external starting poses before the runs below; the S2 table does not use this run.

Cell `S_otyYtETrub`, run requested 15:24:58, output returned 15:27:17 UTC.

```python
for pth in [ROOT + '/framespace_pkg', ROOT + '/framespace_wip']:
    if pth not in sys.path: sys.path.insert(0, pth)
import importlib, framespace as fs, stage0 as s0; importlib.reload(s0.rulers); importlib.reload(s0)
from tmtools import tm_align as tma
import time
CA = {n: z[n + '|CA'][0] if z[n + '|CA'].ndim == 3 else z[n + '|CA'] for n in names}
g = pd.read_csv(ROOT + '/results/v_residue_unit_graded.tsv', sep='\t')
gg = g[g.A.isin(CA) & g.B.isin(CA)]; print(len(gg), gg.src.value_counts().to_dict())
rng = np.random.default_rng(0); smp = gg.iloc[rng.choice(len(gg), 12, replace=False)]
rows = []
for _, rw in smp.iterrows():
    A, B = CA[rw.A], CA[rw.B]; seqA, seqB = 'A' * len(A), 'A' * len(B)
    t0 = time.time(); res = tma(A, B, seqA, seqB); t1 = time.time()
    d0 = s0.d0_tm(len(B)); r = s0.tm_align(A, B, d0); t2 = time.time()
    rows.append((rw.A, rw.B, len(A), len(B), round(res.tm_norm_chain2, 4), round(r['S'] / len(B), 4), round(t1 - t0, 2), round(t2 - t1, 2)))
print(pd.DataFrame(rows, columns=['A', 'B', 'n', 'm', 'TMalign_normB', 'ours_normB', 's_tmalign', 's_ours']).to_string())
```

Output:

```
1379 {'benchmark': 1371, 'homolog': 8}
         A       B    n    m  TMalign_normB  ours_normB  s_tmalign  s_ours
0   4a56_A  1w53_A   93   84         0.4968      0.3690       0.04   13.81
1   4obi_A  1n8u_A   87  100         0.2604      0.2760       0.03    2.21
2   3pf6_A  3t7z_A   60  119         0.3188      0.3240       0.03    2.15
3   3v3l_A  3txs_C   83   91         0.2763      0.2338       0.02    1.08
4   3zdm_A  4npd_A   69   58         0.4583      0.3820       0.01    0.49
5   1ojh_C  3onj_A   53   97         0.3583      0.2954       0.01    0.36
6   1ptq_A  4e1p_A   50   55         0.4231      0.4296       0.00    0.76
7   4exo_A  5hub_A  146   79         0.5473      0.5216       0.03   12.62
8   3txs_C  2p58_A   91   54         0.6949      0.8163       0.01    3.88
9   2p58_A  5cof_A   54  172         0.2366      0.2428       0.02    1.39
10  2pmk_A  3vnb_A  243  149         0.3877      0.3994       0.09   84.20
11  2q9r_A  1ail_A  195   70         0.5638      0.5696       0.02    4.85
```

**Second run: TM-align's alignment scored by our definition, and our optimization with and without TM-align's
superposition as a start.** Columns of the S2 table: TM-align = `TMalign`; TM-align's alignment, rescored =
`TMalign_aln_ours_def`; ours without TM-align start = `ours`; ours with TM-align start = `ours+TMalign_pose`.

Cell `S_otyYtETrub`, run requested 15:30:57, output returned 15:32:33 UTC.

```python
importlib.reload(s0.rulers); importlib.reload(s0)
def aln_pairs(res):
    i = j = 0; ii, jj = [], []
    for a, b in zip(res.seqxA, res.seqyA):
        if a != '-' and b != '-': ii.append(i); jj.append(j)
        i += a != '-'; j += b != '-'
    return np.array(ii), np.array(jj)
def tmalign_pose(res): return res.u.T, -res.u.T @ res.t
_ = s0.tm_align(CA[smp.iloc[0].A][:30], CA[smp.iloc[0].B][:30], 2.0)   # numba 컴파일
rows = []
for _, rw in smp.iterrows():
    A, B = CA[rw.A], CA[rw.B]
    t0 = time.time(); res = tma(A, B, 'A' * len(A), 'A' * len(B)); t1 = time.time()
    d0 = s0.d0_tm(len(B)); r = s0.tm_align(A, B, d0); t2 = time.time()
    U, t = tmalign_pose(res); ii, jj = aln_pairs(res)
    M = np.zeros((len(A), len(B)), bool); M[ii, jj] = True
    S_tmaln = 0.5 * (len(A) + len(B)) - s0.tm_cost(A, B, M, U, t, d0)     # TM-align 의 대응·자세를 우리 정의로 잰 값
    r2 = s0.tm_align(A, B, d0, extra_poses=[(U, t)])
    rows.append((rw.A, rw.B, len(A), len(B), round(res.tm_norm_chain2, 4), round(S_tmaln / len(B), 4), round(r['S'] / len(B), 4), round(r2['S'] / len(B), 4), round(t1 - t0, 2), round(t2 - t1, 2)))
print(pd.DataFrame(rows, columns=['A', 'B', 'n', 'm', 'TMalign', 'TMalign_aln_ours_def', 'ours', 'ours+TMalign_pose', 's_tmalign', 's_ours']).to_string())
```

Output:

```
         A       B    n    m  TMalign  TMalign_aln_ours_def    ours  ours+TMalign_pose  s_tmalign  s_ours
0   4a56_A  1w53_A   93   84   0.4968                0.4968  0.5021             0.5021       0.02    0.33
1   4obi_A  1n8u_A   87  100   0.2604                0.2604  0.2706             0.2706       0.03    0.18
2   3pf6_A  3t7z_A   60  119   0.3188                0.3188  0.3252             0.3252       0.01    0.08
3   3v3l_A  3txs_C   83   91   0.2763                0.2763  0.2408             0.2863       0.02    0.15
4   3zdm_A  4npd_A   69   58   0.4583                0.4583  0.4640             0.4640       0.01    0.10
5   1ojh_C  3onj_A   53   97   0.3583                0.3583  0.3685             0.3685       0.01    0.14
6   1ptq_A  4e1p_A   50   55   0.4231                0.4231  0.4296             0.4296       0.00    4.59
7   4exo_A  5hub_A  146   79   0.5473                0.5473  0.5216             0.5561       0.12    0.44
8   3txs_C  2p58_A   91   54   0.6949                0.6949  0.8163             0.8163       0.05    0.21
9   2p58_A  5cof_A   54  172   0.2366                0.2366  0.2501             0.2501       0.08    0.45
10  2pmk_A  3vnb_A  243  149   0.3877                0.3877  0.4050             0.4051       0.26    0.32
11  2q9r_A  1ail_A  195   70   0.5638                0.5638  0.5197             0.5696       0.02    0.11
```

**C++ TM-align.** Downloaded and compiled for rescoring with `-I`.

Cell `AJGEYYNzgKHE`, run requested 16:17:04, output returned 16:17:08 UTC.

```python
os.makedirs('/content/tmalign_src', exist_ok=True)
r_ = subprocess.run(['wget', '-q', '-O', '/content/tmalign_src/TMalign.cpp', 'https://zhanggroup.org/TM-align/TMalign.cpp'], capture_output=True, text=True)
print(r_.returncode, os.path.getsize('/content/tmalign_src/TMalign.cpp') if os.path.exists('/content/tmalign_src/TMalign.cpp') else None)
print(open('/content/tmalign_src/TMalign.cpp').read()[:600])
tm_build = subprocess.Popen(['g++', '-static', '-O3', '-ffast-math', '-lm', '-o', '/content/tmalign_src/TMalign', '/content/tmalign_src/TMalign.cpp'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
```

Output (first line; the rest is the head of the source file):

```
0 182099
```

Cell `PV7sFW0NgPwF`, run requested 16:18:14, output returned 16:18:18 UTC.

```python
print('build', tm_build.poll(), (tm_build.stdout.read() if tm_build.poll() is not None else '')[-500:], os.path.exists('/content/tmalign_src/TMalign'))
TMB = '/content/tmalign_src/TMalign'; WD = '/content/tm_ind'; os.makedirs(WD, exist_ok=True)
def write_ca_pdb(X_, path_):
    with open(path_, 'w') as f_:
        for i_, (x_, y_, z_) in enumerate(X_):
            f_.write('ATOM  %5d  CA  ALA A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C\n' % (i_ + 1, i_ + 1, x_, y_, z_))
        f_.write('TER\nEND\n')
def aln_strings(ii_, jj_, n_, m_):
    sa, sb = [], []; pi_, pj_ = 0, 0
    for i_, j_ in zip(ii_, jj_):
        while pi_ < i_: sa.append('A'); sb.append('-'); pi_ += 1
        while pj_ < j_: sa.append('-'); sb.append('A'); pj_ += 1
        sa.append('A'); sb.append('A'); pi_ += 1; pj_ += 1
    while pi_ < n_: sa.append('A'); sb.append('-'); pi_ += 1
    while pj_ < m_: sa.append('-'); sb.append('A'); pj_ += 1
    return ''.join(sa), ''.join(sb)
def blocks_gaps(ii_, jj_):
    """연속 짝 조각 수(i, j 모두 1 씩 늘면 같은 조각)와 안쪽 틈(조각 사이) 수."""
    if len(ii_) == 0: return 0, 0
    br = np.sum((np.diff(ii_) != 1) | (np.diff(jj_) != 1)); return int(br + 1), int(br)
def tmalign_cli(pa_, pb_, extra_):
    out_ = subprocess.run([TMB, pa_, pb_] + extra_, capture_output=True, text=True).stdout
    tm2 = [l_ for l_ in out_.splitlines() if l_.startswith('TM-score=') and 'Chain_2' in l_]
    return float(tm2[0].split()[1]) if tm2 else None, out_
```

Output:

```
build 0  True
```

**Third run: rescoring our alignment with C++ TM-align (`-I`), aligned pairs and fragments.** Columns of the S2
table: ours, rescored by TM-align (-I) = `ours_rescored_by_TMalign_I` (rounded to four decimals); aligned pairs =
`k_ours`/`k_tma`; fragments (S2 text) = `blocks_ours`/`blocks_tma`. `TMalign_cli` (compiled C++) agrees with
`tmtools` to four decimals in all 12 pairs.

Cell `2cRORPL3gf6o`, run requested 16:18:34, output returned 16:18:38 UTC.

```python
rows_ind = []
for _, rw in smp.iterrows():
    A_, B_ = CA[rw.A], CA[rw.B]; n_, m_ = len(A_), len(B_)
    pa_, pb_ = f'{WD}/{rw.A}.pdb', f'{WD}/{rw.B}.pdb'; write_ca_pdb(A_, pa_); write_ca_pdb(B_, pb_)
    res_ = tma_orig(A_, B_, 'A' * n_, 'A' * m_); ext_ = [tmalign_pose(res_)]
    d0_ = s0.d0_tm(m_); r_ = s0.tm_align(A_, B_, d0_, extra_poses=ext_)
    oi_, oj_ = np.nonzero(r_['match']); sa_, sb_ = aln_strings(oi_, oj_, n_, m_)
    open(f'{WD}/aln.fa', 'w').write(f'>A\n{sa_}\n>B\n{sb_}\n')
    tm_I, _o = tmalign_cli(pa_, pb_, ['-I', f'{WD}/aln.fa'])
    tm_cli, _o2 = tmalign_cli(pa_, pb_, [])
    ti_, tj_ = aln_pairs(res_)
    rows_ind.append(dict(A=rw.A, B=rw.B, n=n_, m=m_, tmtools=round(res_.tm_norm_chain2, 4), TMalign_cli=tm_cli, ours=round(r_['S'] / m_, 4), ours_rescored_by_TMalign_I=tm_I,
                         k_ours=len(oi_), k_tma=len(ti_), blocks_ours=blocks_gaps(oi_, oj_)[0], blocks_tma=blocks_gaps(ti_, tj_)[0]))
ind_df = pd.DataFrame(rows_ind); print(ind_df.to_string()); print(_o[-900:])
```

Output (table only; the TM-align report printed after it for the last pair is left out):

```
         A       B    n    m  tmtools  TMalign_cli    ours  ours_rescored_by_TMalign_I  k_ours  k_tma  blocks_ours  blocks_tma
0   4a56_A  1w53_A   93   84   0.4968      0.49684  0.5021                     0.50094      82     59            8           8
1   4obi_A  1n8u_A   87  100   0.2604      0.26041  0.2706                     0.27033      49     46           14          14
2   3pf6_A  3t7z_A   60  119   0.3188      0.31877  0.3252                     0.32483      60     53           13          10
3   3v3l_A  3txs_C   83   91   0.2763      0.27634  0.2863                     0.28527      66     40           11          10
4   3zdm_A  4npd_A   69   58   0.4583      0.45829  0.4640                     0.46239      42     35            6           5
5   1ojh_C  3onj_A   53   97   0.3583      0.35831  0.3685                     0.36833      50     43           10           5
6   1ptq_A  4e1p_A   50   55   0.4231      0.42310  0.4296                     0.42756      39     34            5           4
7   4exo_A  5hub_A  146   79   0.5473      0.54732  0.5561                     0.55528      72     69           16          11
8   3txs_C  2p58_A   91   54   0.6949      0.69486  0.8163                     0.81534      54     53            4           6
9   2p58_A  5cof_A   54  172   0.2366      0.23660  0.2501                     0.24999      54     49            8           9
10  2pmk_A  3vnb_A  243  149   0.3877      0.38772  0.4051                     0.40468     127    104           20          18
11  2q9r_A  1ail_A  195   70   0.5638      0.56383  0.5696                     0.56867      62     59           11           9
```

## 2. Reversed chains (16 October 2026)

The eight chains are the first eight of `xprot_structs.npz` in file order. D_TM at fixed d₀ = 4 Å:

Cell `oC2Qk_XKhfGf`, run requested 16:22:54, output returned 16:26:14 UTC.

```python
rev_rows = []
for n_ in names_o[:8]:
    X_ = CA[n_]; Xr_ = X_[::-1].copy(); L_ = len(X_)
    r_ = s0.tm_align(X_, Xr_, 4.0); rself_ = s0.tm_align(X_, X_, 4.0)
    rev_rows.append((n_, L_, round(rself_['D2'], 6), round(r_['D2'], 2), round(np.sqrt(r_['D2']), 2), round(np.sqrt(r_['D2'] / L_), 3), round(r_['S'] / L_, 3)))
print(pd.DataFrame(rev_rows, columns=['name', 'L', 'D2_self', 'D2_rev', 'D_rev', 'D_rev/sqrt(L)', 'S_rev/L']).to_string())
```

Output:

```
     name    L  D2_self  D2_rev  D_rev  D_rev/sqrt(L)  S_rev/L
0  3mx7_A   90      0.0   38.67   6.22          0.655    0.570
1  4obi_A   87      0.0   46.41   6.81          0.730    0.467
2  4ftf_A  109      0.0   56.55   7.52          0.720    0.481
3  4exo_A  146      0.0  101.62  10.08          0.834    0.304
4  3kz7_A  119      0.0   81.23   9.01          0.826    0.317
5  1huf_A  123      0.0   79.58   8.92          0.804    0.353
6  3dso_A   66      0.0   28.63   5.35          0.659    0.566
7  2vri_A  173      0.0  103.60  10.18          0.774    0.401
```

The framed-chain distance d. The cell first set up a run that is not part of this check (the d₀-sensitivity script
of the factorial comparison); only the reversed-chain part is copied. In the code comment, n is the number of
bonds, so x'_i = x_{n−1−i} diag(−1, 1, −1) with 0-based indices is x̄_i = x_{n+1−i} diag(−1, 1, −1),
i = 1, …, n, with 1-based indices, as in S2 Text. The module version was `f728959516a789a936722cb1640ac0b1`
(the 16:19 edit of `triangle_rate` only).

Cell `JBTVYnWkjiau`, run requested 16:32:07, output returned 16:41:04 UTC.

```python
# (2) 뒤집은 사슬: 프레임 기본 판. 뒤집은 사슬의 변 i' 는 원래 변 n−1−i 를 거꾸로 걷는 것: T → −T, U 그대로(같은 펩타이드의 C=O), V → −V. 곧 x'_i = x_{n−1−i} diag(−1, 1, −1)
def bb_(n_, a_):
    v_ = z[n_ + '|' + a_]; return v_[0] if v_.ndim == 3 else v_
rev2 = []
for n_ in names_o[:8]:
    x_ = fs.coordinates_from_backbone(bb_(n_, 'CA'), bb_(n_, 'C'), bb_(n_, 'O'))
    xr_ = np.ascontiguousarray(x_[::-1] @ np.diag([-1.0, 1.0, -1.0]))
    d2s_, _ = fs.sq_distance(torch.as_tensor(x_), torch.as_tensor(x_.copy())); d2r_, _ = fs.sq_distance(torch.as_tensor(x_), torch.as_tensor(xr_))
    ne_ = len(x_); rev2.append((n_, ne_, round(float(d2s_), 6), round(float(d2r_), 2), round(float(fs.length_scaled(float(d2r_) ** 0.5, ne_, ne_)), 3), bool(fs.is_real(xr_).all())))
print(pd.DataFrame(rev2, columns=['name', 'edges', 'd2_self', 'd2_rev', 'd_rev/sqrt(0.75(n+m))', 'rev_is_real']).to_string())
```

Output:

```
     name  edges  d2_self  d2_rev  d_rev/sqrt(0.75(n+m))  rev_is_real
0  3mx7_A     89      0.0   51.33                  0.620         True
1  4obi_A     86      0.0   52.81                  0.640         True
2  4ftf_A    108      0.0   89.18                  0.742         True
3  4exo_A    145      0.0  131.10                  0.776         True
4  3kz7_A    118      0.0   76.53                  0.658         True
5  1huf_A    122      0.0   99.27                  0.737         True
6  3dso_A     65      0.0   36.59                  0.613         True
7  2vri_A    172      0.0  141.23                  0.740         True
```

## 3. Equipment of the computational-cost measurement (22 October 2026)

The runtime started at about 07:28 UTC (`up 20 min` at 07:48). The checks below ran in it before the inputs
were prepared; the measurement inputs were then written to the local directory `/content/ct`, which existed
until the end of the measurement, so the measurement ran in the same runtime.

Cell `lIYdn1woOS1n`, run requested 07:30:41, output returned 07:31:24 UTC.

```python
import subprocess, os
print(subprocess.run("date -u; nvidia-smi --query-gpu=name,driver_version,memory.total,memory.used,utilization.gpu --format=csv,noheader; lscpu | grep -E 'Model name|^CPU\\(s\\)|Thread\\(s\\)|Core\\(s\\)|Socket'; free -g | head -2; nproc; ps -eo pid,etime,args --sort=-%cpu | head -8; ls /content", shell=True, capture_output=True, text=True).stdout)
print('drive', os.path.exists('/content/drive/MyDrive'))
```

Output (the process list is left out; it showed only the Colab kernel and server processes):

```
Thu Oct 22 07:31:25 AM UTC 2026
NVIDIA A100-SXM4-40GB, 580.82.07, 40960 MiB, 0 MiB, 0 %
CPU(s):                                  12
Model name:                              Intel(R) Xeon(R) CPU @ 2.20GHz
Thread(s) per core:                      2
Core(s) per socket:                      6
Socket(s):                               1
               total        used        free      shared  buff/cache   available
Mem:              83           1          78           0           3          81
12
sample_data

drive False
```

Cell `lIYdn1woOS1n`, run requested 07:31:30, output returned 07:32:14 UTC.

```python
import os, json, subprocess, time, shutil
from google.colab import drive
if not os.path.exists('/content/drive/MyDrive'): drive.mount('/content/drive')
R = "/content/drive/MyDrive/atlas_framed_curves"; L = R + "/results/logs"; C = R + "/code"
sh = lambda c: subprocess.run(c, shell=True, capture_output=True).stdout.decode('utf-8', 'replace')
print(sh(f"ls {R}/derived/r38/static | head; wc -l {R}/derived/r38/static/bnb.jsonl; python3 -c 'import torch, numpy; print(torch.__version__, numpy.__version__)'; pip list 2>/dev/null | grep -i -E 'tmtools|h5py|mdtraj'"))
```

Output (torch and numpy versions on the second-to-last line):

```
Mounted at /content/drive
align.jsonl
bnb.jsonl
medoid
table.jsonl
500 /content/drive/MyDrive/atlas_framed_curves/derived/r38/static/bnb.jsonl
2.11.0+cu130 2.1.3
h5py                                  3.14.0
```

Cell `lIYdn1woOS1n`, run requested 07:33:52, output returned 07:34:17 UTC.

```python
T = '/content/ct'; os.makedirs(T, exist_ok=True)
print(sh(f"""pip install -q tmtools 2>&1 | tail -1; python3 -c 'import tmtools; print("tmtools", tmtools.__version__ if hasattr(tmtools,"__version__") else "ok")'
cd {T} && mkdir -p dl && cd dl && curl -sSL -o foldseek-linux-avx2.tar.gz https://github.com/steineggerlab/foldseek/releases/download/10-941cd33/foldseek-linux-avx2.tar.gz && ls -l foldseek-linux-avx2.tar.gz && sha256sum foldseek-linux-avx2.tar.gz && tar xzf foldseek-linux-avx2.tar.gz -C {T} && {T}/foldseek/bin/foldseek version
cd {T} && git clone -q https://github.com/pylelab/USalign.git && cd USalign && git rev-parse HEAD && grep -n 'USalign:' -A2 Makefile | head -5"""))
```

Output (tmtools version, Foldseek archive checksum and commit, US-align commit):

```
     ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 3.1/3.1 MB 95.8 MB/s eta 0:00:00
tmtools 0.3.0
-rw-r--r-- 1 root root 11941946 Oct 22 07:34 foldseek-linux-avx2.tar.gz
e6a3e3860e76b0b7b858736f627fe680d4d7f0fc7298301d8a92b9f711d74ac7  foldseek-linux-avx2.tar.gz
941cd33ff0771cd2e3f144e3293e22a2b87e9fda
23f14ffcf8647f2dc1c0e01043db3a7f3a661c28
12:USalign: USalign.cpp SOIalign.h MMalign.h param_set.h basic_fun.h Kabsch.h NW.h TMalign.h pstream.h se.h NWalign.h BLOSUM.h flexalign.h
13-	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}
14-
```

**Timeline.** Input preparation ended at 08:07:01 and the measurement queue (`r41_compute.py`) ran from 08:07 to
08:42:28 UTC. During that time the notebook ran only status checks (`sleep`, `tail` of the logs, `ls`) and, at
08:09, a read of one rebuilt PDB file to check bond lengths.

Cell `lIYdn1woOS1n`, run requested 08:02:05, output returned 08:09:06 UTC.

```python
time.sleep(420); print(sh("date -u +%T; ls /content/ct/x | wc -l; ls /content/ct/pdb | wc -l; grep -v 'I0000\\|WARNING' /content/ct/prep.log | tail -3 | cut -c1-600; tail -3 /content/ct/run.log 2>/dev/null | cut -c1-200"))
```

Output:

```
08:09:06
570
570
08:02:13 prep 400 570 583 s
08:04:17 prep 500 570 707 s
08:07:01 PREP_DONE {'n_dom': 570, 'seq_mismatch': [], 'no_medoid': [], 'seq_source': {'p3_dom': 479, 'h5': 91}, 'template_edges': 3127, 'template_domains': 30, 'n_pairs': 500, 'ava200': 200, 'ava100': 100, 'ens': [{'dom': '2gauA02', 'ne': 80, 'target': 80}, {'dom': '2fzfA01', 'ne': 149, 'target': 150}, {'dom': '2fzfA01', 'ne': 149, 'target': 250}], 'sec': 871.1}
08:07:45 cpu 200 30 s
08:08:00 cpu 300 45 s
08:08:16 cpu 400 61 s
```

Cell `lIYdn1woOS1n`, run requested 08:39:22, output returned 08:47:22 UTC.

```python
time.sleep(480); print(sh("date -u +%T; grep -v 'I0000\\|WARNING' /content/ct/run.log | tail -6 | cut -c1-200; grep -c CT_ALL_DONE " + R + "/results/logs/r41_compute.log 2>/dev/null"))
```

Output:

```
08:47:23
08:42:24 ens 2fzfA01 149 1 1.09 s
08:42:25 ens 2fzfA01 149 2 1.09 s
08:42:26 ens 2fzfA01 149 3 1.09 s
08:42:26 CT_ENS_DONE
08:42:28 부분 시작 sum
08:42:28 CT_SUM_DONE
1
```

## 4. Decomposition of the unshared part between mean points (21 October 2026)

Manuscript, sequence-identity section: about 89% of the unshared part of the squared distance between mean points is
staggered, one-sided loss of correspondence inside the aligned region, about 7% terminal overhangs and about 4.5%
bonds at the same position turned beyond 120° (a descriptive decomposition, in no record). Inputs:
`derived/r38/c/pibar.jsonl` and `derived/r38/p3/pibar.jsonl` (alignment of the two mean points of each pair, with its
correspondence `best` and score `S`), `derived/r38/c/dom/<domain>.npz` and `derived/r38/p3/dom/<domain>.npz` (mean
point `m`, mean length ratio `sbar`) and `derived/r41/tm_med.jsonl` (TM-score of TM-align on medoids). In the output, C is
the confirmatory set and H the first homolog study; "끝 겹침" is terminal overhang, "구간 초과" the unpaired excess
inside the aligned region, "경첩(K≤0)" bonds facing each other at the same position with a non-positive pair score
(turned beyond 120°) and "밀림(K>0)" bonds facing each other with a positive score. The masses are summed over pairs
after dividing each pair's terms by √n̄. The first cell checks the decomposition on 200 pairs (the rebuilt S agrees with
the stored S, and the parts sum to the unshared mass); the second gives the shares.

Cell `HpplqbJChcyO`, run requested 2026-10-21 09:17:04, output returned 09:17:08 UTC.

```python
import numpy as np, json, sys
sys.path[:0] = [R + '/framespace_pkg']
CC = f'{R}/derived/r38/c'; PB = [json.loads(l) for l in open(CC + '/pibar.jsonl')]
rg = np.random.default_rng(1); pick = rg.choice(len(PB), 200, replace=False)
def proc(M):  # argmax_U tr(U M) (fs.align 규약과 같은 M = Σ y xᵀ)
    W, S_, Vt = np.linalg.svd(M); d = np.sign(np.linalg.det(W @ Vt)); return (Vt.T @ np.diag([1, 1, d]) @ W.T)
def ld(d):
    z = np.load(f'{CC}/dom/{d}.npz'); return z['m'].astype(np.float64), z['sbar'].astype(np.float64)
rel, frac, vneg = [], [], []
for q in pick:
    p = PB[q]; a, b = p['ab']; mA, sA = ld(a); mB, sB = ld(b); pr = np.array(p['best'], int)
    M = np.einsum('kxz,kyz->xy', mB[pr[:, 1]], mA[pr[:, 0]]); U = proc(M)
    K = lambda i, j: float((mA[i] * (U @ mB[j])).sum())
    S2 = sum(K(i, j) for i, j in pr); rel.append(abs(S2 - p['S']) / p['S'])
    nA, nB = (mA ** 2).sum((1, 2)), (mB ** 2).sum((1, 2)); ins = hin = slp = 0.0
    o = np.argsort(pr[:, 0]); pr = pr[o]
    ins += nA[:pr[0, 0]].sum() + nA[pr[-1, 0] + 1:].sum() + nB[:pr[0, 1]].sum() + nB[pr[-1, 1] + 1:].sum()
    for (i1, j1), (i2, j2) in zip(pr[:-1], pr[1:]):
        ga, gb = list(range(i1 + 1, i2)), list(range(j1 + 1, j2)); k = min(len(ga), len(gb))
        for t in range(k):
            v = nA[ga[t]] + nB[gb[t]]
            if K(ga[t], gb[t]) <= 0: hin += v
            else: slp += v
        ins += nA[ga[k:]].sum() + nB[gb[k:]].sum()
    tot = ins + hin + slp; gap = nA.sum() + nB.sum() - nA[pr[:, 0]].sum() - nB[pr[:, 1]].sum()
    frac.append((ins / tot, hin / tot, slp / tot, abs(tot - gap) / gap))
    vA = 0.75 * sA ** 2 - nA; vneg.append(float((vA <= 0).mean()))
fr = np.array(frac)
print('S 재현 상대 차 최대', max(rel), '중앙', np.median(rel))
print('빈칸 질량 몫 중앙 [삽입형, 경첩형(K≤0), 밀림(K>0)]', np.median(fr[:, :3], 0).round(3), '사분위 경첩', np.percentile(fr[:, 1], [25, 75]).round(3), '합 대 빈칸 질량 상대 차 최대', fr[:, 3].max())
print('퍼짐 대용 v = 0.75 s̄² − ‖m‖² ≤ 0 인 변 비율 중앙', np.median(vneg))
```

Output:

```
S 재현 상대 차 최대 5.986152210616383e-05 중앙 4.3328119880274994e-16
빈칸 질량 몫 중앙 [삽입형, 경첩형(K≤0), 밀림(K>0)] [0.951 0.049 0.   ] 사분위 경첩 [0.023 0.072] 합 대 빈칸 질량 상대 차 최대 5.131670388268424e-15
퍼짐 대용 v = 0.75 s̄² − ‖m‖² ≤ 0 인 변 비율 중앙 0.0
```

Cell `HpplqbJChcyO`, run requested 2026-10-21 09:17:42, output returned 09:23:14 UTC.

```python
from scipy.stats import rankdata
tmd = {tuple(json.loads(l)['ab']): json.loads(l).get('tm_long') for l in open(f'{R}/derived/r41/tm_med.jsonl')}
def decomp(base):
    rows = []; cache = {}
    def ld(d):
        if d not in cache:
            if len(cache) > 300: cache.clear()
            z = np.load(f'{base}/dom/{d}.npz'); cache[d] = z['m'].astype(np.float64)
        return cache[d]
    for l in open(base + '/pibar.jsonl'):
        p = json.loads(l); a, b = p['ab']; mA, mB = ld(a), ld(b); pr = np.array(p['best'], int); pr = pr[np.argsort(pr[:, 0])]
        U = proc(np.einsum('kxz,kyz->xy', mB[pr[:, 1]], mA[pr[:, 0]])); nA, nB = (mA ** 2).sum((1, 2)), (mB ** 2).sum((1, 2))
        ov = nA[:pr[0, 0]].sum() + nA[pr[-1, 0] + 1:].sum() + nB[:pr[0, 1]].sum() + nB[pr[-1, 1] + 1:].sum(); ex = hin = slp = 0.0; nh = 0
        for (i1, j1), (i2, j2) in zip(pr[:-1], pr[1:]):
            ga, gb = np.arange(i1 + 1, i2), np.arange(j1 + 1, j2); k = min(len(ga), len(gb))
            if k:
                Kv = np.einsum('kxz,kxz->k', mA[ga[:k]], np.einsum('xy,kyz->kxz', U, mB[gb[:k]])); v = nA[ga[:k]] + nB[gb[:k]]
                hin += v[Kv <= 0].sum(); slp += v[Kv > 0].sum(); nh += int((Kv <= 0).sum())
            ex += nA[ga[k:]].sum() + nB[gb[k:]].sum()
        nb_ = (len(mA) + len(mB)) / 2; s = np.sqrt(nb_)
        rows.append(dict(ab=(a, b), ov=ov / s, ex=ex / s, hin=hin / s, slp=slp / s, gap=(ov + ex + hin + slp) / s, dfrac=abs(len(mA) - len(mB)) / max(len(mA), len(mB)), tm=tmd.get((a, b)), nh=nh))
    return rows
for nm, base in (('C', CC), ('H', f'{R}/derived/r38/p3')):
    rw = [r for r in decomp(base) if r['tm'] is not None]; A_ = lambda k: np.array([r[k] for r in rw], float)
    gap, ov, ex, hin, slp, df, tm = (A_(k) for k in ('gap', 'ov', 'ex', 'hin', 'slp', 'dfrac', 'tm'))
    tot = gap.sum(); print(nm, len(rw), '질량 몫 합계: 끝 겹침 %.3f 구간 초과 %.3f 경첩(K≤0) %.3f 밀림(K>0) %.3f' % (ov.sum() / tot, ex.sum() / tot, hin.sum() / tot, slp.sum() / tot),
          '경첩 몫 쌍별 중앙 %.3f' % np.median(hin / gap), '경첩 변 수 중앙', np.median(A_('nh')))
    Z = np.column_stack([np.ones(len(rw)), rankdata(tm), rankdata(df)]); res = lambda y: rankdata(y) - Z @ np.linalg.lstsq(Z, rankdata(y), rcond=None)[0]
    rg_ = res(gap); print('   TM·dfrac 잔차 뒤: 빈칸 잔차와 각 성분 잔차의 상관', {k: round(float(np.corrcoef(rg_, res(v))[0, 1]), 3) for k, v in (('끝 겹침', ov), ('구간 초과', ex), ('경첩', hin), ('삽입형=끝+초과', ov + ex))})
```

Output:

```
C 4152 질량 몫 합계: 끝 겹침 0.067 구간 초과 0.888 경첩(K≤0) 0.045 밀림(K>0) 0.000 경첩 몫 쌍별 중앙 0.042 경첩 변 수 중앙 3.0
   TM·dfrac 잔차 뒤: 빈칸 잔차와 각 성분 잔차의 상관 {'끝 겹침': 0.065, '구간 초과': 0.938, '경첩': 0.28, '삽입형=끝+초과': 0.992}
H 12003 질량 몫 합계: 끝 겹침 0.162 구간 초과 0.787 경첩(K≤0) 0.052 밀림(K>0) 0.000 경첩 몫 쌍별 중앙 0.040 경첩 변 수 중앙 1.0
   TM·dfrac 잔차 뒤: 빈칸 잔차와 각 성분 잔차의 상관 {'끝 겹침': 0.226, '구간 초과': 0.734, '경첩': 0.156, '삽입형=끝+초과': 0.975}
```

**Re-run on 31 October 2026.** The cell code above, run again on the same input files, printed the same line for the
confirmatory set (C 4152: 0.067, 0.888, 0.045, 0.000; residual correlations 0.065, 0.938, 0.28, 0.992). The first-study line
(H), which the manuscript does not use, was not re-run.


## 5. Where correspondence shifts occur (21 October 2026)

Manuscript, bond-level section: the shifts of the correspondence by a few bonds occurred about equally often within and
between ensembles, twice as often per bond in the 10 bonds at each chain end as in the core, but with most events in
the core. Inputs: `derived/r41/dcomp_gpu/<domain>.json` (per-cell fractions of bonds that are unshared, `gap_frac`, or
shifted, `shift_frac`, in the 10 bonds at each end, `term`, and in the core, between ensembles, `cross`, and within one
ensemble, `within`), written by `code/r41_dcomp.py` (part gpu, round41 g2) and summarised as `gpu` in
`results/r41_dcomp.json`, and the bond counts `ne` in `derived/r38/s0/e7_30/` and `e7_100/`. The cross-minus-within
differences are 0.25–0.31 percentage points; per bond, shifts are 0.124 at the ends against 0.058 in the core; of all
shift events between ensembles, a median 0.387 fall in the 20 end bonds (so about 61% in the core).

Cell `HpplqbJChcyO`, run requested 2026-10-21 14:45:32, output returned 14:45:33 UTC.

```python
g = J('r41_dcomp'); print(g['checks']); print('gpu', json.dumps(g.get('gpu'), ensure_ascii=False)); print(g['reading'])
import glob
G = [json.load(open(f)) for f in glob.glob(f'{R}/derived/r41/dcomp_gpu/*.json')]
rows = [(T, c) for j in G for k, c in j['cells'].items() for T in [k] if k.startswith('MD')]
for part in ('cross', 'within'):
    for x in ('gap_frac_term', 'gap_frac_core', 'shift_frac_term', 'shift_frac_core', 'frac_c2_gt_c1', 'frac_c2_lt_c1'):
        v = np.array([c[part][x] for _, c in rows]); print(part, x, round(float(np.median(v)), 4), [round(float(t), 4) for t in np.percentile(v, [25, 75])])
for T in ('348', '379', '413', '450'):
    sub = [c for k, c in rows if k[2:5] == T]
    print(T, len(sub), {x: round(float(np.median([c['cross'][x] for c in sub])), 4) for x in ('gap_frac_term', 'gap_frac_core', 'shift_frac_term', 'shift_frac_core', 'frac_c2_gt_c1')})
# 교차 − 같은 쪽 (끝/중심의 초과)
for x in ('gap_frac_term', 'gap_frac_core', 'shift_frac_term', 'shift_frac_core'):
    v = np.array([c['cross'][x] - c['within'][x] for _, c in rows]); print('교차−같은쪽', x, round(float(np.median(v)), 4))
```

Output:

```
{'a_vs_A2000_maxrel': 0.0, 'Sb_iter_maxrel': 2.5830321892508642e-15, 'viol_ab': 0, 'viol_bc1': 0, 'n_pairs': 28962000, 'gpu_full_vs_stored_maxrel': 0.0, 'n_gpu': 97}
gpu {"gap_frac_term": {"median": 0.1074965, "q25": 0.057036875, "q75": 0.17311}, "gap_frac_core": {"median": 0.046312112068965516, "q25": 0.027818351841028637, "q75": 0.06292354283707866}, "shift_frac_term": {"median": 0.12364160526088659, "q25": 0.05110920856060537, "q75": 0.2395405262774736}, "shift_frac_core": {"median": 0.05821813723861202, "q25": 0.02146178130849964, "q75": 0.10411769735665727}, "frac_c2_gt_c1": {"median": 0.00025, "q25": 0.0, "q75": 0.00102125}, "n_cells": 464}
{'D2': {'shares': {'a→b': -0.015274405526762671, 'b→c1': 0.2940324418424907, 'c1→c2': 0.7212419636842721}, 'reading': '등록 미끄러짐(끝 재정렬·주기 구조)'}, 'ratio': {'shares': {'a→b': -0.03371513197490173, 'b→c1': 0.24450536170182544, 'c1→c2': 0.7892097702730763}, 'reading': '등록 미끄러짐(끝 재정렬·주기 구조)'}}
cross gap_frac_term 0.1075 [0.057, 0.1731]
cross gap_frac_core 0.0463 [0.0278, 0.0629]
cross shift_frac_term 0.1236 [0.0511, 0.2395]
cross shift_frac_core 0.0582 [0.0215, 0.1041]
cross frac_c2_gt_c1 0.0003 [0.0, 0.001]
cross frac_c2_lt_c1 0.9822 [0.9104, 0.999]
within gap_frac_term 0.1008 [0.0519, 0.1562]
within gap_frac_core 0.0406 [0.0253, 0.054]
within shift_frac_term 0.1059 [0.0452, 0.2139]
within shift_frac_core 0.0475 [0.0195, 0.0829]
within frac_c2_gt_c1 0.0004 [0.0, 0.0011]
within frac_c2_lt_c1 0.9477 [0.832, 0.998]
348 176 {'gap_frac_term': 0.1037, 'gap_frac_core': 0.0441, 'shift_frac_term': 0.1017, 'shift_frac_core': 0.0545, 'frac_c2_gt_c1': 0.0003}
379 138 {'gap_frac_term': 0.1084, 'gap_frac_core': 0.0453, 'shift_frac_term': 0.1353, 'shift_frac_core': 0.0605, 'frac_c2_gt_c1': 0.0003}
413 99 {'gap_frac_term': 0.1082, 'gap_frac_core': 0.0524, 'shift_frac_term': 0.115, 'shift_frac_core': 0.0678, 'frac_c2_gt_c1': 0.0003}
450 51 {'gap_frac_term': 0.1226, 'gap_frac_core': 0.0464, 'shift_frac_term': 0.1508, 'shift_frac_core': 0.0511, 'frac_c2_gt_c1': 0.0003}
교차−같은쪽 gap_frac_term 0.0028
교차−같은쪽 gap_frac_core 0.0025
교차−같은쪽 shift_frac_term 0.0031
교차−같은쪽 shift_frac_core 0.0031
```

Cell `HpplqbJChcyO`, run requested 2026-10-21 14:45:51, output returned 14:45:57 UTC.

```python
NE = {}
for f in glob.glob(f'{R}/derived/r38/s0/e7_30/*.json') + glob.glob(f'{R}/derived/r38/s0/e7_100/*.json'):
    j = json.load(open(f))
    if 'ne' in j: NE[j['dom']] = j['ne']
sh_t, gp_t, nes = [], [], []
for j in G:
    ne = NE.get(j['dom'])
    if not ne: continue
    for k, c in j['cells'].items():
        if not k.startswith('MD'): continue
        x = c['cross']; nt = 20; nc = ne - 20
        st_, sc_ = x['shift_frac_term'] * nt, x['shift_frac_core'] * nc; gt_, gc_ = x['gap_frac_term'] * nt, x['gap_frac_core'] * nc
        if st_ + sc_ > 0: sh_t.append(st_ / (st_ + sc_))
        if gt_ + gc_ > 0: gp_t.append(gt_ / (gt_ + gc_))
        nes.append(ne)
print('도메인 ne 중앙', np.median(nes), '끝 변 몫(20/ne) 중앙', np.median(20 / np.array(nes)))
print('밀림 사건 중 끝 몫 중앙', round(float(np.median(sh_t)), 3), np.percentile(sh_t, [25, 75]).round(3), '빈칸 변 중 끝 몫 중앙', round(float(np.median(gp_t)), 3), np.percentile(gp_t, [25, 75]).round(3))
```

Output:

```
도메인 ne 중앙 89.0 끝 변 몫(20/ne) 중앙 0.2247191011235955
밀림 사건 중 끝 몫 중앙 0.387 [0.199 0.688] 빈칸 변 중 끝 몫 중앙 0.414 [0.255 0.595]
```

**Re-run on 31 October 2026.** The same code on the same input files printed the same values (cross − within
0.0028, 0.0025, 0.0031, 0.0031; 464 cells; median share of shift events in the end bonds 0.387 [0.199, 0.688]).


## 6. Chain breaks (28 October 2026)

Methods: no domain used in the analyses had a Cα–Cα bond longer than 4.2 Å (s > 1.1). For the homolog studies the check
used the frame-averaged length ratio `sbar` stored for each domain (`derived/r38/p3/dom/`, first homolog study, 1,496
domains with `sbar`; `derived/r38/c/dom/`, confirmatory set, 949), and for the 336 calibration domains every bond of the
first 320 K frame of replicate 0 (`derived/mdcath_5T/<domain>.npz`, `320_0_X`). The output tuple is (domains checked,
domains with a long bond, examples).

Cell `gegSq8EQ72mh`, run requested 2026-10-28 05:15:39, output returned 05:16:03 UTC.

```python
# 12. 원고 점검(읽기만): (a) H-풀 12,003 → 11,991 쌍이 된 까닭, (b) 사슬 끊김(긴 결합)이 있는 도메인 수
import numpy as np, glob
p3 = json.load(open(f'{R}/results/r38_p3.json'))
print('r38_p3 keys:', list(p3)[:20])
for k in ('dropouts', 'n_pairs', 'n_rows', 'counts'):
    if k in p3: print(k, p3[k])
def breaks(dirp, thr=1.1):   # 평균 길이 비 s̄ > 1.1 (Cα–Cα > 4.2 Å) 인 결합이 있는 도메인
    n = nb = 0; ex = []
    for f in glob.glob(dirp + '/*.npz'):
        z = np.load(f)
        if 'sbar' not in z.files: continue
        n += 1; s = z['sbar']
        if (s > thr).any(): nb += 1; ex.append((os.path.basename(f)[:-4], int((s > thr).sum()), round(float(s.max()), 2)))
    return n, nb, ex[:8]
for nm, d in (('H-풀', R + '/derived/r38/p3/dom'), ('C', R + '/derived/r38/c/dom')): print(nm, breaks(d))
z = np.load(sorted(glob.glob(R + '/derived/mdcath_5T/*.npz'))[0]); print('줄기 npz keys', z.files[:10])
```

Output:

```
r38_p3 keys: ['written', 'verdict', 'main', 'joint_gate', 'freedman_lane', 'sensitivity', 'dropouts', 'report', 'pending_sensitivity', 'rule']
dropouts {'n_pairs': 14189, 'with_pibar': 12003, 'main_ok': 11991, 'inter_ok': 11372, 'main_fail_minshared': 12, 'inter_fail_minshared': 619}
H-풀 (1496, 0, [])
C (949, 0, [])
줄기 npz keys ['seq', 'native_i', 'native_j', 'native_d', '320_0_Q', '320_0_nframes', '320_0_X', '320_0_rg', '320_0_dssp', '320_1_Q']
```

Cell `46NZaGFB7_0m`, run requested 2026-10-28 05:16:17, output returned 05:22:39 UTC.

```python
# 13. 줄기 336 도메인의 긴 결합(첫 프레임 Cα–Cα > 4.2 Å) 수 (읽기만)
nb = 0; n = 0; mx = []
for f in sorted(glob.glob(R + '/derived/mdcath_5T/*.npz')):
    z = np.load(f); X = z['320_0_X'][0]; n += 1
    ca = X[:, 0] if X.ndim == 3 else X.reshape(-1, 3, 3)[:, 0]
    l = np.linalg.norm(np.diff(ca, axis=0), axis=1); mx.append(float(l.max()))
    if (l > 4.2).any(): nb += 1
print('줄기', n, '도메인, 긴 결합 있는 도메인', nb, '| 최대 Cα–Cα 중앙/최대', round(float(np.median(mx)), 3), round(max(mx), 3), '| X shape', X.shape)
```

Output:

```
줄기 336 도메인, 긴 결합 있는 도메인 0 | 최대 Cα–Cα 중앙/최대 3.984 4.13 | X shape (100, 3, 3)
```

**Re-run on 31 October 2026.** The same computations on the same files printed the same results: H-풀 (the H-pool: 1496, 0, []),
C (949, 0, []) and, for the calibration domains, 336 domains with no long bond (largest Cα–Cα bond: median 3.984 Å,
maximum 4.13 Å). For the calibration domains the re-run read only the array `320_0_X` of each file, several files at a time,
to avoid slow reads of the whole files; the computation on that array was the same.

## 7. framedchain 0.4.1 against the as-run code (28 October 2026)

Methods: on synthetic inputs with fixed seeds, framedchain 0.4.1 reproduced the as-run code exactly. The script is
`code/r44_repro.py` (md5 `a8bd5a23aa5f896d1e267dac160fedc0`); it was run locally (CPU, float64). The first run, at
05:13:56 UTC, imported the as-run modules from the working copy; the paths were then changed to the copies in this
repository (the final version, with the md5 above) and it was run again at 05:15:00 UTC. Both printed the
same table. Command (the local working-copy path is written <atlas>) and printed output of the first run:

```bash
cd <atlas> && PYTHONDONTWRITEBYTECODE=1 timeout 600 uv run --offline -q --no-project --python 3.12 --with numpy --with torch --with scipy --with numba python code/r44_repro.py "$(pwd -W)" 2>&1 | tail -25
```

```
framespace 0.2.0 framespace_pkg | framedchain 0.4.1
distance               max |diff| = 0
align score            max |diff| = 0
align match            max |diff| = 0
bnb S_lb/S_ub          max |diff| = 0
mean_point             max |diff| = 0
ensemble_w2            max |diff| = 0
MeanStart center       max |diff| = 0
MeanStart starts       max |diff| = 0
red_pairs d            max |diff| = 0
red_pairs gap          max |diff| = 0
estats D2              max |diff| = 0
estats ratio           max |diff| = 0
estats D2/sqrt(n)      max |diff| = 0
ALL_EQUAL
```

The second run printed the same lines (the end of its output is shown, followed by the md5 of the final script):

```
50:src = open(os.path.join(A, 'framedchain-paper', 'code', 'r38_s0.py'), encoding='utf-8').read()
estats ratio           max |diff| = 0
estats D2/sqrt(n)      max |diff| = 0
ALL_EQUAL
a8bd5a23aa5f896d1e267dac160fedc0 *code/r44_repro.py
```

**Re-run on 31 October 2026.** The deposited script (same md5) printed the same table, ending in `ALL_EQUAL`.

