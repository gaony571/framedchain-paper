# Interactive checks: code and printed output

Seven results were computed in interactive Colab cells, or by a script whose printed output was not saved to a file:
the comparison of 12 test pairs with TM-align (S2 Text 8.1; section 1), the reversed chains (S2 Text 8.5; section 2),
the equipment check before the computational-cost measurement (S2 Text 7; section 3), the decomposition of the
unshared part between mean points (section 4), where correspondence shifts occur (section 5), the chain-break check
(section 6) and the reproduction of the as-run code by framedchain 0.4.1 (section 7). Three further cells wrote result
files that the paper uses: the violation rates of the factorial comparison behind Fig 2A (S2 Text 8.2; section 8),
the size and survival audit behind the survival row of Table 2 (section 9) and the diagnostic of the window difference
behind the trajectory-window row of Table 2 (section 10). Two more groups of cells wrote and ran short scripts whose
result files the paper cites: the coverage simulations of the first thermophile test (Methods; section 11) and the
diagnostics of the failed truncated-chain check (S2 Text 5; section 12). Two further cells, run after a result and
not specified in advance, are cited in S1 Table: the check of the 450 K version of the temperature conversion without a
missing value, which wrote a result file (section 13), and the description of the accuracy of the 128-start values after
the specification error of the certification (S2 Text 1; section 14), later recomputed by `code/r39_sec6.py`. Section 15
gives the commands and printed output of the negative-type computations of S3 Text 2, run locally with Node.js and
Python. Section 16 gives two cells of 3 November 2026: one, run while the fragments of Fig 1 were chosen, that listed the
eligible stretches, and one, run after the figure was drawn, that bounded its alignment more tightly than the script does. The code and printed output
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
bonds at the same position turned beyond 120° (a descriptive decomposition, not specified in advance). Inputs:
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
ensemble, `within`), written by `code/r41_dcomp.py` (part gpu, round41 g2) and summarized as `gpu` in
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

Methods: no domain or chain checked had a Cα–Cα bond longer than 4.2 Å (s > 1.097; in the two homolog studies, a frame-averaged s̄ above 1.1, 4.21 Å); the sets not checked are listed there and in S1 Table. For the homolog studies the check
used the frame-averaged length ratio `sbar` stored for each domain, with the cut-off s̄ > 1.1 (4.21 Å for l₀ = 3.83 Å; the
code comment rounds it to 4.2 Å) (`derived/r38/p3/dom/`, first homolog study, 1,496
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

**The other sets (31 October 2026).** The sets not covered here (model-error domains, comparison across sequences, natural
insertion–deletion pairs, certified pairs, both thermophile studies, force-field chains and the static chains of the first
two Results sections; the sets still not checked are listed in S1 Table) were checked by `code/r49_chainbreak_fix1.py`, specified in a record before it was run
(`registration/prereg_2026-10-31_addendum12.json` and its correction); see `docs/RESULTS_MAP.md`, Round49.

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

## 8. Violation rates of the factorial comparison (16 October 2026)

S2 Text 8.2 and Fig 2A: the rates of triangle-inequality violations among the 161,700 triples of the 100 structures of
block B, for every normalization, transform and choice of d₀ at every optimization level, and for 1 − TM-score with
TM-align's own alignments. The first cell reads the pair table written by `code/pu_tm_factorial_v2.py`
(`results/pu_tm_factorial_B_pairs.tsv`, 4,950 pairs) and its setup file, applies `triangle_rate` of
`framespace_wip/stage0/rulers.py` (a triple violates when its longest side exceeds the sum of the other two by more than
the tolerance times the longest side, plus 10⁻⁹) at tolerances 10⁻⁹, 10⁻⁶ and 10⁻⁴, to all triples and to those among the
truncated chains (names containing `|`), and writes `results/pu_tm_factorial_B_triangle.tsv`, which
`code/make_figs.py` reads for Fig 2A (rows with tolerance 10⁻⁴). `triangle_rate` was last edited, and reloaded in the
kernel, at 16:19 UTC (section 1); the module's later edit (17:28 UTC) changed `wl_align` only, so the function used is
the one in the deposited `rulers.py`. The comment in the first line reads "block B analysis (main verdict of H-U2); run when
finished", H-U2 being the factorial comparison. The second cell prints the rows of the levels used in the paper
(strongest level for fixed d₀, the next for d₀ of the normalizing length) and the starting-point curve at tolerance
10⁻⁹.

Cell `aC6JZGe-6qBE`, run requested 20:16:56, output returned 20:20:21 UTC.

```python
# 묶음 B 분석 (H-U2 주 판정). 끝나면 돌린다.
pB = pd.read_csv(ROOT + '/results/pu_tm_factorial_B_pairs.tsv', sep='\t'); stB = json.load(open(ROOT + '/results/pu_tm_factorial_B_setup.json'))
blB = stB['block']; biB = {n_: i_ for i_, n_ in enumerate(blB)}; NB_ = len(blB); print(len(pB), 'pairs of', NB_ * (NB_ - 1) // 2)
iaB, ibB = pB.a.map(biB).values, pB.b.map(biB).values
def matB(v):
    M = np.full((NB_, NB_), np.nan); M[iaB, ibB] = v; M[ibB, iaB] = v; np.fill_diagonal(M, 0); return M
nB, mB = pB.n.values.astype(float), pB.m.values.astype(float); Lb_, Ls_, Ll_ = 0.5 * (nB + mB), np.minimum(nB, mB), np.maximum(nB, mB)
cellsB = {}
for lv in ('k1', 'k4', 'k16', 'hi', 'max'):
    S = pB[f'S_fix_{lv}'].values; D2 = np.clip(Lb_ - S, 0, None)
    cellsB[('fix', 'none', 'sqrt', lv)] = np.sqrt(D2); cellsB[('fix', 'none', 'lin', lv)] = D2
    for nm, L in (('avg', Lb_), ('short', Ls_), ('long', Ll_)):
        cellsB[('fix', nm, 'lin', lv)] = 1 - S / L; cellsB[('fix', nm, 'sqrt', lv)] = np.sqrt(np.clip(1 - S / L, 0, None))
for nm, L in (('short', Ls_), ('long', Ll_), ('avg', Lb_)):
    for lv in ('k4', 'hi'):
        S = pB[f'S_{nm}_{lv}'].values; D2 = np.clip(Lb_ - S, 0, None)
        cellsB[(f'var_{nm}', 'none', 'sqrt', lv)] = np.sqrt(D2); cellsB[(f'var_{nm}', 'none', 'lin', lv)] = D2
        cellsB[(f'var_{nm}', nm, 'lin', lv)] = 1 - S / L; cellsB[(f'var_{nm}', nm, 'sqrt', lv)] = np.sqrt(np.clip(1 - S / L, 0, None))
shB = nB <= mB; tsB = np.where(shB, pB.tm_a, pB.tm_b); tlB = np.where(shB, pB.tm_b, pB.tm_a)
for nm, v in (('short', 1 - tsB), ('long', 1 - tlB), ('avg', 1 - 0.5 * (tsB + tlB))): cellsB[('TMalign', nm, 'lin', 'tmalign')] = v
trunc_ = [n_ for n_ in blB if '|' in n_]; tidx_ = np.array([biB[n_] for n_ in trunc_])
rowsB = []
for key, v in cellsB.items():
    M = matB(v)
    if np.isnan(M).any(): continue
    for tol in (1e-9, 1e-6, 1e-4):
        a_ = s0.triangle_rate(M, rel_tol=tol); b_ = s0.triangle_rate(M[np.ix_(tidx_, tidx_)], rel_tol=tol)
        rowsB.append(dict(d0=key[0], norm=key[1], form=key[2], level=key[3], tol=tol, rate_all=a_[0], n_all=a_[1], max_rel=a_[3], n_trunc=b_[1], rate_trunc=b_[0]))
resB = pd.DataFrame(rowsB); resB.to_csv(ROOT + '/results/pu_tm_factorial_B_triangle.tsv', sep='\t', index=False)
pd.set_option('display.width', 250)
print(resB[resB.tol == 1e-4].pivot_table(index=['d0', 'norm', 'form'], columns='level', values='n_all').to_string())
print(resB[(resB.tol == 1e-4) & (resB.d0 == 'fix') & (resB.norm == 'none') & (resB.form == 'sqrt')][['level', 'n_all', 'rate_all', 'max_rel', 'n_trunc']].to_string())
for c_ in ('fix', 'short', 'long', 'avg'):
    if c_ == 'fix': dS = pB.S_fix_max - pB.S_fix_hi; print('fix: max>hi pairs', round(float((dS > 1e-6).mean()), 4), 'max<hi', int((dS < -1e-6).sum()))
    else: dS = pB[f'S_{c_}_hi'] - pB[f'S_{c_}_k4']; print(c_, 'hi>k4', round(float((dS > 1e-6).mean()), 3))
```

Output:

```
4950 pairs of 4950
level                     hi      k1     k16      k4     max  tmalign
d0        norm  form                                                 
TMalign   avg   lin      NaN     NaN     NaN     NaN     NaN     16.0
          long  lin      NaN     NaN     NaN     NaN     NaN      0.0
          short lin      NaN     NaN     NaN     NaN     NaN   3849.0
fix       avg   lin      3.0     6.0     4.0     6.0     3.0      NaN
                sqrt     0.0     0.0     0.0     0.0     0.0      NaN
          long  lin      0.0     0.0     0.0     0.0     0.0      NaN
                sqrt     0.0     0.0     0.0     0.0     0.0      NaN
          none  lin     76.0   334.0   166.0   238.0    42.0      NaN
                sqrt     0.0     0.0     0.0     0.0     0.0      NaN
          short lin   4472.0  4298.0  4432.0  4410.0  4532.0      NaN
                sqrt  3857.0  3884.0  3878.0  3878.0  3863.0      NaN
var_avg   avg   lin     10.0     NaN     NaN    10.0     NaN      NaN
                sqrt     0.0     NaN     NaN     0.0     NaN      NaN
          none  lin      6.0     NaN     NaN    30.0     NaN      NaN
                sqrt     0.0     NaN     NaN     0.0     NaN      NaN
var_long  long  lin      0.0     NaN     NaN     0.0     NaN      NaN
                sqrt     0.0     NaN     NaN     0.0     NaN      NaN
          none  lin      9.0     NaN     NaN    41.0     NaN      NaN
                sqrt     0.0     NaN     NaN     0.0     NaN      NaN
var_short none  lin     78.0     NaN     NaN   188.0     NaN      NaN
                sqrt     0.0     NaN     NaN     0.0     NaN      NaN
          short lin   3854.0     NaN     NaN  3887.0     NaN      NaN
                sqrt  3848.0     NaN     NaN  3885.0     NaN      NaN
   level  n_all  rate_all   max_rel  n_trunc
2     k1      0       0.0 -0.312029        0
26    k4      0       0.0 -0.316110        0
50   k16      0       0.0 -0.326792        0
74    hi      0       0.0 -0.337067        0
98   max      0       0.0 -0.344473        0
fix: max>hi pairs 0.3772 max<hi 525
short hi>k4 0.863
long hi>k4 0.872
avg hi>k4 0.87
```

Cell `ErAfVzOSRR0n`, run requested 20:20:46, output returned 20:20:46 UTC.

```python
r4_ = resB[resB.tol == 1e-4]
show_ = r4_[((r4_.d0 == 'fix') & (r4_.level == 'max')) | ((r4_.d0.str.startswith('var')) & (r4_.level == 'hi')) | (r4_.d0 == 'TMalign')]
print(show_[['d0', 'norm', 'form', 'level', 'n_all', 'rate_all', 'n_trunc', 'rate_trunc', 'max_rel']].to_string())
print('trunc triples', len(tidx_) * (len(tidx_) - 1) * (len(tidx_) - 2) // 6)
r9_ = resB[(resB.tol == 1e-9) & (resB.d0 == 'fix') & (resB.norm == 'none')][['form', 'level', 'n_all']]; print(r9_.to_string())
```

Output:

```
            d0   norm  form    level  n_all  rate_all  n_trunc  rate_trunc   max_rel
98         fix   none  sqrt      max      0  0.000000        0    0.000000 -0.344473
101        fix   none   lin      max     42  0.000260       23    0.000672  0.050752
104        fix    avg   lin      max      3  0.000019        3    0.000088  0.092256
107        fix    avg  sqrt      max      0  0.000000        0    0.000000 -0.344330
110        fix  short   lin      max   4532  0.028027     2280    0.066628  1.000000
113        fix  short  sqrt      max   3863  0.023890     2273    0.066423  1.000000
116        fix   long   lin      max      0  0.000000        0    0.000000 -0.297301
119        fix   long  sqrt      max      0  0.000000        0    0.000000 -0.606343
134  var_short   none  sqrt       hi      0  0.000000        0    0.000000 -0.362079
137  var_short   none   lin       hi     78  0.000482       57    0.001666  0.035075
140  var_short  short   lin       hi   3854  0.023834     2271    0.066365  1.000000
143  var_short  short  sqrt       hi   3848  0.023797     2267    0.066248  1.000000
158   var_long   none  sqrt       hi      0  0.000000        0    0.000000 -0.351051
161   var_long   none   lin       hi      9  0.000056        4    0.000117  0.034999
164   var_long   long   lin       hi      0  0.000000        0    0.000000 -0.323369
167   var_long   long  sqrt       hi      0  0.000000        0    0.000000 -0.622886
182    var_avg   none  sqrt       hi      0  0.000000        0    0.000000 -0.356412
185    var_avg   none   lin       hi      6  0.000037        4    0.000117  0.011273
188    var_avg    avg   lin       hi     10  0.000062       10    0.000292  0.091933
191    var_avg    avg  sqrt       hi      0  0.000000        0    0.000000 -0.344569
194    TMalign  short   lin  tmalign   3849  0.023803     2265    0.066189  1.000000
197    TMalign   long   lin  tmalign      0  0.000000        0    0.000000 -0.280855
200    TMalign    avg   lin  tmalign     16  0.000099       16    0.000468  0.354345
trunc triples 34220
    form level  n_all
0   sqrt    k1      0
3    lin    k1    336
24  sqrt    k4      0
27   lin    k4    242
48  sqrt   k16      0
51   lin   k16    168
72  sqrt    hi      0
75   lin    hi     77
96  sqrt   max      0
99   lin   max     43
```

## 9. Size and survival audit of the unit (18 October 2026)

Table 2, survival row (round37 R1): Δ²/√n of heating from 320 K to 348, 379, 413 and 450 K in aligned coordinates, from
the temperature cells of `code/ladder_temp_single_v2.py` (`derived/ladder_temp_single_v2/<domain>.json`), as superfamily
medians (bootstrap over superfamilies, 2,000 resamples) over all domains (`full`) and over the 149 domains with at least
100 folded frames, pooled over the replicates, at each of the five temperatures (`restricted`; `n_folded` in
`derived/c1_core/<domain>.json`, the folded frames of all replicates of a temperature); in both versions each cell
uses the replicates with at least 40 folded frames (`code/ladder_temp_single_v2.py`), which is the "at least 40 folded
frames per replicate" of Table 2; the
same cell gives the medians by quartile of the number of bonds (`R1_size`). `cath` maps each CATH domain to its
superfamily label and length (from `raw/cath/cath-domain-list.txt`; defined in an earlier cell of the same runtime).
The survival ratios of Table 2 are `ratio` at 379, 413 and 348 K (1.045, 0.946 and 1.586), and the `full` medians are
the constants `full_for_ratio` of `code/r37_trunk.py`. A first run of the cell stopped with a pandas error (a column
named `ne` clashes with `DataFrame.ne`); the cell was edited (column `nedge`) and run again, and the version shown is
the one that produced the output. The next cell wrote the printed dictionary, with a time stamp, to
`results/r37_size_survival.json` (its first line, shown); the rest of that cell wrote `code/r37_trunk.py`, identical to
the deposited file (md5 `61148d6a8c116fe06bf4f3c3ca6f9099`), and started it, and is left out.

Cell `X8J43YO6_TGu`, run requested 2026-10-18 17:09:32, output returned 17:09:37 UTC.

```python
Js = [json.load(open(q)) for q in glob.glob(ROOT + '/derived/ladder_temp_single_v2/*.json')]
nf = {}
for q in glob.glob(ROOT + '/derived/c1_core/*.json'):
    j = json.load(open(q)); nf[j['dom']] = j.get('n_folded', {})
cathH = {d: v[0] for d, v in cath.items()}
rngb = np.random.default_rng(37)
def sfmed(vals, doms, n=2000):
    s = pd.Series(vals).groupby([cathH.get(x, x) for x in doms]).median(); b = [np.median(s.values[rngb.integers(0, len(s), len(s))]) for _ in range(n)]
    return round(float(s.median()), 4), [round(float(x), 4) for x in np.percentile(b, [2.5, 97.5])], int(len(s))
rows = []
for j in Js:
    for T in ('348', '379', '413', '450'):
        if T in j['cross1'] and T in j['floorT']:
            rows.append(dict(dom=j['dom'], T=T, nedge=j['ne'], d2n=(j['cross1'][T] - (j['floor1'] + j['floorT'][T]) / 2) / np.sqrt(j['ne'])))
R1 = pd.DataFrame(rows); OUT37 = {}
sz = {}
for T, g in R1.groupby('T'):
    g = g.assign(q=pd.qcut(g['nedge'], 4, labels=False, duplicates='drop')); med = g.groupby('q')['d2n'].median()
    sz[T] = dict(n=len(g), ne_edges=[int(x) for x in np.quantile(g['nedge'], [0, .25, .5, .75, 1])], quartile_median=[round(float(x), 4) for x in med.values],
                 max_over_min=round(float(med.max() / med.min()), 3) if med.min() > 0 else None, spearman_ne=round(float(g['nedge'].corr(g['d2n'], method='spearman')), 3))
OUT37['R1_size'] = sz
surv = {d for d, v in nf.items() if v and all(v.get(T, 0) >= 100 for T in ('320', '348', '379', '413', '450'))}
sv = {}
for T, g in R1.groupby('T'):
    full = sfmed(g['d2n'].values, g['dom'].values); gs = g[g['dom'].isin(surv)]; res_ = sfmed(gs['d2n'].values, gs['dom'].values) if len(gs) > 5 else None
    sv[T] = dict(full=full, restricted=res_, ratio=round(res_[0] / full[0], 3) if res_ and full[0] else None, n_full=len(g), n_restricted=len(gs))
OUT37['R1_survival'] = dict(n_surviving_domains=len(surv), **sv)
print(json.dumps(OUT37, ensure_ascii=False, indent=1))
```

Output:

```
{
 "R1_size": {
  "348": {
   "n": 314,
   "ne_edges": [
    49,
    70,
    86,
    109,
    149
   ],
   "quartile_median": [
    0.0026,
    0.0041,
    0.004,
    0.002
   ],
   "max_over_min": 2.096,
   "spearman_ne": -0.047
  },
  "379": {
   "n": 278,
   "ne_edges": [
    49,
    69,
    86,
    109,
    149
   ],
   "quartile_median": [
    0.0152,
    0.0192,
    0.0201,
    0.0229
   ],
   "max_over_min": 1.513,
   "spearman_ne": 0.071
  },
  "413": {
   "n": 218,
   "ne_edges": [
    49,
    70,
    86,
    108,
    149
   ],
   "quartile_median": [
    0.0278,
    0.0346,
    0.0359,
    0.0519
   ],
   "max_over_min": 1.872,
   "spearman_ne": 0.195
  },
  "450": {
   "n": 115,
   "ne_edges": [
    49,
    68,
    84,
    108,
    149
   ],
   "quartile_median": [
    0.0494,
    0.0375,
    0.0583,
    0.0831
   ],
   "max_over_min": 2.219,
   "spearman_ne": 0.351
  }
 },
 "R1_survival": {
  "n_surviving_domains": 149,
  "348": {
   "full": [
    0.0029,
    [
     0.0013,
     0.0051
    ],
    68
   ],
   "restricted": [
    0.0046,
    [
     0.0016,
     0.0075
    ],
    47
   ],
   "ratio": 1.586,
   "n_full": 314,
   "n_restricted": 149
  },
  "379": {
   "full": [
    0.0199,
    [
     0.0151,
     0.0259
    ],
    64
   ],
   "restricted": [
    0.0208,
    [
     0.0132,
     0.0256
    ],
    47
   ],
   "ratio": 1.045,
   "n_full": 278,
   "n_restricted": 149
  },
  "413": {
   "full": [
    0.0404,
    [
     0.0299,
     0.0462
    ],
    53
   ],
   "restricted": [
    0.0382,
    [
     0.0275,
     0.0517
    ],
    47
   ],
   "ratio": 0.946,
   "n_full": 218,
   "n_restricted": 148
  },
  "450": {
   "full": [
    0.0657,
    [
     0.0483,
     0.0723
    ],
    41
   ],
   "restricted": [
    0.0657,
    [
     0.0483,
     0.0723
    ],
    41
   ],
   "ratio": 1.0,
   "n_full": 115,
   "n_restricted": 115
  }
 }
}
```

Cell `8qMMN_3A_f7v`, run requested 17:10:11, output returned 17:10:12 UTC (first line; the cell printed `started`).

```python
OUT37['written'] = time.strftime('%F %T'); json.dump(OUT37, open(ROOT + '/results/r37_size_survival.json', 'w'), indent=1, ensure_ascii=False)
```

## 10. Diagnostic of the window difference (21 October 2026)

Results (calibration) and Table 2, trajectory-window row: in the edit distance the ratio of second-half to
whole-trajectory responses was 0.85–1.13 (part `e2` of `code/r39_core.py`), against 1.0–2.1 in the earlier audit with
aligned coordinates, and a diagnostic added after this result traced the difference to the folded-frame criterion and
the domain set rather than the distance. The cell below first added the entry `round39_posthoc_window_diag_2026-10-21`
to the first record (`results/prereg_2026-10-17.json`, with a copy in `results/_backup/2026-10-21/`; its stamp,
12:43:00, is the first output line). The entry gives, in Korean, the question, the design and the reading rule: run
part `e2` of `r39_core.py` unchanged on the aligned-coordinate values stored with the same replicates and frames in
the per-domain files of `r39_D_stem.py` (`derived/r39/D/<domain>.json`; keys `Aall` and `dmu2_A` instead of `D` and
`dmu2_D`; same seed 20261072, same domain set, at least 100 folded frames per replicate), and read the second-half /
whole ratio at +59, +93 and +130 K: "criterion and set" (기준·집합) if its 95% interval contains 1 in at least two of
the three steps, "distance" (거리) if the lower bound is at least 1.2 in at least two, and "not separated"
(가르지 못함) otherwise; +28 K is descriptive. The cell then wrote the modified copy of `r39_core.py` to
`/content/r39_e2_aligned_diag.py` (not kept: the deposited `r39_core.py` up to its `__main__` block with the
replacements shown, followed by a call of `part_e2()`) and ran it, which wrote `results/r39_E2_aligned_diag.json`.
The second output line gives the number of occurrences replaced (3 and 2; the deposited `r39_core.py` has the same
counts in `part_e2`, and the other replaced strings occur in it), and the `R39_E2` line the whole / second-half ratios
of the medians. The second run of the cell added the rule, the reading and the second-half / whole ratios (the
reciprocals of the stored whole / second-half ratios and of their interval bounds) to that file. In this cell `R` is
the project root and `C = R + '/code'`; `b`, used in the second run, is the backup folder set in the first. Times
printed inside the output come from the Colab clock. The character � in the first run's code is in the log as
written; it lies inside a label string of the output file.

Cell `HpplqbJChcyO`, run requested 2026-10-21 12:42:57, output returned 12:43:10 UTC.

```python
pf = f'{R}/results/prereg_2026-10-17.json'; P = json.load(open(pf)); K = 'round39_posthoc_window_diag_2026-10-21'
assert K not in P
P[K] = dict(stamped_utc=time.strftime('%F %T'), nature='사후 진단(판정을 바꾸지 않고 설명을 붙임; 사용자 2026-10-21). E2 결과(D 판 전체/뒤 절반 0.89–1.18, 구간이 1 포함)를 본 뒤 등록.',
  question='round37 맞춘 좌표 판의 "뒤 절반의 온도 차가 전체 창의 1.3–2.1 배" 와 D 판 E2 의 불일치가 거리(D 대 맞춘 좌표) 때문인가, 접힌 프레임 기준(≥ 40 대 ≥ 100)·도메인 집합 때문인가.',
  design='줄기 336 D 파일(r39_D_stem)에 같은 실행·같은 프레임으로 저장된 맞춘 좌표 칸(Aall: 모든 프레임 쌍, Δ²/√변; dmu2_A)으로 r39_core.part_e2 를 그대로 다시 냄(값 키만 D → Aall, dmu2_D → dmu2_A; 같은 seed 20261072, 같은 도메인 집합, 실행당 접힌 프레임 ≥ 100). 출력 results/r39_E2_aligned_diag.json.',
  reading=('맞춘 좌표 판의 뒤 절반/전체 비(= 1/(전체/뒤 절반)) 를 +59·+93·+130 K 에서 본다. 셋 중 둘 이상에서 95% 구간이 1 을 포함하면 "기준·집합" (round37 그 행은 투명성 표로 내림); '
           '셋 중 둘 이상에서 구간 하한 ≥ 1.2 이면 "거리" (실험 2 의 D 압축 기전과 같은 뿌리); 둘 다 아니면 "가르지 못함". +28 K 는 서술.'))
json.dump(P, open(pf + '.tmp', 'w'), ensure_ascii=False, indent=1); os.replace(pf + '.tmp', pf)
b = f"{R}/results/_backup/{time.strftime('%Y-%m-%d')}"; shutil.copy(pf, b + '/' + time.strftime('%H%M%S_') + os.path.basename(pf)); print(P[K]['stamped_utc'])
src = open(f'{C}/r39_core.py').read(); body = src[:src.index("\nif __name__ == '__main__':")]
seg = body[body.index('def part_e2():'):body.index('# ---------- E3')]
seg2 = seg.replace("['D']['d2n']", "['Aall']['d2n']").replace("['dmu2_D']", "['dmu2_A']").replace("RES + '/r39_E2.json'", "RES + '/r39_E2_aligned_diag.json'").replace("log('R39_E2_DONE')", "log('R39_E2_ALIGNED_DIAG_DONE')").replace("'round39 E2 가벼운 판 (D)", "'사후 진단: E2 와 같은 판을 맞�춘 좌표(Aall)로 (D)")
print(seg.count("['D']['d2n']"), seg.count("['dmu2_D']"))
open('/content/r39_e2_aligned_diag.py', 'w').write(body.replace(seg, seg2) + "\npart_e2()\n")
r = subprocess.run(f"cd /content && CUDA_VISIBLE_DEVICES= python3 /content/r39_e2_aligned_diag.py {R} none 2>&1 | grep -v Warning | tail -2 | cut -c1-400", shell=True, capture_output=True, text=True); print(r.stdout)
```

Output:

```
2026-10-21 12:43:00
3 2
12:43:11 R39_E2 {'348': 0.928, '379': 0.917, '413': 0.858, '450': 1.029}
12:43:11 R39_E2_ALIGNED_DIAG_DONE
```

Cell `HpplqbJChcyO`, run requested 12:43:20, output returned 12:43:21 UTC.

```python
fo = f'{R}/results/r39_E2_aligned_diag.json'; d = json.load(open(fo))
d['rule'] = '사후 진단(Drive JSON round39_posthoc_window_diag_2026-10-21): r39_core.part_e2 를 맞춘 좌표(Aall Δ²/√변, dmu2_A)로 그대로 — 같은 도메인 집합·같은 seed·실행당 접힌 프레임 ≥ 100. 뒤 절반/전체 = 1/(전체/뒤 절반).'
res = {}
for T, s in d['steps'].items():
    r = s['ratio_full_over_half']; hf = (1 / r['est'], [1 / r['ci95'][1], 1 / r['ci95'][0]]); res[T] = hf
    print(T, 'n', s['n_dom'], 'Aall d2n 뒤 절반', round(s['d2n_half']['est'], 4), '전체', round(s['d2n_full']['est'], 4), '뒤절반/전체', round(hf[0], 3), [round(x, 3) for x in hf[1]], 'dmu2 전체/뒤절반', round(s['dmu2_ratio_full_over_half']['est'], 3))
ci1 = sum(1 for T in ('379', '413', '450') if res[T][1][0] <= 1 <= res[T][1][1]); lo12 = sum(1 for T in ('379', '413', '450') if res[T][1][0] >= 1.2)
d['reading'] = '기준·집합' if ci1 >= 2 else ('거리' if lo12 >= 2 else '가르지 못함'); d['half_over_full'] = {T: dict(est=v[0], ci95=v[1]) for T, v in res.items()}
json.dump(d, open(fo + '.tmp', 'w'), ensure_ascii=False, indent=1); os.replace(fo + '.tmp', fo); shutil.copy(fo, b + '/' + time.strftime('%H%M%S_') + os.path.basename(fo)); print('읽기', d['reading'], ci1, lo12)
```

Output (per heating step, 348, 379, 413 and 450 K being +28, +59, +93 and +130 K: number of domains, the
aligned-coordinate Δ²/√n over the second half (뒤 절반) and the whole trajectory (전체), the ratio second half /
whole (뒤절반/전체) with its 95% interval, and the whole / second-half ratio of Δμ²₂; last line, the reading (읽기,
"criterion and set") and the numbers of the three steps whose interval contains 1 and whose lower bound is at least
1.2):

```
348 n 264 Aall d2n 뒤 절반 0.0065 전체 0.0061 뒤절반/전체 1.077 [0.488, 2.041] dmu2 전체/뒤절반 0.845
379 n 209 Aall d2n 뒤 절반 0.0256 전체 0.0235 뒤절반/전체 1.091 [0.885, 1.414] dmu2 전체/뒤절반 0.73
413 n 111 Aall d2n 뒤 절반 0.0454 전체 0.0389 뒤절반/전체 1.166 [0.859, 1.533] dmu2 전체/뒤절반 0.912
450 n 47 Aall d2n 뒤 절반 0.0605 전체 0.0623 뒤절반/전체 0.971 [0.916, 1.555] dmu2 전체/뒤절반 0.75
읽기 기준·집합 3 0
```

## 11. Coverage simulations of the first thermophile test (20 October 2026)

Methods, thermophile study: coverage simulations whose output files reached Drive at 12:48:32 and 12:56:47 UTC, before
the decision file (12:57:36), showed that two-sided 95% intervals contained the simulated θ* in 28–70% of runs and the
upper bound alone lay above it in 26–71%, whereas the one-sided lower bound used for the decision lay below it in all
runs. The two-sided coverage is part `sim` of `code/r38_p1.py` (run 12:44:58–12:48:32 UTC; `results/r38_p1_sim.json`,
column `coverage`); the checksum of `r38_p1.py` entered the first record at 04:58:59 UTC (entry
`round36_p1_code_2026-10-20`, stamp 04:58:52), and the script holds the rule that a coverage below 0.90 for θ < 413 K requires the
interval method to be changed before the decision. The same estimator simulation and rule were in the first record
before that, as a row of the freeze checklist in `round36_partial_freeze_7b7c` (19 October 14:56:51 UTC: "P1 estimator
simulation: bias and coverage; below 0.90, change the interval method; bias only described", translated from Korean), and
the same entry fixed the order: the simulation first, its coverage entered in the checklist,
then the decision. The one-sided coverages and a simulation of condition 2 were run by
two short scripts that the cells below wrote to `/content` and started in the background, `/content/p1_sim2.py` and
`/content/p1_sim3.py` (not kept as files; their full text is in the first two cells). Both load the deposited
`r38_p1.py` up to its stage list `ALL1` and call its functions unchanged (`sim_noise`, `sim_once`, `theta`,
`sf_summary`, `cond2`) on the same pairs and replicate noise, with 100 repetitions and 1,000 bootstrap resamples (seeds
20261043 and 20261044). `p1_sim2.py` wrote `results/r38_p1_sim_onesided.json`: for θ = 320, 335, 350 and 380 K, with all
pairs at that θ (`sd` 0) or with a spread of 15 K between pairs (`sd` 15), the share of runs whose bootstrap lower bound
lay at or below θ (`lower_cov`) and whose upper bound lay at or above it (`upper_cov`). `p1_sim3.py` wrote
`results/r38_p1_sim_cond2.json`: the share of runs in which the lower bound of condition 2 exceeded zero
(`frac_lo_gt_0`; `version` 짝지은 판, "paired version"), a false-positive rate when thermophile–mesophile and control
pairs share θ (320 or 350 K) and a power when they differ by 15 or 30 K. No record names these two scripts before they were run. The one-sided coverages carry out the pre-specified rule: the
two-sided coverage had fallen below 0.90, and they chose the changed interval method (a one-sided 97.5% lower bound).
The condition-2 simulation was not specified in advance: the estimator simulation of the records covers the bias and
coverage of θ* (condition 1; part `sim` of `r38_p1.py` simulates the thermophile summaries only), whereas `p1_sim3.py`
estimates the false-positive rate and power of condition 2 and applies no coverage rule; its results were entered in the
record before the decision file was read. The cell run at 12:57:23 entered both summaries, with that of the two-sided simulation, in the entry
`round36_p1_checklist_2026-10-20` of the first record (`sim_onesided`, and `sim_cond2` with the placeholder 돌리는 중,
"running", because the condition-2 file did not yet exist when the cell looked for it); the cell run at 12:57:59
replaced the placeholder with the condition-2 rows (`sim_cond2_added_utc`, 12:58:00 by the Colab clock). S1 Table flags this entry:
it reached Drive 3 s after the decision file. In these cells `R` is the project root, `L = R + '/results/logs'`, `bdir`
the backup folder of the day (`results/_backup/2026-10-20/`) and `sh(c)` runs a shell command and returns the last
1,500 characters of its output; `datetime`, `hashlib` and `shutil` were imported earlier. Times printed inside the
output come from the Colab clock. Cells that only re-ran a log check, and cells about other analyses, are left out.

Cell `H_XE91IsBAsy`, run requested 2026-10-20 12:55:49, output returned 12:55:51 UTC.

```python
SMK4 = r'''import sys, json, time
R = sys.argv[1]
ns = {'__name__': 'p1sim2'}; s = open(R + '/code/r38_p1.py').read(); s = s[:s.index("ALL1 = [") if "ALL1 = [" in s else s.index("if __name__ == '__main__':")]
sys.argv = ['r38_p1', R, 'none', '1']; exec(s, ns)
np = ns['np']; pj = json.load(open(ns['OUT1'] + '/pairs.json')); sp = ns['load_spread'](pj['domains']); nz = [z for z in (ns['sim_noise'](sp, r) for r in pj['HT']) if z is not None]
rng = np.random.default_rng(20261043); key = ns['key_ht']; rows = []
for sdv in (0.0, 15.0):
    for th0 in (320.0, 335.0, 350.0, 380.0):
        lo_ok, hi_ok = [], []
        for _ in range(100):
            med, lo, hi = ns['sim_once'](nz, th0, sdv, rng, B=1000)
            lo_ok.append(True if lo[0] == 'lt' else key(lo) <= (th0, 1)); hi_ok.append(True if hi[0] == 'gt' else (hi[0] != 'lt' and hi[1] >= th0))
        rows.append(dict(theta=th0, sd=sdv, lower_cov=float(np.mean(lo_ok)), upper_cov=float(np.mean(hi_ok))))
        print('P1_SIM2', rows[-1], flush=True)
json.dump(dict(written=time.strftime('%F %T'), rows=rows, note='탐색(판정 전, 실제 θ 값 안 봄): 한쪽 피복률 = P(하한 ≤ θ0), P(상한 ≥ θ0); 반복 100, 부트스트랩 1,000, seed 20261043'), open(R + '/results/r38_p1_sim_onesided.json', 'w'), ensure_ascii=False, indent=1)
print('P1_SIM2_DONE', flush=True)
'''
open('/content/p1_sim2.py', 'w').write(SMK4)
print(sh(f"grep -n \"ALL1 = \\[\\|__main__\" {R}/code/r38_p1.py | head -3; cd /content && CUDA_VISIBLE_DEVICES= setsid nohup nice -n 5 /usr/bin/python3 /content/p1_sim2.py {R} > {L}/r38_p1_sim2.log 2>&1 < /dev/null & sleep 1; echo 띄움"))
```

Output:

```
447:ALL1 = ['pairs', 'mean', 'align', 'spread', 'sim', 'pre', 'stats']
448:if __name__ == '__main__':
띄움
```

Cell `H_XE91IsBAsy`, run requested 12:56:30, output returned 12:56:32 UTC.

```python
SMK5 = r'''import sys, json, time
R = sys.argv[1]
ns = {'__name__': 'p1sim3'}; s = open(R + '/code/r38_p1.py').read(); s = s[:s.index("ALL1 = [")]
sys.argv = ['r38_p1', R, 'none', '1']; exec(s, ns)
np = ns['np']; pj = json.load(open(ns['OUT1'] + '/pairs.json')); sp = ns['load_spread'](pj['domains'])
zH = [z for z in (ns['sim_noise'](sp, r) for r in pj['HT']) if z is not None]; zM = [z for z in (ns['sim_noise'](sp, r) for r in pj['MM']) if z is not None]
b = ns['B_SLOPE']; theta = ns['theta']
def simit(zs, th0, sdv, rng):
    items = {}
    for z in zs:
        th = th0 + (rng.normal(0, sdv) if sdv > 0 else 0.0); er = rng.normal(0, z['se_r'])
        LR = [b * (T - th) + rng.normal(0, s) - er for T, s in zip(z['temps'], z['se_h'])]
        items.setdefault(z['H'], []).append(theta(LR, z['temps'])[:2])
    return items
rng = np.random.default_rng(20261044); rows = []
for thH, thM in ((320.0, 320.0), (350.0, 350.0), (335.0, 320.0), (350.0, 320.0)):
    for sdv in (0.0, 15.0):
        pos, ver = [], []
        for _ in range(100):
            iH, iM = simit(zH, thH, sdv, rng), simit(zM, thM, sdv, rng)
            T_s = {h: ns['sf_summary'](v, 'ht') for h, v in iH.items()}; C_s = {h: ns['sf_summary'](v, 'mm') for h, v in iM.items()}
            c2 = ns['cond2'](T_s, C_s, rng, B=1000); pos.append(c2['passed']); ver.append(c2['version'])
        rows.append(dict(theta_HT=thH, theta_MM=thM, sd=sdv, frac_lo_gt_0=float(np.mean(pos)), version=max(set(ver), key=ver.count)))
        print('P1_SIM3', rows[-1], flush=True)
json.dump(dict(written=time.strftime('%F %T'), n_HT=len(zH), n_MM=len(zM), rows=rows, note='탐색(판정 전, 실제 θ 값 안 봄): 조건 2 의 하한 > 0 비율. θ_HT = θ_MM 이면 거짓 양성률, 다르면 검정력. 반복 100, 부트스트랩 1,000, seed 20261044'), open(R + '/results/r38_p1_sim_cond2.json', 'w'), ensure_ascii=False, indent=1)
print('P1_SIM3_DONE', flush=True)
'''
open('/content/p1_sim3.py', 'w').write(SMK5)
print(sh(f"cd /content && CUDA_VISIBLE_DEVICES= setsid nohup nice -n 5 /usr/bin/python3 /content/p1_sim3.py {R} > {L}/r38_p1_sim3.log 2>&1 < /dev/null & sleep 1; echo 띄움; grep P1_SIM2 {L}/r38_p1_sim2.log"))
```

Output (the one-sided simulation, still running, had printed five of its eight rows):

```
띄움
P1_SIM2 {'theta': 320.0, 'sd': 0.0, 'lower_cov': 1.0, 'upper_cov': 0.71}
P1_SIM2 {'theta': 335.0, 'sd': 0.0, 'lower_cov': 1.0, 'upper_cov': 0.63}
P1_SIM2 {'theta': 350.0, 'sd': 0.0, 'lower_cov': 1.0, 'upper_cov': 0.6}
P1_SIM2 {'theta': 380.0, 'sd': 0.0, 'lower_cov': 1.0, 'upper_cov': 0.31}
P1_SIM2 {'theta': 320.0, 'sd': 15.0, 'lower_cov': 1.0, 'upper_cov': 0.71}
```

Cell `H_XE91IsBAsy`, run requested 12:56:52, output returned 12:56:53 UTC.

```python
print(sh(f"date -u +%T; grep -E 'P1_SIM2|P1_SIM3' {L}/r38_p1_sim2.log {L}/r38_p1_sim3.log | cut -d: -f2- ; grep -E 'Error|Traceback' {L}/r38_p1_sim3.log | head -3"))
```

Output (all eight one-sided rows and the first two condition-2 rows):

```
12:56:53
P1_SIM2 {'theta': 320.0, 'sd': 0.0, 'lower_cov': 1.0, 'upper_cov': 0.71}
P1_SIM2 {'theta': 335.0, 'sd': 0.0, 'lower_cov': 1.0, 'upper_cov': 0.63}
P1_SIM2 {'theta': 350.0, 'sd': 0.0, 'lower_cov': 1.0, 'upper_cov': 0.6}
P1_SIM2 {'theta': 380.0, 'sd': 0.0, 'lower_cov': 1.0, 'upper_cov': 0.31}
P1_SIM2 {'theta': 320.0, 'sd': 15.0, 'lower_cov': 1.0, 'upper_cov': 0.71}
P1_SIM2 {'theta': 335.0, 'sd': 15.0, 'lower_cov': 1.0, 'upper_cov': 0.61}
P1_SIM2 {'theta': 350.0, 'sd': 15.0, 'lower_cov': 1.0, 'upper_cov': 0.59}
P1_SIM2 {'theta': 380.0, 'sd': 15.0, 'lower_cov': 1.0, 'upper_cov': 0.26}
P1_SIM2_DONE
P1_SIM3 {'theta_HT': 320.0, 'theta_MM': 320.0, 'sd': 0.0, 'frac_lo_gt_0': 0.0, 'version': '짝지은 판'}
P1_SIM3 {'theta_HT': 320.0, 'theta_MM': 320.0, 'sd': 15.0, 'frac_lo_gt_0': 0.0, 'version': '짝지은 판'}
```

The same cell, run again: requested 12:57:07, output returned 12:57:08 UTC. Output:

```
12:57:08
P1_SIM2 {'theta': 320.0, 'sd': 0.0, 'lower_cov': 1.0, 'upper_cov': 0.71}
P1_SIM2 {'theta': 335.0, 'sd': 0.0, 'lower_cov': 1.0, 'upper_cov': 0.63}
P1_SIM2 {'theta': 350.0, 'sd': 0.0, 'lower_cov': 1.0, 'upper_cov': 0.6}
P1_SIM2 {'theta': 380.0, 'sd': 0.0, 'lower_cov': 1.0, 'upper_cov': 0.31}
P1_SIM2 {'theta': 320.0, 'sd': 15.0, 'lower_cov': 1.0, 'upper_cov': 0.71}
P1_SIM2 {'theta': 335.0, 'sd': 15.0, 'lower_cov': 1.0, 'upper_cov': 0.61}
P1_SIM2 {'theta': 350.0, 'sd': 15.0, 'lower_cov': 1.0, 'upper_cov': 0.59}
P1_SIM2 {'theta': 380.0, 'sd': 15.0, 'lower_cov': 1.0, 'upper_cov': 0.26}
P1_SIM2_DONE
P1_SIM3 {'theta_HT': 320.0, 'theta_MM': 320.0, 'sd': 0.0, 'frac_lo_gt_0': 0.0, 'version': '짝지은 판'}
P1_SIM3 {'theta_HT': 320.0, 'theta_MM': 320.0, 'sd': 15.0, 'frac_lo_gt_0': 0.0, 'version': '짝지은 판'}
P1_SIM3 {'theta_HT': 350.0, 'theta_MM': 350.0, 'sd': 0.0, 'frac_lo_gt_0': 0.0, 'version': '짝지은 판'}
P1_SIM3 {'theta_HT': 350.0, 'theta_MM': 350.0, 'sd': 15.0, 'frac_lo_gt_0': 0.0, 'version': '짝지은 판'}
P1_SIM3 {'theta_HT': 335.0, 'theta_MM': 320.0, 'sd': 0.0, 'frac_lo_gt_0': 0.94, 'version': '짝지은 판'}
```

Cell `H_XE91IsBAsy`, run requested 12:57:23, output returned 12:57:26 UTC.

```python
print(sh(f"grep -E 'P1_SIM3' {L}/r38_p1_sim3.log | tail -4"))
import os
S1 = json.load(open(f'{R}/results/r38_p1_sim.json')); S2 = json.load(open(f'{R}/results/r38_p1_sim_onesided.json')); S3 = json.load(open(f'{R}/results/r38_p1_sim_cond2.json')) if os.path.exists(f'{R}/results/r38_p1_sim_cond2.json') else None
PR = json.load(open(f'{R}/results/r38_p1_pre.json'))
PJ = R + '/results/prereg_2026-10-17.json'; P = json.load(open(PJ)); shutil.copy(PJ, bdir + '/' + time.strftime('%H%M%S_') + 'prereg_2026-10-17.json')
k = 'round36_p1_checklist_2026-10-20'; assert k not in P
P[k] = dict(stamped_utc=datetime.datetime.now(datetime.UTC).strftime('%F %T'),
    sim=dict(file='results/r38_p1_sim.json', rows=S1['rows'], coverage_min=S1['coverage_min'], se_h_median=S1['se_h_median'], se_r_median=S1['se_r_median'], n_pairs=S1['n_pairs'], n_sf=S1['n_sf']),
    sim_onesided=S2['rows'], sim_cond2=(S3['rows'] if S3 else '돌리는 중'),
    coverage_rule='prereg 7 절: 피복률(θ < 413) < 0.90 이면 구간 방법을 고친 뒤 동결. 양쪽 피복률 최소 0.275 → 고침.',
    interval_fix='조건 1·2 가 쓰는 것은 부트스트랩 2.5% 하한뿐이다. 그 하한의 한쪽 피복률은 θ < 413 의 모든 시나리오(sd 0·15)에서 1.00, θ0 = 320 의 하한 > 320 비율 0, θ0 = 335 에서 0.88–0.98 (조건 2: θ_HT = θ_MM 이면 하한 > 0 비율 0, 15 K 차에서 0.94). 양쪽 피복이 깨지는 것은 위쪽 끝이다(0.26–0.71): 등록된 보수적 순위(상과 아래 중앙값, "> u" 를 u 자리)가 추정량을 정의상 아래로 민다(치우침 −3 ~ −10 K, θ 335–380). 그래서 보고 구간을 "한쪽 97.5% 하한" 으로 바꾸고, 위쪽 끝은 신뢰 한계로 쓰지 않는다. 판정 수치(하한 > 320, 하한 > 0)는 그대로이고 코드는 바꾸지 않는다. 점추정은 시뮬레이션 치우침 표와 함께 서술한다(판정에서 빼지 않음).',
    pre=dict(file='results/r38_p1_pre.json', n_pairs=PR['n_pairs'], n_used=PR['n_used'], n_both_numeric=PR['n_both_numeric'], n_sf=PR['n_sf'], diff_full_minus_gpa=PR['diff_full_minus_gpa'], censor_agree=PR['censor_agree'], note='기록이지 판정·규칙 조정에 쓰지 않음'),
    stats_script_md5=hashlib.md5(open(f'{R}/code/r38_p1.py', 'rb').read()).hexdigest(), next='이 도장 뒤 r38_p1.py stats (판정).')
json.dump(P, open(PJ + '.tmp', 'w'), ensure_ascii=False, indent=1); os.replace(PJ + '.tmp', PJ); print('도장', P[k]['stamped_utc'], P[k]['stats_script_md5'], 'cond2 sim', 'ok' if S3 else '돌리는 중')
```

Output (the last four condition-2 rows printed so far, then the stamp of the entry, the md5 of `r38_p1.py` and the state of the condition-2 file, 돌리는 중 "running"):

```
P1_SIM3 {'theta_HT': 350.0, 'theta_MM': 350.0, 'sd': 15.0, 'frac_lo_gt_0': 0.0, 'version': '짝지은 판'}
P1_SIM3 {'theta_HT': 335.0, 'theta_MM': 320.0, 'sd': 0.0, 'frac_lo_gt_0': 0.94, 'version': '짝지은 판'}
P1_SIM3 {'theta_HT': 335.0, 'theta_MM': 320.0, 'sd': 15.0, 'frac_lo_gt_0': 0.76, 'version': '짝지은 판'}
P1_SIM3 {'theta_HT': 350.0, 'theta_MM': 320.0, 'sd': 0.0, 'frac_lo_gt_0': 1.0, 'version': '짝지은 판'}
도장 2026-10-20 12:57:26 58c1062c13ab8afc72c35d809eef9500 cond2 sim 돌리는 중
```

Cell `H_XE91IsBAsy`, run requested 12:57:51, output returned 12:57:52 UTC.

```python
print(sh(f"date -u +%T; grep -E 'P1_SIM3' {L}/r38_p1_sim3.log | tail -2; grep -v -E 'WARNING|I0000|E0000|W0000' {L}/r38_p1.log | tail -3 | cut -c1-300"))
```

Output (the last condition-2 row and the end of the decision run, which had been started at 12:57:31):

```
12:57:53
P1_SIM3 {'theta_HT': 350.0, 'theta_MM': 320.0, 'sd': 15.0, 'frac_lo_gt_0': 1.0, 'version': '짝지은 판'}
P1_SIM3_DONE
12:57:36 P1_STATS_DONE
12:57:36 부분 끝 stats
12:57:36 P1_DONE stats
```

Cell `H_XE91IsBAsy`, run requested 12:57:59, output returned 12:57:59 UTC.

```python
P = json.load(open(PJ)); S3 = json.load(open(f'{R}/results/r38_p1_sim_cond2.json'))
P['round36_p1_checklist_2026-10-20']['sim_cond2'] = S3['rows']; P['round36_p1_checklist_2026-10-20']['sim_cond2_added_utc'] = datetime.datetime.now(datetime.UTC).strftime('%F %T') + ' (r38_p1.json 을 읽기 전)'
json.dump(P, open(PJ + '.tmp', 'w'), ensure_ascii=False, indent=1); os.replace(PJ + '.tmp', PJ); print('조건 2 시뮬레이션 추가')
J = json.load(open(f'{R}/results/r38_p1.json')); shutil.copy(f'{R}/results/r38_p1.json', bdir + '/' + time.strftime('%H%M%S_') + 'r38_p1.json')
print(json.dumps(J['main'], ensure_ascii=False, indent=0)[:4000])
```

Output (first line, "condition-2 simulation added"; the decision values of `results/r38_p1.json` that the cell printed next are left out):

```
조건 2 시뮬레이션 추가
```

## 12. Diagnostics of the failed truncated-chain check (19 October 2026)

S2 Text 5: diagnostics on four of the domains (run the same day, not specified in advance) located the failure of the
pre-specified truncated-chain check in its assumption that the removed bonds are matched in full-versus-full
comparisons and unmatched in full-versus-truncated ones. The cells below wrote two short scripts to `/content` and ran
them, `/content/s0_trunc_diag.py` and `/content/s0_trunc_exact.py` (not kept as files; their full text is in the
cells). Both load the deposited `code/r38_s0.py` up to its stage list (`ALL = ['pick'`) and use its functions unchanged
(`dom_runs`, `coords`, `MeanStart`, `upairs`, `red_pairs`, `seeded`; `CUT` = 10 removed bonds), for 1ileA03, 2l23A00,
1ffsB00 and 3r1kA02. `s0_trunc_diag.py` compares 300 frame pairs between replicates 0 and 1, full against full (`ff`)
and full against truncated (`ft`), aligned from the mean-point start with the full dynamic program, and wrote
`results/r38_s0_trunc_diag.json`. In its output, one minus `frac_x_cut_edges_matched` under `ff` is the share of the
removed segment's bonds left unmatched between two full frames (0.435 and 0.338 in 1ileA03 and 2l23A00, the "43% and
34%" of S2 Text, and 0.015 in 1ffsB00; `ff_cut_unmatched_frac_x`, the same share of the squared norm, is 0.435, 0.337
and 0.015), and `frac_x_cut_edges_matched` under `ft` is the share of the removed segment's bonds
paired with other bonds of the truncated chain (0.325, 0.085 and 0.003, the "33% and 8%" and "0.3%").
`s0_trunc_exact.py` compares each frame with its own truncation, by the reduced computation and by the 128-start full
computation, against 0.75 Σ s² of the removed bonds; its first run, on 40 frames, stopped at the cell's 300 s time
limit, and it was restarted in the background on 12 frames (the full computation on 4 of them), which wrote
`results/r38_s0_trunc_exact.json`. In its output, `max_abs_d2_minus_cut_reduced` and `max_abs_d2_minus_cut_full128` are
the largest |d² − 0.75 Σ s²| over the frames: at most 1.7 × 10⁻¹³ in three domains, and 0.0018 in 1ileA03 by both
computations. In these cells `R` is the project root, `L = R + '/results/logs'`, `Dd = R + '/derived/r38'` and `sh(c)`
runs a shell command and returns the last 1,500 characters of its output. Times printed inside the output come from the
Colab clock. Progress checks in between that printed only the time, and cells about other analyses, are left out.

Cell `vfgSWJKreH1m`, run requested 2026-10-19 16:38:35, output returned 16:38:36 UTC.

```python
diag=r'''import sys, os, json, time
ROOT=sys.argv[1]
src=open(ROOT+'/code/r38_s0.py').read(); src=src[:src.index("ALL = ['pick'")]
sys.argv=['x',ROOT,'none']; exec(src)
torch.set_num_threads(4)
out={}
for d in ('1ileA03','2l23A00','1ffsB00','3r1kA02'):
    _runs.clear(); R_,L_=dom_runs(d); ne=L_-1; Ef=coords(np.concatenate([X for _,X in R_])); Et=Ef[:,:ne-CUT].contiguous()
    nA=(Ef**2).sum((-1,-2)); rng=seeded('diag',d); res={}
    for tag,EB,ms in (('ff',Ef,MeanStart(Ef,Ef)),('ft',Et,MeanStart(Ef,Et))):
        ia,ib=upairs(0,1,rng,k=300); ia=torch.as_tensor(ia); ib=torch.as_tensor(ib)
        al=fs.align(Ef[ia],EB[ib],starts=ms.starts(ia,ib,n_local=1),iters=10); M=al.match
        nB=(EB**2).sum((-1,-2)); mx=M.any(2); my=M.any(1); cut=torch.arange(ne)>=ne-CUT
        gx_cut=(nA[ia]*(~mx)*cut).sum(1); gx_non=(nA[ia]*(~mx)*(~cut)).sum(1)
        cutB=torch.arange(EB.shape[1])>=ne-CUT if tag=='ff' else torch.zeros(EB.shape[1],dtype=torch.bool)
        gy_cut=(nB[ib]*(~my)*cutB).sum(1); gy_non=(nB[ib]*(~my)*(~cutB)).sum(1)
        pi,ii,jj=torch.nonzero(M,as_tuple=True); shift=float((ii!=jj).float().mean())
        xcut_matched=float((mx[:,ne-CUT:]).float().mean())
        res[tag]=dict(gap_x_cut=float(gx_cut.mean()),gap_x_non=float(gx_non.mean()),gap_y_cut=float(gy_cut.mean()),gap_y_non=float(gy_non.mean()),
                      frac_x_cut_edges_matched=xcut_matched,frac_matches_off_diagonal=shift)
    s2cut=float((Ef[:,ne-CUT:]**2).sum((-1,-2,-3)).mean())
    ff,ft=res['ff'],res['ft']
    res['expect_075_sum_s2_cut']=s2cut
    res['raw_gap_diff']=(ft['gap_x_cut']+ft['gap_x_non']+ft['gap_y_non'])-(ff['gap_x_cut']+ff['gap_x_non']+ff['gap_y_cut']+ff['gap_y_non'])
    res['noncut_gap_diff']=(ft['gap_x_non']+ft['gap_y_non'])-(ff['gap_x_non']+ff['gap_y_non'])
    res['ft_cut_gap_vs_expect']=ft['gap_x_cut']/s2cut
    res['ff_cut_unmatched_frac_x']=ff['gap_x_cut']/s2cut
    out[d]=res; print(d, json.dumps({k:(round(v,4) if isinstance(v,float) else {kk:round(vv,4) for kk,vv in v.items()}) for k,v in res.items()}, ensure_ascii=False), flush=True)
json.dump(out,open(ROOT+'/results/r38_s0_trunc_diag.json','w'),ensure_ascii=False,indent=1); print('DIAG_DONE')
'''
open("/content/s0_trunc_diag.py","w").write(diag)
subprocess.run(f"cd /content && CUDA_VISIBLE_DEVICES= nohup nice -n 5 /usr/bin/python3 /content/s0_trunc_diag.py {R} > {L}/r38_s0_trunc_diag.log 2>&1 &", shell=True); print("띄움", sh("date -u +%T"))
```

Output ("started"):

```
띄움 16:38:37
```

Cell `vfgSWJKreH1m`, run requested 16:39:53, output returned 16:41:08 UTC.

```python
import time; time.sleep(75); print(sh(f"date -u +%T; grep -v 'oneDNN\\|absl\\|cudart_stub\\|WARNING: All' {L}/r38_s0_trunc_diag.log | tail -6 | cut -c1-900; ps -eo etime,pcpu,args | grep '[s]0_trunc_diag' | cut -c1-80"))
```

Output (the four domains of the diagnostic, then `DIAG_DONE`; the last line is the shell command itself, listed by `ps`):

```
16:41:09
1ileA03 {"ff": {"gap_x_cut": 3.2679, "gap_x_non": 8.9956, "gap_y_cut": 3.2345, "gap_y_non": 9.0167, "frac_x_cut_edges_matched": 0.5653, "frac_matches_off_diagonal": 0.1546}, "ft": {"gap_x_cut": 5.0738, "gap_x_non": 10.3779, "gap_y_cut": 0.0, "gap_y_non": 7.9446, "frac_x_cut_edges_matched": 0.3253, "frac_matches_off_diagonal": 0.1363}, "expect_075_sum_s2_cut": 7.5114, "raw_gap_diff": -1.1184, "noncut_gap_diff": 0.3102, "ft_cut_gap_vs_expect": 0.6755, "ff_cut_unmatched_frac_x": 0.4351}
2l23A00 {"ff": {"gap_x_cut": 2.5385, "gap_x_non": 7.8471, "gap_y_cut": 2.5329, "gap_y_non": 7.8505, "frac_x_cut_edges_matched": 0.662, "frac_matches_off_diagonal": 0.1305}, "ft": {"gap_x_cut": 6.8828, "gap_x_non": 8.5009, "gap_y_cut": 0.0, "gap_y_non": 7.8548, "frac_x_cut_edges_matched": 0.0847, "frac_matches_off_diagonal": 0.1184}, "expect_075_sum_s2_cut": 7.5246, "raw_gap_diff": 2.4696, "noncut_gap_diff": 0.6582, "ft_cut_gap_vs_expect": 0.9147, "ff_cut_unmatched_frac_x": 0.3374}
1ffsB00 {"ff": {"gap_x_cut": 0.1098, "gap_x_non": 1.6136, "gap_y_cut": 0.1093, "gap_y_non": 1.6138, "frac_x_cut_edges_matched": 0.9853, "frac_matches_off_diagonal": 0.0154}, "ft": {"gap_x_cut": 7.5157, "gap_x_non": 1.7076, "gap_y_cut": 0.0, "gap_y_non": 1.6865, "frac_x_cut_edges_matched": 0.003, "frac_matches_off_diagonal": 0.0175}, "expect_075_sum_s2_cut": 7.5256, "raw_gap_diff": 7.4633, "noncut_gap_diff": 0.1667, "ft_cut_gap_vs_expect": 0.9987, "ff_cut_unmatched_frac_x": 0.0146}
3r1kA02 {"ff": {"gap_x_cut": 0.602, "gap_x_non": 5.658, "gap_y_cut": 0.6017, "gap_y_non": 5.6491, "frac_x_cut_edges_matched": 0.9197, "frac_matches_off_diagonal": 0.0682}, "ft": {"gap_x_cut": 6.8179, "gap_x_non": 6.3731, "gap_y_cut": 0.0, "gap_y_non": 5.687, "frac_x_cut_edges_matched": 0.09, "frac_matches_off_diagonal": 0.073}, "expect_075_sum_s2_cut": 7.4916, "raw_gap_diff": 6.3673, "noncut_gap_diff": 0.753, "ft_cut_gap_vs_expect": 0.9101, "ff_cut_unmatched_frac_x": 0.0804}
DIAG_DONE
      00:00  0.0 /bin/sh -c date -u +%T; grep -v 'oneDNN\|absl\|cudart_stub\|WAR
```

Cell `vfgSWJKreH1m`, run requested 16:41:42, output returned 16:46:43 UTC.

```python
ex=r'''import sys, json
ROOT=sys.argv[1]
src=open(ROOT+'/code/r38_s0.py').read(); src=src[:src.index("ALL = ['pick'")]
sys.argv=['x',ROOT,'none']; exec(src)
torch.set_num_threads(4); out={}
for d in ('1ileA03','2l23A00','1ffsB00','3r1kA02'):
    _runs.clear(); R_,L_=dom_runs(d); ne=L_-1; Ef=coords(np.concatenate([X for _,X in R_]))[:40]; Et=Ef[:,:ne-CUT].contiguous()
    cut=(Ef[:,ne-CUT:]**2).sum((-1,-2,-3))
    ms=MeanStart(Ef,Et); idx=np.arange(40); dr,gap,S,ku,_=red_pairs(ms,idx,idx,want_gap=True)
    dfull=torch.stack([fs.distance(Ef[k],Et[k])[0] for k in range(40)])
    out[d]=dict(max_abs_d2_minus_cut_reduced=float(np.abs(dr**2-cut.numpy()).max()), max_abs_gap_minus_cut=float(np.abs(gap-cut.numpy()).max()),
                max_abs_d2_minus_cut_full128=float((dfull**2-cut).abs().max()), mean_cut=float(cut.mean()))
    print(d, json.dumps(out[d]), flush=True)
json.dump(out,open(ROOT+'/results/r38_s0_trunc_exact.json','w'),indent=1); print('EXACT_DONE')
'''
open("/content/s0_trunc_exact.py","w").write(ex)
r=subprocess.run(f"cd /content && CUDA_VISIBLE_DEVICES= /usr/bin/python3 /content/s0_trunc_exact.py {R} 2>&1 | grep -v 'oneDNN\\|absl\\|cudart_stub\\|WARNING: All'", shell=True, capture_output=True, text=True, timeout=300); print(r.stdout[-1500:])
```

Output: the error below (traceback left out); the next cell stopped the script and restarted it.

```
TimeoutExpired: Command 'cd /content && CUDA_VISIBLE_DEVICES= /usr/bin/python3 /content/s0_trunc_exact.py /content/drive/MyDrive/atlas_framed_curves 2>&1 | grep -v 'oneDNN\|absl\|cudart_stub\|WARNING: All'' timed out after 300 seconds
```

Cell `vfgSWJKreH1m`, run requested 16:46:57, output returned 16:46:58 UTC.

```python
import signal
for l in sh("ps -eo pid,args").splitlines()[1:]:
    pid,args=l.strip().split(None,1)
    if 's0_trunc_exact.py' in args and 'ps -eo' not in args and int(pid)!=os.getpid(): os.kill(int(pid),signal.SIGTERM); print("멈춤", pid, args[:60])
ex=open("/content/s0_trunc_exact.py").read().replace("[:40]","[:12]").replace("np.arange(40)","np.arange(12)").replace("for k in range(40)","for k in range(4)").replace("(dfull**2-cut).abs().max()","(dfull**2-cut[:4]).abs().max()")
open("/content/s0_trunc_exact.py","w").write(ex)
subprocess.run(f"cd /content && CUDA_VISIBLE_DEVICES= nohup /usr/bin/python3 /content/s0_trunc_exact.py {R} > {L}/r38_s0_trunc_exact.log 2>&1 &", shell=True); print("뒤에서 띄움", sh("date -u +%T"))
```

Output ("stopped" with the process number, then "started in the background"):

```
멈춤 52939 /usr/bin/python3 /content/s0_trunc_exact.py /content/drive/M
뒤에서 띄움 16:46:59
```

Cell `vfgSWJKreH1m`, run requested 16:56:09, output returned 16:56:10 UTC.

```python
print(sh(f"date -u +%T; grep -v 'oneDNN\\|absl\\|cudart_stub\\|WARNING: All' {L}/r38_s0_trunc_exact.log | tail -6 | cut -c1-400; wc -l < {Dd}/s0/bias.jsonl; ps -eo pid,etime,args | grep -E '[r]38_s0.py|[q]_full_gpu2' | cut -c1-120"))
```

Output (the four domains of the identity check, then `EXACT_DONE`; the last three lines are the number of finished lines of the bias check of `r38_s0.py`, which ran at the same time, and the running processes):

```
16:56:11
1ileA03 {"max_abs_d2_minus_cut_reduced": 0.0017919650054656344, "max_abs_gap_minus_cut": 0.02284898063445784, "max_abs_d2_minus_cut_full128": 0.0017919650054656344, "mean_cut": 7.545062078723937}
2l23A00 {"max_abs_d2_minus_cut_reduced": 1.5187850976872141e-13, "max_abs_gap_minus_cut": 1.7763568394002505e-15, "max_abs_d2_minus_cut_full128": 1.5631940186722204e-13, "mean_cut": 7.520251467780422}
1ffsB00 {"max_abs_d2_minus_cut_reduced": 4.529709940470639e-14, "max_abs_gap_minus_cut": 1.7763568394002505e-15, "max_abs_d2_minus_cut_full128": 2.3092638912203256e-14, "mean_cut": 7.4997785403434785}
3r1kA02 {"max_abs_d2_minus_cut_reduced": 1.6697754290362354e-13, "max_abs_gap_minus_cut": 1.7763568394002505e-15, "max_abs_d2_minus_cut_full128": 1.6697754290362354e-13, "mean_cut": 7.497694505237564}
EXACT_DONE
568
  16098    02:28:03 bash /content/q_full_gpu2.sh
  44493       46:32 /usr/bin/python3 /content/drive/MyDrive/atlas_framed_curves/code/r38_s0.py /content/drive/MyDrive/at
```

## 13. The 450 K version of the temperature conversion without a missing value (21 October 2026)

Results (calibration) and S1 Table: the 450 K version of the temperature conversion (E4) is among the pre-specified
predictions of the calibration that were not met. It requires the superfamily-median temperature equivalent at 450 K
to lie within 10 K, which it did in both windows (bias −5.3 K over second halves and −9.3 K over whole trajectories,
first cell below), and the mean absolute error (MAE) of the domains with a 450 K step to be at most 1.5 times the MAE
at 413 K of the same domains. Part `e4` of `code/r39_core.py` computed the latter with
`np.mean`, so over second halves the missing value of a domain whose fitted amplitude was not positive made it missing
(`MAE_413_same_set` `nan` in `results/r39_E4.json`, first cell below), and the 450 K version was marked as not met
without a comparison. The second cell checks that the 31 domains of the 450 K version (72 over whole trajectories) are
all among those of the 413 K version and that the stored rows do not keep the per-domain errors (`'ae' in r413` is
`False`), so the comparison could not be redone from the stored file. Not specified in advance: the third cell wrote a
copy of the deposited `r39_core.py` (md5 `04e57acc…`) up to its `__main__` block, with that one `np.mean` replaced by
`np.nanmean`, the output file renamed `results/r39_E4_nancheck.json` and the end marker renamed, to
`/content/r39_e4_nancheck.py` (not kept; the first output line gives the number of occurrences of the two replaced
strings, 1 and 1, as in the deposited file), followed by a call of `part_e4()`, and ran it. Without the missing value,
the 413 K MAE of the same domains was 41.3 K over second halves (29.3 K over whole trajectories, unchanged), against a
450 K MAE of 118.7 K (114.9 K), so the 450 K version was not met in either window and the last reliable point stayed
+93 K. The same day (12:44 UTC, session log) the check was entered in the transparency table of the
pre-specification document; no record entry names it. In these cells `R` is the project root, `C = R + '/code'`,
`J(name)` reads `results/<name>.json` and `sh(c)` runs a shell command; times printed inside the output come from the
Colab clock. In the first output, `K̂_sf` is the superfamily-median temperature equivalent with its interval, 편향 the
bias, `MAE` the domain MAE, 90%폭 the width of the 90% range of domain errors, 0이하 the share of domains with a
non-positive value, 외삽 the share extrapolated, 평평 the share on a flat curve and 나눔SD the SD between splits.

Cell `HpplqbJChcyO`, run requested 2026-10-21 12:37:08, output returned 12:37:08 UTC.

```python
e4 = J('r39_E4'); print(e4['checks'])
for W in ('half', 'full'):
    print(W, {k: (v if not isinstance(v, dict) else {a: round(b, 1) if isinstance(b, float) else b for a, b in v.items()}) for k, v in e4['verdict'][W].items()})
    for r in e4['rows']:
        if r['window'] == W and r['curve'] == '모양×진폭':
            print('   ', r['T'], 'K', r['K'], 'n', r['n_dom'], r['n_sf'], 'K̂_sf', r['K_sf'], '편향', round(r['bias'], 1), 'MAE', round(r['dom_MAE'], 1), '90%폭', round(r['dom_err_90w'], 1) if r['dom_err_90w'] else None, '0이하', round(r['frac_zero'] or 0, 2), '외삽', round(r['frac_extrap'] or 0, 2), '평평', round(r['frac_flat'] or 0, 2), '나눔SD', round(r['split_sd'] or 0, 1))
sc = J('r39_sec6') if os.path.exists(f'{R}/results/r39_sec6.json') else None; print(sh(f"ls {R}/results | grep -i sec6"))
```

Output:

```
{'subset_all_equals_cell|half': 0.0, 'subset_all_equals_cell|full': 0.0}
half {'MAE_379_413': {'모양×진폭': 45.8, '선형 보간': 61.5, '모집단': 56.1}, 'main_curve': '모양×진폭', 'replaced': False, 'gate_le10K': {'379': False, '413': True}, 'gate_pass': False, 'pass_450': False, 'MAE_413_same_set': nan, 'last_trusted': '+93', 'pred_order': False}
    348 K 28 n 205 55 K̂_sf {'est': 18.107003997355243, 'ci95': [13.122597019779391, 29.820719989139103]} 편향 -9.9 MAE 31.5 90%폭 113.9 0이하 0.38 외삽 0.03 평평 0 나눔SD 19.8
    379 K 59 n 185 52 K̂_sf {'est': 46.23915059457419, 'ci95': [35.43057652149916, 53.31541799949837]} 편향 -12.8 MAE 37.0 90%폭 133.7 0이하 0.12 외삽 0.05 평평 0 나눔SD 17.8
    413 K 93 n 90 38 K̂_sf {'est': 88.36039583018639, 'ci95': [71.83363875066077, 99.0826684486648]} 편향 -4.6 MAE 54.6 90%폭 141.3 0이하 0.02 외삽 0.17 평평 0 나눔SD 17.4
    450 K 130 n 31 20 K̂_sf {'est': 124.73712953958031, 'ci95': [97.6765402106939, 167.89394184113704]} 편향 -5.3 MAE 118.7 90%폭 546.6 0이하 0 외삽 0.7 평평 0 나눔SD 25.2
full {'MAE_379_413': {'모양×진폭': 41.7, '선형 보간': 48.9, '모집단': 48.8}, 'main_curve': '모양×진폭', 'replaced': False, 'gate_le10K': {'379': True, '413': True}, 'gate_pass': True, 'pass_450': False, 'MAE_413_same_set': 29.3360069477329, 'last_trusted': '+93', 'pred_order': False}
    348 K 28 n 262 62 K̂_sf {'est': 17.18134323226584, 'ci95': [12.765231509275136, 25.441141053668478]} 편향 -10.8 MAE 29.8 90%폭 66.1 0이하 0.33 외삽 0.02 평평 0 나눔SD 18.9
    379 K 59 n 244 58 K̂_sf {'est': 56.00096733809558, 'ci95': [48.31333021920029, 61.393768078647724]} 편향 -3.0 MAE 46.1 90%폭 123.1 0이하 0.07 외삽 0.05 평평 0 나눔SD 23.9
    413 K 93 n 173 51 K̂_sf {'est': 85.3279836047711, 'ci95': [78.67082654177088, 99.99201784981929]} 편향 -7.7 MAE 37.2 90%폭 141.3 0이하 0.01 외삽 0.16 평평 0 나눔SD 14.6
    450 K 130 n 72 34 K̂_sf {'est': 120.69139122904228, 'ci95': [112.56903632317227, 146.3795841367003]} 편향 -9.3 MAE 114.9 90%폭 296.4 0이하 0 외삽 0.82 평평 0 나눔SD 60.8
r39_sec6.json
```

Cell `HpplqbJChcyO`, run requested 12:37:24, output returned 12:37:25 UTC.

```python
for W in ('half', 'full'):
    r450 = next(r for r in e4['rows'] if r['window'] == W and r['T'] == '450' and r['curve'] == '모양×진폭'); r413 = next(r for r in e4['rows'] if r['window'] == W and r['T'] == '413' and r['curve'] == '모양×진폭')
    print(W, len(r450['doms']), len(r413['doms']), len(set(r450['doms']) & set(r413['doms'])), 'ae' in r413)
```

Output (window, number of domains of the 450 K and 413 K versions, number shared, whether the stored 413 K row keeps the per-domain errors):

```
half 31 90 31 False
full 72 173 72 False
```

Cell `HpplqbJChcyO`, run requested 12:37:44, output returned 12:37:48 UTC.

```python
src = open(f'{C}/r39_core.py').read(); a = "m413 = np.mean([a for d, a in zip(r413['doms'], r413['ae']) if d in s450])"
print(src.count(a), src.count("RES + '/r39_E4.json'"))
chk = src[:src.index("\nif __name__ == '__main__':")].replace(a, a.replace('np.mean(', 'np.nanmean(')).replace("RES + '/r39_E4.json'", "RES + '/r39_E4_nancheck.json'").replace("log('R39_E4_DONE')", "log('R39_E4_NANCHECK_DONE')")
open('/content/r39_e4_nancheck.py', 'w').write(chk + "\npart_e4()\n")
r = subprocess.run(f"cd /content && CUDA_VISIBLE_DEVICES= python3 /content/r39_e4_nancheck.py {R} none 2>&1 | grep -v Warning | tail -2 | cut -c1-600", shell=True, capture_output=True, text=True); print(r.stdout)
v = json.load(open(f'{R}/results/r39_E4_nancheck.json'))['verdict']; print({W: (round(v[W]['MAE_413_same_set'], 1), v[W]['pass_450'], v[W]['last_trusted']) for W in v})
```

Output (the last line gives, per window, the 413 K MAE of the same domains, whether the 450 K version passed and the last reliable point):

```
1 1
12:37:49 R39_E4 {'half': {'MAE_379_413': {'모양×진폭': np.float64(45.81292734841799), '선형 보간': np.float64(61.527850748759874), '모집단': np.float64(56.07625252154164)}, 'main_curve': '모양×진폭', 'replaced': False, 'gate_le10K': {'379': False, '413': True}, 'gate_pass': False, 'pass_450': False, 'MAE_413_same_set': 41.346274551440565, 'last_trusted': '+93', 'pred_order': False}, 'full': {'MAE_379_413': {'모양×진폭': np.float64(41.67037249202508), '선형 보간': np.float64(48.86740875670304), '모집단': np.float64(48.78136981335953)}, 'main_curve': '모양×진폭'
12:37:49 R39_E4_NANCHECK_DONE

{'half': (41.3, False, '+93'), 'full': (29.3, False, '+93')}
```

## 14. Accuracy of the 128-start values after the specification error of the certification (20 October 2026)

S2 Text 1: the second prediction of the certification was specified (record entry `round36_sec6_2026-10-20`, on Drive
at 17:34:46 UTC) as (S_ub − S_heur)/S_ub ≤ 0.002 in at least 95% of pairs, S_ub being the certified upper bound and
S_heur the 128-start value. After the result of the static comparison (`results/r38_static.json`, written 18:56:22),
the first cell below read it and the per-pair file `derived/r38/static/bnb.jsonl` (500 pairs): all were certified within
60 s, and the specified quantity was at most 0.002 in none. The second cell counted the pairs in which branch-and-bound
found a better alignment than the 128 starts (35) and the relative difference in d this shows (at most 0.0036, the
0.36% of S2 Text 1). A cell that printed the head of `framespace_pkg/framespace/certify.py`, where the tolerance is
rel_tol × 0.75 × max(n, m) and S_ub is the best value found plus the tolerance, is left out. Not specified in advance:
the third cell measured the 128-start value against the best value found, (S_lb − S_heur)/(0.75 max(n, m)) ≤ 0.002,
which held in 499 of 500 pairs (largest shortfall 0.0032 in that unit), and confirmed that S_ub − S_lb equals the
tolerance (0.005 in that unit), so the specified criterion could not hold. The fourth cell added this description to
`code/r39_sec6.py` through the helper `patch` (defined earlier in the runtime: it writes the patched script to Drive
and prints whether its md5 matches the one computed in the cell; 일치, "match") and entered the finding and the change
as `fix1` of the record entry `round36_sec6_2026-10-20` (stamp 19:06:04, on Drive at 19:06:13). Part `stats` of the
patched `r39_sec6.py` (md5 `1b11b675…`, the deposited file) recomputed the description on 21 October
(`bnb_posthoc` in `results/r39_sec6.json`: 0.998, largest shortfall 0.0032, 35 pairs). In these cells `R` is the project
root, `PJ` the path of the first record and `sh(c)` runs a shell command; `json`, `hashlib`, `datetime` and `np` were
imported earlier. In the output, 인증 is certified, 유사도 틈 the specified relative gap, 휴리스틱 the heuristic,
하한 the lower bound and 출처별 by source (H, pairs from the homolog pool; O, other pairs).

Cell `H_XE91IsBAsy`, run requested 2026-10-20 19:04:49, output returned 19:04:50 UTC.

```python
js = json.load(open(f'{R}/results/r38_static.json')); print({k: js[k] for k in ('written', 'n_H', 'n_O', 'n_O_registered')}); print(json.dumps(js['bnb'], ensure_ascii=False))
for k, v in js['by_level_TM'].items(): print(k, v)
bn = [json.loads(l) for l in open(f'{R}/derived/r38/static/bnb.jsonl')]; ok = [r for r in bn if 'err' not in r]
gap = np.array([(r['S_ub'] - r['S_heur']) / r['S_ub'] for r in ok]); sec = np.array([r['sec'] for r in ok]); rel = np.array([r['rel_err_ub'] for r in ok])
print('인증', sum(r['certified'] for r in ok), '/', len(ok), '| ≤60 s 인증', sum(r['certified'] and r['sec'] <= 60 for r in ok), '| 유사도 틈 ≤ 0.002', int((gap <= 0.002).sum()), '중앙', np.median(gap), '최대', gap.max(), '| sec 중앙', np.median(sec), '90%', np.percentile(sec, 90), '최대', sec.max())
print('휴리스틱 > 하한(S_heur ≥ S_lb):', sum(r['S_heur'] >= r['S_lb'] - 1e-9 for r in ok), '| BnB 가 휴리스틱보다 높은 유사도를 찾음(S_lb > S_heur + 1e-6):', sum(r['S_lb'] > r['S_heur'] + 1e-6 for r in ok), '| rel_err_ub 95%', np.percentile(rel, 95))
print(sh(f"ls {R}/results/_backup/2026-10-20 | grep static"))
```

Output:

```
{'written': '2026-10-20 18:56:22', 'n_H': 12003, 'n_O': 4140, 'n_O_registered': 6000}
{"n": 500, "n_err": 0, "n_certified": 500, "rel_err_ub_median": 0.006914233722314057, "rel_err_ub_p95": 0.010296609910040121, "rel_err_ub_max": 0.02352063861786593, "err_ub_max": 0.09098470963927063, "n_heur_below_lb": 35, "by_src": {"H": {"n": 379, "rel_err_ub_max": 0.02352063861786593}, "O": {"n": 121, "rel_err_ub_max": 0.009408713326138376}}}
H|0.2-0.3 {'n': 25, 'median_D2n': 6.378639427013441, 'q25': 5.503906369920392, 'q75': 6.839560659939331}
H|0.3-0.4 {'n': 486, 'median_D2n': 6.505244909011598, 'q25': 5.777509511775944, 'q75': 7.39404595723551}
H|0.4-0.5 {'n': 2464, 'median_D2n': 6.288428953007781, 'q25': 5.488205599129112, 'q75': 6.973841661302108}
H|0.5-0.6 {'n': 4632, 'median_D2n': 5.590587812662261, 'q25': 4.910151104069522, 'q75': 6.239672034223425}
H|0.6-0.7 {'n': 3133, 'median_D2n': 4.823864719573307, 'q25': 4.1277533028219215, 'q75': 5.441764233976805}
H|0.7-0.8 {'n': 1109, 'median_D2n': 4.147222318621939, 'q25': 3.4690663888563615, 'q75': 4.974075858664795}
H|0.8-1.01 {'n': 154, 'median_D2n': 3.504603001587057, 'q25': 2.5768646790221217, 'q75': 4.287050320839409}
T|0-0.2 {'n': 4, 'median_D2n': 6.695594829790144, 'q25': 6.614043440962289, 'q75': 8.366158651414445}
T|0.2-0.3 {'n': 200, 'median_D2n': 7.40071144048486, 'q25': 6.762739514005833, 'q75': 7.971308762186957}
T|0.3-0.4 {'n': 512, 'median_D2n': 6.72019870366984, 'q25': 5.889493648368024, 'q75': 7.443692799063938}
T|0.4-0.5 {'n': 331, 'median_D2n': 5.947636309525806, 'q25': 5.123008629606039, 'q75': 6.744014488903183}
T|0.5-0.6 {'n': 76, 'median_D2n': 4.961284494705178, 'q25': 4.134162021691082, 'q75': 5.63084329901723}
T|0.6-0.7 {'n': 4, 'median_D2n': 4.7131626951998005, 'q25': 4.508588653780282, 'q75': 4.85297266489973}
A|0-0.2 {'n': 47, 'median_D2n': 8.572922197900182, 'q25': 8.01604794101806, 'q75': 9.014888713688684}
A|0.2-0.3 {'n': 502, 'median_D2n': 7.283553878543293, 'q25': 6.659039140717017, 'q75': 8.016710581729201}
A|0.3-0.4 {'n': 354, 'median_D2n': 6.608382155560986, 'q25': 5.809782035311036, 'q75': 7.591764935397661}
A|0.4-0.5 {'n': 41, 'median_D2n': 5.633356497411857, 'q25': 4.880229933885636, 'q75': 6.076350946566754}
A|0.5-0.6 {'n': 1, 'median_D2n': 3.8382887363604246, 'q25': 3.8382887363604246, 'q75': 3.8382887363604246}
C|0-0.2 {'n': 87, 'median_D2n': 8.64789666840621, 'q25': 8.163854323849204, 'q75': 9.296399306698948}
C|0.2-0.3 {'n': 782, 'median_D2n': 7.672491470049747, 'q25': 6.910683720558244, 'q75': 8.308310519987831}
C|0.3-0.4 {'n': 220, 'median_D2n': 6.9446158111397835, 'q25': 6.101141711170594, 'q75': 7.814059046095469}
C|0.4-0.5 {'n': 14, 'median_D2n': 6.254558059298889, 'q25': 5.415721808135067, 'q75': 7.2426849719964395}
X|0-0.2 {'n': 99, 'median_D2n': 8.923061544660964, 'q25': 8.331160899585848, 'q75': 9.503848235608517}
X|0.2-0.3 {'n': 714, 'median_D2n': 7.877417555162422, 'q25': 7.169442430444778, 'q75': 8.529675017134112}
X|0.3-0.4 {'n': 147, 'median_D2n': 7.166634224238272, 'q25': 6.373061812865261, 'q75': 7.731323978310967}
X|0.4-0.5 {'n': 5, 'median_D2n': 6.525692882434912, 'q25': 5.921010709243553, 'q75': 6.699521962317506}
인증 500 / 500 | ≤60 s 인증 500 | 유사도 틈 ≤ 0.002 0 중앙 0.009152300127484444 최대 0.02455260358086823 | sec 중앙 3.4474289417266846 90% 4.848191738128662 최대 7.1466896533966064
휴리스틱 > 하한(S_heur ≥ S_lb): 465 | BnB 가 휴리스틱보다 높은 유사도를 찾음(S_lb > S_heur + 1e-6): 35 | rel_err_ub 95% 0.010296609910040121
185622_r38_static.json
```

Cell `H_XE91IsBAsy`, run requested 19:05:16, output returned 19:05:16 UTC.

```python
g_lb = np.array([(r['S_lb'] - r['S_heur']) / r['S_lb'] for r in ok]); bad = g_lb > 1e-9
Dh = np.array([r['D_heur'] for r in ok]); nx = np.array([r['D_heur'] ** 2 + 2 * r['S_heur'] for r in ok])  # ‖x‖²+‖y‖²
Dlb = np.sqrt(np.maximum(nx - 2 * np.array([r['S_lb'] for r in ok]), 0)); dD = (Dh - Dlb) / Dh
print('BnB 가 더 나은 정렬을 찾은 쌍', int(bad.sum()), '| 그 쌍들의 유사도 차 중앙', np.median(g_lb[bad]), '최대', g_lb.max(), '| D 상대 차(보인 오차) 중앙', np.median(dD[bad]), '최대', dD.max(), '| > 0.002 인 쌍', int((g_lb > 0.002).sum()))
print('출처별 더 나은 정렬', {s: int(sum(1 for r, b in zip(ok, bad) if b and r['src'] == s)) for s in ('H', 'O')}, '| 길이', np.median([max(r['n'], r['m']) for r, b in zip(ok, bad) if b]), '대 전체', np.median([max(r['n'], r['m']) for r in ok]))
slack = np.array([(r['S_ub'] - r['S_lb']) / r['S_ub'] for r in ok]); print('인증 틈 (S_ub − S_lb)/S_ub 중앙', np.median(slack), '최대', slack.max())
```

Output:

```
BnB 가 더 나은 정렬을 찾은 쌍 35 | 그 쌍들의 유사도 차 중앙 0.00014107332171462843 최대 0.007981857891529246 | D 상대 차(보인 오차) 중앙 8.532206978362482e-05 최대 0.0035579760433795374 | > 0.002 인 쌍 4
출처별 더 나은 정렬 {'H': 18, 'O': 17} | 길이 111.0 대 전체 105.0
인증 틈 (S_ub − S_lb)/S_ub 중앙 0.009146143777022067 최대 0.018912056296239977
```

Cell `H_XE91IsBAsy`, run requested 19:05:43, output returned 19:05:44 UTC.

```python
sl = np.array([(r['S_lb'] - r['S_heur']) / (0.75 * max(r['n'], r['m'])) for r in ok]); tolu = np.array([(r['S_ub'] - r['S_lb']) / (0.75 * max(r['n'], r['m'])) for r in ok])
print('S_ub − S_lb (sim_long) 중앙', np.median(tolu), '최대', tolu.max(), '| S_lb − S_heur (sim_long): ≤ 0.002 인 쌍', int((sl <= 0.002).sum()), '/', len(sl), '최대', sl.max(), '| > 0 인 쌍', int((sl > 1e-12).sum()))
print('등록 정의 (S_ub − S_heur)/S_ub ≤ 0.002:', int(((np.array([r['S_ub'] for r in ok]) - np.array([r['S_heur'] for r in ok])) / np.array([r['S_ub'] for r in ok]) <= 0.002).sum()), '/ 500')
```

Output:

```
S_ub − S_lb (sim_long) 중앙 0.0049999999999999975 최대 0.005000000000000074 | S_lb − S_heur (sim_long): ≤ 0.002 인 쌍 499 / 500 최대 0.003244454275038441 | > 0 인 쌍 35
등록 정의 (S_ub − S_heur)/S_ub ≤ 0.002: 0 / 500
```

Cell `H_XE91IsBAsy`, run requested 19:06:03, output returned 19:06:04 UTC.

```python
o1 = "    b = out['bnb']; out['prediction'] = dict(cert_ge95=(b['frac_certified_le60s'] or 0) >= 0.95, heur_ge95=(b['frac_sim_gap_le_0p002'] or 0) >= 0.95)\n"
n1 = o1 + ("    # (고침 1, 2026-10-20 19:1x UTC, 사후 서술) 등록 정의는 S_ub = S_lb + 허용 폭(0.005 sim_long) 이라 정의상 0.002 안일 수 없다. 같은 단위(sim_long = S / (0.75 × 긴 쪽 변 수))로 BnB 가 찾은 최선값과 휴리스틱의 차를 함께 낸다\n"
           "    sl = [(r['S_lb'] - r['S_heur']) / (0.75 * max(r['n'], r['m'])) for r in ok]; out['bnb_posthoc'] = dict(frac_heur_within_0p002_simlong_of_S_lb=float(np.mean([g <= 0.002 for g in sl])), max_shortfall_simlong=float(max(sl)),\n"
           "        n_bnb_better=int(sum(g > 1e-12 for g in sl)), note='사후 서술: 등록 정의(S_ub − S_heur)/S_ub ≤ 0.002 는 S_ub 에 허용 폭 0.005 sim_long 이 늘 더해져 설 수 없음(등록 오류). 판정 값은 그대로 보고')\n")
s_ = open(f'{R}/code/r39_sec6.py').read(); assert s_.count(o1) == 1; md_ = hashlib.md5(s_.replace(o1, n1).encode()).hexdigest(); patch('r39_sec6.py', [(o1, n1)], md_)
print(sh(f"/usr/bin/python3 -m py_compile {R}/code/r39_sec6.py && echo 문법 통과"))
P = json.load(open(PJ)); P['round36_sec6_2026-10-20']['fix1'] = dict(stamped_utc=datetime.datetime.now(datetime.UTC).strftime('%F %T'), script_md5=md_,
    finding='정적 축 BnB 결과(18:56 UTC)를 보고 안 것: fs.bnb 의 허용 폭은 rel_tol × 0.75 × max(n, m) (sim_long 단위)이고 S_ub = S_lb + 허용 폭이다. 그래서 17:34 에 등록한 정의 (S_ub − S_heur)/S_ub ≤ 0.002 는 정의상 0/500 이다(등록 오류).',
    change='판정 값은 등록 정의 그대로 내고, 사후 서술로 (S_lb − S_heur)/(0.75 max(n, m)) ≤ 0.002 의 비율을 더함. 결과를 본 뒤의 고침이라 투명성 표에 적음.')
json.dump(P, open(PJ + '.tmp', 'w'), ensure_ascii=False, indent=1); os.replace(PJ + '.tmp', PJ); print('찍음', md_)
```

Output:

```
r39_sec6.py 1b11b675bfd25b2e34dafebc217b5618 일치
문법 통과
찍음 1b11b675bfd25b2e34dafebc217b5618
```

## 15. Negative type of the edit distance on two-bond chains (18 October 2026; checks of 30 and 31 October)

S3 Text 2: for two-bond chains, d depends only on the angle θ between the two relative rotations,
d = h·min(√6, 4 sin(θ/4)); a positive coefficient of this class function against a character of SO(3) rules out
negative type, and so does a positive eigenvalue of J D J. These computations were not specified in advance. They were
run locally with Node.js on 18 October 2026 from `code/negtype/negtype.js`, `so3fourier.js` and `len2check.js`, in a
scratch directory (written <scratch>) from which the three files were copied to the project at 07:56:38 UTC; their md5
at 07:57:40, `735bb530…`, `acd21303…` and `76e6c6c8…`, are those of the deposited files (`docs/PROVENANCE.md`). The
scripts read no input files, and their random rotations come from the fixed-seed generator of `negtype.js`, so every run
prints the same values. A first search of `negtype.js` (07:54:30, left out) drew structures of one bond only and gave a
distance matrix of zeros; the line that draws the lengths was then changed to at least two bonds, the deposited
version, by the `sed` in the first command below, which reran the search (8 random structures of 2–3 bonds: largest
eigenvalue of J D J zero to rounding). `so3fourier.js` was written at 07:55:21; the `sed` in the second command set its
comparison line (the Frobenius distance between rotation matrices) to the deposited text before it ran. The third
command wrote `len2check.js`, whose text is the deposited file, and ran it. S3 Text 2 takes from these outputs
f̂_5 = +1.15 × 10⁻² (h = 1), the positive coefficients up to l = 14 at l = 5, 7, 10, 12 and 14, the coefficients of the
version without the √6 cap (positive at every even l) and of the Frobenius distance (negative throughout), the agreement
of the closed form with the general computation (1.28 × 10⁻¹⁵) and the largest eigenvalue of J D J (zero to rounding at
N = 60, positive at N = 120 and 200). On 30 October the coefficients up to l = 14 were recomputed independently with
SciPy (fourth command; Haar density (1 − cos θ)/π, χ_l(θ) = sin((2l + 1)θ/2)/sin(θ/2), adaptive quadrature), with the
same values to the digits printed. On 31 October `negtype_eta_len2.js` (fifth command; the file written by the command
is the deposited one, md5 `8fd84e58…`) gave η₊ for random two-bond chains, quoted in S3 Text 2. The local
working-copy path is written <atlas>; the line "Shell cwd was reset …" that the tool appended to some outputs is left
out. In the outputs, 최대 상대 위반 is the largest relative violation, 구조 길이 the lengths of the structures, 양수 계수
the positive coefficients, 상한 포함 / 상한 없음 with and without the √6 cap, 회전 행렬 현 거리 the Frobenius distance
between rotation matrices, and 공식과 실제 D 의 최대 차 the largest difference between the closed form and the general
computation.

Command run 2026-10-18 07:54:38, output returned 07:54:39 UTC (session log).

```bash
cd "<scratch>" && sed -i 's/const n = 1 + Math.floor(rnd()\*maxLen);/const n = 2 + Math.floor(rnd()*(maxLen-1));/' negtype.js && grep -n "const n = 2" negtype.js && timeout 500 node negtype.js search 3 8 40
```

Output:

```
70:function randStruct(maxLen) { const n = 2 + Math.floor(rnd()*(maxLen-1)); const out=[]; for (let i=0;i<n;i++) out.push(scale(randRot(), 0.5)); return out; }
최대 상대 위반 λmax(JDJ)/max D = 1.5937e-15  λmax = 2.3638e-15
구조 길이: 2,2,2,3,2,3,3,2
```

Command run 2026-10-18 07:55:27, output returned 07:55:29 UTC (session log).

```bash
cd "<scratch>" && sed -i 's#const fChord = th => .*#const fChord = th => Math.sqrt(4 * (1 - Math.cos(th)));            // 비교: 회전 행렬 현 거리 ‖R−R′‖_F (유클리드라 음의 유형)#' so3fourier.js && node so3fourier.js
```

Output:

```
D (상한 포함): 양수 계수 l = 5(1.15e-2) 7(5.53e-3) 10(1.92e-3) 12(2.81e-3) 14(9.69e-4) 17(1.21e-3) 19(1.08e-3) 21(7.62e-5) 22(2.95e-4) 24(7.55e-4)
   f̂_1..8 = -4.97e-1 -1.68e-2 -7.80e-3 -2.32e-2 1.15e-2 -1.40e-2 5.53e-3 -3.86e-3
상한 없음: 양수 계수 l = 2(3.71e-2) 4(1.59e-2) 6(8.51e-3) 8(5.26e-3) 10(3.56e-3) 12(2.57e-3) 14(1.94e-3) 16(1.52e-3) 18(1.22e-3) 20(9.99e-4) 22(8.35e-4) 24(7.08e-4)
   f̂_1..8 = -5.56e-1 3.71e-2 -5.49e-2 1.59e-2 -1.92e-2 8.51e-3 -9.65e-3 5.26e-3
회전 행렬 현 거리: 양수 계수 l = 없음
   f̂_1..8 = -4.80e-1 -6.86e-2 -2.29e-2 -1.04e-2 -5.60e-3 -3.36e-3 -2.17e-3 -1.49e-3
```

Command run 2026-10-18 07:55:51, output returned 07:55:53 UTC (session log).

```bash
cd "<scratch>" && cat > len2check.js <<'EOF'
const src = require('fs').readFileSync('negtype.js','utf8').split("const mode =")[0];
eval(src);
function angleBetween(Q1,Q2){ const P = matMul(Q1, transpose(Q2)); const c=(P[0]+P[4]+P[8]-1)/2; return Math.acos(Math.max(-1,Math.min(1,c))); }
const fD = th => 0.5*Math.min(Math.sqrt(6), 4*Math.sin(th/4));
let maxErr=0;
for (let t=0;t<200;t++){ const Q=randRot(), Qp=randRot(); const I=[1,0,0,0,1,0,0,0,1];
  const x=[scale(I,0.5),scale(Q,0.5)], y=[scale(I,0.5),scale(Qp,0.5)];
  const d=D(x,y), f=fD(angleBetween(Q,Qp)); maxErr=Math.max(maxErr,Math.abs(d-f)); }
console.log('공식과 실제 D 의 최대 차:', maxErr.toExponential(2));
// 무작위 길이 2 구조 N 개에서 J D J 의 최대 고유값
for (const N of [60,120,200]) {
  const Qs=[]; for(let i=0;i<N;i++) Qs.push(randRot());
  const Dm=Qs.map(a=>Qs.map(b=>fD(angleBetween(a,b))));
  const rm=Dm.map(r=>r.reduce((p,q)=>p+q,0)/N), all=rm.reduce((p,q)=>p+q,0)/N;
  const C=Dm.map((r,i)=>r.map((v,j)=>v-rm[i]-rm[j]+all));
  const ev=jacobiEig(C).sort((a,b)=>b-a);
  console.log('N='+N+': J D J 고유값 상위 3 =', ev.slice(0,3).map(v=>v.toExponential(3)).join(' '), ' (음의 유형이면 모두 ≤ 0)');
}
EOF
timeout 600 node len2check.js
```

Output:

```
공식과 실제 D 의 최대 차: 1.28e-15
N=60: J D J 고유값 상위 3 = 3.438e-14 -5.522e-2 -6.119e-2  (음의 유형이면 모두 ≤ 0)
N=120: J D J 고유값 상위 3 = 3.980e-2 3.224e-2 1.960e-2  (음의 유형이면 모두 ≤ 0)
N=200: J D J 고유값 상위 3 = 1.276e-1 1.165e-1 1.089e-1  (음의 유형이면 모두 ≤ 0)
```

Command run 2026-10-30 14:51:17, output returned 14:51:48 UTC (session log).

```bash
cd "<scratch>" && PYTHONUTF8=1 uv run --offline --no-project --python 3.12 --with numpy --with scipy python -X utf8 -c "
import numpy as np
from scipy.integrate import quad
# SO(3) class functions: Haar density on angle theta in [0, pi] is (1 - cos theta)/pi; character chi_l = sin((2l+1)t/2)/sin(t/2)
f_cap = lambda t: min(np.sqrt(6), 4*np.sin(t/4))
f_nocap = lambda t: 4*np.sin(t/4)
f_frob = lambda t: np.sqrt(2*(3 - (1 + 2*np.cos(t))))   # ||R - R'||_F = sqrt(6 - 2 tr) ; tr = 1 + 2cos t
chi = lambda l, t: np.sin((2*l+1)*t/2)/np.sin(t/2) if t > 1e-12 else 2*l+1
haar = lambda t: (1 - np.cos(t))/np.pi
print('check norm', quad(haar, 0, np.pi)[0])
for name, f in [('cap', f_cap), ('nocap', f_nocap), ('frob', f_frob)]:
    co = [quad(lambda t: f(t)*chi(l, t)*haar(t), 0, np.pi, limit=400, points=[4*np.arcsin(np.sqrt(6)/4)] if name=='cap' else None)[0] for l in range(1, 15)]
    print(name, ' '.join(f'{l}:{c:+.2e}' for l, c in zip(range(1, 15), co)))
"
```

Output:

```
Installed 2 packages in 3.03s
check norm 1.0
cap 1:-4.97e-01 2:-1.68e-02 3:-7.80e-03 4:-2.32e-02 5:+1.15e-02 6:-1.40e-02 7:+5.53e-03 8:-3.86e-03 9:-1.19e-03 10:+1.92e-03 11:-3.72e-03 12:+2.81e-03 13:-2.57e-03 14:+9.69e-04
nocap 1:-5.56e-01 2:+3.71e-02 3:-5.49e-02 4:+1.59e-02 5:-1.92e-02 6:+8.51e-03 7:-9.65e-03 8:+5.26e-03 9:-5.78e-03 10:+3.56e-03 11:-3.85e-03 12:+2.57e-03 13:-2.74e-03 14:+1.94e-03
frob 1:-4.80e-01 2:-6.86e-02 3:-2.29e-02 4:-1.04e-02 5:-5.60e-03 6:-3.36e-03 7:-2.17e-03 8:-1.49e-03 9:-1.06e-03 10:-7.85e-04 11:-5.96e-04 12:-4.64e-04 13:-3.68e-04 14:-2.97e-04
```

Command run 2026-10-31 08:54:59, output returned 08:55:09 UTC (session log).

```bash
S="<scratch>"; mkdir -p "$S/eta"; cat > "$S/eta/eta_len2.js" <<'EOF'
// eta_+ for random two-bond chains (S3 Text 2), using the helpers of framedchain-paper/code/negtype/negtype.js
const path = process.argv[2];
const src = require('fs').readFileSync(path + '/negtype.js', 'utf8').split("const mode =")[0];
eval(src);
function angleBetween(Q1,Q2){ const P = matMul(Q1, transpose(Q2)); const c=(P[0]+P[4]+P[8]-1)/2; return Math.acos(Math.max(-1,Math.min(1,c))); }
const fD = th => 0.5*Math.min(Math.sqrt(6), 4*Math.sin(th/4));
for (const [N, reps] of [[60,3],[100,5],[200,3],[400,2]]) {
  const out = [];
  for (let r = 0; r < reps; r++) {
    const Qs=[]; for(let i=0;i<N;i++) Qs.push(randRot());
    const Dm=Qs.map(a=>Qs.map(b=>fD(angleBetween(a,b))));
    const rm=Dm.map(r=>r.reduce((p,q)=>p+q,0)/N), all=rm.reduce((p,q)=>p+q,0)/N;
    const C=Dm.map((r,i)=>r.map((v,j)=>v-rm[i]-rm[j]+all));
    const ev=jacobiEig(C); const pos=ev.filter(v=>v>1e-9).reduce((p,q)=>p+q,0), tot=ev.reduce((p,q)=>p+Math.abs(q),0);
    out.push((pos/tot).toExponential(2) + ' (max eig ' + Math.max(...ev).toExponential(2) + ')');
  }
  console.log('N=' + N + ': eta_+ =', out.join(', '));
}
EOF
cd "$S/eta" && timeout 900 node eta_len2.js "<atlas>/framedchain-paper/code/negtype"
```

Output:

```
N=60: eta_+ = 0.00e+0 (max eig 2.14e-14), 0.00e+0 (max eig 2.49e-14), 0.00e+0 (max eig 2.19e-14)
N=100: eta_+ = 2.87e-4 (max eig 2.52e-2), 3.51e-4 (max eig 2.49e-2), 1.63e-4 (max eig 1.19e-2), 1.45e-4 (max eig 1.39e-2), 2.79e-4 (max eig 2.46e-2)
N=200: eta_+ = 6.78e-3 (max eig 1.35e-1), 6.69e-3 (max eig 1.43e-1), 6.76e-3 (max eig 1.40e-1)
N=400: eta_+ = 2.43e-2 (max eig 3.34e-1), 2.43e-2 (max eig 3.25e-1)
```

## 16. Stretch and tighter bound for Fig 1 (3 November 2026)

Not specified in advance. (a) While choosing the fragments of Fig 1, a cell listed every stretch of the whole-chain
alignment of 1AKI and 1LZ1 (255 random starts) with exactly one unmatched bond and four matched bonds on each side; there
is one, at pair index 45, which `make_fig1.py` takes (it keeps the four matched bonds before and three after the
unmatched bond). Lines of the cell as run (its setup, and the alignment of each stretch that printed the last line, omitted):

```python
al = fc.align(xa, xb, n_random=255)
pairs = sorted(zip(*[t.tolist() for t in torch.nonzero(al.match, as_tuple=True)]))
print('whole chains: matched', len(pairs), '| unmatched 1AKI', len(xa) - len(pairs), '| unmatched 1LZ1', len(xb) - len(pairs))
cands = []
for k in range(4, len(pairs) - 4):
    (i0, j0), (i1, j1) = pairs[k - 1], pairs[k]
    run_ok = all(pairs[q + 1][0] - pairs[q][0] == 1 and pairs[q + 1][1] - pairs[q][1] == 1 for q in list(range(k - 4, k - 1)) + list(range(k, k + 3)))
    if run_ok and {(i1 - i0), (j1 - j0)} == {1, 2}:
        cands.append(k)
print('single-bond gaps with 4 matched bonds on each side at pair index:', cands)
```

```
bonds 128 129 | breaks 0 0
whole chains: matched 126 | unmatched 1AKI 2 | unmatched 1LZ1 3
single-bond gaps with 4 matched bonds on each side at pair index: [45]
1AKI 42–49 ATNRNTDG | 1LZ1 42–50 ATNYNAGDR | cols [(0, 0), (1, 1), (2, 2), (3, 3), (-1, 4), (4, 5), (5, 6), (6, 7)] | d² 1.3583 = 0.6190 + 0.7393 | certified True
```

(b) A bound on how far the drawn alignment can be from the optimum, tighter than the script's default tolerance (S
within 0.03, d² within 0.06), run in the same kernel (`fa`, `fb`: the two fragments; `a2`: the drawn alignment). At
tolerance 10⁻⁶ the search stopped at its cell limit (not certified), with upper bound S_ub 4.897544 against the drawn S
4.897116, so the optimal d² lies between 1.357441 and the drawn 1.358296.

```python
c6 = fc.bnb(fa, fb, tol=1e-6)
S_a2 = float(a2.score); nx, ny = float((fa ** 2).sum()), float((fb ** 2).sum())
print('S(drawn)', round(S_a2, 9), '| bnb tol 1e-6: S_lb', round(c6['S_lb'], 9), 'S_ub', round(c6['S_ub'], 9), 'certified', c6['certified'], '| drawn >= S_ub - tol:', S_a2 >= c6['S_ub'] - c6['tol'])
print('d2(drawn)', round(nx + ny - 2 * S_a2, 6), '| d2 lower bound', round(nx + ny - 2 * c6['S_ub'], 6))
```

```
S(drawn) 4.897115906 | bnb tol 1e-6: S_lb 4.897115906 S_ub 4.897543541 certified False | drawn >= S_ub - tol: False
d2(drawn) 1.358296 | d2 lower bound 1.357441
```
