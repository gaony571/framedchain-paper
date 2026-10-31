import matplotlib
"""Figures 2–7 and S1 Fig for the manuscript, drawn only from stored result files (no recomputation).
Output: manuscript_figs/<name>.tif (PLOS: RGB, LZW, 600 dpi) and preview PNG (150 dpi)."""
import os, json, glob, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
R = '/content/drive/MyDrive/atlas_framed_curves'; OUT = R + '/manuscript_figs'; os.makedirs(OUT, exist_ok=True)
J = lambda f: json.load(open(f'{R}/results/{f}'))
plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 8, 'axes.titlesize': 8, 'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
                     'legend.fontsize': 7, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'lines.linewidth': 0.9,
                     'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42, 'savefig.dpi': 600,
                     'mathtext.fontset': 'custom', 'mathtext.rm': 'Liberation Sans', 'mathtext.it': 'Liberation Sans:italic', 'mathtext.bf': 'Liberation Sans:bold'})
OI = dict(black='#000000', orange='#E69F00', sky='#56B4E9', green='#009E73', yellow='#F0E442', blue='#0072B2', red='#D55E00', purple='#CC79A7', grey='#999999')
STEP = {'348': '+28 K', '379': '+59 K', '413': '+93 K', '450': '+130 K'}
SCOL = {'348': OI['sky'], '379': OI['blue'], '413': OI['orange'], '450': OI['red']}

def panel(ax, s, x=-0.14, y=1.04):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom', ha='left')

def save(fig, name):
    from PIL import Image
    fig.savefig(f'{OUT}/_tmp_{name}.png', dpi=600, facecolor='white'); fig.savefig(f'{OUT}/{name}_preview.png', dpi=150, facecolor='white'); plt.close(fig)
    im = Image.open(f'{OUT}/_tmp_{name}.png').convert('RGB'); im.save(f'{OUT}/{name}.tif', compression='tiff_lzw', dpi=(600, 600)); os.remove(f'{OUT}/_tmp_{name}.png')
    im = Image.open(f'{OUT}/{name}.tif'); return dict(name=name, size_px=im.size, size_in=tuple(round(v / 600, 2) for v in im.size), mode=im.mode, dpi=im.info.get('dpi'))

def ci(d, k='est', c='ci95'): return d[k], d[c][0], d[c][1]

# ---------------- Fig 2 ----------------
def fig2():
    t = pd.read_csv(f'{R}/results/pu_tm_factorial_B_triangle.tsv', sep='\t'); t = t[t.tol == 1e-4]
    def rate(d0, norm, form):
        s = t[(t.d0 == d0) & (t.norm == norm) & (t.form == form)]
        if not len(s): return np.nan
        for lev in ('max', 'hi', 'tmalign'):
            if (s.level == lev).any(): return float(s[s.level == lev].rate_all.iloc[0])
    rows = [(r'Unnormalized, fixed $d_0$', 'fix', 'none'), (r'Unnormalized, $d_0(\bar L)$', 'var_avg', 'none'),
            (r'Shorter chain, fixed $d_0$', 'fix', 'short'), (r'Shorter chain, $d_0(L)$', 'var_short', 'short'),
            (r'Mean length, fixed $d_0$', 'fix', 'avg'), (r'Mean length, $d_0(\bar L)$', 'var_avg', 'avg'),
            (r'Longer chain, fixed $d_0$', 'fix', 'long'), (r'Longer chain, $d_0(L)$', 'var_long', 'long')]
    M = np.array([[rate(d0, nm, 'sqrt'), rate(d0, nm, 'lin')] for _, d0, nm in rows]) * 100
    tm = [rate('TMalign', nm, 'lin') * 100 for nm in ('short', 'avg', 'long')]
    fig = plt.figure(figsize=(7.5, 3.4)); gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1], wspace=0.62, left=0.235, right=0.99, bottom=0.17, top=0.9)
    ax = fig.add_subplot(gs[0]); panel(ax, 'A', x=-0.85)
    labs = [r[0] for r in rows] + ['TM-align, shorter chain', 'TM-align, mean length', 'TM-align, longer chain']
    Z = np.vstack([M, np.c_[np.full(3, np.nan), tm]])
    from matplotlib.colors import LogNorm
    show = np.where(Z > 0, Z, np.nan)
    im = ax.imshow(show, cmap='Reds', norm=LogNorm(vmin=1e-3, vmax=10), aspect='auto')
    for i in range(Z.shape[0]):
        for j in range(2):
            v = Z[i, j]
            if np.isnan(v): ax.text(j, i, '—', ha='center', va='center', color=OI['grey'], fontsize=7); continue
            txt = '0' if v == 0 else (f'{v:.2f}' if v >= 0.1 else f'{v:.3f}')
            ax.text(j, i, txt, ha='center', va='center', fontsize=7, color='white' if v > 1 else 'black')
    ax.set_xticks([0, 1]); ax.set_xticklabels(['square root\n' + r'$\sqrt{L - S}$', 'linear\n' + r'$L - S$'])
    ax.set_yticks(range(len(labs))); ax.set_yticklabels(labs)
    ax.axhline(7.5, color='black', lw=0.6); ax.tick_params(length=0)
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.set_title('% of 161,700 triples violating\nthe triangle inequality', fontsize=7.5)
    cb = fig.colorbar(im, ax=ax, fraction=0.06, pad=0.04); cb.set_label('% violated (log scale)', fontsize=6.5); cb.ax.tick_params(labelsize=6)
    # B: chaining
    c1, c2 = J('chaining.json'), J('chaining_v2.json')
    meas = [('1 − TM,\nshorter chain', c2['tms']), ('1 − TM,\nlonger chain', c1['tml']), (r'$D_{TM}/\sqrt{\bar L}$' + ',\n' + r'fixed $d_0$', c2['dfn']), (r'$1 - S/\bar L$' + '\n(linear)', c2['dn_lin'])]
    ax = fig.add_subplot(gs[1]); panel(ax, 'B', x=-0.22)
    x = np.arange(len(meas)); w = 0.38
    for k, (key, col, lab) in enumerate([('triangle_violation', OI['grey'], 'triples violating the triangle inequality'), ('average_link_merge', OI['red'], 'chains whose folds are merged (average linkage)')]):
        v = np.array([m[1][key] for m in meas], float) * 100
        b = ax.bar(x + (k - 0.5) * w, v[:, 0], w, color=col, label=lab, edgecolor='none')
        ax.errorbar(x + (k - 0.5) * w, v[:, 0], yerr=[v[:, 0] - v[:, 1], v[:, 2] - v[:, 0]], fmt='none', ecolor='black', elinewidth=0.6, capsize=1.5)
        for xi, vi, hi_ in zip(x + (k - 0.5) * w, v[:, 0], v[:, 2]): ax.text(xi, hi_ + 2.0, f'{vi:.0f}', ha='center', va='bottom', fontsize=6.5)
    ax.set_xticks(x); ax.set_xticklabels([m[0] for m in meas], fontsize=6.2); ax.set_ylabel('% (150 two-domain chains)'); ax.set_ylim(0, 142); ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.legend(frameon=False, loc='upper center', fontsize=6.3, bbox_to_anchor=(0.5, 1.02))
    return save(fig, 'Fig2'), dict(M=M.round(4).tolist(), tm=np.round(tm, 4).tolist(), chain={m[0]: (m[1]['triangle_violation'], m[1]['average_link_merge']) for m in meas})

# ---------------- Fig 3 ----------------
def fig3():
    E0, E1, E2, E3 = J('r39_E0.json'), J('r39_E1.json'), J('r39_E2.json'), J('r39_E3.json')
    fig = plt.figure(figsize=(7.5, 2.7)); gs = fig.add_gridspec(1, 3, wspace=0.48, left=0.085, right=0.99, bottom=0.19, top=0.9)
    ax = fig.add_subplot(gs[0]); panel(ax, 'A', x=-0.3)
    lab = ['+28 K', '+59 K', '+93 K', '+130 K']; T = ['348', '379', '413', '450']; vals = []
    for i, (l, tt) in enumerate(zip(lab, T)):
        e, lo, hi = ci(E0['cells'][f'온도 {l}|줄기|half']['d2n']); vals.append((e, lo, hi))
        ax.errorbar(i, e, yerr=[[e - lo], [hi - e]], fmt='o', color=SCOL[tt], ms=4, capsize=2, elinewidth=0.8)
    b = E2['R0_band_half']['band']; ax.axhspan(b[0], b[1], color=OI['grey'], alpha=0.3, lw=0); ax.axhline(0, color='black', lw=0.4)
    ax.text(3.35, b[1] + 0.001, 'replicate-only\nband', fontsize=6, ha='right', va='bottom', color='#555555')
    ax.set_xticks(range(4)); ax.set_xticklabels(lab); ax.set_xlabel('heating from 320 K'); ax.set_ylabel(r'$\Delta^2_D/\sqrt{\bar n}$' + ' (superfamily median)'); ax.set_xlim(-0.5, 3.5)
    # B: size dependence
    ax = fig.add_subplot(gs[1]); panel(ax, 'B', x=-0.3); fits = {}
    recs = {tt: [] for tt in T}
    for f in glob.glob(f'{R}/derived/r39/D/*.json'):
        d = json.load(open(f))
        for tt in T:
            c = d['cells'].get(f'MD{tt}-MD320|half')
            if c and c.get('D') and c['D'].get('D2') is not None: recs[tt].append((d['ne'], c['D']['D2']))
    for tt in T:
        a = np.array(recs[tt], float); n, y = a[:, 0], a[:, 1]; pos = y > 0
        ax.scatter(n[pos], y[pos], s=3, color=SCOL[tt], alpha=0.5, lw=0, label=f'{STEP[tt]} (n = {len(n)})')
        al = E1['steps'][tt]['irls']['alpha']; ga = E1['steps'][tt]['irls']['gamma']
        w = n ** (-ga); x = n ** al; cfit = (w * y * x).sum() / (w * x * x).sum(); fits[tt] = (al, ga, cfit, len(n), int((~pos).sum()))
        g = np.linspace(n.min(), n.max(), 50); ax.plot(g, cfit * g ** al, color=SCOL[tt], lw=1)
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xticks([50, 70, 100, 150]); ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter()); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.set_xlabel('number of bonds n'); ax.set_ylabel(r'$\Delta^2_D$' + ' per domain (positive values)')
    ax.legend(frameon=False, fontsize=5.5, loc='lower right', markerscale=2, handletextpad=0.2)
    # C: power
    ax = fig.add_subplot(gs[2]); panel(ax, 'C', x=-0.3); pw = {}
    for tt in T:
        rr = sorted([(int(r['n']), r['p_detect']) for r in E3['rows'] if r['window'] == 'half' and r['method'] == '나' and str(r['T']) == tt and str(r['n']).isdigit()])
        pw[tt] = rr; ax.plot([a for a, _ in rr], [b for _, b in rr], 'o-', color=SCOL[tt], ms=3, label=STEP[tt])
    ax.set_xscale('log'); ax.set_xlabel('number of domains'); ax.set_ylabel('probability of detection'); ax.set_ylim(0, 1.03); ax.axhline(0.8, color=OI['grey'], lw=0.5, ls='--')
    ax.legend(frameon=False, fontsize=6, loc='lower right')
    return save(fig, 'Fig3'), dict(A=vals, R0=b, fits=fits, power=pw)

# ---------------- Fig 4 ----------------
CELL_EN = {'온도 +28 K (뒤 절반)': 'Heating +28 K', '온도 +59 K (뒤 절반)': 'Heating +59 K', '온도 +93 K (뒤 절반)': 'Heating +93 K', '온도 +130 K (뒤 절반)': 'Heating +130 K',
           'R3′ (50–100 ns, 주)': 'Force field, temperature and protocol', '7c aSAMt 320 K': 'aSAMt error', '7c AlphaFlow': 'AlphaFlow error',
           '상동체 TM 0.2-0.3': 'Homologs, TM 0.2–0.3', '상동체 TM 0.3-0.4': 'Homologs, TM 0.3–0.4', '상동체 TM 0.4-0.5': 'Homologs, TM 0.4–0.5', '상동체 TM 0.5-0.6': 'Homologs, TM 0.5–0.6',
           '상동체 TM 0.6-0.7': 'Homologs, TM 0.6–0.7', '상동체 TM 0.7-0.8': 'Homologs, TM 0.7–0.8', '상동체 TM 0.8-1.01': 'Homologs, TM ≥ 0.8', '자연 삽입·결실 (R10a′ 주)': 'Natural insertion–deletion'}
def fig4():
    pts = []
    for f in sorted(glob.glob(f'{R}/derived/r38/b6/*.json')):
        b = json.load(open(f)); D = b.get('D')
        if not D or any(D.get(k) is None for k in ('C', 'B_P', 'B_Q')): continue
        d2 = 2 * D['C'] - D['B_P'] - D['B_Q']; pts.append((b['TM'], d2 / np.sqrt(b['nbar']), b.get('W2'), b['Ha'] == b['Hb']))
    P = pd.DataFrame(pts, columns=['TM', 'd2n', 'W2', 'same'])
    fig = plt.figure(figsize=(7.5, 3.5)); gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.05], wspace=1.05, left=0.075, right=0.975, bottom=0.14, top=0.92)
    ax = fig.add_subplot(gs[0]); panel(ax, 'A', x=-0.22)
    for same, col, lab in [(False, OI['grey'], 'different superfamilies'), (True, OI['blue'], 'same superfamily')]:
        s = P[P.same == same]; ax.scatter(s.TM, s.d2n, s=6, color=col, lw=0, alpha=0.8, label=f'{lab} ({len(s)})')
    ax.set_xlabel('TM-score (first 320 K frames, longer chain)'); ax.set_ylabel(r'$\Delta^2_D/\sqrt{\bar n}$'); ax.set_xlim(0.18, 1.0)
    ax2 = ax.twinx(); w = P.dropna(subset=['W2']); ax2.scatter(w.TM, w.W2, s=6, marker='x', color=OI['red'], lw=0.5, alpha=0.7, label=r'$W_2$' + f' ({len(w)} pairs, TM ≥ 0.5)')
    ax2.set_ylabel(r'$W_2$' + ' (units of ' + r'$d$' + ')', color=OI['red']); ax2.tick_params(axis='y', colors=OI['red']); ax2.spines['right'].set_visible(True); ax2.spines['right'].set_color(OI['red'])
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels(); ax.legend(h1 + h2, l1 + l2, frameon=False, fontsize=6, loc='lower left')
    # B: kappa per cell
    s1 = J('r43_sens.json')['s1']['cells']; order = [k for k in CELL_EN if k in s1]
    ax = fig.add_subplot(gs[1]); panel(ax, 'B', x=-0.98); y = np.arange(len(order))[::-1]; kap = {}
    for yi, k in zip(y, order):
        e, lo, hi = ci(s1[k]['ratio']); kap[CELL_EN[k]] = (e, lo, hi)
        col = OI['red'] if 'AlphaFlow' in CELL_EN[k] else OI['blue'] if 'Homolog' in CELL_EN[k] or 'insertion' in CELL_EN[k] else OI['black']
        ax.errorbar(e, yi, xerr=[[e - lo], [hi - e]], fmt='o', color=col, ms=3, capsize=1.5, elinewidth=0.7)
    ax.set_yticks(y); ax.set_yticklabels([CELL_EN[k] for k in order], fontsize=6.5); ax.set_xscale('log'); ax.set_xlim(0.05, 3)
    ax.axvspan(0.5, 2, color=OI['yellow'], alpha=0.25, lw=0); ax.axvline(1, color='black', lw=0.6, ls='--')
    ax.text(0.07, y[0] + 0.9, 'near', fontsize=6.5, color='#555555'); ax.text(0.93, y[0] + 0.9, 'crossing', fontsize=6.5, color='#555555', ha='right'); ax.text(2.4, y[0] + 0.9, 'far', fontsize=6.5, color='#555555', ha='center')
    ax.set_xlabel(r'$\kappa = \sqrt{(C/\bar w)^2 - 1}$'); ax.set_ylim(-0.8, len(order) + 0.3)
    sp = json.load(open(f'{R}/derived/r44/split.json'))['cells']['7c AlphaFlow 대 MD 320|full']; sh = sp['share_mean']['est']
    ya = y[order.index('7c AlphaFlow')]; yfr = (ya + 0.8) / (len(order) + 1.1); ins = ax.inset_axes([0.8, yfr - 0.018, 0.18, 0.036], transform=ax.transAxes)
    ins.barh([0], [sh], color=OI['red'], height=1); ins.barh([0], [1 - sh], left=[sh], color=OI['orange'], height=1); ins.set_xlim(0, 1); ins.axis('off')
    ins.text(0, 1.2, r'$\kappa^2$' + f' = shift {sh:.0%}\n+ spread {1 - sh:.0%}', ha='left', va='bottom', fontsize=5.5, linespacing=1.0)
    return save(fig, 'Fig4'), dict(n_pairs=len(P), n_w2=int(P.W2.notna().sum()), n_same=int(P.same.sum()), kappa=kap, alphaflow_shift_share=sh)

# ---------------- Fig 5 ----------------
def fig5():
    P2 = J('r38_p2_p2_67_D.json'); E5 = J('r39_E5.json')['cells']
    rows = P2['rows']; rows = sorted(rows, key=lambda r: r['model_dmu'])
    fig = plt.figure(figsize=(7.5, 3.2)); gs = fig.add_gridspec(1, 2, width_ratios=[2.5, 1], wspace=0.32, left=0.09, right=0.99, bottom=0.15, top=0.92)
    ax = fig.add_subplot(gs[0]); panel(ax, 'A', x=-0.08)
    for i, r in enumerate(rows):
        for tt, v in r['resp'].items():
            if v is not None: ax.scatter(i, v, s=7, color=SCOL[tt], lw=0, zorder=2)
        ax.scatter(i, r['model_dmu'], s=12, marker='D', color='black', lw=0, zorder=3)
    ax.set_yscale('symlog', linthresh=1e-3, linscale=0.3); ax.set_xlim(-1, len(rows))
    ax.set_xticks([]); ax.set_xlabel(f'{len(rows)} held-out domains, ordered by model error'); ax.set_ylabel(r'$\Delta\mu^2_{2,D}$' + ' (mean shift per bond)')
    ax.axhline(0, color='black', lw=0.4)
    hs = [Line2D([], [], marker='D', ls='', color='black', ms=4, label='aSAMt error at 320 K')] + [Line2D([], [], marker='o', ls='', color=SCOL[tt], ms=4, label=f'MD heating {STEP[tt]}') for tt in SCOL]
    ax.legend(handles=hs, frameon=False, fontsize=6.3, loc='upper center', ncol=5, bbox_to_anchor=(0.62, -0.075), columnspacing=1.0, handletextpad=0.2)
    ax = fig.add_subplot(gs[1]); panel(ax, 'B', x=-0.32); out = {}
    for j, (key, lab) in enumerate([('7c aSAMt320 대 MD 320 (P2 67)|T=top', 'aSAMt'), ('7c AlphaFlow 대 MD 320 (P2 67)|T=top', 'AlphaFlow')]):
        c = E5[key]
        for dx, k, col, mk in [(-0.22, 'ceil_cos', OI['grey'], 's'), (0, 'cos', OI['red'] if lab == 'AlphaFlow' else 'black', 'o'), (0.22, 'null_cos', OI['grey'], 'v')]:
            e, lo, hi = ci(c[k]); out[f'{lab}.{k}'] = (e, lo, hi)
            ax.errorbar(j + dx, e, yerr=[[e - lo], [hi - e]], fmt=mk, color=col, ms=3.5, capsize=1.5, elinewidth=0.7, mfc=col if k == 'cos' else 'white')
        out[f'{lab}.n_dom_cos'] = c['n_dom_cos']
    ax.axhline(0, color='black', lw=0.4); ax.set_xticks([0, 1]); ax.set_xticklabels(['aSAMt', 'AlphaFlow']); ax.set_xlim(-0.5, 1.5); ax.set_ylim(-0.15, 1.15)
    ax.set_ylabel('cosine with the heating shift')
    ax.legend(handles=[Line2D([], [], marker='s', ls='', mfc='white', color=OI['grey'], label='ceiling (two heating steps)'), Line2D([], [], marker='o', ls='', color='black', label='aSAMt error'), Line2D([], [], marker='o', ls='', color=OI['red'], label='AlphaFlow error'),
                       Line2D([], [], marker='v', ls='', mfc='white', color=OI['grey'], label='replicate-noise null')], frameon=False, fontsize=5.8, loc='center right')
    return save(fig, 'Fig5'), dict(n=len(rows), cos=out)

# ---------------- Fig 6 ----------------
def fig6():
    p1n, pw = J('r38_p1n.json'), J('r44_p1win.json')
    items = [('Q ≥ 0.70 (primary)', p1n['all']['U']['half|0.70']['U2'], 'm'), ('Q ≥ 0.75', p1n['all']['U']['half|0.75']['U2'], 'm'), ('Q ≥ 0.80', p1n['all']['U']['half|0.80']['U2'], 'm'),
             ('50–100 ns', pw['all']['each|w050_100'], 'w'), ('300–400 ns', pw['all']['each|w300_400'], 'w'), ('400–500 ns', pw['all']['each|w400_500'], 'w')]
    fig = plt.figure(figsize=(5.6, 2.9)); gs = fig.add_gridspec(1, 2, width_ratios=[2.2, 1], wspace=0.42, left=0.2, right=0.985, bottom=0.18, top=0.86)
    ax = fig.add_subplot(gs[0]); panel(ax, 'A', x=-0.38); y = np.arange(len(items))[::-1]; vals = {}
    for yi, (lab, d, kind) in zip(y, items):
        if kind == 'm': e, lo, hi, n = d['median'], d['lo'], d['hi'], d['n_sf']
        else: e, lo, hi, n = d['U2'], d['lo025'], d['hi975'], d['n_sf']
        vals[lab] = (e, lo, hi, n)
        ax.errorbar(e, yi, xerr=[[e - lo], [hi - e]], fmt='o', color='black' if yi == y[0] else OI['blue'] if kind == 'm' else OI['orange'], ms=3.5, capsize=1.5, elinewidth=0.8)
        ax.text(0.33, yi, f'{n}', va='center', fontsize=6, color='#555555')
    ax.text(0.33, y[0] + 0.55, 'super-\nfamilies', fontsize=5.5, color='#555555', va='bottom')
    ax.axvline(0, color='black', lw=0.5); ax.set_yticks(y); ax.set_yticklabels([i[0] for i in items]); ax.set_xlim(-0.75, 0.42)
    ax.set_xlabel('U2 (log spread ratio, controlled)'); ax.axhline(2.5, color=OI['grey'], lw=0.4, ls=':')
    ax.text(-0.73, 4.9, 'second half of\neach trajectory', fontsize=5.5, color='#555555', va='center'); ax.text(-0.73, 1.0, 'absolute\nwindows (Q ≥ 0.70)', fontsize=5.5, color='#555555', va='center')
    ax = fig.add_subplot(gs[1]); panel(ax, 'B', x=-0.55); fd = p1n['all']['fold']['379']
    ax.errorbar(0, fd['median'], yerr=[[fd['median'] - fd['lo']], [fd['hi'] - fd['median']]], fmt='o', color=OI['green'], ms=4, capsize=2)
    ax.axhline(0, color='black', lw=0.5); ax.set_xticks([0]); ax.set_xticklabels([f'379 K\n({fd["n_sf"]} superfamilies)']); ax.set_xlim(-0.6, 0.6)
    ax.set_ylabel('paired difference in folded-fraction logit')
    return save(fig, 'Fig6'), dict(U2=vals, fold=(fd['median'], fd['lo'], fd['hi'], fd['n_sf']))

# ---------------- Fig 7 ----------------
def fig7():
    p3, ex2, c, nl, pw, cr, cr2, f1 = J('r38_p3.json'), J('r38_p3_explore2.json'), J('r38_c.json'), J('r38_c_nl.json'), J('r44_pairw.json'), J('r44_corr.json')['by_corr'], J('r44_corr2.json'), J('r45_f1.json')
    g = lambda d: (d['est'], d['lo'], d['hi'])
    m = ex2['real_sf_only']['main']
    items = [('First test, pre-specified interval', g(p3['main']), 'first'), ('First test, superfamily-only interval (exploratory)', (m['est'], *m['ci_sf']), 'first'),
             ('C-1, structural control', g(c['C1']['main']), 'main'), ('C-2, network-model control as well', g(c['C2']['main']), 'main'),
             ('Nonlinear controls (eight covariates)', (nl['C_main']['est']['nl'], *nl['C_main']['ci95']['nl']), 'rob'),
             ('C-1, pairs weighted equally', g(pw['C_C1']['pair_weighted']), 'rob'), ('C-2, pairs weighted equally', g(pw['C_C2']['pair_weighted']), 'rob'),
             ('C-1, TM-align on mean points', g(cr2['tmmean']['full']['C1']), 'corr'), ('C-2, TM-align on mean points', g(cr2['tmmean']['full']['C2']), 'corr'),
             ('C-1, TM-align on medoids', g(cr['tm']['C1']), 'corr'), ('C-2, TM-align on medoids', g(cr['tm']['C2']), 'corr'),
             ('C-1, sequence alignments', g(cr['seq']['C1']), 'corr'), ('C-2, sequence alignments', g(cr['seq']['C2']), 'corr'),
             ('C-2, bond-direction profiles only', g(f1['C2']['tilt']), 'rob')]
    col = dict(first=OI['grey'], main='black', rob=OI['blue'], corr=OI['orange'])
    fig, ax = plt.subplots(figsize=(7.5, 4.0)); fig.subplots_adjust(left=0.37, right=0.8, bottom=0.25, top=0.97); y = np.arange(len(items))[::-1]; vals = {}
    for yi, (lab, (e, lo, hi), grp) in zip(y, items):
        vals[lab] = (e, lo, hi)
        ax.errorbar(e, yi, xerr=[[e - lo], [hi - e]], fmt='o', color=col[grp], ms=3.5, capsize=1.5, elinewidth=0.8, mfc='white' if 'exploratory' in lab else col[grp])
        ax.text(1.02, yi, f'{e:.3f} [{lo:.3f}, {hi:.3f}]'.replace('-', '−'), va='center', fontsize=6.3, transform=ax.get_yaxis_transform())
    ax.axvline(0, color='black', lw=0.5); ax.set_yticks(y); ax.set_yticklabels([i[0] for i in items], fontsize=6.8); ax.set_xlim(-0.1, 0.46)
    ax.set_xlabel('partial Spearman correlation of profile similarity ' + r'$r_{\rm dyn}$' + ' with sequence identity')
    for yy in (11.5, 9.5, 6.5): ax.axhline(yy, color=OI['grey'], lw=0.4, ls=':')
    ax.legend(handles=[Patch(color=col['first'], label='first homolog study (shorter domains)'), Patch(color=col['main'], label='confirmatory set'),
                       Patch(color=col['rob'], label='confirmatory set, further controls'), Patch(color=col['corr'], label='confirmatory set, other correspondences')],
              frameon=False, fontsize=6.3, loc='upper center', bbox_to_anchor=(0.3, -0.13), ncol=2)
    return save(fig, 'Fig7'), vals

# ---------------- S1 Fig ----------------
def s1fig():
    adv = J('r41_synth_adv.json')['generators']
    fams = [('i', 'ℓ2 terms,\ngap cost |a|'), ('ii', 'ℓ1 terms (ERP),\ngap cost |a|'), ('iii', 'planar curves mod SO(2),\ngap cost |a|'), ('iv', 'gap cost 2|a|\n(not 1-Lipschitz)'), ('v', 'dynamic time\nwarping')]
    fig, ax = plt.subplots(figsize=(7.0, 2.9)); fig.subplots_adjust(left=0.08, right=0.99, bottom=0.2, top=0.97); x = np.arange(len(fams)); w = 0.38; out = {}
    for k, (gname, col, lab) in enumerate([('A', OI['blue'], 'A: shrinking and deleting terms'), ('B', OI['orange'], 'B: stretching')]):
        v = []
        for f, _ in fams:
            d = adv[gname][f]; nv = d.get('violations', d.get('violations_certain')); v.append(100 * nv / d['n_ineq']); out[f'{gname}.{f}'] = (nv, d['n_ineq'], d.get('undecided'))
        ax.bar(x + (k - 0.5) * w, v, w, color=col, label=lab)
        for xi, vi in zip(x + (k - 0.5) * w, v): ax.text(xi, vi + 0.8, f'{vi:.1f}' if vi else '0', ha='center', va='bottom', fontsize=6)
    ax.text(2 - w / 2, 4.5, f"A: {adv['A']['iii']['undecided']} undecided\n(equality within\ncertification width)", ha='center', va='bottom', fontsize=5.5, color='#555555')
    ax.set_xticks(x); ax.set_xticklabels([f[1] for f in fams], fontsize=6.5); ax.set_ylabel('% of triangle inequalities violated'); ax.set_ylim(0, 44)
    ax.legend(frameon=False, fontsize=6.5, loc='upper left', bbox_to_anchor=(0.0, 0.88)); ax.axvline(2.5, color=OI['grey'], lw=0.4, ls=':')
    ax.text(1, 41, '1-Lipschitz gap costs (Lemma 1)', ha='center', fontsize=7); ax.text(3.5, 41, 'controls', ha='center', fontsize=7)
    return save(fig, 'S1_Fig'), out

# ---------------- Fig 1 (schematic) ----------------
from matplotlib.patches import FancyArrowPatch, Circle, FancyBboxPatch
SH, UN = OI['blue'], OI['orange']
def glyph(ax, x, y, th, s=1.0, col='black', L=0.5, lw=1.4, zorder=3):
    """2-D stand-in for a bond term x = s R / 2: arrow along f1 (length ∝ s) with a tick along f2."""
    if s <= 1e-6:
        ax.add_patch(Circle((x, y), 0.08, fill=False, ec=col, lw=0.9, zorder=zorder)); return
    dx, dy = np.cos(th) * L * s, np.sin(th) * L * s
    ax.add_patch(FancyArrowPatch((x - dx / 2, y - dy / 2), (x + dx / 2, y + dy / 2), arrowstyle='-|>', mutation_scale=5 * max(s, 0.4), color=col, lw=lw, zorder=zorder))
    px, py = -np.sin(th) * 0.2 * s, np.cos(th) * 0.2 * s
    ax.plot([x, x + px], [y, y + py], color=col, lw=lw * 0.8, zorder=zorder)

def fig1():
    fig = plt.figure(figsize=(7.5, 4.6)); gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.35], height_ratios=[1, 1], wspace=0.12, hspace=0.25, left=0.03, right=0.99, bottom=0.04, top=0.95)
    # A: bond frame
    ax = fig.add_subplot(gs[0, 0]); ax.set_xlim(0, 4.2); ax.set_ylim(0, 2.6); ax.set_aspect('equal'); ax.axis('off'); panel(ax, 'A', x=-0.02, y=0.98)
    ca1, ca2 = np.array([0.6, 0.9]), np.array([3.0, 1.3]); cC, oO, nN = np.array([1.55, 1.55]), np.array([1.45, 2.25]), np.array([2.15, 1.05])
    ax.plot(*zip(ca1, cC, nN, ca2), color='#777777', lw=1.0, zorder=1); ax.plot(*zip(cC, oO), color='#777777', lw=1.0, zorder=1)
    for p, lab, col, r in [(ca1, r'C$\alpha_i$', '#444444', 0.13), (ca2, r'C$\alpha_{i+1}$', '#444444', 0.13), (cC, 'C', '#444444', 0.09), (oO, 'O', OI['red'], 0.09), (nN, 'N', OI['blue'], 0.09)]:
        ax.add_patch(Circle(p, r, color=col, zorder=2)); ax.text(p[0] + (0.18 if lab != 'O' else 0.14), p[1] - 0.05 if lab not in ('C',) else p[1] + 0.12, lab, fontsize=7, va='center')
    f1 = (ca2 - ca1) / np.linalg.norm(ca2 - ca1); co = oO - cC; f2 = co - (co @ f1) * f1; f2 /= np.linalg.norm(f2); o = np.array([3.2, 0.3])
    for v, lab, dxy in [(f1, r'$f_1$', (0.06, -0.14)), (f2, r'$f_2$', (-0.24, 0.02))]:
        ax.add_patch(FancyArrowPatch(tuple(o), tuple(o + 0.6 * v), arrowstyle='-|>', mutation_scale=8, color='black', lw=1.1)); ax.text(*(o + 0.65 * v + np.array(dxy)), lab, fontsize=8)
    ax.add_patch(Circle(tuple(o), 0.07, fill=False, ec='black', lw=0.9)); ax.plot(*o, 'k.', ms=2.5); ax.text(o[0] - 0.3, o[1] - 0.15, r'$f_3$', fontsize=8); ax.text(o[0] + 0.25, o[1] - 0.32, 'bond frame (f₃ out of page)', fontsize=6, color='#555555')
    ax.text(3.2, 2.35, r'$x_i = s_i R_i/2$' + '\n' + r'$R_i = [f_1, f_2, f_3]$' + '\n' + r'$s_i = l_i\,/\,3.83$ Å', fontsize=7.5, va='top', ha='left')
    ax.text(0.15, 0.12, 'absent bond: the zero matrix', fontsize=6.5, color='#555555')
    # B: correspondence
    ax = fig.add_subplot(gs[0, 1]); ax.set_xlim(-1.0, 9.6); ax.set_ylim(-0.3, 2.7); ax.axis('off'); panel(ax, 'B', x=-0.02, y=0.98)
    tx = [0.2, 0.5, 0.9, 1.6, 2.2, 2.4, 2.9]; ty = [0.3, 0.6, 1.0, 2.0, 2.3, 2.9]
    m = {0: 0, 1: 1, 2: 2, 4: 3, 5: 4}  # x index -> y index (order-preserving); x3, x6 and y5 unmatched
    xs = np.arange(7) * 1.2; ys = np.r_[0, 1, 2, 4, 5, 7.3] * 1.2
    yx, yy = 2.1, 0.45
    for i, (xx, th) in enumerate(zip(xs, tx)): glyph(ax, xx, yx, th, col=SH if i in m else UN)
    U = 0.35
    for j, (xx, th) in enumerate(zip(ys, ty)): glyph(ax, xx, yy, th + U, col=SH if j in m.values() else UN)
    for i, j in m.items(): ax.plot([xs[i], ys[j]], [yx - 0.25, yy + 0.25], color=SH, lw=0.6, ls='-', alpha=0.7)
    for i in (3, 6):
        ax.add_patch(Circle((xs[i], yy + 0.05), 0.07, fill=False, ec=UN, lw=0.8)); ax.text(xs[i], yy - 0.27, '0', ha='center', fontsize=6.5, color=UN)
        ax.plot([xs[i], xs[i]], [yx - 0.25, yy + 0.15], color=UN, lw=0.6, ls=':')
    ax.add_patch(Circle((ys[5], yx - 0.0), 0.07, fill=False, ec=UN, lw=0.8)); ax.text(ys[5], yx + 0.2, '0', ha='center', fontsize=6.5, color=UN); ax.plot([ys[5], ys[5]], [yx - 0.1, yy + 0.25], color=UN, lw=0.6, ls=':')
    ax.text(-0.95, yx, r'$x$', fontsize=9, va='center'); ax.text(-0.95, yy, r'$Uy$', fontsize=9, va='center')
    ax.text(4.3, 2.62, 'order-preserving correspondence π and one global rotation U', ha='center', fontsize=7)
    # C: exact split
    ax = fig.add_subplot(gs[1, 0]); ax.set_xlim(0, 4.2); ax.set_ylim(0, 2.6); ax.axis('off'); panel(ax, 'C', x=-0.02, y=0.98)
    a, b = 2.3, 1.2
    ax.add_patch(FancyBboxPatch((0.25, 1.15), a, 0.45, boxstyle='square,pad=0', color=SH)); ax.add_patch(FancyBboxPatch((0.25 + a, 1.15), b, 0.45, boxstyle='square,pad=0', color=UN))
    ax.text(0.25 + (a + b) / 2, 1.85, r'$d(x, y)^2$', ha='center', fontsize=9)
    ax.text(0.25 + a / 2, 0.85, 'shared bonds\n' + r'$\sum_{(i,j)\in\pi}\|x_i - U y_j\|^2$', ha='center', va='top', fontsize=7, color=SH)
    ax.text(0.25 + a + b / 2, 0.85, 'unshared bonds\n' + r'$\sum\|x_i\|^2 + \sum\|y_j\|^2$', ha='center', va='top', fontsize=7, color=UN)
    # D: geodesic
    ax = fig.add_subplot(gs[1, 1]); ax.set_xlim(-1.0, 9.6); ax.set_ylim(-0.4, 2.6); ax.axis('off'); panel(ax, 'D', x=-0.02, y=0.98)
    ts = [0, 0.5, 1]; rows_y = [2.05, 1.15, 0.25]
    for t, yv in zip(ts, rows_y):
        ax.text(-0.95, yv, f't = {t:g}', fontsize=7, va='center')
        k = 0
        for i in range(7):
            if i in m:
                j = m[i]; th = (1 - t) * tx[i] + t * (ty[j] + U); glyph(ax, xs[i] * (1 - t) + ys[j] * t, yv, th, col=SH)
            else:
                glyph(ax, xs[i], yv, tx[i], s=1 - t, col=UN)
        glyph(ax, ys[5], yv, ty[5] + U, s=t, col=UN)
    ax.text(4.3, 2.5, 'shared terms move linearly; unshared terms shrink to or grow from zero', ha='center', fontsize=7)
    return save(fig, 'Fig1'), {}
