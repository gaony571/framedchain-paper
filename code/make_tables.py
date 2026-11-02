"""S2 Table and S3 Table from stored result files (no recomputation). Output: manuscript_figs/S2_Table.xlsx, S3_Table.xlsx (+ .csv)."""
import json, re, numpy as np, pandas as pd
R = '/content/drive/MyDrive/atlas_framed_curves'; OUT = R + '/manuscript_figs'
J = lambda f: json.load(open(f'{R}/results/{f}'))
f3 = lambda v: '' if v is None else f'{v:.3f}'
def eci(d, k='est', c='ci95', small=False):
    if d is None: return ''
    lo, hi = (d[c] if c in d else (d['lo'], d['hi']))
    nd = 4 if (small and abs(d[k]) < 0.01) else 3
    fm = lambda v: f'{v:.{nd}f}' if v else f'{0:.{nd}f}'   # a value in (−0.5, 0) × 10^−nd keeps its sign (−0.000; write() turns - into −); exact zero, also −0.0, prints 0.000
    return f'{fm(d[k])} [{fm(lo)}, {fm(hi)}]'
def eck(d):
    if d is None: return ''
    lo, hi = (d['ci95'] if 'ci95' in d else (d['lo'], d['hi']))
    return f"{d['est']:.4f} [{lo:.4f}, {hi:.4f}]"
def region(k): return 'near' if k < 0.5 else 'crossing' if k <= 2 else 'far'
CELL_EN = {'온도 +28 K (뒤 절반)': 'Heating +28 K (second half)', '온도 +59 K (뒤 절반)': 'Heating +59 K (second half)', '온도 +93 K (뒤 절반)': 'Heating +93 K (second half)',
           '온도 +130 K (뒤 절반)': 'Heating +130 K (second half)', 'R3′ (50–100 ns, 주)': 'Force field, temperature and protocol (ATLAS 300 K vs mdCATH 320 K, 50–100 ns)',
           'R3′ (전체 창, 참고)': 'Force field, temperature and protocol (whole trajectories, reference)', '7c aSAMt 320 K': 'aSAMt error at 320 K', '7c AlphaFlow': 'AlphaFlow error at 320 K',
           '상동체 TM 0.2-0.3': 'Homolog pairs, TM-score 0.2–0.3', '상동체 TM 0.3-0.4': 'Homolog pairs, TM-score 0.3–0.4', '상동체 TM 0.4-0.5': 'Homolog pairs, TM-score 0.4–0.5',
           '상동체 TM 0.5-0.6': 'Homolog pairs, TM-score 0.5–0.6', '상동체 TM 0.6-0.7': 'Homolog pairs, TM-score 0.6–0.7', '상동체 TM 0.7-0.8': 'Homolog pairs, TM-score 0.7–0.8',
           '상동체 TM 0.8-1.01': 'Homolog pairs, TM-score ≥ 0.8', '자연 삽입·결실 (R10a′ 주)': 'Natural insertion–deletion pairs'}

def s2():
    s1 = J('r43_sens.json')['s1']; rows = []
    for k, en in CELL_EN.items():
        c = s1['cells'].get(k)
        if c is None: continue
        rows.append({'Cell': en, 'Units (all)': c['n_units'], 'Median n̄ (all)': c['nbar_median'], 'κ (all)': eck(c['ratio']), 'Regime (all)': region(c['ratio']['est']),
                     'Units (80–150 bonds)': c['n_units_matched'], 'Median n̄ (matched)': c['nbar_median_matched'], 'κ (matched)': eck(c.get('ratio_matched')),
                     'Regime (matched)': region(c['ratio_matched']['est']) if c.get('ratio_matched') else '', 'Δ²_D/√n̄ (all)': eci(c['d2n'], small=True), 'Δ²_D/√n̄ (matched)': eci(c.get('d2n_matched'), small=True)})
    df = pd.DataFrame(rows)
    note = (f"κ = √((C/w̄)² − 1) and Δ²_D/√n̄ for each cell, as superfamily medians of per-unit values with 95% superfamily bootstrap intervals (2,000 resamples), over all units and over "
            f"units with mean bond count n̄ in [80, 150]. Units are domains (heating, model errors), chains (force-field cells) or domain pairs (homolog and insertion–deletion cells). "
            f"Regimes: near κ < 0.5, crossing 0.5–2, far κ > 2. Rank agreement of cells between the two sets: Kendall τ = {s1['kendall_tau_ratio_full_vs_matched']:.2f} over {s1['n_cells_compared']} cells (all cells with size-matched values except the whole-trajectory force-field reference); "
            f"no cell with size-matched values changed regime. Size-matched κ, regime and Δ²_D/√n̄ are given only for cells with at least five units in [80, 150] (13 of the 16 cells; blank otherwise). Values for all units were recomputed here with the same resamples as the matched sets, so their intervals can differ slightly from those of Fig 3A, which come from the calibration analysis (results/r39_E0.json). Source: results/r43_sens.json (s1), a check added after the results (record round43, 22 October 2026 08:54 UTC).")
    return df, note

def s3():
    a, b = J('r46_align_fix1.json'), J('r48_identity_range.json'); rows = []
    E, D, P, S2, AS = a['est'], a['descriptive'], a['paired_diff'], a['seed2'], a['argsort_descriptive']
    lab46 = [('base', 'Original (identity and r_dyn on the edit-distance correspondence)'), ('A', 'A: identity over residues aligned by TM-align on medoids'),
             ('B0', 'B0: precision of the correspondence as covariate (exact pairing)'), ('B2', 'B2: precision as covariate (within two bonds)')]
    for key, lab in lab46:
        rows.append({'Part': 'Alignment accuracy (confirmatory set, 4,152 pairs)', 'Analysis': lab, 'C-1': eci(E[f'{key}_C1']), 'C-2': eci(E[f'{key}_C2']),
                     'C-1, second seed': f"[{S2[f'{key}_C1']['lo']:.3f}, {S2[f'{key}_C1']['hi']:.3f}]", 'C-2, second seed': f"[{S2[f'{key}_C2']['lo']:.3f}, {S2[f'{key}_C2']['hi']:.3f}]",
                     'C-1, ties in sort order': eci(AS.get(f'{key}_C1')), 'C-2, ties in sort order': eci(AS.get(f'{key}_C2'))})
    for key, lab in [('Aprime', "A′: identity at the edit-distance positions, TM-align partners"), ('recall0', 'Recall instead of precision (exact)'), ('recall2', 'Recall instead of precision (within two bonds)')]:
        if f'{key}_C1' in D: rows.append({'Part': 'Alignment accuracy (confirmatory set, 4,152 pairs)', 'Analysis': lab, 'C-1': eci(D[f'{key}_C1']), 'C-2': eci(D[f'{key}_C2'])})
    for key, lab in [('base-A', 'Paired difference, original − A'), ('base-B0', 'Paired difference, original − B0'), ('base-B2', 'Paired difference, original − B2')]:
        rows.append({'Part': 'Alignment accuracy (confirmatory set, 4,152 pairs)', 'Analysis': lab, 'C-1': eci(P[f'{key}_C1']), 'C-2': eci(P[f'{key}_C2'])})
    for key, lab in [('cpl_ident_prec0', 'Partial correlation of precision (exact) with identity'), ('cpl_rdyn_prec0', 'Partial correlation of precision (exact) with r_dyn'),
                     ('cpl_ident_prec2', 'Partial correlation of precision (two bonds) with identity'), ('cpl_rdyn_prec2', 'Partial correlation of precision (two bonds) with r_dyn')]:
        rows.append({'Part': 'Alignment accuracy (confirmatory set, 4,152 pairs)', 'Analysis': lab + ' (structural and network-model controls)', 'C-2': eci(E[key])})   # both controls, as in C-2 (r46_align_fix1.py: Z2 = control, r_anm)
    sets = b['sets']; nsf = b['n_sf']; B = b['est']
    lab48 = [('full', f"All pairs ({sets['full']:,} pairs, {nsf['full']} superfamilies)"), ('S15', f"Identity below 15% ({sets['S15']:,} pairs, {nsf['S15']} superfamilies)"),
             ('SFbase', f"All pairs of the {nsf['S15']} superfamilies with a pair below 15% ({sets['SFbase']:,} pairs)"), ('S10', f"Identity below 10% ({sets['S10']:,} pairs)"),
             ('T15', f"TM-align identity below 15% ({sets['T15']:,} pairs)"), ('A15', 'Below 15%, identity from TM-align (A)'), ('B015', 'Below 15%, precision covariate (B0)')]
    for key, lab in lab48:
        rows.append({'Part': 'Identity range (average ranks, 10,000 resamples)', 'Analysis': lab, 'C-1': eci(B[f'{key}_C1']), 'C-2': eci(B[f'{key}_C2']),
                     'Pairs': B[f'{key}_C1']['n'], 'Superfamilies': B[f'{key}_C1']['n_sf']})
    dd = b['diff_full_minus_S15']; rows.append({'Part': 'Identity range (average ranks, 10,000 resamples)', 'Analysis': 'Paired difference, all pairs − below 15% (descriptive)', 'C-1': eci(dd['C1']), 'C-2': eci(dd['C2'])})
    rr, rv = b['range_restriction'], b['residual_var']
    rows.append({'Part': 'Identity range (average ranks, 10,000 resamples)', 'Analysis': 'Value expected from range restriction alone (mean of 1,000 simulations; reading)',
                 'C-1': f"{rr['C1']['expected_mean']:.3f} ({rr['C1']['reading']})", 'C-2': f"{rr['C2']['expected_mean']:.3f} ({rr['C2']['reading']})"})
    rows.append({'Part': 'Identity range (average ranks, 10,000 resamples)', 'Analysis': 'Residual variance of rank r_dyn, identity < 15% / ≥ 15% (×10⁶)',
                 'C-1': f"{rv['C1']['lt15'] / 1e6:.2f} / {rv['C1']['ge15'] / 1e6:.2f}", 'C-2': f"{rv['C2']['lt15'] / 1e6:.2f} / {rv['C2']['ge15'] / 1e6:.2f}"})
    for i, d in enumerate(b['deciles_C2'], 1):
        rows.append({'Part': 'Identity deciles (all pairs)', 'Analysis': f"Decile {i}: identity {d['ident_lo']:.3f}–{d['ident_hi']:.3f} ({d['n']} pairs)",
                     'C-2': f"{d['mean_rY']:.0f} [{d['ci95'][0]:.0f}, {d['ci95'][1]:.0f}]"})
    for key, lab in [('full', 'All pairs'), ('S15', 'Identity below 15%')]:
        c = b['composition'][key]; rows.append({'Part': 'Composition', 'Analysis': f"{lab}: median observed / composition-expected identity (median difference)",
                                                'C-1': f"{c['obs_median']:.3f} / {c['exp_median']:.3f} ({c['diff_median']:.3f})", 'Pairs': c['n']})
    df = pd.DataFrame(rows).fillna('')
    note = ("Partial Spearman correlations of r_dyn with sequence identity with the structural control (C-1) and with the network-model control as well (C-2), "
            "superfamily-weighted, with 95% superfamily-cluster bootstrap intervals. All rows use average ranks and 10,000 resamples (2,000 for the deciles), so the original estimates and the all-pairs estimates differ slightly from those in the text, which used ranks in sort order and 2,000 resamples. Alignment-accuracy rows: results/r46_align_fix1.json (second seed and "
            "ties ranked in sort order for comparison; the columns with ties in sort order were computed on a different processor from the estimates in the text, where ties were ordered differently, so their point estimates differ from them by up to 5 × 10⁻⁴). Identity-range rows: results/r48_identity_range.json (one bootstrap stream shared by all rows). Decile rows: in the C-2 column, the weighted mean residual rank of r_dyn after both controls (in rank units) within each identity decile of all pairs "
            "(2,000 resamples). Composition-expected identity: Σ_a p_a q_a over the first residues of the shared bonds. Both checks were added after the confirmatory result "
            "(records of 29 October 2026 16:34 UTC and 30 October 2026 07:30 UTC).")
    return df, note

NUMERIC = ('Median n̄ (all)', 'Median n̄ (matched)', 'Pairs', 'Superfamilies')   # S2 and S3 columns stored as numbers in the xlsx (the csv keeps the same text as before)
def fmt(path):
    """Column widths, wrapped text, bold frozen header (the same as fmt() in make_s1.py)."""
    from openpyxl import load_workbook
    from openpyxl.styles import Alignment, Font
    wb = load_workbook(path)
    for ws in wb.worksheets:
        for col in ws.columns:
            n = max(len(str(c.value)) if c.value is not None else 0 for c in col)
            ws.column_dimensions[col[0].column_letter].width = min(max(10, n * 0.9 + 2), 70 if ws.title not in ('Notes', 'Note') else 140)
            for c in col: c.alignment = Alignment(wrap_text=True, vertical='top')
        for c in ws[1]: c.font = Font(bold=True)
        ws.freeze_panes = 'A2'
    wb.save(path)

def write(name, df, note):
    num = [c for c in df.columns if c not in ('Part', 'Analysis', 'Cell')]
    df = df.copy(); df[num] = df[num].apply(lambda col: col.map(lambda v: re.sub(r'(?<![\w])-(?=\d)', '−', v) if isinstance(v, str) else v))
    fmtnum = lambda v: '' if v is None or (isinstance(v, float) and np.isnan(v)) else (f'{round(v, 4):g}' if isinstance(v, (float, np.floating)) else v)
    df[num] = df[num].apply(lambda col: col.map(fmtnum)).astype(object)
    xl = df.copy()
    for c in [c for c in NUMERIC if c in xl.columns]: xl[c] = pd.to_numeric(xl[c].map(lambda v: np.nan if v == '' else v))   # same values as the text; any other text raises
    with pd.ExcelWriter(f'{OUT}/{name}.xlsx') as xw:
        xl.to_excel(xw, index=False, sheet_name=name); pd.DataFrame({'Note': [note]}).to_excel(xw, index=False, sheet_name='Note')
    fmt(f'{OUT}/{name}.xlsx')
    df.to_csv(f'{OUT}/{name}.csv', index=False)
    return df
