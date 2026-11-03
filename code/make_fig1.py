"""Fig 1: framed chains and the edit distance, on two lysozyme fragments of different lengths.

Usage:  python make_fig1.py [OUT_DIR]        (default ./fig1_out)
Output: OUT_DIR/Fig1.tif (RGB, LZW, 600 dpi), Fig1_preview.png (150 dpi), fig1_lysozyme.json (the numbers in the figure).

Steps
1. Download PDB 1AKI (hen egg-white lysozyme) and 1LZ1 (human lysozyme) from RCSB; chain A, first altloc, heavy atoms.
2. Align the whole chains with framedchain (255 random starts) and take the stretch where exactly one bond is unmatched,
   with four matched bonds on each side and no gap in the residue numbers (one stretch: 1AKI 42-49, 1LZ1 42-50).
3. Align the two fragments again (511 random starts), certify the optimum by branch-and-bound, and lay them out as
   aligned columns; each row's term of d^2 is ||x_i - U y_j||^2 (shared) or the squared norm of the unmatched term.
4. Draw both fragments in one orientation (the human fragment after the rotation U) with PyMOL: the long axis of the
   hen fragment vertical, then turned about that axis by the angle (5-degree grid) that keeps every bond of both
   fragments longest when projected onto the page. Bonds are colored by the alignment (blue shared, orange unshared).
5. Draw each term by its first two columns projected onto the page (arrow and tick), with the midpoint of the
   geodesic between the aligned columns and its term of d^2.

Needs numpy, torch, scipy, matplotlib, Pillow, framedchain 0.4.3 and PyMOL (open-source 2.5). Text is Arial (8 pt or
larger) when it is installed (msttcorefonts), otherwise Liberation Sans. FIG1_PYMOL gives the PyMOL command
(default "pymol -cq"); on Colab (Python 3.13, Debian pymol 2.5 built for Python 3.12) use
FIG1_PYMOL="/usr/bin/python3.12 -m pymol -cq". A stand-in `imp` module (removed from Python 3.12, still imported by
pymol 2.5) is put on PYTHONPATH for the PyMOL run.
Run on Colab on 2026-11-04 to draw the figure of the manuscript.
"""
import glob, json, os, subprocess, sys, urllib.request

import numpy as np
import torch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.patheffects as pe
from matplotlib.patches import FancyArrowPatch, Circle
from PIL import Image
import framedchain as fc

OUT = sys.argv[1] if len(sys.argv) > 1 else 'fig1_out'
os.makedirs(OUT, exist_ok=True)
for f in glob.glob('/usr/share/fonts/truetype/msttcorefonts/[Aa]rial*.ttf'):
    fm.fontManager.addfont(f)
FONT = 'Arial' if any(f.name == 'Arial' for f in fm.fontManager.ttflist) else 'Liberation Sans'
plt.rcParams.update({'font.family': FONT, 'mathtext.fontset': 'custom', 'mathtext.rm': FONT, 'mathtext.it': f'{FONT}:italic', 'mathtext.bf': f'{FONT}:bold'})
SH, UN, UN_BOND, UN_TEXT = '#0072B2', '#E69F00', '#F0A000', '#C07800'     # shared, unshared (glyphs, bond rods, numbers)
FS_S, FS_M, FS_H, FS_T = 8, 8.5, 9, 10                                   # small text, titles, headers, total (pt)
AA = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L',
      'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}
TAN10 = np.tan(np.radians(10))      # half of PyMOL's field of view (20 degrees)
rgb = lambda h: [int(h[k:k + 2], 16) / 255 for k in (1, 3, 5)]


# ---------------- structures and alignment ----------------
def read_chain(pid):
    p = f'{OUT}/{pid}.pdb'
    if not os.path.exists(p):
        urllib.request.urlretrieve(f'https://files.rcsb.org/download/{pid}.pdb', p)
    res = {}
    for ln in open(p):
        if ln.startswith('ENDMDL'):
            break
        if ln.startswith('ATOM') and ln[21] == 'A' and ln[16] in ' A' and ln[76:78].strip() != 'H':
            r = int(ln[22:26])
            res.setdefault(r, {'name': ln[17:20], 'atoms': {}})['atoms'][ln[12:16].strip()] = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
    return [r for r in sorted(res) if all(k in res[r]['atoms'] for k in ('N', 'CA', 'C', 'O'))], res

def atoms(res, rr, a):
    return np.array([res[r]['atoms'][a] for r in rr])

def coords(res, rr):
    return fc.coordinates_from_backbone(atoms(res, rr, 'CA'), atoms(res, rr, 'C'), atoms(res, rr, 'O'))

(rsA, resA), (rsB, resB) = read_chain('1AKI'), read_chain('1LZ1')
al = fc.align(coords(resA, rsA), coords(resB, rsB), n_random=255)
pairs = sorted(zip(*[t.tolist() for t in torch.nonzero(al.match, as_tuple=True)]))
win = None
for k in range(4, len(pairs) - 4):
    (i0, j0), (i1, j1) = pairs[k - 1], pairs[k]
    runs = all(pairs[q + 1][0] - pairs[q][0] == 1 and pairs[q + 1][1] - pairs[q][1] == 1 for q in list(range(k - 4, k - 1)) + list(range(k, k + 3)))
    if runs and {i1 - i0, j1 - j0} == {1, 2}:
        ia, ja = pairs[k - 4]; ib, jb = pairs[k + 2]
        ra, rb = rsA[ia:ib + 2], rsB[ja:jb + 2]
        if ra[-1] - ra[0] == len(ra) - 1 and rb[-1] - rb[0] == len(rb) - 1:
            win = (ra, rb); break
rA, rB = win
fa, fb = coords(resA, rA), coords(resB, rB)
d2, a2 = fc.sq_distance(torch.as_tensor(fa), torch.as_tensor(fb), n_random=511)
_, _, cols = fc.aligned_columns(torch.as_tensor(fa), torch.as_tensor(fb), alignment=a2)
U = a2.rotation.numpy()
Ca = np.array([fa[p] if p >= 0 else np.zeros((3, 3)) for p, q in cols])
Cb = np.array([U @ fb[q] if q >= 0 else np.zeros((3, 3)) for p, q in cols])
contrib = [float(((Ca[c] - Cb[c]) ** 2).sum()) for c in range(len(cols))]
shared_rows = [p >= 0 and q >= 0 for p, q in cols]
cert = fc.bnb(fa, fb)
seq = lambda res, rr: ''.join(AA[res[r]['name']] for r in rr)
summary = dict(hen=dict(pdb='1AKI', residues=[rA[0], rA[-1]], sequence=seq(resA, rA)), human=dict(pdb='1LZ1', residues=[rB[0], rB[-1]], sequence=seq(resB, rB)),
               columns=cols, terms=contrib, d2=float(d2), shared=sum(v for v, s in zip(contrib, shared_rows) if s), unshared=sum(v for v, s in zip(contrib, shared_rows) if not s),
               certified=bool(cert['certified']), S_lb=float(cert['S_lb']), S_ub=float(cert['S_ub']), tol=float(cert['tol']), rotation_U=U.tolist())

# ---------------- display orientation ----------------
caA, caB = atoms(resA, rA, 'CA'), atoms(resB, rB, 'CA')
cA, cB = caA.mean(0), caB.mean(0)
_, _, Vt = np.linalg.svd(caA - cA)
Q0 = np.array([Vt[1], Vt[0], np.cross(Vt[1], Vt[0])])               # long axis of the hen fragment vertical
if (Q0 @ (caA[0] - cA))[1] < 0:
    Q0 = np.array([-Vt[1], -Vt[0], np.cross(-Vt[1], -Vt[0])])        # first residue at the top
def turned(th):
    return np.array([[np.cos(th), 0, np.sin(th)], [0, 1, 0], [-np.sin(th), 0, np.cos(th)]]) @ Q0
def min_projected(th):
    Qt = turned(th)
    return min(np.linalg.norm((Qt @ np.diff(caA, axis=0).T)[:2], axis=0).min(), np.linalg.norm((Qt @ U @ np.diff(caB, axis=0).T)[:2], axis=0).min())
ths = np.radians(np.arange(-90, 91, 5))
theta = ths[max(range(len(ths)), key=lambda k: min_projected(ths[k]))]
Q = turned(theta)
summary['view_turn_deg'] = float(np.degrees(theta))
TA = lambda v: (Q @ (np.asarray(v) - cA).T).T
TB = lambda v: (Q @ U @ (np.asarray(v) - cB).T).T
ca2A, ca2B = TA(caA), TB(caB)

# ---------------- PyMOL renderings ----------------
def write_pdb(path, res, rr, T):
    with open(path, 'w') as f:
        k = 1
        for r in rr:
            for a, xyz in res[r]['atoms'].items():
                p = T(xyz[None])[0]
                f.write(f"ATOM  {k:5d} {a:<4s} {res[r]['name']} A{r:4d}    {p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00  0.00          {a[0]:>2s}\n"); k += 1

def bond_rods(ca, shared):
    out = []
    for i in range(len(ca) - 1):
        c = rgb(SH if i in shared else UN_BOND)
        out += [9.0, *map(float, ca[i]), *map(float, ca[i + 1]), 0.36 if i in shared else 0.44, *c, *c]   # CGO CYLINDER
    return out

W_PX, H_PX = 1800, 2800
hy = max(np.abs(ca2A[:, 1]).max(), np.abs(ca2B[:, 1]).max()) + 2.2
hx = max(np.abs(ca2A[:, 0]).max(), np.abs(ca2B[:, 0]).max()) + 2.2
D = max(hy, hx * H_PX / W_PX) / TAN10                                # camera distance: both fragments at one scale

def render(pdb, cgo, out):
    cgo = [round(float(v), 4) for v in cgo]
    script = f'''
from pymol import cmd
cmd.load("{pdb}", "p"); cmd.hide("everything")
cmd.bg_color("white"); cmd.set("ray_opaque_background", 1); cmd.set("antialias", 2); cmd.set("ray_shadows", 0); cmd.set("orthoscopic", 1); cmd.set("field_of_view", 20)
cmd.set("stick_radius", 0.13); cmd.set("sphere_scale", 0.26); cmd.set("stick_transparency", 0.45)
cmd.color("grey70", "elem C"); cmd.color("tv_blue", "elem N"); cmd.color("tv_red", "elem O"); cmd.color("tv_yellow", "elem S")
cmd.show("sticks", "all"); cmd.show("spheres", "(all) and name CA"); cmd.color("grey25", "name CA")
cmd.load_cgo({cgo}, "vb")
cmd.set_view((1,0,0, 0,1,0, 0,0,1, 0,0,{-D:.3f}, 0.000,0.000,0.000, {max(D - 30, 1):.3f},{D + 30:.3f}, 1))
cmd.ray({W_PX}, {H_PX}); cmd.png("{out}")
'''
    open(f'{OUT}/_render.py', 'w').write(script)
    os.makedirs(f'{OUT}/_shim', exist_ok=True)
    open(f'{OUT}/_shim/imp.py', 'w').write('import importlib.util\ndef find_module(name, path=None):\n    spec = importlib.util.find_spec(name)\n'
                                           '    if spec is None:\n        raise ImportError(name)\n    loc = spec.submodule_search_locations\n'
                                           "    return (None, list(loc)[0] if loc else spec.origin, ('', '', 5))\n")
    r = subprocess.run(os.environ.get('FIG1_PYMOL', 'pymol -cq').split() + [f'{OUT}/_render.py'], capture_output=True, text=True, timeout=600,
                       env={**os.environ, 'PYTHONPATH': os.path.abspath(f'{OUT}/_shim')})
    assert r.returncode == 0 and os.path.exists(out), (r.stdout + r.stderr)[-2000:]

shA = {p for p, q in cols if p >= 0 and q >= 0}; shB = {q for p, q in cols if p >= 0 and q >= 0}
write_pdb(f'{OUT}/fragA.pdb', resA, rA, TA); write_pdb(f'{OUT}/fragB.pdb', resB, rB, TB)
render(f'{OUT}/fragA.pdb', bond_rods(ca2A, shA), f'{OUT}/_hen.png')
render(f'{OUT}/fragB.pdb', bond_rods(ca2B, shB), f'{OUT}/_human.png')

def crop_pair(pa, pb, pad=(80, 220, 80, 80)):       # one box size for both images (same scale); pad left, right, top, bottom
    ims = [np.array(Image.open(p).convert('RGB')) for p in (pa, pb)]
    bbs = []
    for im in ims:
        ys, xs = np.where((im < 240).any(-1)); bbs.append((xs.min(), xs.max(), ys.min(), ys.max()))
    wc = max(b[1] - b[0] for b in bbs) + pad[0] + pad[1]; hc = max(b[3] - b[2] for b in bbs) + pad[2] + pad[3]
    out = []
    for im, (x0, x1, y0, y1) in zip(ims, bbs):
        cx = (x0 + x1) / 2 + (pad[1] - pad[0]) / 2; cy = (y0 + y1) / 2 + (pad[3] - pad[2]) / 2
        L_, T_ = int(round(cx - wc / 2)), int(round(cy - hc / 2))
        big = np.pad(im, ((hc, hc), (wc, wc), (0, 0)), constant_values=255)
        out.append((big[T_ + hc:T_ + 2 * hc, L_ + wc:L_ + 2 * wc], (L_, T_)))
    return out, (wc, hc)
((imA, offA), (imB, offB)), (WC, HC) = crop_pair(f'{OUT}/_hen.png', f'{OUT}/_human.png')

def to_px(p):                                        # model (after Q) -> pixel of the orthoscopic identity view
    s = (H_PX / 2) / (D * TAN10); return W_PX / 2 + p[0] * s, H_PX / 2 - p[1] * s

# ---------------- figure ----------------
def glyph(ax, x, y, M, col, L=0.62, lw=1.5):
    """A 3x3 term by its first two columns projected onto the page: arrow (first) and tick (second); zero: open circle."""
    P = Q @ M
    if np.abs(P).sum() < 1e-9:
        ax.add_patch(Circle((x, y), 0.075, fill=False, ec=col, lw=1.0, zorder=3)); return
    a = 2 * P[:2, 0] * L; t = 2 * P[:2, 1] * L * 0.45
    ax.add_patch(FancyArrowPatch((x - a[0] / 2, y - a[1] / 2), (x + a[0] / 2, y + a[1] / 2), arrowstyle='-|>', mutation_scale=6, color=col, lw=lw, zorder=3, shrinkA=0, shrinkB=0))
    ax.plot([x, x + t[0]], [y, y + t[1]], color=col, lw=lw * 0.8, zorder=3, solid_capstyle='round')

X0, XC, DX, WN = 0.0, 1.95, 0.36, 0.80
XG = (XC - 0.42 + XC + DX + WN) / 2                  # centre of the midpoint glyph and its number
X2 = 2 * XG - X0
def columns(ax):
    L = len(cols); dy = -0.16
    for c, (p, q) in enumerate(cols):
        y = L - 1 - c; shared = p >= 0 and q >= 0; col = SH if shared else UN
        glyph(ax, X0, y, Ca[c], col); glyph(ax, XC, y, (Ca[c] + Cb[c]) / 2, col); glyph(ax, X2, y, Cb[c], col)
        ax.text(XC + DX, y + dy, f'{contrib[c]:.3f}', ha='left', va='center', fontsize=FS_S, color=SH if shared else UN_TEXT)
        for A_, B_ in ((X0 + 0.40, XC - 0.42), (XC + DX + WN, X2 - 0.40)):
            ax.plot([A_, B_], [y, y], color=col, lw=0.7, ls='-' if shared else ':', alpha=0.8, zorder=1)
        if p >= 0: ax.text(X0 - 0.55, y, f'{rA[p]}–{rA[p] + 1}', ha='right', va='center', fontsize=FS_S, color='#444444')
        else: ax.text(X0 - 0.22, y - 0.02, '0', ha='right', va='center', fontsize=FS_H, color=UN)
        if q >= 0: ax.text(X2 + 0.55, y, f'{rB[q]}–{rB[q] + 1}', ha='left', va='center', fontsize=FS_S, color='#444444')
        else: ax.text(X2 + 0.22, y - 0.02, '0', ha='left', va='center', fontsize=FS_H, color=UN)
    ax.text(X0, L - 0.15, r'$x$' + '\nhen', ha='center', va='bottom', fontsize=FS_H); ax.text(X2, L - 0.15, r'$Uy$' + '\nhuman', ha='center', va='bottom', fontsize=FS_H)
    ax.text(XC, L - 0.15, 'midpoint\n' + r'$\gamma(1/2)$', ha='center', va='bottom', fontsize=FS_H)
    ax.text(XC + DX, L - 0.32, 'term of ' + r'$d$' + '²', ha='left', va='center', fontsize=FS_S, color='#555555')
    ax.plot([XC + DX - 0.02, XC + DX + WN - 0.08], [-0.45, -0.45], color='#333333', lw=0.6)
    ax.text(XG, -0.95, r'$d(x, y)$' + '² = ' + f'{sum(contrib):.3f}', ha='center', va='center', fontsize=FS_T)
    half = max(XG - (X0 - 1.45), (X2 + 1.45) - XG)
    ax.set_xlim(XG - half, XG + half); ax.set_ylim(-1.5, L + 1.1); ax.set_aspect('equal'); ax.axis('off')

STROKE = [pe.withStroke(linewidth=2.2, foreground='white')]
def structure(ax, im, off, title, ca, rr, mark=()):
    ax.imshow(im, interpolation='lanczos'); ax.axis('off'); ax.set_title(title, fontsize=FS_M, y=-0.02, va='top')
    for k in (0, len(rr) - 1, *mark):
        x, y = to_px(ca[k]); ax.text(x - off[0] + 120, y - off[1], str(rr[k]), fontsize=FS_S, color='#222222', va='center', path_effects=STROKE)

FW, FH = 7.5, 4.9
fig = plt.figure(figsize=(FW, FH))
axm = fig.add_axes([0.25, 0.0, 0.50, 1.0]); columns(axm)
ws = 0.235; hs = ws * FW / FH * HC / WC
tA = f"hen lysozyme (1AKI) {rA[0]}–{rA[-1]}\n{summary['hen']['sequence']}"
tB = f"human lysozyme (1LZ1) {rB[0]}–{rB[-1]}\n{summary['human']['sequence']}, rotated by U"
axL = fig.add_axes([0.125 - ws / 2, 0.2, ws, hs]); structure(axL, imA, offA, tA, ca2A, rA)
axR = fig.add_axes([0.875 - ws / 2, 0.2, ws, hs]); structure(axR, imB, offB, tB, ca2B, rB, mark=tuple(q for p, q in cols if p < 0))
fig.canvas.draw(); rend = fig.canvas.get_renderer(); inv = fig.transFigure.inverted()
bbs = [t.get_window_extent(rend).transformed(inv) for t in axm.texts]
target = (max(b.y1 for b in bbs) + min(b.y0 for b in bbs)) / 2       # vertical centre of the middle block (headers to total)
for a in (axL, axR):                                                 # centre each structure with its title on it
    b = a.get_tightbbox(rend).transformed(inv); pos = a.get_position()
    a.set_position([pos.x0, pos.y0 + target - (b.y0 + b.y1) / 2, pos.width, pos.height])
fig.savefig(f'{OUT}/_Fig1.png', dpi=600, facecolor='white'); fig.savefig(f'{OUT}/Fig1_preview.png', dpi=150, facecolor='white'); plt.close(fig)
Image.open(f'{OUT}/_Fig1.png').convert('RGB').save(f'{OUT}/Fig1.tif', compression='tiff_lzw', dpi=(600, 600)); os.remove(f'{OUT}/_Fig1.png')
json.dump(summary, open(f'{OUT}/fig1_lysozyme.json', 'w'), indent=1)
print(f"1AKI {rA[0]}-{rA[-1]} {summary['hen']['sequence']} | 1LZ1 {rB[0]}-{rB[-1]} {summary['human']['sequence']} | d2 {summary['d2']:.4f} = "
      f"{summary['shared']:.4f} + {summary['unshared']:.4f} | certified {summary['certified']} | view turn {summary['view_turn_deg']:.0f} deg | font {FONT}")
