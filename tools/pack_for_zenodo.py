"""Build the data deposit for the paper from a project root (docs/DATA.md, sections 2 and 3).

Usage:
    python tools/pack_for_zenodo.py ROOT OUT_DIR [--ensembles] [--dry-run | --excluded-only]

Writes OUT_DIR/framedchain_paper_data.tar (uncompressed: most bytes are npz/xtc, which do not compress) and
OUT_DIR/MANIFEST_data.tsv (md5, bytes, path) for every packed file, plus OUT_DIR/EXCLUDED_public.tsv listing the
public inputs matched by PUBLIC below, and the structures and tables taken from ATLAS (CC BY-NC 4.0), by name and size,
with the md5 of the files of up to 64 MiB (PUBLIC_MD5_MAX) that are not simulation data (lists, tables, PDB entries,
structures taken from ATLAS), so that the same files can be identified without redistributing them; simulation files
are identified by name and size within their versioned public releases. Other files under ROOT/raw are inputs of
analyses outside the paper that no script of this repository reads; they are neither packed nor listed. --ensembles
adds the generated aSAMt and AlphaFlow ensembles (large). --dry-run only prints the sizes per group. --excluded-only rewrites only EXCLUDED_public.tsv (for a
deposit whose tar is already uploaded). Nothing under ROOT is modified.
"""
import argparse, fnmatch, hashlib, os, sys, tarfile, time

PUBLIC_MD5_MAX = 64 * 2 ** 20   # public files up to this size get an md5 in EXCLUDED_public.tsv, except simulation data
SIMULATION = ('raw/mdcath', 'raw/atlas', 'raw/models/', 'raw/bioemu_cath/', 'raw/megasim/', 'raw/majewski/', 'raw/mfm_sm/')

# (group, glob relative to ROOT); a file is packed if it matches any pattern of an enabled group
GROUPS = [
    ('results', 'results/*'),
    ('results', 'results/_backup/**'),
    ('results', 'results/logs/**'),
    ('results', 'results/_registration_evidence_2026-10-26/**'),
    ('intermediate', 'raw/asamt_test90_needs.json'),
    ('intermediate', 'raw/mdcath_all_domains_listing.json'),
    ('intermediate', 'derived/r38/*.txt'),
    ('intermediate', 'derived/r38/alphaflow_input.csv'),
    ('intermediate', 'derived/r38/p1n/pairs.json'),
    ('perdomain', 'derived/**/*.json'),
    ('perdomain', 'derived/**/*.jsonl'),
    ('perdomain', 'derived/**/*.tsv'),
    ('ensembles', 'derived/r38/asamt/**'),
    ('ensembles', 'derived/r38/alphaflow/**'),
]
# public inputs: never packed, only listed (licences of some sources do not allow redistribution)
PUBLIC = ['raw/mdcath/**', 'raw/mdcath_c/**', 'raw/mdcath_p1n/**', 'raw/atlas/**', 'raw/atlas_af82/**',
          'raw/models/**', 'raw/bioemu_cath/**', 'raw/megasim/**', 'raw/tsuboyama/**', 'raw/majewski/**',
          'raw/cath/**', 'raw/chaining/**',
          # organism and growth-temperature sources of results/pilot_mdcath_organism_temp_v2.tsv (SIFTS, TEMPURA, Engqvist 2018)
          'raw/sifts/**', 'raw/tempura/**', 'raw/engqvist/**',
          # MFM reference MD converted to npz: read only by the ladder analyses, which the paper does not report
          'raw/mfm_sm/**',
          # taken from ATLAS (CC BY-NC 4.0): listed, not packed
          'raw/atlas_r6/**', 'derived/atlas_pdb/**', 'derived/xprot_structs.npz', 'derived/pu_fold208_ca.npz', 'derived/pu_blockB_ca.npz',
          'results/atlas_info_full.tsv']


def match(rel, pat):
    if pat.endswith('/**'):
        return rel.startswith(pat[:-2])
    if '**/' in pat:
        head, tail = pat.split('**/', 1)
        return rel.startswith(head) and fnmatch.fnmatch(os.path.basename(rel), tail)
    return fnmatch.fnmatch(rel, pat) and rel.count('/') == pat.count('/')


def md5(path, buf=1 << 22):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(buf), b''):
            h.update(b)
    return h.hexdigest()


def want_md5(rel, p):
    return not rel.startswith(SIMULATION) and os.path.getsize(p) <= PUBLIC_MD5_MAX


def walk(root):
    for d, _, fs in os.walk(root):
        for f in fs:
            p = os.path.join(d, f)
            yield os.path.relpath(p, root).replace(os.sep, '/'), p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('root'); ap.add_argument('out')
    ap.add_argument('--ensembles', action='store_true'); ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--excluded-only', action='store_true')
    a = ap.parse_args()
    enabled = {'results', 'intermediate', 'perdomain'} | ({'ensembles'} if a.ensembles else set())
    pack, public, size = [], [], {}
    for rel, p in sorted(walk(a.root)):
        if rel.startswith(('code/', 'framespace_pkg/', 'framespace_wip/')):
            continue                                   # code is archived with the code repository
        if any(match(rel, q) for q in PUBLIC):
            public.append((rel, p)); continue
        g = next((g for g, q in GROUPS if g in enabled and match(rel, q)), None)
        if g:
            pack.append((g, rel, p)); size[g] = size.get(g, 0) + os.path.getsize(p)
    for g in sorted(size):
        print(f'{g:13s} {sum(1 for x in pack if x[0] == g):7d} files {size[g] / 2 ** 30:9.2f} GiB')
    print(f'{"public (list)":13s} {len(public):7d} files {sum(os.path.getsize(p) for _, p in public) / 2 ** 30:9.2f} GiB')
    small = [os.path.getsize(p) for rel, p in public if want_md5(rel, p)]
    print(f'{"  with md5":13s} {len(small):7d} files {sum(small) / 2 ** 30:9.2f} GiB (not simulation data, up to {PUBLIC_MD5_MAX // 2 ** 20} MiB each)')
    if a.dry_run:
        return
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    with open(os.path.join(a.out, 'EXCLUDED_public.tsv'), 'w') as f:
        f.write('md5\tbytes\tpath\n')
        for rel, p in public:
            f.write(f'{md5(p) if want_md5(rel, p) else "-"}\t{os.path.getsize(p)}\t{rel}\n')
    if a.excluded_only:
        print('EXCLUDED_public.tsv only', len(public), 'files', flush=True)
        return
    man = os.path.join(a.out, 'MANIFEST_data.tsv')
    with tarfile.open(os.path.join(a.out, 'framedchain_paper_data.tar'), 'w') as tar, open(man, 'w') as f:
        f.write('md5\tbytes\tgroup\tpath\n')
        for i, (g, rel, p) in enumerate(pack):
            f.write(f'{md5(p)}\t{os.path.getsize(p)}\t{g}\t{rel}\n')
            tar.add(p, arcname=rel)
            if (i + 1) % 1000 == 0:
                print(time.strftime('%H:%M:%S'), i + 1, 'of', len(pack), flush=True)
        f.flush()
        tar.add(man, arcname='MANIFEST_data.tsv')
    print('DONE', len(pack), 'files', round(time.time() - t0), 's', flush=True)


if __name__ == '__main__':
    sys.exit(main())
