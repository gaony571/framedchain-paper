"""Build the data deposit for the paper from a project root (docs/DATA.md, sections 2 and 3).

Usage:
    python tools/pack_for_zenodo.py ROOT OUT_DIR [--ensembles] [--dry-run]

Writes OUT_DIR/framedchain_paper_data.tar (uncompressed: most bytes are npz/xtc, which do not compress) and
OUT_DIR/MANIFEST_data.tsv (md5, bytes, path) for every packed file, plus OUT_DIR/EXCLUDED_public.tsv listing the
public inputs found under ROOT/raw by name, size and md5 so that the same files can be identified without
redistributing them. --ensembles adds the generated aSAMt and AlphaFlow ensembles (large). --dry-run only
prints the sizes per group. Nothing under ROOT is modified.
"""
import argparse, fnmatch, hashlib, os, sys, tarfile, time

# (group, glob relative to ROOT); a file is packed if it matches any pattern of an enabled group
GROUPS = [
    ('results', 'results/*'),
    ('results', 'results/_backup/**'),
    ('results', 'results/logs/**'),
    ('intermediate', 'raw/asamt_test90_needs.json'),
    ('intermediate', 'raw/mdcath_all_domains_listing.json'),
    ('intermediate', 'raw/mfm_sm/*.npz'),
    ('intermediate', 'derived/r38/*.txt'),
    ('intermediate', 'derived/r38/alphaflow_input.csv'),
    ('intermediate', 'derived/r38/p1n/pairs.json'),
    ('intermediate', 'derived/xprot_structs.npz'),
    ('intermediate', 'derived/atlas_pdb/**'),
    ('intermediate', 'derived/pu_fold208_ca.npz'),
    ('intermediate', 'derived/pu_blockB_ca.npz'),
    ('perdomain', 'derived/**/*.json'),
    ('perdomain', 'derived/**/*.jsonl'),
    ('perdomain', 'derived/**/*.tsv'),
    ('ensembles', 'derived/r38/asamt/**'),
    ('ensembles', 'derived/r38/alphaflow/**'),
]
# public inputs: never packed, only listed (licences of some sources do not allow redistribution)
PUBLIC = ['raw/mdcath/**', 'raw/mdcath_c/**', 'raw/mdcath_p1n/**', 'raw/atlas/**', 'raw/atlas_af82/**',
          'raw/models/**', 'raw/bioemu_cath/**', 'raw/megasim/**', 'raw/tsuboyama/**', 'raw/majewski/**',
          'raw/cath/**', 'raw/chaining/**']


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


def walk(root):
    for d, _, fs in os.walk(root):
        for f in fs:
            p = os.path.join(d, f)
            yield os.path.relpath(p, root).replace(os.sep, '/'), p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('root'); ap.add_argument('out')
    ap.add_argument('--ensembles', action='store_true'); ap.add_argument('--dry-run', action='store_true')
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
    if a.dry_run:
        return
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    with open(os.path.join(a.out, 'EXCLUDED_public.tsv'), 'w') as f:
        f.write('md5\tbytes\tpath\n')
        for rel, p in public:
            f.write(f'{md5(p)}\t{os.path.getsize(p)}\t{rel}\n')
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
