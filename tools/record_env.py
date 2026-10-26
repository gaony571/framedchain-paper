"""Record the software environment a script runs in (the base Colab environment was not frozen at run time;
use this when reproducing).

Usage:
    python tools/record_env.py OUT.json [ROOT]

Writes Python and platform details, every installed distribution with its version, the versions of the
packages the analysis imports, CUDA/GPU details if torch can see a GPU, `nvidia-smi` output if available, and,
when ROOT is given, the git commit of ROOT/code/ext/alphaflow and the md5 of the framespace bundles.
"""
import hashlib, json, os, platform, subprocess, sys, time
from importlib import metadata

KEY = ['numpy', 'scipy', 'pandas', 'torch', 'h5py', 'mdtraj', 'pot', 'numba', 'tmtools', 'biopython']


def run(cmd):
    try:
        return subprocess.run(cmd, capture_output=True, text=True, timeout=60).stdout.strip()
    except Exception as e:
        return f'unavailable: {e!r}'


def bundle_md5(d):
    """md5 of the concatenation of the .py files of a package directory in sorted order (as in docs/PROVENANCE.md)."""
    h = hashlib.md5()
    for f in sorted(x for x in os.listdir(d) if x.endswith('.py')):
        with open(os.path.join(d, f), 'rb') as fh:
            h.update(fh.read())
    return h.hexdigest()


def main():
    out = sys.argv[1]; root = sys.argv[2] if len(sys.argv) > 2 else None
    dists = sorted({(d.metadata['Name'] or '').lower(): d.version for d in metadata.distributions()}.items())
    rec = dict(written_utc=time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime()), python=sys.version,
               executable=sys.executable, platform=platform.platform(), machine=platform.machine(),
               key_packages={k: v for k, v in dists if k in KEY}, distributions=dict(dists),
               nvidia_smi=run(['nvidia-smi']))
    try:
        import torch
        rec['torch_cuda'] = dict(version=torch.version.cuda, available=torch.cuda.is_available(),
                                 device=torch.cuda.get_device_name(0) if torch.cuda.is_available() else None)
    except Exception as e:
        rec['torch_cuda'] = f'unavailable: {e!r}'
    if root:
        af = os.path.join(root, 'code', 'ext', 'alphaflow')
        rec['alphaflow_commit'] = run(['git', '-C', af, 'rev-parse', 'HEAD']) if os.path.isdir(af) else None
        for name, d in (('framespace_pkg', 'framespace_pkg/framespace'), ('stage0', 'framespace_wip/stage0')):
            p = os.path.join(root, d)
            rec[f'{name}_bundle_md5'] = bundle_md5(p) if os.path.isdir(p) else None
    with open(out, 'w') as f:
        json.dump(rec, f, indent=1)
    print(json.dumps({k: rec[k] for k in ('python', 'key_packages', 'torch_cuda')}, indent=1, default=str))


if __name__ == '__main__':
    main()
