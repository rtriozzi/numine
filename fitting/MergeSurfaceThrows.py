#!/usr/bin/env python3

"""
  Merge per-throw PROfit surfaces into one file for the Brazil-band code.
  Usage: python3 merge_brazil.py [indir] [outfile]
  Bypasses PROfit's own combination step; preserves each input surface verbatim.
"""
import sys, os, re, glob
import numpy as np
import uproot

indir   = sys.argv[1] if len(sys.argv) > 1 else \
    '/storage/gpfs_data/icarus/plain/user/fpoppiicarus/RiccardoPROfit/NuMI_3P1A_DetSysts_2DFit_v2_3_0'
outpath = sys.argv[2] if len(sys.argv) > 2 else 'combined_brazil_surf.root'
pattern = '3P1A_brazil_*_surf.root'

def idx(p):
    m = re.search(r'3P1A_brazil_(\d+)_surf\.root', os.path.basename(p))
    return int(m.group(1)) if m else -1

def find_th2_key(f):
    th2 = [k for k, c in f.classnames().items() if c.startswith('TH2')]
    if not th2:
        raise RuntimeError('no TH2 found')
    pref = [k for k in th2 if 'surf' in k.lower()]
    return (pref or th2)[0]

files = sorted(glob.glob(os.path.join(indir, pattern)), key=idx)
if not files:
    sys.exit(f'no files matching {pattern} in {indir}')

ref_xe = ref_ye = ref_shape = None
n_written = 0
with uproot.recreate(outpath) as out:
    for path in files:
        with uproot.open(path) as f:
            key = find_th2_key(f)
            vals, xe, ye = f[key].to_numpy(flow=False)   # (nx, ny), x=axis0
        if ref_xe is None:
            ref_xe, ref_ye, ref_shape = xe, ye, vals.shape
            print(f'grid: axis0 {xe[0]:.3g}..{xe[-1]:.3g}  '
                  f'axis1 {ye[0]:.3g}..{ye[-1]:.3g}  shape {vals.shape}')
        else:
            if not (np.allclose(xe, ref_xe) and np.allclose(ye, ref_ye)):
                print(f'  SKIP {os.path.basename(path)}: grid mismatch'); continue
            if vals.shape != ref_shape:
                print(f'  SKIP {os.path.basename(path)}: shape mismatch'); continue
        out[f'brazil_throw_surf_{idx(path)}'] = (vals, xe, ye)
        n_written += 1

print(f'wrote {n_written}/{len(files)} throws to {outpath}')