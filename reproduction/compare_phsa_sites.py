"""Generate upstream PH-SA sites and compare to supplied Cu(100) input."""
import argparse
import json
import os
import sys
import time
from pathlib import Path

for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    os.environ[key] = '4'

import numpy as np
from ase.io import read
from ase.geometry import find_mic

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'reproduction/PH-SA_upstream/PH-SA'))
from phsa.adsorption_sites import SlabAdsorptionsSitesFinder


def plain(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def distances(source, target, cell):
    if not len(source) or not len(target):
        return {'within_0.1_A': 0, 'count': len(source)}
    nearest = []
    for point in source:
        _, length = find_mic(np.asarray(target) - point, cell, pbc=[True, True, False])
        nearest.append(float(length.min()))
    return {'count': len(source), 'within_0.1_A': sum(d <= 0.1 for d in nearest),
            'within_0.0001_A': sum(d <= 0.0001 for d in nearest),
            'median_nearest_A': float(np.median(nearest)), 'max_nearest_A': max(nearest)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--k', type=float, default=1.1)
    parser.add_argument('--co-anchor', choices=['C', 'CO'], default='CO')
    parser.add_argument('--surface-only', action='store_true')
    args = parser.parse_args()
    if not np.isfinite(args.k) or args.k <= 0:
        parser.error('k must be finite and positive')
    args.output_dir.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    atoms = read(ROOT / 'GCMC/GCMC/POSCAR')
    # Preserve upstream all-periodic slab convention; compare in x/y only.
    params = dict(mul=0.8, tol=0.6, k=args.k, radius=5.0,
                  surface_coordination=12, both_surface=False)
    hfinder = SlabAdsorptionsSitesFinder(atoms, **params)
    print('Starting H surface sites', flush=True)
    hsurface = hfinder.get_surface_sites(absorbent=['H'])
    print(f'H surface done: {len(hsurface)}, elapsed {time.monotonic()-start:.1f}s', flush=True)
    hinside = [] if args.surface_only else hfinder.get_inside_sites(absorbent=['H'])
    print(f'H inside done: {len(hinside)}, elapsed {time.monotonic()-start:.1f}s', flush=True)
    cofinder = SlabAdsorptionsSitesFinder(atoms, **params)
    co_adsorbents = ['C'] if args.co_anchor == 'C' else ['C', 'O']
    cosurface = cofinder.get_surface_sites(absorbent=co_adsorbents)
    print(f'CO surface done: {len(cosurface)}, elapsed {time.monotonic()-start:.1f}s', flush=True)
    raw = {'H_surface': hsurface, 'H_inside': hinside, 'CO_surface': cosurface}
    (args.output_dir / 'raw_sites.json').write_text(json.dumps(raw, default=plain, indent=2))
    candidate = {'H': {'positions': [s['position'] for s in hsurface + hinside]},
                 'CO': {'positions': [s['position'] for s in cosurface],
                        'normals': [s['normal'] for s in cosurface]}}
    (args.output_dir / 'candidate_sites.json').write_text(json.dumps(candidate, default=plain, indent=2))
    original = json.loads((ROOT / 'GCMC/GCMC/adsorption_sites.json').read_text())
    report = {'parameters': params, 'adsorbents': {'H': ['H'], 'CO': co_adsorbents},
              'surface_only': args.surface_only,
              'source': json.loads((ROOT / 'reproduction/PH-SA_upstream/SOURCE.json').read_text()),
              'counts': {key: len(value) for key, value in raw.items()},
              'comparison': {}, 'cpu_affinity': sorted(os.sched_getaffinity(0)),
              'limitations': 'Upstream defaults, not author settings. Candidates not approved for GCMC. No DFT or MLIP evaluation.'}
    for species in ['H', 'CO']:
        a = original[species]['positions']; b = candidate[species]['positions']
        report['comparison'][species] = {'original_to_generated': distances(a, b, atoms.cell),
                                           'generated_to_original': distances(b, a, atoms.cell)}
    if cosurface:
        dots = []
        for point, normal in zip(original['CO']['positions'], original['CO']['normals']):
            _, length = find_mic(np.asarray(candidate['CO']['positions']) - point,
                                 atoms.cell, pbc=[True, True, False])
            dots.append(float(np.dot(normal, cosurface[int(length.argmin())]['normal'])))
        report['comparison']['CO']['min_matched_normal_dot'] = min(dots)
    report['elapsed_seconds'] = time.monotonic() - start
    (args.output_dir / 'comparison.json').write_text(json.dumps(report, default=plain, indent=2))
    print(json.dumps(report, default=plain))


if __name__ == '__main__':
    main()
