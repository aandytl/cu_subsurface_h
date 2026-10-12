"""Evaluate supplied labels after explicit user approval; never train a model."""
import argparse
import json
import os
from pathlib import Path
import time
import platform
from importlib.metadata import version


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--model', type=Path, required=True)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--frames', type=int, default=100)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.frames < 1 or args.threads < 1:
        parser.error('frames and threads must be positive')
    for name in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                 'NUMEXPR_NUM_THREADS', 'TF_INTRA_OP_PARALLELISM_THREADS']:
        os.environ[name] = str(args.threads)
    os.environ['TF_INTER_OP_PARALLELISM_THREADS'] = '1'
    os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
    import numpy as np
    from deepmd.infer import DeepPot
    start = time.monotonic()
    print('Loading frozen model on CPU...', flush=True)
    model = DeepPot(str(args.model.resolve(strict=True)))
    model_types = model.get_type_map()
    load_seconds = time.monotonic() - start
    index_start = time.monotonic()
    systems = []
    for coord in sorted(args.data.glob('*/set.*/coord.npy')):
        nframes = np.load(coord, mmap_mode='r').shape[0]
        systems.append((coord.parent, nframes))
    total = sum(n for _, n in systems)
    if total == 0:
        raise ValueError('No dataset frames found')
    # Evenly spaced frames across the complete ordered dataset, without randomness.
    selected = set(np.linspace(0, total - 1, min(args.frames, total), dtype=int).tolist())
    index_seconds = time.monotonic() - index_start
    print('Selected {} of {} frames'.format(len(selected), total), flush=True)
    evaluation_start = time.monotonic()
    rows = []
    dtypes = set()
    energy_sq = force_sq = 0.0
    force_count = offset = 0
    for directory, nframes in systems:
        indices = [i for i in range(nframes) if offset + i in selected]
        offset += nframes
        if not indices:
            continue
        names = (directory.parent / 'type_map.raw').read_text().split()
        raw_types = np.loadtxt(directory.parent / 'type.raw', dtype=int).reshape(-1)
        types = np.array([model_types.index(names[t]) for t in raw_types])
        arrays = {key: np.load(directory / (key + '.npy'), mmap_mode='r')
                  for key in ['coord', 'box', 'energy', 'force']}
        for i in indices:
            energy, force, _ = model.eval(arrays['coord'][i:i+1], arrays['box'][i:i+1], types)
            dtypes.add((str(arrays['coord'].dtype), str(energy.dtype), str(force.dtype)))
            de = (float(energy.reshape(-1)[0]) - float(arrays['energy'][i])) / len(types)
            df = force.reshape(-1) - arrays['force'][i].reshape(-1)
            if not np.isfinite(de) or not np.all(np.isfinite(df)):
                raise ValueError('Non-finite model prediction or label')
            energy_sq += de * de
            force_sq += float(np.dot(df, df))
            force_count += df.size
            rows.append({'system': str(directory), 'frame': i, 'atoms': len(types),
                         'energy_error_eV_per_atom': de,
                         'force_component_rmse_meV_per_angstrom':
                         1000 * float(np.sqrt(np.dot(df, df) / df.size))})
            if len(rows) % 10 == 0:
                print('Evaluated {}/{} frames'.format(len(rows), len(selected)), flush=True)
    report = {'model': str(args.model.resolve()), 'data': str(args.data.resolve()),
              'type_map': model_types, 'available_frames': total, 'evaluated_frames': len(rows),
              'energy_rmse_meV_per_atom': 1000 * (energy_sq / len(rows)) ** 0.5,
              'force_component_rmse_meV_per_angstrom': 1000 * (force_sq / force_count) ** 0.5,
              'elapsed_seconds': time.monotonic() - start, 'threads': args.threads,
              'timing_seconds': {'model_load': load_seconds, 'dataset_index': index_seconds,
                                 'evaluation': time.monotonic() - evaluation_start},
              'python': platform.python_version(),
              'cpu_affinity': sorted(os.sched_getaffinity(0)),
              'packages': {name: version(name) for name in
                           ['numpy', 'scipy', 'ase', 'deepmd-kit', 'tensorflow-cpu']},
              'dtypes_coord_energy_force': sorted(dtypes),
              'sampling': 'evenly spaced in lexicographically ordered dataset', 'frames': rows}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as f:
        json.dump(report, f, indent=2)
    print(json.dumps({k: v for k, v in report.items() if k != 'frames'}, indent=2))


if __name__ == '__main__':
    main()
