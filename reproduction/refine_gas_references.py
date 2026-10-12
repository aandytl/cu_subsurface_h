"""Refine MLIP dimer bond lengths; no DFT or thermodynamic validation."""
import argparse
import json
import os
from pathlib import Path
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output already exists')
    source = json.loads(args.probe.read_text())
    allowed_flag = 'Nominal gas bond has force >0.5 eV/A'
    for gas in source['gas_probes'].values():
        if any(flag != allowed_flag for flag in gas['screening_flags']):
            parser.error('Source gas curve has unresolved screening flags')
    for name in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                 'NUMEXPR_NUM_THREADS', 'TF_INTRA_OP_PARALLELISM_THREADS']:
        os.environ[name] = '4'
    os.environ['TF_INTER_OP_PARALLELISM_THREADS'] = '1'
    os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
    import numpy as np
    from scipy.optimize import minimize_scalar
    from deepmd.infer import DeepPot
    start = time.monotonic()
    model = DeepPot(source['model'])
    names = model.get_type_map()
    results = {}
    for species, gas in source['gas_probes'].items():
        grid = gas['curve']
        center = gas['scanned_minimum']['bond_angstrom']
        spacing = grid[1]['bond_angstrom'] - grid[0]['bond_angstrom']
        bounds = (center - spacing, center + spacing)
        types = [names.index(s) for s in gas['symbols']]
        evaluations = []

        def evaluate(bond):
            positions = np.full((2, 3), 10.0)
            positions[0, 2] -= bond / 2
            positions[1, 2] += bond / 2
            energy, force, _ = model.eval(positions.reshape(1, 6),
                                          np.diag([20.0] * 3).reshape(1, 9), types)
            energy = float(energy.reshape(-1)[0])
            fmax = float(np.max(np.linalg.norm(force.reshape(2, 3), axis=1)))
            if not np.isfinite(energy) or not np.isfinite(fmax):
                raise ValueError('Non-finite gas prediction')
            row = {'bond_angstrom': float(bond), 'energy_eV': energy,
                   'max_force_eV_per_angstrom': fmax}
            evaluations.append(row)
            return row

        opt = minimize_scalar(lambda bond: evaluate(bond)['energy_eV'],
                              bounds=bounds, method='bounded',
                              options={'maxiter': 50, 'xatol': 1e-7})
        final = evaluate(opt.x)
        converged = bool(opt.success and final['max_force_eV_per_angstrom'] <= 0.05
                         and bounds[0] + 1e-5 < opt.x < bounds[1] - 1e-5)
        results[species] = {'final': final, 'bounds_angstrom': bounds,
                            'converged': converged, 'optimizer_message': str(opt.message),
                            'evaluations': evaluations}
        print(species, 'converged:', converged, final, flush=True)
    conditions = source['assumptions']
    mu = {'H': 0.5 * results['H2']['final']['energy_eV'] - conditions['u_rhe_V']
               + conditions['delta_h_eV_per_H'],
          'CO': results['CO']['final']['energy_eV'] + conditions['delta_co_eV_per_CO']}
    report = {'status': 'MLIP gas approximation; workflow pilot only, not author mu',
              'source_probe': str(args.probe.resolve()), 'model': source['model'],
              'assumptions': conditions, 'optimizer': {'method': 'bounded 1D bond search',
                                                       'maxiter_per_gas': 50, 'xatol_A': 1e-7,
                                                       'fmax_eV_per_A': 0.05},
              'gases': results, 'both_converged': all(r['converged'] for r in results.values()),
              'candidate_mu_eV': mu, 'elapsed_seconds': time.monotonic() - start,
              'cpu_affinity': sorted(os.sched_getaffinity(0)),
              'limitation': 'No DFT gas benchmark or adsorbate vibrational thermal correction.'}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as f:
        json.dump(report, f, indent=2)


if __name__ == '__main__':
    main()
