"""Diagnose MLIP gas extrapolation; candidate mu values are NOT author inputs."""
import argparse
import json
import os
from pathlib import Path
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--delta-h', required=True, type=float,
                        help='Assumed gas-minus-adsorbate thermal correction, eV per H')
    parser.add_argument('--delta-co', required=True, type=float,
                        help='Assumed gas-minus-adsorbate thermal correction, eV per CO')
    parser.add_argument('--u-rhe', type=float, default=-1.0)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output already exists; choose a new file')
    for name in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                 'NUMEXPR_NUM_THREADS', 'TF_INTRA_OP_PARALLELISM_THREADS']:
        os.environ[name] = '4'
    os.environ['TF_INTER_OP_PARALLELISM_THREADS'] = '1'
    os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
    import numpy as np
    from deepmd.infer import DeepPot
    if not np.all(np.isfinite([args.delta_h, args.delta_co, args.u_rhe])):
        parser.error('Chemical-potential terms must be finite')
    start = time.monotonic()
    model = DeepPot(str(args.model.resolve(strict=True)))
    names = model.get_type_map()

    def evaluate(symbols, bond, vacuum):
        positions = np.full((2, 3), vacuum / 2)
        positions[0, 2] -= bond / 2
        positions[1, 2] += bond / 2
        box = np.diag([vacuum] * 3).reshape(1, 9)
        energy, force, _ = model.eval(positions.reshape(1, 6), box,
                                      [names.index(s) for s in symbols])
        e = float(energy.reshape(-1)[0])
        f = float(np.max(np.linalg.norm(force.reshape(2, 3), axis=1)))
        if not np.isfinite(e) or not np.isfinite(f):
            raise ValueError('Non-finite gas prediction')
        return {'bond_angstrom': float(bond), 'cell_angstrom': vacuum,
                'energy_eV': e, 'max_force_eV_per_angstrom': f}

    gases = {}
    for species, symbols, nominal, lo, hi in [
            ('H2', ['H', 'H'], 0.74, 0.50, 1.10),
            ('CO', ['C', 'O'], 1.13, 0.90, 1.40)]:
        curve = [evaluate(symbols, bond, 20.0) for bond in np.linspace(lo, hi, 13)]
        nominal_cells = [evaluate(symbols, nominal, cell) for cell in [15.0, 20.0, 25.0]]
        separated = evaluate(symbols, 10.0, 30.0)
        minimum = min(curve, key=lambda row: row['energy_eV'])
        span = max(r['energy_eV'] for r in nominal_cells) - min(r['energy_eV'] for r in nominal_cells)
        flags = []
        if minimum is curve[0] or minimum is curve[-1]:
            flags.append('No internal minimum on the scanned bond interval')
        if abs(minimum['bond_angstrom'] - nominal) > 0.15:
            flags.append('Scanned minimum differs from preset nominal bond by >0.15 A')
        if nominal_cells[1]['max_force_eV_per_angstrom'] > 0.5:
            flags.append('Nominal gas bond has force >0.5 eV/A')
        if minimum['energy_eV'] >= separated['energy_eV'] - 0.1:
            flags.append('No binding stronger than 0.1 eV relative to separated atoms')
        if span > 0.001:
            flags.append('Vacuum size changes nominal energy by >1 meV')
        gases[species] = {'symbols': symbols, 'nominal_bond_angstrom': nominal,
                          'curve': curve, 'nominal_cell_checks': nominal_cells,
                          'separated_atoms': separated, 'scanned_minimum': minimum,
                          'vacuum_energy_span_eV': span, 'screening_flags': flags}
        print(species, 'flags:', flags, flush=True)
    mu_h = 0.5 * gases['H2']['nominal_cell_checks'][1]['energy_eV'] - args.u_rhe + args.delta_h
    mu_co = gases['CO']['nominal_cell_checks'][1]['energy_eV'] + args.delta_co
    report = {'status': 'unvalidated MLIP extrapolation; never auto-start GCMC',
              'model': str(args.model.resolve()), 'type_map': names,
              'rcut_angstrom': float(model.get_rcut()), 'gas_probes': gases,
              'assumptions': {'reference_method': 'isolated molecules evaluated with supplied MLIP',
                              'u_rhe_V': args.u_rhe, 'pH': 7, 'pCO_bar': 1.0, 'p0_bar': 1.0,
                              'sampling_T_K': 300, 'delta_h_eV_per_H': args.delta_h,
                              'delta_co_eV_per_CO': args.delta_co,
                              'nominal_bond_source': 'reproduction presets, not author gas geometries'},
              'candidate_mu_eV': {'H': mu_h, 'CO': mu_co},
              'elapsed_seconds': time.monotonic() - start,
              'cpu_affinity': sorted(os.sched_getaffinity(0)),
              'limitation': 'Passing these screens does not validate gas energies against DFT. '
                            'All supplied dataset systems contain Cu; isolated gases may be out of distribution.'}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as f:
        json.dump(report, f, indent=2)
    print('Unvalidated candidate mu (eV):', report['candidate_mu_eV'], flush=True)


if __name__ == '__main__':
    main()
