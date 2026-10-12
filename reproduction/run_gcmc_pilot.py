"""Run only an explicitly approved 10-proposal approximate Cu(100) pilot."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    reference = json.loads(args.reference.read_text())
    if not reference['both_converged']:
        parser.error('Gas references did not converge; pilot must not start')
    if args.output_dir.exists():
        parser.error('Pilot output directory already exists')
    root = Path(__file__).resolve().parents[1]
    environment = os.environ.copy()
    for name in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                 'NUMEXPR_NUM_THREADS', 'TF_INTRA_OP_PARALLELISM_THREADS']:
        environment[name] = '4'
    environment['TF_INTER_OP_PARALLELISM_THREADS'] = '1'
    environment['CUDA_VISIBLE_DEVICES'] = '-1'
    subprocess.run([sys.executable, str(root / 'GCMC/GCMC/gcmc.py'),
                    str(root / 'GCMC/GCMC/POSCAR'), '--dp_model', reference['model'],
                    '--sites-file', str(root / 'GCMC/GCMC/adsorption_sites.json'),
                    '--mu-co', str(reference['candidate_mu_eV']['CO']),
                    '--mu-h', str(reference['candidate_mu_eV']['H']),
                    '--steps', '10', '--temperature', '300', '--opt-steps', '20',
                    '--seed', '20261010', '--output-dir', str(args.output_dir.resolve())],
                   env=environment, check=True)


if __name__ == '__main__':
    main()
