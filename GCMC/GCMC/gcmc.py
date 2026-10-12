from ase.io import read, write, Trajectory
from ase import Atom, Atoms
from ase.optimize import BFGS
from deepmd.calculator import DP
import numpy as np
import random
import json
import sys
import os
import argparse
from pathlib import Path
from typing import List, Tuple, Optional, Dict
from itertools import combinations
from scipy.spatial import KDTree

kB = 8.617333262e-5

def load_adsorption_sites(json_file: str) -> Tuple[List[np.ndarray], List[np.ndarray], List[np.ndarray]]:
    with open(json_file, "r") as f:
        loaded_data = json.load(f)

    h_positions = [np.array(pos) for pos in loaded_data["H"]["positions"]]
    co_positions = [np.array(pos) for pos in loaded_data["CO"]["positions"]]
    co_normals = [np.array(normal) for normal in loaded_data["CO"]["normals"]]

    return h_positions, co_positions, co_normals

def calculate_distance(atoms: Atoms, point1: np.ndarray, point2: np.ndarray) -> float:
    vectors = atoms.cell
    min_dist = float('inf')

    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            for dz in [-1, 0, 1]:
                offset = dx * vectors[0] + dy * vectors[1] + dz * vectors[2]
                image = point1 + offset
                dist = np.linalg.norm(image - point2)
                if dist < min_dist:
                    min_dist = dist
    return min_dist

def check_distance_pbc(atoms: Atoms, new_pos: np.ndarray, existing_positions: List[np.ndarray], min_dist: float) -> bool:
    for pos in existing_positions:
        if calculate_distance(atoms, new_pos, pos) < min_dist:
            return False
    return True

def sort_atoms(atoms: Atoms) -> Atoms:
    sorted_atoms = Atoms(cell=atoms.cell, pbc=atoms.pbc)

    for atom in atoms:
        if atom.symbol == 'Cu':
            sorted_atoms.append(atom)
    for atom in atoms:
        if atom.symbol == 'C':
            sorted_atoms.append(atom)
    for atom in atoms:
        if atom.symbol == 'O':
            sorted_atoms.append(atom)
    for atom in atoms:
        if atom.symbol == 'H':
            sorted_atoms.append(atom)

    return sorted_atoms

class SiteIndex:
    def __init__(self, sites: List[np.ndarray], normals: Optional[List[np.ndarray]] = None):
        self.sites = np.array(sites)
        self.normals = np.array(normals) if normals is not None else None
        self.kdtree = KDTree(self.sites)

    def find_nearest_in_radius(self, position: np.ndarray,
                               radius: float,
                               max_neighbors: int = 10) -> Tuple[List[np.ndarray], List[np.ndarray]]:
        distances, indices = self.kdtree.query(position, k=min(max_neighbors, len(self.sites)),
                                               distance_upper_bound=radius)

        valid_sites = []
        valid_normals = []

        for i, idx in enumerate(indices):
            if idx < len(self.sites) and distances[i] <= radius:
                valid_sites.append(self.sites[idx])
                if self.normals is not None:
                    valid_normals.append(self.normals[idx])

        return valid_sites, valid_normals

def add_CO(atoms: Atoms,
           surface_sites: List[np.ndarray],
           surface_normals: List[np.ndarray],
           min_distance_C: float = 1.6,
           min_distance_O: float = 1.6) -> Optional[Atoms]:
    existing_C = [a.position.copy() for a in atoms if a.symbol == 'C']
    existing_O = [a.position.copy() for a in atoms if a.symbol == 'O']
    existing_H = [a.position.copy() for a in atoms if a.symbol == 'H']

    all_existing = existing_C + existing_O + existing_H

    shuffled_indices = random.sample(range(len(surface_sites)), min(100, len(surface_sites)))

    for idx in shuffled_indices:
        try:
            C_pos = surface_sites[idx]
            O_pos = C_pos + surface_normals[idx] * 1.13

            valid_C = check_distance_pbc(atoms, C_pos, all_existing, min_distance_C)
            valid_O = check_distance_pbc(atoms, O_pos, all_existing, min_distance_O)

            if valid_C and valid_O:
                new_atoms = atoms.copy()
                new_atoms.append(Atom('C', C_pos))
                new_atoms.append(Atom('O', O_pos))
                return new_atoms
        except IndexError:
            continue
    return None

def add_H(atoms: Atoms,
          surface_sites: List[np.ndarray],
          min_distance: float = 1.6) -> Optional[Atoms]:
    existing_pos = [a.position.copy() for a in atoms if a.symbol in ['C', 'O', 'H']]

    shuffled_sites = random.sample(surface_sites, min(100, len(surface_sites)))

    for site in shuffled_sites:
        if check_distance_pbc(atoms, site, existing_pos, min_distance):
            new_atoms = atoms.copy()
            new_atoms.append(Atom('H', site))
            return new_atoms
    return None

def remove_CO(atoms: Atoms) -> Optional[Atoms]:
    C_indices = [i for i, a in enumerate(atoms) if a.symbol == 'C']
    if not C_indices:
        return None

    co_pairs = []
    for i in C_indices:
        C_pos = atoms[i].position
        for j, a in enumerate(atoms):
            if a.symbol == 'O' and np.linalg.norm(a.position - C_pos) < 1.2:
                co_pairs.append((i, j))
                break

    if not co_pairs:
        return None

    i, j = sorted(random.choice(co_pairs), reverse=True)
    new_atoms = atoms.copy()
    del new_atoms[i]
    del new_atoms[j]
    return new_atoms

def remove_H(atoms: Atoms) -> Optional[Atoms]:
    H_indices = [i for i, a in enumerate(atoms) if a.symbol == 'H']
    if not H_indices:
        return None

    new_atoms = atoms.copy()
    del_idx = random.choice(H_indices)
    del new_atoms[del_idx]
    return new_atoms

def move_CO(atoms: Atoms,
            site_index: SiteIndex,
            search_radius: float = 3.0,
            max_candidates: int = 10,
            min_distance_C: float = 1.6,
            min_distance_O: float = 1.6) -> Atoms:

    co_pairs = []
    for i, a in enumerate(atoms):
        if a.symbol == 'C':
            for j, b in enumerate(atoms):
                if b.symbol == 'O' and np.linalg.norm(b.position - a.position) < 1.2:
                    co_pairs.append((i, j, a.position))
                    break
    if not co_pairs:
        return atoms

    c_idx, o_idx, c_pos = random.choice(co_pairs)

    candidate_sites, candidate_normals = site_index.find_nearest_in_radius(
        c_pos, search_radius, max_neighbors=max_candidates)
    if not candidate_sites:
        return atoms

    existing_pos = [a.position.copy() for i, a in enumerate(atoms)
                    if a.symbol in ['C', 'O', 'H'] and i not in [c_idx, o_idx]]

    candidate_indices = list(range(len(candidate_sites)))
    random.shuffle(candidate_indices)

    for idx in candidate_indices:
        C_pos = candidate_sites[idx]
        O_pos = C_pos + candidate_normals[idx] * 1.13

        valid_C = check_distance_pbc(atoms, C_pos, existing_pos, min_distance_C)
        valid_O = check_distance_pbc(atoms, O_pos, existing_pos, min_distance_O)

        if valid_C and valid_O:
            new_atoms = atoms.copy()
            new_atoms[c_idx].position = C_pos
            new_atoms[o_idx].position = O_pos
            return new_atoms

    return atoms

def move_H(atoms: Atoms,
           site_index: SiteIndex,
           search_radius: float = 3,
           max_candidates: int = 10,
           min_distance: float = 1.6) -> Atoms:

    H_indices = [i for i, a in enumerate(atoms) if a.symbol == 'H']
    if not H_indices:
        return atoms

    h_idx = random.choice(H_indices)
    h_pos = atoms[h_idx].position.copy()

    candidate_sites, _ = site_index.find_nearest_in_radius(
        h_pos, search_radius, max_neighbors=max_candidates)
    if not candidate_sites:
        return atoms

    existing_pos = [a.position.copy() for i, a in enumerate(atoms)
                    if a.symbol in ['C', 'O', 'H'] and i != h_idx]

    random.shuffle(candidate_sites)

    for site in candidate_sites:
        if check_distance_pbc(atoms, site, existing_pos, min_distance):
            new_atoms = atoms.copy()
            new_atoms[h_idx].position = site
            return new_atoms

    return atoms

class GCMC_Simulator:
    def __init__(self,
                 atoms: Atoms,
                 sites_file: str,
                 mu_CO: float,
                 mu_H: float,
                 temperature: float,
                 steps: int = 1000,
                 dp_model: Optional[str] = None,
                 opt_params: Optional[dict] = None,
                 max_H_coverage: float = 1.25,
                 move_search_radius: float = 3.0,
                 move_max_candidates: int = 8):

        self.original = atoms.copy()
        self.atoms = atoms.copy()
        self.h_sites, self.co_sites, self.co_normals = load_adsorption_sites(sites_file)

        self.h_site_index = SiteIndex(self.h_sites)
        self.co_site_index = SiteIndex(self.co_sites, self.co_normals)

        self.mu_CO = mu_CO
        self.mu_H = mu_H
        self.T = temperature
        self.steps = steps
        self.use_dp = dp_model is not None
        self.max_H_coverage = max_H_coverage
        self.move_search_radius = move_search_radius
        self.move_max_candidates = move_max_candidates

        if self.use_dp:
            opt_params = opt_params or {}
            self.calculator = DP(model=dp_model)
            self.opt_steps = opt_params.get('opt_steps', 100)
            self.fmax = opt_params.get('fmax', 0.05)
            self.traj_prefix = opt_params.get('trajectory_prefix', "gcmc")
            os.makedirs("trajectories", exist_ok=True)
            os.makedirs("optimized_structures", exist_ok=True)

        optimized_atoms, energy = self._calculate_energy(self.atoms, -1)
        self.atoms = optimized_atoms
        self.energy = energy
        self.main_traj = Trajectory("gcmc_main.traj", 'w')

    def _calculate_energy(self, atoms: Atoms, step: int) -> tuple[Atoms, float]:
        if not self.use_dp:
            return atoms.copy(), 0.0

        atoms = atoms.copy()
        atoms.calc = self.calculator

        traj_file = f"trajectories/{self.traj_prefix}_step{step}.traj"
        dyn = BFGS(atoms, trajectory=traj_file)
        converged = bool(dyn.run(fmax=self.fmax, steps=self.opt_steps))
        optimization = {'step': step, 'bfgs_steps': int(dyn.nsteps),
                        'converged': converged,
                        'max_force_eV_per_angstrom': float(np.max(np.linalg.norm(atoms.get_forces(), axis=1))),
                        'fmax_target': self.fmax, 'opt_steps_limit': self.opt_steps}
        with open('optimization_log.jsonl', 'a') as f:
            f.write(json.dumps(optimization) + '\n')

        if atoms.calc is None:
            atoms.calc = self.calculator

        energy = float(atoms.get_potential_energy())
        optimized_atoms = atoms.copy()
        write(f"optimized_structures/step_{step}.vasp", sort_atoms(optimized_atoms))

        return optimized_atoms, energy

    def _metropolis(self, delta_G: float) -> bool:
        if delta_G < 0:
            return True
        return bool(random.random() < np.exp(-delta_G / (kB * self.T)))

    def _get_operation(self) -> str:
        n_CO = len([a for a in self.atoms if a.symbol == 'C'])
        n_H = len([a for a in self.atoms if a.symbol == 'H'])
        current_H_coverage = n_H / 144

        if current_H_coverage >= self.max_H_coverage:
            possible_ops = ['add_CO', 'remove_CO', 'move_CO', 'remove_H', 'move_H']
        else:
            if n_CO == 0 and n_H == 0:
                possible_ops = ['add_CO', 'add_H']
            elif n_CO > 0 and n_H == 0:
                possible_ops = ['add_CO', 'remove_CO', 'move_CO', 'add_H']
            elif n_CO == 0 and n_H > 0:
                possible_ops = ['add_H', 'remove_H', 'move_H', 'add_CO']
            else:
                possible_ops = ['add_CO', 'remove_CO', 'move_CO', 'add_H', 'remove_H', 'move_H']

        return random.choice(possible_ops)

    def run_step(self, step: int) -> Dict:
        op = self._get_operation()

        old_atoms = self.atoms.copy()
        old_energy = self.energy

        if op == 'add_CO':
            new_atoms = add_CO(old_atoms, self.co_sites, self.co_normals)
        elif op == 'remove_CO':
            new_atoms = remove_CO(old_atoms)
        elif op == 'add_H':
            new_atoms = add_H(old_atoms, self.h_sites)
        elif op == 'remove_H':
            new_atoms = remove_H(old_atoms)
        elif op == 'move_CO':
            new_atoms = move_CO(old_atoms, self.co_site_index,
                                self.move_search_radius, self.move_max_candidates)
        elif op == 'move_H':
            new_atoms = move_H(old_atoms, self.h_site_index,
                               self.move_search_radius, self.move_max_candidates)
        else:
            new_atoms = None

        if not new_atoms:
            return {'step': step, 'accepted': False, 'operation': op,
                    'proposal_generated': False}

        new_atoms, new_energy = self._calculate_energy(new_atoms, step)

        delta_n_CO = len([a for a in new_atoms if a.symbol == 'C']) - len([a for a in old_atoms if a.symbol == 'C'])
        delta_n_H = len([a for a in new_atoms if a.symbol == 'H']) - len([a for a in old_atoms if a.symbol == 'H'])

        delta_G = float((new_energy - old_energy) - (self.mu_CO * delta_n_CO + self.mu_H * delta_n_H))
        accepted = self._metropolis(delta_G)

        if accepted:
            self.atoms = new_atoms
            self.energy = new_energy
            self.main_traj.write(self.atoms)

        return {
            'step': step,
            'operation': op,
            'accepted': accepted,
            'delta_G': delta_G,
            'energy': new_energy
        }

    def run(self) -> List[Dict]:
        results = []
        for step in range(self.steps):
            result = self.run_step(step)
            result.update({'current_energy_eV': self.energy,
                           'n_CO': sum(a.symbol == 'C' for a in self.atoms),
                           'n_H': sum(a.symbol == 'H' for a in self.atoms)})
            with open('mc_steps.jsonl', 'a') as f:
                f.write(json.dumps(result) + '\n')
            results.append(result)
            write("CONTCAR", sort_atoms(self.atoms))
            if step % 10 == 0:
                print(f"Step {step}/{self.steps}: {result}")
        return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the supplied Cu(100) GCMC implementation. Chemical potentials must be independently verified (eV).")
    parser.add_argument("input_file", type=Path)
    parser.add_argument("--dp_model", required=True, type=Path)
    parser.add_argument("--sites-file", type=Path, default=Path(__file__).resolve().parent / "adsorption_sites.json")
    parser.add_argument("--mu-co", required=True, type=float)
    parser.add_argument("--mu-h", required=True, type=float)
    parser.add_argument("--steps", type=int, default=15000)
    parser.add_argument("--temperature", type=float, default=300)
    parser.add_argument("--opt-steps", type=int, default=400)
    parser.add_argument("--seed", type=int, default=20261010)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    if args.steps <= 0 or args.opt_steps <= 0 or args.temperature <= 0:
        parser.error("steps, opt-steps and temperature must be positive")
    if not all(np.isfinite(x) for x in [args.mu_co, args.mu_h, args.temperature]):
        parser.error("chemical potentials and temperature must be finite")
    input_file, dp_model, sites_file = [str(p.resolve(strict=True)) for p in
                                        [args.input_file, args.dp_model, args.sites_file]]
    args.output_dir.mkdir(parents=True, exist_ok=False)
    random.seed(args.seed)
    np.random.seed(args.seed)
    base_atoms = read(input_file)
    os.chdir(args.output_dir)
    with open("run_config.json", "w") as f:
        json.dump({**vars(args), "input_file": input_file, "dp_model": dp_model,
                   "sites_file": sites_file}, f, default=str, indent=2)

    simulator = GCMC_Simulator(
        atoms=base_atoms,
        sites_file=sites_file,
        mu_CO=args.mu_co,
        mu_H=args.mu_h,
        temperature=args.temperature,
        steps=args.steps,
        dp_model=dp_model,
        max_H_coverage=1,
        move_search_radius=3.0,
        move_max_candidates=100,
        opt_params={'opt_steps': args.opt_steps, 'fmax': 0.05}
    )

    results = simulator.run()
    simulator.main_traj.close()
    with open("results.json", "w") as f:
        json.dump(results, f, indent=2)
    write("CONTCAR", sort_atoms(simulator.atoms))

    accepted = [r for r in results if r['accepted']]
    print(f"Overall Acceptance Rate: {len(accepted)/len(results):.1%}")
