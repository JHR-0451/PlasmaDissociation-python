"""
main.py

This script serves as the main entry point for running molecular simulations using the 
custom-built framework. It orchestrates the entire simulation process, from initializing 
the molecular structure to processing results. 

Key Functionalities:
--------------------
- Loads input parameters from a JSON configuration file.
- Initializes the molecular structure based on specified simulation parameters.
- Runs electronic structure calculations and properties evaluations.
- Handles restart functionality for previously interrupted simulations.
- Outputs molecular data at each timestep for analysis.
- Measures and records the execution time of the simulation.
- Processes and collates the results of the simulation at the end.

Input:
------
- A JSON file (`../inputs.json`) that contains configuration settings for the simulation.
  This file must specify the following parameters:
  - `repeats`: Number of repetitions for the simulation.
  - `cores`: Number of CPU cores to use for calculations.
  - `States`: Number of states to consider in the simulation.
  - `Branches`: Whether to allow for cloning ().
  - `Multiplicity`: Initial spin multiplicity for the molecule.
  - `Timestep`: Time increment for the simulation steps.
  - `Tot_timesteps`: Total number of timesteps to run.
  - `Geom_start`: If reading from a folder of geometries, which to start from.
  - `Spin_flip`: Boolean to indicate if spin-flip is to be used.
  - `method`: The electronic structure package to use.

Output:
-------
- Outputs molecular data to the `output/` directory, including:
  - `xyz.all`: A file containing the coordinates trajectory data.
  - `moementa.all`: A file containing the momenta trajectory data.
  - `forces.all`: A file containing the forces, multiplicity, dissociation flag of the trajectory.
  - `time.out`: A file recording the execution time of the simulation.
  - JSON files representing the molecular structure at each timestep.

Usage:
------
Run the script from the command line. Ensure that the input JSON file is correctly formatted 
and accessible at the specified path. The script will handle initialization, execution, 
and result processing, with proper checks for whether a simulation can be restarted.

Exceptions:
-----------
- The script checks if a simulation has already completed by examining the `xyz.all` file.
  If the run is complete, it will exit and not process results to avoid unnecessary computation.
- It raises a `ValueError` if any expected files or directories are missing.

Author: Ryan Brook
Date: 30/10/2024
"""

import os
import json
import time
import numpy as np
from init import Molecule, initialize_structure, create_empty_molecule, initialize_simulation
import elec
import prop
import output as out
import result
import global_vars as gv

# def initialize_simulation(restart):
#     """Initialize the molecular structure and handle restart logic."""
   
#     time1 =time.time()
#     if restart == 'NO':
#         time2= time.time()
#         molecule1.time[1] = time2-time1
#         molecule1.time[3] += time2-time1
#         molecule1.time[4] = molecule1.timestep/gv.timestep
        
#         molecule1.time[0] = 0
#         return molecule1, molecule2, 1, True
    
#     # If JSON file doesn't exist, initialize a new structure
#     molecule1 = initialize_structure()
#     n = len(molecule1.symbols)
#     molecule2 = create_empty_molecule(n)
#     molecule1 = elec.run_elec_structure(molecule1, Guess=False)
#     time2= time.time()
#     molecule1.time[1] = time2-time1
#     molecule1.time[3] += time2-time1
#     molecule1.time[4] = molecule1.timestep/gv.timestep
#     out.output_molecule(molecule1)
#     molecule1.time[0] = 0
#     return molecule1, molecule2, 1, True

def run_simulation(startstep, molecule1, molecule2):
    """Run the main simulation loop for the specified number of timesteps."""
    n = len(molecule1.symbols)
    guess = True

    for i in range(int(startstep)-1, gv.tot_timesteps + 1):
        time1 = time.time()
        old_coordinates = molecule1.coordinates.copy()

        molecule2 = prop.prop_1(molecule1, molecule2)
        
        molecule2 = elec.run_elec_structure(molecule2, Guess=guess)

        molecule1.time[0] = molecule2.time[0]
        molecule1.time[2] = molecule2.time[2]

        molecule1.elecinfo = molecule2.elecinfo

        molecule1 = prop.prop_2(molecule1, molecule2)

        if gv.remove_atoms == 1:
            molecule1, dissociated = prop.fragments(molecule1)
            molecule1 = prop.prop_diss(molecule1)
            guess = dissociated == 0  # Update guess based on dissociation

        if gv.checks > 0:
            out.run_checks(molecule1,old_coordinates)
        time2= time.time()
        molecule1.time[1] = time2-time1
        molecule1.time[3] += time2-time1


        out.output_molecule(molecule1)

        molecule1.time[0] = 0

def ground_state_prep(molecule1, molecule2):
    """Optional ground-state propagation before the excited-state dynamics (pre-equilibration).

    Controlled by inputs.json -> "Propagation" -> "GS_prep" (absent or "Enabled": false = original behaviour):
        "Enabled": true, "Multiplicity": 1, "Min_steps": int, "Max_steps": int, "Seed": int (optional)
    Each trajectory propagates on the ground state (given multiplicity) for N steps, N drawn uniformly from
    [Min_steps, Max_steps] (seeded per rep from Seed and the rep number, so a restart redraws the same N).
    Propagation uses the code's own prop_1 / elec / prop_2 steps; no fragment detection, atom removal or
    checks run in this phase, and nothing is written to the trajectory outputs. At the switch, the multiplicity
    is set back, the SCF is restarted from scratch (Guess=False, as at the original step 0) at the new geometry,
    and the clock is reset, so the excited-state dynamics start at t = 0 from the pre-equilibrated geometry and
    momenta. The prep phase is logged to output/gs_prep.out.
    """
    with open("../../inputs.json", "r") as f:
        cfg = json.load(f)["Propagation"].get("GS_prep", {})
    if not cfg.get("Enabled", False):
        return molecule1, molecule2
    gs_mult = int(cfg["Multiplicity"])
    nmin, nmax = int(cfg["Min_steps"]), int(cfg["Max_steps"])
    if not 0 <= nmin <= nmax:
        raise ValueError(f"GS_prep: need 0 <= Min_steps <= Max_steps, got {nmin}, {nmax}")
    rep = os.path.basename(os.getcwd())                  # e.g. 'rep-17'
    rep_num = int(rep.split("-")[-1]) if rep.split("-")[-1].isdigit() else 0
    seed = cfg.get("Seed")
    if seed is None:                                     # O3: draw a seed from OS entropy and log it (gs_prep.out)
        seed = int(np.random.SeedSequence().entropy % 2**63)
    rng = np.random.default_rng([int(seed), rep_num])
    nsteps = int(rng.integers(nmin, nmax + 1))
    target_mult = molecule1.multiplicity
    t_start = molecule1.timestep
    AU_FS = 0.02418884

    def ke(mol):
        return float(np.sum(np.sum(mol.momenta**2, axis=1) / (2 * mol.masses)))

    with open("output/gs_prep.out", "w") as log:
        log.write(f"Ground-state prep: {rep}, multiplicity {gs_mult} -> {target_mult}, seed {seed}, "
                  f"steps {nsteps} (drawn from [{nmin}, {nmax}]), timestep {gv.timestep} a.u.\n")
        log.write("step  time_fs  E_scf(Ha)  KE(Ha)  E_total(Ha)\n")
        molecule1.multiplicity = gs_mult
        molecule1 = elec.run_elec_structure(molecule1, Guess=False)
        log.write(f"{0} {0.0:.3f} {molecule1.scf_energy[0]:.10f} {ke(molecule1):.10f} "
                  f"{molecule1.scf_energy[0] + ke(molecule1):.10f}\n")
        for k in range(1, nsteps + 1):
            molecule2 = prop.prop_1(molecule1, molecule2)
            molecule2 = elec.run_elec_structure(molecule2, Guess=True)
            molecule1.elecinfo = molecule2.elecinfo
            molecule1 = prop.prop_2(molecule1, molecule2)
            log.write(f"{k} {(molecule1.timestep - t_start) * AU_FS:.3f} {molecule1.scf_energy[0]:.10f} "
                      f"{ke(molecule1):.10f} {molecule1.scf_energy[0] + ke(molecule1):.10f}\n")
            log.flush()
        # switch to the excited-state dynamics
        molecule1.multiplicity = target_mult
        molecule1.timestep = t_start
        molecule1 = elec.run_elec_structure(molecule1, Guess=False)
        log.write(f"switch: multiplicity {target_mult}, new SCF E = {molecule1.scf_energy[0]:.10f} Ha, "
                  f"KE = {ke(molecule1):.10f} Ha; clock reset to {t_start}\n")
        log.write("geometry at switch (bohr):\n")
        for s, c in zip(molecule1.symbols, molecule1.coordinates):
            log.write(f"{s} {float(c[0])!r} {float(c[1])!r} {float(c[2])!r}\n")
    return molecule1, molecule2


def main():


    gv.load_global_variables()
    # Check basic arguments
    fresh_start = not os.path.exists('output/xyz.all')   # GS prep only for a new trajectory, never on restart

    with open("../../inputs.json", "r") as f:
        gs_enabled = json.load(f)["Propagation"].get("GS_prep", {}).get("Enabled", False)
    # O1: with GS prep on a fresh start, skip the step-0 excited-state SCF (the prep runs its own SCF first and
    # the switch runs a fresh one), so no discarded "stray" SCF is written to f.all.
    molecule1, molecule2, startstep = initialize_simulation(run_scf=not (fresh_start and gs_enabled))

    if fresh_start:
        molecule1, molecule2 = ground_state_prep(molecule1, molecule2)

    out.output_molecule(molecule1)
    
    run_simulation(startstep,molecule1, molecule2)

    result.process_results()

if __name__ == "__main__":
    main()