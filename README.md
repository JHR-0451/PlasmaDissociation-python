# PlasmaDissociation-python (fork)

Fork of [RBrook318/PlasmaDissociation-python](https://github.com/RBrook318/PlasmaDissociation-python):
Python code for simulating molecular dissociation after electron impact, with forces from Q-Chem (or PySCF).
Runs are configured with `run/inputs.json` and submitted with `run/run.py`.
This fork (branch `jhr-developments`) adds the optional input options below; see `CHANGES.md` for the full
list of changes.

## Additional input options

All options below are optional. Every new key has a default that reproduces the original behaviour,
so an existing `inputs.json` runs unchanged. See `CHANGES.md` for what changed in each commit and why.

### `Elec_structure`: Q-Chem symmetry and spin settings

| Key | Default | Meaning |
|---|---|---|
| `Sym_ignore` | `true` | Q-Chem `sym_ignore` in the dynamics (force) calculations |
| `Unrestricted` | `true` | Unrestricted SCF in the dynamics (force) calculations |
| `Sym_tol` | `5` | Q-Chem `sym_tol` in the dynamics (force) calculations |
| `Use_integral_symmetry` | `true` | Q-Chem `symmetry` (integral symmetry) in the force and NAC calculations |
| `Opt_sym_ignore` | `true` | `sym_ignore` for the setup geometry optimisation **and** frequency job |
| `Opt_sym_tol` | `5` | `sym_tol` for the setup geometry optimisation and frequency job |
| `Opt_geom_symflag` | `true` | Q-Chem `geom_opt_symflag` for the setup geometry optimisation |

Note: the NAC (derivative-coupling) input still uses `unrestricted = true` and `sym_ignore = true`
regardless of the keys above; only `Use_integral_symmetry` is passed to it.

### `Molecule_data`: starting geometry and initial conditions

| Key | Default | Meaning |
|---|---|---|
| `Geom_flg: "NIST"` | – | New geometry source: experimental geometry from the NIST CCCBDB (needs internet access at setup) |
| `NIST_CAS` | – | Required with `"NIST"`: CAS number, e.g. `"75-73-0"` (CF4); checksum is validated |
| `NIST_charge` | `0` | Only `0` (neutral) is supported |
| `Init_cond` | `"classical"` | `"classical"`: original sampling (optimised geometry for every repeat, thermal momenta). `"wigner"`: thermal harmonic Wigner sampling of positions **and** momenta per normal mode (includes zero-point motion; each repeat gets its own displaced geometry) |
| `Seed` | none | Integer seed for the initial-condition sampling. If absent, a seed is generated and logged, so every run is reproducible |
| `Random_rotation` | `false` | `true`: each repeat (positions and momenta) is given an independent uniformly random rotation about the centre of mass |

`NIST` with `Molecule` must be a plain formula (e.g. `"CF4"`); the isomer is chosen by `NIST_CAS`.
The NIST request is validated by `run/run.py` before any folders are created or jobs submitted.

### `Propagation`: ground-state pre-equilibration (`GS_prep`)

Before the excited-state dynamics, each repeat can first be propagated on the ground state for a
random number of steps, then switched to the target state.

| Key | Default | Meaning |
|---|---|---|
| `Enabled` | `false` | Turn the pre-equilibration on |
| `Multiplicity` | required if enabled | Multiplicity used during the prep (e.g. `1` for a closed-shell ground state) |
| `Min_steps`, `Max_steps` | required if enabled | Number of prep steps N is drawn uniformly from [Min_steps, Max_steps] per repeat (Min = Max gives a fixed N) |
| `Seed` | none | Seed for drawing N (combined with the repeat number). If absent, one is generated and logged |

At the switch: the multiplicity is set back to `Molecule_data.Multiplicity`, the SCF is restarted from
scratch at the new geometry, and the clock is reset to t = 0. Prep steps are not written to the trajectory
outputs. The prep runs only for a new trajectory, never on a restart.
If `Enabled` is true, a missing `Multiplicity`, `Min_steps` or `Max_steps` stops the run with a `KeyError`.

### Example (`inputs.json` excerpt)

```json
"Molecule_data": {
    "Geom_flg": "NIST", "Molecule": "CF4", "NIST_CAS": "75-73-0", "NIST_charge": 0,
    "Init_cond": "wigner", "Seed": 12345, "Random_rotation": false
},
"Elec_structure": {
    "Sym_ignore": true, "Unrestricted": true, "Sym_tol": 5, "Use_integral_symmetry": true,
    "Opt_sym_ignore": true, "Opt_sym_tol": 5, "Opt_geom_symflag": true
},
"Propagation": {
    "GS_prep": {"Enabled": true, "Multiplicity": 1, "Min_steps": 0, "Max_steps": 400, "Seed": 67890}
}
```
(Other existing keys in each block are unchanged and omitted here.)

### New output files

| File | Written by | Content |
|---|---|---|
| `setup/seed_used.txt` | setup | `Init_cond`, seed used (from inputs or generated), `Random_rotation` |
| `setup/wigner_checks.txt` | setup (`"wigner"` only) | Expected vs sampled kinetic energy, vibrational energy and RMS displacement |
| `repetitions/rep-N/output/gs_prep.out` | each repeat (`GS_prep` only) | Seed, N, energies per prep step, the switch energy and the geometry at the switch |

## Setting up and running on the Leeds HPC (Aire / Calder)

These are the steps used for our runs on Calder. The environment-creation command is reconstructed:
the environment we used has exactly the package versions in `environment.yml`. Module names on Aire
have not been checked here.

1. **Clone into your home folder at exactly `~/PlasmaDissociation-python`.** The generated job scripts
   call `/users/<user>/PlasmaDissociation-python/code/setup.py` and `code/main.py`, so a copy anywhere
   else is not used.
   ```bash
   git clone -b jhr-developments https://github.com/JHR-0451/PlasmaDissociation-python.git ~/PlasmaDissociation-python
   ```
   If that folder already holds another copy of the code, rename it first (e.g.
   `mv ~/PlasmaDissociation-python ~/PlasmaDissociation-python-old`), because only one copy can be used at a time.
2. **Create the conda environment `plasma`** from `environment.yml` (pinned versions, Python 3.12):
   ```bash
   module load calder/cpu && module load miniforge3/26.1.1-3    # Calder
   conda env create -f ~/PlasmaDissociation-python/environment.yml
   ```
   The `tools/` folder does not need to be copied into pyqchem for this code (an older README said so).
3. **Q-Chem:** a packed Q-Chem installation must be at `~/qchem.tar.gz`. Every job unpacks it to local
   temporary storage; the licence comes from the University of Leeds licence server (set in the job scripts).
4. **Run** from the `run/` folder, with the environment active. The job scripts do not activate it
   themselves: Slurm passes on the environment you submit from.
   ```bash
   module load calder/cpu && module load miniforge3/26.1.1-3 && conda activate plasma
   cd ~/PlasmaDissociation-python/run
   # edit inputs.json (Runfolder, repeats, molecule and the options above)
   python run.py
   ```
   `run.py` creates the run folder in scratch (Calder: `/scratch-calder/users/<user>/<Runfolder>`) and
   submits a setup job, then an array job for the repeats that starts once the setup job has succeeded.

**Calder specifics.** `run/run.py` recognises Calder by hostname. Each Slurm job unpacks Q-Chem into its
own temporary folder (`mktemp`, removed on exit), so array tasks cannot overwrite each other's Q-Chem
installation. On both Aire and Calder the job scripts request `--mem=4G` and limit arrays to 100
concurrent tasks (`%100`).
