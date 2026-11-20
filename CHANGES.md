# Changes in this fork (branch `jhr-developments`)

Based on `RBrook318/PlasmaDissociation-python` `main` at `fbae5d8`. One entry per commit, oldest first.
All new input keys are optional and default to the original behaviour (see the README section
"Additional input options"). Old code replaced by the June and August fixes is kept as comments.

## 1. Derivative-coupling input for more than two states (June 2026) – `code/qchem.py`
Q-Chem expects the `$derivative_coupling` block to contain one line listing all states requested by
`CIS_DER_NUMSTATE` (e.g. `1 2 3`). The previous code wrote one line per state pair, which only works for
two states. The block is now written only when `num_states > 1`.

## 2. NAC parser selects valid state pairs (August 2026) – `code/qchem.py` (`readqchemnac`)
The parser now accepts coupling blocks for any state pair (i, j) within the requested states, instead of
the pair (1, `Start_State`). `Start_State` is the initially populated state and does not determine which
coupling blocks Q-Chem prints. The parser stops after reading one value per atom.
Not yet tested: full multistate (three or more states) dynamics; `prop.py` still uses two-state coupling
terms.

## 3. Configurable Q-Chem symmetry and spin settings (October 2026) – `code/global_vars.py`, `code/qchem.py`
`Unrestricted`, `Sym_ignore`, `Sym_tol` and `Use_integral_symmetry` (new, optional, under `Elec_structure`)
replace values that were hard-coded in the force inputs. Integral symmetry is also passed to the NAC
input. This allows controlled tests of how symmetry handling affects the electronic solution.

## 4. Separate symmetry settings for the setup optimisation and frequency job (October 2026) – `code/setup.py`, `code/qchem.py`
`Opt_sym_ignore`, `Opt_sym_tol` and `Opt_geom_symflag` control the setup geometry optimisation.
`sym_ignore`/`sym_tol` are now also passed to the **frequency** job. Before this change the optimisation
used `sym_ignore = true` but the frequency job used Q-Chem's default (`false`). Q-Chem could then report
the normal modes in its own reoriented frame, and they were applied to the un-reoriented geometry. This
gives wrong initial momenta (e.g. a spurious overall rotation). Runs set up with earlier versions of the
code may be affected.

## 5. NIST geometry source (October 2026) – `code/nist_geometry.py` (new), `code/setup.py`, `run/run.py`
New `Geom_flg: "NIST"`: the starting geometry is taken from the NIST CCCBDB experimental geometry for
`NIST_CAS` (neutral molecules only), then optimised as before. `run/run.py` validates the request before
creating folders or submitting jobs.

## 6. Calder support and per-job Q-Chem folder (October 2026) – `run/run.py`
`run/run.py` recognises the Leeds Calder cluster (scratch at `/scratch-calder/users/<user>`). On Calder,
each job unpacks Q-Chem into its own `mktemp` folder, which is removed when the job ends, so concurrent
array tasks cannot collide in a shared temporary folder.

## 7. Job resources (October 2026) – `run/run.py`
Generated Slurm scripts (Aire and Calder) request `--mem=4G` (was 1G) and limit arrays to 100
concurrent tasks (`%100`).

## 8. Thermal Wigner initial conditions (November 2026) – `code/setup.py`
`Molecule_data.Init_cond: "wigner"` samples positions and momenta of each normal mode from the thermal
Wigner distribution of a harmonic oscillator at `Temp` (includes zero-point motion; `Temp = 0` gives
zero-point sampling). Each repeat starts from its own displaced geometry. Optional `Seed`. Check file
`setup/wigner_checks.txt`. `"classical"` (default) keeps the original sampling.

## 9. Ground-state pre-equilibration (November 2026) – `code/main.py`
`Propagation.GS_prep` runs each new trajectory first on the ground state (`Multiplicity`) for N steps,
N drawn per repeat from [`Min_steps`, `Max_steps`]. It then switches to the target multiplicity with a
fresh SCF and resets the clock to t = 0. Logged in `output/gs_prep.out`; the prep steps are not part of
the trajectory output. Off by default.

## 10. Follow-ups to 8 and 9 (November 2026) – `code/init.py`, `code/main.py`, `code/setup.py`
- With `GS_prep` on a new trajectory, the original step-0 excited-state SCF is skipped (its result was
  discarded by the prep anyway).
- New optional `Molecule_data.Random_rotation`: an independent uniformly random rotation of each repeat
  (positions and momenta) about the centre of mass.
- If no seed is given, one is generated and logged (`setup/seed_used.txt`, `gs_prep.out`), so runs are
  reproducible. This seeds the classical sampling too. Results are statistically equivalent to before
  but not bit-identical to an unseeded run.

## 11. Documentation
README section "Additional input options" and this file.
