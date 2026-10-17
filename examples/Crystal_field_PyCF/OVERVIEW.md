# Crystal-field Hamiltonians from PyCF as DBQA benchmarks: overview

This folder adds three rare-earth crystal-field Hamiltonians, built with PyCF ([mikereidnz/pycf](https://github.com/mikereidnz/pycf)), to the qibo-dbqa pipeline and benchmarks double-bracket quantum algorithms (DBQA) on them. This document is the entry point. For each model it says what question was asked, where the Hamiltonian comes from, how it was validated, what came out, and where the work departs from the original plan. Every number here is taken from an executed notebook and is repeated, with more detail, in `RESULTS_Model_1.md`, `RESULTS_Model_2.md` and `RESULTS_Model_3.md`.

## In short

| | Model 1 | Model 2 | Model 3 |
|---|---|---|---|
| System | Ce³⁺:YAG | Er³⁺:LiYF₄ | ¹⁶⁷Er³⁺:Y₂SiO₅, site 2, in a field |
| States → qubits | 4f¹: 14 → 4 | full 4f¹¹: 364 → 9 | ⁴I₁₅/₂ + ⁴I₁₃/₂ × I = 7/2: 240 → 8 |
| Question | how fast DBQA recovers a small crystal-field spectrum; circuit cost | effect of the multiscale structure; compile cost | resolving hyperfine levels inside the full Hamiltonian, across a field sweep, well enough to locate a ZEFOZ point |
| Validation | levels and eigenspaces vs PyCF to 10⁻¹¹ cm⁻¹ | σ = 12.9 cm⁻¹ vs 108 measured levels (paper 14.6) after a free-ion refit | 7 T hyperfine splittings of the source paper reproduced to 0.7 MHz |
| Headline | all 7 levels to < 10⁻³ cm⁻¹ in 15 group-commutator iterations (with durations up to 10²); DB-DOI 1 − F < 10⁻² after 2 steps | two-stage DBQA: ⁴I₁₅/₂ levels to the 0.072 cm⁻¹ restriction floor in 16 + 12 iterations | two-stage DBQA: all 120 Z1 transitions to 2×10⁻⁴ MHz (10⁻¹² of the spectral width); ZEFOZ field to < 0.01 mT; the schedule transfers across the sweep |
| CNOTs per 1st-order Trotter step | 272 | 1.07×10⁶ | 5.6×10⁴ |

**Three lessons carry across the models.**
1. **Duration grids must reach about 1/gap** in the units the Hamiltonian is rescaled to. This was found in Model 3 and changed conclusions of Models 1 and 2 (`Example_11`, `Example_12`): with long enough durations, greedy group-commutator DBQA resolves the closest pair of Model 1 and the crystal-field levels of Model 2, which the original grids (s ≤ 1 or 10) could not.
2. **Isolate, then resolve.** Running DBQA in two stages is much faster than one, and the isolation threshold of the first stage sets the final accuracy.
3. **The group commutator with D = Δ(Ĥ_k) is the most reliable generator** once the duration grid is long enough, in every model and stage. Pauli-Z generators are competitive in small blocks (Model 2); the canonical bracket does not exploit long durations and stalls at a fixed point when the block diagonal is degenerate (Model 3).


## 1. What is where

| Item | Location |
|---|---|
| Interface (PyCF → Qibo, metrics, schedulers, compilation, Model 3 field models and ZEFOZ search) | `src/qibo_dbqa/utils_crystal_field.py` |
| Tests | `tests/test_models_crystal_field.py` |
| Notebooks | `examples/Crystal_field_PyCF/Example_1` … `Example_15` |
| Transcribed tables from the source papers | `examples/Crystal_field_PyCF/data/` |
| Figures | `examples/Crystal_field_PyCF/pdfs_figures/` |
| Pickled results | `examples/Crystal_field_PyCF/simulation_results/` (two large `.npz` caches of Model 2 are not under version control) |

The code follows the pipeline's conventions: modules in `src/qibo_dbqa/`, one notebook that generates data and another that loads it (as for Fig. 4), and pytest tests that run on the numpy and qibojit backends.

| Notebook | Model | Content | Run time |
|---|---|---|---|
| `Example_1_validation_PyCF_to_Qibo_Ce_YLF` | – | PyCF → Qibo conversion checked on PyCF's Ce:LiYF₄ example | seconds |
| `Example_2_DBQA_Ce_YAG` | 1 | matrix-level DBQA | 1.5 min |
| `Example_3_compiled_DBQA_Ce_YAG` | 1 | compiled oracle, DB-DOI | 7 min |
| `Example_4_validation_Er_LiYF4` | 2 | validation, free-ion refit | 20 min |
| `Example_5_DBQA_Er_LiYF4_generate_data` | 2 | DBQA data: full Hamiltonian, then the ground-multiplet block | 35 min |
| `Example_6_DBQA_Er_LiYF4_results` | 2 | analysis, compilation cost | < 1 min |
| `Example_7_compiled_DBQA_Er_LiYF4` | 2 | compiled oracle, symmetry-resolved DB-DOI | 1 h |
| `Example_8_validation_Er_YSO` | 3 | validation, ZEFOZ points, choice of the field sweep | 1 h |
| `Example_9_DBQA_Er_YSO_generate_data` | 3 | two-stage DBQA over a field sweep, cold/warm/replayed starts | 85 min |
| `Example_10_DBQA_Er_YSO_results` | 3 | analysis, ZEFOZ point from DBQA, compilation cost | seconds |
| `Example_11_step_range_Er_LiYF4` | 2 | stage 1 of Model 2 with durations up to 10⁴ (the correction below) | 47 min |
| `Example_12_duration_grids_Ce_YAG_Er_LiYF4` | 1, 2 | Model 1, the Model 2 block and DB-DOI with extended duration grids | 35 min |
| `Example_13_avoided_crossings_Er_YSO` | 3 | avoided crossings (both Hamiltonians), DBQA across two of them, symmetry sectors | 4 min |
| `Example_14_compiled_DBQA_Er_YSO` | 3 (and 2) | compiled oracle, full DB-DOI, two-stage DB-DOI in the 4-qubit blocks | 70 min |
| `Example_15_noisy_DB_DOI` | 1, 2, 3 | depolarizing noise on every CNOT of the DB-DOI circuits; error-rate budgets | 11 min |

## 2. Set-up

PyCF builds a C extension against GSL, NLopt and LAPACKE and does not build on native Windows; everything here ran under WSL2 (Ubuntu) with a conda-forge environment:

```bash
micromamba create -n dbqa-cf -c conda-forge python=3.12 numpy=1.26.4 "scipy>=1.15" gsl nlopt "libblas=*=*openblas" liblapacke openblas cython c-compiler make pkg-config matplotlib jupyterlab
micromamba activate dbqa-cf
CFL_CFLAGS="-I$CONDA_PREFIX/include" CFL_LDLIBS="-L$CONDA_PREFIX/lib -Wl,-rpath,$CONDA_PREFIX/lib -llapacke -lopenblas" LAPACKE_INCLUDE_DIR=$CONDA_PREFIX/include pip install --no-build-isolation /path/to/pycf
pip install qibo==0.2.14 optuna==4.1.0 hyperopt qibojit==0.1.7
pip install --no-deps -e .   # from the root of this repository
```

- **Matrix-element bases.** PyCF does not install its bases (`f1cf`, `f11cf`, `erhfs`, under `examples/*/matel`). Set `PYCF_ROOT` to a PyCF source checkout; the notebooks default to `../../../pycf`, i.e. a checkout next to `qibo-dbqa`. Tests that need PyCF are skipped without it.
- **BLAS threads.** Set `OPENBLAS_NUM_THREADS=1` (the Model 3 notebooks do this before importing numpy). With OpenBLAS's default of one thread per core, diagonalising a 240×240 matrix took 3 s instead of 5 ms on a 32-core machine, and forked worker processes inherit busy-waiting BLAS threads; the matrices here are too small to profit from more threads.
- **Tests.** `pytest -o addopts="" tests/` (pytest-cov is not installed in this environment). The suite has 159 tests and passes on the numpy and qibojit backends in about 2 minutes.
- **Notebooks.** `jupyter nbconvert --to notebook --execute --inplace <notebook>` from `examples/Crystal_field_PyCF`. Every run is deterministic: grid-search scheduling and the choice of generators involve no random sampling.

## 3. Common methods

- **PyCF → Qibo.** `pycf_dense_matrix` assembles $\sum_t c_t T_t$ exactly as PyCF does (complex $c_t$ multiplies $C^k_q$, its conjugate the partner $C^k_{-q}$); it agrees with PyCF's own diagonalisation to $10^{-10}$ cm⁻¹. `embed_in_qubits` pads the $d$ physical states to $2^{\lceil\log_2 d\rceil}$ with uncoupled states above the spectrum. Energies are divided by a model-dependent unit so the Hamiltonian norm is O(1).
- **Metrics.** Levels are grouped with a tolerance fixed from the exact spectrum before any run (`group_levels`, `check_grouping_tolerance`). For the DBQA iterate $\hat H_k$ we report energy errors (diagonal entries, ranked), projector fidelities of each level's eigenspace, gap errors, and the off-diagonal norm (`dbqa_spectrum_report`).
- **Methods.** Canonical (GWW) bracket, single commutator with Pauli-Z diagonal operators, and the group commutator with $D = \Delta(\hat H_k)$. Durations are chosen greedily by grid search on $\|\sigma(\hat H_k)\|$; `spectral_grid_search_step` gives the same losses as `grid_search_step` but diagonalises once per step, which makes 256–512-dimensional runs affordable.
- **Two-stage strategy.** Run DBQA on the full Hamiltonian until the basis states with the lowest diagonal entries span a target block (a multiplet, or a crystal-field level with its hyperfine states) to $1-F<\epsilon$; then run DBQA on that block, rescaled. The isolation $\epsilon$ bounds the final accuracy ("restriction error").
- **Compilation.** `pauli_decomposition` (fast tensor contraction), `PauliCompiledHamiltonian` (CNOT-ladder Trotter–Suzuki circuits), `compiled_trotter_unitary` and `compiled_trotter_depth` (exact unitary and depth of the compiled sequence without building $10^6$-gate circuits), and the ideal oracle `ExactEvolutionHamiltonian`.

## 4. Model 1: Ce³⁺:YAG (4f¹, 4 qubits)

**Question.** How quickly do the DBQA variants of the pipeline recover the low-energy spectrum of a small but realistic crystal-field Hamiltonian, and what does a circuit-level run (DB-DOI) cost?

**Hamiltonian.** Ce³⁺ at the D₂ site of Y₃Al₅O₁₂; 14 |²F_J M_J⟩ states (PyCF `f1cf`) in 16 dimensions. Spin–orbit ζ and nine real B^k_q from Table 10 of L. Hu, M. F. Reid, C. K. Duan, S. Xia and M. Yin, J. Phys.: Condens. Matter 23, 045501 (2011): column [a] ("fitted", the fit of Tanner et al. rotated into Hu's axes) for the main runs, column [c] ("calculated", DV-Xα) as a sensitivity check.

**Validation.** On PyCF's own Ce:LiYF₄ example (`Example_1`): the dense matrix matches PyCF's eigenvalues and eigenvectors to 10⁻¹¹ cm⁻¹; the embedding keeps every Kramers eigenspace (projector fidelity 1 to 12 digits); the padding is uncoupled; the tensors are Wybourne-normalised; the levels match PyCF's example data.

**Results** (`RESULTS_Model_1.md`).
- Matrix level, original grids (durations ≤ 0.3 or 1): the ground multiplet ²F₅/₂ is within 1 cm⁻¹ after 25 group-commutator iterations (32 with Pauli-Z generators); the canonical bracket needs more than 60. The closest pair of levels (94 cm⁻¹) is isolated (pair infidelity 6×10⁻⁵) but not resolved within 60 iterations.
- **Corrected by `Example_12`:** with durations up to 10² the group commutator brings all seven levels, the pair included, within 1 cm⁻¹ after 15 iterations and to < 10⁻³ cm⁻¹ by 60; ²F₅/₂ is within 1 cm⁻¹ after 10. The pair was limited by the grid (it needs s ≈ 10), not by its gap.
- With the other parameter set every method needs 1.4–2.3× more iterations: counts depend on the parameters even at 4 qubits.
- Circuit level: a second-order, 4-step Trotterised oracle costs 2152 CNOTs per call (error < 5×10⁻³). DB-DOI from the best basis state reaches 1 − F < 10⁻² for the ground doublet after 6 steps, 364 oracle calls, about 7.8×10⁵ CNOTs. Oracle calls triple per DB-DOI step. With durations up to 10² (`Example_12`) it reaches 1 − F < 10⁻² after 2 steps (4 oracle calls); each call then needs 32 Trotter steps, and the total falls to 6.9×10⁴ CNOTs, 11 times fewer. Under depolarizing noise after every CNOT (`Example_15`) the best ground-doublet fidelity is 0.96, 0.84 and 0.31 at p = 10⁻⁵, 10⁻⁴ and 10⁻³.

**Departure from the plan.** Hu et al.'s column [c] does not reproduce their own calculated levels (Table 9[c], differences up to 760 cm⁻¹), so the fitted column [a] is the main parameter set and [c] serves only as a sensitivity test.

## 5. Model 2: Er³⁺:LiYF₄ (full 4f¹¹, 9 qubits)

**Question.** How does the multiscale structure of a full 4fᴺ configuration (Coulomb 10⁵, spin–orbit 10³, crystal field 10² cm⁻¹) affect DBQA, and how costly is such a dense Hamiltonian to compile?

**Hamiltonian.** Er³⁺ at the S₄ site of LiYF₄ in the D₂d approximation; all 364 states of 4f¹¹ (PyCF `f11cf`) in 512 dimensions. Crystal field from M. A. Couto dos Santos et al., J. Alloys Compd. 275–277, 435 (1998), Table 2; their Table 1 is transcribed in `data/`.

**Validation and the refit.** The paper's parameters taken literally give σ = 547 cm⁻¹ in PyCF: its Racah E³ follows a convention we could not reproduce (a refit of E³ alone repairs it, the crystal field cannot). Model 2 therefore keeps the paper's crystal field and refits the 15 free-ion parameters to the paper's measured levels: σ = 12.9 cm⁻¹, against the paper's own 14.6 cm⁻¹. The paper's Table 1 has **108** measured levels; its abstract says 109.

**Results** (`RESULTS_Model_2.md`).
- Greedy group-commutator DBQA isolates the ground multiplet ⁴I₁₅/₂ to 1 − F < 10⁻⁶ after 16 iterations.
- **Two-stage strategy.** Restricting a second DBQA stage to that 4-qubit block (Pauli-Z generators) brings every ⁴I₁₅/₂ level within 0.075 cm⁻¹ of the exact value after 16 + 45 iterations. Switching at 1 − F < 10⁻³ instead leaves 7–46 cm⁻¹ of restriction error: the isolation threshold, not stage 2, sets the final accuracy.
- The slowest mode is a coupling between two levels of the same irrep (Γ₇, 0 ↔ 5.98 cm⁻¹), not the closest pair (0.36 cm⁻¹).
- Compilation: 93 184 of 262 144 Pauli strings, about 1.1×10⁶ CNOTs per first-order Trotter step (2.1×10⁶ at second order), against 36 for a 9-qubit XXZ step; an oracle accurate to 1.4×10⁻³ costs 3.4×10⁷ CNOTs per call.
- DB-DOI trap: the D₂d Hamiltonian splits into four sectors of 91 states that a Z-field D never leaves; the reference state must be chosen per sector. From the right sector, five steps (4.2×10⁹ CNOTs) take the energy from 907 to 98 cm⁻¹ but leave only 5 % on the ground doublet.

**Correction made after Model 3.** `RESULTS_Model_2.md` first concluded that greedy DBQA on the full Hamiltonian cannot resolve crystal-field levels (55–76 cm⁻¹ errors after 150 iterations) because the large free-ion couplings dominate the objective. Model 3 showed the actual cause: the duration grid stopped at s = 10, while resolving a gap g (in the rescaled units) needs durations of about 1/g, here up to 10³–10⁴. `Example_11` reruns that stage with s ≤ 10⁴: every ⁴I₁₅/₂ level then comes within 0.19 cm⁻¹ in 150 iterations. `Example_12` does the same for stage 2: with durations up to 10³ the group commutator (not Pauli-Z) is fastest inside the block and reaches the 0.072 cm⁻¹ restriction floor after 16 + 12 iterations, and the same-irrep coupling stops being the slowest mode. DB-DOI improves in energy (72 against 98 cm⁻¹ after 5 steps) but still leaves 94 % of the state off the ground doublet. The two-stage strategy is still the fastest route, but the earlier explanations were wrong, and `RESULTS_Model_2.md` §3, §4, §6 and §8 now say so.

**Departures from the plan.** The free-ion refit (above), and the count of 108 rather than 109 measured levels.

## 6. Model 3: ¹⁶⁷Er³⁺:Y₂SiO₅, site 2, with hyperfine structure (8 qubits)

**Question.** Can DBQA resolve the hyperfine levels of the lowest crystal-field level Z1 inside the full 240-state Hamiltonian H(**B**), across a magnetic-field sweep, accurately enough to locate a zero-first-order-Zeeman (ZEFOZ) point, and how much carries over from one field to the next?

**Hamiltonian.** Er³⁺ at site 2 of Y₂SiO₅ (C₁, all 27 crystal-field parameters independent), ¹⁶⁷Er nuclear spin I = 7/2, magnetic-dipole and quadrupole hyperfine terms, electronic and nuclear Zeeman terms. Parameters: N. L. Jobbitt et al., Phys. Rev. B 104, 155121 (2021), arXiv:2206.09080, Tables III–IV. Basis: PyCF's `erhfs`, the intermediate-coupled ⁴I₁₅/₂ and ⁴I₁₃/₂ multiplets × |m_I⟩, 240 states in 256 dimensions. The scales span 6 500 cm⁻¹ (multiplets) to below 10⁻⁴ cm⁻¹ (hyperfine avoided crossings). The Hamiltonian is a single connected block, so there are no symmetry sectors.

**Validation** (`Example_8`, `RESULTS_Model_3.md` §2).
- Transcription: the paper's crystal-field strengths S₂, S₄, S₆ are reproduced exactly; Table II is transcribed in `data/`.
- **Hyperfine model: the paper's predicted 7 T hyperfine splittings (Table VI, 14 values) are reproduced to 0.7 MHz**, the rounding of the table. This works only without the nuclear Zeeman term, so the paper's model omits it; Model 3 includes it.
- Full 4f¹¹ basis (`f11cf`): crystal-field splittings inside every multiplet agree with Table II to 0.6 cm⁻¹ rms, and the g values of all 59 levels to 0.07. The multiplets, however, are rigidly shifted by −25 to +37 cm⁻¹, for both sites. No refit or truncation explains it; it is a convention difference of the bundled basis. Model 3 does not use this basis.
- ZEFOZ reference: K. Matsuura et al., Phys. Rev. B 113, 085421 (2026), arXiv:2412.10126 (first author verified), computed ZEFOZ points from a measured spin Hamiltonian (S.-J. Wang et al.'s tensors). Our search on that spin Hamiltonian returns their five best site-2 points in the same order. Their Cartesian field has the wrong sign of B_b; their polar form is right.
- **The crystal-field model has a different ZEFOZ landscape**: its best points lie at 90–150 mT near the D1–D2 plane, against 280–670 mT at θ ≈ −40° for the measured spin Hamiltonian. The fitted quadrupole tensor is several times smaller than the measured one. So Model 3 benchmarks DBQA on a realistic hyperfine Hamiltonian, but does not predict the crystal's ZEFOZ points.

**Results** (`Example_9`, `Example_10`, `RESULTS_Model_3.md` §3–5). Sweep: 21 fields, 120–180 mT, through the crystal-field model's best ZEFOZ point (|6⟩–|7⟩ at 149.9 mT). ν₆,₇ changes by only 0.083 MHz over the sweep.
- **Stage A** (8 qubits, group commutator, durations up to 10²) isolates Z1 to 1 − F < 10⁻⁶ in 81 iterations and to 10⁻⁹ in 105, at every field.
- **Isolation sets accuracy:** the Z1 transition frequencies inherit 300, 0.3, 2×10⁻⁴ and 3×10⁻⁵ MHz of error at 1 − F = 10⁻³, 10⁻⁶, 10⁻⁹ and 10⁻¹².
- **Stage B** (4 qubits): the group commutator reaches that floor (1.7×10⁻⁴ MHz from 10⁻⁹) in 60 iterations. The canonical bracket stalls at a fixed point, because after stage A the hyperfine structure sits in couplings between nearly degenerate diagonal entries. Pauli-Z is slow (50–85 MHz).
- **ZEFOZ from DBQA:** the turning point of ν₆,₇ along the sweep is located within 0.2 mT (10⁻⁶ isolation) and better than 0.01 mT (10⁻⁹).
- **Across the sweep:** warm-starting from the previous field's rotation needs no iteration to 10⁻⁶ and 8 to 10⁻⁹. Replaying the previous field's schedule of durations reaches 10⁻¹² with no new step, so on hardware the durations would only have to be found once; the circuit length stays that of a cold start.
- **Compilation:** 5 978 Pauli strings, 5.6×10⁴ CNOTs per first-order Trotter step, 19 times cheaper than Model 2; the 4-qubit stage-B block costs 1 018.
- **Avoided crossings** (`Example_13`): the crystal-field model and the measured spin Hamiltonian agree on where the 77 matched crossings of the Z1 levels are (fields within a few per cent, median ratio 0.976) but not on their gaps (crystal-field gaps typically 10 times smaller; 11 below 0.01 MHz). Across a 1.12 MHz and a 93.6 MHz crossing, two-stage DBQA resolves the pair at every field (gap errors below 6×10⁻⁴ MHz, pair infidelity down to 10⁻¹⁴); the canonical bracket and Pauli-Z isolate the pair but do not resolve it; warm starts make every field after the first nearly free.
- **Circuits** (`Example_14`): one second-order Trotter step (1.1×10⁵ CNOTs) is an accurate oracle. Full 8-qubit DB-DOI keeps the state in ⁴I₁₅/₂ but moves it only 15 cm⁻¹ towards Z1 in 5 steps (1.3×10⁷ CNOTs). Two-stage DB-DOI inside the 4-qubit Z1 block converges slowly: infidelity 0.98 → 0.26 in 8 steps (3 280 oracle calls).
- **Noise** (`Example_15`): with a depolarizing channel after every CNOT, the two-stage block circuits lose their signal already at p = 10⁻⁵ (8.5×10⁸ CNOTs for Model 3 at k = 8). Keeping half of the noiseless signal needs p ≲ 10⁻⁵ for Model 1 (with the extended grid) and 10⁻⁸–10⁻¹⁰ for Models 2 and 3: circuit-level DBQA of these Hamiltonians is a fault-tolerant-era task.

**Departures from the plan.**
- *ZEFOZ comparison.* The plan was to compare DBQA-located ZEFOZ points with Matsuura et al. Their points are reproduced from their spin Hamiltonian, but the crystal-field model of Jobbitt et al. does not have them. DBQA was therefore benchmarked at a ZEFOZ point of the crystal-field model itself.
- *Symmetry-resolved reference states*, planned after Model 2, are moot for Model 3, which has no symmetry sectors (C₁ plus a field). They matter again for the two-stage DB-DOI of Model 2: its multiplet block keeps the four sectors of the full Hamiltonian (`Example_14`).
- *Warm start.* Planned as reusing the previous rotation; on hardware that keeps every earlier step in the circuit. We added the schedule replay, which turned out to transfer exactly.
- *Model 2.* Model 3 exposed the duration-grid limit and led to the Model 2 correction above.
- *Paper details.* The paper gives 224 states for the two-multiplet hyperfine basis; PyCF's basis, from the authors' group, has 240 = (16 + 14) × 8.

## 7. Deviations from the original plan, collected

| Model | Planned | Done instead | Why | Where |
|---|---|---|---|---|
| 1 | Hu et al. parameters | column [a] (fitted) as main set; [c] only as sensitivity test | [c] does not reproduce the paper's own Table 9[c] levels (up to 760 cm⁻¹ off) | `RESULTS_Model_1.md` §1 |
| 2 | Couto dos Santos et al. parameters as published | paper's crystal field + free-ion parameters refitted in PyCF (σ 12.9 vs 14.6 cm⁻¹) | the paper's E³ convention cannot be reproduced (σ = 547 cm⁻¹ as printed) | `RESULTS_Model_2.md` §2 |
| 2 | 109 measured levels (abstract) | 108 | Table 1 lists 108; the conclusion says 108 | `RESULTS_Model_2.md` §2.1 |
| 1 | conclusion: the closest pair limits convergence | corrected: the duration grid did; with s ≤ 10² all levels converge in 15 iterations | found while running Model 3 | `RESULTS_Model_1.md` §5, §7; `Example_12` |
| 2 | conclusion: greedy DBQA cannot resolve crystal-field levels; Pauli-Z best in stage 2 | corrected: both were grid limits | found while running Model 3 | `RESULTS_Model_2.md` §3, §4, §8; `Example_11`, `Example_12` |
| 3 | compare DBQA ZEFOZ points with Matsuura et al. | Matsuura et al. reproduced from their spin Hamiltonian; DBQA benchmarked at the crystal-field model's own ZEFOZ point | the crystal-field parameters do not reproduce the measured ZEFOZ landscape | `RESULTS_Model_3.md` §2.5–2.6 |
| 3 | full 4f¹¹ validation against Table II | done; multiplets shifted rigidly, unexplained | convention of the bundled `f11cf` basis | `RESULTS_Model_3.md` §2.2 |
| 3 | symmetry-resolved reference states | not applicable | one connected block | §6 |
| 3 | DBQA near avoided crossings with pair metrics; compare crossing fields with Matsuura et al. | done (`Example_13`); crossing fields compared with the spin Hamiltonian Matsuura et al. use, since they give crossings only in figures | §6 |

## 8. Limitations and open items

- **Circuit level.** DB-DOI was run for all three models with Z-field generators and greedy duration grids. It removes the large energy scales quickly but resolves small ones slowly or not at all (Model 2 and 3 full circuits; the two-stage blocks reach about 70 % fidelity after 8 steps). The group commutator with D = Δ(Ĥ_k), which does most of the work at matrix level, is not directly available on hardware.
- **Noise.** A depolarizing channel after every CNOT, simulated for the 4-qubit circuits (`Example_15`); the 8- and 9-qubit circuits are covered only by the global-depolarizing estimate. No coherent errors, readout errors or connectivity constraints were modelled.
- **Duration grids.** The full 8-qubit DB-DOI of Model 3 used durations ≤ 1, as in `Example_7`; repeating it with the extended grid was not done because the run already peaks near 19 GB of memory.
- **Scheduling.** Only greedy grid search was used. The duration grids matter more than any other setting we varied.
- **Physical accuracy of Model 3.** The ZEFOZ points and crossing gaps of Model 3 are those of Jobbitt et al.'s parameters, not of the crystal. A refit including the measured site-2 A and Q tensors would make the model predictive.
- **Basis encoding.** All models use a binary map from basis index to qubits. The Pauli content (and so the CNOT cost) of Models 2 and 3 might drop with another map; this was not tested.

## 9. Relation to previous work

Quantum algorithms have been applied to crystal-field Hamiltonians before: I. Makushin and E. Baibekov ran the variational quantum eigensolver on lowest-multiplet crystal-field and Zeeman Hamiltonians of Yb³⁺ and Er³⁺ on IBM hardware (arXiv:2005.03712). In searches on 17 October 2026 we found no application of double-bracket algorithms to crystal-field or hyperfine Hamiltonians. We do not claim more than that; a proper literature review should precede any novelty statement in a paper.

## 10. Data not in the repository

The two large Model 2 caches, `Er_LiYF4_DBQA_checkpoints.npz` (17 MB, the iterated matrices of Example 5) and `Er_LiYF4_trotter_steps.npz` (33 MB, the compiled Trotter steps of Example 7), are not under version control. Examples 5 and 7 regenerate them when they are run; every other notebook runs from the files in `simulation_results/`.
