# Model 1: DBQA for Ce³⁺:YAG (4f¹): results

All numbers below come from the executed notebooks in this folder:
- `Example_1_validation_PyCF_to_Qibo_Ce_YLF.ipynb`: conversion check
- `Example_2_DBQA_Ce_YAG.ipynb`: DBQA on the matrix
- `Example_3_compiled_DBQA_Ce_YAG.ipynb`: compiled circuits

The data are in `simulation_results/` and the figures in `pdfs_figures/`. Every run is deterministic, because grid-search scheduling and the generator selection involve no random sampling.

## 1. Model

**System and basis.** Ce³⁺ sits at the D₂ site of Y₃Al₅O₁₂. The 4f¹ configuration has 14 |²F_J M_J⟩ states, taken from PyCF's bundled `f1cf` basis. They are embedded in 16 dimensions (4 qubits).

**Padding.** The two padding states are uncoupled from the physical states. They sit on the diagonal a quarter of the physical spectral width (1012 cm⁻¹) above the highest level.

**Energy unit.** Energies enter Qibo in units of 1000 cm⁻¹.

**Parameters.** These are from Table 10 of Hu, Reid, Duan, Xia and Yin, J. Phys.: Condens. Matter 23, 045501 (2011), in their axis convention. The main runs use column [a], labelled "fitted": the fit of Tanner et al. (2007), rotated by π/2 about z into Hu's axes. Column [c], their DV-Xα calculation, is the "calculated" set and is used as a sensitivity check.

| | ζ | B²₀ | B²₂ | B⁴₀ | B⁴₂ | B⁴₄ | B⁶₀ | B⁶₂ | B⁶₄ | B⁶₆ |
|---|---|---|---|---|---|---|---|---|---|---|
| fitted [a] | 647 | −465 | 96 | −3739 | 380 | 1602 | 901 | −307 | 2136 | −246 |
| calculated [c] | 752 | −275 | 449 | −2154 | −74 | 1038 | 842 | 79 | 1613 | −740 |

**Exact spectrum (fitted set).** Seven Kramers doublets, in cm⁻¹ above the ground level:
- ²F₅/₂: 0, 321.9, 733.7
- ²F₇/₂: 2250.4, 2492.8, 2586.8, 4046.9
- Padding: 5058.6

The smallest gap is 94.0 cm⁻¹, between doublets 4 and 5. The largest splitting inside a doublet is 3×10⁻¹² cm⁻¹. Levels are grouped with a tolerance of 0.01 cm⁻¹, fixed from the exact spectrum before any DBQA run.

**Note on the calculated set.** Column [c] does not reproduce the calculated 4f levels that Hu et al. list in Table 9[c]; the differences reach 760 cm⁻¹. No sign or normalisation change removes this. Separately, the Wybourne normalisation check below shows the conventions on our side are standard. The likely cause is that their effective Hamiltonian was truncated to ζ and the nine B^k_q. This does not affect the benchmark, whose reference is the exact spectrum of the matrix given to DBQA.

## 2. Validation of the PyCF → Qibo conversion (Example 1, Ce:YLF)

| Check | Result |
|---|---|
| Dense matrix vs PyCF eigenvalues | max deviation 7×10⁻¹² cm⁻¹ |
| Dense matrix vs PyCF eigenvectors | max ‖Hz − zw‖ = 8×10⁻¹² cm⁻¹ |
| Eigenvalues after qubit embedding | max deviation 6×10⁻¹² cm⁻¹ |
| Doublet eigenspaces after embedding | projector fidelity 1.000000000000 for all 7 doublets |
| Padding coupling to physical states | exactly 0 |
| Wybourne normalisation (‖C^k_q + C^k_−q‖² / ‖C^k_0‖²) | 2.000000 for every q > 0, k = 2, 4, 6 |
| Ce:YLF levels vs PyCF example data | 0, 212, 413, 2214, 2311, 2429, 3157 cm⁻¹; measured: 0, 216, –, 2216, 2313, 2429, 3158 |

## 3. DBQA on the matrix (Example 2)

**Methods.** Three double-bracket rotations from `qibo_dbqa` were run for 60 iterations each. In every case the step was chosen by grid search to minimise ‖σ(Ĥ_k)‖.
- **canonical:** the GWW bracket [Δ(Ĥ_k), Ĥ_k].
- **Pauli-Z:** a single commutator with D_k chosen each step from all 15 products of Pauli Z (`select_best_dbr_generator`).
- **group commutator:** D_k = Δ(Ĥ_k).

**Metrics.** For each level we record:
- the energy error, from the diagonal of Ĥ_k;
- the projector fidelity of each Kramers doublet;
- the error of every gap between consecutive levels;
- the off-diagonal norm ‖σ(Ĥ_k)‖.

### 3.1 Iterations until the error stays below a threshold (fitted set)

| Quantity | Threshold | canonical | Pauli-Z | group commutator |
|---|---|---|---|---|
| ²F₅/₂ max \|ΔE\| | 10 cm⁻¹ | 45 | 14 | 14 |
| | 1 cm⁻¹ | > 60 | 32 | 25 |
| | 0.1 cm⁻¹ | > 60 | > 60 | 43 |
| ²F₅/₂ max 1 − F | 10⁻² | > 60 | 24 | 21 |
| | 10⁻³ | > 60 | 46 | 36 |
| | 10⁻⁴ | > 60 | > 60 | 54 |
| All 7 levels max \|ΔE\| | 10 cm⁻¹ | > 60 | 47 | > 60 |
| | 1 cm⁻¹ | > 60 | > 60 | > 60 |
| All 7 levels max 1 − F | 10⁻² | > 60 | > 60 | > 60 |

### 3.2 Off-diagonal norm, 0 → 60 iterations

| Method | ‖σ(Ĥ₀)‖ → ‖σ(Ĥ₆₀)‖ | Reduction factor |
|---|---|---|
| canonical | 2.392 → 0.142 | 16.9 |
| Pauli-Z | 2.392 → 0.047 | 50.6 |
| group commutator | 2.392 → 0.083 | 29.0 |

### 3.3 Errors after 60 iterations, per doublet (levels 0–6)

| Method | Quantity | Doublet 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|---|
| canonical | \|ΔE\| (cm⁻¹) | 1.7 | 3.6 | 5.3 | 1.5 | 26 | 28 | 0.009 |
| | 1 − F | 5.4e−3 | 1.8e−2 | 1.3e−2 | 5.0e−3 | 0.29 | 0.29 | 2.2e−6 |
| Pauli-Z | \|ΔE\| (cm⁻¹) | 0.14 | 0.20 | 0.20 | 0.49 | 3.7 | 4.2 | 5e−4 |
| | 1 − F | 2.6e−4 | 5.5e−4 | 3.0e−4 | 1.8e−3 | 4.4e−2 | 4.2e−2 | 1.3e−7 |
| group commutator | \|ΔE\| (cm⁻¹) | 0.012 | 0.011 | 0.001 | 0.032 | 24 | 24 | 2e−4 |
| | 1 − F | 3.8e−5 | 4.1e−5 | 2.9e−6 | 1.2e−4 | 0.26 | 0.26 | 5.8e−8 |

Errors between consecutive levels after 60 iterations, in cm⁻¹ (gap k means levels k and k+1):

| Method | gap 0 | gap 1 | gap 2 | gap 3 | gap 4 | gap 5 |
|---|---|---|---|---|---|---|
| canonical | 1.9 | 9.0 | 6.9 | 25 | 54 | 28 |
| Pauli-Z | 0.045 | 0.15 | 0.54 | 3.1 | 7.5 | 4.0 |
| group commutator | 0.023 | 0.010 | 0.034 | 24 | 48 | 24 |

### 3.4 The closest pair (doublets 4 and 5, gap 94 cm⁻¹)

| Method | Individual 1 − F | Pair-subspace 1 − F | Gap error |
|---|---|---|---|
| canonical | 0.29, 0.29 | 2.5×10⁻³ | 53.8 cm⁻¹ |
| Pauli-Z | 0.044, 0.042 | 8.9×10⁻⁴ | 7.5 cm⁻¹ |
| group commutator | 0.26, 0.26 | 6.0×10⁻⁵ | 48.4 cm⁻¹ |

The group commutator isolates the pair from the rest of the spectrum best: its pair-subspace infidelity is 6×10⁻⁵. It has not resolved the pair internally, since the two doublets are still mixed and the gap error is 48 cm⁻¹. Pauli-Z resolves the pair best but isolates it less well. This is the case the two-state subspace fidelity and gap error were proposed for, and they separate the two failure modes.

### 3.5 Relation between fidelity and energy error

Diagonal entries of Ĥ_k are expectation values in the approximate eigenvectors, so the energy error is second order in the eigenvector error. Over all iterates and all ²F₅/₂ levels of the three methods, the bound |ΔE| ≤ c·(1 − F)·Δ holds with c = 0.22 at most (median c = 0.08). Here Δ is the largest distance from the level to any other physical level (figure `Ce_YAG_energy_vs_fidelity.pdf`).

### 3.6 Sensitivity to the parameter set (calculated set)

The calculated set has levels at 0, 202.5, 577.5, 2591.8, 2743.6, 2961.0 and 3737.1 cm⁻¹, with a smallest gap of 151.8 cm⁻¹.

| Quantity | Threshold | canonical | Pauli-Z | group commutator |
|---|---|---|---|---|
| ²F₅/₂ max \|ΔE\| | 10 cm⁻¹ | > 60 | 25 | 20 |
| | 1 cm⁻¹ | > 60 | > 60 | 58 |
| ²F₅/₂ max 1 − F | 10⁻² | > 60 | 56 | 44 |
| All 7 levels max \|ΔE\| | 10 cm⁻¹ | > 60 | 46 | > 60 |

For the ground multiplet the ranking of the methods is the same for both sets: group commutator, then Pauli-Z, then canonical. With the calculated set, though, every method needs 1.4–2.3× more iterations for the ground multiplet. Iteration counts therefore depend strongly on the specific parameters, even for this 4-qubit model.

## 4. Compiled circuits (Example 3)

### 4.1 Evolution oracle

**Pauli content.** The embedded Hamiltonian is real. It has 71 non-zero Pauli strings out of 256; by weight 0 to 4 there are 1, 6, 15, 26 and 23 strings. The identity coefficient, 0.410, only contributes a global phase and is not compiled.

**Compilation.** `PauliCompiledHamiltonian` compiles a Trotter–Suzuki step from Pauli rotations, each a CNOT ladder costing 2(w − 1) CNOTs for a weight-w string. The table gives the operator-norm error against e^{−it(Ĥ − c_I)}. CNOT counts are without cancellation between neighbouring rotations.

| Order | Steps | Error, t = 0.1 | t = 0.3 | t = 0.55 | CNOT | Depth |
|---|---|---|---|---|---|---|
| 1 | 1 | 1.3e−2 | 0.11 | 0.35 | 272 | 437 |
| 1 | 2 | 6.3e−3 | 5.5e−2 | 0.17 | 544 | 872 |
| 1 | 4 | 3.1e−3 | 2.7e−2 | 8.3e−2 | 1088 | 1742 |
| 1 | 8 | 1.6e−3 | 1.4e−2 | 4.2e−2 | 2176 | 3482 |
| 2 | 1 | 5.1e−4 | 1.3e−2 | 7.9e−2 | 538 | 867 |
| 2 | 2 | 1.3e−4 | 3.3e−3 | 1.9e−2 | 1076 | 1734 |
| **2** | **4** | **3.2e−5** | **8.3e−4** | **4.7e−3** | **2152** | **3468** |
| 2 | 8 | 7.9e−6 | 2.1e−4 | 1.2e−3 | 4304 | 6936 |

**Choice.** We use second order with 4 steps, which keeps the error below 5×10⁻³ up to the largest duration DB-DOI uses (√0.3 = 0.55).

**Comparison with direct synthesis.** A generic 4-qubit unitary needs at least 61 CNOTs (Shende–Markov–Bullock bound). Direct synthesis of each oracle call would therefore be far cheaper for this small model; Trotterisation is the compilation that scales to Models 2 and 3.

### 4.2 DB-DOI ground-state preparation

**Set-up.** We use `GroupCommutatorIteration` with the reduced group commutator, the same algorithm as the Fig. 4 and Quantinuum examples. D_k is a magnetic field Σ b_i Z_i, with the signs b_i = ±1 and the step s_k ∈ [0.001, 0.3] chosen each step to minimise the energy. The reference state is |0101⟩, the basis state with the lowest diagonal energy.

**How the table was built.** The ideal oracle (`ExactEvolutionHamiltonian`) gives the convergence. For k ≤ 3 the compiled circuits were built explicitly with the order-2, 4-step oracle. Their CNOT count equals queries × 2152 exactly, because the magnetic fields and state preparation contain no two-qubit gates. For k = 4 to 7 the CNOT counts follow from that identity; those depths were not built.

| k | E_k − E₀ (cm⁻¹) | 1 − F (ground doublet) | Witness | Oracle queries | CNOT | Depth |
|---|---|---|---|---|---|---|
| 0 | 19.67 | 3.50e−2 | 0.939 | 0 | 0 | 1 |
| 1 | 11.25 | 2.83e−2 | 0.965 | 1 | 2152 | 3475 |
| 2 | 8.38 | 2.18e−2 | 0.974 | 4 | 8608 | 13897 |
| 3 | 6.53 | 1.90e−2 | 0.980 | 13 | 27976 | 45163 |
| 4 | 5.04 | 1.44e−2 | 0.984 | 40 | 86080 | not built |
| 5 | 3.83 | 1.05e−2 | 0.988 | 121 | 260392 | not built |
| 6 | 2.95 | 7.75e−3 | 0.991 | 364 | 783328 | not built |
| 7 | 2.28 | 6.68e−3 | 0.993 | 1093 | 2352136 | not built |

The witness is 1 − (E − E₀)/(E₁ − E₀), as in `utils_reporting.fidelity_witness`.

**Trotter error.** It is negligible here. For the explicit compiled circuits at k = 1, 2, 3 the energies are 11.244, 8.377 and 6.524 cm⁻¹, against 11.247, 8.377 and 6.529 with the ideal oracle. The infidelities agree to the digits shown.

**Cost.** Oracle queries follow q_{k+1} = 3q_k + 1, because the evolution under Ĥ_k is the evolution under Ĥ conjugated by the circuit prepared so far. The CNOT count therefore grows by a factor of 3 per step, while after the first step the energy error falls by only a factor of about 1.3 per step.

## 5. Answer: how quickly does DBQA recover the low-energy spectrum?

**Matrix level.** Ground-multiplet energies are recovered to 1 cm⁻¹ after 25 iterations with the group commutator and 32 with Pauli-Z generators, on the fitted set. The doublet eigenspaces reach 1 − F < 10⁻² after 21 and 24 iterations. The canonical bracket with greedy steps zig-zags and needs more than 60 iterations.

**What limits convergence (corrected 2026-10-17, see §7).** We concluded that convergence is limited by the smallest gap: the 94 cm⁻¹ pair in ²F₇/₂ is isolated within 60 iterations but not resolved. The limit was the duration grid (s ≤ 1, while resolving a gap of 0.094 in units of 1000 cm⁻¹ needs s ≈ 10). With durations up to 10², the group commutator resolves all seven levels, the pair included, within 1 cm⁻¹ after 15 iterations.

**Circuit level.** DB-DOI from the best basis state reaches 1 − F < 10⁻² for the ground doublet after 6 steps (364 oracle queries, about 7.8×10⁵ CNOTs with the Trotterised oracle).

**Tolerances proposed for Models 2 and 3**, to be fixed before their DBQA runs:
- **Grouping.** Group levels with a tolerance above the largest intra-level splitting and below the smallest gap claimed to be resolved, checked with `check_grouping_tolerance` as here.
- **Energy.** Use an energy tolerance of 1% of the smallest resolved gap g_min (≈ 1 cm⁻¹ for Model 1).
- **Fidelity.** Use a projector-fidelity tolerance of 1 − F ≤ 10⁻² per level. By §3.5 this implies |ΔE| ≲ 0.2·(1 − F)·Δ, so for wide spectra (Model 2) the energy tolerance is the binding one.
- **Close pairs.** For any pair closer than about 2× the smallest gap here, report the pair-subspace fidelity and the gap error as well. Count the pair as resolved only when the gap error is below 10% of the gap.
- **Duration grid (added 2026-10-17).** Let the grid reach about 1/g for the smallest gap g to be resolved, in the units the Hamiltonian is rescaled to (§7).
- **Iteration budget.** Budget at least 60 iterations per strategy and compare canonical, Pauli-Z and group-commutator rotations. Section 3.6 shows the counts change by up to 2.3× between parameter sets of the same ion.

## 6. Limitations

- **Scheduling.** Steps are chosen greedily by grid search: 300 points in [10⁻⁴, 0.3], or [10⁻³, 1] for the group commutator. Other schedulers (hyperopt, polynomial) were not compared.
- **Diagonal operators in DB-DOI.** These are restricted to ±1 magnetic fields. The gradient-optimised D of `utils_gci_optimization` was not used.
- **Compiled gate counts.** CNOT counts carry no peephole optimisation. Depths beyond k = 3 were not computed, because those circuits have more than 10⁵ gates.
- **Noise.** None in §4; depolarizing noise after every CNOT is simulated in `Example_15` (§7).
- **Single benchmark Hamiltonian.** Experimental validation of the Ce:YAG parameters is outside the scope: the benchmark is the exact spectrum of the PyCF Hamiltonian.

## 7. Addendum (2026-10-17): duration grids (Example 12)

Models 2 and 3 showed that a double-bracket step resolves a gap g (in the rescaled units) only if the duration grid reaches about 1/g. Here the closest pair is 94 cm⁻¹ = 0.094 in units of 1000 cm⁻¹, while the grids of §3 stopped at 0.3 (canonical, Pauli-Z) and 1 (group commutator). `Example_12_duration_grids_Ce_YAG_Er_LiYF4.ipynb` repeats §3 with `run_dbqa_history` on the original grids, which reproduces the iteration counts and pair metrics of §3.1 and §3.4, and on 300 log-spaced durations in [10⁻⁴, 10²]:

| Grid | Method | ²F₅/₂ < 1 cm⁻¹ | < 0.1 cm⁻¹ | all levels < 1 cm⁻¹ | final max \|ΔE\|, all levels | pair gap error | pair 1 − F | largest duration used |
|---|---|---|---|---|---|---|---|---|
| original | canonical | > 60 | > 60 | > 60 | 27.6 cm⁻¹ | 53.8 cm⁻¹ | 2.5×10⁻³ | 0.19 |
| original | Pauli-Z | 32 | > 60 | > 60 | 4.2 cm⁻¹ | 7.5 cm⁻¹ | 8.9×10⁻⁴ | 0.30 |
| original | group commutator | 25 | 43 | > 60 | 24.2 cm⁻¹ | 48.4 cm⁻¹ | 6.0×10⁻⁵ | 1.00 |
| ≤ 10² | canonical | > 60 | > 60 | > 60 | 27.6 cm⁻¹ | 53.8 cm⁻¹ | 2.5×10⁻³ | 0.20 |
| ≤ 10² | Pauli-Z | 23 | 49 | > 60 | 1.8 cm⁻¹ | 3.5 cm⁻¹ | 2.8×10⁻⁴ | 0.94 |
| ≤ 10² | **group commutator** | **10** | **12** | **15** | **< 10⁻³ cm⁻¹** | **< 10⁻³ cm⁻¹** | **2×10⁻¹⁵** | 75.8 |

With long enough durations the group commutator diagonalises the whole Hamiltonian in 15 iterations; the closest pair is not a bottleneck. The canonical bracket does not use the longer durations (its greedy optimum stays below 0.2) and is unchanged.

**DB-DOI.** The same DB-DOI as §4.2 (ideal oracle) with 40 log-spaced durations in [10⁻³, 10²] instead of 20 in [10⁻³, 0.3]:

| Step k | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|
| E_k − E₀, original grid (cm⁻¹) | 11.25 | 8.38 | 6.53 | 5.04 | 3.83 | 2.95 | 2.28 |
| E_k − E₀, durations ≤ 10² | 11.26 | 5.09 | 2.84 | 1.72 | 1.17 | 0.71 | 0.44 |
| 1 − F, original grid | 2.8×10⁻² | 2.2×10⁻² | 1.9×10⁻² | 1.4×10⁻² | 1.1×10⁻² | 7.8×10⁻³ | 6.7×10⁻³ |
| 1 − F, durations ≤ 10² | 2.8×10⁻² | **7.6×10⁻³** | 5.0×10⁻³ | 3.2×10⁻³ | 2.4×10⁻³ | 1.7×10⁻³ | 1.1×10⁻³ |

The ground doublet reaches 1 − F < 10⁻² after 2 steps (4 oracle calls) instead of 6 (364 calls). Longer evolutions need more Trotter steps per call: at the longest one used (s = 7.0, t = 2.65) the oracle needs 32 second-order steps (17 216 CNOTs per call, error 3.7×10⁻³) instead of 4 (2 152). The net cost of 1 − F < 10⁻² still falls from 7.8×10⁵ to **6.9×10⁴ CNOTs**, 11 times less (`Example_15`).

**Noise (Example 15).** With a depolarizing channel of strength p after every CNOT, the best fidelity with the ground doublet over the DB-DOI steps is 0.958 (p = 10⁻⁵), 0.84 (p = 10⁻⁴) and 0.31 (p = 10⁻³), always at k = 1 (original grid; noiseless 0.993 at k = 7). Keeping half of the noiseless signal at the step that reaches 1 − F < 10⁻² requires p ≲ 9×10⁻⁷ on the original grid and p ≲ 10⁻⁵ with the extended grid.

## Reproduction

```bash
export PYCF_ROOT=/path/to/pycf
cd examples/Crystal_field_PyCF
jupyter nbconvert --to notebook --execute --inplace Example_1_validation_PyCF_to_Qibo_Ce_YLF.ipynb Example_2_DBQA_Ce_YAG.ipynb Example_3_compiled_DBQA_Ce_YAG.ipynb Example_12_duration_grids_Ce_YAG_Er_LiYF4.ipynb
cd ../.. && pytest -o addopts="" tests/
```

On this laptop Example 2 runs in about 1.5 minutes and Example 3 in about 7 minutes. The test suite passes in full: 110 tests on the numpy and qibojit backends, 4 minutes.
