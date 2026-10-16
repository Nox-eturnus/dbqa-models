# Model 2: DBQA for Er³⁺:LiYF₄ (full 4f¹¹): results

All numbers below come from the executed notebooks in this folder:
- `Example_4_validation_Er_LiYF4.ipynb`: Hamiltonian and validation
- `Example_5_DBQA_Er_LiYF4_generate_data.ipynb`: DBQA runs
- `Example_6_DBQA_Er_LiYF4_results.ipynb`: analysis and compilation cost

The data are in `simulation_results/` and the figures in `pdfs_figures/`. Every run is deterministic, because grid-search scheduling and the generator selection involve no random sampling.

## 1. Model

**System and basis.** Er³⁺ sits at the S₄ site of LiYF₄, treated in the D₂d approximation as in the source paper. The full 4f¹¹ configuration has 364 states in the |SLJM⟩ basis, taken from PyCF's bundled `f11cf` basis. They are embedded in 512 dimensions (9 qubits).

**Padding.** The 148 padding states are uncoupled from the physical states.

**Energy unit.** Energies enter Qibo in units of 20 000 cm⁻¹. The traceless part of the embedded Hamiltonian then has norm 3.07, comparable to Model 1. In the PyCF basis, 5.6 % of the matrix elements are non-zero.

**Hamiltonian.** Free ion (Eⁿ, α, β, γ, Tⁱ, ζ, Mᵏ, Pᵏ) plus the D₂d crystal field: B²₀ C²₀ + B⁴₀ C⁴₀ + B⁴₄ (C⁴₄ + C⁴₋₄) + B⁶₀ C⁶₀ + B⁶₄ (C⁶₄ + C⁶₋₄). It follows M. A. Couto dos Santos, E. Antic-Fidancev, J. Y. Gesland, J. C. Krupa, M. Lemaître-Blaise and P. Porcher, J. Alloys Compd. 275–277, 435–441 (1998), doi:10.1016/S0925-8388(98)00363-6. Their Table 1 (levels) is transcribed in `data/Er_LiYF4_CoutoDosSantos1998_Table1.csv`.

## 2. Validation of the Hamiltonian (Example 4)

### 2.1 The paper's data

**Number of measured levels.** Table 1 lists 134 Kramers doublets up to 44 000 cm⁻¹, of which **108** have an experimental energy. The abstract and Table 2 say 109; the conclusion says 108.

**Transcription check.** The paper's own calculated levels against its experimental ones give σ = 14.63 cm⁻¹, with σ = √(ΣΔ²/(N − p)), N = 108 and p = 21 free parameters. This reproduces the paper's stated σ = 14.6 cm⁻¹, so the transcription and the definition of σ are confirmed.

### 2.2 The paper's parameters fail in PyCF

**Conversion.** The Racah parameters were converted to PyCF's Slater–Condon Fᵏ with the standard relations: E¹ = (70F₂ + 231F₄ + 2002F₆)/9, E² = (F₂ − 3F₄ + 7F₆)/9, E³ = (5F₂ + 6F₄ − 91F₆)/3, where F₂ = F²/225, F₄ = F⁴/1089 and F₆ = 25F⁶/184041. M², M⁴, P⁴ and P⁶ follow from M⁰ and P² with the paper's fixed ratios 0.56, 0.31, 0.50 and 0.10.

**Dense matrix.** The dense matrix agrees with PyCF's own diagonalisation to 4×10⁻¹⁰ cm⁻¹.

**Mismatch.** With the paper's parameters as printed, the PyCF levels differ from the paper's calculated levels by up to 1602 cm⁻¹, and from experiment by σ = 547 cm⁻¹.

**Diagnosis.** We refitted groups of parameters to the paper's own calculated levels:

| Fit to the paper's E_calc | rms | Outcome |
|---|---|---|
| E¹, E², E³ only | 15.5 cm⁻¹ | E¹ 6607 → 6615 and E² 32.99 → 32.92 (unchanged); E³ 664.9 → 619.3 |
| crystal field only | 515 cm⁻¹ | no repair; B⁴₀ is pushed to −116 cm⁻¹ |

So the paper's E³ is not the standard Racah/Carnall E³. As a cross-check, Carnall's Er:LaF₃ Slater parameters (97483, 67904, 54010) give E³ = 624, also about 40 cm⁻¹ below the paper's value. The convention used by the authors' fitting program cannot be recovered from the paper, so we do not guess it.

### 2.3 The Model 2 Hamiltonian

We keep the paper's crystal field, whose convention matches PyCF's (checked in Example 1), and refit the 15 free-ion parameters to the 108 experimental levels with the M/P ratios fixed. σ always counts the paper's 21 parameters.

| Variant | σ (cm⁻¹) | rms (cm⁻¹) | max \|Δ\| (cm⁻¹) |
|---|---|---|---|
| paper parameters as given | 546.89 | 490.85 | 1677.8 |
| refit E¹, E², E³ only | 17.08 | 15.33 | 38.4 |
| **refit free ion, paper crystal field (Model 2)** | **12.91** | 11.59 | 38.3 |
| refit all 20 | 12.75 | 11.45 | 39.6 |
| `ER_LIYF4_FREE_ION_REFIT` (rounded to 0.01) | 12.92 | 11.60 | 37.7 |

**Result.** The Model 2 Hamiltonian reproduces the experimental levels better than the paper's own fit (12.9 against 14.6 cm⁻¹).

**The crystal field is consistent.** Freeing the crystal field as well gains only 0.16 cm⁻¹ in σ. In that fit, B²₀, B⁴₀, B⁴₄ and B⁶₄ move by at most 3 %, and only the small B⁶₀ (−43.6 → −14.0 cm⁻¹) changes appreciably.

**Free-ion parameters.** The refitted Slater parameters are F² = 97 315, F⁴ = 67 960 and F⁶ = 53 525 cm⁻¹, typical Er³⁺ values. The largest relative changes from the paper are in poorly determined parameters: γ (+8 %), T² (−18 %) and P² (−10 %).

| Parameter | Paper | Model 2 |
|---|---|---|
| E¹, E², E³ | 6607.22, 32.99, 664.86 | 6583.09, 32.91, 625.12 |
| α, β, γ | 18.29, −589.84, 1755.61 | 18.08, −596.69, 1892.51 |
| T², T³, T⁴ | 455.12, 38.43, 83.67 | 372.53, 41.74, 71.19 |
| T⁶, T⁷, T⁸ | −348.26, 243.92, 388.37 | −353.31, 234.42, 382.67 |
| ζ, M⁰, P² | 2378.00, 4.10, 801.82 | 2375.21, 3.76, 721.45 |
| B²₀, B⁴₀, B⁴₄, B⁶₀, B⁶₄ | 323.07, −753.83, −1052.24, −43.60, −623.59 | unchanged |

### 2.4 Level structure and groupings

**Groupings.** All three were fixed before any DBQA run.
- **Levels:** 182 Kramers doublets, grouping tolerance 10⁻³ cm⁻¹. The largest splitting inside a doublet is 3×10⁻¹⁰ cm⁻¹ and the smallest gap between doublets is 0.03 cm⁻¹.
- **Manifolds:** 35 clusters separated by more than 300 cm⁻¹. That threshold is larger than every gap inside the ground multiplet (at most 203 cm⁻¹), so ⁴I₁₅/₂ is one manifold. The smallest gap between manifolds is 313 cm⁻¹. Manifolds separate the free-ion scale (which multiplet) from the crystal-field scale (which level inside it).
- **Clusters:** levels closer than 1 cm⁻¹ are merged, giving 180 physical clusters.

**Ground multiplet ⁴I₁₅/₂.** Its levels are at 0, 0.36, 5.98, 44.91, 247.70, 283.76, 320.20 and 339.21 cm⁻¹. The whole spectrum spans about 100 000 cm⁻¹, so the crystal-field splittings are 10⁻³–10⁻⁶ of the total energy scale.

## 3. DBQA on the full Hamiltonian (stage 1)

**Methods.** Each method ran for 150 iterations, with the step at every iteration chosen by grid search to minimise ‖σ(Ĥ_k)‖.
- **canonical:** the GWW bracket.
- **Pauli-Z:** a single commutator with the best single-qubit Zᵢ, either sign. These are the magnetic-field generators of DB-DOI; all 511 Z products are too many at 512 dimensions.
- **group commutator:** D_k = Δ(Ĥ_k).

**Scheduler.** `spectral_grid_search_step` minimises the same loss as `grid_search_step`; the tests check that it gives identical losses. It diagonalises once per step, so it is affordable at this size.

**Grids.** 100 log-spaced durations in [10⁻⁴, 10]; 30 per sign for Pauli-Z. The three methods ran in parallel in 35 minutes.

### 3.1 Iterations until the error stays below a threshold

| Quantity | Threshold | canonical | Pauli-Z | group commutator |
|---|---|---|---|---|
| ground manifold 1 − F | 10⁻² | 1 | 5 | 1 |
| | 10⁻⁶ | > 150 | > 150 | **26** |
| all 35 manifolds, max 1 − F | 10⁻² | > 150 | > 150 | > 150 |
| ⁴I₁₅/₂ levels, max \|ΔE\| | 10 cm⁻¹ | > 150 | > 150 | > 150 |
| ⁴I₁₅/₂ levels, max 1 − F | 10⁻² | > 150 | > 150 | > 150 |
| ⁴I₁₅/₂ clusters, max 1 − F | 10⁻² | > 150 | > 150 | > 150 |
| all 182 levels, max \|ΔE\| | 100 cm⁻¹ | > 150 | > 150 | > 150 |

### 3.2 State after 150 iterations

| Method | ‖σ‖ (start 9.89) | 1 − F, first 6 manifolds | ⁴I₁₅/₂ \|ΔE\| range | ⁴I₁₅/₂ 1 − F range |
|---|---|---|---|---|
| canonical | 0.82 | 1.4e−5 … 0.20 | 81–120 cm⁻¹ | 0.90–1.0 |
| Pauli-Z | 2.16 | 4.1e−4 … 0.37 | 82–130 cm⁻¹ | 0.65–0.98 |
| group commutator | 0.074 | ≤ 6.9e−7 | 55–76 cm⁻¹ | 0.73–1.0 |

**Free-ion scale.** The group commutator isolates the six lowest multiplets to 1 − F ≤ 7×10⁻⁷ and the ground multiplet to 10⁻⁶ in 26 iterations. The higher manifolds, which are closer together, are still mixed after 150 iterations (minimum over all 35 manifolds of 1 − F ≈ 0.15).

**Crystal-field scale.** No method resolves the crystal-field levels inside ⁴I₁₅/₂, even where the multiplet itself is isolated. The energy errors stay at 55–130 cm⁻¹, comparable to the whole 339 cm⁻¹ splitting, and the individual level fidelities stay near 0.

**Why.** The greedy objective ‖σ(Ĥ_k)‖ is dominated by the Coulomb and spin–orbit couplings, which are 10²–10³ times larger than the crystal-field couplings. Every step is spent on them, while the convergence rate for a coupling across a gap g scales with g². This is the multiscale effect the proposal asked about.

## 4. Two-stage strategy: DBQA inside the isolated ground multiplet (stage 2)

**Switch.** The group-commutator stage 1 is stopped at the first stored iteration (30) at which the ground manifold is isolated to 1 − F < 10⁻⁶.

**Block.** The 16 basis states with the lowest diagonal entries of Ĥ₃₀ then span ⁴I₁₅/₂. Their 16×16 block is a 4-qubit Hamiltonian, the same size as Model 1. Its eigenvalues differ from the exact ground-multiplet levels of the full Hamiltonian by at most 0.018 cm⁻¹; this is the error the restriction itself introduces.

**Runs.** The Model 1 strategies ran on the block for 100 iterations: canonical, Pauli-Z with all 15 Z products, and group commutator. The block was in units of 100 cm⁻¹, so its norm (1.8) is comparable to Model 1's and the Model 1 step grids apply. Errors are measured against the block eigenvalues.

| Stage 2 method | max \|ΔE\| < 1 cm⁻¹ | < 0.1 cm⁻¹ | max 1 − F < 10⁻² | After 100: max \|ΔE\| | After 100: max 1 − F |
|---|---|---|---|---|---|
| canonical | 76 | > 100 | > 100 | 0.39 cm⁻¹ | 6.6e−2 |
| **Pauli-Z** | **10** | **45** | **53** | **2.9e−3 cm⁻¹** | **4.8e−4** |
| group commutator | 27 | > 100 | > 100 | 0.31 cm⁻¹ | 5.2e−2 |

**The two-stage strategy resolves the crystal field.** After 30 full iterations, Pauli-Z generators on the block bring every ⁴I₁₅/₂ level within 1 cm⁻¹ in 10 more iterations and within 0.1 cm⁻¹ in 45. In 150 full-Hamiltonian iterations, no method got below 55 cm⁻¹.

**The slowest mode is a symmetry-allowed coupling, not the closest pair.** The residual errors come in equal and opposite pairs: after 100 Pauli-Z steps, +0.0029 cm⁻¹ on level 0 and −0.0029 cm⁻¹ on level 2 (5.98 cm⁻¹), and likewise for levels 1 and 3 (44.9 cm⁻¹). This is the second-order signature of a remaining coupling between two levels. In D₂d the crystal field couples only levels of the same irrep, so levels 0 and 2 share an irrep, and so do levels 1 and 3. This is consistent with the paper's assignment Γ₇, Γ₆, Γ₇, Γ₆ for the four lowest levels; our model has the same level ordering as the paper's calculation. The closest pair (0 and 0.36 cm⁻¹) has different symmetry, never mixes, and is not the bottleneck.

**The other block methods.** The canonical and group-commutator rotations leave the 0↔2 and 1↔3 couplings largely unresolved after 100 iterations (0.3–0.4 cm⁻¹).

## 5. Compilation cost (Example 6)

**Pauli content.** The embedded Model 2 Hamiltonian has 93 184 non-zero Pauli strings out of 4⁹ = 262 144. Most have weight 6–8:

| Weight | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|---|---|
| Strings | 1 | 17 | 156 | 926 | 3794 | 10804 | 21204 | 27530 | 21309 | 7443 |

**CNOT counts per Trotter step** (CNOT-ladder compilation, no gate cancellation):

| Hamiltonian | Qubits | Pauli terms | CNOTs per Trotter step |
|---|---|---|---|
| **Er:LiYF₄ (Model 2)** | 9 | 93 184 | **1 073 026** (1st order), 2 146 036 (2nd order) |
| Ce:YAG (Model 1) | 4 | 71 | 272 (1st order) |
| XXZ chain, `XXZ_compilation_line` | 9 | 24 | 36 (2nd order) |
| XXZ chain | 4 | 9 | 12 (2nd order) |

**Comparison with the spin chains.** A second-order Trotter step of the rare-earth Hamiltonian costs about 6×10⁴ times more CNOTs than the second-order 9-qubit XXZ step used in the Fig. 4 examples.

**Why it is dense.** In the |SLJM⟩ basis with a binary state-to-qubit map, the free-ion operators couple states whose qubit labels differ in many bits. A different state-to-qubit map might reduce the cost; this was not tested.

**DB-DOI recursion.** On top of that, the circuit-level DBQA recursion triples the number of oracle queries per step (Model 1, §4.2). With Trotterised evolution in this basis, circuit-level DBQA for Model 2 is far beyond near-term hardware.

## 6. Answers

**How does the multiscale structure affect DBQA convergence?**
- **Greedy DBQA works scale by scale.** It isolates the free-ion multiplets quickly: the ground multiplet to 1 − F < 10⁻⁶ in 26 group-commutator iterations. It then stalls on the crystal-field scale; after 150 iterations the ⁴I₁₅/₂ levels are still 55–130 cm⁻¹ off.
- **A multiscale-aware two-stage strategy fixes this.** Isolate the multiplet, then run DBQA in the 4-qubit multiplet block. Every ⁴I₁₅/₂ level then comes within 0.1 cm⁻¹ in 30 + 45 iterations.
- **Symmetry sets the remaining convergence.** The slowest mode is the coupling between same-irrep levels (Γ₇ 0 ↔ 5.98 cm⁻¹), not the smallest gap.

**How costly is the dense rare-earth Hamiltonian to compile?** About 2×10⁶ CNOTs per second-order Trotter step at 9 qubits (10⁶ at first order), against 36 per second-order step for the 9-qubit XXZ chain. 93 184 of the 262 144 possible Pauli strings are present, mostly of weight 6–8.

**Implications for Model 3.** The ¹⁶⁷Er:Y₂SiO₅ hyperfine model has a much smaller energy scale under much larger ones. Two things follow:
- the two-stage, multiplet-restricted strategy should be the default;
- tolerances should be set relative to the gaps between same-symmetry levels, not to the smallest gap overall.

## 7. Limitations

- **Free-ion parameters.** The Model 2 Hamiltonian uses refitted free-ion parameters (§2.3), not the paper's values, because the paper's E³ convention could not be reproduced. The crystal field is the paper's.
- **Scheduling.** Stage 1 uses greedy grid-search scheduling with 100 durations (60 for Pauli-Z). The Pauli-Z generators are restricted to single-qubit Zᵢ. Other schedulers, higher-weight Z products and gradient-optimised D were not tried at 9 qubits.
- **Stage-2 switch.** Stage 2 starts only from the group-commutator run; the canonical and Pauli-Z stage 1 runs do not isolate the ground multiplet to 10⁻⁶ within 150 iterations. The switch uses the stored checkpoints, every 10 iterations, so it happens at iteration 30 although isolation is reached at 26.
- **Not simulated.** Compiled circuits were not simulated for Model 2; only the gate counts were computed. No noise was simulated.

## Reproduction

```bash
export PYCF_ROOT=/path/to/pycf
cd examples/Crystal_field_PyCF
jupyter nbconvert --to notebook --execute --inplace Example_4_validation_Er_LiYF4.ipynb Example_5_DBQA_Er_LiYF4_generate_data.ipynb Example_6_DBQA_Er_LiYF4_results.ipynb
```

Run times on this laptop are about 20 minutes for Example 4 (the fits), 35 minutes for Example 5 and under a minute for Example 6. Example 5 also writes the iterated matrices to `simulation_results/Er_LiYF4_DBQA_checkpoints.npz` (17 MB), which is not under version control.
