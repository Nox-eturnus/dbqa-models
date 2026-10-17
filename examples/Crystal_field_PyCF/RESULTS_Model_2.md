# Model 2: DBQA for Er³⁺:LiYF₄ (full 4f¹¹): results

All numbers below come from the executed notebooks in this folder:
- `Example_4_validation_Er_LiYF4.ipynb`: Hamiltonian and validation
- `Example_5_DBQA_Er_LiYF4_generate_data.ipynb`: DBQA runs
- `Example_6_DBQA_Er_LiYF4_results.ipynb`: analysis and compilation cost
- `Example_7_compiled_DBQA_Er_LiYF4.ipynb`: compiled oracle, Trotter error and circuit-level DBQA (DB-DOI)

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

**Free-ion scale.** The group commutator isolates the six lowest multiplets to 1 − F ≤ 7×10⁻⁷. The ground multiplet first reaches 1 − F < 10⁻⁶ at iteration 16 and stays there from iteration 26. The higher manifolds, which are closer together, are still mixed after 150 iterations (minimum over all 35 manifolds of 1 − F ≈ 0.15).

**Crystal-field scale.** No method resolves the crystal-field levels inside ⁴I₁₅/₂, even where the multiplet itself is isolated. The energy errors stay at 55–130 cm⁻¹, comparable to the whole 339 cm⁻¹ splitting, and the individual level fidelities stay near 0.

**Why (corrected 2026-10-17, see §8).** We first attributed this to the greedy objective ‖σ(Ĥ_k)‖ being dominated by the Coulomb and spin–orbit couplings, which are 10²–10³ times larger than the crystal-field couplings. The direct cause is the duration grid. A group-commutator step removes a coupling across a gap g (in the rescaled units) only if its duration reaches about 1/g, and the crystal-field gaps of ⁴I₁₅/₂ are 2×10⁻⁵–10⁻² in units of 20 000 cm⁻¹, while the grid stopped at s = 10. With durations up to 10⁴ (`Example_11`), the same group-commutator stage 1 brings every ⁴I₁₅/₂ level within 0.19 cm⁻¹ in 150 iterations. The multiscale structure does set how long the steps must be and how many are needed, but it does not stop greedy DBQA.

## 4. Two-stage strategy: DBQA inside the isolated ground multiplet (stage 2)

**Switch.** Stage 1 stores the iterated Hamiltonian at the first iteration at which the ground manifold reaches 1 − F < ε, for ε = 10⁻³ and ε = 10⁻⁶. Stage 2 starts from every stage-1 method that reaches ε within 150 iterations. The group commutator first reaches 10⁻⁶ at iteration 16. The table in §3.1 counts 26 because it requires the infidelity to stay below the threshold from then on.

**Block.** At the switch, the 16 basis states with the lowest diagonal entries of Ĥ_k span ⁴I₁₅/₂. Their 16×16 block is a 4-qubit Hamiltonian, the same size as Model 1. The **restriction error** is the largest difference between the block eigenvalues and the exact ground-multiplet levels of the full Hamiltonian. It is the price of switching at that point.

**Runs.** The Model 1 strategies ran on the block for 100 iterations: canonical, Pauli-Z with all 15 Z products, and group commutator. The block was in units of 100 cm⁻¹, so its norm (1.8) is comparable to Model 1's and the Model 1 step grids apply.

**Errors.** Stage-2 errors are measured against the block eigenvalues. The error against the exact levels is at most the stage-2 error plus the restriction error, shown as "vs exact ≤".

| Stage 1 | ε | Switch k | Restriction error | Stage 2 | its to max\|ΔE\| < 1 cm⁻¹ | < 0.1 cm⁻¹ | max 1 − F < 10⁻² | Final max\|ΔE\| (block) | vs exact ≤ |
|---|---|---|---|---|---|---|---|---|---|
| canonical | 10⁻³ | 1 | 12 cm⁻¹ | canonical | 74 | > 100 | > 100 | 0.17 cm⁻¹ | 13 cm⁻¹ |
| | | | | Pauli-Z | 11 | 31 | 67 | 2.4e−3 cm⁻¹ | 12 cm⁻¹ |
| | | | | group commutator | 27 | > 100 | > 100 | 0.13 cm⁻¹ | 13 cm⁻¹ |
| canonical | 10⁻⁶ | not reached in 150 | | | | | | | |
| Pauli-Z | 10⁻³ | 24 | 46 cm⁻¹ | canonical | 98 | > 100 | > 100 | 0.91 cm⁻¹ | 47 cm⁻¹ |
| | | | | Pauli-Z | 12 | 52 | 63 | 6.5e−3 cm⁻¹ | 46 cm⁻¹ |
| | | | | group commutator | 31 | > 100 | > 100 | 0.38 cm⁻¹ | 46 cm⁻¹ |
| Pauli-Z | 10⁻⁶ | not reached in 150 | | | | | | | |
| group commutator | 10⁻³ | 1 | 7.5 cm⁻¹ | canonical | 80 | > 100 | > 100 | 0.48 cm⁻¹ | 7.9 cm⁻¹ |
| | | | | Pauli-Z | 11 | 51 | 62 | 6.2e−3 cm⁻¹ | 7.5 cm⁻¹ |
| | | | | group commutator | 29 | > 100 | > 100 | 0.41 cm⁻¹ | 7.9 cm⁻¹ |
| **group commutator** | **10⁻⁶** | **16** | **0.072 cm⁻¹** | canonical | 76 | > 100 | > 100 | 0.39 cm⁻¹ | 0.46 cm⁻¹ |
| | | | | **Pauli-Z** | **10** | **45** | **53** | **2.9e−3 cm⁻¹** | **0.075 cm⁻¹** |
| | | | | group commutator | 27 | > 100 | > 100 | 0.31 cm⁻¹ | 0.38 cm⁻¹ |

**Stage 2 always works on the block.** From every switching point, Pauli-Z generators bring the block levels within 1 cm⁻¹ in 10–12 iterations and within 0.1 cm⁻¹ in 31–52. The canonical and group-commutator rotations need 74–98 and 27–31 iterations to reach 1 cm⁻¹.

**How well stage 1 isolates the multiplet sets the final accuracy.** Switching at 1 − F < 10⁻³ is quick (1–24 iterations), but the block is then 7.5–46 cm⁻¹ off the true levels, so the final error is dominated by the restriction. Only the group commutator reaches 10⁻⁶ (16 iterations), and there the two-stage strategy gives **every ⁴I₁₅/₂ level within 0.075 cm⁻¹ of the exact value after 16 + 45 iterations** (16 + 12 with the group commutator and longer durations, §8). In 150 full-Hamiltonian iterations no method got below 55 cm⁻¹. The canonical and Pauli-Z stage 1 runs never reach 10⁻⁶ within 150 iterations (they end at 1.4×10⁻⁵ and 4.1×10⁻⁴), so with them the two-stage strategy is limited to errors of order 10 cm⁻¹.

**The slowest mode is a symmetry-allowed coupling, not the closest pair** (on the original duration grids; with durations up to 10³ this coupling is resolved in about 12 iterations, §8). The residual errors come in equal and opposite pairs. After 100 Pauli-Z steps (group-commutator stage 1, ε = 10⁻⁶) they are +0.0029 cm⁻¹ on level 0 and −0.0029 cm⁻¹ on level 2 (5.98 cm⁻¹), and likewise, much smaller, for levels 1 and 3. This is the second-order signature of a remaining coupling between two levels. In D₂d the crystal field couples only levels of the same irrep, so levels 0 and 2 share an irrep, and so do levels 1 and 3. This is consistent with the paper's assignment Γ₇, Γ₆, Γ₇, Γ₆ for the four lowest levels; our model has the same level ordering as the paper's calculation. The closest pair (0 and 0.36 cm⁻¹) has different symmetry, never mixes, and is not the bottleneck.

**The other block methods.** The canonical and group-commutator rotations leave the 0↔2 coupling largely unresolved after 100 iterations (0.3–0.4 cm⁻¹).

## 5. Compilation cost and compiled circuits (Examples 6 and 7)

### 5.1 Pauli content and CNOTs per Trotter step (Example 6)

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

### 5.2 Compiled evolution oracle (Example 7)

**Method.** At 9 qubits a compiled Trotter step has about 10⁶ CNOTs, so the circuits are not built gate by gate; several million gate objects would not fit in memory. Two tools replace them:
- **Exact unitary.** `compiled_trotter_unitary` computes the exact unitary of the compiled gate sequence, term by term. Each CNOT ladder equals e^{−iθP} exactly, and the tests check the whole sequence against the explicit circuits at 3 qubits.
- **Exact depth.** `compiled_trotter_depth` computes the exact depth by scheduling the same gate sequence. It is checked against `Circuit.depth` in the tests.

**Trotter error.** The table gives the error against e^{−it(Ĥ − c_I)} for the longest oracle call DB-DOI makes (t = 1, i.e. s ≤ 1).

| Order | Steps n | δt | Error at t = 1 | CNOTs per call | Depth per call |
|---|---|---|---|---|---|
| 1 | 1 | 1.0 | 0.97 | 1 073 026 | 1 323 115 |
| 1 | 4 | 0.25 | 0.26 | 4 292 104 | 5 292 439 |
| 1 | 16 | 0.0625 | 6.4e−2 | 17 168 416 | 21 169 735 |
| 2 | 1 | 1.0 | 0.44 | 2 146 036 | 2 646 213 |
| 2 | 2 | 0.5 | 8.6e−2 | 4 292 072 | 5 292 426 |
| 2 | 4 | 0.25 | 2.2e−2 | 8 584 144 | 10 584 852 |
| 2 | 8 | 0.125 | 5.4e−3 | 17 168 288 | 21 169 704 |
| **2** | **16** | **0.0625** | **1.4e−3** | **34 336 576** | **42 339 408** |

**Choice.** The oracle is second order with 16 steps: the fewest steps that keep the error at or below 5×10⁻³, the criterion used for Model 1. Each call costs 3.4×10⁷ CNOTs at depth 4.2×10⁷. For Model 1 the same accuracy took 2152 CNOTs.

### 5.3 Circuit-level DBQA: DB-DOI ground-state preparation (Example 7)

**Symmetry sectors block DB-DOI.** In the |SLJM⟩ basis the D₂d Hamiltonian splits into four sectors of 91 states each, plus the 148 uncoupled padding states. The free ion conserves M_J and the crystal field changes it by 4. A magnetic-field D also preserves the sectors, so DB-DOI never leaves the sector of its reference state. The ground Kramers doublet lies in one pair of sectors (one partner in each).

The basis state with the lowest diagonal energy overall (|157⟩) lies in the other pair, whose lowest level is the 0.36 cm⁻¹ doublet; it has zero overlap with the ground doublet. A first run from it showed exactly this: the ground-doublet fidelity stayed at 0.

Model 1 does not have this problem, because each of its two sectors contains one partner of every Kramers doublet.

**Remedy.** We find the sectors from the sparsity pattern of Ĥ alone (connected components, no diagonalisation), merge Kramers-partner sectors, and run DB-DOI from the lowest-diagonal basis state of each sector.

**Run.** The diagonal operator was a magnetic field chosen at each step among ±Σᵢ Zᵢ and ±Zᵢ, with the step on a grid s ∈ [10⁻³, 1]; all chosen steps are inside the grid. The run used the ideal oracle, and the CNOT counts follow from the number of queries.

| Reference | Sector's lowest level | k | E_k − E₀ | 1 − F (sector's lowest level) | 1 − F (ground doublet) | 1 − F (⁴I₁₅/₂) | Queries | CNOTs |
|---|---|---|---|---|---|---|---|---|
| \|238⟩ | ground doublet | 0 | 907.0 cm⁻¹ | 0.951 | 0.951 | 2.9e−2 | 0 | 0 |
| | | 1 | 121.2 cm⁻¹ | 0.950 | 0.950 | 1.2e−3 | 1 | 3.4e7 |
| | | 3 | 103.3 cm⁻¹ | 0.949 | 0.949 | 4.4e−4 | 13 | 4.5e8 |
| | | 5 | 97.7 cm⁻¹ | 0.948 | 0.948 | 2.1e−4 | 121 | 4.2e9 |
| \|157⟩ | 0.36 cm⁻¹ doublet | 0 | 893.4 cm⁻¹ | 0.768 | 1 | 2.9e−2 | 0 | 0 |
| | | 1 | 113.5 cm⁻¹ | 0.761 | 1 | 1.6e−3 | 1 | 3.4e7 |
| | | 3 | 89.7 cm⁻¹ | 0.759 | 1 | 4.9e−4 | 13 | 4.5e8 |
| | | 5 | 86.2 cm⁻¹ | 0.757 | 1 | 2.6e−4 | 121 | 4.2e9 |

**Depths.** These are upper bounds: the depth of the ideal-oracle circuit with every oracle gate replaced by the oracle depth. At k = 5 the bound is 5.1×10⁹.

**The circuit level shows the same multiscale behaviour as the matrix level.**
- **Free-ion scale.** DB-DOI removes the free-ion part of the energy quickly: 907 → 121 cm⁻¹ in one step, and the state's leakage out of the ground multiplet falls from 3×10⁻² to 2×10⁻⁴.
- **Crystal-field scale.** It does not resolve the crystal-field levels: after 5 steps the state is still spread over ⁴I₁₅/₂, with only 5 % weight on the ground doublet, 98 cm⁻¹ above it.
- **Cost.** These 5 steps already need 121 oracle calls, 4.2×10⁹ CNOTs.

**Compiled replay.** Replaying the first two steps from |238⟩ with the compiled oracle (`CompiledUnitaryHamiltonian`) gives 120.92 and 107.52 cm⁻¹, against 121.21 and 107.55 cm⁻¹ with the ideal oracle. The ground-doublet infidelities agree to four digits. The Trotter error of the chosen compilation is therefore negligible for this algorithm, as for Model 1.

**DB-DOI recursion.** The number of oracle calls triples per step. With Trotterised evolution in this basis, circuit-level DBQA for Model 2 is far beyond near-term hardware.

## 6. Answers

**How does the multiscale structure affect DBQA convergence?**
- **Greedy DBQA works scale by scale.** It isolates the free-ion multiplets quickly: the ground multiplet first reaches 1 − F < 10⁻⁶ after 16 group-commutator iterations (and stays there from 26). With durations s ≤ 10 it then stalls on the crystal-field scale; after 150 iterations the ⁴I₁₅/₂ levels are still 55–130 cm⁻¹ off. With s ≤ 10⁴ it does not stall: the group commutator reaches 0.19 cm⁻¹ for every ⁴I₁₅/₂ level in 150 iterations (§8).
- **A multiscale-aware two-stage strategy fixes this.** Isolate the multiplet well (1 − F < 10⁻⁶, reached only by the group commutator, at iteration 16), then run DBQA in the 4-qubit multiplet block. Every ⁴I₁₅/₂ level then comes within 0.075 cm⁻¹ of the exact value after 16 + 45 iterations (16 + 12 with longer durations, §8). Switching at a looser isolation (10⁻³) leaves 7–46 cm⁻¹ errors from the restriction alone.
- **Symmetry sets the remaining convergence on short duration grids.** The slowest mode is the coupling between same-irrep levels (Γ₇ 0 ↔ 5.98 cm⁻¹), not the smallest gap; with durations up to 10³ it is no longer slow (§8).

**How costly is the dense rare-earth Hamiltonian to compile?**
- **Per Trotter step.** About 2×10⁶ CNOTs per second-order step at 9 qubits (10⁶ at first order), against 36 for the 9-qubit XXZ chain. 93 184 of the 262 144 possible Pauli strings are present, mostly of weight 6–8.
- **Per oracle call.** An oracle accurate to 1.4×10⁻³ for the durations DB-DOI uses costs 3.4×10⁷ CNOTs at depth 4.2×10⁷ per call.
- **For DB-DOI.** Five DB-DOI steps cost 4.2×10⁹ CNOTs. They take the state into the ground multiplet but do not resolve its crystal-field levels.
- **A trap for state preparation.** Circuit-level state preparation must start in the right symmetry sector: from the lowest-diagonal basis state overall it can never reach the ground doublet.

**Implications for Model 3.** The ¹⁶⁷Er:Y₂SiO₅ hyperfine model has a much smaller energy scale under much larger ones. Two things follow:
- the two-stage, multiplet-restricted strategy should be the default, and the quality of the stage-1 isolation should be reported, since it bounds the final accuracy;
- reference states for circuit-level preparation must be chosen per symmetry sector;
- tolerances should be set relative to the gaps between same-symmetry levels, not to the smallest gap overall.

## 7. Limitations

- **Free-ion parameters.** The Model 2 Hamiltonian uses refitted free-ion parameters (§2.3), not the paper's values, because the paper's E³ convention could not be reproduced. The crystal field is the paper's.
- **Scheduling.** Stage 1 uses greedy grid-search scheduling with 100 durations (60 for Pauli-Z). The Pauli-Z generators are restricted to single-qubit Zᵢ. Other schedulers, higher-weight Z products and gradient-optimised D were not tried at 9 qubits.
- **Stage-2 switch.** The canonical and Pauli-Z stage 1 runs do not isolate the ground multiplet to 10⁻⁶ within 150 iterations, so for them stage 2 is only available at the 10⁻³ switch. Running stage 1 longer, or with other generators, might reach it.
- **Compiled circuits.** These are simulated through the exact unitary of the compiled gate sequence, not gate by gate. Depths are exact per oracle call and upper bounds for the full DB-DOI circuits.
- **DB-DOI scope.** DB-DOI was run for 5 steps (121 oracle calls) with magnetic-field generators. Two of those steps were replayed with the compiled oracle. A two-stage, multiplet-restricted variant was added later (`Example_14`): inside the ⁴I₁₅/₂ block, from a reference in the ground-doublet sector, the ground-doublet infidelity falls from 0.95 to 0.31 in 8 steps (3 280 oracle calls, 4.0×10⁷ CNOTs).
- **Noise.** None in §5; `Example_15` simulates the two-stage block circuits with depolarizing noise after every CNOT: their signal is lost already at p = 10⁻⁵, and keeping half of it at k = 8 needs p ≲ 2×10⁻⁸ (`RESULTS_Model_3.md` §4.2).

## 8. Addendum (2026-10-17): the duration grid (Examples 11 and 12)

Model 3 showed that the duration grid has to reach the inverse of the smallest gap to be resolved, in the units the Hamiltonian is rescaled to. `Example_11_step_range_Er_LiYF4.ipynb` repeats the group-commutator stage 1 of Example 5 with 160 durations in [10⁻⁴, 10⁴] instead of 100 in [10⁻⁴, 10] (47 minutes):

| Durations | max \|ΔE\|, ⁴I₁₅/₂ levels, k = 25 | k = 50 | k = 100 | k = 150 | max 1 − F (levels), k = 150 | 1 − F (⁴I₁₅/₂ manifold), k = 150 |
|---|---|---|---|---|---|---|
| s ≤ 10 (Example 5) | 117 cm⁻¹ | 115 | 98 | 76 | 1.0 | 3×10⁻¹⁰ |
| s ≤ 10⁴ | 117 cm⁻¹ | 37 | 11 | **0.19** | 3×10⁻² | 9×10⁻¹¹ |

Forty of the 150 chosen durations exceed 10 (median 1.7, maximum 2×10³). So greedy DBQA on the full 9-qubit Hamiltonian does resolve the crystal-field levels when the grid allows long enough steps. The two-stage strategy remains faster and more accurate: 0.075 cm⁻¹ after 16 + 45 iterations, against 0.19 cm⁻¹ after 150. The remaining infidelity at k = 150 (3 %, energy errors 0.19 cm⁻¹) sits on levels 0 and 2 (0 and 5.98 cm⁻¹), the same-irrep coupling identified in §4 as the slowest mode; the closest pair (0 and 0.36 cm⁻¹) is resolved to 8×10⁻⁴. This revises the explanation in §3.2 and the first answer in §6.

**Stage 2 and DB-DOI (Example 12).** The same criterion applies inside the 16-state block (units of 100 cm⁻¹): the slowest coupling, 0 ↔ 5.98 cm⁻¹, is 0.06, so it needs durations of about 17, while §4 used grids ending at 0.3 and 1. `Example_12_duration_grids_Ce_YAG_Er_LiYF4.ipynb` recomputes the stage-1 switch (iteration 16, 1 − F < 10⁻⁶) and reruns stage 2 on the original grids, which reproduces §4, and with 300 log-spaced durations in [10⁻⁴, 10³]. Errors are against the block eigenvalues:

| Grid | Method | max \|ΔE\| < 1 cm⁻¹ | < 0.1 | < 0.01 | after 100 iterations | max 1 − F after 100 |
|---|---|---|---|---|---|---|
| original | canonical | 76 | > 100 | > 100 | 0.39 cm⁻¹ | 6.5×10⁻² |
| original | Pauli-Z | 10 | 45 | 82 | 2.9×10⁻³ cm⁻¹ | 4.9×10⁻⁴ |
| original | group commutator | 27 | > 100 | > 100 | 0.31 cm⁻¹ | 5.2×10⁻² |
| ≤ 10³ | canonical | 25 | > 100 | > 100 | 0.18 cm⁻¹ | 3.0×10⁻² |
| ≤ 10³ | Pauli-Z | 7 | 13 | 25 | 4×10⁻¹¹ cm⁻¹ | 5×10⁻¹² |
| ≤ 10³ | **group commutator** | **7** | **12** | **12** | **2×10⁻¹⁰ cm⁻¹** | **3×10⁻¹⁴** |

So with long enough durations the group commutator also wins inside the block, and the same-irrep coupling of §4 stops being slow: every ⁴I₁₅/₂ level is within the restriction error (0.072 cm⁻¹) of the exact value after **16 + 12** iterations, against 16 + 45 with Pauli-Z on the original grid. The statements in §4 and §6 about which generator wins in stage 2 and about the slowest mode hold for the original grids only.

DB-DOI from |238⟩ (§5.3, ideal oracle) with 40 log-spaced durations in [10⁻³, 10²] instead of 20 in [10⁻³, 1]: after 5 steps the energy is 72 cm⁻¹ above the ground state (98 cm⁻¹ on the original grid), but the ground-doublet infidelity is still 0.94 (0.95). Two of the five chosen durations sit at the edge of the grid (10²), so the DB-DOI result is still grid-limited; the crystal-field scale is not resolved at circuit level within 5 steps either way.

## Reproduction

```bash
export PYCF_ROOT=/path/to/pycf
cd examples/Crystal_field_PyCF
jupyter nbconvert --to notebook --execute --inplace Example_4_validation_Er_LiYF4.ipynb Example_5_DBQA_Er_LiYF4_generate_data.ipynb Example_6_DBQA_Er_LiYF4_results.ipynb Example_7_compiled_DBQA_Er_LiYF4.ipynb Example_11_step_range_Er_LiYF4.ipynb Example_12_duration_grids_Ce_YAG_Er_LiYF4.ipynb
```

Approximate run times on this laptop, each running alone:

| Notebook | Run time |
|---|---|
| Example 4 (the fits) | about 20 minutes |
| Example 5 | about 35 minutes |
| Example 6 | under a minute |
| Example 7 | about an hour, of which about 20 minutes compute the compiled Trotter steps |
| Example 11 (addendum, §8) | about 50 minutes |
| Example 12 (addendum, §8; also Model 1) | about 35 minutes |

Two outputs are large and not under version control:
- `simulation_results/Er_LiYF4_DBQA_checkpoints.npz` (17 MB), the iterated matrices written by Example 5;
- `simulation_results/Er_LiYF4_trotter_steps.npz` (33 MB), the cached compiled Trotter steps of Example 7.
