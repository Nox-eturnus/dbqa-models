# Model 3: DBQA for ¹⁶⁷Er³⁺:Y₂SiO₅ (site 2) with hyperfine structure, over a field sweep: results

All numbers below come from the executed notebooks in this folder:
- `Example_8_validation_Er_YSO.ipynb`: Hamiltonian, validation, ZEFOZ points, choice of the field sweep
- `Example_9_DBQA_Er_YSO_generate_data.ipynb`: DBQA runs over the sweep
- `Example_10_DBQA_Er_YSO_results.ipynb`: analysis and compilation cost
- `Example_13_avoided_crossings_Er_YSO.ipynb`: avoided crossings and DBQA across them
- `Example_14_compiled_DBQA_Er_YSO.ipynb`: compiled oracle, DB-DOI and two-stage DB-DOI
- `Example_15_noisy_DB_DOI.ipynb`: DB-DOI under depolarizing noise (also Models 1 and 2)

The data are in `simulation_results/` and the figures in `pdfs_figures/`. Every run is deterministic, because grid-search scheduling, the generator selection and the ZEFOZ search involve no random sampling.

## 1. Model

**System.** Er³⁺ at crystallographic site 2 of Y₂SiO₅ (C₁ point symmetry), isotope ¹⁶⁷Er with nuclear spin I = 7/2, in a magnetic field **B**. Axes x, y, z are the optical extinction axes D1, D2 and the crystal b axis.

**Source.** N. L. Jobbitt, J.-P. R. Wells, M. F. Reid, S. P. Horvath, P. Goldner and A. Ferrier, *Prediction of the optical polarization and high field hyperfine structure via a parametrized crystal-field model for the low symmetry centers in Er³⁺ doped Y₂SiO₅*, Phys. Rev. B 104, 155121 (2021), arXiv:2206.09080. We use their Table III ("This study", site 2) and the free-ion parameters they fixed (Table IV, Carnall's Er³⁺:LaF₃ values). M², M⁴, P⁴ and P⁶ follow from M⁰ and P² with the ratios of PyCF's Er:YSO example (0.56, 0.31, 0.50, 0.10). The constants are in `utils_crystal_field.py` (`ER_YSO_SITE2_JOBBITT2021`, `ER_YSO_SITE2_HYPERFINE_JOBBITT2021`).

**Hamiltonian.**
H(**B**) = H_FI + Σ B^k_q C^(k)_q + a_l H_HF + a_Q H_Q + μ_B **B**·(**L** + g_s **S**) − g_n μ_N **B**·**I**,
with all 27 crystal-field parameters independent (complex B^k_q for q > 0), a_l = 0.005389 cm⁻¹ and a_Q = 0.0240 cm⁻¹. The nuclear Zeeman term uses g_n = −0.1618 (Matsuura et al., Table I).

**Basis.** PyCF's `erhfs` basis: the intermediate-coupled ⁴I₁₅/₂ and ⁴I₁₃/₂ multiplets times |m_I⟩, (16 + 14) × 8 = 240 states, as in the paper's hyperfine fit. They are embedded in 256 dimensions (8 qubits) with 16 uncoupled padding states. `hyperfine_field_model` returns H₀ and the field operators, so that H(**B**) = H₀ + Σ_a B_a V_a.

**Energy scales.** 6 500 cm⁻¹ between the two multiplets; 27–140 cm⁻¹ between crystal-field levels; about 0.2 cm⁻¹ (6 GHz) of hyperfine structure inside each level; hyperfine gaps of 0.03 cm⁻¹ (800 MHz) and less.

## 2. Validation of the Hamiltonian (Example 8)

### 2.1 Transcription

The crystal-field strengths S_k = [(|B^k_0|² + 2Σ_{q>0}|B^k_q|²)/(2k+1)]^{1/2} computed from our transcription are 363.1, 653.3 and 151.5 cm⁻¹, exactly the values of Table III. Table II (59 levels, 57 measured energies, and 25, 24 and 21 measured g values along D1, D2 and b) is transcribed in `data/Er_YSO_site2_Jobbitt2021_Table2.csv`; the counts match those stated in the paper.

### 2.2 The full 4f¹¹ configuration (PyCF `f11cf`) against Table II

| Quantity | Result |
|---|---|
| paper fit, E_calc − E_exp (57 levels) | rms 12.7 cm⁻¹, max 26 cm⁻¹ |
| PyCF − paper E_calc (59 levels) | rms 21.9 cm⁻¹, max 37.7 cm⁻¹ |
| the same after removing one shift per multiplet | **rms 0.59 cm⁻¹** |
| g values, PyCF − paper calc (59 levels × 3 axes) | **rms 0.07, max 0.22** |
| Z1 g values (D1, D2, b) | PyCF 14.57, 2.13, 3.33; paper 14.65, 2.21, 3.13; experiment 15.31, 2.66, 2.84 |

**The multiplets are shifted.** The crystal-field splittings inside every multiplet agree with the paper to 0.6 cm⁻¹ rms, but the multiplets are displaced rigidly by +33 (⁴I₁₅/₂), +8 (⁴I₁₃/₂), +26, +37, +4, −4, −25, −4, −4, −5, −25 and −15 cm⁻¹ (up to ⁴G₁₁/₂).

**What it is not.**
- A refit of E_avg, F², F⁴, F⁶ and ζ to the paper's calculated levels leaves rms 12.1 cm⁻¹; adding α, β, γ, M⁰ or P² one at a time leaves 9.3–12.0 cm⁻¹. A refit of all 16 free-ion parameters reaches 0.8 cm⁻¹ only with unphysical values (F⁴ and F⁶ +40 %, γ = −8 900 cm⁻¹), so it explains nothing.
- Carnall's own M and P ratios (0.56, 0.38, 0.75, 0.50): rms 21.4 cm⁻¹.
- Truncating the crystal field to the 216 free-ion states below 40 000 cm⁻¹, as the paper says it did: rms 21.1 cm⁻¹.
- A site-2 transcription error: the site-1 parameters of the same table give the same multiplet pattern against Table I (+32, +7, +24, +34, …, −21 cm⁻¹).

**Conclusion.** Some convention of the free-ion operators in PyCF's bundled `f11cf` basis differs from the program behind Tables I and II; we could not identify which from the paper. The shifts are rigid, and the g values, which probe the crystal-field wave functions, agree to 0.07. Model 3 uses the two-multiplet `erhfs` basis, which reproduces the paper's hyperfine predictions exactly (§2.3).

### 2.3 The 240-state hyperfine Hamiltonian against Table VI

**Basis size.** The paper says its two-multiplet hyperfine basis has 224 states; PyCF's `erhfs` basis, which the authors' group distributes, has 240 = (16 + 14) × 8. The test suite checks its ordering on PyCF's state labels.

**Electronic levels.** Relative to Z1, the `erhfs` levels of ⁴I₁₅/₂ are 29.5, 67.8, 127.6, 171.0, 313.3, 346.8 and 409.7 cm⁻¹ (paper: 30, 67, 130, 171, 314, 348, 409); inside the two multiplets the splittings agree with the paper to 2.1 cm⁻¹ rms. Y1 lies 6 493 cm⁻¹ above Z1 (paper 6 516; the f11cf shift of §2.2).

**High-field hyperfine splittings (7 T along D1).**

| | Our model, no nuclear Zeeman | Table VI (Jobbitt et al.) | with nuclear Zeeman |
|---|---|---|---|
| Z1 splittings (MHz) | 897, 881, 865, 849, 833, 818, 802 | 897, 881, 865, 849, 833, 817, 801 | 888, 873, …, 793 |
| Y1 splittings (MHz) | 928, 912, 896, 880, 863, 847, 831 | 928, 912, 895, 879, 863, 847, 831 | 920, 904, …, 823 |

All 14 splittings agree to within 0.7 MHz, the rounding of the table. The paper's model therefore has no nuclear Zeeman term (with it every splitting drops by 8 MHz, g_n μ_N × 7 T). Model 3 keeps the term: it is physical and enters the ZEFOZ condition.

### 2.4 Zero field: the crystal-field model against the measured spin Hamiltonian

Z1 hyperfine levels at zero field (MHz above the lowest):

| Model | Levels |
|---|---|
| crystal-field model (Jobbitt et al.) | 0, 0, 840, 840, 1647, 1657, 2249, 2774, 3292, 3429, 4133, 4134, 4901, 4901, 5658, 5658 |
| spin Hamiltonian (S.-J. Wang et al., as used by Matsuura et al.) | 0, 0, 867, 872, 1669, 1727, 2404, 2609, 3148, 3525, 4108, 4198, 4931, 4941, 5725, 5725 |

The overall hyperfine scale agrees to about 1–3 %, but individual levels differ by up to 165 MHz. The pattern of near-degenerate pairs differs as well; it is set by the quadrupole tensor, which the crystal-field fit gets several times too small for site 2 (Jobbitt et al., Table V: their Q_bb = 1.69 MHz against 25.8 MHz measured). The paper's site-2 hyperfine fit used no excited-state data.

### 2.5 ZEFOZ points

**Method.** `zefoz_newton` solves ∇_B ν_ij = 0 by Newton's method with the exact gradient (Hellmann–Feynman) and exact Hessian (second-order perturbation theory) of `transition_derivatives`. `zefoz_search` runs it from many starting fields for every transition and merges ±**B** (time reversal).

**Matsuura et al. reproduced.** From the spin Hamiltonian of their Table II, our Newton search reproduces their five best site-2 points (their Table V): frequencies to 0.1 MHz, |S₁| below 10⁻¹² MHz/mT after refinement. Their text prints the optimal field as (∓378.9, ±73.2, ∓502.3) mT. That has the wrong sign of the b component relative to their polar angles (θ, φ) = (∓37.5383°, −10.9417°): with the Cartesian field as printed, ν₁₄,₁₅ = 783.7 MHz and |S₁| = 0.025 MHz/mT, so it is not a ZEFOZ point. The polar form reproduces their table.

**The crystal-field model at Matsuura et al.'s fields.** At their five fields the crystal-field model gives frequencies 9–62 MHz lower and gradients of 0.015–0.033 MHz/mT: none is a ZEFOZ point of Model 3.

**Search.** Like Matsuura et al. (their Appendix D), we search all 120 transitions of Z1 from many starting fields with |**B**| ≤ 3 T: their 6³ grid of components ±5, ±15, ±25 mT, plus 26 cube directions at 0.1, 0.3, 0.6, 1 and 2 T (346 starts), Newton capped at 60 iterations. Points are ranked by the largest curvature ‖S₂‖, which limits the coherence time once S₁ = 0 (their Eq. (5)).

| Hamiltonian | ZEFOZ points found | Five with the smallest ‖S₂‖ |
|---|---|---|
| spin Hamiltonian (Wang et al.) | 396 | \|14⟩–\|15⟩ 633.5 mT, \|6⟩–\|7⟩ 668.5 mT, \|13⟩–\|14⟩ 413.7 mT, \|13⟩–\|15⟩ 510.5 mT, \|12⟩–\|13⟩ 283.9 mT; all at θ ≈ −38° to −42°, ‖S₂‖ = 1.7–5.0×10⁻⁴ MHz/mT² |
| crystal-field model (Jobbitt et al.) | 258 | \|6⟩–\|7⟩ 149.9 mT, \|14⟩–\|15⟩ 106.1 mT, \|13⟩–\|14⟩ 100.1 mT, \|12⟩–\|13⟩ 92.9 mT, \|5⟩–\|6⟩ 103.2 mT; all at θ ≈ −77° to −82°, ‖S₂‖ = 5.3–11×10⁻⁴ MHz/mT² |

**The search reproduces Matsuura et al.** On their spin Hamiltonian it returns exactly the five transitions of their Table V, in the same order, at their fields. (With Newton capped at 40 iterations it missed the second one, so the cap matters.)

**The crystal-field model predicts a different ZEFOZ landscape.** Its best points are at 90–150 mT, near the D1–D2 plane (θ ≈ −80°), against 280–670 mT at θ ≈ −40° for the measured spin Hamiltonian. For each of Matsuura et al.'s five transitions, the nearest crystal-field ZEFOZ point of the same transition lies 215–571 mT away, at frequencies 7–53 MHz lower, with 1.8–3.4 times larger curvature. The likely cause is the quadrupole tensor (§2.4). **Jobbitt et al.'s site-2 parameters reproduce the high-field hyperfine splittings they were tested on (Table VI) but are not accurate enough to predict ZEFOZ points.** This is a statement about the parameters; a crystal-field fit including the measured site-2 A and Q tensors would be needed.

### 2.6 The field sweep for DBQA

DBQA is benchmarked on the crystal-field model itself, against the exact spectrum of the matrix it is given. The sweep was fixed by a rule set before the search: the crystal-field ZEFOZ point with the smallest ‖S₂‖ among those with 0.1 ≤ |**B**| ≤ 1 T. That is the transition |6⟩–|7⟩ (766.8 MHz) at **B*** = (−148.2, 8.3, 20.8) mT, |**B***| = 149.9 mT. The field is varied along its direction over 21 values in [0.8, 1.2]|**B***| (119.9–179.8 mT).

Along the sweep the curvature of ν₆,₇ is 1.4×10⁻⁴ MHz/mT², so ν₆,₇ varies by only **0.083 MHz** over 60 mT, while the hyperfine levels of Z1 span 37 GHz and the smallest gap inside Z1 is 762 MHz. Z1 stays at least 28.3 cm⁻¹ below Z2. Locating the turning point therefore needs transition frequencies accurate to well below 10⁻² MHz, i.e. 3×10⁻⁷ cm⁻¹ on a Hamiltonian whose spectrum spans 6 900 cm⁻¹ (a relative accuracy of 5×10⁻¹¹).

## 3. DBQA over the field sweep (Examples 9 and 10)

### 3.1 Strategy

**Stage A** (8 qubits, units of 5 000 cm⁻¹): group-commutator DBQA, D = Δ(Ĥ_k), until the 16 basis states with the lowest diagonal entries span Z1 to 1 − F < ε, for ε = 10⁻³, 10⁻⁶, 10⁻⁹ and 10⁻¹². **Stage B** (4 qubits): DBQA on the captured 16 × 16 block, shifted to zero trace and divided by its spectral norm, 60 iterations.

**Duration grids.** The Z1–Z2 gap is 5.4×10⁻³ in units of 5 000 cm⁻¹, so stage A uses 100 log-spaced durations in [10⁻³, 10²]. Stage B uses 120 in [10⁻³, 10³] (40 per sign for Pauli-Z). In a first test with s ≤ 10, stage A stalled at 1 − F = 0.15 after 300 iterations. Reaching durations of order 1/gap is what lets greedy DBQA cross the scales (this also revised Model 2, see `RESULTS_Model_2.md` §8).

**Three starts per field point (stage A).**
1. *cold*: from H(b);
2. *warm*: from U†_prev H(b) U_prev, with U_prev the accumulated stage-A rotation of the previous field point;
3. *replay*: the previous field point's schedule of durations applied to H(b), with D = Δ(Ĥ_k) recomputed, then grid search.

### 3.2 Stage A: isolating Z1

Mean iterations to 1 − F < ε over field points 2–21 (no run failed within 250 iterations):

| ε | cold | warm | replay (new steps + replayed steps) |
|---|---|---|---|
| 10⁻³ | 54 | 0 | 0 + 150 |
| 10⁻⁶ | 81 | 0 | 0 + 150 |
| 10⁻⁹ | 105 | 8 | 0 + 150 |
| 10⁻¹² | 150 | 63 | 0 + 150 |

- **Cold starts are almost field independent**: 54, 79–81 and 104–108 iterations to 10⁻³, 10⁻⁶ and 10⁻⁹ at every one of the 21 fields; 128–247 to 10⁻¹².
- **Warm starts** are already isolated to 10⁻⁶ at the next field (3 mT away) and need 7–9 iterations to 10⁻⁹.
- **The schedule transfers.** Replaying the previous field's durations on the new Hamiltonian reaches 10⁻¹² with no further step at every field. On a quantum computer this is the useful form: the warm start keeps every earlier step in the circuit (U_prev is built from evolutions under the previous Hamiltonians), whereas the replay has the circuit of a cold start and saves only the classical search for the durations.

### 3.3 Isolation sets the accuracy of the transition frequencies

Largest error of the 120 Z1 transition frequencies when the captured block is diagonalised exactly (all 21 fields, cold starts):

| ε | median | max |
|---|---|---|
| 10⁻³ | 300 MHz | 320 MHz |
| 10⁻⁶ | 0.28 MHz | 0.41 MHz |
| 10⁻⁹ | 1.5×10⁻⁴ MHz | 2.8×10⁻⁴ MHz |
| 10⁻¹² | 2.0×10⁻⁵ MHz | 4.2×10⁻⁵ MHz |

Down to ε = 10⁻⁹ the restriction error falls about in proportion to ε, as expected for a second-order shift (coupling²/gap ∝ (1 − F) × gap); below that it levels off near 2×10⁻⁵ MHz.

### 3.4 Stage B: resolving the hyperfine levels

Largest Z1 transition-frequency error from the DBQA diagonal at |**B**| = 149.9 mT:

| Block from | Method | after 20 iterations | after 60 | first below 0.01 MHz |
|---|---|---|---|---|
| 1 − F < 10⁻⁶ | group commutator | 0.35 MHz | 0.36 MHz | never (restriction floor 0.28 MHz) |
| 1 − F < 10⁻⁹ | group commutator | 0.042 MHz | **1.7×10⁻⁴ MHz** | 21 |
| 1 − F < 10⁻¹² | group commutator | 0.042 MHz | **2.6×10⁻⁵ MHz** | 21 |
| any | Pauli-Z | 570 MHz | 53–85 MHz | never |
| any | canonical | 4.7 GHz | 4.6 GHz | never |

- **The group commutator** reaches the restriction floor of §3.3 within 60 iterations.
- **The canonical bracket fails at a fixed point.** After stage A, the diagonal of the block holds two tight clusters (0–200 MHz and 12.7–12.9 GHz above the lowest entry), while its eigenvalues span 0–5.7 GHz and 31.7–37.3 GHz: the whole hyperfine structure sits in couplings between nearly degenerate diagonal entries. The canonical generator [Δ(Ĥ), Ĥ] scales those couplings by the diagonal differences, so after separating the two Zeeman branches (off-diagonal norm 3.12 → 0.39 in 3 steps) it stops moving. The group commutator, built from e^{isĤ}, is not affected.
- **Pauli-Z** with the 15 Z products as candidate D makes slow progress (53–85 MHz after 60 iterations).
- **Warm-started stage B** (group commutator, the previous field's block rotation applied first) needs 0–3 iterations to 0.01 MHz at fields 2–21, against 21–23 from a cold start.

### 3.5 Locating the ZEFOZ point from DBQA

A parabola fitted to ν₆,₇(b) over the 21 fields has its vertex at 152.21 mT for the exact frequencies. (The true turning point is at 149.86 mT; the 2.4 mT offset is the cubic term over a ±20 % window and is the same for every curve below.) From the DBQA estimates (group commutator, final stage-B iterate, cold stage A):

| Stage-A isolation | max error of ν₆,₇ over the sweep | vertex | error |
|---|---|---|---|
| 1 − F < 10⁻⁶ | 1.2×10⁻² MHz | 152.39 mT | +0.18 mT |
| 1 − F < 10⁻⁹ | 1.6×10⁻⁵ MHz | 152.20 mT | < 0.01 mT |
| 1 − F < 10⁻¹² | 2.3×10⁻⁵ MHz | 152.20 mT | < 0.01 mT |

DBQA locates the ZEFOZ turning point to 0.2 mT from 10⁻⁶ isolation and to better than 0.01 mT from 10⁻⁹, although the frequency changes by only 0.083 MHz over the sweep.

### 3.6 Avoided crossings (Example 13)

**No symmetry sectors.** At zero field and at **B*** the 240-state Hamiltonian is one connected block (10.8 % and 11.5 % of the matrix elements non-zero), so the reference-state trap of Model 2's DB-DOI cannot occur. `test_model3_single_sector` checks this.

**Where the crossings are.** `avoided_crossings` locates the interior minima of every gap between neighbouring Z1 levels along D1, D2 and b (801 fields in 0–200 mT, refined), for the crystal-field model and for the spin Hamiltonian of Wang et al. Matsuura et al. show avoided crossings only in figures (site 1), so both Hamiltonians are computed here.

| | Crystal-field model | Spin Hamiltonian |
|---|---|---|
| crossings found | 80 | 77 |
| matched (same level pair, field within 25 %) | 77 | 77 |
| crossing field, CF / SH | median 0.976, range 0.91–1.14 | |
| gap, CF / SH | median 0.087 | |
| gaps below 0.01 MHz | 11 | 0 |

The two models agree on **where** the crossings are: the fields are set by the g and A tensors, which the crystal-field fit reproduces. They disagree on **how wide** they are. Gaps come from the anisotropic couplings and the quadrupole term, and the crystal-field gaps are typically 10 times smaller, some below 10 kHz where the spin Hamiltonian has 0.4–19 MHz. Along D1 the first family of crossings sits at 3.7–4.2 mT in both models.

**DBQA across two crossings** (rule fixed in advance, along D1): the narrowest crystal-field crossing with a gap of at least 1 MHz, levels 5–6, 1.12 MHz at 3.949 mT; and the one closest to 100 MHz, levels 8–9, 93.6 MHz at 23.04 mT. Each was swept over 11 fields where the gap is below three times its minimum (3.931–3.967 mT and 21.72–24.34 mT). Stage A ran to 1 − F < 10⁻⁹; stage B ran 100 iterations with durations up to 10⁴, since the narrow gap is 10⁻⁴ of the block norm.

| | narrow (1.12 MHz) | wide (93.6 MHz) |
|---|---|---|
| stage A, cold start (first field) | 113 iterations | 106 |
| stage A, warm start (other fields) | 0 at every field | 0–4 |
| stage B (group commutator), iterations to resolve the pair (gap error < 10 % of the gap), cold | 33–64 | 15–43 |
| the same, warm-started from the previous field | 0–2 | 0–1 |
| at the centre, group commutator: gap error / pair 1 − F / individual 1 − F | 6×10⁻⁵ MHz / 10⁻¹⁴ / 2×10⁻⁷ | 1.4×10⁻⁵ MHz / 7×10⁻¹³ / 4×10⁻¹¹ |
| at the centre, canonical | not resolved (gap error 10 MHz, 900 %) | resolved at 92 (gap error 3 %) |
| at the centre, Pauli-Z | not resolved (56 %) | not resolved (27 %) |

At every field of both sweeps the final gap error of the group commutator is below 6×10⁻⁴ MHz. The pair-subspace fidelity converges much faster than the individual fidelities, so the pair is isolated before it is resolved internally, as in Model 1. The canonical bracket and Pauli-Z generators separate the pair from the rest (pair 1 − F between 3×10⁻² and 10⁻⁴) but do not resolve it internally.

## 4. Compilation cost (Example 10)

| Hamiltonian | Qubits | Pauli strings | CNOTs per first-order Trotter step |
|---|---|---|---|
| **¹⁶⁷Er:Y₂SiO₅, H(**B***), Model 3** | 8 | 5 978 of 65 536 | **55 530** (111 048 at second order) |
| its field-independent part H₀ | 8 | 5 718 | |
| field operators V_x, V_y, V_z | 8 | 256, 256, 143 | |
| stage-B block (Z1, captured at 10⁻¹²) | 4 | 251 | 1 018 |
| Er:LiYF₄, Model 2 | 9 | 93 184 | 1 073 026 |
| Ce:YAG, Model 1 | 4 | 71 | 272 |

By weight 0–8 the Model 3 strings number 1, 18, 107, 383, 820, 1 228, 1 524, 1 379 and 518. The two-multiplet basis makes Model 3 19 times cheaper per Trotter step than the full 4f¹¹ configuration of Model 2, despite the lower symmetry and the hyperfine terms. The stage-B block, which has no structure left after the stage-A rotation, is dense (251 of 256 strings).

### 4.1 Compiled oracle and circuit-level DBQA (Example 14)

**Oracle.** For H(**B***) in units of 5 000 cm⁻¹, the Trotter error at t = 1 (the longest evolution for durations s ≤ 1) is 2.2×10⁻² for one first-order step and 2.7×10⁻³ for one second-order step. A single second-order step (111 048 CNOTs, depth 141 465) therefore meets the 5×10⁻³ criterion used for Models 1 and 2. Model 2 needed 16 steps (3.4×10⁷ CNOTs) for the same accuracy.

**DB-DOI on the full 8-qubit Hamiltonian** (ideal oracle, D ∈ {±ΣZᵢ, ±Zᵢ}, 20 durations in [10⁻³, 1] as in `Example_7`), from the lowest-diagonal basis state |120⟩:

| k | E_k − E₀ (cm⁻¹) | 1 − F (lowest state) | weight outside Z1 | weight outside ⁴I₁₅/₂ | oracle calls | CNOTs |
|---|---|---|---|---|---|---|
| 0 | 143.7 | 0.9989 | 0.930 | 1.3×10⁻⁴ | 0 | 0 |
| 1 | 140.4 | 0.9988 | 0.925 | 1.1×10⁻⁴ | 1 | 1.1×10⁵ |
| 3 | 134.2 | 0.9987 | 0.916 | 9.0×10⁻⁵ | 13 | 1.4×10⁶ |
| 5 | 128.3 | 0.9986 | 0.906 | 8.4×10⁻⁵ | 121 | 1.3×10⁷ |

The reference already lies in ⁴I₁₅/₂, and DB-DOI keeps it there, but it barely moves towards Z1: 15 cm⁻¹ in 5 steps. The compiled replay of the first two steps agrees with the ideal oracle to 0.34 cm⁻¹ and in the infidelity to four digits. The longer duration grid of `Example_12` was not used for this run, which already peaks near 19 GB of memory (deep-copied circuits of up to 121 oracle calls); for Model 2 it improved the energy but not the ground-state weight.

**Two-stage DB-DOI.** Circuit-level preparation inside the 4-qubit block isolated by matrix-level DBQA: the Z1 block at **B*** (1 − F < 10⁻⁹) for Model 3, and the ⁴I₁₅/₂ block (1 − F < 10⁻⁶) for Model 2. It assumes the block Hamiltonian, which is known classically after stage A, as the oracle. D are the 16 fields Σ bᵢZᵢ, and the durations are 40 log-spaced values in [10⁻³, 10²]. The Model 2 block keeps the four symmetry sectors of the full Hamiltonian (the stage-A rotation is built from evolutions under H and under diagonal operators). Its reference must therefore be taken in the sector of the ground doublet: the lowest-diagonal state overall converges to the 0.39 cm⁻¹ doublet instead, the trap of `RESULTS_Model_2.md` §5.3.

| Block | target, gap above it | 1 − F at k = 0, 2, 4, 6, 8 | oracle (error ≤ 5×10⁻³ at t = 10) | CNOTs at k = 3 (explicit circuit) |
|---|---|---|---|---|
| Model 3, Z1 | lowest hyperfine state, 851 MHz | 0.98, 0.90, 0.70, 0.54, 0.26 | 128 second-order steps, 2.6×10⁵ CNOTs per call | 3.4×10⁶ |
| Model 2, ⁴I₁₅/₂ | ground doublet, 0.39 cm⁻¹ | 0.95, 0.92, 0.82, 0.60, 0.31 | 64 steps, 1.2×10⁴ CNOTs per call | 1.6×10⁵ |

Both converge, slowly: after 8 steps (3 280 oracle calls) about 30 % of the state is still outside the target, and most steps choose the longest duration on the grid. Long evolutions also make the oracle expensive (128 Trotter steps at t = 10 for Model 3). The explicit compiled circuits reproduce the ideal-oracle energies and infidelities to three or four digits.

### 4.2 Noise (Example 15)

Every CNOT is followed by a two-qubit depolarizing channel of strength p (single-qubit gates noiseless). The 4-qubit DB-DOI circuits of Models 1, 2 and 3 are simulated exactly, through 256×256 superoperators of the noisy oracle calls composed in the order of the DB-DOI recursion. This reproduces qibo's noisy density-matrix simulation of the explicit circuits to all printed digits and handles circuits of up to 8.5×10⁸ CNOTs; building them as noisy qibo circuits exhausted 24 GB of memory.

| Circuit (4 qubits) | best fidelity with the target over the DB-DOI steps: noiseless | p = 10⁻⁵ | p = 10⁻⁴ | p = 10⁻³ |
|---|---|---|---|---|
| Model 1, Example 3 (2 152 CNOTs per call) | 0.993 (k = 7) | 0.958 (k = 1) | 0.84 (k = 1) | 0.31 (k = 1) |
| Model 1, durations ≤ 10² (17 216 per call) | 0.997 (k = 4) | 0.87 (k = 1) | 0.37 (k = 1) | 0.13 |
| Model 2, two-stage ⁴I₁₅/₂ block (12 160 per call) | 0.69 (k = 8) | 0.15 | 0.13 | 0.125 |
| Model 3, two-stage Z1 block (259 840 per call) | 0.74 (k = 8) | 0.066 | 0.0625 | 0.0625 |

The fully mixed state has fidelity g/16 = 0.125 for a doublet and 0.0625 for a single state. The global-depolarizing estimate F ≈ F₀(1 − p)^N + (1 − (1 − p)^N) g/16 is conservative: for the first Model 1 step (2 152 CNOTs) it gives 0.95, 0.81 and 0.22 at p = 10⁻⁵, 10⁻⁴ and 10⁻³, against 0.958, 0.84 and 0.31 simulated. It then gives the error rate that keeps half of the noiseless signal, p ≈ ln 2 / N:

| Circuit | CNOTs | p for (1 − p)^N = 1/2 |
|---|---|---|
| Model 1, Example 3, k = 6 (1 − F < 10⁻²) | 7.8×10⁵ | 9×10⁻⁷ |
| Model 1, durations ≤ 10², k = 2 (1 − F < 10⁻²) | 6.9×10⁴ | 1.0×10⁻⁵ |
| Model 3, full 8-qubit DB-DOI, k = 5 | 1.3×10⁷ | 5×10⁻⁸ |
| Model 2, two-stage block, k = 8 | 4.0×10⁷ | 2×10⁻⁸ |
| Model 3, two-stage block, k = 8 | 8.5×10⁸ | 8×10⁻¹⁰ |
| Model 2, full 9-qubit DB-DOI, k = 5 | 4.2×10⁹ | 2×10⁻¹⁰ |

With CNOT-ladder Trotter oracles, circuit-level DBQA for these Hamiltonians needs error rates of 10⁻⁵ (Model 1) to 10⁻¹⁰ (Model 2), i.e. fault-tolerant hardware; at two-qubit error rates around 10⁻³ only the first DB-DOI step of Model 1 retains a useful signal.

## 5. Answers

**Can DBQA resolve the hyperfine structure inside a full crystal-field Hamiltonian?** Yes, with two stages and long enough durations. Stage A isolates the 16 hyperfine states of Z1 from the other 240 to 1 − F < 10⁻⁹ in about 105 group-commutator iterations, at every field of the sweep. Stage B then gives every one of the 120 Z1 transition frequencies to 2×10⁻⁴ MHz (2.6×10⁻⁵ MHz from a 10⁻¹² isolation) within 60 iterations. That is 7×10⁻⁹ cm⁻¹, or 10⁻¹² of the spectral width.

**Accurately enough to locate a ZEFOZ point?** Yes: the turning point of ν₆,₇ along the sweep comes out within 0.2 mT (10⁻⁶ isolation) or 0.01 mT (10⁻⁹), for a transition whose frequency changes by 0.08 MHz across 60 mT.

**What carries over between fields?** Everything. The rotation of the previous field isolates the next field's Z1 to 10⁻⁶ immediately, and its schedule of durations, replayed on the new Hamiltonian, reaches 10⁻¹² with no extra step. For a field sweep on hardware, the durations need to be found once.

**What sets the cost?**
- *Durations must reach about 1/gap* in the rescaled units: with s ≤ 10 stage A stalls. This also corrects the Model 2 conclusion.
- *The isolation threshold sets the accuracy*: the restriction error is about 300 MHz, 0.3 MHz, 2×10⁻⁴ MHz and 3×10⁻⁵ MHz at 1 − F = 10⁻³, 10⁻⁶, 10⁻⁹ and 10⁻¹².
- *The generator matters inside the block*: the group commutator works, the canonical bracket stops at a fixed point because the block diagonal is degenerate, and Pauli-Z is slow.

**What does the physics comparison say?** The crystal-field model of Jobbitt et al. reproduces their 7 T hyperfine predictions exactly in PyCF, but its ZEFOZ points differ from those of the measured spin Hamiltonian (90–150 mT near the D1–D2 plane, against 280–670 mT at θ ≈ −40°). DBQA on Model 3 is therefore a benchmark of the algorithm on a realistic hyperfine Hamiltonian, not a prediction of where the ZEFOZ points of the real crystal are.

## 6. Limitations

- **Parameters.** Jobbitt et al.'s site-2 hyperfine parameters do not reproduce the measured Z1 quadrupole tensor (§2.4), so the ZEFOZ points of Model 3 are not those of the crystal (§2.5). The `f11cf` multiplet offsets (§2.2) are unexplained; Model 3 does not use that basis.
- **Effective basis.** As in the paper, the hyperfine model is restricted to ⁴I₁₅/₂ and ⁴I₁₃/₂; J-mixing with higher multiplets enters only through the intermediate-coupled basis.
- **Sweeps.** DBQA was run along one line through one ZEFOZ point (21 fields) and across two avoided crossings (11 fields each). The transferability of the schedule was tested over 3 mT steps; larger steps were not tried.
- **Circuit level.** DB-DOI (§4.1) uses Z-field generators; the full 8-qubit run used durations ≤ 1, and the two-stage run assumes the classically known block Hamiltonian as its oracle. Stage A of the matrix-level DBQA uses D = Δ(Ĥ_k), which is not directly available on hardware. Noise is a depolarizing channel after every CNOT, simulated for the 4-qubit circuits only (§4.2).
- **Scheduling.** Greedy grid search only (100 durations in stage A, 120 in stage B).

## Reproduction

```bash
export PYCF_ROOT=/path/to/pycf
cd examples/Crystal_field_PyCF
jupyter nbconvert --to notebook --execute --inplace Example_8_validation_Er_YSO.ipynb Example_9_DBQA_Er_YSO_generate_data.ipynb Example_10_DBQA_Er_YSO_results.ipynb Example_13_avoided_crossings_Er_YSO.ipynb Example_14_compiled_DBQA_Er_YSO.ipynb Example_15_noisy_DB_DOI.ipynb
```

Run times on a 32-core laptop under WSL2: Example 8 about an hour (the ZEFOZ search on 24 processes), Example 9 84 minutes with another notebook running alongside (cold starts on 8 processes, then the sequential warm and replay chains), Example 10 a few seconds, Example 13 about 4 minutes, Example 14 about 70 minutes (peak memory near 19 GB, from the full 8-qubit DB-DOI), Example 15 about 11 minutes. Example 15 needs the pickle of Example 14. Example 9 reads the sweep chosen by Example 8 from `simulation_results/Er_YSO_validation.pickle`.
