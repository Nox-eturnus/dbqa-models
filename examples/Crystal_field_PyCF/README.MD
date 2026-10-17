DBQA for rare-earth crystal-field Hamiltonians built with PyCF (https://github.com/mikereidnz/pycf).

The interface lives in `src/qibo_dbqa/utils_crystal_field.py` and is tested in `tests/test_models_crystal_field.py`. It covers:
- conversion of PyCF tensors to a dense Hamiltonian and its embedding into qubits, with uncoupled padding states;
- spectral metrics for the DBQA iterates: level grouping, energy errors, projector fidelities and gap errors;
- compilation of the evolution oracle. `PauliCompiledHamiltonian` is a Trotterized CNOT circuit; `ExactEvolutionHamiltonian` is an ideal oracle.

Notebooks:
- `Example_1_validation_PyCF_to_Qibo_Ce_YLF.ipynb` checks the PyCF-to-Qibo conversion on the Ce:YLF example bundled with PyCF. It shows that the eigenvalues and Kramers-doublet eigenspaces survive the 4-qubit embedding, that the padding is uncoupled, and that the crystal-field tensors are Wybourne-normalised.
- `Example_2_DBQA_Ce_YAG.ipynb` is Model 1, Ce3+:YAG (4f1, 14 states in 16 dimensions), run at matrix level. Parameters come from Table 10 of Hu, Reid, Duan, Xia and Yin, J. Phys.: Condens. Matter 23, 045501 (2011). The notebook compares canonical, Pauli-Z and group-commutator DBQA. It reports the tolerance tables, the closest-pair analysis and the relation between fidelity and energy error, and repeats the benchmark with the second parameter set.
- `Example_3_compiled_DBQA_Ce_YAG.ipynb` is Model 1 compiled into circuits. It runs DB-DOI ground-state preparation with `GroupCommutatorIteration` and magnetic-field diagonal operators. The evolution oracle is a Trotterized Pauli-rotation circuit. The notebook reports CNOT counts, circuit depth and the effect of the Trotter error.

- `Example_4_validation_Er_LiYF4.ipynb` builds and validates Model 2, Er3+:LiYF4 with the full 4f11 configuration: 364 states in 512 dimensions, i.e. 9 qubits. The parametric Hamiltonian is that of Couto dos Santos et al., J. Alloys Compd. 275-277, 435 (1998), in the D2d approximation.
  - The paper's Racah E3 does not follow the standard Racah/Carnall convention, so its parameters fail in PyCF: σ = 547 cm^-1.
  - The notebook diagnoses this and builds the Model 2 Hamiltonian from the paper's crystal field plus a refit of the free-ion parameters to the 108 measured levels. That gives σ = 12.9 cm^-1, against the paper's 14.6.
- `Example_5_DBQA_Er_LiYF4_generate_data.ipynb` runs DBQA on Model 2, which takes about 35 minutes. It first runs canonical, Pauli-Z and group-commutator DBQA on the full Hamiltonian. It then runs a second stage restricted to the isolated ground multiplet.
- `Example_6_DBQA_Er_LiYF4_results.ipynb` analyses the Model 2 data. It covers convergence of the free-ion manifolds and of the crystal-field levels, the two-stage strategy from every stage-1 method at two isolation thresholds, and the compilation cost against Model 1 and the XXZ chain.
- `Example_7_compiled_DBQA_Er_LiYF4.ipynb` is Model 2 compiled into circuits, and takes about an hour.
  - It reports the Trotter error, exact CNOT count and exact depth of the compiled oracle. These come from the exact unitary of the compiled gate sequence, since the circuits are too large to build gate by gate.
  - It runs symmetry-resolved DB-DOI ground-state preparation, with a compiled-oracle replay.
- `Example_8_validation_Er_YSO.ipynb` builds and validates Model 3: 167Er3+ at site 2 of Y2SiO5, the two lowest multiplets times the nuclear spin I = 7/2 (240 states in 256 dimensions, 8 qubits) in a magnetic field. The parameters are those of Jobbitt et al., Phys. Rev. B 104, 155121 (2021).
  - It checks the parameters against the paper's crystal-field strengths, levels, g values and 7 T hyperfine splittings (Tables II, III and VI).
  - It reproduces the site-2 ZEFOZ points of Matsuura et al., Phys. Rev. B 113, 085421 (2026), from their spin Hamiltonian. It then searches the ZEFOZ points of the crystal-field model and chooses the field sweep for DBQA.
- `Example_9_DBQA_Er_YSO_generate_data.ipynb` runs two-stage DBQA over that field sweep. Stage A isolates the 16 hyperfine states of the lowest crystal-field level Z1 on 8 qubits; stage B resolves them on 4 qubits. Cold, warm and replayed starts are compared.
- `Example_10_DBQA_Er_YSO_results.ipynb` analyses the Model 3 data. It covers the isolation, the accuracy of the hyperfine transition frequencies, the location of the ZEFOZ point from DBQA, and the compilation cost.
- `Example_11_step_range_Er_LiYF4.ipynb` repeats the group-commutator stage 1 of Model 2 with durations up to 10^4 instead of 10. Every level of the ground multiplet then comes within 0.19 cm^-1, which corrects the Model 2 conclusion that greedy DBQA cannot resolve crystal-field levels (`RESULTS_Model_2.md` §8).
- `Example_12_duration_grids_Ce_YAG_Er_LiYF4.ipynb` checks the remaining Model 1 and Model 2 runs against the same duration-grid criterion: Model 1 at matrix level, the Model 2 multiplet block, and DB-DOI for both (`RESULTS_Model_1.md` §7, `RESULTS_Model_2.md` §8).
- `Example_13_avoided_crossings_Er_YSO.ipynb` locates the avoided crossings of the Z1 hyperfine levels in the crystal-field model and in the measured spin Hamiltonian, checks that Model 3 has no symmetry sectors, and runs two-stage DBQA across two crossings with pair-subspace fidelities and gap errors.
- `Example_14_compiled_DBQA_Er_YSO.ipynb` compiles Model 3 into circuits (Trotter error, CNOT count, depth) and runs DB-DOI on the full 8-qubit Hamiltonian. It also runs two-stage DB-DOI, i.e. circuit-level state preparation inside the 4-qubit blocks isolated by matrix-level DBQA, for Models 2 and 3.
- `Example_15_noisy_DB_DOI.ipynb` simulates the compiled DB-DOI circuits of Models 1, 2 and 3 with a depolarizing channel after every CNOT and derives error-rate budgets for the full circuits.

`data/Er_LiYF4_CoutoDosSantos1998_Table1.csv` is a transcription of Table 1 of that paper: experimental and calculated energies of the lowest 134 Kramers doublets. `data/Er_YSO_site2_Jobbitt2021_Table2.csv` transcribes Table II of Jobbitt et al.: calculated and experimental energies and g values of site 2.

`OVERVIEW.md` is the entry point for all three models. `RESULTS_Model_1.md`, `RESULTS_Model_2.md` and `RESULTS_Model_3.md` summarise the results. `pdfs_figures` contains the figures and `simulation_results` the pickled data. The iterated Model 2 matrices (`Er_LiYF4_DBQA_checkpoints.npz`) and the cached compiled Trotter steps (`Er_LiYF4_trotter_steps.npz`) are large and are not under version control.

Setup:
PyCF builds a C extension against GSL, NLopt and LAPACKE and does not support native Windows. On Windows, the work was carried out under WSL2 with a conda-forge environment:

```bash
micromamba create -n dbqa-cf -c conda-forge python=3.12 numpy=1.26.4 "scipy>=1.15" gsl nlopt "libblas=*=*openblas" liblapacke openblas cython c-compiler make pkg-config matplotlib jupyterlab
CFL_CFLAGS="-I$CONDA_PREFIX/include" CFL_LDLIBS="-L$CONDA_PREFIX/lib -Wl,-rpath,$CONDA_PREFIX/lib -llapacke -lopenblas" LAPACKE_INCLUDE_DIR=$CONDA_PREFIX/include pip install --no-build-isolation /path/to/pycf
pip install qibo==0.2.14 optuna==4.1.0 hyperopt qibojit==0.1.7
pip install --no-deps -e .   # from the root of this repository
```

The matrix-element bases (for example `examples/ceylf/matel/f1cf`) are not installed with PyCF. The notebooks and tests look for them in the PyCF source tree given by `PYCF_ROOT`. The notebooks default to `../../../pycf`, a PyCF checkout next to `qibo-dbqa`. The tests that need PyCF are skipped when it or `PYCF_ROOT` is missing.

Energies are passed to Qibo in units of 1000 cm^-1 (Model 1), 20 000 cm^-1 (Model 2) or 5000 cm^-1 (Model 3), so the Hamiltonian norm is O(1). The DBR duration grid must still reach about the inverse of the smallest gap to be resolved in those units; see `OVERVIEW.md`.

Set `OPENBLAS_NUM_THREADS=1` (the Model 3 notebooks do so before importing numpy). With one OpenBLAS thread per core, a 240x240 diagonalisation took 3 s instead of 5 ms on a 32-core machine, and forked worker processes inherit busy-waiting BLAS threads.
