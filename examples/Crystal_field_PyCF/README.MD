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
- `Example_5_DBQA_Er_LiYF4_generate_data.ipynb` runs DBQA on Model 2, which takes about 1.5 hours. It first runs canonical, Pauli-Z and group-commutator DBQA on the full Hamiltonian. It then runs a second stage restricted to the isolated ground multiplet.
- `Example_6_DBQA_Er_LiYF4_results.ipynb` analyses the Model 2 data. It covers convergence of the free-ion manifolds and of the crystal-field levels, the two-stage strategy, and the compilation cost against Model 1 and the XXZ chain.

`data/Er_LiYF4_CoutoDosSantos1998_Table1.csv` is a transcription of Table 1 of that paper: experimental and calculated energies of the lowest 134 Kramers doublets.

`RESULTS_Model_1.md` and `RESULTS_Model_2.md` summarise the results. `pdfs_figures` contains the figures and `simulation_results` the pickled data. The iterated Model 2 matrices (`Er_LiYF4_DBQA_checkpoints.npz`) are large and are not under version control.

Setup:
PyCF builds a C extension against GSL, NLopt and LAPACKE and does not support native Windows. On Windows we run under WSL2 with a conda-forge environment:

```bash
micromamba create -n dbqa-cf -c conda-forge python=3.12 numpy=1.26.4 "scipy>=1.15" gsl nlopt "libblas=*=*openblas" liblapacke openblas cython c-compiler make pkg-config matplotlib jupyterlab
CFL_CFLAGS="-I$CONDA_PREFIX/include" CFL_LDLIBS="-L$CONDA_PREFIX/lib -Wl,-rpath,$CONDA_PREFIX/lib -llapacke -lopenblas" LAPACKE_INCLUDE_DIR=$CONDA_PREFIX/include pip install --no-build-isolation /path/to/pycf
pip install qibo==0.2.14 optuna==4.1.0 hyperopt qibojit==0.1.7
pip install --no-deps -e /path/to/qibo-dbqa
```

The matrix-element bases (for example `examples/ceylf/matel/f1cf`) are not installed with PyCF. The notebooks and tests look for them in the PyCF source tree given by `PYCF_ROOT`. The notebooks default to `../../../pycf`, a PyCF checkout next to `qibo-dbqa`. The tests that need PyCF are skipped when it or `PYCF_ROOT` is missing.

Energies are passed to Qibo in units of 1000 cm^-1, so the Hamiltonian norm is O(1) and the default DBR step ranges apply.
