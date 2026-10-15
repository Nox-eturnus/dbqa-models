DBQA for rare-earth crystal-field Hamiltonians built with PyCF (https://github.com/mikereidnz/pycf).

The interface lives in `src/qibo_dbqa/utils_crystal_field.py` and is tested in `tests/test_models_crystal_field.py`.

Notebooks:
- `Example_1_validation_PyCF_to_Qibo_Ce_YLF.ipynb` checks the PyCF-to-Qibo conversion on the Ce:YLF example bundled with PyCF: eigenvalues and Kramers-doublet eigenspaces survive the 4-qubit embedding, and the padding states are uncoupled.
- `Example_2_DBQA_Ce_YAG.ipynb` is Model 1, Ce3+:YAG (4f1, 14 states in 16 dimensions). Parameters are from Table 10 of Hu, Reid, Duan, Xia and Yin, J. Phys.: Condens. Matter 23, 045501 (2011). It compares canonical, Pauli-Z and group-commutator DBQA.

Setup:
PyCF builds a C extension against GSL, NLopt and LAPACKE and does not support native Windows. On Windows we run under WSL2 with a conda-forge environment:

```bash
micromamba create -n dbqa-cf -c conda-forge python=3.12 numpy=1.26.4 "scipy>=1.15" gsl nlopt "libblas=*=*openblas" liblapacke openblas cython c-compiler make pkg-config matplotlib jupyterlab
CFL_CFLAGS="-I$CONDA_PREFIX/include" CFL_LDLIBS="-L$CONDA_PREFIX/lib -Wl,-rpath,$CONDA_PREFIX/lib -llapacke -lopenblas" LAPACKE_INCLUDE_DIR=$CONDA_PREFIX/include pip install --no-build-isolation /path/to/pycf
pip install qibo==0.2.14 optuna==4.1.0 hyperopt qibojit==0.1.7
pip install --no-deps -e /path/to/qibo-dbqa
```

The matrix-element bases (for example `examples/ceylf/matel/f1cf`) are not installed with PyCF. The notebooks look for them in the PyCF source tree given by `PYCF_ROOT`, which defaults to `../../../pycf`, a PyCF checkout next to `qibo-dbqa`.

Energies are passed to Qibo in units of 1000 cm^-1, so the Hamiltonian norm is O(1) and the default DBR step ranges apply.

The folder `simulation_results` contains the pickled DBQA histories.
