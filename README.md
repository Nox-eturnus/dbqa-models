# Crystal-field Hamiltonians as DBQA benchmarks

Rare-earth ions in crystals are among the best-characterised quantum systems we have. Their energy levels are measured to fractions of a wavenumber, their magnetic and hyperfine structure is known in detail, and they are serious candidates for quantum memories. They also make demanding test problems: a single Hamiltonian mixes interactions whose strengths differ by six orders of magnitude, from the Coulomb repulsion between the 4f electrons down to the coupling of the electron to its own nucleus.

This repository builds three such Hamiltonians, from simple to hard, and runs them through the double-bracket quantum algorithm pipeline of [qiboteam/qibo-dbqa](https://github.com/qiboteam/qibo-dbqa), on which it is based. Each Hamiltonian is built with [PyCF](https://github.com/mikereidnz/pycf) from parameters published for a real crystal, checked against that paper before anything else is done, and then diagonalised with the pipeline at matrix level and as compiled circuits.

The detailed write-up is [`examples/Crystal_field_PyCF/OVERVIEW.md`](examples/Crystal_field_PyCF/OVERVIEW.md), and every number quoted below is traced to an executed notebook in the three `RESULTS_Model_*.md` files next to it.

## The three models

| | Model 1 | Model 2 | Model 3 |
|---|---|---|---|
| System | Ce³⁺ in YAG | Er³⁺ in LiYF₄ | ¹⁶⁷Er³⁺ in Y₂SiO₅ (site 2), in a magnetic field |
| What is included | the whole 4f¹ configuration | the whole 4f¹¹ configuration | the two lowest multiplets, with the nuclear spin |
| States → qubits | 14 → 4 | 364 → 9 | 240 → 8 |
| Source of the parameters | Hu *et al.*, J. Phys.: Condens. Matter **23**, 045501 (2011) | Couto dos Santos *et al.*, J. Alloys Compd. **275–277**, 435 (1998) | Jobbitt *et al.*, Phys. Rev. B **104**, 155121 (2021) |

### Model 1: Ce³⁺:YAG

Cerium has a single 4f electron, so its Hamiltonian has only fourteen states and fits on four qubits. It is the warm-up: small enough to check every step by hand, yet already a genuine crystal-field problem with seven Kramers doublets spread over 4 000 cm⁻¹.

I first checked the conversion from PyCF to Qibo on PyCF's own Ce:LiYF₄ example, where the energies and eigenvectors agree to 10⁻¹¹ cm⁻¹. The group-commutator iteration then brings the ground multiplet within 1 cm⁻¹ of the exact levels in 25 iterations. Once the durations it may choose are long enough (more on this below), all seven levels are within 1 cm⁻¹ after 15 iterations and below 10⁻³ cm⁻¹ well before 60, including the two levels only 94 cm⁻¹ apart that first looked like the bottleneck. As a circuit, the evolution under this Hamiltonian costs 2 152 CNOTs per call, and DB-DOI prepares the ground doublet to 99 % fidelity with 6.9×10⁴ CNOTs.

### Model 2: Er³⁺:LiYF₄

Erbium has eleven 4f electrons, and keeping all of them gives 364 states on nine qubits. Here the multiscale structure appears in full: the multiplets are separated by thousands of wavenumbers, while the crystal field splits each of them by a few hundred.

The published parameters did not reproduce the paper's own energy levels in PyCF, because one of the free-ion parameters follows a convention I could not recover. I kept the paper's crystal field and refitted the free-ion part to the 108 measured levels, which matches experiment slightly better than the original fit (σ = 12.9 against 14.6 cm⁻¹).

The most useful result is a two-stage strategy: let the algorithm first isolate the ground multiplet, then diagonalise that 16-state block on its own. This brings every level of the ground multiplet to within 0.075 cm⁻¹ of the exact value after 16 + 12 iterations. The price of the full Hamiltonian shows at circuit level: about 10⁶ CNOTs per Trotter step, against a few dozen for the spin chains the pipeline was first tested on.

### Model 3: ¹⁶⁷Er³⁺:Y₂SiO₅ with hyperfine structure

The last model is the one closest to current experiments. ¹⁶⁷Er in Y₂SiO₅ is studied as a telecom-band quantum memory, and the quantities of interest there are its hyperfine levels in a magnetic field, in particular the "ZEFOZ" fields at which a hyperfine transition becomes insensitive to magnetic noise. The model keeps the two lowest multiplets of erbium together with its nuclear spin I = 7/2, so it spans everything from 6 500 cm⁻¹ between the multiplets down to avoided crossings a few MHz wide.

The checks came first. The Hamiltonian reproduces the paper's predicted hyperfine splittings at 7 T to within 0.7 MHz. I then compared it with the measured spin Hamiltonian used by Matsuura *et al.* [Phys. Rev. B **113**, 085421 (2026)] to locate ZEFOZ points. From their spin Hamiltonian I recover their five best points in the same order. The crystal-field model, on the other hand, places its avoided crossings at the right fields but gets their widths and its ZEFOZ points wrong. Its quadrupole coupling is several times too small, so it is a realistic benchmark rather than a prediction for the real crystal.

With durations matched to the gaps, the two-stage strategy resolves the hyperfine structure inside the full 240-state Hamiltonian. Every transition frequency of the lowest level comes out to 2×10⁻⁴ MHz, about 10⁻¹² of the spectral width. That is precise enough to locate a ZEFOZ point to better than 0.01 mT, and to resolve an avoided crossing only 1.1 MHz wide. Along a sweep of the magnetic field, the work done at one field carries over almost completely to the next. At circuit level the picture is more sober: preparing these states with DB-DOI is slow, and under realistic gate noise the deeper circuits lose their signal entirely.

## What the three models taught

- **The step durations have to match the gaps.** A double-bracket step can only resolve a gap g if it is allowed to run for a time of about 1/g, measured in the units the Hamiltonian is scaled to. I first ran Models 1 and 2 with durations capped at 1 or 10, and they seemed unable to resolve close or small splittings. The cap was the problem, not the algorithm. `Example_11` and `Example_12` repeat those runs with longer durations, and the results files record which conclusions changed.
- **Isolate first, then resolve.** Splitting the problem into a block-isolation stage and a block-diagonalisation stage is much faster than attacking everything at once, and how well the block is isolated sets the final accuracy.
- **The group commutator is the most reliable generator.** Once the durations are long enough, it does best in every model and at every stage. The canonical bracket can get stuck when the diagonal of a block is degenerate.
- **Circuits are the hard part.** With Trotterised oracles and DB-DOI, keeping half of the signal needs two-qubit error rates of 10⁻⁵ for Model 1 and 10⁻⁸ to 10⁻¹⁰ for the larger models: a task for fault-tolerant hardware.

Where the work departs from the published parameters or from the original plan, the overview says so openly in a table of its own.

## Repository layout

The layout follows qibo-dbqa: library code in `src/`, one folder of example notebooks per topic in `examples/`, and tests in `tests/`.

```
src/qibo_dbqa/                 the qibo-dbqa package
src/qibo_dbqa/utils_crystal_field.py   everything added here: PyCF interface, metrics, schedulers,
                                       compilation, field-dependent models, ZEFOZ and crossing searches
src/db_doi/                    DB-DOI circuits (from qibo-dbqa)
examples/Crystal_field_PyCF/   notebooks Example_1 ... Example_15, OVERVIEW.md, RESULTS_Model_1-3.md,
                               transcribed tables (data/), figures (pdfs_figures/), results (simulation_results/)
tests/test_models_crystal_field.py     tests of the crystal-field code
```

| Notebooks | Model | What they do |
|---|---|---|
| `Example_1` | – | checks the PyCF → Qibo conversion |
| `Example_2`, `Example_3` | 1 | matrix-level DBQA; compiled circuits and DB-DOI |
| `Example_4` to `Example_7` | 2 | validation and refit; DBQA data and results; compiled circuits |
| `Example_8` to `Example_10` | 3 | validation and ZEFOZ search; DBQA over a field sweep; results |
| `Example_11`, `Example_12` | 1, 2 | the runs repeated with longer durations |
| `Example_13` | 3 | avoided crossings, and DBQA across them |
| `Example_14`, `Example_15` | 3 (and 1, 2) | compiled circuits and DB-DOI; the same circuits under noise |

The notebooks of qibo-dbqa itself (the Fig. 4 and Quantinuum examples, which some notebooks here refer to) are not included; they are in the original repository.

## Installation

PyCF builds a C extension that does not compile on native Windows, so everything here was run under WSL2 (Ubuntu) with a conda-forge environment:

```bash
micromamba create -n dbqa-cf -c conda-forge python=3.12 numpy=1.26.4 "scipy>=1.15" gsl nlopt "libblas=*=*openblas" liblapacke openblas cython c-compiler make pkg-config matplotlib jupyterlab
micromamba activate dbqa-cf
CFL_CFLAGS="-I$CONDA_PREFIX/include" CFL_LDLIBS="-L$CONDA_PREFIX/lib -Wl,-rpath,$CONDA_PREFIX/lib -llapacke -lopenblas" LAPACKE_INCLUDE_DIR=$CONDA_PREFIX/include pip install --no-build-isolation /path/to/pycf
pip install qibo==0.2.14 optuna==4.1.0 hyperopt qibojit==0.1.7
pip install --no-deps -e .
```

PyCF does not install its matrix-element files, so point `PYCF_ROOT` at a PyCF source checkout; the notebooks look for one next to this repository by default. Set `OPENBLAS_NUM_THREADS=1`: the matrices here are small, and OpenBLAS's default of one thread per core slowed them down by up to 600 times.

Run the tests with `pytest -o addopts="" tests/` (159 tests, about two minutes; the crystal-field tests are skipped without PyCF) and a notebook with `jupyter nbconvert --to notebook --execute --inplace <notebook>` from `examples/Crystal_field_PyCF`. Every run is deterministic. Run times range from seconds to about an hour and are listed in the overview; `Example_14` needs about 19 GB of memory.

## Credits and license

The double-bracket pipeline (`src/qibo_dbqa` apart from `utils_crystal_field.py`, `src/db_doi`, the package configuration and the original tests) is taken unchanged from [qiboteam/qibo-dbqa](https://github.com/qiboteam/qibo-dbqa) at commit `5aab199` and is distributed under its Apache License 2.0 (see `LICENSE`). The crystal-field code, notebooks, tests and results were added on top of it. The Hamiltonians are built with M. F. Reid's [PyCF](https://github.com/mikereidnz/pycf) and its bundled matrix-element files.

## References

1. L. Hu, M. F. Reid, C. K. Duan, S. Xia and M. Yin, J. Phys.: Condens. Matter **23**, 045501 (2011).
2. M. A. Couto dos Santos, E. Antic-Fidancev, J. Y. Gesland, J. C. Krupa, M. Lemaître-Blaise and P. Porcher, J. Alloys Compd. **275–277**, 435 (1998), doi:10.1016/S0925-8388(98)00363-6.
3. N. L. Jobbitt, J.-P. R. Wells, M. F. Reid, S. P. Horvath, P. Goldner and A. Ferrier, Phys. Rev. B **104**, 155121 (2021), arXiv:2206.09080.
4. K. Matsuura, S. Yasui, R. Kaji, H. Sasakura, T. Tawara and S. Adachi, Phys. Rev. B **113**, 085421 (2026), arXiv:2412.10126.
5. qibo-dbqa, [github.com/qiboteam/qibo-dbqa](https://github.com/qiboteam/qibo-dbqa).
6. PyCF, [github.com/mikereidnz/pycf](https://github.com/mikereidnz/pycf).
