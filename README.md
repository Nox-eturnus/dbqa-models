# Crystal-field Hamiltonians as DBQA benchmarks

Rare-earth ions in crystals are among the most thoroughly characterised quantum systems. Their energy levels are known to fractions of a wavenumber, their magnetic and hyperfine structure has been measured in detail, and several of them are under active study as quantum memories. These properties also make them demanding benchmarks: a single crystal-field Hamiltonian combines interactions whose strengths span six orders of magnitude, from the Coulomb repulsion between the 4f electrons to the hyperfine coupling between the electrons and the nucleus.

This repository constructs three such Hamiltonians of increasing complexity and studies them with the double-bracket quantum algorithm pipeline of [qiboteam/qibo-dbqa](https://github.com/qiboteam/qibo-dbqa), on which it is based. Each Hamiltonian is built with [PyCF](https://github.com/mikereidnz/pycf) from parameters published for a specific crystal, validated against the energy levels and spectroscopic data of the source paper, and then diagonalised with the pipeline, both at matrix level and as compiled circuits.

A detailed account is given in [`examples/Crystal_field_PyCF/OVERVIEW.md`](examples/Crystal_field_PyCF/OVERVIEW.md). Every number quoted below is traced to an executed notebook in the three `RESULTS_Model_*.md` files in the same folder.

## The three models

| | Model 1 | Model 2 | Model 3 |
|---|---|---|---|
| System | Ce³⁺ in YAG | Er³⁺ in LiYF₄ | ¹⁶⁷Er³⁺ in Y₂SiO₅ (site 2), in a magnetic field |
| Configuration space | the full 4f¹ configuration | the full 4f¹¹ configuration | the two lowest multiplets, including the nuclear spin |
| States → qubits | 14 → 4 | 364 → 9 | 240 → 8 |
| Source of the parameters | Hu *et al.*, J. Phys.: Condens. Matter **23**, 045501 (2011) | Couto dos Santos *et al.*, J. Alloys Compd. **275–277**, 435 (1998) | Jobbitt *et al.*, Phys. Rev. B **104**, 155121 (2021) |

### Model 1: Ce³⁺:YAG

Ce³⁺ carries a single 4f electron, so its Hamiltonian has fourteen states and is embedded in four qubits. The model serves as a controlled starting point: it is small enough for every step to be verified directly, yet it is a complete crystal-field problem, with seven Kramers doublets spread over 4 000 cm⁻¹.

The conversion from PyCF to Qibo was first validated on the Ce:LiYF₄ example distributed with PyCF, where energies and eigenvectors agree to 10⁻¹¹ cm⁻¹. The group-commutator iteration then reaches the ground multiplet to within 1 cm⁻¹ after 25 iterations. With durations long enough to resolve the smallest gap (see below), all seven levels lie within 1 cm⁻¹ of the exact spectrum after 15 iterations and within 10⁻³ cm⁻¹ well before 60. This includes the pair of levels separated by only 94 cm⁻¹, which initially appeared to limit convergence. At circuit level, the evolution oracle costs 2 152 CNOTs per call, and DB-DOI prepares the ground doublet with 99 % fidelity using 6.9×10⁴ CNOTs.

### Model 2: Er³⁺:LiYF₄

Er³⁺ has eleven 4f electrons, and the complete configuration comprises 364 states, embedded in nine qubits. The multiscale structure of rare-earth spectra is fully present here: the multiplets are separated by thousands of wavenumbers, whereas the crystal field splits each of them by a few hundred.

The published parameters do not reproduce the paper's calculated levels in PyCF, because one of the free-ion parameters follows a convention that could not be identified. The Model 2 Hamiltonian therefore retains the published crystal field and refits the free-ion parameters to the 108 measured levels. This yields a slightly better agreement with experiment than the original fit (σ = 12.9 against 14.6 cm⁻¹).

The central result is a two-stage strategy, in which the ground multiplet is first isolated and the resulting 16-state block is then diagonalised separately. This places every level of the ground multiplet within 0.075 cm⁻¹ of its exact value after 16 + 12 iterations. The cost of the full configuration becomes apparent at circuit level: a single Trotter step requires about 10⁶ CNOTs, compared with a few dozen for the spin chains on which the pipeline was originally tested.

### Model 3: ¹⁶⁷Er³⁺:Y₂SiO₅ with hyperfine structure

The third model is the closest to current experiments. ¹⁶⁷Er in Y₂SiO₅ is studied as a telecom-band quantum memory, and the relevant quantities are its hyperfine levels in a magnetic field, in particular the zero-first-order-Zeeman (ZEFOZ) fields at which a hyperfine transition becomes insensitive to magnetic noise. The model retains the two lowest multiplets of Er³⁺ together with the nuclear spin I = 7/2, and therefore spans energy scales from 6 500 cm⁻¹ between the multiplets down to avoided crossings a few MHz wide.

The Hamiltonian reproduces the hyperfine splittings at 7 T predicted in the source paper to within 0.7 MHz. It was further compared with the measured spin Hamiltonian employed by Matsuura *et al.* [Phys. Rev. B **113**, 085421 (2026)] to locate ZEFOZ points; the same search applied to that spin Hamiltonian recovers their five optimal points in the reported order. The crystal-field model places its avoided crossings at the correct fields, but neither their widths nor its ZEFOZ points agree with the measured spin Hamiltonian. Its quadrupole coupling is several times too small. Model 3 is therefore a realistic benchmark, not a prediction for the physical crystal.

With durations matched to the relevant gaps, the two-stage strategy resolves the hyperfine structure within the full 240-state Hamiltonian. All transition frequencies of the lowest crystal-field level are obtained to 2×10⁻⁴ MHz, about 10⁻¹² of the spectral width. This precision suffices to locate a ZEFOZ point to better than 0.01 mT and to resolve an avoided crossing of 1.1 MHz. Along a magnetic-field sweep, the rotation obtained at one field carries over almost entirely to the next. At circuit level the outlook is less favourable: DB-DOI state preparation converges slowly, and under realistic gate noise the deeper circuits lose their signal entirely.

## Principal findings

- **Step durations must match the gaps.** A double-bracket step resolves a gap g only if its duration reaches about 1/g, in the units to which the Hamiltonian is rescaled. Models 1 and 2 were first run with durations capped at 1 or 10, which made close and small splittings appear unresolvable. The limitation lay in the duration grid rather than in the algorithm. `Example_11` and `Example_12` repeat these runs with longer durations, and the results files record which conclusions changed.
- **Isolation followed by resolution.** Separating block isolation from block diagonalisation is considerably faster than treating the full Hamiltonian at once, and the quality of the isolation sets the final accuracy.
- **The group commutator is the most reliable generator.** With sufficiently long durations, it performs best in every model and at every stage. The canonical bracket can stall when the diagonal of a block is degenerate.
- **Circuit-level costs are prohibitive for current hardware.** With Trotterised oracles and DB-DOI, retaining half of the noiseless signal requires two-qubit error rates of 10⁻⁵ for Model 1 and 10⁻⁸ to 10⁻¹⁰ for Models 2 and 3, which places these computations in the fault-tolerant regime.

Departures from the published parameters and from the original plan are listed in a dedicated table in the overview.

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
