"""Testing the PyCF -> Qibo crystal-field interface used for DBQA"""

import os
from functools import reduce

import numpy as np
import pytest
from qibo import Circuit
from scipy.linalg import expm

from qibo_dbqa.double_bracket import (
    DoubleBracketGeneratorType,
    DoubleBracketIteration,
)
from qibo_dbqa.utils_crystal_field import *

NSTEPS = 3
TOL = 1e-8


@pytest.fixture(scope="module")
def f1cf():
    """The PyCF 4f1 basis; PyCF tests are skipped without PyCF or PYCF_ROOT."""
    import_sljm = pytest.importorskip("pycf.import_sljm")
    if os.environ.get("PYCF_ROOT") is None:  # pragma: no cover
        pytest.skip("PYCF_ROOT is not set")
    return import_sljm.ImportSLJM(pycf_matel_path("f1cf", "ceylf"))


ER_LIYF4_COUTO1998_RACAH = [
    ER_LIYF4_COUTO1998["E1"],
    ER_LIYF4_COUTO1998["E2"],
    ER_LIYF4_COUTO1998["E3"],
]


def ce_yag_coeff(parameter_set="fitted"):
    return dict(CE_YAG_HU2011[parameter_set], EAVG=0.0)


def test_dense_matrix_matches_pycf(f1cf):
    pycf_cfl = pytest.importorskip("pycf.cfl")
    tensors = [getattr(f1cf, name) for name in CE_YLF_PYCF_EXAMPLE]
    h = pycf_cfl.Hamiltonian(tensors)
    h.set_coeff(CE_YLF_PYCF_EXAMPLE)
    w, z = h.diag()
    matrix = pycf_dense_matrix(tensors, CE_YLF_PYCF_EXAMPLE)
    assert np.allclose(matrix, matrix.conj().T)
    assert np.allclose(np.linalg.eigvalsh(matrix), np.sort(w), atol=TOL)
    assert np.allclose(matrix @ z, z * w, atol=TOL)


@pytest.mark.parametrize("k", [2, 4, 6])
def test_wybourne_normalisation(f1cf, k):
    """Each Hermitian combination C^k_q + C^k_{-q} (q>0) must carry twice the HS norm of C^k_0."""
    norm_0 = np.linalg.norm(pycf_dense_matrix([getattr(f1cf, f"C{k}0")], {f"C{k}0": 1.0}))
    for q in range(1, k + 1):
        name = f"C{k}{q}"
        if hasattr(f1cf, name):
            norm_q = np.linalg.norm(pycf_dense_matrix([getattr(f1cf, name)], {name: 1.0}))
            assert np.isclose(norm_q**2, 2 * norm_0**2)


def test_rotation_invariance_ce_yag(f1cf):
    """Rotating about z by angle a multiplies B^k_q by exp(i q a) and leaves the spectrum unchanged."""
    coeff = ce_yag_coeff()
    tensors = [f1cf.EAVG] + [getattr(f1cf, name) for name in CE_YAG_TENSORS]
    energies = np.linalg.eigvalsh(pycf_dense_matrix(tensors, coeff))
    angle = 0.3
    rotated = {
        name: value * np.exp(1j * int(name[2]) * angle) if name.startswith("C") else value
        for name, value in coeff.items()
    }
    assert np.allclose(
        np.linalg.eigvalsh(pycf_dense_matrix(tensors, rotated)), energies, atol=TOL
    )


@pytest.mark.parametrize("parameter_set", ["fitted", "calculated"])
def test_embedding_ce_yag(backend, f1cf, parameter_set):
    tensors = [f1cf.EAVG] + [getattr(f1cf, name) for name in CE_YAG_TENSORS]
    coeff = ce_yag_coeff(parameter_set)
    physical = np.linalg.eigvalsh(pycf_dense_matrix(tensors, coeff)) / 1000.0
    h, n_physical = crystal_field_hamiltonian(tensors, coeff, backend=backend)
    matrix = backend.to_numpy(h.matrix)
    assert h.nqubits == 4 and n_physical == 14
    assert np.allclose(matrix[:n_physical, n_physical:], 0)
    energies = np.linalg.eigvalsh(matrix)
    assert np.allclose(energies[:n_physical], physical, atol=TOL)
    assert np.all(energies[n_physical:] > physical[-1])
    groups = group_levels(energies, tol=1e-5)
    assert [len(g) for g in groups] == [2] * 7 + [2]


def test_embedding_padding_offset():
    matrix = np.diag([0.0, 1.0, 2.0])
    padded, nqubits, n_physical = embed_in_qubits(matrix, padding_offset=0.5)
    assert nqubits == 2 and n_physical == 3
    assert np.isclose(padded[3, 3], 2.5)
    padded, _, _ = embed_in_qubits(matrix)
    assert np.isclose(padded[3, 3], 2.0 + 0.25 * 2.0)


def test_group_levels_and_tolerance():
    energies = np.array([0.0, 1e-9, 1.0, 1.0 + 1e-9, 1.5, 3.0])
    groups = group_levels(energies, tol=1e-6)
    assert groups == [[0, 1], [2, 3], [4], [5]]
    max_inside, min_between = check_grouping_tolerance(energies, tol=1e-6)
    assert np.isclose(max_inside, 1e-9) and np.isclose(min_between, 0.5)


def test_subspace_fidelity_and_gap_error():
    identity = np.eye(4)
    rotation = np.array([[np.cos(0.3), -np.sin(0.3)], [np.sin(0.3), np.cos(0.3)]])
    mixed = identity[:, :2] @ rotation
    assert np.isclose(subspace_fidelity(identity[:, :2], mixed), 1.0)
    assert np.isclose(subspace_fidelity(identity[:, :1], mixed[:, :1]), np.cos(0.3) ** 2)
    assert np.isclose(gap_error(np.array([0.0, 1.1]), np.array([0.0, 1.0]), 0, 1), 0.1)


def test_pauli_decomposition(backend):
    matrix = random_hermitian_matrix(3)
    decomposition = pauli_decomposition(matrix)
    rebuilt = sum(
        c * reduce(np.kron, [PAULI_MATRICES[p] for p in s]) for s, c in decomposition.items()
    )
    assert np.allclose(rebuilt, matrix)
    h = PauliCompiledHamiltonian(matrix, backend=backend)
    assert np.allclose(backend.to_numpy(h.matrix), matrix)


@pytest.mark.parametrize("string", ["XIZ", "YYI", "ZXY", "IIZ", "XYZ"])
def test_pauli_rotation(backend, string):
    theta = 0.37
    circuit = add_pauli_rotation(Circuit(3), string, theta)
    pauli = reduce(np.kron, [PAULI_MATRICES[p] for p in string])
    assert np.allclose(circuit.unitary(backend=backend), expm(-1j * theta * pauli))
    assert two_qubit_gate_count(circuit) == 2 * (pauli_weight(string) - 1)


@pytest.mark.parametrize("ts_order,steps,tol", [(1, 40, 2e-2), (2, 4, 1e-3)])
def test_trotter_evolution(backend, ts_order, steps, tol):
    matrix = random_hermitian_matrix(3)
    h = PauliCompiledHamiltonian(matrix, trotter_steps=steps, ts_order=ts_order, backend=backend)
    dt = 0.2
    traceless = matrix - h.identity_coefficient * np.eye(len(matrix))
    error = np.linalg.norm(h.circuit(dt).unitary(backend=backend) - expm(-1j * dt * traceless), 2)
    assert error < tol
    exact = ExactEvolutionHamiltonian(matrix, backend=backend)
    assert np.allclose(exact.circuit(dt).unitary(backend=backend), expm(-1j * dt * matrix))


def test_basis_state_circuit(backend):
    state = backend.to_numpy(basis_state_circuit(4, 5)().state())
    assert np.isclose(abs(state[5]), 1.0)


@pytest.mark.parametrize(
    "mode",
    [
        DoubleBracketGeneratorType.canonical,
        DoubleBracketGeneratorType.single_commutator,
        DoubleBracketGeneratorType.group_commutator,
    ],
)
def test_spectral_grid_search_step(backend, mode):
    """The spectral grid search must reproduce the losses of `DoubleBracketIteration.loss`."""
    from qibo.hamiltonians import Hamiltonian

    matrix = random_hermitian_matrix(3)
    dbi = DoubleBracketIteration(Hamiltonian(3, backend.cast(matrix), backend=backend), mode=mode)
    d = np.diag(np.linspace(-1.0, 1.0, 8)) if mode is not DoubleBracketGeneratorType.canonical else None
    space = np.linspace(1e-3, 0.5, 15)
    step, loss = spectral_grid_search_step(dbi, space, d=d)
    losses = [dbi.loss(s, d=d) for s in space]
    assert np.isclose(loss, np.min(losses))
    assert np.isclose(step, space[int(np.argmin(losses))])


@pytest.mark.parametrize("method", ["canonical", "Pauli-Z", "group commutator"])
def test_run_dbqa_history(method):
    matrix = random_hermitian_matrix(3)
    energies = np.linalg.eigvalsh(matrix)
    groups = group_levels(energies, 1e-8)
    d_list = [np.diag(np.linspace(-1.0, 1.0, 8)), np.diag(np.arange(8) % 2 * 2.0 - 1)]
    run = run_dbqa_history(
        matrix, method, dict(levels=groups), 4, np.linspace(1e-3, 0.5, 20),
        d_list=d_list if method == "Pauli-Z" else None, checkpoints=(0, 4),
    )
    norms = [r["levels"]["off_diagonal_norm"] for r in run["history"]]
    assert len(run["history"]) == 5 and len(run["steps"]) == 4
    assert norms[-1] < norms[0]
    assert set(run["checkpoints"]) == {0, 4}
    assert np.allclose(np.linalg.eigvalsh(run["checkpoints"][4]), energies)


def test_racah_slater_round_trip():
    f = racah_to_slater(*ER_LIYF4_COUTO1998_RACAH)
    assert np.allclose(slater_to_racah(*f), ER_LIYF4_COUTO1998_RACAH)
    # Carnall's Er3+:LaF3 Slater parameters give E^1, E^2, E^3 close to 6602, 33.1, 624
    assert np.allclose(slater_to_racah(97483, 67904, 54010), [6602, 33.06, 624.3], atol=1.0)


def test_trotter_cnot_count():
    decomposition = {"III": 1.0, "XIZ": 0.5, "YYZ": 0.2, "IIZ": 0.1}
    assert trotter_cnot_count(decomposition, ts_order=1) == 2 + 4 + 0
    assert trotter_cnot_count(decomposition, ts_order=2) == 2 * 6 - 0


ER_TABLE = os.path.join(
    os.path.dirname(__file__),
    "..",
    "examples",
    "Crystal_field_PyCF",
    "data",
    "Er_LiYF4_CoutoDosSantos1998_Table1.csv",
)


@pytest.fixture(scope="module")
def f11cf():
    """The PyCF 4f11 basis bundled with the Er:YSO example."""
    import_sljm = pytest.importorskip("pycf.import_sljm")
    if os.environ.get("PYCF_ROOT") is None:  # pragma: no cover
        pytest.skip("PYCF_ROOT is not set")
    return import_sljm.ImportSLJM(pycf_matel_path("f11cf", "eryso"))


def er_table():
    import csv

    with open(ER_TABLE) as file:
        return list(csv.DictReader(line for line in file if not line.startswith("#")))


def test_er_table_transcription():
    rows = er_table()
    e_exp = [float(r["E_exp"]) for r in rows if r["E_exp"]]
    e_calc = np.array([float(r["E_calc"]) for r in rows])
    assert len(rows) == 134 and len(e_exp) == 108
    assert np.all(np.diff(e_calc) >= 0)
    # the paper's own fit: sigma = 14.6 cm^-1 for 21 parameters
    residuals = np.array([float(r["E_calc"]) - float(r["E_exp"]) for r in rows if r["E_exp"]])
    assert np.isclose(np.sqrt(np.sum(residuals**2) / (len(residuals) - 21)), 14.6, atol=0.05)


@pytest.mark.parametrize("parameters,sigma_max", [(ER_LIYF4_FREE_ION_REFIT, 13.0), (ER_LIYF4_COUTO1998, None)])
def test_er_liyf4_levels(f11cf, parameters, sigma_max):
    coeff = er_liyf4_pycf_coeff(parameters)
    matrix = pycf_dense_matrix([getattr(f11cf, name) for name in coeff], coeff)
    assert matrix.shape == (364, 364) and np.allclose(matrix, matrix.conj().T)
    energies = np.linalg.eigvalsh(matrix)
    assert np.allclose(energies[::2], energies[1::2], atol=1e-6)  # Kramers doublets
    rows = er_table()
    mask = np.array([r["E_exp"] != "" for r in rows])
    e_exp = np.array([float(r["E_exp"]) for r in rows if r["E_exp"]])
    levels = energies[::2][: len(rows)][mask]
    residuals = levels - e_exp + np.mean(e_exp - levels)
    sigma = np.sqrt(np.sum(residuals**2) / (len(residuals) - 21))
    if sigma_max is None:
        # the paper's E^3 convention differs from PyCF's: its parameters taken literally fail by far
        assert sigma > 100
    else:
        assert sigma < sigma_max


def test_er_liyf4_embedding(backend, f11cf):
    coeff = er_liyf4_pycf_coeff(ER_LIYF4_FREE_ION_REFIT)
    h, n_physical = crystal_field_hamiltonian(
        [getattr(f11cf, name) for name in coeff], coeff, energy_unit=20000.0, backend=backend
    )
    matrix = backend.to_numpy(h.matrix)
    assert h.nqubits == 9 and n_physical == 364
    assert np.allclose(matrix[:n_physical, n_physical:], 0)


def random_hermitian_matrix(nqubits, seed=10):
    rng = np.random.default_rng(seed)
    a = rng.normal(size=(2**nqubits, 2**nqubits)) + 1j * rng.normal(size=(2**nqubits, 2**nqubits))
    return (a + a.conj().T) / 4


def test_dbqa_report_ce_yag(backend, f1cf):
    tensors = [f1cf.EAVG] + [getattr(f1cf, name) for name in CE_YAG_TENSORS]
    h, _ = crystal_field_hamiltonian(tensors, ce_yag_coeff(), backend=backend)
    dbi = DoubleBracketIteration(h, mode=DoubleBracketGeneratorType.canonical)
    energies = np.linalg.eigvalsh(backend.to_numpy(h.matrix))
    groups = group_levels(energies, tol=1e-5)
    initial = dbqa_spectrum_report(dbi, groups, energies)
    for _ in range(NSTEPS):
        dbi(step=dbi.choose_step(step_max=0.2))
    final = dbqa_spectrum_report(dbi, groups, energies)
    assert final["off_diagonal_norm"] < initial["off_diagonal_norm"]
    assert np.all(final["level_fidelity"] <= 1 + TOL)
    assert np.max(final["level_energy_error"]) < np.max(initial["level_energy_error"])
