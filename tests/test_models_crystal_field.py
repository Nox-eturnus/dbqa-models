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
