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


@pytest.mark.parametrize("string", ["XIZ", "YYI", "ZXY", "IIZ", "XYZ"])
def test_apply_pauli_rotation(string):
    theta = 0.37
    pauli = reduce(np.kron, [PAULI_MATRICES[p] for p in string])
    rng = np.random.default_rng(3)
    vector = rng.normal(size=8) + 1j * rng.normal(size=8)
    assert np.allclose(apply_pauli_rotation(np.eye(8, dtype=complex), string, theta), expm(-1j * theta * pauli))
    assert np.allclose(apply_pauli_rotation(vector, string, theta), expm(-1j * theta * pauli) @ vector)


@pytest.mark.parametrize("ts_order,steps", [(1, 1), (1, 3), (2, 1), (2, 2)])
def test_compiled_trotter_unitary_and_depth(backend, ts_order, steps):
    """The matrix-free unitary and depth must match the explicit compiled circuit."""
    matrix = random_hermitian_matrix(3)
    h = PauliCompiledHamiltonian(matrix, trotter_steps=steps, ts_order=ts_order, backend=backend)
    dt = 0.3
    circuit = h.circuit(dt)
    step = compiled_trotter_unitary(h.decomposition, dt / steps, ts_order)
    assert np.allclose(circuit.unitary(backend=backend), np.linalg.matrix_power(step, steps))
    assert compiled_trotter_depth(h.decomposition, ts_order, steps) == circuit.depth
    oracle = CompiledUnitaryHamiltonian(matrix, trotter_steps=steps, ts_order=ts_order, backend=backend)
    assert np.allclose(oracle.circuit(dt).unitary(backend=backend), circuit.unitary(backend=backend))


def test_dense_symbolic_hamiltonian(backend):
    matrix = random_hermitian_matrix(3)
    rng = np.random.default_rng(5)
    state = rng.normal(size=8) + 1j * rng.normal(size=8)
    state /= np.linalg.norm(state)
    for h in (ExactEvolutionHamiltonian(matrix, backend=backend), PauliCompiledHamiltonian(matrix, backend=backend)):
        assert np.allclose(backend.to_numpy(h.matrix), matrix)
        assert np.isclose(h.expectation(backend.cast(state)), np.real(state.conj() @ matrix @ state))


def test_run_dbqa_history_capture():
    matrix = random_hermitian_matrix(3)
    energies = np.linalg.eigvalsh(matrix)
    groups = group_levels(energies, 1e-8)
    run = run_dbqa_history(
        matrix, "group commutator", dict(levels=groups), 30, np.linspace(1e-3, 1.0, 50),
        capture=("levels", 0, (1e-1, 1e-3)),
    )
    for threshold, (k, captured_matrix) in run["captured"].items():
        infidelities = [1 - r["levels"]["level_fidelity"][0] for r in run["history"]]
        assert infidelities[k] < threshold and all(f >= threshold for f in infidelities[:k])
        assert np.allclose(np.linalg.eigvalsh(captured_matrix), energies)


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


ER_YSO_TABLE = os.path.join(
    os.path.dirname(__file__), "..", "examples", "Crystal_field_PyCF", "data", "Er_YSO_site2_Jobbitt2021_Table2.csv"
)
WANG = ER_YSO_SITE2_WANG_SPIN_HAMILTONIAN


@pytest.fixture(scope="module")
def erhfs():
    """The PyCF hyperfine basis (4I15/2 and 4I13/2 times I = 7/2, 240 states) of the Er:YSO example."""
    import_sljm = pytest.importorskip("pycf.import_sljm")
    if os.environ.get("PYCF_ROOT") is None:  # pragma: no cover
        pytest.skip("PYCF_ROOT is not set")
    return import_sljm.ImportSLJM(pycf_matel_path("erhfs", "eryso"))


@pytest.fixture(scope="module")
def er_yso_model(erhfs):
    return hyperfine_field_model(erhfs, er_yso_pycf_coeff(hyperfine=ER_YSO_SITE2_HYPERFINE_JOBBITT2021))


def er_yso_table():
    import csv

    with open(ER_YSO_TABLE) as file:
        return list(csv.DictReader(line for line in file if not line.startswith("#")))


def wang_spin_hamiltonian():
    return spin_hamiltonian_model(WANG["A"], WANG["Q"], WANG["g"])


def test_er_yso_transcription():
    """Table III of Jobbitt et al. lists S2, S4, S6 = 363.1, 653.3, 151.5 cm^-1 for site 2; Table II has 57
    experimental levels and 25, 24 and 21 g values along D1, D2 and b (their Sec. III C)."""
    coeff = er_yso_pycf_coeff()
    strengths = [crystal_field_strength(coeff, k) for k in (2, 4, 6)]
    assert np.allclose(strengths, [363.1, 653.3, 151.5], atol=0.05)
    rows = er_yso_table()
    assert len(rows) == 59 and sum(r["E_exp"] != "" for r in rows) == 57
    assert [sum(r[c] != "" for r in rows) for c in ("gD1_exp", "gD2_exp", "gb_exp")] == [25, 24, 21]


def test_nuclear_spin_matrices():
    ix, iy, iz = nuclear_spin_matrices(3.5)
    assert np.allclose(ix @ iy - iy @ ix, 1j * iz)
    assert np.allclose(ix @ ix + iy @ iy + iz @ iz, 3.5 * 4.5 * np.eye(8))
    assert np.allclose(np.diag(iz), np.arange(3.5, -4, -1))


def test_polar_field():
    field = polar_field(0.63352, -37.5383, -10.9417)
    assert np.isclose(np.linalg.norm(field), 0.63352)
    # Matsuura et al. print (378.9, 73.2, 502.3) mT up to signs
    assert np.allclose(np.abs(field) * 1e3, [378.9, 73.2, 502.3], atol=0.2)


def test_transition_derivatives_finite_differences():
    h0, operators = wang_spin_hamiltonian()
    field = np.array([0.05, -0.02, 0.03])

    def frequency(b):
        energies = np.linalg.eigvalsh(linear_field_hamiltonian(h0, operators, b))
        return energies[9] - energies[4]

    nu, gradient, hessian = transition_derivatives(h0, operators, field, 4, 9)
    assert np.isclose(nu, frequency(field))
    unit, h = np.eye(3), 1e-6
    numeric_gradient = [(frequency(field + h * e) - frequency(field - h * e)) / (2 * h) for e in unit]
    assert np.allclose(gradient, numeric_gradient, rtol=1e-5)
    h = 1e-4
    numeric_hessian = [
        [
            (
                frequency(field + h * (a + b))
                - frequency(field + h * (a - b))
                - frequency(field - h * (a - b))
                + frequency(field - h * (a + b))
            )
            / (4 * h**2)
            for b in unit
        ]
        for a in unit
    ]
    assert np.allclose(hessian, numeric_hessian, atol=1e-3 * np.max(np.abs(hessian)))


@pytest.mark.parametrize("zefoz", ER_YSO_SITE2_MATSUURA2026_ZEFOZ)
def test_matsuura_zefoz_points(zefoz):
    """The site-2 ZEFOZ points of Table V of Matsuura et al. follow from the spin Hamiltonian of their Table II."""
    i, j, frequency, magnitude, theta, phi, _ = zefoz
    h0, operators = wang_spin_hamiltonian()
    result = zefoz_newton(h0, operators, polar_field(magnitude, theta, phi), i, j)
    assert result["converged"]
    assert np.isclose(result["frequency"], frequency, atol=0.1)
    assert np.linalg.norm(result["field"] - polar_field(magnitude, theta, phi)) < 2e-3
    assert np.linalg.norm(result["gradient"]) < 1e-6


def test_zefoz_search_time_reversal():
    h0, operators = wang_spin_hamiltonian()
    start = polar_field(0.633, -37.53, -10.94)
    points = zefoz_search(h0, operators, [start, -start, 1.0001 * start], [(14, 15)])
    assert len(points) == 1 and points[0]["field"][2] >= 0 and points[0]["transition"] == (14, 15)


def test_run_dbqa_history_unitary():
    matrix = random_hermitian_matrix(3)
    groups = group_levels(np.linalg.eigvalsh(matrix), 1e-8)
    run = run_dbqa_history(
        matrix, "group commutator", dict(levels=groups), 30, np.linspace(1e-3, 1.0, 50),
        checkpoints=(30,), capture=("levels", 0, (1e-1, 1e-3)), track_unitary=True,
    )
    unitary = run["unitary"]
    assert np.allclose(unitary.conj().T @ unitary, np.eye(8))
    assert np.allclose(unitary.conj().T @ matrix @ unitary, run["checkpoints"][30])
    for threshold, (_, captured_matrix) in run["captured"].items():
        u = run["captured_unitary"][threshold]
        assert np.allclose(u.conj().T @ matrix @ u, captured_matrix)
    stopped = run_dbqa_history(
        matrix, "group commutator", dict(levels=groups), 30, np.linspace(1e-3, 1.0, 50),
        capture=("levels", 0, (1e-1,)), stop_when_captured=True,
    )
    assert len(stopped["steps"]) == stopped["captured"][1e-1][0]


def test_erhfs_basis_ordering(erhfs, er_yso_model):
    """`hyperfine_field_model` assumes electronic-major ordering with m_I = 7/2, ..., -7/2. Check it on the state
    labels, on the commutation of the electronic operators with I, and on the conservation of M_J + m_I by the
    magnetic-dipole hyperfine operator."""
    import re

    with open(pycf_matel_path("erhfs", "eryso") + ".st_") as file:
        lines = file.read().splitlines()[3:]
    labels = np.array(
        [re.match(r"\[4I\s*(\d+)\s*(-?\d+),\s*(-?\d+)\)", line).groups() for line in lines], dtype=int
    )
    assert labels.shape == (240, 3)
    assert np.all(labels[:, 2] == 7 - 2 * (np.arange(240) % 8))
    assert np.all(labels.reshape(30, 8, 3)[:, :, :2] == labels[::8, None, :2])
    i_ops = [np.kron(np.eye(30), i_a) for i_a in nuclear_spin_matrices(3.5)]
    for v in er_yso_model["electronic"]:
        assert all(np.allclose(v @ i_a, i_a @ v) for i_a in i_ops)
    hyp = pycf_dense_matrix([erhfs.HYP], {"HYP": 1.0})
    f_z = np.diag((labels[:, 1] + labels[:, 2]) / 2.0)
    assert np.allclose(hyp @ f_z, f_z @ hyp)
    assert np.allclose(er_yso_model["h0"], er_yso_model["h0"].conj().T)


def test_er_yso_high_field_hyperfine(er_yso_model):
    """Table VI of Jobbitt et al.: hyperfine splittings (MHz) of Z1 and Y1 of site 2 at 7 T along D1 from their
    crystal-field model, reproduced without the nuclear Zeeman term."""
    operators = field_operators(er_yso_model, nuclear_zeeman=False)
    energies = np.linalg.eigvalsh(linear_field_hamiltonian(er_yso_model["h0"], operators, [7.0, 0.0, 0.0]))
    z1 = np.diff(energies[:8]) * CM_TO_MHZ
    y1 = np.diff(energies[energies > 5000][:8]) * CM_TO_MHZ
    assert np.allclose(z1, [897, 881, 865, 849, 833, 817, 801], atol=1.5)
    assert np.allclose(y1, [928, 912, 895, 879, 863, 847, 831], atol=1.5)


def test_er_yso_g_values(f11cf):
    """Z1 and Y1 g values along D1, D2 and b in the full 4f11 basis, against Table II of Jobbitt et al."""
    coeff = er_yso_pycf_coeff()
    _, vectors = np.linalg.eigh(pycf_dense_matrix([getattr(f11cf, name) for name in coeff], coeff))
    rows = {r["state"]: r for r in er_yso_table()}
    for state, index in (("Z1", 0), ("Y1", 16)):
        u = vectors[:, index : index + 2]
        for axis, column in zip("XYZ", ("gD1_calc", "gD2_calc", "gb_calc")):
            moment = u.conj().T @ pycf_dense_matrix([getattr(f11cf, "MAG" + axis)], {"MAG" + axis: 1.0}) @ u
            assert np.isclose(np.ptp(np.linalg.eigvalsh(moment)), float(rows[state][column]), atol=0.25)


def test_er_yso_levels(f11cf):
    """The crystal-field splittings inside each multiplet reproduce Table II; the multiplet positions do not
    (see `Example_8_validation_Er_YSO.ipynb`)."""
    coeff = er_yso_pycf_coeff()
    levels = np.linalg.eigvalsh(pycf_dense_matrix([getattr(f11cf, name) for name in coeff], coeff))[::2]
    rows = er_yso_table()
    residuals = levels[: len(rows)] - np.array([float(r["E_calc"]) for r in rows])
    multiplets = [r["multiplet"] for r in rows]
    for name in dict.fromkeys(multiplets):
        assert np.ptp(residuals[[m == name for m in multiplets]]) < 4.5
    assert np.max(np.abs(residuals)) > 30


def test_avoided_crossings_two_level():
    """H(b) = c X + (b - b0) Z has one avoided crossing at b0 with gap 2c."""
    c, b0 = 0.05, 0.3
    x, z = PAULI_MATRICES["X"], PAULI_MATRICES["Z"]
    crossings = avoided_crossings(c * x - b0 * z, [z, 0 * z, 0 * z], [1.0, 0.0, 0.0], np.linspace(0, 1, 101), 2)
    assert len(crossings) == 1
    assert crossings[0]["levels"] == (0, 1)
    assert np.isclose(crossings[0]["field"], b0, atol=1e-6)
    assert np.isclose(crossings[0]["gap"], 2 * c, atol=1e-9)


def test_model3_single_sector(er_yso_model):
    """C1 crystal field, hyperfine and Zeeman terms couple every state: H(B) has no symmetry sectors, so the
    DB-DOI reference-state problem of Model 2 does not arise for Model 3."""
    from scipy.sparse.csgraph import connected_components

    for field in ([0.0, 0.0, 0.0], [-0.1482, 0.0083, 0.0208]):
        matrix = linear_field_hamiltonian(er_yso_model["h0"], field_operators(er_yso_model), field)
        n_components, _ = connected_components(np.abs(matrix) > 1e-12, directed=False)
        assert n_components == 1
