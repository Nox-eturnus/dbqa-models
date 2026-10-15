"""Testing the PyCF -> Qibo crystal-field interface used for DBQA"""

import os

import numpy as np
import pytest

pycf_cfl = pytest.importorskip("pycf.cfl")
import_sljm = pytest.importorskip("pycf.import_sljm")

from qibo_dbqa.double_bracket import (
    DoubleBracketGeneratorType,
    DoubleBracketIteration,
)
from qibo_dbqa.utils_crystal_field import *

NSTEPS = 3
TOL = 1e-8

if os.environ.get("PYCF_ROOT") is None:  # pragma: no cover
    pytest.skip("PYCF_ROOT is not set", allow_module_level=True)


@pytest.fixture(scope="module")
def f1cf():
    return import_sljm.ImportSLJM(pycf_matel_path("f1cf", "ceylf"))


def ce_yag_coeff(parameter_set="fitted"):
    return dict(CE_YAG_HU2011[parameter_set], EAVG=0.0)


def test_dense_matrix_matches_pycf(f1cf):
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
