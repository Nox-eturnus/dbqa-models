import os
from typing import Optional

import numpy as np
from qibo.backends import _check_backend
from qibo.hamiltonians import Hamiltonian

CE_YAG_TENSORS = ["ZETA", "C20", "C22", "C40", "C42", "C44", "C60", "C62", "C64", "C66"]
"""PyCF ``f1cf`` tensors needed for Ce3+ at the D2 site of YAG (EAVG is added separately)."""

CE_YAG_HU2011 = {
    "fitted": dict(
        ZETA=647.0,
        C20=-465.0,
        C22=96.0,
        C40=-3739.0,
        C42=380.0,
        C44=1602.0,
        C60=901.0,
        C62=-307.0,
        C64=2136.0,
        C66=-246.0,
    ),
    "calculated": dict(
        ZETA=752.0,
        C20=-275.0,
        C22=449.0,
        C40=-2154.0,
        C42=-74.0,
        C44=1038.0,
        C60=842.0,
        C62=79.0,
        C64=1613.0,
        C66=-740.0,
    ),
}
"""4f spin-orbit and crystal-field parameters (cm^-1) of Ce3+:YAG from Table 10 of
L. Hu, M. F. Reid, C. K. Duan, S. Xia and M. Yin, J. Phys.: Condens. Matter 23, 045501 (2011).

``fitted`` is column [a] (fit of Tanner et al., J. Phys.: Condens. Matter 19, 216213 (2007),
rotated by pi/2 about z to the axis choice of Hu et al.); ``calculated`` is column [c] (DV-X alpha).
Both sets use the same axis convention, Wybourne normalisation, B^k_q = B^k_{-q} real (D2 site)."""

CE_YLF_PYCF_EXAMPLE = {
    "EAVG": 1035.1277,
    "ZETA": 625.6990,
    "C20": 297.8906,
    "C40": -1328.1522,
    "C44": -1282.4766,
    "C60": -191.5100,
    "C64": -1743.1424 + 692.8662j,
}
"""Ce3+:LiYF4 parameters (cm^-1) used in ``examples/ceylf/pyfit_example.py`` of PyCF."""


def pycf_matel_path(name: str, example: str, pycf_root: Optional[str] = None):
    """Path prefix of a matrix-element basis bundled with the PyCF examples, e.g. ``ceylf/matel/f1cf``.

    The PyCF source tree is found from `pycf_root` or from the ``PYCF_ROOT`` environment variable.
    """
    if pycf_root is None:
        pycf_root = os.environ.get("PYCF_ROOT")
    if pycf_root is None:
        raise ValueError(
            "Set `pycf_root` or the PYCF_ROOT environment variable to the PyCF source tree."
        )
    path = os.path.join(pycf_root, "examples", example, "matel", name)
    if not os.path.exists(path + ".txt"):
        raise FileNotFoundError(f"No PyCF matrix-element file found at {path}.txt")
    return path


def pycf_dense_matrix(tensors: list, coeff: dict):
    r"""Dense Hermitian matrix of the PyCF Hamiltonian $\sum_t c_t T_t$ in the PyCF basis.

    PyCF stores only the upper triangle of each tensor and completes the summed matrix
    hermitian-ly, so a complex coefficient $c_t$ multiplies $T^k_q$ and its conjugate partner
    $T^k_{-q}$ receives $c_t^*$. We follow the same assembly here.

    Args:
        tensors (list): `pycf.cfl.Tensor` objects (e.g. attributes of `ImportSLJM`).
        coeff (dict): coefficients keyed by tensor name, in cm^-1; tensors missing from `coeff` are skipped.
    """
    upper = sum(
        coeff[t.name] * np.triu(t.get_matel()) for t in tensors if t.name in coeff
    )
    strictly_upper = np.triu(upper, 1)
    return strictly_upper + strictly_upper.conj().T + np.diag(np.real(np.diag(upper)))


def embed_in_qubits(matrix: np.array, padding_offset: Optional[float] = None):
    r"""Embed a $d\times d$ Hermitian matrix into $2^L$ dimensions with $L=\lceil\log_2 d\rceil$.

    The $2^L-d$ padding states are uncoupled from the physical states and placed on the
    diagonal at $E_{\max} + $ `padding_offset`, i.e. moderately above the highest physical energy.

    Args:
        matrix (np.array): Hermitian matrix of the physical states.
        padding_offset (float, optional): distance of the padding level above the highest physical
            eigenvalue. Defaults to a quarter of the physical spectral width.

    Returns:
        (np.array, int, int): padded matrix, number of qubits, number of physical states.
    """
    n_physical = matrix.shape[0]
    nqubits = int(np.ceil(np.log2(n_physical)))
    energies = np.linalg.eigvalsh(matrix)
    if padding_offset is None:
        padding_offset = 0.25 * (energies[-1] - energies[0])
    padded = np.zeros((2**nqubits, 2**nqubits), dtype=complex)
    padded[:n_physical, :n_physical] = matrix
    for i in range(n_physical, 2**nqubits):
        padded[i, i] = energies[-1] + padding_offset
    return padded, nqubits, n_physical


def crystal_field_hamiltonian(
    tensors: list,
    coeff: dict,
    energy_unit: float = 1000.0,
    padding_offset: Optional[float] = None,
    backend=None,
):
    r"""Qibo `Hamiltonian` of a PyCF crystal-field model embedded in qubits.

    The matrix is expressed in units of `energy_unit` (cm^-1), so that its norm is O(1) like the
    spin-chain models used elsewhere in the pipeline and the default DBR step ranges apply.

    Returns:
        (`qibo.hamiltonians.Hamiltonian`, int): the embedded Hamiltonian and the number of physical states.
    """
    backend = _check_backend(backend)
    matrix = pycf_dense_matrix(tensors, coeff)
    if padding_offset is not None:
        padding_offset = padding_offset * energy_unit
    padded, nqubits, n_physical = embed_in_qubits(matrix, padding_offset)
    return (
        Hamiltonian(nqubits, backend.cast(padded / energy_unit), backend=backend),
        n_physical,
    )


def group_levels(energies: np.array, tol: float):
    """Group sorted `energies` into levels: consecutive values closer than `tol` share a level.

    Returns:
        list: lists of indices into the sorted `energies`, one list per level.
    """
    groups = [[0]]
    for i in range(1, len(energies)):
        if energies[i] - energies[i - 1] < tol:
            groups[-1].append(i)
        else:
            groups.append([i])
    return groups


def check_grouping_tolerance(energies: np.array, tol: float):
    """Largest splitting inside a level and smallest gap between levels for the grouping `tol`.

    The tolerance is acceptable when the first is smaller than `tol` and `tol` is smaller than
    the second, i.e. no splitting claimed to be resolved is merged.
    """
    groups = group_levels(energies, tol)
    max_inside = max(energies[g[-1]] - energies[g[0]] for g in groups)
    min_between = min(
        energies[groups[k + 1][0]] - energies[groups[k][-1]] for k in range(len(groups) - 1)
    )
    return max_inside, min_between


def subspace_fidelity(vectors_a: np.array, vectors_b: np.array):
    r"""Projector fidelity $\mathrm{Tr}(P_a P_b)/d$ between the spans of two sets of $d$ orthonormal columns."""
    overlap = vectors_a.conj().T @ vectors_b
    return np.real(np.trace(overlap @ overlap.conj().T)) / vectors_a.shape[1]


def gap_error(energies_approx: np.array, energies_exact: np.array, i: int, j: int):
    r"""$|\Delta E_{\rm DBQA} - \Delta E_{\rm exact}|$ for the gap between levels `i` and `j`."""
    return np.abs(
        (energies_approx[j] - energies_approx[i]) - (energies_exact[j] - energies_exact[i])
    )


def dbqa_spectrum_report(dbi_object, groups: list, energies_exact: np.array):
    r"""Energy errors and level fidelities of the current DBQA iterate.

    After $k$ rotations the pipeline stores $\hat H_k = \hat U^\dagger \hat H_0 \hat U$, so the
    approximate eigenvalues are the diagonal entries of $\hat H_k$ and the approximate eigenvectors
    are $\hat U|i\rangle$. The overlap $\langle v|\hat U|i\rangle$ with an exact eigenvector $|v\rangle$
    of $\hat H_0$ equals $\langle \hat U^\dagger v|i\rangle$, where $\hat U^\dagger|v\rangle$ is an exact
    eigenvector of $\hat H_k$. Fidelities are therefore computed from $\hat H_k$ alone, without
    tracking $\hat U$.

    Computational basis states are assigned to levels by ranking the diagonal of $\hat H_k$.

    Args:
        dbi_object (`DoubleBracketIteration`): the DBQA object, in the same energy units as `energies_exact`.
        groups (list): level grouping of the sorted exact spectrum, see `group_levels`.
        energies_exact (np.array): sorted exact eigenvalues of the embedded Hamiltonian.

    Returns:
        dict: ``energy_error`` (per state, sorted), ``level_energy_error`` (max per level),
        ``level_fidelity`` (projector fidelity per level) and ``off_diagonal_norm``.
    """
    h_matrix = dbi_object.backend.to_numpy(dbi_object.h.matrix)
    diagonal = np.real(np.diag(h_matrix))
    ranking = np.argsort(diagonal)
    energy_error = diagonal[ranking] - energies_exact
    _, vectors = np.linalg.eigh(h_matrix)
    identity = np.eye(len(diagonal))
    level_fidelity = [
        subspace_fidelity(vectors[:, g], identity[:, ranking[g]]) for g in groups
    ]
    level_energy_error = [np.max(np.abs(energy_error[g])) for g in groups]
    return dict(
        energy_error=energy_error,
        level_energy_error=np.array(level_energy_error),
        level_fidelity=np.array(level_fidelity),
        off_diagonal_norm=dbi_object.off_diagonal_norm,
    )
