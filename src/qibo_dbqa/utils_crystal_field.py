import os
from functools import reduce
from itertools import product
from typing import Optional

import numpy as np
from qibo import Circuit, gates, symbols
from qibo.backends import _check_backend
from qibo.hamiltonians import Hamiltonian, SymbolicHamiltonian
from scipy.linalg import expm

PAULI_MATRICES = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}

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
    level_energies = np.array([np.mean(energies_exact[g]) for g in groups])
    level_energies_dbqa = np.array([np.mean(diagonal[ranking[g]]) for g in groups])
    level_gap_error = [
        gap_error(level_energies_dbqa, level_energies, k, k + 1)
        for k in range(len(groups) - 1)
    ]
    return dict(
        energy_error=energy_error,
        level_energy_error=np.array(level_energy_error),
        level_fidelity=np.array(level_fidelity),
        level_gap_error=np.array(level_gap_error),
        off_diagonal_norm=dbi_object.off_diagonal_norm,
    )


def pauli_decomposition(matrix: np.array, tol: float = 1e-12):
    r"""Coefficients $c_P = \mathrm{Tr}(P H)/2^L$ of a Hermitian matrix in the Pauli basis.

    Pauli strings follow the Qibo ordering, i.e. the first character acts on qubit 0, the most
    significant bit of the computational basis index.

    Returns:
        dict: ``{"XIZY": c_P}`` for all strings with :math:`|c_P| >` `tol`.
    """
    nqubits = int(np.log2(matrix.shape[0]))
    decomposition = {}
    for string in product("IXYZ", repeat=nqubits):
        pauli = reduce(np.kron, [PAULI_MATRICES[c] for c in string])
        coefficient = np.real(np.trace(pauli @ matrix)) / 2**nqubits
        if abs(coefficient) > tol:
            decomposition["".join(string)] = coefficient
    return decomposition


def pauli_weight(string: str):
    """Number of qubits on which a Pauli string acts non-trivially."""
    return sum(c != "I" for c in string)


def add_pauli_rotation(circuit: Circuit, string: str, theta: float):
    r"""Append $e^{-i\theta P}$ for the Pauli string `string`, compiled into a CNOT ladder.

    A weight-$w$ string costs $2(w-1)$ CNOTs, one RZ and single-qubit basis changes.
    """
    support = [q for q, c in enumerate(string) if c != "I"]
    if len(support) == 0:
        return circuit
    basis_change = []
    for q in support:
        if string[q] == "X":
            basis_change.append(gates.H(q))
        elif string[q] == "Y":
            basis_change.append(gates.RX(q, np.pi / 2))
    ladder = [gates.CNOT(support[i], support[i + 1]) for i in range(len(support) - 1)]
    circuit.add(basis_change)
    circuit.add(ladder)
    circuit.add(gates.RZ(support[-1], 2 * theta))
    circuit.add(ladder[::-1])
    circuit.add([g.dagger() for g in basis_change])
    return circuit


def pauli_symbolic_form(decomposition: dict):
    """Symbolic Qibo form of a Pauli decomposition."""
    form = 0
    for string, coefficient in decomposition.items():
        term = coefficient
        for q, c in enumerate(string):
            if c != "I":
                term *= getattr(symbols, c)(q)
        form += term
    return form


class PauliCompiledHamiltonian(SymbolicHamiltonian):
    r"""Dense Hamiltonian written in the Pauli basis whose evolution is compiled into CNOTs.

    `circuit(dt)` returns a Trotter-Suzuki approximation of $e^{-i\,dt\,H}$ (order `ts_order`,
    `trotter_steps` steps) built from Pauli rotations, see `add_pauli_rotation`.
    The identity component `identity_coefficient` only contributes a global phase and is not compiled.
    This plays the role of `XXZ_compilation_line` for the crystal-field Hamiltonians.
    """

    def __init__(self, matrix, trotter_steps: int = 1, ts_order: int = 2, tol: float = 1e-12, backend=None):
        self.decomposition = pauli_decomposition(matrix, tol)
        nqubits = int(np.log2(matrix.shape[0]))
        super().__init__(pauli_symbolic_form(self.decomposition), nqubits=nqubits, backend=backend)
        self.trotter_steps = trotter_steps
        self.ts_order = ts_order
        self.pauli_terms = [
            (string, coefficient)
            for string, coefficient in self.decomposition.items()
            if pauli_weight(string) > 0
        ]
        self.identity_coefficient = self.decomposition.get("I" * nqubits, 0.0)

    def trotter_step(self, dt: float):
        circuit = Circuit(self.nqubits)
        if self.ts_order == 1:
            for string, coefficient in self.pauli_terms:
                add_pauli_rotation(circuit, string, coefficient * dt)
        elif self.ts_order == 2:
            for string, coefficient in self.pauli_terms[:-1]:
                add_pauli_rotation(circuit, string, coefficient * dt / 2)
            string, coefficient = self.pauli_terms[-1]
            add_pauli_rotation(circuit, string, coefficient * dt)
            for string, coefficient in self.pauli_terms[-2::-1]:
                add_pauli_rotation(circuit, string, coefficient * dt / 2)
        else:
            raise ValueError(f"Trotter-Suzuki order {self.ts_order} is not supported.")
        return circuit

    def circuit(self, dt, accelerators=None):
        return reduce(
            Circuit.__add__,
            [self.trotter_step(dt / self.trotter_steps) for _ in range(self.trotter_steps)],
        )


class ExactEvolutionHamiltonian(SymbolicHamiltonian):
    r"""Same Hamiltonian as `PauliCompiledHamiltonian` but `circuit(dt)` is the exact $e^{-i\,dt\,H}$
    as a single `Unitary` gate, i.e. an ideal evolution oracle used to separate DBQA convergence
    from Trotter errors."""

    def __init__(self, matrix, tol: float = 1e-12, backend=None):
        self.dense_matrix = np.array(matrix, dtype=complex)
        nqubits = int(np.log2(matrix.shape[0]))
        super().__init__(
            pauli_symbolic_form(pauli_decomposition(matrix, tol)), nqubits=nqubits, backend=backend
        )

    def circuit(self, dt, accelerators=None):
        circuit = Circuit(self.nqubits)
        circuit.add(gates.Unitary(expm(-1j * dt * self.dense_matrix), *range(self.nqubits)))
        return circuit


def two_qubit_gate_count(circuit: Circuit):
    """Number of gates acting on exactly two qubits."""
    return sum(1 for gate in circuit.queue if len(gate.qubits) == 2)


def basis_state_circuit(nqubits: int, index: int):
    r"""Circuit preparing the computational basis state $|$`index`$\rangle$ from $|0\ldots0\rangle$."""
    circuit = Circuit(nqubits)
    bits = format(index, f"0{nqubits}b")
    circuit.add([gates.X(q) for q, b in enumerate(bits) if b == "1"])
    return circuit
