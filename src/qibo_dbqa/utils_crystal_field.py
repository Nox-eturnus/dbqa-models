import os
from functools import reduce
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


ER_LIYF4_COUTO1998 = dict(
    E0=35214.19,
    E1=6607.22,
    E2=32.99,
    E3=664.86,
    ALPHA=18.29,
    BETA=-589.84,
    GAMMA=1755.61,
    T2=455.12,
    T3=38.43,
    T4=83.67,
    T6=-348.26,
    T7=243.92,
    T8=388.37,
    ZETA=2378.00,
    M0=4.10,
    P2=801.82,
    C20=323.07,
    C40=-753.83,
    C44=-1052.24,
    C60=-43.60,
    C64=-623.59,
)
"""Free-ion and crystal-field parameters (cm^-1) of Er3+:LiYF4 from Table 2 of M. A. Couto dos Santos,
E. Antic-Fidancev, J. Y. Gesland, J. C. Krupa, M. Lemaitre-Blaise and P. Porcher, J. Alloys Compd. 275-277,
435-441 (1998), D2d approximation of the S4 site. The Coulomb interaction is given by the Racah parameters
E^0..E^3; `racah_to_slater` converts them to the Slater-Condon F^k used by PyCF. M2, M4, P4 and P6 follow from
M0 and P2 through the fixed ratios `MP_RATIOS` (Table 2 lists 2.30, 1.27, 400.90 and 80.18).

These parameters do not reproduce the paper's levels in PyCF: the paper's E^3 follows a different convention
(see `Example_4_validation_Er_LiYF4.ipynb`). `ER_LIYF4_FREE_ION_REFIT` keeps the crystal field and refits the
free-ion parameters to the experimental levels."""

ER_LIYF4_FREE_ION_REFIT = dict(
    ER_LIYF4_COUTO1998,
    E1=6583.09,
    E2=32.91,
    E3=625.12,
    ALPHA=18.08,
    BETA=-596.69,
    GAMMA=1892.51,
    T2=372.53,
    T3=41.74,
    T4=71.19,
    T6=-353.31,
    T7=234.42,
    T8=382.67,
    ZETA=2375.21,
    M0=3.76,
    P2=721.45,
)
"""Model 2 Hamiltonian: crystal field of Couto dos Santos et al. unchanged, free-ion parameters refitted in PyCF
conventions to their 108 experimental levels (sigma = 12.9 cm^-1 for 21 parameters, paper: 14.6 cm^-1),
see `Example_4_validation_Er_LiYF4.ipynb`."""

MP_RATIOS = dict(M2=0.56, M4=0.31, P4=0.50, P6=0.10)
"""Ratios M^2/M^0, M^4/M^0, P^4/P^2 and P^6/P^2 kept fixed by Couto dos Santos et al. (from Carnall et al.)."""

MU_B = 0.46686447783
"""Bohr magneton in cm^-1/T (CODATA 2018)."""

MU_N = 2.542623413e-4
"""Nuclear magneton in cm^-1/T (CODATA 2018)."""

CM_TO_MHZ = 29979.2458
"""1 cm^-1 in MHz."""

ER167_NUCLEAR_SPIN = 3.5
ER167_GN = -0.1618
"""Nuclear spin and nuclear g-factor of 167Er (Table I of Matsuura et al., Phys. Rev. B 113, 085421 (2026))."""

ER_YSO_CARNALL_FIXED = dict(
    ALPHA=17.79,
    BETA=-582.10,
    GAMMA=1800.00,
    T2=400.00,
    T3=43.00,
    T4=73.00,
    T6=-271.00,
    T7=308.00,
    T8=299.00,
    M0=3.86,
    P2=594.00,
)
"""Free-ion parameters (cm^-1) held fixed by Jobbitt et al. (their Table IV, from Carnall et al., Er3+:LaF3)."""

ER_YSO_SITE2_JOBBITT2021 = dict(
    ER_YSO_CARNALL_FIXED,
    EAVG=35507.5,
    F2=96121.9,
    F4=67722.4,
    F6=53241.2,
    ZETA=2362.3,
    C20=389.0,
    C21=-325.7 - 95.8j,
    C22=-368.5 + 53.7j,
    C40=17.2,
    C41=-378.7 - 519.5j,
    C42=-72.0 - 146.0j,
    C43=-890.8 + 570.4j,
    C44=-198.7 - 567.9j,
    C60=73.4,
    C61=-37.5 + 49.9j,
    C62=135.5 + 60.6j,
    C63=-166.7 + 131.8j,
    C64=227.2 + 47.6j,
    C65=119.5 + 64.3j,
    C66=37.6 - 41.3j,
)
"""Free-ion and crystal-field parameters (cm^-1) of Er3+ at site 2 of Y2SiO5 from Table III ("This study") of
N. L. Jobbitt, J.-P. R. Wells, M. F. Reid, S. P. Horvath, P. Goldner and A. Ferrier, Phys. Rev. B 104, 155121 (2021),
arXiv:2206.09080. Axes x, y, z = D1, D2, b; C1 site, Wybourne normalisation. The other magnetically inequivalent
orientation follows by changing the sign of every odd-q parameter. M2, M4, P4 and P6 follow from M0 and P2 through
`MP_RATIOS`."""

ER_YSO_SITE2_HYPERFINE_JOBBITT2021 = dict(HYP=0.005389, EQHYP=0.0240)
"""Magnetic-dipole (a_l) and electric-quadrupole (a_Q) hyperfine parameters (cm^-1) of 167Er3+ at site 2,
Table III of Jobbitt et al.; they multiply PyCF's ``HYP`` and ``EQHYP`` tensors of the ``erhfs`` basis."""

ER_YSO_SITE2_WANG_SPIN_HAMILTONIAN = dict(
    A=np.array([[-1570.0, 224.0, -131.0], [224.0, -17.0, -15.0], [-131.0, -15.0, 141.0]]),
    Q=np.array([[-9.8, -21.0, -0.4], [-21.0, -16.0, -12.4], [-0.4, -12.4, 25.8]]),
    g=np.array([[14.44, -1.76, 2.35], [-1.76, 1.91, -0.46], [2.35, -0.46, 1.424]]),
)
"""Effective spin Hamiltonian of the site-2 ground doublet 4I15/2 Z1 (A and Q in MHz, frame D1, D2, b) used by
K. Matsuura, S. Yasui, R. Kaji, H. Sasakura, T. Tawara and S. Adachi, Phys. Rev. B 113, 085421 (2026),
arXiv:2412.10126 (their Table II, parameters of S.-J. Wang et al.)."""

ER_YSO_SITE2_MATSUURA2026_ZEFOZ = [
    # (i, j, frequency in MHz, |B| in T, theta, phi in degrees, T2 in s)
    (14, 15, 796.3, 0.633, -37.53, -10.94, 95.3),
    (6, 7, 797.6, 0.668, -38.16, -12.59, 89.5),
    (13, 14, 796.0, 0.413, -40.22, -11.85, 54.3),
    (13, 15, 1592.2, 0.510, -38.88, -11.42, 36.0),
    (12, 13, 798.0, 0.283, -42.23, -12.58, 30.4),
]
"""The five site-2 ZEFOZ transitions with the longest coherence time in Table V of Matsuura et al. (orientation 1).
Levels are numbered from the lowest of the 16 Z1 hyperfine states; the field is given by `polar_field`. The
Cartesian field of their optimal point printed in the text, (-378.9, 73.2, -502.3) mT, has the sign of the b
component reversed with respect to their polar angles; the polar form reproduces their Table V."""


def racah_to_slater(e1: float, e2: float, e3: float):
    r"""Slater-Condon parameters $F^2, F^4, F^6$ from the Racah parameters $E^1, E^2, E^3$ of an $f^N$ configuration.

    Uses $E^1 = (70F_2 + 231F_4 + 2002F_6)/9$, $E^2 = (F_2 - 3F_4 + 7F_6)/9$, $E^3 = (5F_2 + 6F_4 - 91F_6)/3$
    with the reduced parameters $F_2 = F^2/225$, $F_4 = F^4/1089$, $F_6 = 25F^6/184041$ (Racah; Wybourne,
    *Spectroscopic Properties of Rare Earths*, 1965).
    """
    racah = np.array(
        [[70 / 9, 231 / 9, 2002 / 9], [1 / 9, -3 / 9, 7 / 9], [5 / 3, 6 / 3, -91 / 3]]
    )
    reduced = np.linalg.solve(racah, np.array([e1, e2, e3]))
    return reduced * np.array([225.0, 1089.0, 184041.0 / 25.0])


def slater_to_racah(f2: float, f4: float, f6: float):
    r"""Racah parameters $E^1, E^2, E^3$ from the Slater-Condon $F^2, F^4, F^6$; inverse of `racah_to_slater`."""
    reduced = np.array([f2 / 225.0, f4 / 1089.0, 25.0 * f6 / 184041.0])
    return np.array(
        [
            (70 * reduced[0] + 231 * reduced[1] + 2002 * reduced[2]) / 9,
            (reduced[0] - 3 * reduced[1] + 7 * reduced[2]) / 9,
            (5 * reduced[0] + 6 * reduced[1] - 91 * reduced[2]) / 3,
        ]
    )


def er_liyf4_pycf_coeff(parameters: dict = None):
    """PyCF coefficients (cm^-1) for the `f11cf` basis from the Couto dos Santos et al. parameters.

    E^0 only shifts the whole configuration and is replaced by ``EAVG = 0``; energies are compared relative to
    a fitted constant offset.
    """
    if parameters is None:
        parameters = ER_LIYF4_COUTO1998
    f2, f4, f6 = racah_to_slater(parameters["E1"], parameters["E2"], parameters["E3"])
    coeff = {name: value for name, value in parameters.items() if name not in ("E0", "E1", "E2", "E3")}
    coeff.update(EAVG=0.0, F2=f2, F4=f4, F6=f6)
    coeff.update(
        M2=MP_RATIOS["M2"] * parameters["M0"],
        M4=MP_RATIOS["M4"] * parameters["M0"],
        P4=MP_RATIOS["P4"] * parameters["P2"],
        P6=MP_RATIOS["P6"] * parameters["P2"],
    )
    return coeff


def er_yso_pycf_coeff(parameters: dict = None, hyperfine: Optional[dict] = None):
    """PyCF coefficients (cm^-1) for Er3+:Y2SiO5: `parameters` (default site 2 of Jobbitt et al.) with M2, M4, P4
    and P6 from `MP_RATIOS`, plus the hyperfine coefficients `hyperfine` if given (``erhfs`` basis only)."""
    if parameters is None:
        parameters = ER_YSO_SITE2_JOBBITT2021
    coeff = dict(parameters)
    coeff.update(
        M2=MP_RATIOS["M2"] * parameters["M0"],
        M4=MP_RATIOS["M4"] * parameters["M0"],
        P4=MP_RATIOS["P4"] * parameters["P2"],
        P6=MP_RATIOS["P6"] * parameters["P2"],
    )
    if hyperfine is not None:
        coeff.update(hyperfine)
    return coeff


def crystal_field_strength(coeff: dict, k: int):
    r"""Crystal-field strength $S_k = [(|B^k_0|^2 + 2\sum_{q>0}|B^k_q|^2)/(2k+1)]^{1/2}$ (Eq. (5) of Jobbitt et al.)."""
    total = abs(coeff.get(f"C{k}0", 0.0)) ** 2 + 2 * sum(
        abs(coeff.get(f"C{k}{q}", 0.0)) ** 2 for q in range(1, k + 1)
    )
    return np.sqrt(total / (2 * k + 1))


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


def spectral_grid_search_step(dbi_object, space: np.array, d: Optional[np.array] = None, mode=None):
    r"""Grid search of the DBR duration with the same result as `grid_search_step` but faster for large matrices.

    Within one step the generator is fixed, so the rotation is diagonalised once and every duration in
    `space` only costs a few matrix products:
    canonical and single-commutator rotations $e^{-sW}$ use the eigenbasis of $iW$, the group commutator
    uses the eigenbasis of $\hat H_k$ and the diagonal $D_k$.
    The off-diagonal norm of $U^\dagger \hat H_k U$ is minimised, as in `DoubleBracketIteration.loss`.

    Returns:
        (float, float): optimal duration and the corresponding off-diagonal norm.
    """
    from .double_bracket import DoubleBracketGeneratorType  # pylint: disable=import-outside-toplevel

    if mode is None:
        mode = dbi_object.mode
    backend = dbi_object.backend
    h = backend.to_numpy(dbi_object.h.matrix)
    if d is None:
        d = np.diag(np.diag(h))
    d = backend.to_numpy(backend.cast(d))
    if mode in (DoubleBracketGeneratorType.canonical, DoubleBracketGeneratorType.single_commutator):
        if mode is DoubleBracketGeneratorType.canonical:
            d = np.diag(np.diag(h))
        w = d @ h - h @ d
        a, v = np.linalg.eigh(1j * w)

        def unitary(s):
            return (v * np.exp(1j * s * a)) @ v.conj().T

    elif mode is DoubleBracketGeneratorType.group_commutator:
        e, v = np.linalg.eigh(h)
        d_diag = np.real(np.diag(d))

        def unitary(s):
            h_exp_plus = (v * np.exp(-1j * s * e)) @ v.conj().T
            h_exp_minus = (v * np.exp(1j * s * e)) @ v.conj().T
            return (h_exp_plus * np.exp(1j * s * d_diag)) @ (h_exp_minus * np.exp(-1j * s * d_diag))

    else:  # pragma: no cover
        raise ValueError(f"Mode {mode} is not supported by spectral_grid_search_step.")
    # ||sigma(U^dag H U)||^2 = ||H||^2 - sum_i |(U^dag H U)_ii|^2, the Frobenius norm being unitarily invariant
    norm_squared = np.linalg.norm(h) ** 2
    losses = []
    for s in space:
        u = unitary(s)
        diagonal = np.real(np.sum(u.conj() * (h @ u), axis=0))
        losses.append(np.sqrt(max(norm_squared - np.sum(diagonal**2), 0.0)))
    index = int(np.argmin(losses))
    return space[index], losses[index]


def run_dbqa_history(
    matrix: np.array,
    method: str,
    group_sets: dict,
    nsteps: int,
    space: np.array,
    d_list: Optional[list] = None,
    checkpoints: tuple = (),
    backend: str = "numpy",
    capture: Optional[tuple] = None,
    track_unitary: bool = False,
    stop_when_captured: bool = False,
):
    """Run `nsteps` DBQA iterations on `matrix` with spectral grid-search scheduling and record the spectrum report.

    Args:
        matrix (np.array): embedded Hamiltonian (dimensionless).
        method (str): ``"canonical"``, ``"group commutator"`` (with D = Delta(H_k)) or ``"Pauli-Z"``
            (single commutator with the best D from `d_list`; negative durations cover the sign of D).
        group_sets (dict): named level groupings, e.g. ``{"doublets": ..., "manifolds": ...}``; a report is
            recorded for each of them at every iteration.
        space (np.array): durations tried at each step.
        checkpoints (tuple): iterations at which the iterated matrix is stored.
        backend (str): name of the Qibo backend, set explicitly so that the function can run in worker processes.
        capture (tuple, optional): ``(grouping, index, thresholds)``; for every threshold, the first iterated
            matrix with ``1 - level_fidelity[index] < threshold`` in that grouping is stored.
        track_unitary (bool): also accumulate the total rotation $U$, $\\hat H_k = U^\\dagger \\hat H_0 U$, e.g. to
            warm-start the next point of a parameter sweep from $U^\\dagger \\hat H_0' U$.
        stop_when_captured (bool): stop as soon as every threshold of `capture` has been reached.

    Returns:
        dict: ``history`` (list of reports per grouping), ``steps``, ``generators``, ``checkpoints`` and
        ``captured`` (``{threshold: (iteration, matrix)}``); with `track_unitary` also ``unitary`` (after the
        last iteration) and ``captured_unitary`` (``{threshold: U}``).
    """
    from qibo.backends import construct_backend  # pylint: disable=import-outside-toplevel

    from .double_bracket import (  # pylint: disable=import-outside-toplevel
        DoubleBracketGeneratorType,
        DoubleBracketIteration,
    )

    backend = construct_backend(backend)

    nqubits = int(np.log2(matrix.shape[0]))
    energies = np.linalg.eigvalsh(matrix)
    mode = {
        "canonical": DoubleBracketGeneratorType.canonical,
        "group commutator": DoubleBracketGeneratorType.group_commutator,
        "Pauli-Z": DoubleBracketGeneratorType.single_commutator,
    }[method]
    dbi = DoubleBracketIteration(Hamiltonian(nqubits, backend.cast(np.array(matrix)), backend=backend), mode=mode)

    def report():
        return {name: dbqa_spectrum_report(dbi, groups, energies) for name, groups in group_sets.items()}

    captured, captured_unitary = {}, {}
    unitary = np.eye(matrix.shape[0], dtype=complex)

    def capture_matrix(k, reports):
        if capture is None:
            return
        grouping, index, thresholds = capture
        infidelity = 1 - reports[grouping]["level_fidelity"][index]
        for threshold in thresholds:
            if threshold not in captured and infidelity < threshold:
                captured[threshold] = (k, backend.to_numpy(dbi.h.matrix))
                if track_unitary:
                    captured_unitary[threshold] = unitary.copy()

    def all_captured():
        return stop_when_captured and capture is not None and len(captured) == len(capture[2])

    history, steps, generators, stored = [report()], [], [], {}
    capture_matrix(0, history[0])
    if 0 in checkpoints:
        stored[0] = backend.to_numpy(dbi.h.matrix)
    for k in range(1, nsteps + 1):
        if all_captured():
            break
        if method == "Pauli-Z":
            signed_space = np.concatenate([-space[::-1], space])
            best = None
            for j, d in enumerate(d_list):
                step, loss = spectral_grid_search_step(dbi, signed_space, d=d)
                if best is None or loss < best[1]:
                    best = (step, loss, j)
            step, _, j = best
            operator = dbi(step=step, d=d_list[j])
            generators.append(j)
        else:
            step, _ = spectral_grid_search_step(dbi, space)
            operator = dbi(step=step)
        if track_unitary:
            unitary = unitary @ backend.to_numpy(operator)
        steps.append(step)
        history.append(report())
        capture_matrix(k, history[-1])
        if k in checkpoints:
            stored[k] = backend.to_numpy(dbi.h.matrix)
    result = dict(history=history, steps=steps, generators=generators, checkpoints=stored, captured=captured)
    if track_unitary:
        result.update(unitary=unitary, captured_unitary=captured_unitary)
    return result


def pauli_decomposition(matrix: np.array, tol: float = 1e-12):
    r"""Coefficients $c_P = \mathrm{Tr}(P H)/2^L$ of a Hermitian matrix in the Pauli basis.

    Pauli strings follow the Qibo ordering, i.e. the first character acts on qubit 0, the most
    significant bit of the computational basis index.

    Returns:
        dict: ``{"XIZY": c_P}`` for all strings with :math:`|c_P| >` `tol`.
    """
    nqubits = int(np.log2(matrix.shape[0]))
    # interleave row and column bits so that axis q indexes (a_q, b_q) of qubit q
    tensor = np.asarray(matrix, dtype=complex).reshape([2] * (2 * nqubits))
    tensor = tensor.transpose([i for q in range(nqubits) for i in (q, nqubits + q)])
    tensor = tensor.reshape([4] * nqubits)
    # kernel[p, 2a + b] = P_p[b, a] / 2, so that contracting every axis gives Tr(P H) / 2^L
    kernel = np.array([PAULI_MATRICES[c].T.reshape(4) for c in "IXYZ"]) / 2
    for axis in range(nqubits):
        tensor = np.moveaxis(np.tensordot(kernel, tensor, axes=([1], [axis])), 0, axis)
    coefficients = np.real(tensor).reshape(-1)
    decomposition = {}
    for index in np.nonzero(np.abs(coefficients) > tol)[0]:
        string = "".join("IXYZ"[int(d)] for d in np.base_repr(index, 4).zfill(nqubits))
        decomposition[string] = coefficients[index]
    return decomposition


def trotter_cnot_count(decomposition: dict, ts_order: int = 1):
    """CNOTs of one Trotter-Suzuki step compiled with `add_pauli_rotation` (no gate cancellations)."""
    cnots = [2 * (pauli_weight(s) - 1) for s in decomposition if pauli_weight(s) > 0]
    if ts_order == 1:
        return sum(cnots)
    return 2 * sum(cnots) - cnots[-1]


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


class DenseSymbolicHamiltonian(SymbolicHamiltonian):
    """`SymbolicHamiltonian` whose matrix and expectation values come from a stored dense matrix.

    `GroupCommutatorIteration` expects a `SymbolicHamiltonian` with a custom `circuit(dt)`; building the
    symbolic form of a dense crystal-field Hamiltonian (up to 4^L Pauli strings) is too slow at 9 qubits,
    so the subclasses below only keep the dense matrix.
    """

    def __init__(self, matrix, backend=None):
        self.dense_matrix = np.array(matrix, dtype=complex)
        nqubits = int(np.log2(self.dense_matrix.shape[0]))
        super().__init__(symbols.I(0), nqubits=nqubits, backend=backend)

    @property
    def matrix(self):
        return self.backend.cast(self.dense_matrix)

    def expectation(self, state, normalize=False):
        state = self.backend.to_numpy(state)
        value = np.real(np.conj(state) @ self.dense_matrix @ state)
        if normalize:
            value = value / np.real(np.conj(state) @ state)
        return value


class PauliCompiledHamiltonian(DenseSymbolicHamiltonian):
    r"""Dense Hamiltonian written in the Pauli basis whose evolution is compiled into CNOTs.

    `circuit(dt)` returns a Trotter-Suzuki approximation of $e^{-i\,dt\,H}$ (order `ts_order`,
    `trotter_steps` steps) built from Pauli rotations, see `add_pauli_rotation`.
    The identity component `identity_coefficient` only contributes a global phase and is not compiled.
    This plays the role of `XXZ_compilation_line` for the crystal-field Hamiltonians.
    """

    def __init__(self, matrix, trotter_steps: int = 1, ts_order: int = 2, tol: float = 1e-12, backend=None):
        super().__init__(matrix, backend=backend)
        self.decomposition = pauli_decomposition(matrix, tol)
        self.trotter_steps = trotter_steps
        self.ts_order = ts_order
        self.pauli_terms = [
            (string, coefficient)
            for string, coefficient in self.decomposition.items()
            if pauli_weight(string) > 0
        ]
        self.identity_coefficient = self.decomposition.get("I" * self.nqubits, 0.0)

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


class ExactEvolutionHamiltonian(DenseSymbolicHamiltonian):
    r"""Same Hamiltonian as `PauliCompiledHamiltonian` but `circuit(dt)` is the exact $e^{-i\,dt\,H}$
    as a single `Unitary` gate, i.e. an ideal evolution oracle used to separate DBQA convergence
    from Trotter errors. Unitaries are cached by duration."""

    def __init__(self, matrix, backend=None):
        super().__init__(matrix, backend=backend)
        self._cache = {}

    def unitary(self, dt):
        key = round(float(dt), 15)
        if key not in self._cache:
            self._cache[key] = expm(-1j * dt * self.dense_matrix)
        return self._cache[key]

    def circuit(self, dt, accelerators=None):
        circuit = Circuit(self.nqubits)
        circuit.add(gates.Unitary(self.unitary(dt), *range(self.nqubits)))
        return circuit


class CompiledUnitaryHamiltonian(DenseSymbolicHamiltonian):
    r"""Evolution oracle given by the exact unitary of the compiled Trotter circuit, as one `Unitary` gate.

    `circuit(dt)` applies the same operator as `PauliCompiledHamiltonian(matrix, trotter_steps, ts_order).circuit(dt)`
    up to the global phase of the identity component, computed with `compiled_trotter_unitary` instead of
    simulating millions of gates. Unitaries are cached by duration.
    """

    def __init__(self, matrix, trotter_steps: int = 1, ts_order: int = 2, tol: float = 1e-12, backend=None):
        super().__init__(matrix, backend=backend)
        self.decomposition = pauli_decomposition(matrix, tol)
        self.trotter_steps = trotter_steps
        self.ts_order = ts_order
        self._cache = {}

    def unitary(self, dt):
        key = round(float(dt), 15)
        if key not in self._cache:
            step = compiled_trotter_unitary(self.decomposition, dt / self.trotter_steps, self.ts_order)
            self._cache[key] = np.linalg.matrix_power(step, self.trotter_steps)
        return self._cache[key]

    def circuit(self, dt, accelerators=None):
        circuit = Circuit(self.nqubits)
        circuit.add(gates.Unitary(self.unitary(dt), *range(self.nqubits)))
        return circuit


def apply_pauli_rotation(matrix: np.array, string: str, theta: float):
    r"""Return $e^{-i\theta P}$ `matrix` for the Pauli string `string` (qubit 0 = most significant bit).

    Matrix-free: $P|b\rangle = i^{n_Y} (-1)^{|b \wedge z|} |b \oplus x\rangle$ with $x$ the X/Y positions and
    $z$ the Z/Y positions, so $e^{-i\theta P} = \cos\theta - i \sin\theta\, P$ costs one permutation and one
    phase per term. This is exactly the operator implemented by `add_pauli_rotation`.
    """
    nqubits = len(string)
    index = np.arange(2**nqubits)
    x_mask = z_mask = 0
    n_y = 0
    for q, c in enumerate(string):
        bit = 1 << (nqubits - 1 - q)
        if c in "XY":
            x_mask |= bit
        if c in "ZY":
            z_mask |= bit
        n_y += c == "Y"
    masked = index & z_mask
    parity = np.zeros(len(index), dtype=int)
    for q in range(nqubits):
        parity ^= (masked >> q) & 1
    phase = (1j**n_y) * (1 - 2 * parity)
    pauli_matrix = np.empty_like(matrix)
    pauli_matrix[index ^ x_mask] = phase[:, None] * matrix if matrix.ndim == 2 else phase * matrix
    return np.cos(theta) * matrix - 1j * np.sin(theta) * pauli_matrix


def compiled_trotter_unitary(decomposition: dict, dt: float, ts_order: int = 1, matrix: np.array = None):
    r"""Exact unitary of one Trotter-Suzuki step compiled with `add_pauli_rotation`, without building the circuit.

    The terms are applied in the order of `PauliCompiledHamiltonian.trotter_step`; the identity term is omitted
    (global phase). If `matrix` is given, the step is applied to it instead of to the identity.
    """
    terms = [(s, c) for s, c in decomposition.items() if pauli_weight(s) > 0]
    nqubits = len(terms[0][0])
    result = np.eye(2**nqubits, dtype=complex) if matrix is None else np.array(matrix, dtype=complex)
    if ts_order == 1:
        sequence = [(s, c * dt) for s, c in terms]
    elif ts_order == 2:
        sequence = (
            [(s, c * dt / 2) for s, c in terms[:-1]]
            + [(terms[-1][0], terms[-1][1] * dt)]
            + [(s, c * dt / 2) for s, c in terms[-2::-1]]
        )
    else:
        raise ValueError(f"Trotter-Suzuki order {ts_order} is not supported.")
    for string, theta in sequence:
        result = apply_pauli_rotation(result, string, theta)
    return result


def compiled_trotter_depth(decomposition: dict, ts_order: int = 1, trotter_steps: int = 1):
    """Circuit depth of `PauliCompiledHamiltonian(...).circuit(dt)` computed without building the circuit.

    Gates are scheduled as soon as all their qubits are free (as for `Circuit.depth`), using the gate sequence
    of `add_pauli_rotation`: basis changes, CNOT ladder, RZ, reversed ladder, inverse basis changes.
    """
    terms = [s for s in decomposition if pauli_weight(s) > 0]
    nqubits = len(terms[0])
    if ts_order == 1:
        sequence = terms
    elif ts_order == 2:
        sequence = terms[:-1] + [terms[-1]] + terms[-2::-1]
    else:
        raise ValueError(f"Trotter-Suzuki order {ts_order} is not supported.")
    free = [0] * nqubits

    def schedule(*qubits):
        moment = max(free[q] for q in qubits)
        for q in qubits:
            free[q] = moment + 1

    for _ in range(trotter_steps):
        for string in sequence:
            support = [q for q, c in enumerate(string) if c != "I"]
            changed = [q for q in support if string[q] in "XY"]
            for q in changed:
                schedule(q)
            for i in range(len(support) - 1):
                schedule(support[i], support[i + 1])
            schedule(support[-1])
            for i in reversed(range(len(support) - 1)):
                schedule(support[i], support[i + 1])
            for q in changed:
                schedule(q)
    return max(free)


def two_qubit_gate_count(circuit: Circuit):
    """Number of gates acting on exactly two qubits."""
    return sum(1 for gate in circuit.queue if len(gate.qubits) == 2)


def basis_state_circuit(nqubits: int, index: int):
    r"""Circuit preparing the computational basis state $|$`index`$\rangle$ from $|0\ldots0\rangle$."""
    circuit = Circuit(nqubits)
    bits = format(index, f"0{nqubits}b")
    circuit.add([gates.X(q) for q, b in enumerate(bits) if b == "1"])
    return circuit


def nuclear_spin_matrices(spin: float):
    r"""Spin matrices $I_x, I_y, I_z$ for nuclear spin `spin`, in the basis $m_I = I, I-1, \dots, -I$."""
    m = np.arange(spin, -spin - 1, -1)
    raising = np.diag(np.sqrt(spin * (spin + 1) - m[1:] * (m[1:] + 1)), 1)
    return [(raising + raising.T) / 2, (raising - raising.T) / 2j, np.diag(m).astype(complex)]


def hyperfine_field_model(basis, coeff: dict, nuclear_spin: float = ER167_NUCLEAR_SPIN, gn: float = ER167_GN):
    r"""Field-independent Hamiltonian and field operators of a PyCF hyperfine basis such as ``erhfs``.

    The Hamiltonian is linear in the magnetic field $\vec B$ (in tesla):
    $H(\vec B) = H_0 + \sum_a B_a (V^{\rm el}_a + V^{\rm nuc}_a)$, with the electronic Zeeman operators
    $V^{\rm el}_a = \mu_B (L + g_s S)_a$ (PyCF ``MAGX``, ``MAGY``, ``MAGZ``) and the nuclear Zeeman operators
    $V^{\rm nuc}_a = -g_n \mu_N I_a$. The basis is ordered electronic state first, then $m_I = I, \dots, -I$, as in
    PyCF's ``erhfs`` files; `tests/test_models_crystal_field.py` checks this ordering.

    Args:
        basis: `ImportSLJM` object of the hyperfine basis.
        coeff (dict): PyCF coefficients in cm^-1, including ``HYP`` and ``EQHYP``.

    Returns:
        dict: ``h0`` (cm^-1), ``electronic`` and ``nuclear`` (lists of three matrices in cm^-1/T, axes x, y, z).
    """
    tensors = [getattr(basis, name) for name in coeff]
    h0 = pycf_dense_matrix(tensors, coeff)
    n_nuclear = int(round(2 * nuclear_spin + 1))
    n_electronic = h0.shape[0] // n_nuclear
    electronic = [
        MU_B * pycf_dense_matrix([getattr(basis, "MAG" + axis)], {"MAG" + axis: 1.0}) for axis in "XYZ"
    ]
    nuclear = [
        -gn * MU_N * np.kron(np.eye(n_electronic), i_a) for i_a in nuclear_spin_matrices(nuclear_spin)
    ]
    return dict(h0=h0, electronic=electronic, nuclear=nuclear)


def field_operators(model: dict, nuclear_zeeman: bool = True):
    r"""The three operators $\partial H/\partial B_a$ of a `hyperfine_field_model`, with or without the nuclear Zeeman term."""
    if not nuclear_zeeman:
        return list(model["electronic"])
    return [v_el + v_nuc for v_el, v_nuc in zip(model["electronic"], model["nuclear"])]


def linear_field_hamiltonian(h0: np.array, operators: list, field):
    r"""$H(\vec B) = H_0 + \sum_a B_a V_a$."""
    return h0 + sum(b * v for b, v in zip(field, operators))


def polar_field(magnitude: float, theta: float, phi: float):
    r"""Field vector in the frame (D1, D2, b) from $|B|$, the polar angle $\theta$ from b and the azimuth $\phi$ from
    D1 (degrees), the convention of Matsuura et al."""
    theta, phi = np.radians(theta), np.radians(phi)
    return magnitude * np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])


def spin_hamiltonian_model(
    a_tensor: np.array,
    q_tensor: np.array,
    g_tensor: np.array,
    nuclear_spin: float = ER167_NUCLEAR_SPIN,
    gn: float = ER167_GN,
):
    r"""Effective spin Hamiltonian $I\cdot A\cdot S + I\cdot Q\cdot I + \mu_B \vec B\cdot g\cdot S - g_n\mu_N \vec B\cdot I$
    of an effective spin $S=1/2$ in MHz (field in tesla), in the basis $|m_S\rangle\otimes|m_I\rangle$.

    Returns:
        (np.array, list): $H_0$ and the three field operators, as for `linear_field_hamiltonian`.
    """
    spin_half = nuclear_spin_matrices(0.5)
    nuclear = nuclear_spin_matrices(nuclear_spin)
    n_nuclear = len(nuclear[0])
    s_ops = [np.kron(s_a, np.eye(n_nuclear)) for s_a in spin_half]
    i_ops = [np.kron(np.eye(2), i_a) for i_a in nuclear]
    h0 = sum(a_tensor[a, b] * i_ops[a] @ s_ops[b] for a in range(3) for b in range(3))
    h0 = h0 + sum(q_tensor[a, b] * i_ops[a] @ i_ops[b] for a in range(3) for b in range(3))
    mu_b, mu_n = MU_B * CM_TO_MHZ, MU_N * CM_TO_MHZ
    operators = [
        mu_b * sum(g_tensor[a, b] * s_ops[b] for b in range(3)) - gn * mu_n * i_ops[a] for a in range(3)
    ]
    return h0, operators


def transition_derivatives(h0: np.array, operators: list, field, i: int, j: int):
    r"""Frequency $\nu = E_j - E_i$ of the transition between the eigenstates `i` and `j` (sorted from the lowest)
    of $H(\vec B)$, its gradient $S_1 = \nabla_B \nu$ (Hellmann-Feynman) and its Hessian $S_2$ (second-order
    perturbation theory, $\partial_a\partial_b E_k = 2\sum_{l\neq k}\mathrm{Re}(V_{a,kl}V_{b,lk})/(E_k-E_l)$).

    Both derivatives are exact for non-degenerate levels.

    Returns:
        (float, np.array, np.array): $\nu$, $S_1$ and $S_2$ in units of `h0`, `h0`/T and `h0`/T^2.
    """
    energies, vectors = np.linalg.eigh(linear_field_hamiltonian(h0, operators, field))

    def derivatives(k):
        # only row k of each operator in the eigenbasis is needed: V_{a,kl} = <k|V_a|l>
        rows = [(vectors[:, k].conj() @ v) @ vectors for v in operators]
        gradient = np.array([np.real(row[k]) for row in rows])
        denominator = energies[k] - energies
        denominator[k] = np.inf
        hessian = np.array(
            [[2 * np.sum(np.real(ra * rb.conj()) / denominator) for rb in rows] for ra in rows]
        )
        return gradient, hessian

    gradient_i, hessian_i = derivatives(i)
    gradient_j, hessian_j = derivatives(j)
    return energies[j] - energies[i], gradient_j - gradient_i, hessian_j - hessian_i


def zefoz_newton(
    h0: np.array,
    operators: list,
    field,
    i: int,
    j: int,
    max_iter: int = 100,
    tol: float = 1e-10,
    max_step: float = 0.02,
):
    r"""Zero-first-order-Zeeman (ZEFOZ) point $\nabla_B \nu = 0$ of the transition `i`-`j` near `field`, found by
    Newton's method on the exact gradient and Hessian of `transition_derivatives`. Steps are capped at `max_step` (T).

    Returns:
        dict: ``field``, ``frequency``, ``gradient``, ``hessian``, ``iterations`` and ``converged``.
    """
    field = np.array(field, dtype=float)
    converged = False
    for iteration in range(1, max_iter + 1):
        _, gradient, hessian = transition_derivatives(h0, operators, field, i, j)
        step = np.linalg.solve(hessian, gradient)
        norm = np.linalg.norm(step)
        if norm > max_step:
            step *= max_step / norm
        field = field - step
        if norm < tol:
            converged = True
            break
    frequency, gradient, hessian = transition_derivatives(h0, operators, field, i, j)
    return dict(
        field=field,
        frequency=frequency,
        gradient=gradient,
        hessian=hessian,
        iterations=iteration,
        converged=converged,
    )


def zefoz_search(
    h0: np.array,
    operators: list,
    starts: list,
    transitions: list,
    max_field: Optional[float] = None,
    distinct: float = 1e-5,
    **newton_options,
):
    r"""ZEFOZ points of several transitions, from Newton searches (`zefoz_newton`) started at every field of `starts`.

    Converged points with $|\vec B|$ above `max_field` are dropped. Since $H(-\vec B)$ is the time reverse of
    $H(\vec B)$, $\vec B$ and $-\vec B$ give the same point; it is stored once, with $B_b \geq 0$, and points closer
    than `distinct` (T) are merged.

    Returns:
        list: one dict per point, as returned by `zefoz_newton`, with the key ``transition`` added.
    """
    points = []
    for i, j in transitions:
        found = []
        for start in starts:
            result = zefoz_newton(h0, operators, start, i, j, **newton_options)
            if not result["converged"]:
                continue
            field = result["field"] if result["field"][2] >= 0 else -result["field"]
            if max_field is not None and np.linalg.norm(field) > max_field:
                continue
            if any(np.linalg.norm(field - point["field"]) < distinct for point in found):
                continue
            found.append(dict(result, field=field, transition=(i, j)))
        points += found
    return points
