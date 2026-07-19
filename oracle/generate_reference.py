#!/usr/bin/env python3
"""Generate a small, deterministic QuSpin reference corpus.

This program is an oracle only.  It is run while preparing the private Minos
verification campaign and is never a dependency of the Julia package.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import quspin
from quspin.basis import (
    basis_int_to_python_int,
    bitwise_and,
    bitwise_leftshift,
    bitwise_not,
    bitwise_or,
    bitwise_rightshift,
    bitwise_xor,
    coherent_state,
    get_basis_type,
    photon_Hspace_dim,
    python_int_to_basis_int,
    spin_basis_1d,
    uint1024,
    uint16384,
    uint256,
    uint4096,
)
from quspin.operators import (
    anti_commutator,
    commutator,
    exp_op,
    hamiltonian,
    quantum_LinearOperator,
    quantum_operator,
)
from quspin.tools.evolution import ED_state_vs_time, ExpmMultiplyParallel, evolve
from quspin.tools.Floquet import Floquet, Floquet_t_vec
from quspin.tools.block_tools import block_diag_hamiltonian
from quspin.tools.lanczos import (
    FTLM_static_iteration,
    LTLM_static_iteration,
    expm_lanczos,
    lanczos_full,
    lin_comb_Q_T,
)
from quspin.tools.misc import (
    KL_div,
    array_to_ints,
    ints_to_array,
    matvec,
    mean_level_spacing,
    project_op,
)
from quspin.tools.measurements import diag_ensemble, ent_entropy, obs_vs_time


def canonical_float(value: float, significant_digits: int = 12) -> float:
    """Remove platform-specific BLAS tails from serialized observations."""

    return float(f"{float(value):.{significant_digits}g}")


def spin_basis_case() -> dict:
    basis = spin_basis_1d(L=4, Nup=2, pauli=False)
    return {
        "id": "spin_basis_1d_L4_Nup2",
        "target_symbol": "SpinBasis1D",
        "parameters": {"L": 4, "nup": 2, "pauli": False},
        "dimension": int(basis.Ns),
        # Ordering is deliberately not contractual.  The verifier compares
        # this as a set because a Julia-native implementation may enumerate
        # the basis in a different, equally valid order.
        "states_as_set": sorted(int(x) for x in basis.states),
    }


def xxz_spectrum_case() -> dict:
    """A four-site open XXZ chain in the two-up-spin sector.

    Site labels in the serialized contract are one-based for the Julia API.
    QuSpin receives the corresponding zero-based labels below.
    """

    L = 4
    basis = spin_basis_1d(L=L, Nup=2, pauli=False)
    jxy = np.sqrt(2.0)
    jzz = 1.0
    hz = 1.0 / np.sqrt(3.0)

    j_zz = [[jzz, i, i + 1] for i in range(L - 1)]
    j_xy = [[jxy / 2.0, i, i + 1] for i in range(L - 1)]
    h_z = [[hz, i] for i in range(L)]
    static = [["+-", j_xy], ["-+", j_xy], ["zz", j_zz], ["z", h_z]]

    H = hamiltonian(
        static,
        [],
        basis=basis,
        dtype=np.float64,
        check_herm=False,
        check_symm=False,
        check_pcon=False,
    )
    dense = np.asarray(H.toarray())
    eigenvalues = np.linalg.eigvalsh(dense)

    return {
        "id": "xxz_open_L4_Nup2",
        "target_symbol": "Hamiltonian",
        "parameters": {
            "L": L,
            "nup": 2,
            "pauli": False,
            "site_indexing": "one_based",
            "jxy": canonical_float(jxy),
            "jzz": jzz,
            "hz": canonical_float(hz),
            "boundary": "open",
        },
        "dimension": int(basis.Ns),
        "spectrum": [canonical_float(x) for x in eigenvalues],
        "trace": canonical_float(np.trace(dense)),
        "frobenius_norm": canonical_float(np.linalg.norm(dense)),
    }


def _dtype_contract(dtype: type) -> str:
    dtype = np.dtype(dtype)
    known = {
        np.dtype(np.uint32): "UInt32",
        np.dtype(np.uint64): "UInt64",
        np.dtype(uint256): "UInt256",
        np.dtype(uint1024): "UInt1024",
        np.dtype(uint4096): "UInt4096",
        np.dtype(uint16384): "UInt16384",
    }
    return known[dtype]


def basis_integer_case() -> dict:
    values = [
        0,
        2**32 - 1,
        2**32,
        2**64 - 1,
        2**64,
        2**255,
        2**256,
    ]
    round_trips = []
    for value in values:
        encoded = python_int_to_basis_int(value)
        round_trips.append(
            {
                "value": str(value),
                "dtype": _dtype_contract(encoded.dtype),
                "decoded": str(basis_int_to_python_int(encoded)),
            }
        )

    type_queries = [
        (10, 5, 2),
        (32, None, 2),
        (64, None, 2),
        (100, 1, 2),
        (300, None, 2),
        (5000, None, 2),
    ]
    return {
        "id": "basis_integer_boundaries",
        "target_symbol": "FixedUInt",
        "round_trips": round_trips,
        "type_queries": [
            {
                "N": N,
                "Np": Np,
                "sps": sps,
                "dtype": _dtype_contract(get_basis_type(N, Np, sps)),
            }
            for N, Np, sps in type_queries
        ],
    }


def bitwise_case() -> dict:
    left = np.array([0x0F, 0xF0, 0xAA], dtype=np.uint32)
    right = np.array([0x33, 0x55, 0xFF], dtype=np.uint32)
    return {
        "id": "basis_bitwise_uint32",
        "target_symbol": "bitwise_and",
        "left": left.tolist(),
        "right": right.tolist(),
        "and": bitwise_and(left, right).tolist(),
        "or": bitwise_or(left, right).tolist(),
        "xor": bitwise_xor(left, right).tolist(),
        "not_uint32": bitwise_not(
            np.array([0, 2**32 - 1], dtype=np.uint32)
        ).tolist(),
        "leftshift": bitwise_leftshift(
            np.array([1, 2, 3], dtype=np.uint32),
            np.array([1, 2, 3], dtype=np.uint32),
        ).tolist(),
        "rightshift": bitwise_rightshift(
            np.array([8, 16, 32], dtype=np.uint32),
            np.array([1, 2, 3], dtype=np.uint32),
        ).tolist(),
    }


def photon_helpers_case() -> dict:
    real_state = coherent_state(0.5, 5, dtype=np.float64)
    complex_state = coherent_state(0.5 + 0.25j, 5, dtype=np.complex128)
    return {
        "id": "photon_helpers",
        "target_symbol": "coherent_state",
        "coherent_real": [float(value) for value in real_state],
        "coherent_complex": [
            [float(value.real), float(value.imag)] for value in complex_state
        ],
        "dimensions": [
            {
                "N": N,
                "Ntot": Ntot,
                "Nph": Nph,
                "dimension": int(photon_Hspace_dim(N, Ntot, Nph)),
            }
            for N, Ntot, Nph in [
                (4, None, 3),
                (4, 0, None),
                (4, 1, None),
                (4, 2, None),
                (4, 4, None),
                (8, 4, None),
            ]
        ],
    }


def misc_tools_case() -> dict:
    encoded = np.array([0, 1, 2, 3, 10], dtype=np.uint32)
    binary = ints_to_array(encoded, N=4)
    A = np.array([[1.0, 2.0], [3.0, 4.0]])
    v = np.array([2.0, -1.0])
    P = np.array([[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]])
    return {
        "id": "misc_tools",
        "target_symbol": "array_to_ints",
        "binary": binary.tolist(),
        "round_trip": [
            str(basis_int_to_python_int(value))
            for value in array_to_ints(binary)
        ],
        "kl_divergences": [
            float(KL_div([0.25, 0.75], [0.5, 0.5])),
            float(KL_div([0.1, 0.2, 0.3, 0.4], [0.4, 0.3, 0.2, 0.1])),
        ],
        "mean_level_spacings": [
            float(mean_level_spacing([0.0, 1.0, 3.0, 6.0])),
            float(mean_level_spacing([0.0, 1.0, 2.0, 4.0, 8.0])),
        ],
        "matvec": matvec(A, v, a=2.0).tolist(),
        "project_down": project_op(
            np.diag([1.0, 2.0, 3.0]), P
        )["Proj_Obs"].tolist(),
        "project_up": project_op(
            np.diag([4.0, 5.0]), P
        )["Proj_Obs"].tolist(),
    }


def operator_krylov_evolution_case() -> dict:
    A = np.array(
        [
            [1.0, 0.2, 0.0, 0.0],
            [0.2, 2.0, 0.3, 0.0],
            [0.0, 0.3, 3.0, 0.4],
            [0.0, 0.0, 0.4, 4.0],
        ]
    )
    B = np.diag([2.0, -1.0, 0.5, 3.0])
    v0 = np.array([1.0, 2.0, -1.0, 0.5])
    E, V, Q_T = lanczos_full(A, v0, 3)
    exact_E, exact_V = np.linalg.eigh(A)
    psi = np.array([1.0, 0.0, 0.0, 0.0], dtype=np.complex128)
    psi_t = ED_state_vs_time(psi, exact_E, exact_V, [0.0, 0.25, 1.0])
    ode_times = np.array([0.0, 0.2, 1.0])
    ode_states = evolve(
        np.array([1.0 + 0.0j]),
        0.0,
        ode_times,
        lambda time, state, frequency: -1j * frequency * state,
        f_params=(2.0,),
    )
    beta = np.array([0.0, 0.5, 2.0])
    ftlm, ftlm_identity = FTLM_static_iteration(
        {"A": A, "B": B}, E, V, Q_T, beta=beta
    )
    ltlm, ltlm_identity = LTLM_static_iteration(
        {"A": A, "B": B}, E, V, Q_T, beta=beta
    )

    return {
        "id": "operator_krylov_evolution",
        "target_symbol": "lanczos_full",
        "commutator": commutator(A, B).tolist(),
        "anti_commutator": anti_commutator(A, B).tolist(),
        "lanczos_E": E.tolist(),
        "lanczos_V": V.tolist(),
        "lanczos_Q_T": Q_T.tolist(),
        "linear_combination": lin_comb_Q_T(
            [1.0, -2.0, 0.5], Q_T
        ).tolist(),
        "expm_lanczos": expm_lanczos(E, V, Q_T, a=-0.25).tolist(),
        "evolved_last": [
            [float(value.real), float(value.imag)]
            for value in psi_t[:, -1]
        ],
        "ode_evolution": [
            [float(value.real), float(value.imag)]
            for value in ode_states[0, :]
        ],
        "ftlm": {
            "identity": ftlm_identity.tolist(),
            "A": ftlm["A"].tolist(),
            "B": ftlm["B"].tolist(),
        },
        "ltlm": {
            "identity": ltlm_identity.tolist(),
            "A": ltlm["A"].tolist(),
            "B": ltlm["B"].tolist(),
        },
    }


def expm_multiply_parallel_case() -> dict:
    generator = np.array([[0.0, -1.0], [1.0, 0.0]])
    operator = ExpmMultiplyParallel(generator, a=0.25)
    real_action = operator.dot(np.array([1.0, 0.0]))
    operator.set_a(-0.5j)
    complex_action = operator.dot(np.array([1.0, 0.0], dtype=np.complex128))
    return {
        "id": "expm_multiply_parallel",
        "target_symbol": "ExpmMultiplyParallel",
        "real_action": [float(value) for value in real_action],
        "complex_action": [
            [float(value.real), float(value.imag)]
            for value in complex_action
        ],
    }


def floquet_time_vector_case() -> dict:
    times = Floquet_t_vec(2 * np.pi, 2, len_T=4, N_up=1, N_down=1)
    return {
        "id": "floquet_time_vector",
        "target_symbol": "FloquetTimeVector",
        "N": int(times.N),
        "T": float(times.T),
        "len_T": int(times.len_T),
        "vals": times.vals.tolist(),
        "strobo_inds_zero_based": times.strobo.inds.tolist(),
        "strobo_vals": times.strobo.vals.tolist(),
        "up_vals": times.up.vals.tolist(),
        "constant_vals": times.const.vals.tolist(),
        "down_vals": times.down.vals.tolist(),
        "coordinates_zero_based": [
            int(value) for value in times.get_coordinates(6)
        ],
    }


def floquet_case() -> dict:
    basis = spin_basis_1d(1, pauli=False)
    H = hamiltonian(
        [["z", [[2.0, 0]]]],
        [],
        basis=basis,
        dtype=np.float64,
        check_herm=False,
        check_symm=False,
        check_pcon=False,
    )
    floquet = Floquet(
        {"H": H, "T": 0.5},
        HF=True,
        UF=True,
        thetaF=True,
        VF=True,
    )
    return {
        "id": "floquet_static",
        "target_symbol": "Floquet",
        "T": float(floquet.T),
        "EF": floquet.EF.tolist(),
        "UF": [
            [[float(value.real), float(value.imag)] for value in row]
            for row in floquet.UF
        ],
        "thetaF": [
            [float(value.real), float(value.imag)]
            for value in floquet.thetaF
        ],
    }


def measurements_case() -> dict:
    basis = spin_basis_1d(4)
    ghz = np.zeros(basis.Ns, dtype=np.complex128)
    ghz[basis.index(0)] = 1 / np.sqrt(2)
    ghz[basis.index(15)] = 1 / np.sqrt(2)
    entropy = ent_entropy(ghz, basis=basis, DM="both")

    A = np.array(
        [
            [1.0, 0.2, 0.0, 0.0],
            [0.2, 2.0, 0.3, 0.0],
            [0.0, 0.3, 3.0, 0.4],
            [0.0, 0.0, 0.4, 4.0],
        ]
    )
    B = np.diag([2.0, -1.0, 0.5, 3.0])
    E, V = np.linalg.eigh(A)
    observations = obs_vs_time(
        (np.array([1.0, 0.0, 0.0, 0.0], dtype=np.complex128), E, V),
        np.array([0.0, 0.25, 1.0]),
        {"A": A, "B": B},
        return_state=True,
        enforce_pure=True,
    )

    return {
        "id": "measurements",
        "target_symbol": "ent_entropy",
        "entropy": float(entropy["Sent"]),
        "rdm_A": entropy["DM_chain_subsys"].real.tolist(),
        "rdm_B": entropy["DM_other_subsys"].real.tolist(),
        "expectations": {
            key: [
                [float(value.real), float(value.imag)]
                for value in observations[key]
            ]
            for key in ("A", "B")
        },
    }


def exp_op_case() -> dict:
    O = np.array([[0.0, 1.0], [-1.0, 0.0]])
    vector = np.array([1.0, 0.0])
    observable = np.diag([2.0, 3.0])
    operator = exp_op(O, a=0.5)
    grid_operator = exp_op(
        O,
        a=0.25,
        start=0.0,
        stop=1.0,
        num=3,
        endpoint=True,
    )
    return {
        "id": "exp_op",
        "target_symbol": "ExpOp",
        "matrix": operator.get_mat(dense=True).tolist(),
        "dot": operator.dot(vector).tolist(),
        "rdot": operator.rdot(vector).tolist(),
        "sandwich": operator.sandwich(observable).tolist(),
        "grid": grid_operator.grid.tolist(),
        "step": float(grid_operator.step),
        "grid_dot": grid_operator.dot(vector).tolist(),
        "grid_sandwich": grid_operator.sandwich(observable).tolist(),
    }


def diagonal_ensemble_case() -> dict:
    eigenvectors = np.array([[1.0, 1.0], [1.0, -1.0]]) / np.sqrt(2.0)
    energies = np.array([-1.0, 1.0])
    state = np.array([1.0, 0.0])
    observable = np.diag([1.0, -1.0])
    result = diag_ensemble(
        2,
        state,
        energies,
        eigenvectors,
        density=False,
        rho_d=True,
        Obs=observable,
        Sd_Renyi=True,
    )
    return {
        "id": "diagonal_ensemble",
        "target_symbol": "diag_ensemble",
        "rho_d": result["rho_d"].tolist(),
        "Obs_pure": float(result["Obs_pure"]),
        "Sd_pure": float(result["Sd_pure"]),
    }


def block_tools_case() -> dict:
    blocks = [{"Nup": sector} for sector in range(3)]
    static = [["z", [[1.0, 0], [0.5, 1]]]]
    projector, block_hamiltonian = block_diag_hamiltonian(
        blocks,
        static,
        [],
        spin_basis_1d,
        (2,),
        np.complex128,
        basis_kwargs={"pauli": False},
        check_symm=False,
        check_herm=False,
        check_pcon=False,
    )
    matrix = block_hamiltonian.toarray()
    return {
        "id": "block_tools_spin_sectors",
        "target_symbol": "block_diag_hamiltonian",
        "projector_shape": list(projector.shape),
        "block_matrix": [
            [[float(value.real), float(value.imag)] for value in row]
            for row in matrix
        ],
    }


def quantum_linear_operator_case() -> dict:
    basis = spin_basis_1d(3, Nup=1, pauli=False)
    static = [
        ["+-", [[0.65, 0, 1], [-0.2, 1, 2]]],
        ["-+", [[0.65, 0, 1], [-0.2, 1, 2]]],
        ["z", [[0.17, 0], [0.31, 1], [-0.23, 2]]],
    ]
    operator = quantum_LinearOperator(
        static,
        basis=basis,
        dtype=np.complex128,
        check_herm=False,
        check_symm=False,
        check_pcon=False,
    )
    identity = np.eye(operator.Ns, dtype=np.complex128)
    matrix = np.column_stack(
        [operator.dot(identity[:, index]) for index in range(operator.Ns)]
    )
    return {
        "id": "quantum_linear_operator",
        "target_symbol": "QuantumLinearOperator",
        "matrix": matrix.real.tolist(),
        "spectrum": np.linalg.eigvalsh(matrix).tolist(),
    }


def quantum_operator_case() -> dict:
    basis = spin_basis_1d(2, pauli=False)
    input_dict = {
        "x": [
            ["+-", [[1.0, 0, 1]]],
            ["-+", [[1.0, 0, 1]]],
        ],
        "z": [["z", [[1.0, 0]]]],
    }
    operator = quantum_operator(
        input_dict,
        basis=basis,
        dtype=np.float64,
        check_herm=False,
        check_symm=False,
        check_pcon=False,
    )
    pars = {"x": 0.7, "z": 0.3}
    return {
        "id": "quantum_operator",
        "target_symbol": "QuantumOperator",
        "matrix": operator.toarray(pars).tolist(),
        "spectrum": operator.eigvalsh(pars).tolist(),
        "trace": float(operator.trace(pars)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("reference.json"),
    )
    args = parser.parse_args()

    payload = {
        "source": {
            "package": "quspin",
            "reported_version": str(quspin.__version__),
        },
        "target": {
            "package": "QuSpin.jl",
            "language": "Julia",
            "contract": "Julia-native API; Python is an offline oracle only",
        },
        "cases": [
            basis_integer_case(),
            bitwise_case(),
            block_tools_case(),
            diagonal_ensemble_case(),
            exp_op_case(),
            expm_multiply_parallel_case(),
            floquet_case(),
            floquet_time_vector_case(),
            misc_tools_case(),
            measurements_case(),
            operator_krylov_evolution_case(),
            photon_helpers_case(),
            quantum_linear_operator_case(),
            quantum_operator_case(),
            spin_basis_case(),
            xxz_spectrum_case(),
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(args.output)


if __name__ == "__main__":
    main()
