#!/usr/bin/env python3
"""Warm-operation timing baseline for the pinned Python QuSpin implementation."""

from __future__ import annotations

import argparse
import csv
import gc
import math
import platform
import statistics
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional

import numpy as np
import quspin
from quspin.basis import (
    spin_basis_1d,
    spin_basis_general,
    spinful_fermion_basis_1d,
)
from quspin.operators import exp_op, hamiltonian, quantum_LinearOperator
from quspin.tools.lanczos import lanczos_full
from quspin.tools.measurements import diag_ensemble
from scipy import sparse
from scipy.sparse.linalg import eigsh


SAMPLES = 15
TARGET_SAMPLE_SECONDS = 0.08
MAX_ITERATIONS = 1_000
SINK: Optional[object] = None


@dataclass(frozen=True)
class Case:
    name: str
    category: str
    comparison: str
    storage: str
    parameters: str
    function: Callable[[], object]
    validator: Optional[Callable[[object], bool]] = None
    note: str = ""


def xxz_static(length: int):
    jxy = math.sqrt(2.0)
    jzz = 1.0
    hz = 1.0 / math.sqrt(3.0)
    return [
        ["+-", [[jxy / 2.0, site, site + 1] for site in range(length - 1)]],
        ["-+", [[jxy / 2.0, site, site + 1] for site in range(length - 1)]],
        ["zz", [[jzz, site, site + 1] for site in range(length - 1)]],
        ["z", [[hz, site] for site in range(length)]],
    ]


def xxz_hamiltonian(length: int, nup: int, static_fmt: str):
    basis = spin_basis_1d(L=length, Nup=nup, pauli=False)
    static = xxz_static(length)
    operator = hamiltonian(
        static,
        [],
        basis=basis,
        dtype=np.float64,
        check_herm=False,
        check_symm=False,
        check_pcon=False,
        static_fmt=static_fmt,
    )
    return basis, operator


def deterministic_state(size: int) -> np.ndarray:
    indices = np.arange(1, size + 1, dtype=np.float64)
    state = np.sin(0.17 * indices) + 1j * np.cos(0.11 * indices)
    return state / np.linalg.norm(state)


def deterministic_state_batch(size: int, count: int) -> np.ndarray:
    indices = np.arange(1, size + 1, dtype=np.float64)[:, None]
    columns = np.arange(1, count + 1, dtype=np.float64)[None, :]
    states = (
        np.sin(0.17 * indices + columns)
        + 1j * np.cos(0.11 * indices - columns)
    )
    return states / np.linalg.norm(states, axis=0, keepdims=True)


def general_spin_maps(lx: int, ly: int) -> tuple[np.ndarray, np.ndarray]:
    site = lambda x, y: x + lx * y
    translation_x = np.asarray(
        [site((x + 1) % lx, y) for y in range(ly) for x in range(lx)]
    )
    translation_y = np.asarray(
        [site(x, (y + 1) % ly) for y in range(ly) for x in range(lx)]
    )
    return translation_x, translation_y


def general_spin_basis():
    translation_x, translation_y = general_spin_maps(4, 3)
    return spin_basis_general(
        N=12,
        Nup=6,
        pauli=False,
        kxblock=(translation_x, 1),
        kyblock=(translation_y, 1),
        block_order=["kxblock", "kyblock"],
    )


def general_spin_static():
    lx, ly = 4, 3
    site = lambda x, y: x + lx * y
    bonds = [
        [site(x, y), site((x + 1) % lx, y)]
        for y in range(ly)
        for x in range(lx)
    ]
    bonds += [
        [site(x, y), site(x, (y + 1) % ly)]
        for y in range(ly)
        for x in range(lx)
    ]
    return [
        ["+-", [[0.5, left, right] for left, right in bonds]],
        ["-+", [[0.5, left, right] for left, right in bonds]],
        ["zz", [[1.0, left, right] for left, right in bonds]],
    ]


def general_spin_hamiltonian(basis, static):
    return hamiltonian(
        static,
        [],
        basis=basis,
        dtype=np.complex128,
        static_fmt="csc",
        check_herm=False,
        check_symm=False,
        check_pcon=False,
    )


def spin_one_hamiltonian(length: int):
    basis = spin_basis_1d(
        L=length,
        S="1",
        Nup=length,
        pauli=False,
    )
    bonds = [[site, (site + 1) % length] for site in range(length)]
    static = [
        ["+-", [[0.5, left, right] for left, right in bonds]],
        ["-+", [[0.5, left, right] for left, right in bonds]],
        ["zz", [[1.0, left, right] for left, right in bonds]],
    ]
    return hamiltonian(
        static,
        [],
        basis=basis,
        dtype=np.float64,
        static_fmt="csc",
        check_herm=False,
        check_symm=False,
        check_pcon=False,
    )


def validate_hamiltonian_fingerprint(
    operator, expected_dimension: int, expected_trace: float, expected_norm: float
) -> bool:
    matrix = operator.toarray()
    return (
        matrix.shape == (expected_dimension, expected_dimension)
        and np.allclose(matrix, matrix.T.conj(), rtol=0.0, atol=2e-12)
        and np.isclose(np.trace(matrix).real, expected_trace, rtol=1e-12, atol=1e-12)
        and np.isclose(np.linalg.norm(matrix), expected_norm, rtol=1e-12, atol=1e-12)
    )


def validate_batch_entropy(result: object) -> bool:
    entropy = np.asarray(result["Sent_A"])
    return (
        entropy.shape == (8,)
        and np.all(np.isfinite(entropy))
        and np.isclose(entropy.sum(), 7.360709128824066, rtol=1e-12, atol=1e-12)
    )


def validate_batch_evolution(result: object) -> bool:
    states = np.asarray(result)
    return states.shape == (70, 4, 9) and np.allclose(
        np.sum(np.abs(states) ** 2, axis=0),
        1.0,
        rtol=0.0,
        atol=2e-10,
    )


def spinful_static(length: int):
    bonds = [(site, site + 1) for site in range(length - 1)]
    return [
        ["+-|", [[-1.0, left, right] for left, right in bonds]],
        ["-+|", [[1.0, left, right] for left, right in bonds]],
        ["|+-", [[-1.0, left, right] for left, right in bonds]],
        ["|-+", [[1.0, left, right] for left, right in bonds]],
        ["n|n", [[2.0, site, site] for site in range(length)]],
    ]


def spinful_hamiltonian(length: int):
    particles = length // 2
    basis = spinful_fermion_basis_1d(
        L=length,
        Nf=(particles, particles),
    )
    operator = hamiltonian(
        spinful_static(length),
        [],
        basis=basis,
        dtype=np.float64,
        static_fmt="csc",
        check_herm=False,
        check_symm=False,
        check_pcon=False,
    )
    return basis, operator


def make_cases() -> list[Case]:
    basis_12, hamiltonian_12_csr = xxz_hamiltonian(12, 6, "csr")
    _, hamiltonian_12_dense = xxz_hamiltonian(12, 6, "dense")
    _, hamiltonian_12_dia = xxz_hamiltonian(12, 6, "dia")
    hamiltonian_12_csc = hamiltonian_12_csr.tocsc()
    matrix_12_dense = hamiltonian_12_dense.toarray()
    vector_12 = deterministic_state(basis_12.Ns)
    linear_static_12 = xxz_static(12)
    linear_operator_12 = quantum_LinearOperator(
        linear_static_12,
        basis=basis_12,
        dtype=np.complex128,
        check_herm=False,
        check_symm=False,
        check_pcon=False,
    )
    basis_10, hamiltonian_10_csr = xxz_hamiltonian(10, 5, "csr")
    _, hamiltonian_10_dense = xxz_hamiltonian(10, 5, "dense")
    hamiltonian_10_csc = hamiltonian_10_csr.tocsc()
    matrix_10_dense = hamiltonian_10_dense.toarray()
    vector_10 = deterministic_state(basis_10.Ns)
    eigsh_seed_10 = np.real(vector_10)
    eigsh_seed_10 /= np.linalg.norm(eigsh_seed_10)
    basis_8, hamiltonian_8_csr = xxz_hamiltonian(8, 4, "csr")
    vector_8 = deterministic_state(basis_8.Ns)
    batch_vector_8 = deterministic_state_batch(basis_8.Ns, 4)
    general_basis_12 = general_spin_basis()
    general_static_12 = general_spin_static()
    full_basis_12 = spin_basis_1d(L=12, pauli=False)
    full_state_12 = deterministic_state(full_basis_12.Ns)
    full_basis_10 = spin_basis_1d(L=10, pauli=False)
    batch_state_10 = deterministic_state_batch(full_basis_10.Ns, 8)
    times = np.linspace(0.0, 1.0, 9)
    spinful_basis_6, spinful_hamiltonian_6 = spinful_hamiltonian(6)
    spinful_state_6 = deterministic_state(spinful_basis_6.Ns)
    spinful_qlo_6 = quantum_LinearOperator(
        spinful_static(6),
        basis=spinful_basis_6,
        dtype=np.complex128,
        check_herm=False,
        check_symm=False,
        check_pcon=False,
    )
    exp_dimension = 256
    exp_matrix = sparse.diags(
        (
            np.full(exp_dimension - 1, 0.15),
            np.linspace(-1.0, 1.0, exp_dimension),
            np.full(exp_dimension - 1, 0.15),
        ),
        offsets=(-1, 0, 1),
        format="csc",
    )
    exp_state = deterministic_state(exp_dimension)
    exp_operator = exp_op(exp_matrix, a=-0.2j)
    exp_grid_dimension = 128
    exp_grid_matrix = sparse.diags(
        (
            np.full(exp_grid_dimension - 1, 0.15),
            np.linspace(-1.0, 1.0, exp_grid_dimension),
            np.full(exp_grid_dimension - 1, 0.15),
        ),
        offsets=(-1, 0, 1),
        format="csc",
    )
    exp_grid_state = deterministic_state(exp_grid_dimension)
    exp_grid_operator = exp_op(
        exp_grid_matrix,
        a=-0.2j,
        start=0.0,
        stop=1.0,
        num=9,
    )
    ensemble_dimension = 256
    ensemble_indices = np.arange(ensemble_dimension, dtype=np.float64)
    ensemble_observable = (
        np.sin(0.013 * ensemble_indices[:, None] + 0.017 * ensemble_indices[None, :])
        + 1j
        * np.cos(0.019 * ensemble_indices[:, None] - 0.023 * ensemble_indices[None, :])
    )
    ensemble_observable = ensemble_observable + ensemble_observable.T.conj()
    # diag_ensemble's legacy fluctuation guard evaluates ``bool(Obs)``.
    # NumPy arrays intentionally reject that operation, whereas the documented
    # QuSpin Hamiltonian observable protocol supplies an unambiguous truth
    # value and preserves the dense matrix action used by this benchmark.
    ensemble_observable_operator = hamiltonian(
        [ensemble_observable],
        [],
        dtype=np.complex128,
        static_fmt="dense",
    )
    ensemble_state = deterministic_state(ensemble_dimension)
    ensemble_energies = np.arange(1, ensemble_dimension + 1, dtype=np.float64)
    ensemble_vectors = np.eye(ensemble_dimension, dtype=np.complex128)
    symmetry_basis_14 = spin_basis_1d(
        L=14,
        Nup=7,
        pauli=False,
        kblock=0,
    )
    symmetry_static_14 = [
        [
            "+-",
            [[0.5, site, (site + 1) % 14] for site in range(14)],
        ],
        [
            "-+",
            [[0.5, site, (site + 1) % 14] for site in range(14)],
        ],
    ]

    return [
        Case(
            "spin_basis_construction",
            "api",
            "storage_independent",
            "n/a",
            "L=16;nup=8;dimension=12870",
            lambda: spin_basis_1d(L=16, Nup=8, pauli=False),
        ),
        Case(
            "symmetry_basis_construction",
            "integration",
            "storage_independent",
            "projector",
            "L=14;nup=7;kblock=3;parent_dimension=3432",
            lambda: spin_basis_1d(
                L=14,
                Nup=7,
                pauli=False,
                kblock=3,
            ),
        ),
        Case(
            "symmetry_hamiltonian_construction_sparse",
            "integration",
            "controlled",
            "csc",
            f"L=14;nup=7;kblock=0;dimension={symmetry_basis_14.Ns}",
            lambda: hamiltonian(
                symmetry_static_14,
                [],
                basis=symmetry_basis_14,
                dtype=np.float64,
                static_fmt="csc",
                check_herm=False,
                check_symm=False,
                check_pcon=False,
            ),
        ),
        Case(
            "general_2d_basis_construction",
            "integration",
            "storage_independent",
            "projector",
            "Lx=4;Ly=3;Nup=6;kx=1;ky=1;dimension=75",
            general_spin_basis,
            validator=lambda basis: basis.Ns == 75,
            note="Two commuting translation maps through spin_basis_general.",
        ),
        Case(
            "general_2d_hamiltonian_construction_sparse",
            "integration",
            "controlled",
            "csc",
            "Lx=4;Ly=3;Nup=6;kx=1;ky=1;dimension=75",
            lambda: general_spin_hamiltonian(
                general_basis_12, general_static_12
            ),
            validator=lambda operator: validate_hamiltonian_fingerprint(
                operator, 75, -41.0, 18.30300521772313
            ),
            note="Periodic 2D XXZ assembly in a two-map symmetry sector.",
        ),
        Case(
            "higher_spin_hamiltonian_construction_sparse",
            "integration",
            "controlled",
            "csc",
            "L=8;S=1;Nup=8;dimension=1107;periodic=true",
            lambda: spin_one_hamiltonian(8),
            validator=lambda operator: validate_hamiltonian_fingerprint(
                operator, 1107, -816.0, 112.7829774389737
            ),
            note="Spin-one XXZ construction with exact angular-momentum factors.",
        ),
        Case(
            "xxz_hamiltonian_construction_dense",
            "integration",
            "controlled",
            "dense",
            "L=10;nup=5;dimension=252;open=true",
            lambda: xxz_hamiltonian(10, 5, "dense")[1],
        ),
        Case(
            "xxz_hamiltonian_construction_sparse",
            "integration",
            "controlled",
            "csc",
            "L=10;nup=5;dimension=252;open=true",
            lambda: xxz_hamiltonian(10, 5, "csc")[1],
        ),
        Case(
            "hamiltonian_matvec_current_storage",
            "api",
            "current_backend",
            "csr",
            "L=12;nup=6;dimension=924",
            lambda: hamiltonian_12_csr.dot(vector_12),
        ),
        Case(
            "matrix_matvec_dense",
            "kernel",
            "controlled",
            "dense",
            "L=12;nup=6;dimension=924",
            lambda: matrix_12_dense @ vector_12,
        ),
        Case(
            "matrix_matvec_sparse_csc",
            "kernel",
            "controlled",
            "csc",
            "L=12;nup=6;dimension=924",
            lambda: hamiltonian_12_csc @ vector_12,
        ),
        Case(
            "matrix_matvec_sparse_csr",
            "kernel",
            "controlled",
            "csr",
            "L=12;nup=6;dimension=924",
            lambda: hamiltonian_12_csr.static @ vector_12,
        ),
        Case(
            "matrix_matvec_sparse_dia",
            "kernel",
            "controlled",
            "dia",
            "L=12;nup=6;dimension=924",
            lambda: hamiltonian_12_dia.static @ vector_12,
        ),
        Case(
            "quantum_linear_operator_construction",
            "integration",
            "controlled",
            "matrix_free",
            "L=12;nup=6;dimension=924",
            lambda: quantum_LinearOperator(
                linear_static_12,
                basis=basis_12,
                dtype=np.complex128,
                check_herm=False,
                check_symm=False,
                check_pcon=False,
            ),
        ),
        Case(
            "quantum_linear_operator_matvec",
            "kernel",
            "controlled",
            "matrix_free",
            "L=12;nup=6;dimension=924",
            lambda: linear_operator_12.dot(vector_12),
        ),
        Case(
            "full_eigenspectrum_dense",
            "kernel",
            "controlled",
            "dense",
            "L=10;nup=5;dimension=252",
            lambda: np.linalg.eigvalsh(matrix_10_dense),
        ),
        Case(
            "lanczos_decomposition_dense",
            "integration",
            "controlled",
            "dense",
            "L=10;nup=5;dimension=252;m=32",
            lambda: lanczos_full(matrix_10_dense, vector_10, 32),
        ),
        Case(
            "lanczos_decomposition_sparse_csc",
            "integration",
            "controlled",
            "csc",
            "L=10;nup=5;dimension=252;m=32",
            lambda: lanczos_full(hamiltonian_10_csc, vector_10, 32),
        ),
        Case(
            "partial_eigenspectrum_sparse_csc",
            "integration",
            "controlled",
            "csc",
            "L=10;nup=5;dimension=252;k=4;which=SA",
            lambda: eigsh(
                hamiltonian_10_csc,
                k=4,
                which="SA",
                v0=eigsh_seed_10,
                tol=1e-10,
            ),
        ),
        Case(
            "static_time_evolution_current_storage",
            "integration",
            "current_backend",
            "csr",
            "L=8;nup=4;dimension=70;times=9;tmax=1",
            lambda: hamiltonian_8_csr.evolve(vector_8, 0.0, times),
        ),
        Case(
            "static_batch_time_evolution_current_storage",
            "integration",
            "current_backend",
            "csr_krylov",
            "L=8;nup=4;dimension=70;states=4;times=9;tmax=1",
            lambda: hamiltonian_8_csr.evolve(batch_vector_8, 0.0, times),
            validator=validate_batch_evolution,
            note="Batched column-state evolution across one shared time grid.",
        ),
        Case(
            "entanglement_entropy",
            "integration",
            "storage_independent",
            "state_vector",
            "L=12;dimension=4096;subsystem=6",
            lambda: full_basis_12.ent_entropy(
                full_state_12, sub_sys_A=range(6)
            ),
        ),
        Case(
            "batched_entanglement_entropy",
            "integration",
            "storage_independent",
            "state_matrix",
            "L=10;dimension=1024;states=8;subsystem=5",
            lambda: full_basis_10.ent_entropy(
                batch_state_10,
                sub_sys_A=range(5),
                density=False,
                enforce_pure=True,
            ),
            validator=validate_batch_entropy,
            note="Eight pure states evaluated through the batched entropy API.",
        ),
        Case(
            "spinful_hamiltonian_construction_sparse",
            "integration",
            "controlled",
            "csc",
            "L=6;Nf=(3,3);dimension=400",
            lambda: spinful_hamiltonian(6)[1],
        ),
        Case(
            "spinful_quantum_linear_operator_matvec",
            "kernel",
            "controlled",
            "matrix_free",
            "L=6;Nf=(3,3);dimension=400",
            lambda: spinful_qlo_6.dot(spinful_state_6),
        ),
        Case(
            "expop_sparse_vector_action",
            "kernel",
            "controlled",
            "csc_expm_action",
            "dimension=256;a=-0.2im",
            lambda: exp_operator.dot(exp_state),
        ),
        Case(
            "expop_sparse_grid_action",
            "integration",
            "controlled",
            "csc_expm_action",
            "dimension=128;times=9;a=-0.2im",
            lambda: exp_grid_operator.dot(exp_grid_state),
        ),
        Case(
            "diag_ensemble_quantum_fluctuation",
            "integration",
            "controlled",
            "dense",
            "dimension=256;delta_t=true;delta_q=true",
            lambda: diag_ensemble(
                1,
                ensemble_state,
                ensemble_energies,
                ensemble_vectors,
                density=False,
                Obs=ensemble_observable_operator,
                delta_t_Obs=True,
                delta_q_Obs=True,
            ),
        ),
    ]


def percentile(sorted_values: list[float], probability: float) -> float:
    position = probability * (len(sorted_values) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return sorted_values[lower]
    weight = position - lower
    return (
        sorted_values[lower] * (1.0 - weight)
        + sorted_values[upper] * weight
    )


def time_case(case: Case) -> dict[str, object]:
    global SINK
    SINK = case.function()
    if case.validator is None:
        validation = "smoke"
    else:
        if not case.validator(SINK):
            raise AssertionError(f"validation failed for benchmark {case.name}")
        validation = "passed"
    for _ in range(2):
        SINK = case.function()

    started = time.perf_counter_ns()
    SINK = case.function()
    estimate = max((time.perf_counter_ns() - started) / 1.0e9, 1.0e-9)
    iterations = max(
        1,
        min(MAX_ITERATIONS, math.ceil(TARGET_SAMPLE_SECONDS / estimate)),
    )

    timings: list[float] = []
    for _ in range(SAMPLES):
        gc.collect()
        started = time.perf_counter_ns()
        for _ in range(iterations):
            SINK = case.function()
        elapsed = (time.perf_counter_ns() - started) / 1.0e9
        timings.append(elapsed / iterations)

    ordered = sorted(timings)
    mean = statistics.fmean(timings)
    stdev = statistics.stdev(timings) if len(timings) > 1 else 0.0
    return {
        "language": "python",
        "benchmark": case.name,
        "category": case.category,
        "comparison": case.comparison,
        "storage": case.storage,
        "supported": "true",
        "note": case.note,
        "parameters": case.parameters,
        "validation": validation,
        "samples": len(timings),
        "iterations_per_sample": iterations,
        "median_seconds": statistics.median(timings),
        "mean_seconds": mean,
        "stdev_seconds": stdev,
        "min_seconds": ordered[0],
        "p05_seconds": percentile(ordered, 0.05),
        "p25_seconds": percentile(ordered, 0.25),
        "p75_seconds": percentile(ordered, 0.75),
        "p95_seconds": percentile(ordered, 0.95),
        "max_seconds": ordered[-1],
        "runtime": (
            f"Python {platform.python_version()}; "
            f"QuSpin {getattr(quspin, '__version__', 'unknown')}"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    rows = [time_case(case) for case in make_cases()]
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
