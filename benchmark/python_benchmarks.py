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
from quspin.basis import spin_basis_1d
from quspin.operators import hamiltonian
from quspin.tools.lanczos import lanczos_full
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


def xxz_hamiltonian(length: int, nup: int, static_fmt: str):
    basis = spin_basis_1d(L=length, Nup=nup, pauli=False)
    jxy = math.sqrt(2.0)
    jzz = 1.0
    hz = 1.0 / math.sqrt(3.0)
    static = [
        ["+-", [[jxy / 2.0, site, site + 1] for site in range(length - 1)]],
        ["-+", [[jxy / 2.0, site, site + 1] for site in range(length - 1)]],
        ["zz", [[jzz, site, site + 1] for site in range(length - 1)]],
        ["z", [[hz, site] for site in range(length)]],
    ]
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


def make_cases() -> list[Case]:
    basis_12, hamiltonian_12_csr = xxz_hamiltonian(12, 6, "csr")
    _, hamiltonian_12_dense = xxz_hamiltonian(12, 6, "dense")
    hamiltonian_12_csc = hamiltonian_12_csr.tocsc()
    matrix_12_dense = hamiltonian_12_dense.toarray()
    vector_12 = deterministic_state(basis_12.Ns)
    basis_10, hamiltonian_10_csr = xxz_hamiltonian(10, 5, "csr")
    _, hamiltonian_10_dense = xxz_hamiltonian(10, 5, "dense")
    hamiltonian_10_csc = hamiltonian_10_csr.tocsc()
    matrix_10_dense = hamiltonian_10_dense.toarray()
    vector_10 = deterministic_state(basis_10.Ns)
    eigsh_seed_10 = np.real(vector_10)
    eigsh_seed_10 /= np.linalg.norm(eigsh_seed_10)
    basis_8, hamiltonian_8_csr = xxz_hamiltonian(8, 4, "csr")
    vector_8 = deterministic_state(basis_8.Ns)
    full_basis_12 = spin_basis_1d(L=12, pauli=False)
    full_state_12 = deterministic_state(full_basis_12.Ns)
    times = np.linspace(0.0, 1.0, 9)

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
            "entanglement_entropy",
            "integration",
            "storage_independent",
            "state_vector",
            "L=12;dimension=4096;subsystem=6",
            lambda: full_basis_12.ent_entropy(
                full_state_12, sub_sys_A=range(6)
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
    for _ in range(3):
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
        "note": "",
        "parameters": case.parameters,
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
