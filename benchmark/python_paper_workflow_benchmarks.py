#!/usr/bin/env python3
"""Paper-scale end-to-end ED workflows for the pinned Python QuSpin baseline.

The cases deliberately time complete user operations (basis, Hamiltonian, and
solver/observable) rather than isolated matrix kernels.  Correctness is checked
before timing.  Use ``--case`` to run one heavy case in an isolated CI step.
"""

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
from typing import Callable

import numpy as np
import quspin
from quspin.basis import (
    spin_basis_1d,
    spinful_fermion_basis_1d,
    spinless_fermion_basis_1d,
)
from quspin.operators import hamiltonian
from quspin.tools.lanczos import expm_lanczos, lanczos_full
from scipy.linalg import expm
from scipy.sparse.linalg import eigsh


SINK: object | None = None


@dataclass(frozen=True)
class PaperCase:
    case_id: str
    family_id: str
    name: str
    parameters: str
    operation: Callable[[], dict[str, object]]


def uniform_seed(size: int) -> np.ndarray:
    """Basis-order-independent physical seed shared with the Julia cases."""
    return np.full(size, 1.0 / math.sqrt(size), dtype=np.float64)


def spin_exchange_static(
    length: int,
    *,
    delta: float = 1.0,
    periodic: bool = False,
):
    bonds = [(site, site + 1) for site in range(length - 1)]
    if periodic:
        bonds.append((length - 1, 0))
    return [
        ["+-", [[0.5, left, right] for left, right in bonds]],
        ["-+", [[0.5, left, right] for left, right in bonds]],
        ["zz", [[delta, left, right] for left, right in bonds]],
    ]


def make_hamiltonian(static, basis, *, dtype=np.float64):
    return hamiltonian(
        static,
        [],
        basis=basis,
        dtype=dtype,
        static_fmt="csc",
        check_herm=False,
        check_symm=False,
        check_pcon=False,
    )


def mbl_shift_invert() -> dict[str, object]:
    length = 14
    basis = spin_basis_1d(L=length, Nup=length // 2, pauli=False)
    fields = [
        2.13,
        -1.77,
        0.31,
        3.24,
        -2.63,
        0.82,
        1.46,
        -3.17,
        2.71,
        -0.54,
        1.09,
        -2.28,
        0.67,
        2.94,
    ]
    static = spin_exchange_static(length, periodic=True)
    static.append(["z", [[fields[site], site] for site in range(length)]])
    operator = make_hamiltonian(static, basis)
    seed = uniform_seed(basis.Ns)
    values, vectors = eigsh(
        operator.static,
        k=6,
        sigma=0.0,
        which="LM",
        v0=seed,
        ncv=32,
        maxiter=5_000,
        tol=1e-9,
    )
    residual = np.linalg.norm(
        operator.static @ vectors - vectors * values[np.newaxis, :]
    )
    return {
        "valid": residual < 2e-7 and np.all(np.isfinite(values)),
        "residual": residual,
        "value": values,
    }


def xxz_lanczos_quench() -> dict[str, object]:
    length = 16
    basis = spin_basis_1d(L=length, Nup=length // 2, pauli=False)
    operator = make_hamiltonian(
        spin_exchange_static(length, delta=0.8, periodic=False),
        basis,
    )
    neel = sum(1 << site for site in range(0, length, 2))
    matches = np.flatnonzero(basis.states == neel)
    if matches.size != 1:
        raise RuntimeError("Néel product state is absent from the basis")
    initial = np.zeros(basis.Ns, dtype=np.complex128)
    initial[matches[0]] = 1.0
    energies, vectors, krylov_vectors = lanczos_full(
        operator,
        initial,
        80,
        full_ortho=True,
    )
    evolved = expm_lanczos(
        energies,
        vectors,
        krylov_vectors,
        a=-0.7j,
    )
    return {
        "valid": abs(np.linalg.norm(evolved) - 1.0) < 2e-9,
        "residual": abs(np.linalg.norm(evolved) - 1.0),
        "value": evolved,
    }


def floquet_heating() -> dict[str, object]:
    length = 9
    basis = spin_basis_1d(L=length)
    bonds = [(site, (site + 1) % length) for site in range(length)]
    first = make_hamiltonian(
        [["zz", [[0.9, left, right] for left, right in bonds]]],
        basis,
    )
    second = make_hamiltonian(
        [["x", [[0.73, site] for site in range(length)]]],
        basis,
    )
    dt_first, dt_second = 0.17, 0.23
    unitary = expm(-1j * dt_second * second.toarray()) @ expm(
        -1j * dt_first * first.toarray()
    )
    phases = np.linalg.eigvals(unitary)
    residual = np.linalg.norm(
        unitary.conj().T @ unitary - np.eye(basis.Ns)
    ) / basis.Ns
    return {
        "valid": residual < 3e-11
        and np.max(np.abs(np.abs(phases) - 1.0)) < 3e-10,
        "residual": residual,
        "value": phases,
    }


def spinful_hubbard() -> dict[str, object]:
    length = 8
    basis = spinful_fermion_basis_1d(
        L=length,
        Nf=(length // 2, length // 2),
    )
    bonds = [(site, site + 1) for site in range(length - 1)]
    static = [
        ["+-|", [[-1.0, left, right] for left, right in bonds]],
        ["-+|", [[1.0, left, right] for left, right in bonds]],
        ["|+-", [[-1.0, left, right] for left, right in bonds]],
        ["|-+", [[1.0, left, right] for left, right in bonds]],
        ["n|n", [[4.0, site, site] for site in range(length)]],
    ]
    operator = make_hamiltonian(static, basis)
    seed = uniform_seed(basis.Ns)
    values, vectors = eigsh(
        operator.static,
        k=6,
        which="SA",
        v0=seed,
        ncv=32,
        maxiter=5_000,
        tol=1e-9,
    )
    residual = np.linalg.norm(
        operator.static @ vectors - vectors * values[np.newaxis, :]
    )
    return {
        "valid": residual < 2e-7,
        "residual": residual,
        "value": values,
    }


def interacting_ssh() -> dict[str, object]:
    length = 16
    particles = length // 2
    basis = spinless_fermion_basis_1d(L=length, Nf=particles)
    bonds = [(site, site + 1) for site in range(length - 1)]
    hopping = [
        0.6 if site % 2 == 0 else 1.0 for site in range(length - 1)
    ]
    static = [
        [
            "+-",
            [
                [-hopping[index], left, right]
                for index, (left, right) in enumerate(bonds)
            ],
        ],
        [
            "-+",
            [
                [hopping[index], left, right]
                for index, (left, right) in enumerate(bonds)
            ],
        ],
        ["nn", [[2.0, left, right] for left, right in bonds]],
    ]
    operator = make_hamiltonian(static, basis)
    seed = uniform_seed(basis.Ns)
    values, vectors = eigsh(
        operator.static,
        k=6,
        which="SA",
        v0=seed,
        ncv=32,
        maxiter=5_000,
        tol=1e-9,
    )
    residual = np.linalg.norm(
        operator.static @ vectors - vectors * values[np.newaxis, :]
    )
    return {
        "valid": residual < 2e-7,
        "residual": residual,
        "value": values,
    }


def translation_sector_xxz() -> dict[str, object]:
    length = 18
    basis = spin_basis_1d(
        L=length,
        Nup=length // 2,
        pauli=False,
        kblock=0,
    )
    operator = make_hamiltonian(
        spin_exchange_static(length, delta=0.9, periodic=True),
        basis,
    )
    seed = uniform_seed(basis.Ns)
    values, vectors = eigsh(
        operator.static,
        k=4,
        which="SA",
        v0=seed,
        ncv=24,
        maxiter=5_000,
        tol=1e-9,
    )
    residual = np.linalg.norm(
        operator.static @ vectors - vectors * values[np.newaxis, :]
    )
    return {
        "valid": residual < 2e-7,
        "residual": residual,
        "value": values,
    }


def make_cases() -> list[PaperCase]:
    return [
        PaperCase(
            "paper_mbl_shift_invert_l14",
            "19021",
            "MBL mid-spectrum shift-invert",
            "L=14;Nup=7;dimension=3432;k=6;sigma=0;ncv=32",
            mbl_shift_invert,
        ),
        PaperCase(
            "paper_xxz_lanczos_quench_l16",
            "3831",
            "XXZ Lanczos quench",
            "L=16;Nup=8;dimension=12870;m=80;t=0.7",
            xxz_lanczos_quench,
        ),
        PaperCase(
            "paper_floquet_heating_l9",
            "1069",
            "Floquet heating full unitary",
            "L=9;dimension=512;steps=2",
            floquet_heating,
        ),
        PaperCase(
            "paper_spinful_hubbard_l8",
            "3127",
            "Spinful Hubbard low-energy spectrum",
            "L=8;Nf=(4,4);dimension=4900;k=6",
            spinful_hubbard,
        ),
        PaperCase(
            "paper_interacting_ssh_l16",
            "19000",
            "Interacting SSH low-energy spectrum",
            "L=16;Nf=8;dimension=12870;k=6;t1=0.6;t2=1;V=2",
            interacting_ssh,
        ),
        PaperCase(
            "paper_translation_xxz_l18",
            "3831",
            "Translation-sector XXZ spectrum",
            "L=18;Nup=9;kblock=0;parent_dimension=48620;k=4",
            translation_sector_xxz,
        ),
    ]


def percentile(values: list[float], probability: float) -> float:
    ordered = np.sort(np.asarray(values, dtype=np.float64))
    return float(np.quantile(ordered, probability))


def time_case(case: PaperCase, samples: int) -> dict[str, object]:
    global SINK
    preflight = case.operation()
    if not bool(preflight["valid"]):
        raise RuntimeError(
            f"{case.case_id} preflight failed: residual={preflight['residual']}"
        )
    SINK = case.operation()
    timings: list[float] = []
    for _ in range(samples):
        gc.collect()
        started = time.perf_counter()
        SINK = case.operation()
        timings.append(time.perf_counter() - started)
    return {
        "language": "python",
        "suite": "paper",
        "case_id": case.case_id,
        "family_id": case.family_id,
        "benchmark": case.name,
        "category": "workflow",
        "comparison": "end_to_end",
        "storage": "current_backend",
        "supported": "true",
        "note": "Full basis + Hamiltonian + solver/observable pipeline.",
        "parameters": case.parameters,
        "validation": "passed",
        "samples": samples,
        "iterations_per_sample": 1,
        "median_seconds": statistics.median(timings),
        "mean_seconds": statistics.fmean(timings),
        "stdev_seconds": statistics.stdev(timings) if samples > 1 else 0.0,
        "min_seconds": min(timings),
        "p05_seconds": percentile(timings, 0.05),
        "p25_seconds": percentile(timings, 0.25),
        "p75_seconds": percentile(timings, 0.75),
        "p95_seconds": percentile(timings, 0.95),
        "max_seconds": max(timings),
        "median_allocated_bytes": "",
        "runtime": (
            f"Python {platform.python_version()}; "
            f"QuSpin {quspin.__version__}"
        ),
        "raw_samples_seconds": ";".join(f"{value:.12g}" for value in timings),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--case", action="append", dest="case_ids")
    parser.add_argument("--samples", type=int, default=5)
    args = parser.parse_args()
    if args.samples < 1:
        parser.error("--samples must be positive")
    cases = make_cases()
    if args.case_ids:
        requested = set(args.case_ids)
        known = {case.case_id for case in cases}
        missing = requested - known
        if missing:
            parser.error(f"unknown case(s): {', '.join(sorted(missing))}")
        cases = [case for case in cases if case.case_id in requested]
    rows = [time_case(case, args.samples) for case in cases]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
