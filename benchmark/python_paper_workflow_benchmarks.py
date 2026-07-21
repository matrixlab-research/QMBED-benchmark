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
    boson_basis_1d,
    spin_basis_1d,
    spinful_fermion_basis_1d,
    spinless_fermion_basis_1d,
    spinless_fermion_basis_general,
)
from quspin.basis.user import op_sig_32, pre_check_state_sig_32, user_basis
from quspin.operators import hamiltonian
from quspin.tools.lanczos import expm_lanczos, lanczos_full
from numba import carray, cfunc, uint32
from scipy.linalg import expm
from scipy.sparse.linalg import eigsh


SINK: object | None = None


@cfunc(op_sig_32, locals={"bit": uint32})
def pxp_operator(op_struct_ptr, op_str, site, length, args):
    """Spin flip callback for the pinned Python QuSpin ``user_basis``."""
    op_struct = carray(op_struct_ptr, 1)[0]
    site = length - site - 1
    bit = 1 << site
    if op_str == 120:  # ord("x")
        op_struct.state ^= bit
        return 0
    op_struct.matrix_ele = 0
    return -1


@cfunc(
    pre_check_state_sig_32,
    locals={"mask": uint32, "shifted_left": uint32, "shifted_right": uint32},
)
def periodic_blockade(state, length, args):
    """Keep bit strings without adjacent excitations, including the wrap bond."""
    mask = 0xFFFFFFFF >> (32 - length)
    shifted_left = ((state << 1) & mask) | (state >> (length - 1))
    shifted_right = (state >> 1) | ((state << (length - 1)) & mask)
    return ((shifted_left | shifted_right) & state) == 0


PXP_OPERATOR_ARGS = np.array([], dtype=np.uint32)


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


def deterministic_seed(size: int) -> np.ndarray:
    values = np.sin(np.arange(1, size + 1, dtype=np.float64))
    return values / np.linalg.norm(values)


def krylov_spectrum(
    operator,
    source: np.ndarray,
    frequencies: np.ndarray,
    *,
    reference_energy: float,
    broadening: float,
    steps: int,
) -> np.ndarray:
    source_norm = np.linalg.norm(source)
    if source_norm == 0:
        return np.zeros_like(frequencies)
    energies, vectors, _ = lanczos_full(
        operator,
        source,
        min(steps, source.size - 1),
        full_ortho=True,
    )
    weights = source_norm**2 * np.abs(vectors[0, :]) ** 2
    offsets = (
        frequencies[:, np.newaxis]
        + reference_energy
        - energies[np.newaxis, :]
    )
    return np.sum(
        weights[np.newaxis, :]
        * broadening
        / (math.pi * (offsets**2 + broadening**2)),
        axis=1,
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


def tfim_fidelity_scan() -> dict[str, object]:
    length = 16
    basis = spin_basis_1d(L=length, pauli=True)
    bonds = [(site, (site + 1) % length) for site in range(length)]
    fields = (0.8, 0.9, 1.0, 1.1, 1.2)
    spaces: list[np.ndarray] = []
    residuals: list[float] = []
    for field in fields:
        operator = make_hamiltonian(
            [
                ["zz", [[-1.0, left, right] for left, right in bonds]],
                ["x", [[-field, site] for site in range(length)]],
            ],
            basis,
        )
        values, vectors = eigsh(
            operator.static,
            k=2,
            which="SA",
            v0=deterministic_seed(basis.Ns),
            ncv=20,
            maxiter=5_000,
            tol=1e-9,
        )
        order = np.argsort(values)
        spaces.append(vectors[:, order])
        residuals.append(
            float(
                np.linalg.norm(
                    operator.static @ vectors
                    - vectors * values[np.newaxis, :]
                )
            )
        )
    fidelities = np.asarray(
        [
            np.min(np.linalg.svd(left.conj().T @ right, compute_uv=False))
            for left, right in zip(spaces[:-1], spaces[1:])
        ]
    )
    residual = max(residuals)
    return {
        "valid": residual < 3e-7
        and np.all(np.isfinite(fidelities))
        and np.all((fidelities > 0) & (fidelities <= 1.0 + 1e-12))
        and int(np.argmin(fidelities)) in (1, 2),
        "residual": residual,
        "value": fidelities,
    }


def pxp_revival() -> dict[str, object]:
    length = 24
    basis = user_basis(
        np.uint32,
        length,
        {"op": pxp_operator, "op_args": PXP_OPERATOR_ARGS},
        allowed_ops={"x"},
        sps=2,
        pre_check_state=(periodic_blockade, None),
        # Python QuSpin requires a strict upper bound, not the exact final Ns.
        Ns_block_est=120_000,
    )
    operator = make_hamiltonian(
        [["x", [[1.0, site] for site in range(length)]]],
        basis,
    )
    neel = sum(1 << site for site in range(0, length, 2))
    matches = np.flatnonzero(basis.states == neel)
    if matches.size != 1:
        raise RuntimeError("staggered PXP state is absent from the constrained basis")
    initial = np.zeros(basis.Ns, dtype=np.complex128)
    initial[matches[0]] = 1.0
    energies, vectors, krylov_vectors = lanczos_full(
        operator,
        initial,
        100,
        full_ortho=True,
    )
    times = (0.0, 2.4, 4.8, 7.2, 9.6)
    evolved = [
        expm_lanczos(energies, vectors, krylov_vectors, a=-1j * time)
        for time in times
    ]
    norm_error = max(abs(np.linalg.norm(state) - 1.0) for state in evolved)
    fidelities = np.asarray([abs(np.vdot(initial, state)) ** 2 for state in evolved])
    return {
        "valid": basis.Ns == 103_682
        and norm_error < 5e-8
        and fidelities[2] > fidelities[1],
        "residual": norm_error,
        "value": fidelities,
    }


def bose_hubbard_mott_quench() -> dict[str, object]:
    length = 11
    hopping = 0.1
    interaction = 1.0
    basis = boson_basis_1d(L=length, Nb=length, sps=3)
    bonds = [(site, site + 1) for site in range(length - 1)]
    operator = make_hamiltonian(
        [
            ["+-", [[-hopping, left, right] for left, right in bonds]],
            ["-+", [[-hopping, left, right] for left, right in bonds]],
            ["nn", [[0.5 * interaction, site, site] for site in range(length)]],
            ["n", [[-0.5 * interaction, site] for site in range(length)]],
        ],
        basis,
    )
    initial = np.zeros(basis.Ns, dtype=np.complex128)
    initial[basis.index("1" * length)] = 1.0
    energies, vectors, krylov_vectors = lanczos_full(
        operator,
        initial,
        100,
        full_ortho=True,
    )
    times = (0.0, 25.0, 50.0, 100.0, 200.0)
    evolved = [
        expm_lanczos(energies, vectors, krylov_vectors, a=-1j * time)
        for time in times
    ]
    norm_error = max(abs(np.linalg.norm(state) - 1.0) for state in evolved)
    returns = np.asarray([abs(np.vdot(initial, state)) ** 2 for state in evolved])
    return {
        "valid": norm_error < 5e-8 and np.min(returns[1:]) < 0.99,
        "residual": norm_error,
        "value": returns,
    }


def spinful_hubbard_current_quench() -> dict[str, object]:
    length = 10
    basis = spinful_fermion_basis_1d(L=length, Nf=(5, 5))
    bonds = [(site, site + 1) for site in range(length - 1)]
    kinetic = [
        ["+-|", [[-1.0, left, right] for left, right in bonds]],
        ["-+|", [[1.0, left, right] for left, right in bonds]],
        ["|+-", [[-1.0, left, right] for left, right in bonds]],
        ["|-+", [[1.0, left, right] for left, right in bonds]],
    ]
    interaction = ["n|n", [[8.0, site, site] for site in range(length)]]
    unbiased = make_hamiltonian(kinetic + [interaction], basis)
    bias = [
        [
            "n|",
            [[-1.5 if site < length // 2 else 1.5, site] for site in range(length)],
        ],
        [
            "|n",
            [[-1.5 if site < length // 2 else 1.5, site] for site in range(length)],
        ],
    ]
    biased = make_hamiltonian(kinetic + [interaction] + bias, basis)
    ground_energy, initial = eigsh(
        biased.static,
        k=1,
        which="SA",
        v0=deterministic_seed(basis.Ns),
        ncv=20,
        maxiter=5_000,
        tol=1e-9,
    )
    initial = initial[:, 0]
    ground_residual = np.linalg.norm(
        biased.static @ initial - ground_energy[0] * initial
    )
    center = length // 2 - 1
    forward = make_hamiltonian(
        [
            ["+-|", [[1.0, center, center + 1]]],
            ["|+-", [[1.0, center, center + 1]]],
        ],
        basis,
        dtype=np.complex128,
    ).static
    current = -1j * (forward - forward.getH())
    energies, vectors, krylov_vectors = lanczos_full(
        unbiased,
        initial,
        100,
        full_ortho=True,
    )
    evolved = [
        expm_lanczos(energies, vectors, krylov_vectors, a=-1j * time)
        for time in (0.0, 0.5, 1.0, 1.5, 2.0)
    ]
    norm_error = max(abs(np.linalg.norm(state) - 1.0) for state in evolved)
    currents = np.asarray(
        [float(np.real(np.vdot(state, current @ state))) for state in evolved]
    )
    residual = max(float(ground_residual), norm_error)
    return {
        "valid": residual < 3e-7 and np.max(np.abs(currents)) > 1e-3,
        "residual": residual,
        "value": currents,
    }


def conb_dynamical_structure_factor() -> dict[str, object]:
    length = 16
    basis = spin_basis_1d(L=length, pauli=False)
    nearest = [(site, (site + 1) % length) for site in range(length)]
    next_nearest = [(site, (site + 2) % length) for site in range(length)]
    transverse_field = 3.21 * 0.0578838 * 7.0 / 2.88
    operator = make_hamiltonian(
        [
            ["zz", [[-1.0, left, right] for left, right in nearest]],
            ["xx", [[-0.205, left, right] for left, right in nearest]],
            ["yy", [[-0.205, left, right] for left, right in nearest]],
            ["zz", [[0.135, left, right] for left, right in next_nearest]],
            ["xx", [[0.003, left, right] for left, right in next_nearest]],
            ["yy", [[0.003, left, right] for left, right in next_nearest]],
            ["x", [[-transverse_field, site] for site in range(length)]],
        ],
        basis,
    )
    ground_energy, ground = eigsh(
        operator.static,
        k=1,
        which="SA",
        v0=deterministic_seed(basis.Ns),
        ncv=20,
        maxiter=5_000,
        tol=1e-9,
    )
    ground = ground[:, 0]
    spin_q = make_hamiltonian(
        [["z", [[(-1.0) ** site, site] for site in range(length)]]],
        basis,
    )
    source = spin_q.static @ ground
    frequencies = np.linspace(0.0, 4.0, 81)
    spectrum = krylov_spectrum(
        operator,
        source,
        frequencies,
        reference_energy=float(ground_energy[0]),
        broadening=0.05,
        steps=100,
    )
    residual = np.linalg.norm(
        operator.static @ ground - ground_energy[0] * ground
    )
    return {
        "valid": residual < 3e-7
        and np.all(np.isfinite(spectrum))
        and np.min(spectrum) >= -1e-12
        and np.max(spectrum) > 1e-3,
        "residual": residual,
        "value": spectrum,
    }


def triangular_bonds(length_x: int, length_y: int) -> list[tuple[int, int]]:
    bonds: set[tuple[int, int]] = set()
    for y in range(length_y):
        for x in range(length_x):
            site = x + length_x * y
            for next_x, next_y in (
                ((x + 1) % length_x, y),
                (x, (y + 1) % length_y),
                ((x + 1) % length_x, (y + 1) % length_y),
            ):
                neighbor = next_x + length_x * next_y
                bonds.add(tuple(sorted((site, neighbor))))
    return sorted(bonds)


def particle_addition_spectrum() -> dict[str, object]:
    length_x, length_y = 6, 3
    length = length_x * length_y
    source_basis = spinless_fermion_basis_general(length, Nf=6)
    target_basis = spinless_fermion_basis_general(length, Nf=7)
    bonds = triangular_bonds(length_x, length_y)

    def interacting_model(basis):
        return make_hamiltonian(
            [
                ["+-", [[-1.0, left, right] for left, right in bonds]],
                ["-+", [[1.0, left, right] for left, right in bonds]],
                ["nn", [[2.0, left, right] for left, right in bonds]],
            ],
            basis,
        )

    source_hamiltonian = interacting_model(source_basis)
    target_hamiltonian = interacting_model(target_basis)
    source_energy, source_state = eigsh(
        source_hamiltonian.static,
        k=1,
        which="SA",
        v0=deterministic_seed(source_basis.Ns),
        ncv=20,
        maxiter=5_000,
        tol=1e-9,
    )
    source_state = source_state[:, 0]
    transition = target_basis.Op_shift_sector(
        source_basis,
        [["+", [length // 2], 1.0]],
        source_state,
    )
    frequencies = np.linspace(-4.0, 12.0, 81)
    spectrum = krylov_spectrum(
        target_hamiltonian,
        transition,
        frequencies,
        reference_energy=float(source_energy[0]),
        broadening=0.1,
        steps=100,
    )
    residual = np.linalg.norm(
        source_hamiltonian.static @ source_state
        - source_energy[0] * source_state
    )
    return {
        "valid": residual < 3e-7
        and np.linalg.norm(transition) > 1e-6
        and np.all(np.isfinite(spectrum))
        and np.max(spectrum) > 1e-4,
        "residual": residual,
        "value": spectrum,
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
        PaperCase(
            "paper_tfim_fidelity_l16",
            "30690",
            "TFIM degenerate-subspace fidelity scan",
            "L=16;dimension=65536;fields=0.8:0.1:1.2;k=2",
            tfim_fidelity_scan,
        ),
        PaperCase(
            "paper_pxp_revival_l24",
            "PXP",
            "PXP constrained-state revival",
            "L=24;periodic=true;dimension=103682;m=100;times=5",
            pxp_revival,
        ),
        PaperCase(
            "paper_bose_hubbard_quench_l11",
            "BHM",
            "Bose-Hubbard Mott quench",
            "L=11;Nb=11;sps=3;dimension=25653;J/U=0.1;m=100;times=5",
            bose_hubbard_mott_quench,
        ),
        PaperCase(
            "paper_hubbard_current_l10",
            "32600",
            "Spinful-Hubbard current quench",
            "L=10;Nf=(5,5);dimension=63504;U/t=8;m=100;times=5",
            spinful_hubbard_current_quench,
        ),
        PaperCase(
            "paper_conb_dsf_l16",
            "CoNb2O6",
            "CoNb2O6 dynamical structure factor",
            "L=16;dimension=65536;B=7T;m=100;frequencies=81",
            conb_dynamical_structure_factor,
        ),
        PaperCase(
            "paper_particle_addition_6x3",
            "FQAH",
            "Particle-addition spectrum (6x3 proxy)",
            "Lx=6;Ly=3;Nf=6->7;dimensions=18564->31824;m=100;frequencies=81",
            particle_addition_spectrum,
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
