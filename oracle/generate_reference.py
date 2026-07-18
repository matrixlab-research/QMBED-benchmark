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
from quspin.basis import spin_basis_1d
from quspin.operators import hamiltonian


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
        "cases": [spin_basis_case(), xxz_spectrum_case()],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(args.output)


if __name__ == "__main__":
    main()
