#!/usr/bin/env python3
"""Verify bidirectional QuSpin Python/Julia operator-archive interchange."""

from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from quspin.operators import load_zip, quantum_operator, save_zip
from scipy import sparse


def components():
    dense = np.array(
        [
            [1, 2j, 0, 0],
            [-2j, 3, 0, 0],
            [0, 0, 4, 0.5],
            [0, 0, 0.5, 5],
        ],
        dtype=np.complex128,
    )
    csr = sparse.csr_matrix(
        np.array(
            [
                [0, 0.25, 0, 0],
                [0.25, 0, 0, 0],
                [0, 0, 0, -0.75j],
                [0, 0, 0.75j, 0],
            ],
            dtype=np.complex128,
        )
    )
    return dense, csr


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", default=".")
    args = parser.parse_args()

    dense, csr = components()
    operator = quantum_operator(
        {"dense": [dense], "sparse": [csr]},
        basis=None,
        dtype=np.complex128,
        copy=False,
        check_symm=False,
        check_herm=False,
        check_pcon=False,
    )

    helper = Path(__file__).with_suffix(".jl")
    with tempfile.TemporaryDirectory() as directory:
        python_archive = Path(directory) / "python.zip"
        julia_archive = Path(directory) / "julia.zip"
        save_zip(str(python_archive), operator, save_basis=False)
        subprocess.run(
            [
                "julia",
                "--startup-file=no",
                f"--project={args.project}",
                str(helper),
                str(python_archive),
                str(julia_archive),
            ],
            check=True,
        )

        from_julia = load_zip(str(julia_archive))
        np.testing.assert_array_equal(
            from_julia._quantum_operator["dense"],
            dense,
        )
        julia_sparse = from_julia._quantum_operator["sparse"]
        if not sparse.isspmatrix_csc(julia_sparse):
            raise AssertionError(
                "Python did not decode Julia's sparse entry as CSC"
            )
        np.testing.assert_array_equal(julia_sparse.toarray(), csr.toarray())

    print("Python CSR -> Julia CSC and Julia CSC -> Python CSC: OK")


if __name__ == "__main__":
    main()
