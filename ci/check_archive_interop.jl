#!/usr/bin/env julia

using QuSpin
using SparseArrays

length(ARGS) == 2 ||
    error("usage: check_archive_interop.jl PYTHON_ARCHIVE JULIA_ARCHIVE")
python_archive, julia_archive = ARGS

expected_dense = ComplexF64[
    1 2im 0 0
    -2im 3 0 0
    0 0 4 0.5
    0 0 0.5 5
]
expected_sparse = sparse(ComplexF64[
    0 0.25 0 0
    0.25 0 0 0
    0 0 0 -0.75im
    0 0 0.75im 0
])

from_python = load_zip(python_archive)
from_python.components["dense"] == expected_dense ||
    error("Julia failed to decode Python's dense NPZ entry")
from_python.components["sparse"] isa SparseMatrixCSC ||
    error("Julia did not decode Python's CSR entry as a sparse matrix")
from_python.components["sparse"] == expected_sparse ||
    error("Julia failed to decode Python's CSR indices or values")

operator = QuantumOperator(
    SpinBasis1D(2),
    Dict(
        "dense" => expected_dense,
        "sparse" => expected_sparse,
    );
    matrix_formats=Dict("dense" => :dense, "sparse" => :csc),
)
save_zip(
    julia_archive,
    operator;
    save_basis=false,
    format=:python,
)
