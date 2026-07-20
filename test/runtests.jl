using Test
using QuSpin
using LinearAlgebra
using SparseArrays

include("timing.jl")

@testset "QuSpin full API verification" begin
    timed_include("full_api/basis/integer_helpers.jl", "api")
    timed_include("full_api/basis/photon_helpers.jl", "api")
    timed_include("full_api/basis/spin_basis_1d.jl", "api")
    timed_include("full_api/basis/discrete_bases.jl", "mixed")
    timed_include("full_api/basis/symmetry_sectors.jl", "integration")
    timed_include("full_api/basis/composite_bases.jl", "integration")
    timed_include("full_api/basis/user_basis.jl", "integration")
    timed_include("full_api/operators/algebra.jl", "mixed")
    timed_include("full_api/operators/exp_op.jl", "integration")
    timed_include("full_api/operators/hamiltonian.jl", "integration")
    timed_include("full_api/operators/hamiltonian_constructor.jl", "integration")
    timed_include("full_api/operators/archive.jl", "integration")
    timed_include("full_api/operators/quantum_linear_operator.jl", "integration")
    timed_include("full_api/operators/quantum_operator.jl", "integration")
    timed_include("full_api/operators/storage_formats.jl", "mixed")
    timed_include("full_api/operators/matrix_free_regression.jl", "integration")
    timed_include("full_api/tools/evolution.jl", "integration")
    timed_include("full_api/tools/expm_multiply_parallel.jl", "mixed")
    timed_include("full_api/tools/floquet_time_vector.jl", "api")
    timed_include("full_api/tools/floquet.jl", "integration")
    timed_include("full_api/tools/block_tools.jl", "integration")
    timed_include("full_api/tools/lanczos.jl", "integration")
    timed_include("full_api/tools/diag_ensemble.jl", "integration")
    timed_include("full_api/tools/measurements.jl", "integration")
    timed_include("full_api/tools/misc.jl", "mixed")
    timed_include("full_api/tools/sparse_krylov_paths.jl", "integration")
    timed_include("full_api/integration/general_basis_hamiltonians.jl", "integration")
    timed_include("full_api/integration/paper_workflows_sparse.jl", "integration")
    timed_include("full_api/integration/ed_workflow_catalog.jl", "integration")
    timed_include("full_api/integration/algorithm_regressions.jl", "integration")
end
