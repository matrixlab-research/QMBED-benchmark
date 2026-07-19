using Test
using QuSpin
using LinearAlgebra

@testset "QuSpin full API verification" begin
    include("full_api/basis/integer_helpers.jl")
    include("full_api/basis/photon_helpers.jl")
    include("full_api/basis/spin_basis_1d.jl")
    include("full_api/basis/discrete_bases.jl")
    include("full_api/basis/composite_bases.jl")
    include("full_api/basis/user_basis.jl")
    include("full_api/operators/algebra.jl")
    include("full_api/operators/exp_op.jl")
    include("full_api/operators/hamiltonian.jl")
    include("full_api/operators/hamiltonian_constructor.jl")
    include("full_api/operators/archive.jl")
    include("full_api/operators/quantum_linear_operator.jl")
    include("full_api/operators/quantum_operator.jl")
    include("full_api/tools/evolution.jl")
    include("full_api/tools/expm_multiply_parallel.jl")
    include("full_api/tools/floquet_time_vector.jl")
    include("full_api/tools/floquet.jl")
    include("full_api/tools/block_tools.jl")
    include("full_api/tools/lanczos.jl")
    include("full_api/tools/diag_ensemble.jl")
    include("full_api/tools/measurements.jl")
    include("full_api/tools/misc.jl")
end
