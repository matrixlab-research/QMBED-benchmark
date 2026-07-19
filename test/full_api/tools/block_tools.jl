@testset "symmetry block decomposition and evolution" begin
    terms = [OperatorTerm("z", [(1.0, 1), (0.5, 2)])]
    blocks = [Dict(:nup => sector) for sector in 0:2]
    P, block_H = block_diag_hamiltonian(
        blocks,
        terms,
        Any[],
        SpinBasis1D,
        (2,),
        ComplexF64;
        basis_kwargs=Dict(:pauli => false),
    )
    full_basis = SpinBasis1D(2; pauli=false)
    full_H = Hamiltonian(full_basis, terms)
    expected = ComplexF64[
        -0.75 0 0 0
        0 0.25 0 0
        0 0 -0.25 0
        0 0 0 0.75
    ]
    @test P' * Matrix(full_H) * P ≈ expected atol=2e-16
    @test block_H ≈ expected atol=2e-16

    operator = BlockOps(
        blocks,
        terms,
        Any[],
        SpinBasis1D,
        (2,),
        ComplexF64;
        basis_kwargs=Dict(:pauli => false),
        compute_all_blocks=true,
    )
    @test operator.dtype == ComplexF64
    @test operator.save_previous_data
    @test length(operator.basis_dict) == 3
    @test length(operator.H_dict) == 3
    @test length(operator.P_dict) == 3
    @test operator.static == terms
    @test isempty(operator.dynamic)

    psi = normalize(ComplexF64[1, 2, 3, 4])
    times = [0.0, 0.2, 1.0]
    evolved = evolve(operator, psi, 0.0, times)
    expected_evolved = evolve(full_H, psi, 0.0, times)
    @test evolved ≈ expected_evolved atol=5e-15
    @test collect(evolve(operator, psi, 0.0, times; iterate=true)) ≈
        collect(eachcol(expected_evolved)) atol=5e-15

    @test block_expm(operator, psi) ≈
        exp((-im) .* Matrix(full_H)) * psi atol=5e-15
    exponential_grid = block_expm(
        operator,
        psi;
        start=0.0,
        stop=1.0,
        num=3,
    )
    @test exponential_grid[:, 1] ≈ psi atol=3e-16
    @test exponential_grid[:, 2] ≈
        exp((-0.5im) .* Matrix(full_H)) * psi atol=5e-15
    @test exponential_grid[:, 3] ≈
        exp((-im) .* Matrix(full_H)) * psi atol=5e-15
end
