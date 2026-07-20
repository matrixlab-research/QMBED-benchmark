@testset "held-out sparse Krylov evolution paths" begin
    L = 9
    basis = SpinBasis1D(L; nup=4, pauli=false)
    bonds = [(site, site + 1) for site in 1:(L - 1)]
    terms = [
        OperatorTerm("zz", [(0.57, pair...) for pair in bonds]),
        OperatorTerm("+-", [(0.22, pair...) for pair in bonds]),
        OperatorTerm("-+", [(0.22, pair...) for pair in bonds]),
        OperatorTerm("z", [(0.04(-1)^site, site) for site in 1:L]),
    ]
    H = Hamiltonian(basis, terms; static_fmt=:csc)
    state = normalize(ComplexF64[
        cos(0.12index) + im * sin(0.31index)
        for index in 1:length(basis)
    ])
    times = [0.0, 0.11, 0.46, 1.23]
    evolved = evolve(H, state, 0.0, times; tol=2e-13, krylov_dim=24)
    exact = exp(-im * times[end] * Matrix(H)) * state
    @test evolved[:, end] ≈ exact atol=8e-11
    @test all(
        isapprox(norm(column), 1.0; atol=8e-12)
        for column in eachcol(evolved)
    )
    @test H.data isa SparseMatrixCSC

    csr = tocsr(H)
    dia = DIAMatrix(H.data)
    for storage in (H.data, csr, dia)
        exponential = ExpmMultiplyParallel(storage, -0.63im)
        @test exponential.A === storage
        @test exponential * state ≈
            exp(-0.63im * Matrix(H)) * state atol=8e-11
    end

    small_basis = SpinBasis1D(5; nup=2, pauli=false)
    small_terms = [
        OperatorTerm("zz", [(0.4, site, site + 1) for site in 1:4]),
        OperatorTerm("+-", [(0.3, site, site + 1) for site in 1:4]),
        OperatorTerm("-+", [(0.3, site, site + 1) for site in 1:4]),
    ]
    first_step = Hamiltonian(small_basis, small_terms; static_fmt=:csc)
    second_step = Hamiltonian(
        small_basis,
        [OperatorTerm("z", [(0.21(-1)^site, site) for site in 1:5])];
        static_fmt=:csc,
    )
    durations = [0.37, 0.19]
    floquet = Floquet(
        Dict(
            :H_list => [first_step, second_step],
            :dt_list => durations,
            :T => sum(durations),
        );
        UF=true,
        VF=true,
        thetaF=true,
    )
    exact_unitary =
        exp(-im * durations[2] * Matrix(second_step)) *
        exp(-im * durations[1] * Matrix(first_step))
    @test floquet.UF ≈ exact_unitary atol=8e-11
    @test floquet.UF' * floquet.UF ≈ I atol=9e-11
    @test abs.(floquet.thetaF) ≈ ones(length(small_basis)) atol=9e-11
end
