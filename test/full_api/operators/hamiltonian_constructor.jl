@testset "Hamiltonian constructor and dynamic evaluation" begin
    basis = SpinBasis1D(1; pauli=false)
    X = ComplexF64[0 0.5; 0.5 0]
    Z = ComplexF64[-0.5 0; 0 0.5]
    drive = (time, frequency) -> cos(frequency * time)
    H = Hamiltonian(
        Any[Any["z", [(1.0, 1)]]],
        Any[Any[X, drive, (2.0,)]];
        basis,
        dtype=ComplexF64,
    )
    @test Matrix(H) == Z
    @test toarray(H; time=0.0) ≈ Z + X atol=2e-16
    @test toarray(H; time=π / 4) ≈ Z atol=3e-16
    @test diagonal(H; time=0.0) == [-0.5, 0.5]
    @test apply(H, ComplexF64[1, 0]; time=0.0) ≈
        (Z + X) * ComplexF64[1, 0] atol=2e-16
    @test eigh(H; time=0.0)[1] ≈ [-sqrt(0.5), sqrt(0.5)] atol=3e-16

    constant = Hamiltonian(
        Any[Z],
        Any[Any[X, (time,) -> 1.0, ()]];
        basis,
        dtype=ComplexF64,
    )
    psi = ComplexF64[1, 0]
    times = [0.0, 0.2, 1.0]
    evolved = evolve(constant, psi, 0.0, times; max_step=0.001)
    exact = reduce(
        hcat,
        (exp((-im * time) .* (Z + X)) * psi for time in times),
    )
    @test evolved ≈ exact atol=8e-12

    sparse_H = Hamiltonian(
        Any[Any["z", [(1.0, 1)]]],
        Any[Any[X, drive, (2.0,)]];
        basis,
        dtype=ComplexF64,
        static_fmt=:csc,
        dynamic_fmt=:csc,
    )
    @test sparse_H.data isa SparseMatrixCSC
    @test all(first(term) isa SparseMatrixCSC for term in sparse_H.dynamic_terms)
    @test tocsc(sparse_H; time=0.0) isa SparseMatrixCSC
    @test Matrix(tocsc(sparse_H; time=0.0)) ≈ Z + X atol=2e-16
    csr_H = Hamiltonian(
        Any[Any["z", [(1.0, 1)]]],
        Any[];
        basis,
        static_fmt=:csr,
    )
    @test csr_H.data isa SparseMatrixCSR
end
