using LinearAlgebra
using QuSpin

@testset "dynamic Hamiltonian end-to-end semantics" begin
    basis = SpinBasis1D(2; pauli=false)
    X1 = operator_matrix(basis, "x", [(1.0, 1)])
    Z2 = operator_matrix(basis, "z", [(1.0, 2)])
    exchange = operator_matrix(basis, "+-", [(1.0, 1, 2)]) +
        operator_matrix(basis, "-+", [(1.0, 1, 2)])
    drive(time, frequency, phase) = sin(frequency * time + phase)
    H = Hamiltonian(
        Any[0.4Z2],
        Any[Any[0.7X1, drive, (1.3, 0.2)]];
        basis,
        dtype=ComplexF64,
        check_herm=false,
        check_symm=false,
        check_pcon=false,
    )
    K = Hamiltonian(
        Any[0.2exchange],
        Any[Any[0.1Z2, (time,) -> time, ()]];
        basis,
        dtype=ComplexF64,
        check_herm=false,
        check_symm=false,
        check_pcon=false,
    )
    rotation = exp(-0.17im .* exchange)

    for time in range(0.0, 1.0; length=7)
        Ht = toarray(H; time)
        Kt = toarray(K; time)
        @test toarray(rotate_by(H, rotation); time) ≈ rotation' * Ht * rotation atol=3e-14
        @test toarray(H + K; time) ≈ Ht + Kt atol=3e-14
        @test toarray(H * K; time) ≈ Ht * Kt atol=4e-14
        @test toarray(commutator(H, K); time) ≈ Ht * Kt - Kt * Ht atol=5e-14
    end

    initial_states = ComplexF64[1 0; 0 1; 0 0; 0 0]
    times = [0.0, 0.15, 0.4]
    batch = evolve(
        H,
        initial_states,
        0.0,
        times;
        eom=:SE,
        max_step=0.01,
        rtol=1e-10,
        atol=1e-12,
    )
    @test size(batch) == (4, 2, 3)
    @test batch[:, 1, :] ≈ evolve(
        H,
        initial_states[:, 1],
        0.0,
        times;
        eom=:SE,
        max_step=0.01,
        rtol=1e-10,
        atol=1e-12,
    ) atol=5e-10
end

@testset "expanded forms and dynamic matrix arithmetic" begin
    drive(time, frequency) = cos(frequency * time)
    basis = SpinBasis1D(2; pauli=false)
    static = Any[Any["xy", [(0.6, 1, 2)]]]
    dynamic = Any[Any["yx", [(-0.4, 1, 2)], drive, (0.8,)]]
    expanded_static, expanded_dynamic =
        expanded_form(basis, static, dynamic)
    original = Hamiltonian(
        static,
        dynamic;
        basis,
        dtype=ComplexF64,
        check_herm=false,
    )
    expanded = Hamiltonian(
        expanded_static,
        expanded_dynamic;
        basis,
        dtype=ComplexF64,
        check_herm=false,
    )
    matrix = ComplexF64[
        0.2 0.1 0 0
        -0.3 0.4 0.2im 0
        0 -0.2im -0.1 0.5
        0 0 -0.4 0.3
    ]
    for time in (0.0, 0.31, 1.1)
        dense = toarray(original; time)
        @test toarray(expanded; time) ≈ dense atol=5e-14
        @test toarray(original + matrix; time) ≈ dense + matrix atol=5e-14
        @test toarray(matrix * original; time) ≈ matrix * dense atol=8e-14
        @test toarray(original^2; time) ≈ dense^2 atol=8e-14
    end
end

@testset "symmetry-changing operator integration" begin
    source = SpinBasis1D(6; nup=3, pauli=false, kblock=0)
    target = SpinBasis1D(6; nup=3, pauli=false, kblock=2)
    op_list = [
        ("z", [site], cis(2π * 2 * (site - 1) / 6))
        for site in 1:6
    ]
    initial = normalize(ComplexF64.(1:length(source)))
    shifted = op_shift_sector(target, source, op_list, initial)

    full = SpinBasis1D(6; pauli=false)
    full_operator = sum(
        operator_matrix(full, opstring, [(coefficient, sites...)])
        for (opstring, sites, coefficient) in op_list
    )
    explicit =
        projection_matrix(target, ComplexF64; sparse=true)' *
        full_operator *
        projection_matrix(source, ComplexF64; sparse=true) *
        initial
    @test shifted ≈ explicit atol=4e-14
    @test norm(shifted) ≈ norm(explicit) atol=4e-14
end
