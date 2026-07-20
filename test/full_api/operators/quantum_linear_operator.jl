@testset "QuantumLinearOperator full protocol" begin
    basis = SpinBasis1D(3; nup=1, pauli=false)
    terms = [
        OperatorTerm("+-", [(0.65, 1, 2), (-0.2, 2, 3)]),
        OperatorTerm("-+", [(0.65, 1, 2), (-0.2, 2, 3)]),
        OperatorTerm("z", [(0.17, 1), (0.31, 2), (-0.23, 3)]),
    ]
    diagonal_shift = [0.11, -0.07, 0.29]
    operator = QuantumLinearOperator(
        basis,
        terms;
        diagonal=diagonal_shift,
    )
    base = Matrix(Hamiltonian(basis, terms))
    matrix = base + Diagonal(diagonal_shift)
    oracle_operator = QuantumLinearOperator(basis, terms)
    @test eigvals(Hermitian(Matrix(oracle_operator))) ≈ [
        -0.6126648705143887,
        -0.3006894850870097,
        0.7883543556013984,
    ] atol=2e-15

    @test isquantum_LinearOperator(operator)
    @test operator.Ns == 3
    @test operator.shape == (3, 3) == operator.get_shape
    @test operator.ndim == 2
    @test operator.dtype == Float64
    @test operator.basis === basis
    @test operator.static_list == terms
    @test operator.diagonal == diagonal_shift
    @test Matrix(operator) == matrix
    @test Matrix(operator.H) == matrix'
    @test Matrix(operator.T) == transpose(matrix)

    vector = normalize(ComplexF64[1.0, -0.4im, 0.7 + 0.2im])
    @test operator * vector ≈ matrix * vector atol=3e-16
    @test apply(operator, vector) ≈ matrix * vector atol=3e-16
    out = zeros(ComplexF64, 3)
    @test apply(operator, vector; out, a=-0.5) === out
    @test out ≈ -0.5matrix * vector atol=3e-16
    @test right_apply(operator, vector) == vec(transpose(vector) * matrix)
    @test Matrix(conj(operator)) == conj(matrix)
    @test Matrix(transpose(operator)) == transpose(matrix)
    @test Matrix(adjoint(operator)) == matrix'
    @test copy(operator) !== operator

    values, vectors = eigsh(operator; k=2, which=:SA)
    @test matrix * vectors ≈ vectors * Diagonal(values) atol=2e-14
    expectation = dot(vector, matrix * vector)
    @test expt_value(operator, vector) ≈ expectation
    @test matrix_ele(operator, vector, vector) ≈ expectation
    @test quant_fluct(operator, vector) ≈
        dot(vector, matrix^2 * vector) - expectation^2 atol=3e-15

    set_diagonal!(operator, zeros(3))
    @test operator.diagonal == zeros(3)
    @test Matrix(operator) == base
end
