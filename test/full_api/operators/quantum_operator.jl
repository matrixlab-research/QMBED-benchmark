@testset "QuantumOperator full protocol and oracle" begin
    basis = SpinBasis1D(2; pauli=false)
    x_terms = [
        OperatorTerm("+-", [(1.0, 1, 2)]),
        OperatorTerm("-+", [(1.0, 1, 2)]),
    ]
    z_terms = [OperatorTerm("z", [(1.0, 1)])]
    operator = QuantumOperator(
        basis,
        Dict(:x => x_terms, :z => z_terms),
    )
    pars = Dict(:x => 0.7, :z => 0.3)
    expected = 0.7Matrix(Hamiltonian(basis, x_terms)) +
        0.3Matrix(Hamiltonian(basis, z_terms))

    @test eigvals(operator; pars) ≈ [
        -0.7158910531638175,
        -0.15,
        0.15,
        0.7158910531638175,
    ] atol=2e-15
    @test isquantum_operator(operator)
    @test operator.Ns == 4
    @test operator.shape == (4, 4) == operator.get_shape
    @test operator.ndim == 2
    @test operator.dtype == Float64
    @test operator.basis === basis
    @test operator.is_dense
    @test get_operators(operator, :x) == Matrix(Hamiltonian(basis, x_terms))
    @test toarray(operator; pars) == expected
    @test todense(operator; pars) == expected
    @test Matrix(tocsc(operator; pars)) == expected
    @test Matrix(tocsr(operator; pars)) == expected
    @test diagonal(operator; pars) == diag(expected)
    @test tr(operator; pars) == 0.0
    @test toarray(operator.H; pars) == expected'
    @test toarray(operator.T; pars) == transpose(expected)

    vector = normalize(ComplexF64[1.0, -0.3im, 0.4 + 0.2im, -0.7])
    @test apply(operator, vector; pars) == expected * vector
    out = zeros(ComplexF64, 4)
    @test apply(operator, vector; pars, out, a=0.25) === out
    @test out ≈ 0.25expected * vector atol=3e-16
    @test right_apply(operator, vector; pars) ==
        vec(transpose(vector) * expected)
    @test toarray(conj(operator); pars) == conj(expected)
    @test toarray(transpose(operator); pars) == transpose(expected)
    @test toarray(adjoint(operator); pars) == expected'
    @test copy(operator) !== operator
    @test astype(operator, ComplexF64).dtype == ComplexF64

    values, vectors = eigh(operator; pars)
    @test expected * vectors ≈ vectors * Diagonal(values) atol=2e-14
    selected, selected_vectors = eigsh(operator; pars, k=2, which=:SA)
    @test selected == values[1:2]
    @test size(selected_vectors) == (4, 2)
    expectation = dot(vector, expected * vector)
    @test expt_value(operator, vector; pars) ≈ expectation
    @test matrix_ele(operator, vector, vector; pars) ≈ expectation
    @test quant_fluct(operator, vector; pars) ≈
        dot(vector, expected^2 * vector) - expectation^2 atol=3e-15
    @test Matrix(tohamiltonian(operator; pars)) == expected
    @test Matrix(aslinearoperator(operator; pars)) == expected
    @test update_matrix_formats!(operator, Dict(:x => :dense)) === operator
end
