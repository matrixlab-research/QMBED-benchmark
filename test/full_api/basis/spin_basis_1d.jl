@testset "SpinBasis1D public contract" begin
    basis = SpinBasis1D(4; nup=2, pauli=false)

    @test length(basis) == 6
    @test basis.N == basis.L == 4
    @test basis.Ns == 6
    @test basis.sps == 2
    @test basis.dtype === UInt64
    @test basis.states == states(basis)
    @test basis.blocks[:nup] == 2
    @test "z" in basis.operators
    @test isempty(basis.noncommuting_bits)
    @test !isempty(basis.description)
    @test Set(states(basis)) == Set(UInt64[3, 5, 6, 9, 10, 12])

    for state in states(basis)
        @test state_at(basis, state_index(basis, state)) == state
        formatted = int_to_state(basis, state)
        @test state_to_int(basis, formatted) == state
    end
    @test int_to_state(basis, 5; bracket_notation=false) == "0101"
    @test operator_matrix(basis, "z", [(1.0, 1)]) ==
        Matrix(Hamiltonian(basis, [OperatorTerm("z", [(1.0, 1)])]))
    output = zeros(ComplexF64, basis.Ns, basis.Ns)
    inplace_op!(output, basis, "z", [(1.0, 1)])
    @test output == operator_matrix(basis, "z", [(1.0, 1)])

    projector = projection_matrix(basis)
    @test size(projector) == (16, 6)
    @test projector' * projector == Matrix(I, 6, 6)
    vector = collect(1.0:6.0)
    @test get_vec(basis, vector) == projector * vector
    @test project_from(basis, vector) == projector * vector
    static, dynamic = expanded_form(basis, [:static], [:dynamic])
    @test static == [:static]
    @test dynamic == [:dynamic]
end
