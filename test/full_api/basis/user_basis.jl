@testset "UserBasis callbacks and state filtering" begin
    operators = Dict(
        "n" => [0.0 0.0; 0.0 1.0],
        "+" => [0.0 0.0; 1.0 0.0],
        "x" => ((state, site) -> (xor(state, UInt64(1) << (site - 1)), 1.0)),
    )
    basis = UserBasis(UInt64, 2, operators)
    @test basis.N == 2
    @test basis.Ns == 4
    @test basis.sps == 2
    @test basis.dtype == UInt64
    @test basis.description == "user-defined finite basis"
    @test Set(basis.operators) == Set(("n", "+", "x"))
    @test Set(basis.states) == Set(UInt64[0, 1, 2, 3])
    @test isempty(basis.noncommuting_bits)
    @test state_at(basis, state_index(basis, 3)) == 3
    @test state_to_int(basis, int_to_state(basis, 2)) == 2

    number = operator_matrix(basis, "n", [(1.0, 1)])
    @test diag(number) == ComplexF64[0, 1, 0, 1]
    flip = operator_matrix(basis, "x", [(1.0, 1)])
    @test flip == ComplexF64[
        0 1 0 0
        1 0 0 0
        0 0 0 1
        0 0 1 0
    ]
    output = zeros(ComplexF64, 4, 4)
    inplace_op!(output, basis, "x", [(1.0, 1)])
    @test output == flip

    even = UserBasis(
        UInt64,
        2,
        operators;
        pre_check_state=state -> iseven(count_ones(state)),
        allowed_ops=("n", "x"),
        parity=1,
    )
    @test even.states == UInt64[0, 3]
    @test even.blocks[:parity] == 1
    @test even.operators == ("n", "x")
    projector = projection_matrix(even)
    @test size(projector) == (4, 2)
    @test projector' * projector ≈ Matrix(I, 2, 2)
    state = ComplexF64[inv(sqrt(2)), inv(sqrt(2))]
    @test ent_entropy(even, state; sub_sys_A=[1])["Sent_A"] ≈ log(2) atol=3e-16
    @test make_basis!(even) === even
    @test make_basis_blocks(even) == [1:2]
    @test normalization(even, 0) == 1
    @test get_amp(even, 0) == 1
    @test representative(even, 3) == 3
end
