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

    constrained(state) =
        count_ones(state) == 2 &&
        iszero(state & (state << 1))
    serial = UserBasis(
        UInt64,
        9,
        Dict("n" => operators["n"]);
        pre_check_state=constrained,
        parallel=false,
    )
    threaded = UserBasis(
        UInt64,
        9,
        Dict("n" => operators["n"]);
        pre_check_state=constrained,
        parallel=true,
    )
    @test states(threaded) == states(serial)
    @test threaded.blocks[:parallel]

    deferred = UserBasis(
        UInt64,
        5,
        Dict("n" => operators["n"]);
        pre_check_state=state -> count_ones(state) == 2,
        _make_basis=false,
    )
    @test length(deferred) == 1
    @test deferred.blocks[:made_basis] == false
    @test make_basis!(deferred) === deferred
    @test length(deferred) == 10
    @test deferred.blocks[:made_basis] == true
end
@testset "UserBasis callback symmetry reduction" begin
    N = 5
    reflect(state, N, args) = begin
        output = zero(UInt64)
        for site in 0:(N - 1)
            output |= ((UInt64(state) >> site) & 1) << (N - site - 1)
        end
        output
    end
    basis = UserBasis(
        UInt64,
        N,
        Dict(
            "x" => ComplexF64[0 1; 1 0],
            "z" => ComplexF64[1 0; 0 -1],
        );
        parity=(reflect, 2, 0, ()),
        block_order=[:parity],
    )
    P = basis.base.symmetry.projector
    @test P' * P ≈ Matrix{ComplexF64}(I, length(basis), length(basis)) atol=4e-14
    @test basis.blocks[:parity] == 0
    @test basis.blocks[:parity_period] == 2
end

@testset "UserBasis cross-sector callback action" begin
    flip(state, site) =
        (xor(state, UInt64(1) << (site - 1)), 1.0)
    source = UserBasis(
        UInt64,
        3,
        Dict("x" => flip);
        states=UInt64[0, 2, 4, 6],
    )
    target = UserBasis(
        UInt64,
        3,
        Dict("x" => flip);
        states=UInt64[1, 3, 5, 7],
    )
    state = normalize(ComplexF64[1, 2, 3im, -0.5])
    shifted = op_shift_sector(
        target,
        source,
        [("x", [1], -0.4)],
        state,
    )
    @test shifted ≈ -0.4state atol=2e-16

    matrix_elements, bras, kets = op_bra_ket(
        source,
        "x",
        [1],
        -0.4,
        ComplexF64,
        UInt64[0, 2, 4],
    )
    @test matrix_elements == fill(-0.4 + 0im, 3)
    @test bras == UInt64[1, 3, 5]
    @test kets == UInt64[0, 2, 4]
end

@testset "UserBasis branching direct-CSC semantics" begin
    first_local = ComplexF64[
        1.0 0.5im 0.0
        -0.25 1.0 0.75
        0.5 0.0 -1.0im
    ]
    second_local = ComplexF64[
        0.0 1.0 0.25im
        0.5 0.0 -0.75
        1.0im 0.25 0.0
    ]
    cycle_factor = 0.5 - 0.25im
    cycle = function (state, site)
        weight = UInt64(3)^(site - 1)
        old = Int((state ÷ weight) % UInt64(3))
        new = mod(old + 1, 3)
        return UInt64(
            Int128(state) + Int128(new - old) * Int128(weight),
        ), cycle_factor
    end
    cycle_local = zeros(ComplexF64, 3, 3)
    for old in 0:2
        cycle_local[mod(old + 1, 3) + 1, old + 1] = cycle_factor
    end

    selected_states = UInt64[0, 1, 2, 3, 4, 6, 7, 8]
    basis = UserBasis(
        UInt64,
        2,
        Dict('a' => first_local, 'b' => second_local, 'c' => cycle);
        sps=3,
        states=selected_states,
        allowed_ops=('a', 'b', 'c'),
    )
    coefficient_ab = 0.7 + 0.2im
    coefficient_ac = -0.4 + 0.1im
    assembled = operator_matrix(
        basis,
        "ab",
        [(coefficient_ab, 1, 1)];
        sparse=true,
    ) + operator_matrix(
        basis,
        "ac",
        [(coefficient_ac, 1, 1)];
        sparse=true,
    )

    local_reference = coefficient_ab * first_local * second_local +
        coefficient_ac * first_local * cycle_local
    state_lookup = Dict(state => index for (index, state) in enumerate(selected_states))
    expected = zeros(ComplexF64, length(basis), length(basis))
    for (column, state) in enumerate(selected_states)
        old = Int(state % UInt64(3))
        spectator = state - UInt64(old)
        for new in 0:2
            row = get(state_lookup, spectator + UInt64(new), 0)
            iszero(row) && continue
            expected[row, column] += local_reference[new + 1, old + 1]
        end
    end

    @test assembled isa SparseMatrixCSC{ComplexF64,Int}
    @test Matrix(assembled) ≈ expected atol=3e-15
    @test all(
        issorted(rowvals(assembled)[nzrange(assembled, column)])
        for column in axes(assembled, 2)
    )
end
