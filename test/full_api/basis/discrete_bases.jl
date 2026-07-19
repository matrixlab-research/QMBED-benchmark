@testset "boson and fermion basis oracles" begin
    bosons = BosonBasis1D(3; Nb=2)
    @test bosons.N == 3
    @test bosons.Ns == 6
    @test bosons.sps == 3
    @test Set(Int.(bosons.states)) == Set([2, 4, 6, 10, 12, 18])
    @test all(sum(row) == 2 for row in eachrow(bosons.occupations))
    @test bosons.description == "boson lattice basis"
    @test bosons.dtype == UInt64
    @test bosons.operators == ("I", "+", "-", "n", "z")
    @test isempty(bosons.noncommuting_bits)
    @test state_index(bosons, 10) == 4
    @test state_at(bosons, 4) == 10
    @test state_to_int(bosons, int_to_state(bosons, 10)) == 10

    projector = projection_matrix(bosons, ComplexF64)
    @test size(projector) == (27, 6)
    @test projector' * projector ≈ Matrix{ComplexF64}(I, 6, 6)
    state = normalize(ComplexF64[1, 2, 3, 4, 5, 6])
    @test project_from(bosons, state) == projector * state
    @test get_vec(bosons, state) == projector * state
    @test tr(partial_trace(bosons, state; sub_sys_A=[1])) ≈ 1 atol=4e-16
    @test ent_entropy(bosons, state; sub_sys_A=[1])["Sent_A"] >= 0
    @test diag(operator_matrix(bosons, "n", [(1.0, 1)])) ==
        ComplexF64.(bosons.occupations[:, 1])

    spinless = SpinlessFermionBasis1D(4; Nf=2)
    @test spinless.Ns == 6
    @test Set(Int.(spinless.states)) == Set([3, 5, 6, 9, 10, 12])
    @test all(sum(row) == 2 for row in eachrow(spinless.occupations))
    @test !isempty(spinless.noncommuting_bits)
    hopping = operator_matrix(
        spinless,
        "+-",
        [(1.0, 1, 2), (1.0, 2, 1)],
    )
    @test ishermitian(hopping)

    unrestricted = SpinlessFermionBasis1D(2)
    creation = operator_matrix(unrestricted, "+", [(1.0, 2)])
    @test creation[state_index(unrestricted, 3), state_index(unrestricted, 1)] == -1

    spinful = SpinfulFermionBasis1D(2; Nf=(1, 1))
    @test spinful.Ns == 4
    @test all(
        count(digit -> digit & 1 == 1, row) == 1 &&
        count(digit -> digit & 2 == 2, row) == 1
        for row in eachrow(spinful.occupations)
    )
end
