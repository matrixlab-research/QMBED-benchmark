@testset "held-out physical symmetry sectors" begin
    L = 6
    bonds = [(site, mod1(site + 1, L)) for site in 1:L]
    terms = [
        OperatorTerm("zz", [(0.73, pair...) for pair in bonds]),
        OperatorTerm("+-", [(-0.31, pair...) for pair in bonds]),
        OperatorTerm("-+", [(-0.31, pair...) for pair in bonds]),
    ]
    full_basis = SpinBasis1D(L; nup=3, pauli=false)
    full_values = sort(real.(eigvals(Hamiltonian(full_basis, terms))))

    momentum_bases = [
        SpinBasis1D(L; nup=3, pauli=false, kblock=momentum)
        for momentum in 0:(L - 1)
    ]
    @test length.(momentum_bases) == [4, 3, 3, 4, 3, 3]
    @test sum(length, momentum_bases) == length(full_basis)
    @test sort(vcat([
        real.(eigvals(Hamiltonian(basis, terms)))
        for basis in momentum_bases
    ]...)) ≈ full_values atol=4e-12
    @test sort(real.(eigvals(Hamiltonian(momentum_bases[3], terms)))) ≈
        [-0.675, -0.331106984109, 0.641106984109] atol=8e-12

    parity_bases = [
        SpinBasis1D(L; nup=3, pauli=false, pblock=sector)
        for sector in (-1, 1)
    ]
    @test sum(length, parity_bases) == length(full_basis)
    @test sort(vcat([
        real.(eigvals(Hamiltonian(basis, terms)))
        for basis in parity_bases
    ]...)) ≈ full_values atol=4e-12

    selected = SpinBasis1D(L; nup=3, pauli=false, kblock=2)
    @test selected.symmetry.projector isa SparseMatrixCSC
    @test selected.symmetry.projector' * selected.symmetry.projector ≈
        I atol=5e-13
    @test projection_matrix(selected, ComplexF64)' *
          projection_matrix(selected, ComplexF64) ≈ I atol=5e-13
    @test_throws ArgumentError projection_matrix(selected, Float64)
    @test_throws ArgumentError SpinBasis1D(
        L;
        nup=3,
        kblock=1,
        pblock=1,
    )

    fermion_L = 5
    fermion_bonds = [
        (site, mod1(site + 1, fermion_L))
        for site in 1:fermion_L
    ]
    plus = [(-0.83, pair...) for pair in fermion_bonds]
    minus = [(0.83, pair...) for pair in fermion_bonds]
    fermion_terms = [
        OperatorTerm("+-", plus),
        OperatorTerm("-+", minus),
        OperatorTerm("n", [(0.19, site) for site in 1:fermion_L]),
    ]
    full_fermion = SpinlessFermionBasis1D(fermion_L; Nf=2)
    fermion_values = sort(real.(eigvals(Hamiltonian(
        full_fermion,
        fermion_terms,
    ))))
    fermion_blocks = [
        SpinlessFermionBasis1D(fermion_L; Nf=2, kblock=momentum)
        for momentum in 0:(fermion_L - 1)
    ]
    @test length.(fermion_blocks) == fill(2, fermion_L)
    @test sort(vcat([
        real.(eigvals(Hamiltonian(basis, fermion_terms)))
        for basis in fermion_blocks
    ]...)) ≈ fermion_values atol=5e-12
    fermion_projector = projection_matrix(fermion_blocks[2], ComplexF64)
    @test size(fermion_projector) ==
        (2^fermion_L, length(fermion_blocks[2]))
    @test fermion_projector' * fermion_projector ≈ I atol=5e-13

    boson_terms = [
        OperatorTerm("+-", plus),
        OperatorTerm("-+", plus),
        OperatorTerm("n", [(0.19, site) for site in 1:fermion_L]),
    ]
    full_boson = BosonBasis1D(fermion_L; Nb=2, sps=3)
    boson_values = sort(real.(eigvals(Hamiltonian(
        full_boson,
        boson_terms,
    ))))
    boson_blocks = [
        BosonBasis1D(
            fermion_L;
            Nb=2,
            sps=3,
            kblock=momentum,
        )
        for momentum in 0:(fermion_L - 1)
    ]
    @test length.(boson_blocks) == fill(3, fermion_L)
    @test sort(vcat([
        real.(eigvals(Hamiltonian(basis, boson_terms)))
        for basis in boson_blocks
    ]...)) ≈ boson_values atol=5e-12
    boson_projector = projection_matrix(boson_blocks[2], ComplexF64)
    @test size(boson_projector) ==
        (3^fermion_L, length(boson_blocks[2]))
    @test boson_projector' * boson_projector ≈ I atol=5e-13
    reduced_state = normalize(ComplexF64.(1:length(boson_blocks[2])))
    full_boson_state = boson_projector * reduced_state
    @test partial_trace(
        boson_blocks[2],
        reduced_state;
        sub_sys_A=[1, 3],
    ) ≈ partial_trace(
        BosonBasis1D(fermion_L; sps=3),
        full_boson_state;
        sub_sys_A=[1, 3],
    ) atol=7e-13

    spinful_terms = [
        OperatorTerm("+-|", plus),
        OperatorTerm("-+|", minus),
        OperatorTerm("|+-", plus),
        OperatorTerm("|-+", minus),
    ]
    full_spinful = SpinfulFermionBasis1D(fermion_L; Nf=(1, 1))
    spinful_values = sort(real.(eigvals(Hamiltonian(
        full_spinful,
        spinful_terms,
    ))))
    spinful_momentum = [
        SpinfulFermionBasis1D(
            fermion_L;
            Nf=(1, 1),
            kblock=momentum,
        )
        for momentum in 0:(fermion_L - 1)
    ]
    @test length.(spinful_momentum) == fill(5, fermion_L)
    @test sort(vcat([
        real.(eigvals(Hamiltonian(basis, spinful_terms)))
        for basis in spinful_momentum
    ]...)) ≈ spinful_values atol=7e-12
    spinful_projector = projection_matrix(spinful_momentum[2], ComplexF64)
    @test size(spinful_projector) ==
        (4^fermion_L, length(spinful_momentum[2]))
    @test spinful_projector' * spinful_projector ≈ I atol=5e-13

    exchange_bases = [
        SpinfulFermionBasis1D(fermion_L; Nf=(1, 1), sblock=sector)
        for sector in (-1, 1)
    ]
    @test sum(length, exchange_bases) == length(full_spinful)
    @test sort(vcat([
        real.(eigvals(Hamiltonian(basis, spinful_terms)))
        for basis in exchange_bases
    ]...)) ≈ spinful_values atol=7e-12
end
