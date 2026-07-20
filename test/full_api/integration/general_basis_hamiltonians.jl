@testset "held-out general-basis Hamiltonians" begin
    single_particle = SpinlessFermionBasis1D(4; Nf=1)
    bonds = [(1, 2), (2, 3), (3, 4)]
    hopping = [
        OperatorTerm("+-", [(-0.6, pair...) for pair in bonds]),
        OperatorTerm("-+", [(0.6, pair...) for pair in bonds]),
    ]
    fermion_H = Hamiltonian(single_particle, hopping; static_fmt=:csc)
    analytic = sort([
        -1.2cos(mode * π / 5)
        for mode in 1:4
    ])
    @test sort(real.(eigvals(fermion_H))) ≈ analytic atol=2e-14
    @test ishermitian(fermion_H)

    bosons = BosonBasis1D(3; Nb=2, sps=4)
    boson_terms = [
        OperatorTerm("+-", [(-0.4, 1, 2), (-0.7, 2, 3)]),
        OperatorTerm("-+", [(-0.4, 1, 2), (-0.7, 2, 3)]),
        OperatorTerm("n", [(0.2, 1), (-0.1, 2), (0.3, 3)]),
    ]
    boson_H = Hamiltonian(bosons, boson_terms; static_fmt=:csc)
    @test size(boson_H) == (6, 6)
    @test ishermitian(boson_H)
    @test tr(Matrix(boson_H)) ≈ 1.6 atol=3e-15
    @test norm(Matrix(boson_H)) ≈ 2.71293199325 atol=8e-11
    @test sort(real.(eigvals(boson_H))) ≈ [
        -1.48167314497,
        -0.517122144247,
        0.176285571762,
        0.447428856476,
        1.14083657248,
        1.83424428849,
    ] atol=8e-11

    spinful = SpinfulFermionBasis1D(3; Nf=(1, 1))
    spinful_bonds = [(1, 2), (2, 3)]
    spinful_terms = [
        OperatorTerm("+-|", [(-0.9, pair...) for pair in spinful_bonds]),
        OperatorTerm("-+|", [(0.9, pair...) for pair in spinful_bonds]),
        OperatorTerm("|+-", [(-0.9, pair...) for pair in spinful_bonds]),
        OperatorTerm("|-+", [(0.9, pair...) for pair in spinful_bonds]),
        OperatorTerm("n|n", [(1.7, site, site) for site in 1:3]),
    ]
    hubbard = Hamiltonian(spinful, spinful_terms; static_fmt=:csc)
    @test size(hubbard) == (9, 9)
    @test ishermitian(hubbard)
    @test sort(real.(eigvals(hubbard))) ≈ [
        -2.07100792258,
        -1.27279220614,
        -0.680522786501,
        0.0,
        0.393747310799,
        1.27279220614,
        1.7,
        2.3805227865,
        3.37726061178,
    ] atol=8e-11

    shift = ComplexF64[
        0 1 0
        0 0 1
        1 0 0
    ]
    qutrits = UserBasis(UInt64, 2, Dict('s' => shift); sps=3)
    custom = Hamiltonian(
        qutrits,
        [
            OperatorTerm("s", [(1.0, 1), (0.4, 2)]),
            OperatorTerm("s", [(1.0, 1), (0.4, 2)]),
        ];
        static_fmt=:csc,
        check_herm=false,
    )
    @test size(custom) == (9, 9)
    @test nnz(custom.data) == 18

    @test check_hermitian(spinful, spinful_terms)
    @test check_pcon(spinful, spinful_terms)
    @test !check_pcon(spinful, [OperatorTerm("+|", [(1.0, 1)])])
end
