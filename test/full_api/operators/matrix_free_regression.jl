@testset "held-out matrix-free operator regression" begin
    L = 10
    basis = SpinBasis1D(L; nup=5, pauli=false)
    bonds = [(site, site + 1) for site in 1:(L - 1)]
    terms = [
        OperatorTerm("zz", [(0.61, pair...) for pair in bonds]),
        OperatorTerm("+-", [(-0.27, pair...) for pair in bonds]),
        OperatorTerm("-+", [(-0.27, pair...) for pair in bonds]),
        OperatorTerm("z", [(0.03site - 0.11, site) for site in 1:L]),
    ]
    operator = QuantumLinearOperator(basis, terms)
    reference = Hamiltonian(basis, terms; static_fmt=:csc)
    vector = normalize(ComplexF64[
        sin(0.13index) + im * cos(0.29index)
        for index in 1:length(basis)
    ])

    @test operator.explicit_data === nothing
    @test operator * vector ≈ reference * vector atol=4e-14
    @test apply(operator, vector) ≈ reference * vector atol=4e-14
    @test operator.explicit_data === nothing
    values, vectors = eigsh(
        operator;
        k=3,
        which=:SA,
        tol=1e-11,
        maxiter=2_000,
    )
    @test operator * vectors ≈ vectors * Diagonal(values) atol=8e-10
    @test values ≈ eigsh(
        reference;
        k=3,
        which=:SA,
        tol=1e-11,
        maxiter=2_000,
        return_eigenvectors=false,
    ) atol=3e-10
    @test operator.explicit_data === nothing

    diagonal_shift = collect(range(-0.2, 0.3; length=length(basis)))
    set_diagonal!(operator, diagonal_shift)
    @test operator * vector ≈
        reference * vector + diagonal_shift .* vector atol=4e-14

    periodic_terms = [
        OperatorTerm("zz", [
            (0.4, site, mod1(site + 1, 6))
            for site in 1:6
        ]),
    ]
    reduced_basis = SpinBasis1D(
        6;
        nup=3,
        pauli=false,
        kblock=2,
    )
    reduced_operator = QuantumLinearOperator(reduced_basis, periodic_terms)
    reduced_matrix = Hamiltonian(
        reduced_basis,
        periodic_terms;
        static_fmt=:csc,
    )
    reduced_vector = normalize(ComplexF64.(1:length(reduced_basis)))
    @test reduced_operator * reduced_vector ≈
        reduced_matrix * reduced_vector atol=3e-14
    @test reduced_operator.explicit_data === nothing

    bosons = BosonBasis1D(3; Nb=2, sps=4)
    boson_terms = [
        OperatorTerm("+-", [(-0.4, 1, 2), (-0.7, 2, 3)]),
        OperatorTerm("-+", [(-0.4, 1, 2), (-0.7, 2, 3)]),
        OperatorTerm("n", [(0.2, 1), (-0.1, 2), (0.3, 3)]),
    ]
    boson_operator = QuantumLinearOperator(bosons, boson_terms)
    boson_matrix = Hamiltonian(bosons, boson_terms; static_fmt=:csc)
    boson_vector = normalize(ComplexF64.(1:length(bosons)))
    @test boson_operator * boson_vector ≈
        boson_matrix * boson_vector atol=3e-14
    @test boson_operator.explicit_data === nothing

    spinful = SpinfulFermionBasis1D(2; Nf=(1, 1))
    spinful_terms = [
        OperatorTerm("+-|", [(-0.8, 1, 2)]),
        OperatorTerm("-+|", [(0.8, 1, 2)]),
        OperatorTerm("|+-", [(-0.8, 1, 2)]),
        OperatorTerm("|-+", [(0.8, 1, 2)]),
        OperatorTerm("n|n", [(1.3, 1, 1), (1.3, 2, 2)]),
    ]
    spinful_operator = QuantumLinearOperator(spinful, spinful_terms)
    spinful_matrix = Hamiltonian(spinful, spinful_terms; static_fmt=:csc)
    spinful_vector = normalize(ComplexF64.(1:length(spinful)))
    @test spinful_operator * spinful_vector ≈
        spinful_matrix * spinful_vector atol=3e-14
    @test spinful_operator.explicit_data === nothing

    flip = ComplexF64[0 1; 1 0]
    custom = UserBasis(UInt64, 2, Dict('x' => flip))
    custom_terms = [
        OperatorTerm("x", [(0.7, 1), (-0.2, 2)]),
    ]
    custom_operator = QuantumLinearOperator(custom, custom_terms)
    custom_matrix = Hamiltonian(custom, custom_terms; static_fmt=:csc)
    custom_vector = normalize(ComplexF64[1, 2im, -0.5, 0.3im])
    @test custom_operator * custom_vector ≈
        custom_matrix * custom_vector atol=3e-14
    @test custom_operator.explicit_data === nothing

    @test_throws ArgumentError QuantumLinearOperator(
        SpinBasis1D(4),
        [OperatorTerm("+", [(1.0, 1)])],
    )

    oracle_basis = SpinBasis1D(5; nup=2, pauli=false)
    oracle_terms = [
        OperatorTerm("zz", [(0.41, site, site + 1) for site in 1:4]),
        OperatorTerm("+-", [(-0.23, site, site + 1) for site in 1:4]),
        OperatorTerm("-+", [(-0.23, site, site + 1) for site in 1:4]),
        OperatorTerm("z", [(0.07(site - 3), site) for site in 1:5]),
    ]
    oracle_vector = ComplexF64[
        sin(0.21(Int(state) + 1)) + im * cos(0.16(Int(state) + 1))
        for state in oracle_basis.states
    ]
    oracle_vector ./= norm(oracle_vector)
    expected_by_state = Dict(
        3 => -0.0665209031487 - 0.0405988505859im,
        5 => -0.276924462538 - 0.142130992486im,
        6 => -0.136981227391 - 0.0355895966618im,
        9 => -0.146144312938 + 0.0422814411285im,
        10 => -0.194100864953 + 0.0969435651939im,
        12 => 0.0231782618918 + 0.0697964089521im,
        17 => -0.00776181931881 + 0.0703519825611im,
        18 => 0.0860051640002 + 0.186468132062im,
        20 => 0.101599306154 + 0.165695567682im,
        24 => -0.0408895873313 - 0.0139569027972im,
    )
    expected_action = ComplexF64[
        expected_by_state[Int(state)]
        for state in oracle_basis.states
    ]
    oracle_operator = QuantumLinearOperator(oracle_basis, oracle_terms)
    @test oracle_operator * oracle_vector ≈ expected_action atol=8e-12
    @test oracle_operator.explicit_data === nothing
end
