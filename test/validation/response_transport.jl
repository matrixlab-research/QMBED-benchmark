@testset "LKM validation: response and transport" begin
    @testset "#32600 spinful-Hubbard current quench" begin
        L = 6
        basis = SpinfulFermionBasis1D(L; Nf=(L ÷ 2, L ÷ 2))
        kinetic = _lkm_fermion_hopping_terms(L; spinful=true)
        interaction = OperatorTerm(
            "n|n",
            [(4.0, site, site) for site in 1:L],
        )
        H = Hamiltonian(basis, vcat(kinetic, [interaction]); static_fmt=:csc)
        bias_terms = [
            OperatorTerm(
                "n|",
                [(site <= L ÷ 2 ? -1.5 : 1.5, site) for site in 1:L],
            ),
            OperatorTerm(
                "|n",
                [(site <= L ÷ 2 ? -1.5 : 1.5, site) for site in 1:L],
            ),
        ]
        biased = Hamiltonian(
            basis,
            vcat(kinetic, [interaction], bias_terms);
            static_fmt=:csc,
        )
        _, initial = _lkm_ground_state(biased)

        imbalance = Matrix(Hamiltonian(basis, [
            OperatorTerm(
                "n|",
                [(site <= L ÷ 2 ? 1.0 : -1.0, site) for site in 1:L],
            ),
            OperatorTerm(
                "|n",
                [(site <= L ÷ 2 ? 1.0 : -1.0, site) for site in 1:L],
            ),
        ]))
        center = L ÷ 2
        forward =
            operator_matrix(basis, "+-|", [(1.0, center, center + 1)]; sparse=true) +
            operator_matrix(basis, "|+-", [(1.0, center, center + 1)]; sparse=true)
        current = -im .* (forward - forward')
        @test ishermitian(current)

        times = collect(0.0:0.1:1.5)
        evolved = evolve(H, initial, 0.0, times; tol=1e-10, krylov_dim=80)
        imbalance_values = [real(dot(state, imbalance * state)) for state in eachcol(evolved)]
        currents = [real(dot(state, current * state)) for state in eachcol(evolved)]
        @test abs(first(imbalance_values)) > 1
        @test maximum(abs.(currents)) > 0.05
        @test abs(last(imbalance_values)) < abs(first(imbalance_values))
    end

    @testset "#31254 strong-coupling doublon survival" begin
        L = 6
        basis = SpinfulFermionBasis1D(L; Nf=(1, 1))
        kinetic = _lkm_fermion_hopping_terms(L; spinful=true)
        function hubbard(interaction)
            return Hamiltonian(
                basis,
                vcat(
                    kinetic,
                    [OperatorTerm(
                        "n|n",
                        [(interaction, site, site) for site in 1:L],
                    )],
                );
                static_fmt=:csc,
            )
        end
        doublon = Matrix(Hamiltonian(
            basis,
            [OperatorTerm("n|n", [(1.0, site, site) for site in 1:L])],
        ))
        center = L ÷ 2
        encoded = UInt64(3) << (2 * (center - 1))
        initial = zeros(ComplexF64, length(basis))
        initial[state_index(basis, encoded)] = 1
        times = collect(0.0:0.1:2.0)
        free = evolve(hubbard(0.0), initial, 0.0, times; tol=1e-11, krylov_dim=50)
        bound = evolve(hubbard(10.0), initial, 0.0, times; tol=1e-11, krylov_dim=50)
        free_survival = real(dot(free[:, end], doublon * free[:, end]))
        bound_survival = real(dot(bound[:, end], doublon * bound[:, end]))
        @test bound_survival > free_survival + 0.2
        @test bound_survival > 0.5
    end

    @testset "#1848 dynamical spin structure factor" begin
        L = 10
        basis = SpinBasis1D(L; nup=L ÷ 2, pauli=false)
        H = Hamiltonian(
            basis,
            _lkm_spin_exchange_terms(L);
            static_fmt=:csc,
        )
        ground_energy, ground = _lkm_ground_state(H)
        momentum = 2pi / L
        spin_q = Matrix(Hamiltonian(
            basis,
            [OperatorTerm(
                "z",
                [(cis(momentum * (site - 1)), site) for site in 1:L],
            )];
            static_fmt=:csc,
            check_herm=false,
        ))
        frequencies = collect(0.0:0.04:4.0)
        broadening = 0.08
        spectrum = spectral_function(
            H,
            ground,
            spin_q,
            frequencies;
            reference_energy=ground_energy,
            broadening,
            method=:krylov,
            krylov_dim=180,
        )
        sum_rule = norm(spin_q * ground)^2
        integrated = _lkm_trapezoid(
            spectrum,
            frequencies[2] - frequencies[1],
        )
        @test minimum(spectrum) > -2e-10
        @test integrated ≈ sum_rule rtol=0.12
        @test frequencies[argmax(spectrum)] > 0
    end

    @testset "sector-changing particle-addition spectrum" begin
        L = 8
        function interacting_chain(basis)
            return Hamiltonian(
                basis,
                vcat(
                    _lkm_fermion_hopping_terms(L),
                    [OperatorTerm(
                        "nn",
                        [(2.0, site, site + 1) for site in 1:(L - 1)],
                    )],
                );
                static_fmt=:csc,
            )
        end
        source_basis = SpinlessFermionBasis1D(L; Nf=3)
        target_basis = SpinlessFermionBasis1D(L; Nf=4)
        source_hamiltonian = interacting_chain(source_basis)
        target_hamiltonian = interacting_chain(target_basis)
        source_energy, source_state = _lkm_ground_state(source_hamiltonian)
        transition = op_shift_sector(
            target_basis,
            source_basis,
            [("+", [L ÷ 2], 1.0)],
            Matrix{ComplexF64}(I, length(source_basis), length(source_basis)),
        )
        frequencies = collect(-2.0:0.04:7.0)
        spectrum = spectral_function(
            target_hamiltonian,
            source_state,
            transition,
            frequencies;
            reference_energy=source_energy,
            broadening=0.1,
            method=:lehmann,
        )
        expected_weight = norm(transition * source_state)^2
        integrated = _lkm_trapezoid(
            spectrum,
            frequencies[2] - frequencies[1],
        )
        @test integrated ≈ expected_weight rtol=0.08
        @test maximum(spectrum) > 0.1
    end

    @testset "small-system OTOC growth" begin
        L = 7
        basis = SpinBasis1D(L)
        H = Matrix(_lkm_tfim(basis, 1.05))
        decomposition = eigen(Hermitian(H))
        ground = decomposition.vectors[:, 1]
        W = operator_matrix(basis, "z", [(1.0, 1)])
        V = operator_matrix(basis, "z", [(1.0, L)])
        values = Float64[]
        for time in (0.0, 0.3, 0.8, 1.4)
            unitary = decomposition.vectors *
                Diagonal(exp.(-im .* decomposition.values .* time)) *
                decomposition.vectors'
            Wt = unitary' * W * unitary
            push!(values, norm((Wt * V - V * Wt) * ground)^2)
        end
        @test first(values) < 1e-20
        @test last(values) > values[2]
        @test maximum(values) > 1e-4
    end
end
