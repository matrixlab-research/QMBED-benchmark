@testset "LKM validation: constrained and truncated dynamics" begin
    @testset "PXP many-body-scar revival" begin
        L = 14
        constrained = constraint_states(
            L;
            prefix_allowed=(occupations, site) ->
                site == 1 || occupations[site - 1] + occupations[site] <= 1,
            state_allowed=occupations -> occupations[1] + occupations[end] <= 1,
        )
        basis = UserBasis(
            UInt64,
            L,
            Dict(
                'x' => ComplexF64[0 1; 1 0],
                'n' => ComplexF64[0 0; 0 1],
            );
            states=constrained,
            allowed_ops=('x', 'n'),
        )
        @test length(basis) == 843
        H = Hamiltonian(
            basis,
            [OperatorTerm("x", [(1.0, site) for site in 1:L])];
            static_fmt=:csc,
        )
        z2_state = sum(UInt64(1) << (site - 1) for site in 1:2:L)
        initial = zeros(ComplexF64, length(basis))
        initial[state_index(basis, z2_state)] = 1
        times = collect(0.0:0.2:8.0)
        evolved = evolve(H, initial, 0.0, times; tol=1e-11, krylov_dim=60)
        fidelity = [abs2(dot(initial, state)) for state in eachcol(evolved)]
        @test all(isapprox(norm(state), 1; atol=2e-10) for state in eachcol(evolved))
        @test maximum(fidelity[times .> 3.0]) > 0.2
        @test minimum(fidelity) < 0.05
    end

    @testset "Bose-Hubbard Mott quench" begin
        L = 7
        basis = BosonBasis1D(L; Nb=L, sps=3)
        hopping_plus = OperatorTerm(
            "+-",
            [(-1.0, site, site + 1) for site in 1:(L - 1)],
        )
        hopping_minus = OperatorTerm(
            "-+",
            [(-1.0, site, site + 1) for site in 1:(L - 1)],
        )
        interaction = OperatorTerm(
            "nn",
            [(4.0, site, site) for site in 1:L],
        )
        chemical = OperatorTerm("n", [(-4.0, site) for site in 1:L])
        H = Hamiltonian(
            basis,
            [hopping_plus, hopping_minus, interaction, chemical];
            static_fmt=:csc,
        )
        mott = sum(UInt64(3)^(site - 1) for site in 1:L)
        initial = zeros(ComplexF64, length(basis))
        initial[state_index(basis, mott)] = 1
        number_operator = Hamiltonian(
            basis,
            [OperatorTerm("n", [(1.0, site) for site in 1:L])];
            static_fmt=:csc,
        )
        times = collect(0.0:0.1:2.0)
        evolved = evolve(H, initial, 0.0, times; tol=1e-11, krylov_dim=70)
        @test all(isapprox(norm(state), 1; atol=3e-10) for state in eachcol(evolved))
        @test all(
            isapprox(real(expt_value(number_operator, state)), L; atol=2e-10)
            for state in eachcol(evolved)
        )
        @test abs2(dot(initial, evolved[:, end])) < 0.95
    end

    @testset "Loschmidt return after an Ising quench" begin
        basis = SpinBasis1D(10)
        _, initial = _lkm_ground_state(_lkm_tfim(basis, 0.25))
        final_hamiltonian = _lkm_tfim(basis, 1.8)
        times = collect(0.0:0.1:3.0)
        evolved = evolve(
            final_hamiltonian,
            initial,
            0.0,
            times;
            tol=1e-11,
            krylov_dim=70,
        )
        echo = [abs2(dot(initial, state)) for state in eachcol(evolved)]
        return_rate = -log.(max.(echo, eps())) ./ basis.L
        @test first(echo) ≈ 1 atol=2e-12
        @test minimum(echo) < 0.2
        @test maximum(return_rate) > 0.1
    end

    @testset "finite-size Dicke photon-cutoff convergence" begin
        function dicke_ground(cutoff)
            spins = 3
            basis = PhotonBasis(SpinBasis1D, spins; Nph=cutoff, pauli=false)
            coupling = 0.35 / sqrt(spins)
            terms = OperatorTerm[
                OperatorTerm("|n", [(1.0, 1)]),
                OperatorTerm("z|", [(1.0, site) for site in 1:spins]),
                OperatorTerm(
                    "x|+",
                    [(coupling, site, 1) for site in 1:spins],
                ),
                OperatorTerm(
                    "x|-",
                    [(coupling, site, 1) for site in 1:spins],
                ),
            ]
            H = Hamiltonian(basis, terms; static_fmt=:csc)
            energy, state = _lkm_ground_state(H)
            photons = Hamiltonian(
                basis,
                [OperatorTerm("|n", [(1.0, 1)])];
                static_fmt=:csc,
            )
            return energy, real(expt_value(photons, state))
        end

        energy4, photons4 = dicke_ground(4)
        energy6, photons6 = dicke_ground(6)
        @test abs(energy6 - energy4) < 2e-4
        @test abs(photons6 - photons4) < 2e-3
        @test photons6 > 0
    end
end
