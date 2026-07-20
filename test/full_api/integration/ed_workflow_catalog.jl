using Statistics

function _edwf_spin_exchange_terms(
    L;
    J1=1.0,
    J2=0.0,
    delta=1.0,
    periodic=false,
)
    nearest = [(site, site + 1) for site in 1:(L - 1)]
    periodic && push!(nearest, (L, 1))
    next_nearest = [(site, site + 2) for site in 1:(L - 2)]
    periodic && append!(next_nearest, [(L - 1, 1), (L, 2)])
    xy = [
        (J1 / 2, left, right) for (left, right) in nearest
    ]
    append!(
        xy,
        [(J2 / 2, left, right) for (left, right) in next_nearest],
    )
    zz = [(delta * J1, left, right) for (left, right) in nearest]
    append!(
        zz,
        [(delta * J2, left, right) for (left, right) in next_nearest],
    )
    return [
        OperatorTerm("+-", xy),
        OperatorTerm("-+", xy),
        OperatorTerm("zz", zz),
    ]
end

function _edwf_fermion_hopping_terms(bonds)
    return [
        OperatorTerm(
            "+-",
            [(-amplitude, left, right) for (amplitude, left, right) in bonds],
        ),
        OperatorTerm(
            "-+",
            [(amplitude, left, right) for (amplitude, left, right) in bonds],
        ),
    ]
end

function _edwf_sorted_spectrum(H)
    return sort(real.(eigvals(Hermitian(Matrix(H)))))
end

function _edwf_flux_ring(L, particles, flux; hopping=1.0, interaction=0.0)
    basis = SpinlessFermionBasis1D(L; Nf=particles)
    phase = cis(flux / L)
    plus = Tuple{ComplexF64,Int,Int}[]
    minus = Tuple{ComplexF64,Int,Int}[]
    for site in 1:L
        neighbor = mod1(site + 1, L)
        push!(plus, (-hopping * phase, site, neighbor))
        push!(minus, (hopping * conj(phase), site, neighbor))
    end
    terms = OperatorTerm[
        OperatorTerm("+-", plus),
        OperatorTerm("-+", minus),
    ]
    if !iszero(interaction)
        push!(
            terms,
            OperatorTerm(
                "nn",
                [
                    (interaction, site, mod1(site + 1, L))
                    for site in 1:L
                ],
            ),
        )
    end
    return basis, Hamiltonian(
        basis,
        terms;
        static_fmt=:csc,
        check_herm=false,
    )
end

@testset "workflow-derived exact-diagonalization user scenarios" begin
    @testset "01 #3831 Heisenberg dimer singlet and triplet" begin
        basis = SpinBasis1D(2; pauli=false)
        H = Hamiltonian(
            basis,
            _edwf_spin_exchange_terms(2);
            static_fmt=:csc,
        )
        @test _edwf_sorted_spectrum(H) ≈ [-0.75, 0.25, 0.25, 0.25] atol=2e-14
    end

    @testset "02 #3831 symmetry-resolved XXZ low-energy spectrum" begin
        L = 6
        terms = _edwf_spin_exchange_terms(L; delta=0.7, periodic=true)
        full = Hamiltonian(
            SpinBasis1D(L; pauli=false),
            terms;
            static_fmt=:csc,
        )
        sector_basis = SpinBasis1D(
            L;
            nup=L ÷ 2,
            pauli=false,
            kblock=0,
        )
        sector = Hamiltonian(sector_basis, terms; static_fmt=:csc)
        sector_values = _edwf_sorted_spectrum(sector)
        full_values = _edwf_sorted_spectrum(full)
        @test length(sector_basis) < length(full.basis)
        @test all(
            minimum(abs.(full_values .- value)) < 2e-12
            for value in sector_values
        )
        selected = first(eigsh(sector; k=2, which=:SA, tol=1e-11))
        @test sort(selected) ≈ sector_values[1:2] atol=2e-10
    end

    @testset "03 #23806 frustrated J1-J2 spin-gap kernel" begin
        L = 8
        singlet_basis = SpinBasis1D(L; nup=L ÷ 2, pauli=false)
        triplet_basis = SpinBasis1D(L; nup=L ÷ 2 + 1, pauli=false)
        terms = _edwf_spin_exchange_terms(
            L;
            J1=1.0,
            J2=0.5,
            periodic=true,
        )
        singlet = Hamiltonian(singlet_basis, terms; static_fmt=:csc)
        triplet = Hamiltonian(triplet_basis, terms; static_fmt=:csc)
        E0 = minimum(_edwf_sorted_spectrum(singlet))
        E1 = minimum(_edwf_sorted_spectrum(triplet))
        @test E1 - E0 > 0.1
        @test E1 - E0 < 1.0
    end

    @testset "04 #23806 transverse-field Ising phase-scan point" begin
        L = 7
        basis = SpinBasis1D(L)
        bonds = [(1.0, site, mod1(site + 1, L)) for site in 1:L]
        function ising(g)
            return Hamiltonian(
                basis,
                [
                    OperatorTerm("zz", [(-value, i, j) for (value, i, j) in bonds]),
                    OperatorTerm("x", [(-g, site) for site in 1:L]),
                ];
                static_fmt=:csc,
            )
        end
        ordered = _edwf_sorted_spectrum(ising(0.35))
        polarized = _edwf_sorted_spectrum(ising(2.0))
        @test ordered[2] - ordered[1] < polarized[2] - polarized[1]
        @test ishermitian(ising(1.0))
    end

    @testset "05 #11058 magnetic susceptibility of a Heisenberg dimer" begin
        beta = 1.3
        basis = SpinBasis1D(2; pauli=false)
        H = Hamiltonian(
            basis,
            _edwf_spin_exchange_terms(2);
            static_fmt=:dense,
        )
        Mz = Matrix(
            Hamiltonian(
                basis,
                [OperatorTerm("z", [(1.0, 1), (1.0, 2)])],
            ),
        )
        spectrum = eigen(Hermitian(Matrix(H)))
        weights = exp.(-beta .* spectrum.values)
        thermal_m2 = sum(
            weights[index] *
            real(
                dot(
                    spectrum.vectors[:, index],
                    Mz * Mz * spectrum.vectors[:, index],
                ),
            )
            for index in eachindex(weights)
        ) / sum(weights)
        expected = 2beta * exp(-beta / 4) /
            (exp(3beta / 4) + 3exp(-beta / 4))
        @test beta * thermal_m2 ≈ expected atol=3e-14
    end

    @testset "06 #12290 variational energy certified by ED" begin
        L = 6
        basis = SpinBasis1D(L; nup=L ÷ 2, pauli=false)
        H = Hamiltonian(
            basis,
            _edwf_spin_exchange_terms(L; J2=0.35);
            static_fmt=:csc,
        )
        neel = sum(UInt64(1) << (site - 1) for site in 1:2:L)
        trial = zeros(ComplexF64, length(basis))
        trial[state_index(basis, neel)] = 1
        variational_energy = real(expt_value(H, trial))
        exact_energy = minimum(_edwf_sorted_spectrum(H))
        @test variational_energy >= exact_energy - 2e-13
        @test variational_energy - exact_energy > 0.2
    end

    @testset "07 #1069 step-driven Floquet quasienergies" begin
        basis = SpinBasis1D(4)
        Hzz = Hamiltonian(
            basis,
            [
                OperatorTerm(
                    "zz",
                    [(0.7, site, mod1(site + 1, 4)) for site in 1:4],
                ),
            ];
            static_fmt=:csc,
        )
        Hx = Hamiltonian(
            basis,
            [OperatorTerm("x", [(0.45, site) for site in 1:4])];
            static_fmt=:csc,
        )
        floquet = Floquet(
            Dict(:H_list => [Hzz, Hx], :dt_list => [0.23, 0.31]);
            UF=true,
            thetaF=true,
        )
        identity_matrix = Matrix{ComplexF64}(I, length(basis), length(basis))
        @test floquet.UF' * floquet.UF ≈ identity_matrix atol=2e-12
        @test all(abs.(abs.(floquet.thetaF) .- 1) .< 2e-12)
    end

    @testset "08 #1069 high-frequency effective Hamiltonian limit" begin
        X = ComplexF64[0 1; 1 0]
        Z = ComplexF64[1 0; 0 -1]
        dt = 0.005
        floquet = Floquet(
            Dict(:H_list => [X, Z], :dt_list => [dt, dt]);
            HF=true,
            UF=true,
        )
        average_H = (X + Z) / 2
        @test norm(floquet.HF - average_H) < 0.01
        @test floquet.UF ≈ exp(-im * dt .* Z) * exp(-im * dt .* X) atol=3e-15
    end

    @testset "09 #19021 Anderson-localization inverse participation ratio" begin
        L = 10
        basis = SpinlessFermionBasis1D(L; Nf=1)
        bonds = [(1.0, site, site + 1) for site in 1:(L - 1)]
        hopping = _edwf_fermion_hopping_terms(bonds)
        clean = Hamiltonian(basis, hopping; static_fmt=:csc)
        disordered = Hamiltonian(
            basis,
            vcat(
                hopping,
                [
                    OperatorTerm(
                        "n",
                        [(20.0 * (-1)^site + 0.31site, site) for site in 1:L],
                    ),
                ],
            );
            static_fmt=:csc,
        )
        clean_vectors = eigen(Hermitian(Matrix(clean))).vectors
        localized_vectors = eigen(Hermitian(Matrix(disordered))).vectors
        clean_ipr = mean(sum(abs2.(clean_vectors) .^ 2; dims=1))
        localized_ipr = mean(sum(abs2.(localized_vectors) .^ 2; dims=1))
        @test localized_ipr > 2clean_ipr
    end

    @testset "10 #19021 MBL mid-spectrum level statistics" begin
        L = 8
        basis = SpinBasis1D(L; nup=L ÷ 2, pauli=false)
        fields = [2.1, -1.7, 0.3, 3.2, -2.6, 0.8, 1.4, -3.1]
        terms = vcat(
            _edwf_spin_exchange_terms(L; periodic=true),
            [OperatorTerm("z", [(fields[site], site) for site in 1:L])],
        )
        H = Hamiltonian(basis, terms; static_fmt=:csc)
        levels = _edwf_sorted_spectrum(H)
        spacings = diff(levels)
        ratios = [
            min(spacings[index], spacings[index + 1]) /
            max(spacings[index], spacings[index + 1])
            for index in 1:(length(spacings) - 1)
            if max(spacings[index], spacings[index + 1]) > 1e-10
        ]
        @test 0 < mean(ratios) < 1
        selected = first(
            eigsh(H; k=4, sigma=mean(levels), which=:LM, tol=1e-10),
        )
        @test all(
            minimum(abs.(levels .- value)) < 2e-9 for value in selected
        )
    end

    @testset "11 #19000 SSH topological edge modes" begin
        L = 8
        basis = SpinlessFermionBasis1D(L; Nf=1)
        bonds = [
            (isodd(site) ? 0.25 : 1.0, site, site + 1)
            for site in 1:(L - 1)
        ]
        H = Hamiltonian(
            basis,
            _edwf_fermion_hopping_terms(bonds);
            static_fmt=:csc,
        )
        decomposition = eigen(Hermitian(Matrix(H)))
        edge_indices = sortperm(abs.(decomposition.values))[1:2]
        @test maximum(abs.(decomposition.values[edge_indices])) < 0.01
        edge_weight = mean(
            sum(
                abs2(decomposition.vectors[row, index]) *
                (basis.occupations[row, 1] + basis.occupations[row, L])
                for row in 1:length(basis)
            )
            for index in edge_indices
        )
        @test edge_weight > 0.8
    end

    @testset "12 #19000 interacting extended SSH many-body gap" begin
        L = 6
        basis = SpinlessFermionBasis1D(L; Nf=3)
        bonds = [
            (isodd(site) ? 0.45 : 1.0, site, site + 1)
            for site in 1:(L - 1)
        ]
        kinetic = _edwf_fermion_hopping_terms(bonds)
        free = Hamiltonian(basis, kinetic; static_fmt=:csc)
        interacting = Hamiltonian(
            basis,
            vcat(
                kinetic,
                [
                    OperatorTerm(
                        "nn",
                        [(2.0, site, site + 1) for site in 1:(L - 1)],
                    ),
                ],
            );
            static_fmt=:csc,
        )
        free_levels = _edwf_sorted_spectrum(free)
        interacting_levels = _edwf_sorted_spectrum(interacting)
        @test interacting_levels[1] != free_levels[1]
        @test interacting_levels[2] - interacting_levels[1] > 1e-3
    end

    @testset "13 #16282 chiral spectral pairing in the SSH model" begin
        L = 7
        basis = SpinlessFermionBasis1D(L; Nf=1)
        bonds = [
            (isodd(site) ? 0.6 : 1.2, site, site + 1)
            for site in 1:(L - 1)
        ]
        H = Hamiltonian(
            basis,
            _edwf_fermion_hopping_terms(bonds);
            static_fmt=:csc,
        )
        levels = _edwf_sorted_spectrum(H)
        @test levels ≈ -reverse(levels) atol=3e-14
        @test minimum(abs.(levels)) < 2e-14
    end

    @testset "14 #16282 flux periodicity as a topological diagnostic" begin
        _, zero_flux = _edwf_flux_ring(7, 1, 0.0)
        _, flux_quantum = _edwf_flux_ring(7, 1, 2π)
        @test _edwf_sorted_spectrum(zero_flux) ≈
            _edwf_sorted_spectrum(flux_quantum) atol=3e-13
    end

    @testset "15 #2880 few-electron quantum-dot shell filling" begin
        basis = SpinlessFermionBasis1D(4; Nf=2)
        orbital_energies = [0.0, 0.7, 1.9, 3.1]
        interactions = [
            (0.2 + 0.1abs(left - right), left, right)
            for left in 1:4 for right in (left + 1):4
        ]
        H = Hamiltonian(
            basis,
            [
                OperatorTerm(
                    "n",
                    [(orbital_energies[site], site) for site in 1:4],
                ),
                OperatorTerm("nn", interactions),
            ];
            static_fmt=:csc,
        )
        @test minimum(_edwf_sorted_spectrum(H)) ≈ 1.0 atol=2e-14
        @test count(value -> !iszero(value), nonzeros(H.data)) <= length(basis)
    end

    @testset "16 #2880 two-electron Hubbard-dot singlet energy" begin
        hopping = [
            OperatorTerm("+-|", [(-1.0, 1, 2)]),
            OperatorTerm("-+|", [(1.0, 1, 2)]),
            OperatorTerm("|+-", [(-1.0, 1, 2)]),
            OperatorTerm("|-+", [(1.0, 1, 2)]),
        ]
        U = 3.0
        basis = SpinfulFermionBasis1D(2; Nf=(1, 1))
        H = Hamiltonian(
            basis,
            vcat(
                hopping,
                [OperatorTerm("n|n", [(U, 1, 1), (U, 2, 2)])],
            );
            static_fmt=:csc,
        )
        expected = (U - sqrt(U^2 + 16)) / 2
        @test minimum(_edwf_sorted_spectrum(H)) ≈ expected atol=3e-14
    end

    @testset "17 #3127 Cu-O cluster orbital-hole fractions" begin
        basis = SpinfulFermionBasis1D(3; Nf=(1, 1))
        bonds = [(1, 2), (1, 3)]
        terms = [
            OperatorTerm("+-|", [(-0.4, pair...) for pair in bonds]),
            OperatorTerm("-+|", [(0.4, pair...) for pair in bonds]),
            OperatorTerm("|+-", [(-0.4, pair...) for pair in bonds]),
            OperatorTerm("|-+", [(0.4, pair...) for pair in bonds]),
            OperatorTerm("n|", [(4.0, 2), (4.0, 3)]),
            OperatorTerm("|n", [(4.0, 2), (4.0, 3)]),
            OperatorTerm("n|n", [(0.5, site, site) for site in 1:3]),
        ]
        H = Hamiltonian(basis, terms; static_fmt=:csc)
        ground = eigen(Hermitian(Matrix(H))).vectors[:, 1]
        site_occupations = [
            real(
                expt_value(
                    Hamiltonian(
                        basis,
                        [
                            OperatorTerm("n|", [(1.0, site)]),
                            OperatorTerm("|n", [(1.0, site)]),
                        ],
                    ),
                    ground,
                ),
            )
            for site in 1:3
        ]
        @test sum(site_occupations) ≈ 2 atol=2e-13
        @test site_occupations[1] > 1.8
    end

    @testset "18 #3127 Hubbard-cluster charge-gap extraction" begin
        function cluster_energy(particles)
            basis = SpinlessFermionBasis1D(5; Nf=particles)
            kinetic = _edwf_fermion_hopping_terms(
                [(0.35, site, site + 1) for site in 1:4],
            )
            H = Hamiltonian(
                basis,
                vcat(
                    kinetic,
                    [
                        OperatorTerm(
                            "nn",
                            [(3.0, site, site + 1) for site in 1:4],
                        ),
                    ],
                );
                static_fmt=:csc,
            )
            return minimum(_edwf_sorted_spectrum(H))
        end
        charge_gap = cluster_energy(3) + cluster_energy(1) - 2cluster_energy(2)
        @test charge_gap > 0
    end

    @testset "19 #12454 FQHE-style interacting flux spectral flow" begin
        _, at_zero = _edwf_flux_ring(6, 2, 0.0; interaction=4.0)
        _, at_quantum = _edwf_flux_ring(6, 2, 2π; interaction=4.0)
        spectrum_zero = _edwf_sorted_spectrum(at_zero)
        spectrum_quantum = _edwf_sorted_spectrum(at_quantum)
        @test spectrum_zero ≈ spectrum_quantum atol=3e-13
        @test spectrum_zero[2] - spectrum_zero[1] < 3e-13
        @test spectrum_zero[3] - spectrum_zero[2] > 1e-4
    end

    @testset "20 #14428 collective-emission bright and dark modes" begin
        emitters = 4
        decay_rate = 0.6
        basis = SpinBasis1D(emitters; nup=1, pauli=false)
        collective_decay = Hamiltonian(
            basis,
            [
                OperatorTerm(
                    "+-",
                    [
                        (-im * decay_rate / 2, left, right)
                        for left in 1:emitters for right in 1:emitters
                    ],
                ),
            ];
            static_fmt=:csc,
            check_herm=false,
        )
        values = eigvals(Matrix(collective_decay))
        @test count(value -> abs(value) < 2e-13, values) == emitters - 1
        @test minimum(abs.(values .+ im * emitters * decay_rate / 2)) < 2e-13
    end

    @testset "21 #26341 FCI-style many-body entanglement diagnostic" begin
        basis, H = _edwf_flux_ring(6, 3, π / 3; interaction=5.0)
        ground = eigen(Hermitian(Matrix(H))).vectors[:, 1]
        entropy = ent_entropy(
            basis,
            ground;
            sub_sys_A=[1, 2, 3],
            return_rdm=:A,
        )
        @test entropy["Sent_A"] > 0
        @test tr(entropy["rdm_A"]) ≈ 1 atol=3e-13
    end
end
