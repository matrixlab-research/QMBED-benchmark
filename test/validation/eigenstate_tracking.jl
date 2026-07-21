@testset "LKM validation: eigenstate continuation" begin
    @testset "#30690 TFIM fidelity susceptibility" begin
        L = 10
        basis = SpinBasis1D(L)
        fields = collect(0.75:0.1:1.25)
        spaces = Matrix{ComplexF64}[]
        gaps = Float64[]
        for field in fields
            values, vectors = eigsh(
                _lkm_tfim(basis, field);
                k=2,
                which=:SA,
                tol=1e-10,
            )
            order = sortperm(values)
            push!(spaces, vectors[:, order[1:1]])
            push!(gaps, values[order[2]] - values[order[1]])
        end
        tracked = track_eigenspaces(spaces)
        @test all(0 .< tracked.fidelities .<= 1)
        @test argmin(tracked.fidelities) in 2:5
        @test minimum(gaps) < maximum(gaps)
        @test subspace_fidelity(
            tracked.spaces[3],
            tracked.spaces[4],
        ) ≈ tracked.fidelities[3] atol=3e-12
    end

    @testset "#21745 Kitaev-chain parity-sector gap" begin
        L = 8
        function parity_hamiltonian(parity, chemical_potential)
            sectors = parity === :even ? collect(0:2:L) : collect(1:2:L)
            basis = SpinlessFermionBasis1D(L; Nf=sectors)
            terms = vcat(
                _lkm_fermion_hopping_terms(L),
                [
                    OperatorTerm(
                        "n",
                        [(-chemical_potential, site) for site in 1:L],
                    ),
                    OperatorTerm(
                        "++",
                        [(1.0, site, site + 1) for site in 1:(L - 1)],
                    ),
                    OperatorTerm(
                        "--",
                        [(-1.0, site, site + 1) for site in 1:(L - 1)],
                    ),
                ],
            )
            return Hamiltonian(
                basis,
                terms;
                static_fmt=:csc,
                check_pcon=false,
            )
        end

        function parity_splitting(chemical_potential)
            even_energy = first(eigsh(
                parity_hamiltonian(:even, chemical_potential);
                k=1,
                which=:SA,
                return_eigenvectors=false,
                tol=1e-11,
            ))
            odd_energy = first(eigsh(
                parity_hamiltonian(:odd, chemical_potential);
                k=1,
                which=:SA,
                return_eigenvectors=false,
                tol=1e-11,
            ))
            return abs(even_energy - odd_energy)
        end

        topological = parity_splitting(0.0)
        trivial = parity_splitting(4.0)
        @test topological < 1e-8
        @test trivial > 0.5
    end
end
