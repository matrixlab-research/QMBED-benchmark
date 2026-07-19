@testset "held-out sparse workflows from published spin-chain studies" begin
    @testset "MBL shift-invert exact diagonalization" begin
        # Pal & Huse, Phys. Rev. B 82, 174411 (2010), arXiv:1003.2613.
        L = 8
        basis = SpinBasis1D(L; nup=L ÷ 2, pauli=false)
        fields = [1.9, -2.6, 0.3, 1.1, -0.8, 2.2, -1.5, 0.7]
        bonds = [(1.0, site, mod1(site + 1, L)) for site in 1:L]
        H = Hamiltonian(
            basis,
            [
                OperatorTerm("+-", [(0.5, i, j) for (_, i, j) in bonds]),
                OperatorTerm("-+", [(0.5, i, j) for (_, i, j) in bonds]),
                OperatorTerm("zz", bonds),
                OperatorTerm("z", [(fields[site], site) for site in 1:L]),
            ];
            static_fmt=:csc,
        )
        seed = normalize(collect(1.0:length(basis)))
        selected, vectors = eigsh(
            H;
            k=3,
            sigma=0.0,
            which=:LM,
            v0=seed,
            tol=1e-12,
            maxiter=1_500,
        )
        exact = eigvals(Hermitian(Matrix(H)))
        expected = sort(exact; by=abs)[1:3]
        @test H.data isa SparseMatrixCSC
        @test nnz(H.data) < length(basis)^2 ÷ 5
        @test sort(selected; by=abs) ≈ expected atol=2e-10
        @test norm(Matrix(H) * vectors - vectors * Diagonal(selected)) < 2e-9
    end

    @testset "Lanczos quench and imbalance" begin
        # Sparse Krylov time evolution is the computational kernel used by the
        # QuSpin paper, SciPost Phys. 2, 003 (2017), arXiv:1610.03042.
        L = 8
        basis = SpinBasis1D(L; nup=L ÷ 2, pauli=false)
        bonds = [(1.0, site, mod1(site + 1, L)) for site in 1:L]
        H = Hamiltonian(
            basis,
            [
                OperatorTerm("+-", [(0.5, i, j) for (_, i, j) in bonds]),
                OperatorTerm("-+", [(0.5, i, j) for (_, i, j) in bonds]),
                OperatorTerm("zz", bonds),
            ];
            static_fmt=:csc,
        )
        neel = sum(UInt64(1) << (site - 1) for site in 1:2:L)
        psi0 = zeros(ComplexF64, length(basis))
        psi0[state_index(basis, neel)] = 1
        E, V, Q_T = lanczos_full(H, psi0, 30; full_ortho=true)
        evolved = expm_lanczos(E, V, Q_T; a=-0.31im)
        exact = exp((-0.31im) .* Matrix(H)) * psi0
        @test evolved ≈ exact atol=2e-12
        @test norm(evolved) ≈ 1.0 atol=2e-13
    end

    @testset "Floquet heating step drive" begin
        # Periodically driven transverse-field Ising workflow from
        # arXiv:1610.03042, evaluated here at a distinct held-out size.
        L = 5
        J, g, h, Omega = 1.0, 0.809, 0.9045, 4.5
        period = 2π / Omega
        basis = SpinBasis1D(L)
        bonds = [(J, site, mod1(site + 1, L)) for site in 1:L]
        z_field = [(h, site) for site in 1:L]
        x_field = [(g, site) for site in 1:L]
        drive(time, frequency) = sign(cos(frequency * time))
        H = Hamiltonian(
            Any[
                Any["zz", bonds],
                Any["z", z_field],
                Any["x", x_field],
            ],
            Any[
                Any["zz", bonds, drive, (Omega,)],
                Any["z", z_field, drive, (Omega,)],
                Any[
                    "x",
                    [(-value, site) for (value, site) in x_field],
                    drive,
                    (Omega,),
                ],
            ];
            basis,
            dtype=Float64,
            static_fmt=:csc,
            dynamic_fmt=:csc,
        )
        sample_times = [eps(Float64), period / 2 + eps(Float64)]
        matrices = [tocsc(H; time) for time in sample_times]
        spectrum = Floquet(
            Dict(
                :H => H,
                :t_list => sample_times,
                :dt_list => [period / 2, period / 2],
            );
            UF=true,
        )
        expected = exp((-im * period / 2) .* Matrix(matrices[2])) *
            exp((-im * period / 2) .* Matrix(matrices[1]))
        identity_matrix = Matrix{ComplexF64}(I, length(basis), length(basis))
        @test all(matrix isa SparseMatrixCSC for matrix in matrices)
        @test spectrum.UF ≈ expected atol=2e-13
        @test spectrum.UF' * spectrum.UF ≈ identity_matrix atol=4e-13
    end
end
