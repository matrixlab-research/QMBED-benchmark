@testset "Second-audit algorithm regressions" begin
    @testset "continuous dynamic Floquet preserves time ordering" begin
        basis = SpinBasis1D(1; pauli=false)
        X = ComplexF64[0 0.5; 0.5 0]
        Z = ComplexF64[-0.5 0; 0 0.5]
        drive = (time, frequency) -> cos(frequency * time)
        H = Hamiltonian(
            Any[Z],
            Any[Any[X, drive, (1.0,)]];
            basis,
            dtype=ComplexF64,
        )
        period = 2π
        floquet = Floquet(Dict(:H => H, :T => period); UF=true)

        steps = 4_000
        dt = period / steps
        reference = Matrix{ComplexF64}(I, 2, 2)
        for step in 1:steps
            midpoint = (step - 0.5) * dt
            reference =
                exp((-im * dt) .* (Z + cos(midpoint) .* X)) *
                reference
        end
        @test floquet.UF ≈ reference atol=3e-7 rtol=3e-7
        @test floquet.UF' * floquet.UF ≈ I atol=2e-10
    end

    @testset "ExpOp acts without changing sparse semantics" begin
        dimension = 64
        diagonal = collect(range(-1.0, 1.0; length=dimension))
        off_diagonal = fill(0.15, dimension - 1)
        matrix = spdiagm(
            -1 => off_diagonal,
            0 => diagonal,
            1 => off_diagonal,
        )
        state = normalize(ComplexF64[
            sin(0.17index) + im * cos(0.11index)
            for index in 1:dimension
        ])
        exponential = ExpOp(matrix; a=-0.3im)
        @test exponential * state ≈
            exp((-0.3im) .* Matrix(matrix)) * state atol=2e-11 rtol=2e-11

        grid = ExpOp(
            matrix;
            a=-0.3im,
            start=0.0,
            stop=1.0,
            num=7,
        )
        values = grid * state
        @test size(values) == (dimension, 7)
        @test values[:, 1] ≈ state atol=2e-13
        @test values[:, end] ≈ exponential * state atol=2e-11 rtol=2e-11

        lazy_grid = ExpOp(
            matrix;
            a=-0.3im,
            start=0.0,
            stop=1.0,
            num=7,
            iterate=true,
        )
        @test reduce(hcat, collect(lazy_grid * state)) ≈ values atol=2e-11 rtol=2e-11
    end

    @testset "discrete Hamiltonian remains sparse and agrees with QLO" begin
        L = 6
        basis = SpinfulFermionBasis1D(L; Nf=(3, 3))
        bonds = [(site, site + 1) for site in 1:(L - 1)]
        terms = [
            OperatorTerm("+-|", [(-1.0, i, j) for (i, j) in bonds]),
            OperatorTerm("-+|", [(1.0, i, j) for (i, j) in bonds]),
            OperatorTerm("|+-", [(-1.0, i, j) for (i, j) in bonds]),
            OperatorTerm("|-+", [(1.0, i, j) for (i, j) in bonds]),
            OperatorTerm("n|n", [(2.0, site, site) for site in 1:L]),
        ]
        sparse_H = Hamiltonian(
            basis,
            terms;
            static_fmt=:csc,
            check_symm=false,
            check_herm=false,
            check_pcon=false,
        )
        qlo = QuantumLinearOperator(
            basis,
            terms;
            check_symm=false,
            check_herm=false,
            check_pcon=false,
        )
        @test sparse_H.data isa SparseMatrixCSC
        state = normalize(ComplexF64.(1:length(basis)))
        @test qlo * state ≈ sparse_H * state atol=3e-12 rtol=3e-12
        @test qlo.H * state ≈ Matrix(sparse_H)' * state atol=3e-12 rtol=3e-12
        @test right_apply(qlo, state) ≈
            vec(transpose(state) * Matrix(sparse_H)) atol=3e-12 rtol=3e-12
    end

    @testset "diagonal ensemble contraction matches dense square" begin
        dimension = 24
        raw = ComplexF64[
            sin(0.13row + 0.07column) +
            im * cos(0.05row - 0.11column)
            for row in 1:dimension, column in 1:dimension
        ]
        observable = Hermitian(raw + raw')
        state = normalize(ComplexF64[
            sin(0.19index) + im * cos(0.23index)
            for index in 1:dimension
        ])
        energies = collect(1.0:dimension)
        eigenvectors = Matrix{ComplexF64}(I, dimension, dimension)
        result = diag_ensemble(
            1,
            state,
            energies,
            eigenvectors;
            density=false,
            Obs=Matrix(observable),
            delta_t_Obs=true,
            delta_q_Obs=true,
        )
        @test isfinite(result["delta_t_Obs_pure"])
        @test isfinite(result["delta_q_Obs_pure"])
    end

    @testset "adaptive generic evolution and user particle enumeration" begin
        calls = Ref(0)
        derivative! = function (destination, time, state, rate)
            calls[] += 1
            @. destination = rate * state
            return destination
        end
        evolved = evolve(
            ComplexF64[1],
            0.0,
            [2.0],
            derivative!;
            f_params=(-0.7im,),
            max_step=0.5,
            rtol=1e-10,
            atol=1e-12,
        )
        @test only(evolved) ≈ exp(-1.4im) atol=2e-10
        @test calls[] > 0

        flip = (state, site) -> (
            xor(UInt64(state), UInt64(1) << (site - 1)),
            1.0,
        )
        next_state = function (state, counter, N, arguments)
            prefix = (state | (state - UInt64(1))) + UInt64(1)
            return prefix | (
                (
                    (
                        (prefix & (-prefix)) ÷
                        (state & (-state))
                    ) >> 1
                ) - UInt64(1)
            )
        end
        basis = UserBasis(
            UInt64,
            24,
            Dict('x' => flip);
            pcon_dict=Dict(
                :Np => 2,
                :next_state => next_state,
                :next_state_args => (),
                :get_Ns_pcon => (N, Np) -> binomial(N, Np),
                :get_s0_pcon =>
                    (N, Np) -> (UInt64(1) << Np) - UInt64(1),
            ),
        )
        @test length(basis) == binomial(24, 2)
        @test all(count_ones(state) == 2 for state in states(basis))
    end
end
