@testset "LKM validation: matrix-free open-system dynamics" begin
    @testset "collective-emission bright-state decay" begin
        emitters = 4
        decay = 0.3
        basis = SpinBasis1D(emitters; pauli=false)
        lowering = sum(
            operator_matrix(
                basis,
                "-",
                [(sqrt(decay), site)];
                sparse=true,
            )
            for site in 1:emitters
        )
        generator = LindbladGenerator(
            spzeros(ComplexF64, length(basis), length(basis)),
            [lowering],
        )
        bright = zeros(ComplexF64, length(basis))
        for site in 1:emitters
            bright[state_index(basis, UInt64(1) << (site - 1))] =
                inv(sqrt(emitters))
        end
        initial = bright * bright'
        excitation = sum(
            operator_matrix(basis, "z", [(1.0, site)])
            for site in 1:emitters
        ) + (emitters / 2) .* Matrix{Float64}(I, length(basis), length(basis))
        times = [0.0, 0.25, 0.75, 1.25]
        evolved = evolve(
            generator,
            initial,
            0.0,
            times;
            max_step=0.02,
            rtol=2e-9,
            atol=1e-11,
        )
        for (index, time) in pairs(times)
            rho = @view evolved[:, :, index]
            population = real(tr(excitation * rho))
            @test population ≈ exp(-emitters * decay * time) atol=3e-7
            @test tr(rho) ≈ 1 atol=3e-9
            @test rho ≈ rho' atol=3e-9
        end
    end
end
