@testset "complete-eigensystem evolution oracle" begin
    A = [
        1.0 0.2 0.0 0.0
        0.2 2.0 0.3 0.0
        0.0 0.3 3.0 0.4
        0.0 0.0 0.4 4.0
    ]
    E, V = eigen(Symmetric(A))
    psi = ComplexF64[1, 0, 0, 0]
    psi_t = ed_state_vs_time(psi, E, V, [0.0, 0.25, 1.0])
    expected_last = ComplexF64[
        0.5356826557616654 - 0.8227743574000482im,
        -0.18743595940496593 - 0.014632856019069006im,
        0.011020113215605741 + 0.024628676702582557im,
        0.002085443746596086 - 0.002779582829544703im,
    ]
    @test psi_t[:, 1] ≈ psi atol=5e-16
    @test psi_t[:, end] ≈ expected_last atol=8e-16
    @test all(isapprox(norm(column), 1.0; atol=5e-15) for column in eachcol(psi_t))

    rho = psi * psi'
    rho_t = ed_state_vs_time(rho, E, V, [0.0, 0.25, 1.0])
    @test rho_t[:, :, 1] ≈ rho atol=8e-16
    @test all(isapprox(tr(rho_t[:, :, index]), 1.0; atol=5e-15) for index in axes(rho_t, 3))
    @test all(isapprox(rho_t[:, :, index], rho_t[:, :, index]'; atol=5e-15) for index in axes(rho_t, 3))

    ode = (time, state, frequency) -> -im * frequency .* state
    ode_times = [0.0, 0.2, 1.0]
    ode_states = evolve(
        ComplexF64[1],
        0.0,
        ode_times,
        ode;
        f_params=(2.0,),
        max_step=0.002,
    )
    @test vec(ode_states) ≈ ComplexF64[
        1.0 + 0.0im,
        0.9210609940028851 - 0.3894183423086505im,
        -0.4161468365471424 - 0.9092974268256817im,
    ] atol=1e-11
    @test collect(
        evolve(
            ComplexF64[1],
            0.0,
            ode_times,
            ode;
            f_params=(2.0,),
            max_step=0.002,
            iterate=true,
        ),
    ) ≈ collect(eachcol(ode_states)) atol=2e-12
end
