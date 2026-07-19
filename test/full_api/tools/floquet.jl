@testset "Floquet oracle and protocol" begin
    H = ComplexF64[1 0; 0 -1]
    floquet = Floquet(
        Dict(:H => H, :T => 0.5);
        HF=true,
        UF=true,
        thetaF=true,
        VF=true,
    )

    @test floquet.T == 0.5
    @test floquet.EF ≈ [-1.0, 1.0] atol=2e-16
    @test floquet.UF ≈ ComplexF64[
        0.8775825618903728-0.479425538604203im 0
        0 0.8775825618903728+0.479425538604203im
    ] atol=4e-16
    @test floquet.thetaF ≈ ComplexF64[
        0.8775825618903728+0.479425538604203im,
        0.8775825618903728-0.479425538604203im,
    ] atol=4e-16
    @test floquet.HF ≈ H atol=2e-15
    @test floquet.VF' * floquet.VF ≈ Matrix{ComplexF64}(I, 2, 2) atol=2e-16

    X = ComplexF64[0 1; 1 0]
    expected = exp((-0.3im) .* X) * exp((-0.2im) .* H)
    stepped = Floquet(
        Dict(:H_list => [H, X], :dt_list => [0.2, 0.3]);
        UF=true,
    )
    @test stepped.UF ≈ expected atol=4e-16
    @test stepped.T == 0.5
    @test stepped.HF === nothing
    @test stepped.thetaF === nothing
    @test stepped.VF === nothing
end
