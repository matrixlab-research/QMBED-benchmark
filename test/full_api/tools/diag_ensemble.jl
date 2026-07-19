@testset "diagonal ensemble oracle and protocol" begin
    V2 = [1.0 1.0; 1.0 -1.0] ./ sqrt(2)
    E2 = [-1.0, 1.0]
    psi = [1.0, 0.0]
    observable = Diagonal([1.0, -1.0])

    result = diag_ensemble(
        2,
        psi,
        E2,
        V2;
        density=false,
        rho_d=true,
        Obs=observable,
        delta_t_Obs=true,
        delta_q_Obs=true,
        Sd_Renyi=true,
    )
    @test result["rho_d"] ≈ [0.5, 0.5] atol=3e-16
    @test result["Obs_pure"] ≈ -2.2371143170757376e-17 atol=3e-16
    @test result["delta_t_Obs_pure"] ≈ 0.7071067811865476 atol=4e-16
    @test result["delta_q_Obs_pure"] ≈ 0.7071067811865476 atol=4e-16
    @test result["Sd_pure"] ≈ 0.6931471805599454 atol=3e-16

    density_result = diag_ensemble(
        2,
        psi,
        E2,
        V2;
        density=true,
        rho_d=true,
        Obs=observable,
        Sd_Renyi=true,
    )
    @test density_result["Obs_pure"] ≈ -1.1185571585378688e-17 atol=3e-16
    @test density_result["Sd_pure"] ≈ 0.3465735902799727 atol=3e-16

    dm_result = diag_ensemble(
        2,
        Diagonal([0.75, 0.25]),
        E2,
        Matrix(I, 2, 2);
        density=false,
        rho_d=true,
        Obs=observable,
        Sd_Renyi=true,
    )
    @test dm_result["rho_d"] ≈ [0.75, 0.25] atol=3e-16
    @test dm_result["Obs_DM"] ≈ 0.5 atol=3e-16
    @test dm_result["Sd_DM"] ≈ 0.5623351446188083 atol=3e-16
end
