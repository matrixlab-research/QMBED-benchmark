@testset "Lanczos oracle, reconstruction, and integration" begin
    A = [
        1.0 0.2 0.0 0.0
        0.2 2.0 0.3 0.0
        0.0 0.3 3.0 0.4
        0.0 0.0 0.4 4.0
    ]
    v0 = [1.0, 2.0, -1.0, 0.5]
    E, V, Q_T = lanczos_full(A, v0, 3)
    @test E ≈ [
        1.0981001507517356,
        1.9892033688829862,
        3.8667165876454432,
    ] atol=5e-15
    @test Q_T ≈ [
        0.4 0.8 -0.4 0.2
        -0.6350283337363477 -0.00460165459229314 -0.28530258472212694 0.7178581163976099
        0.6324241608687684 -0.3665783021968704 0.22191773407160648 0.6453003551931594
    ] atol=3e-15
    @test Q_T * Q_T' ≈ Matrix(I, 3, 3) atol=3e-15
    @test lin_comb_Q_T([1.0, -2.0, 0.5], Q_T) ≈ [
        1.9862687479070797,
        0.6259141580861511,
        0.2815640364800571,
        -0.91306605519864,
    ] atol=3e-15
    @test expm_lanczos(E, V, Q_T; a=-0.25) ≈ [
        0.28370463815222685,
        0.49009609392027137,
        -0.2301972310236773,
        0.09222361279644782,
    ] atol=5e-15

    E_iter, V_iter, Q_iter = lanczos_iter(A, v0, 3)
    @test E_iter ≈ E atol=1e-14
    @test abs.(V_iter) ≈ abs.(V) atol=1e-14
    @test reduce(vcat, permutedims.(collect(Q_iter))) ≈ Q_T atol=1e-14

    B = Diagonal([2.0, -1.0, 0.5, 3.0])
    beta = [0.0, 0.5, 2.0]
    ftlm, ftlm_identity = ftlm_static_iteration(
        Dict("A" => A, "B" => B),
        E,
        V,
        Q_T;
        beta,
    )
    ltlm, ltlm_identity = ltlm_static_iteration(
        Dict("A" => A, "B" => B),
        E,
        V,
        Q_T;
        beta,
    )
    @test ftlm_identity ≈ [
        0.9999999999999998,
        0.38217846291319085,
        0.02560632611811718,
    ] atol=7e-16
    @test ftlm["A"] ≈ [
        1.952,
        0.724811818302737,
        0.04323175995830908,
    ] atol=2e-15
    @test ftlm["B"] ≈ [
        -0.12000000000000034,
        -0.02325674028781951,
        0.01521990440418417,
    ] atol=1e-15
    @test ltlm_identity ≈ ftlm_identity atol=7e-16
    @test ltlm["A"] ≈ [
        1.9519999999999995,
        0.724811818302737,
        0.04323175995830904,
    ] atol=2e-15
    @test ltlm["B"] ≈ [
        -0.12000000000000059,
        -0.02720657100056404,
        0.01224408682892336,
    ] atol=1e-15

    basis = SpinBasis1D(3; nup=1, pauli=false)
    H = Hamiltonian(
        basis,
        [
            OperatorTerm("+-", [(0.5, 1, 2), (0.5, 2, 3)]),
            OperatorTerm("-+", [(0.5, 1, 2), (0.5, 2, 3)]),
            OperatorTerm("z", [(0.3, 1)]),
        ],
    )
    EH, VH, QH = lanczos_full(H, [1.0, 0.5, -0.25], 2)
    @test length(EH) == 2
    @test size(VH) == (2, 2)
    @test size(QH) == (2, 3)
    @test QH * QH' ≈ Matrix(I, 2, 2) atol=3e-15
    H_ftlm, H_identity = ftlm_static_iteration(
        Dict("H" => H),
        EH,
        VH,
        QH;
        beta=[0.0, 1.0],
    )
    @test length(H_ftlm["H"]) == 2
    @test length(H_identity) == 2
end
