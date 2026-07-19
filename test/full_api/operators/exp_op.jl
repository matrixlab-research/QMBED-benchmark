@testset "ExpOp full protocol and oracle" begin
    O = [0.0 1.0; -1.0 0.0]
    vector = [1.0, 0.0]
    observable = Diagonal([2.0, 3.0])
    expO = ExpOp(O; a=0.5)

    @test isexp_op(expO)
    @test expO.Ns == 2
    @test expO.get_shape == (2, 2)
    @test expO.ndim == 2
    @test expO.O == O
    @test expO.a == 0.5
    @test expO.grid === nothing
    @test expO.step === nothing
    @test expO.iterate == false
    @test get_mat(expO) ≈ [
        0.8775825618903728 0.479425538604203
        -0.47942553860420306 0.8775825618903728
    ] atol=8e-16
    @test apply(expO, vector) ≈
        [0.8775825618903728, -0.479425538604203] atol=8e-16
    @test right_apply(expO, vector) ≈
        [0.8775825618903728, 0.479425538604203] atol=8e-16
    @test sandwich(expO, observable) ≈ [
        2.2298488470659295 0.4207354924039482
        0.4207354924039482 2.7701511529340697
    ] atol=2e-15

    @test get_mat(expO.H) ≈ get_mat(expO)' atol=8e-16
    @test get_mat(expO.T) ≈ transpose(get_mat(expO)) atol=8e-16
    @test get_mat(conj(expO)) ≈ conj(get_mat(expO)) atol=8e-16
    copied = copy(expO)
    @test copied !== expO
    @test get_mat(copied) == get_mat(expO)

    set_a!(copied, 0.25)
    set_grid!(copied, 0.0, 1.0; num=3, endpoint=true)
    @test copied.grid == [0.0, 0.5, 1.0]
    @test copied.step == 0.5
    expected_grid = [
        1.0 0.992197667229329 0.9689124217106447
        0.0 -0.12467473338522769 -0.24740395925452294
    ]
    @test apply(copied, vector) ≈ expected_grid atol=8e-16
    @test size(sandwich(copied, observable)) == (2, 2, 3)

    set_iterate!(copied, true)
    @test collect(apply(copied, vector)) ≈
        collect(eachcol(expected_grid)) atol=8e-16
    unset_grid!(copied)
    @test copied.grid === nothing
    @test copied.step === nothing
    @test !copied.iterate
    @test_throws ArgumentError set_iterate!(copied, true)

    basis = SpinBasis1D(2; pauli=false)
    H = Hamiltonian(basis, [OperatorTerm("z", [(1.0, 1)])])
    expH = ExpOp(H; a=-0.2im)
    @test apply(expH, ones(ComplexF64, length(basis))) ≈
        exp(-0.2im * Matrix(H)) * ones(ComplexF64, length(basis)) atol=2e-15
end
