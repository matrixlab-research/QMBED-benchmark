@testset "ExpmMultiplyParallel oracle and protocol" begin
    generator = [0.0 -1.0; 1.0 0.0]
    operator = ExpmMultiplyParallel(generator, 0.25)

    @test operator.A === generator
    @test operator.a == 0.25
    @test apply(operator, [1.0, 0.0]) ≈
        [0.9689124217106447, 0.24740395925452294] atol=4e-16

    batch = [1.0 0.0; 0.0 1.0]
    @test operator * batch ≈ [
        0.9689124217106447 -0.24740395925452294
        0.24740395925452294 0.9689124217106447
    ] atol=4e-16

    output = ComplexF64[1, 0]
    set_a!(operator, -0.5im)
    @test operator.a == -0.5im
    @test apply(operator, output; overwrite_v=true) === output
    @test output ≈
        ComplexF64[1.1276259652063807, -0.5210953054937474im] atol=8e-16

    @test_throws ArgumentError ExpmMultiplyParallel(ones(2, 3))
    @test_throws DimensionMismatch apply(operator, ones(3))
    @test_throws DimensionMismatch apply(operator, ones(2); work_array=zeros(3))
end
