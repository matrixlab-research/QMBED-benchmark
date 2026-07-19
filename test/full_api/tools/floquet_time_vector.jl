@testset "FloquetTimeVector oracle and protocol" begin
    times = FloquetTimeVector(2π, 2; len_T=4, N_up=1, N_down=1)

    @test times.N == 4
    @test times.len_T == 4
    @test times.T ≈ 1.0 atol=2e-16
    @test times.vals ≈ collect(-1.0:0.25:3.0) atol=2e-16
    @test times.len == 17
    @test times.shape == (17,)
    @test times.dt ≈ 0.25 atol=2e-16
    @test times.i ≈ -1.0 atol=2e-16
    @test times.f ≈ 3.0 atol=2e-16
    @test times.tot ≈ 4.0 atol=2e-16
    @test times.strobo.inds == [1, 5, 9, 13, 17]
    @test times.strobo.vals ≈ [-1.0, 0.0, 1.0, 2.0, 3.0] atol=2e-16
    @test times.up.vals ≈ [-1.0, -0.75, -0.5, -0.25] atol=2e-16
    @test times.constant.vals ≈ collect(0.0:0.25:2.0) atol=2e-16
    @test times.down.vals ≈ collect(2.25:0.25:3.0) atol=2e-16
    @test get_coordinates(times, 7) == (2, 3)

    @test collect(times) == times.vals
    @test times[7] == 0.5
    @test times * 2 == times.vals * 2
    @test times / 2 == times.vals / 2
    @test_throws BoundsError get_coordinates(times, 0)
    @test_throws ArgumentError FloquetTimeVector(0.0, 2)
    @test_throws ArgumentError FloquetTimeVector(2π, 0)
end
