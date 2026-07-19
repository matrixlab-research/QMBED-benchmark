@testset "tensor and photon basis protocols" begin
    spin = SpinBasis1D(1; pauli=false)
    boson = BosonBasis1D(1; sps=3)
    tensor = TensorBasis(spin, boson)

    @test tensor.N == (1, 1)
    @test tensor.Ns == 6
    @test tensor.sps == (2, 3)
    @test tensor.basis_left === spin
    @test tensor.basis_right === boson
    @test tensor.blocks[:left] == spin.blocks
    @test tensor.blocks[:right] == boson.blocks
    @test size(projection_matrix(tensor)) == (6, 6)
    @test projection_matrix(tensor) ≈ Matrix(I, 6, 6)
    @test state_index(tensor, 1, 1) == 1
    @test state_index(tensor, 2, 2) == 5

    combined = operator_matrix(tensor, "z|n", [(1.0, 1, 1)])
    expected = kron(
        operator_matrix(spin, "z", [(1.0, 1)]),
        operator_matrix(boson, "n", [(1.0, 1)]),
    )
    @test combined ≈ expected atol=2e-16
    output = zeros(ComplexF64, 6, 6)
    inplace_op!(output, tensor, "z|n", [(1.0, 1, 1)])
    @test output == combined

    state = zeros(ComplexF64, 6)
    state[state_index(tensor, 1, 1)] = inv(sqrt(2))
    state[state_index(tensor, 2, 2)] = inv(sqrt(2))
    rho_left, rho_right = partial_trace(tensor, state; return_rdm=:both)
    @test rho_left ≈ Diagonal([0.5, 0.5]) atol=3e-16
    @test rho_right ≈ Diagonal([0.5, 0.5, 0.0]) atol=3e-16
    @test ent_entropy(tensor, state)["Sent_A"] ≈ log(2) atol=3e-16

    photon = PhotonBasis(SpinBasis1D, 1; Nph=2, pauli=false)
    @test photon.N == 2
    @test photon.Ns == 6
    @test photon.chain_N == 1
    @test photon.chain_Ns == 2
    @test photon.particle_N == 1
    @test photon.particle_Ns == 2
    @test photon.particle_sps == 2
    @test photon.photon_Ns == 3
    @test photon.photon_sps == 3
    @test photon.sps == (2, 3)
    @test photon.blocks[:Nph] == 2
    @test operator_matrix(photon, "z|n", [(1.0, 1, 1)]) ≈ expected atol=2e-16
end
