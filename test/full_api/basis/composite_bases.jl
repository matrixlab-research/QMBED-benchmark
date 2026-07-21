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
    product = zeros(ComplexF64, 6)
    product[state_index(tensor, 1, 1)] = 1
    batch = hcat(state, product)
    batch_result = ent_entropy(
        tensor,
        batch;
        enforce_pure=true,
        return_rdm=:both,
    )
    @test batch_result["Sent_A"] ≈ [log(2), 0.0] atol=4e-14
    @test size(batch_result["rdm_A"]) == (2, 2, 2)

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
    photon_batch = ent_entropy(
        photon,
        batch;
        enforce_pure=true,
        return_rdm=:A,
    )
    @test photon_batch["Sent_A"] ≈ [log(2), 0.0] atol=4e-14
end
@testset "three-factor tensor and conserving photon basis" begin
    left = SpinBasis1D(1; pauli=false)
    middle = BosonBasis1D(1; sps=3)
    right = SpinBasis1D(1; pauli=false)
    tensor = TensorBasis(left, middle, right)
    matrix = operator_matrix(tensor, "x|n|z", [(0.4, 1, 1, 1)])
    @test matrix ≈ kron(
        operator_matrix(left, "x", [(0.4, 1)]),
        operator_matrix(middle, "n", [(1.0, 1)]),
        operator_matrix(right, "z", [(1.0, 1)]),
    ) atol=3e-15

    photon = PhotonBasis(SpinBasis1D, 4; Ntot=2, pauli=false)
    @test length(photon) == photon_Hspace_dim(4, 2, nothing) == 11
    total = sum(
        operator_matrix(photon, "z|", [(1.0, site)])
        for site in 1:4
    ) + operator_matrix(photon, "|n", [(1.0, 1)])
    @test total ≈ zeros(size(total)) atol=4e-15
end
