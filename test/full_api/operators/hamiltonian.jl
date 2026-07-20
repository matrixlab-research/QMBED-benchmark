using LinearAlgebra

@testset "XXZ chain: basis -> operator -> spectrum" begin
    basis = SpinBasis1D(4; nup=2, pauli=false)
    jxy = sqrt(2.0)
    jzz = 1.0
    hz = inv(sqrt(3.0))

    terms = [
        OperatorTerm("+-", [(jxy / 2, i, i + 1) for i in 1:3]),
        OperatorTerm("-+", [(jxy / 2, i, i + 1) for i in 1:3]),
        OperatorTerm("zz", [(jzz, i, i + 1) for i in 1:3]),
        OperatorTerm("z", [(hz, i) for i in 1:4]),
    ]
    H = Hamiltonian(basis, terms)

    expected_spectrum = [
        -2.06671263224485,
        -1.116025403784439,
        -0.25,
        0.13427005520587385,
        0.6160254037844387,
        1.1824425770389764,
    ]

    @test sort(eigvals(H)) ≈ expected_spectrum atol=1e-12 rtol=1e-12
    @test tr(Matrix(H)) ≈ -1.5 atol=1e-12
    @test norm(Matrix(H)) ≈ 2.715695122800054 atol=1e-12
end

@testset "Hamiltonian static protocol" begin
    basis = SpinBasis1D(3; nup=1, pauli=false)
    H = Hamiltonian(
        basis,
        [
            OperatorTerm("+-", [(0.7, 1, 2), (0.4, 2, 3)]),
            OperatorTerm("-+", [(0.7, 1, 2), (0.4, 2, 3)]),
            OperatorTerm("z", [(0.2, 1), (-0.3, 3)]),
        ],
    )
    matrix = Matrix(H)
    @test ishamiltonian(H)
    @test H.Ns == 3
    @test H.shape == (3, 3) == H.get_shape
    @test H.ndim == 2
    @test H.dtype == Float64
    @test H.basis === basis
    @test H.static == matrix
    @test isempty(H.dynamic)
    @test H.is_dense
    @test H.nbytes >= sizeof(matrix)
    @test Matrix(H.H) == matrix'
    @test Matrix(H.T) == transpose(matrix)

    @test as_dense_format(H) === H
    @test as_dense_format(H; copy=true) !== H
    @test Matrix(as_sparse_format(H)) == matrix
    @test aslinearoperator(H) === H
    @test check_is_dense(H)
    @test astype(H, ComplexF64).dtype == ComplexF64
    @test Matrix(conj(H)) == conj(matrix)
    @test Matrix(transpose(H)) == transpose(matrix)
    @test Matrix(adjoint(H)) == matrix'
    @test copy(H) !== H
    @test diagonal(H) == diag(matrix)

    vector = normalize(ComplexF64[1, 2im, -0.5])
    @test apply(H, vector) == matrix * vector
    out = zeros(ComplexF64, 3)
    @test apply(H, vector; out, a=-0.3) === out
    @test out ≈ -0.3matrix * vector atol=3e-16
    @test right_apply(H, vector) == vec(transpose(vector) * matrix)
    values, vectors = eigh(H)
    @test matrix * vectors ≈ vectors * Diagonal(values) atol=2e-14
    selected, selected_vectors = eigsh(H; k=2, which=:SA)
    @test selected ≈ values[1:2] atol=2e-13
    @test size(selected_vectors) == (3, 2)
    @test eigvals(H) ≈ values atol=1e-14

    expected = dot(vector, matrix * vector)
    @test expt_value(H, vector) ≈ expected
    @test matrix_ele(H, vector, vector) ≈ expected
    @test quant_fluct(H, vector) ≈
        dot(vector, matrix^2 * vector) - expected^2 atol=2e-15
    @test Matrix(project_to(H, Matrix{Float64}(I, 3, 3))) == matrix
    rotation = ExpOp(H; a=-0.17im)
    @test Matrix(rotate_by(H, rotation)) ≈
        get_mat(rotation)' * matrix * get_mat(rotation) atol=2e-14

    evolved = evolve(H, vector, 0.2, [0.2, 0.5, 1.3])
    @test evolved[:, 1] ≈ vector atol=2e-15
    @test all(isapprox(norm(column), 1.0; atol=2e-14) for column in eachcol(evolved))
    @test toarray(H) == matrix
    @test todense(H) == matrix
    @test Matrix(tocsc(H)) == matrix
    @test Matrix(tocsr(H)) == matrix
    @test tr(H) == tr(matrix)
    sparse_H = as_sparse_format(H; static_fmt=:csc)
    @test sparse_H.data isa SparseMatrixCSC
    @test !sparse_H.is_dense
    @test Matrix(sparse_H) == matrix
    @test H.data isa Matrix
    @test update_matrix_formats!(H, :csc, Dict()) === H
    @test H.data isa SparseMatrixCSC
    @test !check_is_dense(H)
    @test update_matrix_formats!(H, :dense, Dict()) === H
    @test H.data isa Matrix
    @test update_matrix_formats!(H, :csr, Dict()) === H
    @test H.data isa SparseMatrixCSR
end
