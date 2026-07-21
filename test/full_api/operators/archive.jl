@testset "QuantumOperator native archive round trip" begin
    basis = SpinBasis1D(2)
    operator = QuantumOperator(
        basis,
        Dict(
            "z" => [OperatorTerm("z", [(1.0, 1)])],
            "x" => [
                OperatorTerm("+", [(0.5, 2)]),
                OperatorTerm("-", [(0.5, 2)]),
            ],
        ),
    )
    mktempdir() do directory
        archive = joinpath(directory, "operator.zip")
        save_zip(archive, operator)
        restored = load_zip(archive)
        @test restored.basis == basis
        @test restored.components == operator.components
        @test toarray(restored; pars=Dict("z" => 0.3, "x" => 0.7)) ==
            toarray(operator; pars=Dict("z" => 0.3, "x" => 0.7))

        archive_without_basis = joinpath(directory, "operator-no-basis.zip")
        save_zip(archive_without_basis, operator; save_basis=false)
        inferred = load_zip(archive_without_basis)
        @test inferred.basis.L == 2
        @test inferred.basis.nup === nothing
        @test inferred.components == operator.components
    end
end

@testset "QuantumOperator Python-compatible archive round trip" begin
    basis = SpinBasis1D(2)
    dense_component = ComplexF64[
        1 2im 0 0
        -2im 3 0 0
        0 0 4 0.5
        0 0 0.5 5
    ]
    sparse_component = sparse(ComplexF64[
        0 0.25 0 0
        0.25 0 0 0
        0 0 0 -0.75im
        0 0 0.75im 0
    ])
    operator = QuantumOperator(
        basis,
        Dict(
            "dense" => dense_component,
            "sparse" => sparse_component,
        );
        matrix_formats=Dict("dense" => :dense, "sparse" => :csc),
    )
    mktempdir() do directory
        archive = joinpath(directory, "python-compatible.zip")
        save_zip(
            archive,
            operator;
            save_basis=false,
            format=:python,
        )
        restored = load_zip(archive)
        @test restored.basis == basis
        @test restored.components["dense"] == dense_component
        @test restored.components["sparse"] isa SparseMatrixCSC
        @test restored.components["sparse"] == sparse_component
        @test toarray(
            restored;
            pars=Dict("dense" => 0.3, "sparse" => 0.7),
        ) == toarray(
            operator;
            pars=Dict("dense" => 0.3, "sparse" => 0.7),
        )
    end
end
