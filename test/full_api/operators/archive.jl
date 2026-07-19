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
