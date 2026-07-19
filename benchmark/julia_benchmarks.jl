#!/usr/bin/env julia

using LinearAlgebra
using QuSpin
using SparseArrays

BLAS.set_num_threads(1)

const SAMPLES = 15
const TARGET_SAMPLE_SECONDS = 0.08
const MAX_ITERATIONS = 1_000
const SINK = Ref{Any}(nothing)

function xxz_hamiltonian(
    length::Int,
    nup::Int;
    static_fmt=:dense,
)
    basis = SpinBasis1D(length; nup, pauli=false)
    jxy = sqrt(2.0)
    jzz = 1.0
    hz = inv(sqrt(3.0))
    terms = [
        OperatorTerm(
            "+-",
            [(jxy / 2.0, site, site + 1) for site in 1:(length - 1)],
        ),
        OperatorTerm(
            "-+",
            [(jxy / 2.0, site, site + 1) for site in 1:(length - 1)],
        ),
        OperatorTerm(
            "zz",
            [(jzz, site, site + 1) for site in 1:(length - 1)],
        ),
        OperatorTerm("z", [(hz, site) for site in 1:length]),
    ]
    return basis, Hamiltonian(basis, terms; static_fmt)
end

function deterministic_state(size::Int)
    state = ComplexF64[
        sin(0.17 * index) + im * cos(0.11 * index) for index in 1:size
    ]
    return state / norm(state)
end

struct BenchmarkCase
    name::String
    category::String
    comparison::String
    storage::String
    parameters::String
    function_value::Any
    supported::Bool
    note::String
end

function make_cases()
    basis_12, hamiltonian_12_dense = xxz_hamiltonian(12, 6; static_fmt=:dense)
    _, hamiltonian_12_csc = xxz_hamiltonian(12, 6; static_fmt=:csc)
    matrix_12_dense = Matrix(hamiltonian_12_dense)
    matrix_12_csc = hamiltonian_12_csc.data
    vector_12 = deterministic_state(length(basis_12))
    basis_10, hamiltonian_10_dense = xxz_hamiltonian(10, 5; static_fmt=:dense)
    _, hamiltonian_10_csc = xxz_hamiltonian(10, 5; static_fmt=:csc)
    matrix_10_dense = Matrix(hamiltonian_10_dense)
    matrix_10_csc = hamiltonian_10_csc.data
    vector_10 = deterministic_state(length(basis_10))
    eigsh_seed_10 = normalize(real.(vector_10))
    basis_8, hamiltonian_8 = xxz_hamiltonian(8, 4; static_fmt=:csc)
    vector_8 = deterministic_state(length(basis_8))
    full_basis_12 = SpinBasis1D(12; pauli=false)
    full_state_12 = deterministic_state(length(full_basis_12))
    times = collect(range(0.0, 1.0; length=9))

    return [
        BenchmarkCase(
            "spin_basis_construction",
            "api",
            "storage_independent",
            "n/a",
            "L=16;nup=8;dimension=12870",
            () -> SpinBasis1D(16; nup=8, pauli=false),
            true,
            "",
        ),
        BenchmarkCase(
            "xxz_hamiltonian_construction_dense",
            "integration",
            "controlled",
            "dense",
            "L=10;nup=5;dimension=252;open=true",
            () -> last(xxz_hamiltonian(10, 5; static_fmt=:dense)),
            true,
            "",
        ),
        BenchmarkCase(
            "xxz_hamiltonian_construction_sparse",
            "integration",
            "controlled",
            "csc",
            "L=10;nup=5;dimension=252;open=true",
            () -> last(xxz_hamiltonian(10, 5; static_fmt=:csc)),
            true,
            "Native CSC triplet assembly; no intermediate dense Hamiltonian.",
        ),
        BenchmarkCase(
            "hamiltonian_matvec_current_storage",
            "api",
            "current_backend",
            "csc",
            "L=12;nup=6;dimension=924",
            () -> hamiltonian_12_csc * vector_12,
            true,
            "",
        ),
        BenchmarkCase(
            "matrix_matvec_dense",
            "kernel",
            "controlled",
            "dense",
            "L=12;nup=6;dimension=924",
            () -> matrix_12_dense * vector_12,
            true,
            "",
        ),
        BenchmarkCase(
            "matrix_matvec_sparse_csc",
            "kernel",
            "controlled",
            "csc",
            "L=12;nup=6;dimension=924",
            () -> matrix_12_csc * vector_12,
            true,
            "",
        ),
        BenchmarkCase(
            "full_eigenspectrum_dense",
            "kernel",
            "controlled",
            "dense",
            "L=10;nup=5;dimension=252",
            () -> eigvals(Hermitian(matrix_10_dense)),
            true,
            "",
        ),
        BenchmarkCase(
            "lanczos_decomposition_dense",
            "integration",
            "controlled",
            "dense",
            "L=10;nup=5;dimension=252;m=32",
            () -> lanczos_full(matrix_10_dense, vector_10, 32),
            true,
            "",
        ),
        BenchmarkCase(
            "lanczos_decomposition_sparse_csc",
            "integration",
            "controlled",
            "csc",
            "L=10;nup=5;dimension=252;m=32",
            () -> lanczos_full(matrix_10_csc, vector_10, 32),
            true,
            "",
        ),
        BenchmarkCase(
            "partial_eigenspectrum_sparse_csc",
            "integration",
            "controlled",
            "csc",
            "L=10;nup=5;dimension=252;k=4;which=SA",
            () -> eigsh(
                hamiltonian_10_csc;
                k=4,
                which=:SA,
                v0=eigsh_seed_10,
                tol=1e-10,
            ),
            true,
            "ARPACK iterative eigensolve on the native CSC Hamiltonian.",
        ),
        BenchmarkCase(
            "static_time_evolution_current_storage",
            "integration",
            "current_backend",
            "csc_input_dense_solver",
            "L=8;nup=4;dimension=70;times=9;tmax=1",
            () -> evolve(hamiltonian_8, vector_8, 0.0, times),
            true,
            "The Hamiltonian is CSC; exact evolution currently uses a dense full eigendecomposition.",
        ),
        BenchmarkCase(
            "entanglement_entropy",
            "integration",
            "storage_independent",
            "state_vector",
            "L=12;dimension=4096;subsystem=6",
            () -> ent_entropy(
                full_basis_12,
                full_state_12;
                sub_sys_A=collect(1:6),
            ),
            true,
            "",
        ),
    ]
end

function percentile(sorted_values::Vector{Float64}, probability::Float64)
    position = probability * (length(sorted_values) - 1) + 1
    lower = floor(Int, position)
    upper = ceil(Int, position)
    lower == upper && return sorted_values[lower]
    weight = position - lower
    return sorted_values[lower] * (1 - weight) + sorted_values[upper] * weight
end

function median_value(values::Vector{Float64})
    ordered = sort(values)
    midpoint = length(ordered) ÷ 2
    return isodd(length(ordered)) ?
        ordered[midpoint + 1] :
        (ordered[midpoint] + ordered[midpoint + 1]) / 2
end

function time_case(case::BenchmarkCase)
    if !case.supported
        return (
            language="julia",
            benchmark=case.name,
            category=case.category,
            comparison=case.comparison,
            storage=case.storage,
            supported="false",
            note=case.note,
            parameters=case.parameters,
            samples="",
            iterations_per_sample="",
            median_seconds="",
            mean_seconds="",
            stdev_seconds="",
            min_seconds="",
            p05_seconds="",
            p25_seconds="",
            p75_seconds="",
            p95_seconds="",
            max_seconds="",
            median_allocated_bytes="",
            runtime="Julia $(VERSION); QuSpin candidate",
        )
    end

    for _ in 1:3
        SINK[] = case.function_value()
    end

    started_ns = time_ns()
    SINK[] = case.function_value()
    estimate = max((time_ns() - started_ns) / 1.0e9, 1.0e-9)
    iterations = clamp(
        ceil(Int, TARGET_SAMPLE_SECONDS / estimate),
        1,
        MAX_ITERATIONS,
    )

    timings = Float64[]
    allocated_bytes = Int[]
    for _ in 1:SAMPLES
        GC.gc()
        measurement = @timed begin
            for _ in 1:iterations
                SINK[] = case.function_value()
            end
        end
        push!(timings, measurement.time / iterations)
        push!(allocated_bytes, measurement.bytes ÷ iterations)
    end

    ordered = sort(timings)
    mean_value = sum(timings) / length(timings)
    stdev_value = sqrt(
        sum((value - mean_value)^2 for value in timings) /
        (length(timings) - 1),
    )
    return (
        language="julia",
        benchmark=case.name,
        category=case.category,
        comparison=case.comparison,
        storage=case.storage,
        supported="true",
        note=case.note,
        parameters=case.parameters,
        samples=length(timings),
        iterations_per_sample=iterations,
        median_seconds=median_value(timings),
        mean_seconds=mean_value,
        stdev_seconds=stdev_value,
        min_seconds=ordered[1],
        p05_seconds=percentile(ordered, 0.05),
        p25_seconds=percentile(ordered, 0.25),
        p75_seconds=percentile(ordered, 0.75),
        p95_seconds=percentile(ordered, 0.95),
        max_seconds=ordered[end],
        median_allocated_bytes=round(Int, median_value(Float64.(allocated_bytes))),
        runtime="Julia $(VERSION); QuSpin candidate",
    )
end

function csv_cell(value)
    text = string(value)
    return occursin(r"[\",\n]", text) ?
        "\"" * replace(text, "\"" => "\"\"") * "\"" :
        text
end

function main()
    output_index = findfirst(==("--output"), ARGS)
    output_index === nothing &&
        error("usage: julia_benchmarks.jl --output PATH")
    output_index < length(ARGS) ||
        error("--output requires a path")
    output = ARGS[output_index + 1]
    mkpath(dirname(output))
    rows = Any[time_case(case) for case in make_cases()]
    columns = propertynames(first(rows))
    open(output, "w") do io
        println(io, join(columns, ","))
        for row in rows
            println(
                io,
                join((csv_cell(getproperty(row, column)) for column in columns), ","),
            )
        end
    end
end

main()
