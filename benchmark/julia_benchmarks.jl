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
    return basis, Hamiltonian(
        basis,
        terms;
        static_fmt,
        check_symm=false,
        check_herm=false,
        check_pcon=false,
    )
end

function deterministic_state(size::Int)
    state = ComplexF64[
        sin(0.17 * index) + im * cos(0.11 * index) for index in 1:size
    ]
    return state / norm(state)
end

function deterministic_state_batch(size::Int, count::Int)
    states = ComplexF64[
        sin(0.17 * row + column) + im * cos(0.11 * row - column)
        for row in 1:size, column in 1:count
    ]
    return states ./ sqrt.(sum(abs2, states; dims=1))
end

function general_spin_maps(lx::Int, ly::Int)
    site(x, y) = x + lx * y
    translation_x = [
        site(mod(x + 1, lx), y)
        for y in 0:(ly - 1) for x in 0:(lx - 1)
    ]
    translation_y = [
        site(x, mod(y + 1, ly))
        for y in 0:(ly - 1) for x in 0:(lx - 1)
    ]
    return translation_x, translation_y
end

function general_spin_basis()
    translation_x, translation_y = general_spin_maps(4, 3)
    return SpinBasisGeneral(
        12;
        Nup=6,
        pauli=false,
        kxblock=(translation_x, 1),
        kyblock=(translation_y, 1),
        block_order=[:kxblock, :kyblock],
    )
end

function general_spin_terms()
    lx, ly = 4, 3
    site(x, y) = x + lx * y
    bonds = [
        (site(x, y) + 1, site(mod(x + 1, lx), y) + 1)
        for y in 0:(ly - 1) for x in 0:(lx - 1)
    ]
    append!(
        bonds,
        [
            (site(x, y) + 1, site(x, mod(y + 1, ly)) + 1)
            for y in 0:(ly - 1) for x in 0:(lx - 1)
        ],
    )
    return [
        OperatorTerm("+-", [(0.5, left, right) for (left, right) in bonds]),
        OperatorTerm("-+", [(0.5, left, right) for (left, right) in bonds]),
        OperatorTerm("zz", [(1.0, left, right) for (left, right) in bonds]),
    ]
end

function general_spin_hamiltonian(basis, terms)
    return Hamiltonian(
        basis,
        terms;
        static_fmt=:csc,
        check_herm=false,
        check_symm=false,
        check_pcon=false,
    )
end

function spin_one_hamiltonian(length::Int)
    basis = SpinBasis1D(
        length;
        S="1",
        Nup=length,
        pauli=false,
    )
    bonds = [(site, mod1(site + 1, length)) for site in 1:length]
    terms = [
        OperatorTerm("+-", [(0.5, left, right) for (left, right) in bonds]),
        OperatorTerm("-+", [(0.5, left, right) for (left, right) in bonds]),
        OperatorTerm("zz", [(1.0, left, right) for (left, right) in bonds]),
    ]
    return Hamiltonian(
        basis,
        terms;
        static_fmt=:csc,
        check_herm=false,
        check_symm=false,
        check_pcon=false,
    )
end

function validate_hamiltonian_fingerprint(
    operator,
    expected_dimension::Int,
    expected_trace::Float64,
    expected_norm::Float64,
)
    matrix = Matrix(operator)
    return size(matrix) == (expected_dimension, expected_dimension) &&
        ishermitian(matrix) &&
        isapprox(real(tr(matrix)), expected_trace; rtol=1e-12, atol=1e-12) &&
        isapprox(norm(matrix), expected_norm; rtol=1e-12, atol=1e-12)
end

function validate_batch_entropy(result)
    entropy = result["Sent_A"]
    return size(entropy) == (8,) &&
        all(isfinite, entropy) &&
        isapprox(sum(entropy), 7.360709128824066; rtol=1e-12, atol=1e-12)
end

function validate_batch_evolution(states)
    return size(states) == (70, 4, 9) &&
        maximum(abs.(sum(abs2, states; dims=1) .- 1)) < 2e-10
end

function spinful_terms(length::Int)
    bonds = [(site, site + 1) for site in 1:(length - 1)]
    return [
        OperatorTerm("+-|", [(-1.0, left, right) for (left, right) in bonds]),
        OperatorTerm("-+|", [(1.0, left, right) for (left, right) in bonds]),
        OperatorTerm("|+-", [(-1.0, left, right) for (left, right) in bonds]),
        OperatorTerm("|-+", [(1.0, left, right) for (left, right) in bonds]),
        OperatorTerm("n|n", [(2.0, site, site) for site in 1:length]),
    ]
end

function spinful_hamiltonian(length::Int)
    particles = length ÷ 2
    basis = SpinfulFermionBasis1D(
        length;
        Nf=(particles, particles),
    )
    operator = Hamiltonian(
        basis,
        spinful_terms(length);
        static_fmt=:csc,
        check_herm=false,
        check_symm=false,
        check_pcon=false,
    )
    return basis, operator
end

function spinless_hamiltonian(length::Int)
    particles = length ÷ 2
    bonds = [(site, site + 1) for site in 1:(length - 1)]
    hopping = [isodd(site) ? 0.6 : 1.0 for site in 1:(length - 1)]
    basis = SpinlessFermionBasis1D(length; Nf=particles)
    terms = [
        OperatorTerm(
            "+-",
            [
                (-hopping[index], left, right)
                for (index, (left, right)) in enumerate(bonds)
            ],
        ),
        OperatorTerm(
            "-+",
            [
                (hopping[index], left, right)
                for (index, (left, right)) in enumerate(bonds)
            ],
        ),
        OperatorTerm(
            "nn",
            [(2.0, left, right) for (left, right) in bonds],
        ),
    ]
    operator = Hamiltonian(
        basis,
        terms;
        static_fmt=:csc,
        check_herm=false,
        check_symm=false,
        check_pcon=false,
    )
    return basis, operator
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
    validator::Any
end

BenchmarkCase(
    name::String,
    category::String,
    comparison::String,
    storage::String,
    parameters::String,
    function_value,
    supported::Bool,
    note::String,
) = BenchmarkCase(
    name,
    category,
    comparison,
    storage,
    parameters,
    function_value,
    supported,
    note,
    nothing,
)

function make_cases()
    basis_12, hamiltonian_12_dense = xxz_hamiltonian(12, 6; static_fmt=:dense)
    _, hamiltonian_12_csc = xxz_hamiltonian(12, 6; static_fmt=:csc)
    _, hamiltonian_12_csr = xxz_hamiltonian(12, 6; static_fmt=:csr)
    _, hamiltonian_12_dia = xxz_hamiltonian(12, 6; static_fmt=:dia)
    matrix_12_dense = Matrix(hamiltonian_12_dense)
    matrix_12_csc = hamiltonian_12_csc.data
    matrix_12_csr = hamiltonian_12_csr.data
    matrix_12_dia = hamiltonian_12_dia.data
    vector_12 = deterministic_state(length(basis_12))
    linear_operator_12 = QuantumLinearOperator(
        basis_12,
        hamiltonian_12_csc.terms,
        check_symm=false,
        check_herm=false,
        check_pcon=false,
    )
    basis_10, hamiltonian_10_dense = xxz_hamiltonian(10, 5; static_fmt=:dense)
    _, hamiltonian_10_csc = xxz_hamiltonian(10, 5; static_fmt=:csc)
    matrix_10_dense = Matrix(hamiltonian_10_dense)
    matrix_10_csc = hamiltonian_10_csc.data
    vector_10 = deterministic_state(length(basis_10))
    eigsh_seed_10 = normalize(real.(vector_10))
    basis_8, hamiltonian_8 = xxz_hamiltonian(8, 4; static_fmt=:csc)
    vector_8 = deterministic_state(length(basis_8))
    batch_vector_8 = deterministic_state_batch(length(basis_8), 4)
    general_basis_12 = general_spin_basis()
    general_terms_12 = general_spin_terms()
    full_basis_12 = SpinBasis1D(12; pauli=false)
    full_state_12 = deterministic_state(length(full_basis_12))
    full_basis_10 = SpinBasis1D(10; pauli=false)
    batch_state_10 = deterministic_state_batch(length(full_basis_10), 8)
    times = collect(range(0.0, 1.0; length=9))
    spinful_basis_6, spinful_hamiltonian_6 = spinful_hamiltonian(6)
    spinful_state_6 = deterministic_state(length(spinful_basis_6))
    spinful_qlo_6 = QuantumLinearOperator(
        spinful_basis_6,
        spinful_terms(6);
        check_herm=false,
        check_symm=false,
        check_pcon=false,
    )
    exp_dimension = 256
    exp_matrix = spdiagm(
        -1 => fill(0.15, exp_dimension - 1),
        0 => collect(range(-1.0, 1.0; length=exp_dimension)),
        1 => fill(0.15, exp_dimension - 1),
    )
    exp_state = deterministic_state(exp_dimension)
    exp_operator = ExpOp(exp_matrix; a=-0.2im)
    exp_grid_dimension = 128
    exp_grid_matrix = spdiagm(
        -1 => fill(0.15, exp_grid_dimension - 1),
        0 => collect(range(-1.0, 1.0; length=exp_grid_dimension)),
        1 => fill(0.15, exp_grid_dimension - 1),
    )
    exp_grid_state = deterministic_state(exp_grid_dimension)
    exp_grid_operator = ExpOp(
        exp_grid_matrix;
        a=-0.2im,
        start=0.0,
        stop=1.0,
        num=9,
    )
    ensemble_dimension = 256
    ensemble_observable = ComplexF64[
        sin(0.013row + 0.017column) +
        im * cos(0.019row - 0.023column)
        for row in 0:(ensemble_dimension - 1),
            column in 0:(ensemble_dimension - 1)
    ]
    ensemble_observable = ensemble_observable + ensemble_observable'
    ensemble_state = deterministic_state(ensemble_dimension)
    ensemble_energies = collect(1.0:ensemble_dimension)
    ensemble_vectors = Matrix{ComplexF64}(I, ensemble_dimension, ensemble_dimension)
    symmetry_basis_14 =
        SpinBasis1D(14; nup=7, pauli=false, kblock=0)
    symmetry_terms_14 = [
        OperatorTerm(
            "+-",
            [(0.5, site, mod1(site + 1, 14)) for site in 1:14],
        ),
        OperatorTerm(
            "-+",
            [(0.5, site, mod1(site + 1, 14)) for site in 1:14],
        ),
    ]

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
            "symmetry_basis_construction",
            "integration",
            "storage_independent",
            "projector",
            "L=14;nup=7;kblock=3;parent_dimension=3432",
            () -> SpinBasis1D(
                14;
                nup=7,
                pauli=false,
                kblock=3,
            ),
            true,
            "Constructs an orthonormal translation-sector projector.",
        ),
        BenchmarkCase(
            "symmetry_hamiltonian_construction_sparse",
            "integration",
            "controlled",
            "csc",
            "L=14;nup=7;kblock=0;dimension=$(length(symmetry_basis_14))",
            () -> Hamiltonian(
                symmetry_basis_14,
                symmetry_terms_14;
                static_fmt=:csc,
                check_herm=false,
                check_symm=false,
                check_pcon=false,
            ),
            true,
            "Direct reduced-sector triplet assembly without a parent Hamiltonian.",
        ),
        BenchmarkCase(
            "general_2d_basis_construction",
            "integration",
            "storage_independent",
            "projector",
            "Lx=4;Ly=3;Nup=6;kx=1;ky=1;dimension=75",
            general_spin_basis,
            true,
            "Two commuting translation maps through SpinBasisGeneral.",
            basis -> length(basis) == 75,
        ),
        BenchmarkCase(
            "general_2d_hamiltonian_construction_sparse",
            "integration",
            "controlled",
            "csc",
            "Lx=4;Ly=3;Nup=6;kx=1;ky=1;dimension=75",
            () -> general_spin_hamiltonian(general_basis_12, general_terms_12),
            true,
            "Periodic 2D XXZ assembly in a two-map symmetry sector.",
            operator -> validate_hamiltonian_fingerprint(
                operator,
                75,
                -41.0,
                18.30300521772313,
            ),
        ),
        BenchmarkCase(
            "higher_spin_hamiltonian_construction_sparse",
            "integration",
            "controlled",
            "csc",
            "L=8;S=1;Nup=8;dimension=1107;periodic=true",
            () -> spin_one_hamiltonian(8),
            true,
            "Spin-one XXZ construction with exact angular-momentum factors.",
            operator -> validate_hamiltonian_fingerprint(
                operator,
                1107,
                -816.0,
                112.7829774389737,
            ),
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
            "matrix_matvec_sparse_csr",
            "kernel",
            "controlled",
            "csr",
            "L=12;nup=6;dimension=924",
            () -> matrix_12_csr * vector_12,
            true,
            "",
        ),
        BenchmarkCase(
            "matrix_matvec_sparse_dia",
            "kernel",
            "controlled",
            "dia",
            "L=12;nup=6;dimension=924",
            () -> matrix_12_dia * vector_12,
            true,
            "",
        ),
        BenchmarkCase(
            "quantum_linear_operator_construction",
            "integration",
            "controlled",
            "matrix_free",
            "L=12;nup=6;dimension=924",
            () -> QuantumLinearOperator(
                basis_12,
                hamiltonian_12_csc.terms,
                check_symm=false,
                check_herm=false,
                check_pcon=false,
            ),
            true,
            "Off-diagonal terms stay matrix-free; one diagonal vector is " *
            "precomputed to accelerate repeated actions.",
        ),
        BenchmarkCase(
            "quantum_linear_operator_matvec",
            "kernel",
            "controlled",
            "matrix_free",
            "L=12;nup=6;dimension=924",
            () -> linear_operator_12 * vector_12,
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
            "csc_krylov",
            "L=8;nup=4;dimension=70;times=9;tmax=1",
            () -> evolve(hamiltonian_8, vector_8, 0.0, times),
            true,
            "Native CSC Arnoldi exponential action; no full eigendecomposition.",
        ),
        BenchmarkCase(
            "static_batch_time_evolution_current_storage",
            "integration",
            "current_backend",
            "csc_krylov",
            "L=8;nup=4;dimension=70;states=4;times=9;tmax=1",
            () -> evolve(hamiltonian_8, batch_vector_8, 0.0, times),
            true,
            "Batched column-state evolution across one shared time grid.",
            validate_batch_evolution,
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
        BenchmarkCase(
            "batched_entanglement_entropy",
            "integration",
            "storage_independent",
            "state_matrix",
            "L=10;dimension=1024;states=8;subsystem=5",
            () -> ent_entropy(
                full_basis_10,
                batch_state_10;
                sub_sys_A=collect(1:5),
                density=false,
                enforce_pure=true,
            ),
            true,
            "Eight pure states evaluated through the batched entropy API.",
            validate_batch_entropy,
        ),
        BenchmarkCase(
            "spinful_hamiltonian_construction_sparse",
            "integration",
            "controlled",
            "csc",
            "L=6;Nf=(3,3);dimension=400",
            () -> last(spinful_hamiltonian(6)),
            true,
            "Native sparse assembly for a spinful fermion Hamiltonian.",
        ),
        BenchmarkCase(
            "spinless_hamiltonian_construction_sparse",
            "integration",
            "controlled",
            "csc",
            "L=16;Nf=8;dimension=12870;open=true",
            () -> last(spinless_hamiltonian(16)),
            true,
            "Column-wise transition-to-CSC assembly for an interacting " *
            "spinless-fermion Hamiltonian.",
            operator -> size(operator) == (12870, 12870) &&
                ishermitian(operator.data) && nnz(operator.data) > 12870,
        ),
        BenchmarkCase(
            "spinful_quantum_linear_operator_matvec",
            "kernel",
            "controlled",
            "matrix_free",
            "L=6;Nf=(3,3);dimension=400",
            () -> spinful_qlo_6 * spinful_state_6,
            true,
            "Encoded-state matrix-free action without per-column occupation copies.",
        ),
        BenchmarkCase(
            "expop_sparse_vector_action",
            "kernel",
            "controlled",
            "csc_expm_action",
            "dimension=256;a=-0.2im",
            () -> exp_operator * exp_state,
            true,
            "Sparse exponential action; no dense matrix exponential.",
        ),
        BenchmarkCase(
            "expop_sparse_grid_action",
            "integration",
            "controlled",
            "csc_expm_action",
            "dimension=128;times=9;a=-0.2im",
            () -> exp_grid_operator * exp_grid_state,
            true,
            "One Krylov analysis reused across the time grid.",
        ),
        BenchmarkCase(
            "diag_ensemble_quantum_fluctuation",
            "integration",
            "controlled",
            "dense",
            "dimension=256;delta_t=true;delta_q=true",
            () -> diag_ensemble(
                1,
                ensemble_state,
                ensemble_energies,
                ensemble_vectors;
                density=false,
                Obs=ensemble_observable,
                delta_t_Obs=true,
                delta_q_Obs=true,
            ),
            true,
            "Diagonal of the squared observable uses an O(n^2) contraction.",
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
            validation="unsupported",
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

    SINK[] = case.function_value()
    validation = if case.validator === nothing
        "smoke"
    else
        case.validator(SINK[]) ||
            error("validation failed for benchmark $(case.name)")
        "passed"
    end
    for _ in 1:2
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
        validation,
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
