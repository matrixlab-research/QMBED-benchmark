#!/usr/bin/env julia

using LinearAlgebra
using QuSpin
using SparseArrays

BLAS.set_num_threads(1)

const DEFAULT_SAMPLES = 5
const ITERATIONS_PER_SAMPLE = 1
const SINK = Ref{Any}(nothing)

struct PaperCase
    case_id::String
    family_id::String
    benchmark::String
    parameters::String
    setup::Any
    run::Any
    validate::Any
end

uniform_seed(dimension::Integer) =
    fill(inv(sqrt(Float64(dimension))), dimension)

periodic_bonds(length::Integer, coefficient::Real=1.0) = [
    (coefficient, site, mod1(site + 1, length))
    for site in 1:length
]

function xxz_terms(
    length::Integer;
    Jxy::Real=1.0,
    Jz::Real=1.0,
    fields=nothing,
)
    bonds = periodic_bonds(length)
    terms = OperatorTerm[
        OperatorTerm(
            "+-",
            [(Jxy / 2, left, right) for (_, left, right) in bonds],
        ),
        OperatorTerm(
            "-+",
            [(Jxy / 2, left, right) for (_, left, right) in bonds],
        ),
        OperatorTerm(
            "zz",
            [(Jz, left, right) for (_, left, right) in bonds],
        ),
    ]
    fields === nothing ||
        push!(
            terms,
            OperatorTerm(
                "z",
                [(fields[site], site) for site in 1:length],
            ),
        )
    return terms
end

function sparse_hamiltonian(basis, terms)
    return Hamiltonian(
        basis,
        terms;
        static_fmt=:csc,
        check_symm=false,
        check_herm=false,
        check_pcon=false,
    )
end

function eigensystem_residual(H, values, vectors)
    residual = H * vectors - vectors * Diagonal(values)
    return norm(residual)
end

function validate_eigensystem(
    H,
    result;
    tolerance::Real=2e-7,
)
    values, vectors = result
    length(values) == size(vectors, 2) ||
        error("eigensolver returned inconsistent eigenpairs")
    all(isfinite, values) ||
        error("eigensolver returned non-finite eigenvalues")
    residual = eigensystem_residual(H, values, vectors)
    residual <= tolerance ||
        error("eigensystem residual $residual exceeds $tolerance")
    return "residual=$(residual)"
end

setup_mbl_shift_invert() = nothing

function run_mbl_shift_invert(_)
    L = 14
    basis = SpinBasis1D(L; nup=7, pauli=false)
    fields = [
        2.13,
        -1.77,
        0.31,
        3.24,
        -2.63,
        0.82,
        1.46,
        -3.17,
        2.71,
        -0.54,
        1.09,
        -2.28,
        0.67,
        2.94,
    ]
    H = sparse_hamiltonian(basis, xxz_terms(L; fields))
    seed = uniform_seed(length(basis))
    values, vectors = eigsh(
        H;
        k=6,
        sigma=0.0,
        which=:LM,
        v0=seed,
        ncv=32,
        tol=1e-9,
        maxiter=5_000,
    )
    residual = eigensystem_residual(H, values, vectors)
    return (
        valid=residual < 2e-7 && all(isfinite, values),
        residual=residual,
        value=values,
    )
end

validate_mbl_shift_invert(_, result) = result.valid

setup_xxz_lanczos_quench() = nothing

function run_xxz_lanczos_quench(_)
    L = 16
    basis = SpinBasis1D(L; nup=8, pauli=false)
    # The paired Python case uses open bonds.
    H = sparse_hamiltonian(
        basis,
        [
            OperatorTerm(
                "+-",
                [(0.5, site, site + 1) for site in 1:(L - 1)],
            ),
            OperatorTerm(
                "-+",
                [(0.5, site, site + 1) for site in 1:(L - 1)],
            ),
            OperatorTerm(
                "zz",
                [(0.8, site, site + 1) for site in 1:(L - 1)],
            ),
        ],
    )
    neel = sum(UInt64(1) << (site - 1) for site in 1:2:L)
    psi0 = zeros(ComplexF64, length(basis))
    psi0[state_index(basis, neel)] = 1
    energies, vectors, basis_vectors = lanczos_full(
        H,
        psi0,
        80;
        full_ortho=true,
    )
    evolved = expm_lanczos(
        energies,
        vectors,
        basis_vectors;
        a=-0.7im,
    )
    return (
        valid=abs(norm(evolved) - 1) < 2e-9,
        residual=abs(norm(evolved) - 1),
        value=evolved,
    )
end

validate_xxz_lanczos_quench(_, result) = result.valid

setup_floquet_heating() = nothing

function run_floquet_heating(_)
    L = 9
    basis = SpinBasis1D(L; pauli=true)
    Hzz = sparse_hamiltonian(
        basis,
        [
            OperatorTerm(
                "zz",
                [(0.9, site, mod1(site + 1, L)) for site in 1:L],
            ),
        ],
    )
    Hx = sparse_hamiltonian(
        basis,
        [OperatorTerm("x", [(0.73, site) for site in 1:L])],
    )
    floquet = Floquet(
        Dict(:H_list => [Hzz, Hx], :dt_list => [0.17, 0.23]);
        UF=true,
        thetaF=true,
    )
    identity_matrix = Matrix{ComplexF64}(
        I,
        length(basis),
        length(basis),
    )
    unitarity_error =
        norm(floquet.UF' * floquet.UF - identity_matrix) / length(basis)
    phase_error = maximum(abs.(abs.(floquet.thetaF) .- 1))
    return (
        valid=unitarity_error < 3e-11 && phase_error < 3e-10,
        residual=unitarity_error,
        value=floquet.thetaF,
    )
end

validate_floquet_heating(_, result) = result.valid

setup_spinful_hubbard() = nothing

function run_spinful_hubbard(_)
    L = 8
    interaction = 4.0
    bonds = [(site, site + 1) for site in 1:(L - 1)]
    basis = SpinfulFermionBasis1D(L; Nf=(4, 4))
    terms = [
        OperatorTerm(
            "+-|",
            [(-1.0, left, right) for (left, right) in bonds],
        ),
        OperatorTerm(
            "-+|",
            [(1.0, left, right) for (left, right) in bonds],
        ),
        OperatorTerm(
            "|+-",
            [(-1.0, left, right) for (left, right) in bonds],
        ),
        OperatorTerm(
            "|-+",
            [(1.0, left, right) for (left, right) in bonds],
        ),
        OperatorTerm(
            "n|n",
            [(interaction, site, site) for site in 1:L],
        ),
    ]
    H = sparse_hamiltonian(basis, terms)
    seed = uniform_seed(length(basis))
    values, vectors = eigsh(
        H;
        k=6,
        which=:SA,
        v0=seed,
        ncv=32,
        tol=1e-9,
        maxiter=5_000,
    )
    residual = eigensystem_residual(H, values, vectors)
    return (
        valid=residual < 2e-7,
        residual=residual,
        value=values,
    )
end

validate_spinful_hubbard(_, result) = result.valid

setup_interacting_ssh() = nothing

function run_interacting_ssh(_)
    L = 16
    particles = L ÷ 2
    bonds = [(site, site + 1) for site in 1:(L - 1)]
    hopping = [isodd(site) ? 0.6 : 1.0 for site in 1:(L - 1)]
    basis = SpinlessFermionBasis1D(L; Nf=particles)
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
    H = sparse_hamiltonian(basis, terms)
    seed = uniform_seed(length(basis))
    values, vectors = eigsh(
        H;
        k=6,
        which=:SA,
        v0=seed,
        ncv=32,
        tol=1e-9,
        maxiter=5_000,
    )
    residual = eigensystem_residual(H, values, vectors)
    return (
        valid=residual < 2e-7,
        residual=residual,
        value=values,
    )
end

validate_interacting_ssh(_, result) = result.valid

setup_translation_xxz() = nothing

function run_translation_xxz(_)
    L = 18
    basis = SpinBasis1D(
        L;
        nup=9,
        pauli=false,
        kblock=0,
    )
    H = sparse_hamiltonian(
        basis,
        xxz_terms(L; Jxy=1.0, Jz=0.9),
    )
    seed = uniform_seed(length(basis))
    values, vectors = eigsh(
        H;
        k=4,
        which=:SA,
        v0=seed,
        ncv=24,
        tol=1e-9,
        maxiter=5_000,
    )
    residual = eigensystem_residual(H, values, vectors)
    return (
        valid=residual < 2e-7,
        residual=residual,
        value=values,
    )
end

validate_translation_xxz(_, result) = result.valid

function make_cases()
    return PaperCase[
        PaperCase(
            "paper_mbl_shift_invert_l14",
            "19021",
            "MBL mid-spectrum shift-invert",
            "L=14;Nup=7;dimension=3432;k=6;sigma=0;ncv=32",
            setup_mbl_shift_invert,
            run_mbl_shift_invert,
            validate_mbl_shift_invert,
        ),
        PaperCase(
            "paper_xxz_lanczos_quench_l16",
            "3831",
            "XXZ Lanczos quench",
            "L=16;Nup=8;dimension=12870;m=80;t=0.7",
            setup_xxz_lanczos_quench,
            run_xxz_lanczos_quench,
            validate_xxz_lanczos_quench,
        ),
        PaperCase(
            "paper_floquet_heating_l9",
            "1069",
            "Floquet heating full unitary",
            "L=9;dimension=512;steps=2",
            setup_floquet_heating,
            run_floquet_heating,
            validate_floquet_heating,
        ),
        PaperCase(
            "paper_spinful_hubbard_l8",
            "3127",
            "Spinful Hubbard low-energy spectrum",
            "L=8;Nf=(4,4);dimension=4900;k=6",
            setup_spinful_hubbard,
            run_spinful_hubbard,
            validate_spinful_hubbard,
        ),
        PaperCase(
            "paper_interacting_ssh_l16",
            "19000",
            "Interacting SSH low-energy spectrum",
            "L=16;Nf=8;dimension=12870;k=6;t1=0.6;t2=1;V=2",
            setup_interacting_ssh,
            run_interacting_ssh,
            validate_interacting_ssh,
        ),
        PaperCase(
            "paper_translation_xxz_l18",
            "3831",
            "Translation-sector XXZ spectrum",
            "L=18;Nup=9;kblock=0;parent_dimension=48620;k=4",
            setup_translation_xxz,
            run_translation_xxz,
            validate_translation_xxz,
        ),
    ]
end

function percentile(sorted_values::Vector{Float64}, probability::Float64)
    position = probability * (length(sorted_values) - 1) + 1
    lower = floor(Int, position)
    upper = ceil(Int, position)
    lower == upper && return sorted_values[lower]
    weight = position - lower
    return sorted_values[lower] * (1 - weight) +
        sorted_values[upper] * weight
end

function median_value(values::AbstractVector{<:Real})
    ordered = sort(collect(values))
    midpoint = length(ordered) ÷ 2
    return isodd(length(ordered)) ?
        ordered[midpoint + 1] :
        (ordered[midpoint] + ordered[midpoint + 1]) / 2
end

function time_case(case::PaperCase, samples::Integer=DEFAULT_SAMPLES)
    context = case.setup()

    preflight = case.run(context)
    case.validate(context, preflight) ||
        error(
            "$(case.case_id) preflight failed: " *
            "residual=$(preflight.residual)",
        )
    SINK[] = preflight

    # Exactly one warmup follows the correctness preflight.
    SINK[] = case.run(context)

    timings = Float64[]
    allocated_bytes = Int[]
    for _ in 1:samples
        GC.gc()
        measurement = @timed begin
            SINK[] = case.run(context)
        end
        push!(timings, measurement.time)
        push!(allocated_bytes, measurement.bytes)
    end
    ordered = sort(timings)
    return (
        language="julia",
        suite="paper",
        case_id=case.case_id,
        family_id=case.family_id,
        benchmark=case.benchmark,
        category="workflow",
        comparison="end_to_end",
        storage="current_backend",
        supported="true",
        note="Full basis + Hamiltonian + solver/observable pipeline.",
        parameters=case.parameters,
        validation="passed",
        samples=samples,
        iterations_per_sample=ITERATIONS_PER_SAMPLE,
        median_seconds=median_value(timings),
        mean_seconds=sum(timings) / length(timings),
        stdev_seconds=length(timings) > 1 ?
            sqrt(
                sum(
                    (value - sum(timings) / length(timings))^2
                    for value in timings
                ) / (length(timings) - 1),
            ) :
            0.0,
        min_seconds=ordered[1],
        p05_seconds=percentile(ordered, 0.05),
        p25_seconds=percentile(ordered, 0.25),
        p75_seconds=percentile(ordered, 0.75),
        p95_seconds=percentile(ordered, 0.95),
        max_seconds=ordered[end],
        median_allocated_bytes=round(
            Int,
            median_value(allocated_bytes),
        ),
        runtime="Julia $(VERSION); QuSpin candidate",
        raw_samples_seconds=join(timings, ";"),
    )
end

function csv_cell(value)
    text = string(value)
    return occursin(r"[\",\n]", text) ?
        "\"" * replace(text, "\"" => "\"\"") * "\"" :
        text
end

function parse_cli(arguments)
    output = nothing
    selected = String[]
    samples = DEFAULT_SAMPLES
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        if argument == "--output"
            index < length(arguments) ||
                error("--output requires a path")
            output = arguments[index + 1]
            index += 2
        elseif argument == "--case"
            index < length(arguments) ||
                error("--case requires a case id")
            append!(
                selected,
                filter(!isempty, split(arguments[index + 1], ",")),
            )
            index += 2
        elseif argument == "--samples"
            index < length(arguments) ||
                error("--samples requires a positive integer")
            samples = parse(Int, arguments[index + 1])
            samples >= 1 || error("--samples must be positive")
            index += 2
        elseif argument in ("-h", "--help")
            println(
                "usage: julia_paper_workflow_benchmarks.jl ",
                "--output PATH [--case CASE_ID[,CASE_ID...]] ",
                "[--samples N]",
            )
            exit()
        else
            error("unknown argument: $argument")
        end
    end
    output === nothing && error("--output PATH is required")
    return String(output), selected, samples
end

function select_cases(cases, selected)
    isempty(selected) && return cases
    available = Set(case.case_id for case in cases)
    unknown = [case_id for case_id in selected if case_id ∉ available]
    isempty(unknown) ||
        error(
            "unknown case(s): $(join(unknown, ", ")); available: " *
            join(sort!(collect(available)), ", "),
        )
    requested = Set(selected)
    return [case for case in cases if case.case_id in requested]
end

function write_csv(output, rows)
    mkpath(dirname(abspath(output)))
    columns = propertynames(first(rows))
    open(output, "w") do io
        println(io, join(columns, ","))
        for row in rows
            println(
                io,
                join(
                    (
                        csv_cell(getproperty(row, column))
                        for column in columns
                    ),
                    ",",
                ),
            )
        end
    end
end

function main(arguments=ARGS)
    output, selected, samples = parse_cli(arguments)
    cases = select_cases(make_cases(), selected)
    isempty(cases) && error("no benchmark cases selected")
    rows = Any[]
    for case in cases
        println(stderr, "running $(case.case_id)")
        push!(rows, time_case(case, samples))
    end
    write_csv(output, rows)
    return rows
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__) &&
   get(ENV, "QUSPIN_BENCHMARK_LIBRARY_MODE", "") != "1"
    main()
end
