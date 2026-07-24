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

function deterministic_seed(dimension::Integer)
    values = sin.(Float64.(1:dimension))
    return values ./ norm(values)
end

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

setup_tfim_fidelity_scan() = nothing

function run_tfim_fidelity_scan(_)
    L = 16
    basis = SpinBasis1D(L; pauli=true)
    bonds = [(site, mod1(site + 1, L)) for site in 1:L]
    fields = (0.8, 0.9, 1.0, 1.1, 1.2)
    spaces = Matrix{ComplexF64}[]
    residuals = Float64[]
    for field in fields
        H = sparse_hamiltonian(
            basis,
            [
                OperatorTerm(
                    "zz",
                    [(-1.0, left, right) for (left, right) in bonds],
                ),
                OperatorTerm("x", [(-field, site) for site in 1:L]),
            ],
        )
        values, vectors = eigsh(
            H;
            k=2,
            which=:SA,
            v0=deterministic_seed(length(basis)),
            ncv=20,
            tol=1e-9,
            maxiter=5_000,
        )
        order = sortperm(values)
        push!(spaces, ComplexF64.(vectors[:, order]))
        push!(residuals, eigensystem_residual(H, values, vectors))
    end
    tracked = track_eigenspaces(spaces)
    residual = maximum(residuals)
    return (
        valid=residual < 3e-7 &&
            all(isfinite, tracked.fidelities) &&
            all(0 .< tracked.fidelities .<= 1 + 1e-12) &&
            argmin(tracked.fidelities) in 2:3,
        residual=residual,
        value=tracked.fidelities,
    )
end

validate_tfim_fidelity_scan(_, result) = result.valid

setup_pxp_revival() = nothing

function run_pxp_revival(_)
    L = 24
    constrained = constraint_states(
        L;
        prefix_allowed=(occupations, site) ->
            site == 1 || occupations[site - 1] + occupations[site] <= 1,
        state_allowed=occupations -> occupations[1] + occupations[end] <= 1,
    )
    basis = UserBasis(
        UInt64,
        L,
        Dict('x' => ComplexF64[0 1; 1 0]);
        states=constrained,
        allowed_ops=('x',),
    )
    H = sparse_hamiltonian(
        basis,
        [OperatorTerm("x", [(1.0, site) for site in 1:L])],
    )
    neel = sum(UInt64(1) << (site - 1) for site in 1:2:L)
    initial = zeros(ComplexF64, length(basis))
    initial[state_index(basis, neel)] = 1
    energies, vectors, basis_vectors = lanczos_full(
        H,
        initial,
        100;
        full_ortho=true,
    )
    times = (0.0, 2.4, 4.8, 7.2, 9.6)
    evolved = [
        expm_lanczos(energies, vectors, basis_vectors; a=-im * time)
        for time in times
    ]
    norm_error = maximum(abs(norm(state) - 1) for state in evolved)
    fidelities = [abs2(dot(initial, state)) for state in evolved]
    return (
        valid=length(basis) == 103_682 &&
            norm_error < 5e-8 &&
            fidelities[3] > fidelities[2],
        residual=norm_error,
        value=fidelities,
    )
end

validate_pxp_revival(_, result) = result.valid

setup_bose_hubbard_mott_quench() = nothing

function run_bose_hubbard_mott_quench(_)
    L = 11
    hopping = 0.1
    interaction = 1.0
    basis = BosonBasis1D(L; Nb=L, sps=3)
    bonds = [(site, site + 1) for site in 1:(L - 1)]
    H = sparse_hamiltonian(
        basis,
        [
            OperatorTerm(
                "+-",
                [(-hopping, left, right) for (left, right) in bonds],
            ),
            OperatorTerm(
                "-+",
                [(-hopping, left, right) for (left, right) in bonds],
            ),
            OperatorTerm(
                "nn",
                [(0.5 * interaction, site, site) for site in 1:L],
            ),
            OperatorTerm("n", [(-0.5 * interaction, site) for site in 1:L]),
        ],
    )
    mott = sum(UInt64(3)^(site - 1) for site in 1:L)
    initial = zeros(ComplexF64, length(basis))
    initial[state_index(basis, mott)] = 1
    energies, vectors, basis_vectors = lanczos_full(
        H,
        initial,
        100;
        full_ortho=true,
    )
    times = (0.0, 25.0, 50.0, 100.0, 200.0)
    evolved = [
        expm_lanczos(energies, vectors, basis_vectors; a=-im * time)
        for time in times
    ]
    norm_error = maximum(abs(norm(state) - 1) for state in evolved)
    returns = [abs2(dot(initial, state)) for state in evolved]
    return (
        valid=norm_error < 5e-8 && minimum(returns[2:end]) < 0.99,
        residual=norm_error,
        value=returns,
    )
end

validate_bose_hubbard_mott_quench(_, result) = result.valid

setup_spinful_hubbard_current_quench() = nothing

function run_spinful_hubbard_current_quench(_)
    L = 10
    basis = SpinfulFermionBasis1D(L; Nf=(5, 5))
    bonds = [(site, site + 1) for site in 1:(L - 1)]
    kinetic = OperatorTerm[
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
    ]
    interaction = OperatorTerm(
        "n|n",
        [(8.0, site, site) for site in 1:L],
    )
    unbiased = sparse_hamiltonian(basis, vcat(kinetic, [interaction]))
    bias = OperatorTerm[
        OperatorTerm(
            "n|",
            [(site <= L ÷ 2 ? -1.5 : 1.5, site) for site in 1:L],
        ),
        OperatorTerm(
            "|n",
            [(site <= L ÷ 2 ? -1.5 : 1.5, site) for site in 1:L],
        ),
    ]
    biased = sparse_hamiltonian(
        basis,
        vcat(kinetic, [interaction], bias),
    )
    values, vectors = eigsh(
        biased;
        k=1,
        which=:SA,
        v0=deterministic_seed(length(basis)),
        ncv=20,
        tol=1e-9,
        maxiter=5_000,
    )
    initial = vectors[:, 1]
    ground_residual = norm(biased * initial - values[1] * initial)
    center = L ÷ 2
    forward =
        operator_matrix(
            basis,
            "+-|",
            [(1.0, center, center + 1)];
            sparse=true,
        ) +
        operator_matrix(
            basis,
            "|+-",
            [(1.0, center, center + 1)];
            sparse=true,
        )
    current = -im .* (forward - forward')
    energies, lanczos_vectors, basis_vectors = lanczos_full(
        unbiased,
        initial,
        100;
        full_ortho=true,
    )
    evolved = [
        expm_lanczos(
            energies,
            lanczos_vectors,
            basis_vectors;
            a=-im * time,
        )
        for time in (0.0, 0.5, 1.0, 1.5, 2.0)
    ]
    norm_error = maximum(abs(norm(state) - 1) for state in evolved)
    currents = [real(dot(state, current * state)) for state in evolved]
    residual = max(ground_residual, norm_error)
    return (
        valid=residual < 3e-7 && maximum(abs.(currents)) > 1e-3,
        residual=residual,
        value=currents,
    )
end

validate_spinful_hubbard_current_quench(_, result) = result.valid

setup_conb_dynamical_structure_factor() = nothing

function run_conb_dynamical_structure_factor(_)
    L = 16
    basis = SpinBasis1D(L; pauli=false)
    nearest = [(site, mod1(site + 1, L)) for site in 1:L]
    next_nearest = [(site, mod1(site + 2, L)) for site in 1:L]
    transverse_field = 3.21 * 0.0578838 * 7.0 / 2.88
    H = sparse_hamiltonian(
        basis,
        [
            OperatorTerm(
                "zz",
                [(-1.0, left, right) for (left, right) in nearest],
            ),
            OperatorTerm(
                "xx",
                [(-0.205, left, right) for (left, right) in nearest],
            ),
            OperatorTerm(
                "yy",
                [(-0.205, left, right) for (left, right) in nearest],
            ),
            OperatorTerm(
                "zz",
                [(0.135, left, right) for (left, right) in next_nearest],
            ),
            OperatorTerm(
                "xx",
                [(0.003, left, right) for (left, right) in next_nearest],
            ),
            OperatorTerm(
                "yy",
                [(0.003, left, right) for (left, right) in next_nearest],
            ),
            OperatorTerm(
                "x",
                [(-transverse_field, site) for site in 1:L],
            ),
        ],
    )
    values, vectors = eigsh(
        H;
        k=1,
        which=:SA,
        v0=ComplexF64.(deterministic_seed(length(basis))),
        ncv=20,
        tol=1e-9,
        maxiter=5_000,
    )
    ground = vectors[:, 1]
    spin_q = operator_matrix(
        basis,
        "z",
        [((-1.0)^(site - 1), site) for site in 1:L];
        sparse=true,
    )
    frequencies = collect(range(0.0, 4.0; length=81))
    spectrum = spectral_function(
        H,
        ground,
        spin_q,
        frequencies;
        reference_energy=values[1],
        broadening=0.05,
        method=:krylov,
        krylov_dim=100,
    )
    residual = norm(H * ground - values[1] * ground)
    return (
        valid=residual < 3e-7 &&
            all(isfinite, spectrum) &&
            minimum(spectrum) >= -1e-12 &&
            maximum(spectrum) > 1e-3,
        residual=residual,
        value=spectrum,
    )
end

validate_conb_dynamical_structure_factor(_, result) = result.valid

setup_particle_addition_spectrum() = nothing

function triangular_bonds(length_x::Integer, length_y::Integer)
    bonds = Set{Tuple{Int,Int}}()
    for y in 0:(length_y - 1), x in 0:(length_x - 1)
        site = x + length_x * y + 1
        for (next_x, next_y) in (
            (mod(x + 1, length_x), y),
            (x, mod(y + 1, length_y)),
            (mod(x + 1, length_x), mod(y + 1, length_y)),
        )
            neighbor = next_x + length_x * next_y + 1
            push!(bonds, minmax(site, neighbor))
        end
    end
    return sort!(collect(bonds))
end

function run_particle_addition_spectrum(_)
    length_x, length_y = 6, 3
    L = length_x * length_y
    source_basis = SpinlessFermionBasisGeneral(L; Nf=6)
    target_basis = SpinlessFermionBasisGeneral(L; Nf=7)
    bonds = triangular_bonds(length_x, length_y)
    function interacting_model(basis)
        return sparse_hamiltonian(
            basis,
            [
                OperatorTerm(
                    "+-",
                    [(-1.0, left, right) for (left, right) in bonds],
                ),
                OperatorTerm(
                    "-+",
                    [(1.0, left, right) for (left, right) in bonds],
                ),
                OperatorTerm(
                    "nn",
                    [(2.0, left, right) for (left, right) in bonds],
                ),
            ],
        )
    end
    source_hamiltonian = interacting_model(source_basis)
    target_hamiltonian = interacting_model(target_basis)
    values, vectors = eigsh(
        source_hamiltonian;
        k=1,
        which=:SA,
        v0=deterministic_seed(length(source_basis)),
        ncv=20,
        tol=1e-9,
        maxiter=5_000,
    )
    source_state = vectors[:, 1]
    transition = op_shift_sector(
        target_basis,
        source_basis,
        [("+", [L ÷ 2 + 1], 1.0)],
        source_state,
    )
    identity_target = sparse(I, length(target_basis), length(target_basis))
    frequencies = collect(range(-4.0, 12.0; length=81))
    spectrum = spectral_function(
        target_hamiltonian,
        transition,
        identity_target,
        frequencies;
        reference_energy=values[1],
        broadening=0.1,
        method=:krylov,
        krylov_dim=100,
    )
    residual = norm(source_hamiltonian * source_state - values[1] * source_state)
    return (
        valid=residual < 3e-7 &&
            norm(transition) > 1e-6 &&
            all(isfinite, spectrum) &&
            maximum(spectrum) > 1e-4,
        residual=residual,
        value=spectrum,
    )
end

validate_particle_addition_spectrum(_, result) = result.valid

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
        PaperCase(
            "paper_tfim_fidelity_l16",
            "30690",
            "TFIM degenerate-subspace fidelity scan",
            "L=16;dimension=65536;fields=0.8:0.1:1.2;k=2",
            setup_tfim_fidelity_scan,
            run_tfim_fidelity_scan,
            validate_tfim_fidelity_scan,
        ),
        PaperCase(
            "paper_pxp_revival_l24",
            "PXP",
            "PXP constrained-state revival",
            "L=24;periodic=true;dimension=103682;m=100;times=5",
            setup_pxp_revival,
            run_pxp_revival,
            validate_pxp_revival,
        ),
        PaperCase(
            "paper_bose_hubbard_quench_l11",
            "BHM",
            "Bose-Hubbard Mott quench",
            "L=11;Nb=11;sps=3;dimension=25653;J/U=0.1;m=100;times=5",
            setup_bose_hubbard_mott_quench,
            run_bose_hubbard_mott_quench,
            validate_bose_hubbard_mott_quench,
        ),
        PaperCase(
            "paper_hubbard_current_l10",
            "32600",
            "Spinful-Hubbard current quench",
            "L=10;Nf=(5,5);dimension=63504;U/t=8;m=100;times=5",
            setup_spinful_hubbard_current_quench,
            run_spinful_hubbard_current_quench,
            validate_spinful_hubbard_current_quench,
        ),
        PaperCase(
            "paper_conb_dsf_l16",
            "CoNb2O6",
            "CoNb2O6 dynamical structure factor",
            "L=16;dimension=65536;B=7T;m=100;frequencies=81",
            setup_conb_dynamical_structure_factor,
            run_conb_dynamical_structure_factor,
            validate_conb_dynamical_structure_factor,
        ),
        PaperCase(
            "paper_particle_addition_6x3",
            "FQAH",
            "Particle-addition spectrum (6x3 proxy)",
            "Lx=6;Ly=3;Nf=6->7;dimensions=18564->31824;m=100;frequencies=81",
            setup_particle_addition_spectrum,
            run_particle_addition_spectrum,
            validate_particle_addition_spectrum,
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
   get(ENV, "QMBED_BENCHMARK_LIBRARY_MODE", "") != "1"
    main()
end
