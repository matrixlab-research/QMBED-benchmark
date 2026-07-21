#!/usr/bin/env julia

using LinearAlgebra
using QuSpin
using SparseArrays
using Statistics

BLAS.set_num_threads(1)

const SUITE = "small_ed_workflows"
const SAMPLES = 9
const WARMUPS = 3
const TARGET_SAMPLE_SECONDS = 0.12
const MAX_ITERATIONS = 50_000
const SINK = Ref{Any}(nothing)

struct WorkflowCase
    case_id::String
    family_id::String
    benchmark::String
    parameters::String
    setup::Any
    run::Any
    validate::Any
end

function spin_exchange_terms(
    L;
    J1=1.0,
    J2=0.0,
    delta=1.0,
    periodic=false,
)
    nearest = [(site, site + 1) for site in 1:(L - 1)]
    periodic && push!(nearest, (L, 1))
    next_nearest = [(site, site + 2) for site in 1:(L - 2)]
    periodic && append!(next_nearest, [(L - 1, 1), (L, 2)])
    xy = [(J1 / 2, left, right) for (left, right) in nearest]
    append!(
        xy,
        [(J2 / 2, left, right) for (left, right) in next_nearest],
    )
    zz = [(delta * J1, left, right) for (left, right) in nearest]
    append!(
        zz,
        [(delta * J2, left, right) for (left, right) in next_nearest],
    )
    return [
        OperatorTerm("+-", xy),
        OperatorTerm("-+", xy),
        OperatorTerm("zz", zz),
    ]
end

function fermion_hopping_terms(bonds)
    return [
        OperatorTerm(
            "+-",
            [(-amplitude, left, right) for (amplitude, left, right) in bonds],
        ),
        OperatorTerm(
            "-+",
            [(amplitude, left, right) for (amplitude, left, right) in bonds],
        ),
    ]
end

sorted_spectrum(H) = sort(real.(eigvals(Hermitian(Matrix(H)))))

function flux_ring(L, particles, flux; hopping=1.0, interaction=0.0)
    basis = SpinlessFermionBasis1D(L; Nf=particles)
    phase = cis(flux / L)
    plus = Tuple{ComplexF64,Int,Int}[]
    minus = Tuple{ComplexF64,Int,Int}[]
    for site in 1:L
        neighbor = mod1(site + 1, L)
        push!(plus, (-hopping * phase, site, neighbor))
        push!(minus, (hopping * conj(phase), site, neighbor))
    end
    terms = OperatorTerm[
        OperatorTerm("+-", plus),
        OperatorTerm("-+", minus),
    ]
    if !iszero(interaction)
        push!(
            terms,
            OperatorTerm(
                "nn",
                [
                    (interaction, site, mod1(site + 1, L))
                    for site in 1:L
                ],
            ),
        )
    end
    return basis, Hamiltonian(
        basis,
        terms;
        static_fmt=:csc,
        check_herm=false,
    )
end

function make_cases()
    cases = WorkflowCase[]

    push!(
        cases,
        WorkflowCase(
            "01",
            "3831",
            "heisenberg_dimer_singlet_triplet",
            "L=2;full_basis=true;storage=csc",
            () -> nothing,
            _ -> begin
                basis = SpinBasis1D(2; pauli=false)
                H = Hamiltonian(
                    basis,
                    spin_exchange_terms(2);
                    static_fmt=:csc,
                )
                sorted_spectrum(H)
            end,
            values -> isapprox(
                values,
                [-0.75, 0.25, 0.25, 0.25];
                atol=2e-14,
                rtol=0,
            ),
        ),
    )

    push!(
        cases,
        WorkflowCase(
            "02",
            "3831",
            "symmetry_resolved_xxz_low_energy",
            "L=6;nup=3;kblock=0;delta=0.7;periodic=true;k=2",
            () -> nothing,
            _ -> begin
                L = 6
                terms =
                    spin_exchange_terms(L; delta=0.7, periodic=true)
                full = Hamiltonian(
                    SpinBasis1D(L; pauli=false),
                    terms;
                    static_fmt=:csc,
                )
                sector_basis = SpinBasis1D(
                    L;
                    nup=L ÷ 2,
                    pauli=false,
                    kblock=0,
                )
                sector =
                    Hamiltonian(sector_basis, terms; static_fmt=:csc)
                sector_values = sorted_spectrum(sector)
                full_values = sorted_spectrum(full)
                seed = normalize(collect(1.0:length(sector_basis)))
                selected = first(
                    eigsh(
                        sector;
                        k=2,
                        which=:SA,
                        v0=seed,
                        tol=1e-11,
                        maxiter=1_000,
                    ),
                )
                (
                    full_values=full_values,
                    sector_values=sector_values,
                    selected=selected,
                    full_dimension=length(full.basis),
                    sector_dimension=length(sector_basis),
                )
            end,
            result ->
                result.sector_dimension < result.full_dimension &&
                all(
                    minimum(abs.(result.full_values .- value)) < 2e-12
                    for value in result.sector_values
                ) &&
                isapprox(
                    sort(result.selected),
                    result.sector_values[1:2];
                    atol=2e-10,
                    rtol=0,
                ),
        ),
    )

    push!(
        cases,
        WorkflowCase(
            "03",
            "23806",
            "frustrated_j1_j2_spin_gap",
            "L=8;J1=1;J2=0.5;periodic=true;nup=4,5",
            () -> nothing,
            _ -> begin
                L = 8
                terms = spin_exchange_terms(
                    L;
                    J1=1.0,
                    J2=0.5,
                    periodic=true,
                )
                singlet = Hamiltonian(
                    SpinBasis1D(L; nup=L ÷ 2, pauli=false),
                    terms;
                    static_fmt=:csc,
                )
                triplet = Hamiltonian(
                    SpinBasis1D(L; nup=L ÷ 2 + 1, pauli=false),
                    terms;
                    static_fmt=:csc,
                )
                E0 = minimum(sorted_spectrum(singlet))
                E1 = minimum(sorted_spectrum(triplet))
                (E0=E0, E1=E1, gap=E1 - E0)
            end,
            result -> 0.1 < result.gap < 1.0,
        ),
    )

    push!(
        cases,
        WorkflowCase(
            "04",
            "23806",
            "transverse_field_ising_phase_scan",
            "L=7;periodic=true;g=0.35,1.0,2.0",
            () -> nothing,
            _ -> begin
                L = 7
                basis = SpinBasis1D(L)
                bonds =
                    [(1.0, site, mod1(site + 1, L)) for site in 1:L]
                ising = function (g)
                    Hamiltonian(
                        basis,
                        [
                            OperatorTerm(
                                "zz",
                                [
                                    (-value, left, right)
                                    for (value, left, right) in bonds
                                ],
                            ),
                            OperatorTerm(
                                "x",
                                [(-g, site) for site in 1:L],
                            ),
                        ];
                        static_fmt=:csc,
                    )
                end
                ordered = sorted_spectrum(ising(0.35))
                polarized = sorted_spectrum(ising(2.0))
                (
                    ordered_gap=ordered[2] - ordered[1],
                    polarized_gap=polarized[2] - polarized[1],
                    hermitian=ishermitian(ising(1.0)),
                )
            end,
            result ->
                result.hermitian &&
                result.ordered_gap < result.polarized_gap,
        ),
    )

    push!(
        cases,
        WorkflowCase(
            "05",
            "11058",
            "heisenberg_dimer_susceptibility",
            "L=2;beta=1.3;storage=dense",
            () -> 1.3,
            beta -> begin
                basis = SpinBasis1D(2; pauli=false)
                H = Hamiltonian(
                    basis,
                    spin_exchange_terms(2);
                    static_fmt=:dense,
                )
                Mz = Matrix(
                    Hamiltonian(
                        basis,
                        [
                            OperatorTerm(
                                "z",
                                [(1.0, 1), (1.0, 2)],
                            ),
                        ],
                    ),
                )
                decomposition = eigen(Hermitian(Matrix(H)))
                weights = exp.(-beta .* decomposition.values)
                thermal_m2 = sum(
                    weights[index] *
                    real(
                        dot(
                            decomposition.vectors[:, index],
                            Mz * Mz * decomposition.vectors[:, index],
                        ),
                    )
                    for index in eachindex(weights)
                ) / sum(weights)
                expected = 2beta * exp(-beta / 4) /
                    (exp(3beta / 4) + 3exp(-beta / 4))
                (susceptibility=beta * thermal_m2, expected=expected)
            end,
            result -> isapprox(
                result.susceptibility,
                result.expected;
                atol=3e-14,
                rtol=0,
            ),
        ),
    )

    push!(
        cases,
        WorkflowCase(
            "06",
            "12290",
            "variational_energy_certified_by_ed",
            "L=6;nup=3;J2=0.35;trial=neel",
            () -> nothing,
            _ -> begin
                L = 6
                basis = SpinBasis1D(L; nup=L ÷ 2, pauli=false)
                H = Hamiltonian(
                    basis,
                    spin_exchange_terms(L; J2=0.35);
                    static_fmt=:csc,
                )
                neel =
                    sum(UInt64(1) << (site - 1) for site in 1:2:L)
                trial = zeros(ComplexF64, length(basis))
                trial[state_index(basis, neel)] = 1
                variational = real(expt_value(H, trial))
                exact = minimum(sorted_spectrum(H))
                (variational=variational, exact=exact)
            end,
            result ->
                result.variational >= result.exact - 2e-13 &&
                result.variational - result.exact > 0.2,
        ),
    )

    push!(
        cases,
        WorkflowCase(
            "07",
            "1069",
            "step_driven_floquet_quasienergies",
            "L=4;dt=0.23,0.31;full_dimension=16",
            () -> nothing,
            _ -> begin
                basis = SpinBasis1D(4)
                Hzz = Hamiltonian(
                    basis,
                    [
                        OperatorTerm(
                            "zz",
                            [
                                (0.7, site, mod1(site + 1, 4))
                                for site in 1:4
                            ],
                        ),
                    ];
                    static_fmt=:csc,
                )
                Hx = Hamiltonian(
                    basis,
                    [OperatorTerm("x", [(0.45, site) for site in 1:4])];
                    static_fmt=:csc,
                )
                floquet = Floquet(
                    Dict(
                        :H_list => [Hzz, Hx],
                        :dt_list => [0.23, 0.31],
                    );
                    UF=true,
                    thetaF=true,
                )
                identity_matrix =
                    Matrix{ComplexF64}(I, length(basis), length(basis))
                (
                    unitarity=norm(
                        floquet.UF' * floquet.UF - identity_matrix,
                    ),
                    phase_error=maximum(
                        abs.(abs.(floquet.thetaF) .- 1),
                    ),
                )
            end,
            result ->
                result.unitarity < 2e-12 &&
                result.phase_error < 2e-12,
        ),
    )

    push!(
        cases,
        WorkflowCase(
            "08",
            "1069",
            "high_frequency_effective_hamiltonian",
            "dimension=2;dt=0.005;steps=X,Z",
            () -> (
                X=ComplexF64[0 1; 1 0],
                Z=ComplexF64[1 0; 0 -1],
                dt=0.005,
            ),
            fixture -> begin
                floquet = Floquet(
                    Dict(
                        :H_list => [fixture.X, fixture.Z],
                        :dt_list => [fixture.dt, fixture.dt],
                    );
                    HF=true,
                    UF=true,
                )
                (
                    effective_error=norm(
                        floquet.HF - (fixture.X + fixture.Z) / 2,
                    ),
                    unitary_error=norm(
                        floquet.UF -
                        exp(-im * fixture.dt .* fixture.Z) *
                        exp(-im * fixture.dt .* fixture.X),
                    ),
                )
            end,
            result ->
                result.effective_error < 0.01 &&
                result.unitary_error < 3e-15,
        ),
    )

    push!(
        cases,
        WorkflowCase(
            "09",
            "19021",
            "anderson_localization_ipr",
            "L=10;Nf=1;disorder=20*(-1)^site+0.31site",
            () -> nothing,
            _ -> begin
                L = 10
                basis = SpinlessFermionBasis1D(L; Nf=1)
                hopping = fermion_hopping_terms(
                    [(1.0, site, site + 1) for site in 1:(L - 1)],
                )
                clean =
                    Hamiltonian(basis, hopping; static_fmt=:csc)
                disordered = Hamiltonian(
                    basis,
                    vcat(
                        hopping,
                        [
                            OperatorTerm(
                                "n",
                                [
                                    (
                                        20.0 * (-1)^site + 0.31site,
                                        site,
                                    )
                                    for site in 1:L
                                ],
                            ),
                        ],
                    );
                    static_fmt=:csc,
                )
                clean_vectors =
                    eigen(Hermitian(Matrix(clean))).vectors
                localized_vectors =
                    eigen(Hermitian(Matrix(disordered))).vectors
                clean_ipr =
                    mean(sum(abs2.(clean_vectors) .^ 2; dims=1))
                localized_ipr =
                    mean(sum(abs2.(localized_vectors) .^ 2; dims=1))
                (clean_ipr=clean_ipr, localized_ipr=localized_ipr)
            end,
            result -> result.localized_ipr > 2result.clean_ipr,
        ),
    )

    push!(
        cases,
        WorkflowCase(
            "10",
            "19021",
            "mbl_mid_spectrum_level_statistics",
            "L=8;nup=4;periodic=true;k=4;sigma=mean(E)",
            () -> [
                2.1,
                -1.7,
                0.3,
                3.2,
                -2.6,
                0.8,
                1.4,
                -3.1,
            ],
            fields -> begin
                L = 8
                basis = SpinBasis1D(L; nup=L ÷ 2, pauli=false)
                terms = vcat(
                    spin_exchange_terms(L; periodic=true),
                    [
                        OperatorTerm(
                            "z",
                            [(fields[site], site) for site in 1:L],
                        ),
                    ],
                )
                H = Hamiltonian(basis, terms; static_fmt=:csc)
                levels = sorted_spectrum(H)
                spacings = diff(levels)
                ratios = [
                    min(spacings[index], spacings[index + 1]) /
                    max(spacings[index], spacings[index + 1])
                    for index in 1:(length(spacings) - 1)
                    if max(
                        spacings[index],
                        spacings[index + 1],
                    ) > 1e-10
                ]
                seed = normalize(collect(1.0:length(basis)))
                selected = first(
                    eigsh(
                        H;
                        k=4,
                        sigma=mean(levels),
                        which=:LM,
                        v0=seed,
                        tol=1e-10,
                        maxiter=1_500,
                    ),
                )
                (
                    mean_ratio=mean(ratios),
                    selected_distance=maximum(
                        minimum(abs.(levels .- value))
                        for value in selected
                    ),
                )
            end,
            result ->
                0 < result.mean_ratio < 1 &&
                result.selected_distance < 2e-9,
        ),
    )

    push!(
        cases,
        WorkflowCase(
            "11",
            "19000",
            "ssh_topological_edge_modes",
            "L=8;Nf=1;t1=0.25;t2=1.0;open=true",
            () -> nothing,
            _ -> begin
                L = 8
                basis = SpinlessFermionBasis1D(L; Nf=1)
                bonds = [
                    (isodd(site) ? 0.25 : 1.0, site, site + 1)
                    for site in 1:(L - 1)
                ]
                H = Hamiltonian(
                    basis,
                    fermion_hopping_terms(bonds);
                    static_fmt=:csc,
                )
                decomposition = eigen(Hermitian(Matrix(H)))
                edge_indices =
                    sortperm(abs.(decomposition.values))[1:2]
                edge_weight = mean(
                    sum(
                        abs2(decomposition.vectors[row, index]) *
                        (
                            basis.occupations[row, 1] +
                            basis.occupations[row, L]
                        )
                        for row in 1:length(basis)
                    )
                    for index in edge_indices
                )
                (
                    edge_energy=maximum(
                        abs.(decomposition.values[edge_indices]),
                    ),
                    edge_weight=edge_weight,
                )
            end,
            result ->
                result.edge_energy < 0.01 &&
                result.edge_weight > 0.8,
        ),
    )

    push!(
        cases,
        WorkflowCase(
            "12",
            "19000",
            "interacting_extended_ssh_gap",
            "L=6;Nf=3;t1=0.45;t2=1.0;V=2.0",
            () -> nothing,
            _ -> begin
                L = 6
                basis = SpinlessFermionBasis1D(L; Nf=3)
                bonds = [
                    (isodd(site) ? 0.45 : 1.0, site, site + 1)
                    for site in 1:(L - 1)
                ]
                kinetic = fermion_hopping_terms(bonds)
                free =
                    Hamiltonian(basis, kinetic; static_fmt=:csc)
                interacting = Hamiltonian(
                    basis,
                    vcat(
                        kinetic,
                        [
                            OperatorTerm(
                                "nn",
                                [
                                    (2.0, site, site + 1)
                                    for site in 1:(L - 1)
                                ],
                            ),
                        ],
                    );
                    static_fmt=:csc,
                )
                free_levels = sorted_spectrum(free)
                interacting_levels = sorted_spectrum(interacting)
                (
                    ground_shift=interacting_levels[1] - free_levels[1],
                    gap=interacting_levels[2] - interacting_levels[1],
                )
            end,
            result ->
                !iszero(result.ground_shift) &&
                result.gap > 1e-3,
        ),
    )

    push!(
        cases,
        WorkflowCase(
            "13",
            "16282",
            "ssh_chiral_spectral_pairing",
            "L=7;Nf=1;t1=0.6;t2=1.2;open=true",
            () -> nothing,
            _ -> begin
                L = 7
                basis = SpinlessFermionBasis1D(L; Nf=1)
                bonds = [
                    (isodd(site) ? 0.6 : 1.2, site, site + 1)
                    for site in 1:(L - 1)
                ]
                H = Hamiltonian(
                    basis,
                    fermion_hopping_terms(bonds);
                    static_fmt=:csc,
                )
                sorted_spectrum(H)
            end,
            levels ->
                isapprox(
                    levels,
                    -reverse(levels);
                    atol=3e-14,
                    rtol=0,
                ) &&
                minimum(abs.(levels)) < 2e-14,
        ),
    )

    push!(
        cases,
        WorkflowCase(
            "14",
            "16282",
            "flux_periodicity_topological_diagnostic",
            "L=7;Nf=1;flux=0,2pi",
            () -> nothing,
            _ -> begin
                _, zero_flux = flux_ring(7, 1, 0.0)
                _, flux_quantum = flux_ring(7, 1, 2π)
                (
                    zero=sorted_spectrum(zero_flux),
                    quantum=sorted_spectrum(flux_quantum),
                )
            end,
            result -> isapprox(
                result.zero,
                result.quantum;
                atol=3e-13,
                rtol=0,
            ),
        ),
    )

    push!(
        cases,
        WorkflowCase(
            "15",
            "2880",
            "few_electron_quantum_dot_shell_filling",
            "orbitals=4;Nf=2;density_interaction=true",
            () -> nothing,
            _ -> begin
                basis = SpinlessFermionBasis1D(4; Nf=2)
                orbital_energies = [0.0, 0.7, 1.9, 3.1]
                interactions = [
                    (0.2 + 0.1abs(left - right), left, right)
                    for left in 1:4 for right in (left + 1):4
                ]
                H = Hamiltonian(
                    basis,
                    [
                        OperatorTerm(
                            "n",
                            [
                                (orbital_energies[site], site)
                                for site in 1:4
                            ],
                        ),
                        OperatorTerm("nn", interactions),
                    ];
                    static_fmt=:csc,
                )
                (
                    ground=minimum(sorted_spectrum(H)),
                    nonzeros=count(
                        value -> !iszero(value),
                        nonzeros(H.data),
                    ),
                    dimension=length(basis),
                )
            end,
            result ->
                isapprox(result.ground, 1.0; atol=2e-14, rtol=0) &&
                result.nonzeros <= result.dimension,
        ),
    )

    push!(
        cases,
        WorkflowCase(
            "16",
            "2880",
            "two_electron_hubbard_dot_singlet",
            "sites=2;Nf=(1,1);t=1;U=3",
            () -> 3.0,
            U -> begin
                hopping = [
                    OperatorTerm("+-|", [(-1.0, 1, 2)]),
                    OperatorTerm("-+|", [(1.0, 1, 2)]),
                    OperatorTerm("|+-", [(-1.0, 1, 2)]),
                    OperatorTerm("|-+", [(1.0, 1, 2)]),
                ]
                basis = SpinfulFermionBasis1D(2; Nf=(1, 1))
                H = Hamiltonian(
                    basis,
                    vcat(
                        hopping,
                        [
                            OperatorTerm(
                                "n|n",
                                [(U, 1, 1), (U, 2, 2)],
                            ),
                        ],
                    );
                    static_fmt=:csc,
                )
                (
                    ground=minimum(sorted_spectrum(H)),
                    expected=(U - sqrt(U^2 + 16)) / 2,
                )
            end,
            result -> isapprox(
                result.ground,
                result.expected;
                atol=3e-14,
                rtol=0,
            ),
        ),
    )

    push!(
        cases,
        WorkflowCase(
            "17",
            "3127",
            "cu_o_cluster_orbital_hole_fractions",
            "sites=3;Nf=(1,1);hybridization=0.4;charge_transfer=4",
            () -> nothing,
            _ -> begin
                basis = SpinfulFermionBasis1D(3; Nf=(1, 1))
                bonds = [(1, 2), (1, 3)]
                terms = [
                    OperatorTerm(
                        "+-|",
                        [(-0.4, pair...) for pair in bonds],
                    ),
                    OperatorTerm(
                        "-+|",
                        [(0.4, pair...) for pair in bonds],
                    ),
                    OperatorTerm(
                        "|+-",
                        [(-0.4, pair...) for pair in bonds],
                    ),
                    OperatorTerm(
                        "|-+",
                        [(0.4, pair...) for pair in bonds],
                    ),
                    OperatorTerm("n|", [(4.0, 2), (4.0, 3)]),
                    OperatorTerm("|n", [(4.0, 2), (4.0, 3)]),
                    OperatorTerm(
                        "n|n",
                        [(0.5, site, site) for site in 1:3],
                    ),
                ]
                H = Hamiltonian(basis, terms; static_fmt=:csc)
                ground =
                    eigen(Hermitian(Matrix(H))).vectors[:, 1]
                occupations = [
                    real(
                        expt_value(
                            Hamiltonian(
                                basis,
                                [
                                    OperatorTerm(
                                        "n|",
                                        [(1.0, site)],
                                    ),
                                    OperatorTerm(
                                        "|n",
                                        [(1.0, site)],
                                    ),
                                ],
                            ),
                            ground,
                        ),
                    )
                    for site in 1:3
                ]
                occupations
            end,
            occupations ->
                isapprox(sum(occupations), 2; atol=2e-13, rtol=0) &&
                occupations[1] > 1.8,
        ),
    )

    push!(
        cases,
        WorkflowCase(
            "18",
            "3127",
            "hubbard_cluster_charge_gap",
            "L=5;Nf=1,2,3;t=0.35;V=3",
            () -> nothing,
            _ -> begin
                cluster_energy = function (particles)
                    basis = SpinlessFermionBasis1D(5; Nf=particles)
                    kinetic = fermion_hopping_terms(
                        [
                            (0.35, site, site + 1)
                            for site in 1:4
                        ],
                    )
                    H = Hamiltonian(
                        basis,
                        vcat(
                            kinetic,
                            [
                                OperatorTerm(
                                    "nn",
                                    [
                                        (3.0, site, site + 1)
                                        for site in 1:4
                                    ],
                                ),
                            ],
                        );
                        static_fmt=:csc,
                    )
                    minimum(sorted_spectrum(H))
                end
                E1 = cluster_energy(1)
                E2 = cluster_energy(2)
                E3 = cluster_energy(3)
                (
                    E1=E1,
                    E2=E2,
                    E3=E3,
                    gap=E3 + E1 - 2 * E2,
                )
            end,
            result -> result.gap > 0,
        ),
    )

    push!(
        cases,
        WorkflowCase(
            "19",
            "12454",
            "fqhe_flux_spectral_flow",
            "L=6;Nf=2;flux=0,2pi;interaction=4",
            () -> nothing,
            _ -> begin
                _, at_zero =
                    flux_ring(6, 2, 0.0; interaction=4.0)
                _, at_quantum =
                    flux_ring(6, 2, 2π; interaction=4.0)
                (
                    zero=sorted_spectrum(at_zero),
                    quantum=sorted_spectrum(at_quantum),
                )
            end,
            result ->
                isapprox(
                    result.zero,
                    result.quantum;
                    atol=3e-13,
                    rtol=0,
                ) &&
                result.zero[2] - result.zero[1] < 3e-13 &&
                result.zero[3] - result.zero[2] > 1e-4,
        ),
    )

    push!(
        cases,
        WorkflowCase(
            "20",
            "14428",
            "collective_emission_bright_dark_modes",
            "emitters=4;nup=1;decay_rate=0.6",
            () -> 0.6,
            decay_rate -> begin
                emitters = 4
                basis =
                    SpinBasis1D(emitters; nup=1, pauli=false)
                collective_decay = Hamiltonian(
                    basis,
                    [
                        OperatorTerm(
                            "+-",
                            [
                                (
                                    -im * decay_rate / 2,
                                    left,
                                    right,
                                )
                                for left in 1:emitters
                                for right in 1:emitters
                            ],
                        ),
                    ];
                    static_fmt=:csc,
                    check_herm=false,
                )
                eigvals(Matrix(collective_decay))
            end,
            values ->
                count(value -> abs(value) < 2e-13, values) == 3 &&
                minimum(abs.(values .+ 1.2im)) < 2e-13,
        ),
    )

    push!(
        cases,
        WorkflowCase(
            "21",
            "26341",
            "fci_many_body_entanglement",
            "L=6;Nf=3;flux=pi/3;interaction=5;subsystem=3",
            () -> nothing,
            _ -> begin
                basis, H =
                    flux_ring(6, 3, π / 3; interaction=5.0)
                ground =
                    eigen(Hermitian(Matrix(H))).vectors[:, 1]
                entropy = ent_entropy(
                    basis,
                    ground;
                    sub_sys_A=[1, 2, 3],
                    return_rdm=:A,
                )
                (
                    entropy=entropy["Sent_A"],
                    trace=real(tr(entropy["rdm_A"])),
                )
            end,
            result ->
                result.entropy > 0 &&
                isapprox(result.trace, 1; atol=3e-13, rtol=0),
        ),
    )

    push!(
        cases,
        WorkflowCase(
            "22",
            "PXP",
            "constrained_prefix_state_generation",
            "L=20;periodic_blockade=true;algorithm=prefix_pruning",
            () -> nothing,
            _ -> constraint_states(
                20;
                prefix_allowed=(occupations, site) ->
                    site == 1 ||
                    occupations[site - 1] + occupations[site] <= 1,
                state_allowed=occupations ->
                    occupations[1] + occupations[end] <= 1,
            ),
            states -> length(states) == 15_127 && issorted(states),
        ),
    )

    push!(
        cases,
        WorkflowCase(
            "23",
            "PXP",
            "constrained_full_space_filter",
            "L=20;periodic_blockade=true;algorithm=full_space_predicate",
            () -> Dict('x' => ComplexF64[0 1; 1 0]),
            op_dict -> UserBasis(
                UInt64,
                20,
                op_dict;
                pre_check_state=state -> all(
                    Int((state >> (site - 1)) & 1) +
                    Int((state >> mod(site, 20)) & 1) <= 1
                    for site in 1:20
                ),
                allowed_ops=('x',),
            ),
            basis -> length(basis) == 15_127,
        ),
    )

    return cases
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
    ordered = sort(Float64.(values))
    midpoint = length(ordered) ÷ 2
    return isodd(length(ordered)) ?
        ordered[midpoint + 1] :
        (ordered[midpoint] + ordered[midpoint + 1]) / 2
end

function benchmark_case(case::WorkflowCase, fixture)
    preflight = case.run(fixture)
    case.validate(preflight) ||
        error("correctness preflight returned false for case $(case.case_id)")
    SINK[] = preflight

    warm_times = Float64[]
    for _ in 1:WARMUPS
        started_ns = time_ns()
        SINK[] = case.run(fixture)
        push!(warm_times, (time_ns() - started_ns) / 1.0e9)
    end
    estimate = max(median_value(warm_times), 1.0e-9)
    iterations = clamp(
        ceil(Int, TARGET_SAMPLE_SECONDS / estimate),
        1,
        MAX_ITERATIONS,
    )
    GC.gc()
    calibration_started_ns = time_ns()
    for _ in 1:iterations
        SINK[] = case.run(fixture)
    end
    calibration_seconds =
        max((time_ns() - calibration_started_ns) / 1.0e9, 1.0e-9)
    iterations = clamp(
        round(
            Int,
            iterations * TARGET_SAMPLE_SECONDS / calibration_seconds,
        ),
        1,
        MAX_ITERATIONS,
    )

    timings = Float64[]
    allocated_bytes = Int[]
    for _ in 1:SAMPLES
        GC.gc()
        measurement = @timed begin
            for _ in 1:iterations
                SINK[] = case.run(fixture)
            end
        end
        push!(timings, measurement.time / iterations)
        push!(allocated_bytes, measurement.bytes ÷ iterations)
    end

    ordered = sort(timings)
    return (
        language="julia",
        suite=SUITE,
        case_id=case.case_id,
        family_id=case.family_id,
        benchmark=case.benchmark,
        category="workflow",
        comparison="coverage_only",
        storage="current_backend",
        supported="true",
        note="Small end-to-end ED coverage timing; no cross-language claim.",
        parameters=case.parameters,
        validation="passed",
        samples=length(timings),
        iterations_per_sample=iterations,
        median_seconds=median_value(timings),
        p25_seconds=percentile(ordered, 0.25),
        p75_seconds=percentile(ordered, 0.75),
        median_allocated_bytes=round(
            Int,
            median_value(allocated_bytes),
        ),
        raw_samples_seconds=join(timings, ";"),
    )
end

function failed_row(case::WorkflowCase, exception)
    return (
        language="julia",
        suite=SUITE,
        case_id=case.case_id,
        family_id=case.family_id,
        benchmark=case.benchmark,
        category="workflow",
        comparison="coverage_only",
        storage="current_backend",
        supported="true",
        note="Small end-to-end ED coverage timing; no cross-language claim.",
        parameters=case.parameters,
        validation="failed: " * sprint(showerror, exception),
        samples=0,
        iterations_per_sample=0,
        median_seconds="",
        p25_seconds="",
        p75_seconds="",
        median_allocated_bytes="",
        raw_samples_seconds="",
    )
end

function csv_cell(value)
    text = string(value)
    return occursin(r"[\",\n]", text) ?
        "\"" * replace(text, "\"" => "\"\"") * "\"" :
        text
end

function parse_output(args)
    output_index = findfirst(==("--output"), args)
    output_index === nothing &&
        error("usage: julia_ed_workflow_benchmarks.jl --output PATH")
    output_index < length(args) ||
        error("--output requires a path")
    return args[output_index + 1]
end

function write_row(io, row, columns)
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
    flush(io)
end

function main()
    output = parse_output(ARGS)
    mkpath(dirname(output))
    cases = make_cases()
    length(cases) == 23 ||
        error("expected 23 ED workflow cases, found $(length(cases))")

    columns = propertynames(failed_row(first(cases), ErrorException("")))
    open(output, "w") do io
        println(io, join(columns, ","))
        for case in cases
            fixture = nothing
            try
                fixture = case.setup()
                row = benchmark_case(case, fixture)
                write_row(io, row, columns)
            catch exception
                write_row(io, failed_row(case, exception), columns)
                rethrow()
            finally
                fixture = nothing
                SINK[] = nothing
                GC.gc()
            end
        end
    end
end

main()
