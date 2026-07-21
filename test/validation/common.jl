using Statistics

function _lkm_spin_exchange_terms(L; periodic=true, delta=1.0)
    bonds = [(site, site + 1) for site in 1:(L - 1)]
    periodic && push!(bonds, (L, 1))
    xy = [(0.5, left, right) for (left, right) in bonds]
    zz = [(delta, left, right) for (left, right) in bonds]
    return [
        OperatorTerm("+-", xy),
        OperatorTerm("-+", xy),
        OperatorTerm("zz", zz),
    ]
end

function _lkm_tfim(basis, field)
    L = basis.L
    return Hamiltonian(
        basis,
        [
            OperatorTerm(
                "zz",
                [(-1.0, site, mod1(site + 1, L)) for site in 1:L],
            ),
            OperatorTerm("x", [(-field, site) for site in 1:L]),
        ];
        static_fmt=:csc,
    )
end

function _lkm_fermion_hopping_terms(L; hopping=1.0, spinful=false)
    bonds = [(site, site + 1) for site in 1:(L - 1)]
    if spinful
        return [
            OperatorTerm("+-|", [(-hopping, bond...) for bond in bonds]),
            OperatorTerm("-+|", [(hopping, bond...) for bond in bonds]),
            OperatorTerm("|+-", [(-hopping, bond...) for bond in bonds]),
            OperatorTerm("|-+", [(hopping, bond...) for bond in bonds]),
        ]
    end
    return [
        OperatorTerm("+-", [(-hopping, bond...) for bond in bonds]),
        OperatorTerm("-+", [(hopping, bond...) for bond in bonds]),
    ]
end

function _lkm_ground_state(H)
    values, vectors = eigsh(H; k=2, which=:SA, tol=1e-11)
    index = argmin(values)
    return values[index], vectors[:, index]
end

function _lkm_trapezoid(values, spacing)
    length(values) >= 2 || return zero(eltype(values))
    return spacing * (sum(values) - (first(values) + last(values)) / 2)
end
