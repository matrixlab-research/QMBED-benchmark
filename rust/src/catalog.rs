use crate::Invariant;

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum Capability {
    SpinBasis,
    BosonBasis,
    SpinlessFermionBasis,
    SpinfulFermionBasis,
    UserBasis,
    SymmetrySector,
    CscAssembly,
    PartialEigensolver,
    ShiftInvert,
    KrylovEvolution,
    Floquet,
    SubspaceTracking,
    KrylovSpectrum,
    CrossSectorOperator,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum WorkflowOrigin {
    OriginalBaseline,
    LkmExtension,
}

#[derive(Clone, Debug, PartialEq)]
pub struct WorkflowCase {
    pub case_id: &'static str,
    pub family_id: &'static str,
    pub name: &'static str,
    pub parameters: &'static str,
    pub origin: WorkflowOrigin,
    pub capabilities: &'static [Capability],
    pub invariants: &'static [Invariant],
}

const MBL_CAPS: &[Capability] = &[
    Capability::SpinBasis,
    Capability::CscAssembly,
    Capability::PartialEigensolver,
    Capability::ShiftInvert,
];
const XXZ_QUENCH_CAPS: &[Capability] = &[
    Capability::SpinBasis,
    Capability::CscAssembly,
    Capability::KrylovEvolution,
];
const FLOQUET_CAPS: &[Capability] = &[
    Capability::SpinBasis,
    Capability::CscAssembly,
    Capability::Floquet,
];
const SPINFUL_SPECTRUM_CAPS: &[Capability] = &[
    Capability::SpinfulFermionBasis,
    Capability::CscAssembly,
    Capability::PartialEigensolver,
];
const SPINLESS_SPECTRUM_CAPS: &[Capability] = &[
    Capability::SpinlessFermionBasis,
    Capability::CscAssembly,
    Capability::PartialEigensolver,
];
const TRANSLATION_CAPS: &[Capability] = &[
    Capability::SpinBasis,
    Capability::SymmetrySector,
    Capability::CscAssembly,
    Capability::PartialEigensolver,
];
const FIDELITY_CAPS: &[Capability] = &[
    Capability::SpinBasis,
    Capability::CscAssembly,
    Capability::PartialEigensolver,
    Capability::SubspaceTracking,
];
const PXP_CAPS: &[Capability] = &[
    Capability::UserBasis,
    Capability::CscAssembly,
    Capability::KrylovEvolution,
];
const BOSE_QUENCH_CAPS: &[Capability] = &[
    Capability::BosonBasis,
    Capability::CscAssembly,
    Capability::KrylovEvolution,
];
const HUBBARD_CURRENT_CAPS: &[Capability] = &[
    Capability::SpinfulFermionBasis,
    Capability::CscAssembly,
    Capability::PartialEigensolver,
    Capability::KrylovEvolution,
];
const RESPONSE_CAPS: &[Capability] = &[
    Capability::SpinBasis,
    Capability::CscAssembly,
    Capability::PartialEigensolver,
    Capability::KrylovSpectrum,
];
const PARTICLE_ADDITION_CAPS: &[Capability] = &[
    Capability::SpinlessFermionBasis,
    Capability::CscAssembly,
    Capability::PartialEigensolver,
    Capability::CrossSectorOperator,
    Capability::KrylovSpectrum,
];

const RESIDUAL_2E7: &[Invariant] = &[Invariant::AtMost("residual", 2.0e-7)];
const XXZ_QUENCH_INVARIANTS: &[Invariant] = &[Invariant::AtMost("norm_error", 2.0e-9)];
const FLOQUET_INVARIANTS: &[Invariant] = &[
    Invariant::AtMost("unitarity_error", 3.0e-11),
    Invariant::AtMost("phase_modulus_error", 3.0e-10),
];
const TFIM_INVARIANTS: &[Invariant] = &[
    Invariant::AtMost("residual", 3.0e-7),
    Invariant::Between("minimum_fidelity", 0.0, 1.0 + 1.0e-12),
    Invariant::Between("minimum_fidelity_interval", 2.0, 3.0),
];
const PXP_INVARIANTS: &[Invariant] = &[
    Invariant::Equal("basis_dimension", 103_682.0),
    Invariant::AtMost("norm_error", 5.0e-8),
    Invariant::AtLeast("revival_gain", 0.0),
];
const BOSE_INVARIANTS: &[Invariant] = &[
    Invariant::AtMost("norm_error", 5.0e-8),
    Invariant::AtMost("minimum_return_after_t0", 0.99),
];
const CURRENT_INVARIANTS: &[Invariant] = &[
    Invariant::AtMost("residual", 3.0e-7),
    Invariant::AtLeast("maximum_absolute_current", 1.0e-3),
];
const CONB_INVARIANTS: &[Invariant] = &[
    Invariant::AtMost("residual", 3.0e-7),
    Invariant::AtLeast("minimum_spectrum", -1.0e-12),
    Invariant::AtLeast("maximum_spectrum", 1.0e-3),
    Invariant::Equal("finite_fraction", 1.0),
];
const PARTICLE_ADDITION_INVARIANTS: &[Invariant] = &[
    Invariant::AtMost("residual", 3.0e-7),
    Invariant::AtLeast("transition_norm", 1.0e-6),
    Invariant::AtLeast("maximum_spectrum", 1.0e-4),
    Invariant::Equal("finite_fraction", 1.0),
];

/// The same twelve medium-size paper workflows used by the Python/Julia suite.
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn paper_workflows() -> Vec<WorkflowCase> {
    vec![
        WorkflowCase {
            case_id: "paper_mbl_shift_invert_l14",
            family_id: "19021",
            name: "MBL mid-spectrum shift-invert",
            parameters: "L=14;Nup=7;dimension=3432;k=6;sigma=0;ncv=32",
            origin: WorkflowOrigin::OriginalBaseline,
            capabilities: MBL_CAPS,
            invariants: RESIDUAL_2E7,
        },
        WorkflowCase {
            case_id: "paper_xxz_lanczos_quench_l16",
            family_id: "3831",
            name: "XXZ Lanczos quench",
            parameters: "L=16;Nup=8;dimension=12870;m=80;t=0.7",
            origin: WorkflowOrigin::OriginalBaseline,
            capabilities: XXZ_QUENCH_CAPS,
            invariants: XXZ_QUENCH_INVARIANTS,
        },
        WorkflowCase {
            case_id: "paper_floquet_heating_l9",
            family_id: "1069",
            name: "Floquet heating full unitary",
            parameters: "L=9;dimension=512;steps=2",
            origin: WorkflowOrigin::OriginalBaseline,
            capabilities: FLOQUET_CAPS,
            invariants: FLOQUET_INVARIANTS,
        },
        WorkflowCase {
            case_id: "paper_spinful_hubbard_l8",
            family_id: "3127",
            name: "Spinful Hubbard low-energy spectrum",
            parameters: "L=8;Nf=(4,4);dimension=4900;k=6",
            origin: WorkflowOrigin::OriginalBaseline,
            capabilities: SPINFUL_SPECTRUM_CAPS,
            invariants: RESIDUAL_2E7,
        },
        WorkflowCase {
            case_id: "paper_interacting_ssh_l16",
            family_id: "19000",
            name: "Interacting SSH low-energy spectrum",
            parameters: "L=16;Nf=8;dimension=12870;k=6;t1=0.6;t2=1;V=2",
            origin: WorkflowOrigin::OriginalBaseline,
            capabilities: SPINLESS_SPECTRUM_CAPS,
            invariants: RESIDUAL_2E7,
        },
        WorkflowCase {
            case_id: "paper_translation_xxz_l18",
            family_id: "3831",
            name: "Translation-sector XXZ spectrum",
            parameters: "L=18;Nup=9;kblock=0;parent_dimension=48620;k=4",
            origin: WorkflowOrigin::OriginalBaseline,
            capabilities: TRANSLATION_CAPS,
            invariants: RESIDUAL_2E7,
        },
        WorkflowCase {
            case_id: "paper_tfim_fidelity_l16",
            family_id: "30690",
            name: "TFIM degenerate-subspace fidelity scan",
            parameters: "L=16;dimension=65536;fields=0.8:0.1:1.2;k=2",
            origin: WorkflowOrigin::LkmExtension,
            capabilities: FIDELITY_CAPS,
            invariants: TFIM_INVARIANTS,
        },
        WorkflowCase {
            case_id: "paper_pxp_revival_l24",
            family_id: "PXP",
            name: "PXP constrained-state revival",
            parameters: "L=24;periodic=true;dimension=103682;m=100;times=5",
            origin: WorkflowOrigin::LkmExtension,
            capabilities: PXP_CAPS,
            invariants: PXP_INVARIANTS,
        },
        WorkflowCase {
            case_id: "paper_bose_hubbard_quench_l11",
            family_id: "BHM",
            name: "Bose-Hubbard Mott quench",
            parameters: "L=11;Nb=11;sps=3;dimension=25653;J/U=0.1;m=100;times=5",
            origin: WorkflowOrigin::LkmExtension,
            capabilities: BOSE_QUENCH_CAPS,
            invariants: BOSE_INVARIANTS,
        },
        WorkflowCase {
            case_id: "paper_hubbard_current_l10",
            family_id: "32600",
            name: "Spinful-Hubbard current quench",
            parameters: "L=10;Nf=(5,5);dimension=63504;U/t=8;m=100;times=5",
            origin: WorkflowOrigin::LkmExtension,
            capabilities: HUBBARD_CURRENT_CAPS,
            invariants: CURRENT_INVARIANTS,
        },
        WorkflowCase {
            case_id: "paper_conb_dsf_l16",
            family_id: "CoNb2O6",
            name: "CoNb2O6 dynamical structure factor",
            parameters: "L=16;dimension=65536;B=7T;m=100;frequencies=81",
            origin: WorkflowOrigin::LkmExtension,
            capabilities: RESPONSE_CAPS,
            invariants: CONB_INVARIANTS,
        },
        WorkflowCase {
            case_id: "paper_particle_addition_6x3",
            family_id: "FQAH",
            name: "Particle-addition spectrum (6x3 proxy)",
            parameters: "Lx=6;Ly=3;Nf=6->7;dimensions=18564->31824;m=100;frequencies=81",
            origin: WorkflowOrigin::LkmExtension,
            capabilities: PARTICLE_ADDITION_CAPS,
            invariants: PARTICLE_ADDITION_INVARIANTS,
        },
    ]
}
