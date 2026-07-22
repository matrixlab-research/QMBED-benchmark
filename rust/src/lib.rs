//! Language-neutral verification contracts for a future Rust `QuSpin` package.
//!
//! This crate intentionally does not implement exact diagonalization. It owns
//! the private test-side adapter boundary, the paper-workflow catalog, and the
//! timing protocol. A candidate package supplies an adapter that implements
//! [`WorkflowBackend`] and, for API-level tests, [`QuSpinApi`].

pub mod api;
pub mod benchmark;
pub mod candidate;
mod candidate_workflows;
pub mod catalog;
pub mod observation;

pub use api::{
    AssemblyChecks, BasisHandle, BasisSpec, ComplexCoefficient, Coupling, EigshOptions,
    EvolutionOptions, HamiltonianOptions, LinearOperator, MatrixFormat, OperatorTerm, QuSpinApi,
    SpectrumOptions, SpectrumTarget,
};
pub use benchmark::{benchmark_suite, BenchmarkOptions, BenchmarkRow, WorkflowBackend};
pub use catalog::{paper_workflows, Capability, WorkflowCase, WorkflowOrigin};
pub use observation::{validate_observation, Invariant, Observation, ValidationError};
