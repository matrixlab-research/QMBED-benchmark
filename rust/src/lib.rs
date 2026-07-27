//! Language-neutral verification contracts for the Rust QMBED package.
//!
//! This crate intentionally does not implement exact diagonalization. It owns
//! the independent test-side adapter boundary, the paper-workflow catalog, and the
//! timing protocol. A candidate package supplies an adapter that implements
//! [`WorkflowBackend`] and, for API-level tests, [`QmbedApi`].

pub mod ad_workflows;
pub mod api;
pub mod benchmark;
pub mod candidate;
mod candidate_workflows;
pub mod catalog;
pub mod observation;

pub use api::{
    AssemblyChecks, BasisHandle, BasisSpec, ComplexCoefficient, Coupling, EigshOptions,
    EvolutionOptions, HamiltonianOptions, LinearOperator, MatrixFormat, OperatorSpec, QmbedApi,
    SpectrumOptions, SpectrumTarget,
};
pub use benchmark::{benchmark_suite, BenchmarkOptions, BenchmarkRow, WorkflowBackend};
pub use catalog::{paper_workflows, Capability, WorkflowCase, WorkflowOrigin};
pub use observation::{validate_observation, Invariant, Observation, ValidationError};
