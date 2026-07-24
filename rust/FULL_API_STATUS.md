# QMBED full-package candidate status

Date: 2026-07-23  
Pinned candidate: `85a3f5d44b2f5b5bce58e0c35b3e6d31148d43ce`

## Completion boundary

The frozen source denominator is Python QuSpin 1.0.1: 64 public objects, 282
non-constructor methods, and 180 documented attributes. Rust intentionally
maps these into fewer mathematical interfaces. Completion means that every
scientific behavior family has a public Rust boundary and passes public
properties, independent numerical oracles, composed workflows, and relevant
structural scale checks. It does not mean one Rust method for every Python
spelling.

| Scientific family | Rust boundary | Evidence state |
|---|---|---|
| State integers and helpers | `WideState`, `U256`–`U16384`, `BigUint` round trip, bit helpers | complete candidate |
| Built-in bases | spin/higher-spin, boson, spinless/spinful fermion, fixed and union sectors | complete candidate |
| General and composite bases | finite-map sectors, projectors, tensor, photon, deferred/branching user bases | complete candidate |
| Local and cross-sector action | universal transitions, bra/ket tables, streamed sector shift, CSC/CSR/DIA/matrix-free assembly | complete candidate |
| Operator families | static/dynamic `Hamiltonian`, `QuantumOperator`, `QuantumLinearOperator`, sparse algebra and transforms | complete candidate |
| Exponential and low-level action | `ExpOp` grid/iterator/right action/transforms, matvec/matmat plans | complete candidate |
| Solvers and evolution | complex Hermitian dense/partial spectra, all targets, Lanczos, cached shift-invert, batches, density, FTLM/LTLM | complete candidate |
| Dynamics and measurements | Floquet, block tools, response, arbitrary-site pure/mixed entropy, ensembles, statistics | complete candidate |
| Persistence and workflows | safe named dense/sparse archives, state tracking, Lindblad generator, twelve paper workflows | complete candidate |

## Active evidence

- Public crate: all targets, strict Clippy, rustdoc warnings, and twelve ignored
  release-mode paper workflows.
- Independent oracle gate: seven families covering capabilities outside
  the original 23-symbol workflow adapter.
- Independent scale gate: four tests sweeping several sizes and asserting sparse,
  sector-sized, or transition-sized storage rather than full-square fallback.
- Workflow gate: all twelve medium-size literature-derived cases must satisfy
  physical invariants before their timings are accepted.
- Performance continuity: the current five-sample local record uses one
  streamed transition/assembler contract and one coarse numerical-backend
  boundary. Real CSC shift-invert stays real through factorization and
  Lanczos, reducing the MBL median from 0.403 s to 0.114 s while the summed
  twelve-workflow median falls from 8.99 s to 8.79 s.

## Language-boundary exclusions

Python object identity, pickle execution, NumPy dtype aliases, warning text,
and Python-specific mutable iterator state are not Rust scientific
capabilities. Their behavior is represented by typed Rust values, safe named
archives, `Result` errors, and explicit iterators/plans. Cross-language timing
on a common hosted runner remains useful performance evidence, but is not a
semantic completion requirement.
