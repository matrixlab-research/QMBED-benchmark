# quspin — Motivation

> Generated from frozen SpecBundle `sha256:a9093788a8b2687689113102daaef231e58d33f6188b74927e9b60c097ee1721` — do not hand-edit.
> Read together with CONTRACT.md (the frozen boundary) and TESTS.md
> (worked examples). This file says WHAT each function is; HOW to
> implement it is entirely yours.

## `error::QuSpinError`

**Why it exists.** Give callers stable, actionable failure categories without parsing strings or catching panics.

**Definition.** Every recoverable public failure maps to one stable category; variants may carry additional context while preserving the category.

**Identities.**
- Valid inputs do not return an error
- Invalid user input never requires a panic
- Non-convergence carries the final iteration count and residual

**Stated properties** (the gate probes these on randomized inputs — they must hold, not just on the examples):
- (compatibility) The listed categories remain distinguishable by pattern matching.
- (safety) Invalid public input returns QuSpinError rather than panicking.
- (diagnostic) Non-convergence exposes iteration count and final residual as structured data.
- (ergonomics) Errors implement the standard Error and Display traits and preserve nested dependency causes when available.

## `basis::Basis`

**Why it exists.** Provide the one abstraction shared by basis enumeration, sparse assembly, and matrix-free action.

**Definition.** A basis is a finite bijection between indices 0 through len-1 and physical states, plus the action of a parsed local operator on one state.

**Identities.**
- index(state(i)) = i for every valid index i
- state(index(s)) = s for every represented state s
- apply_local returns None exactly when the local action annihilates the state or leaves the target sector

**Stated properties** (the gate probes these on randomized inputs — they must hold, not just on the examples):
- (bijection) State and index are mutual inverses on the represented sector.
- (error) Out-of-range indices and unrepresented states return structured errors and never panic.
- (invariance) Physical results are invariant under basis enumeration order.
- (linearity) The amplitude returned by local action is the exact coefficient of the destination state before the coupling coefficient is applied.

## `basis::SpinBasis1D`

**Why it exists.** Represent spin chains in the full space or fixed magnetization and lattice-symmetry sectors.

**Definition.** For spin one-half with up-spin count n, the unsymmetrized sector contains all site bitstrings of length L and Hamming weight n; optional momentum and parity select symmetry-adapted representatives.

**Identities.**
- Without symmetry reduction, len = binomial(L,n)
- The full spin-one-half space has len = 2^L
- pauli=false uses spin operators with Sz eigenvalues plus or minus one half

**Stated properties** (the gate probes these on randomized inputs — they must hold, not just on the examples):
- (combinatorial) The fixed-up-spin unsymmetrized dimension is binomial(L,n).
- (symmetry) A symmetry-resolved basis spans only states compatible with every requested commuting symmetry.
- (normalization) The Pauli and spin-operator conventions differ by the documented factor of two on single-site generators.
- (error) Impossible magnetization or incompatible symmetry labels return an error rather than an empty successful basis.

## `basis::BosonBasis1D`

**Why it exists.** Represent lattice bosons with a finite local occupation cutoff and optional total-number sector.

**Definition.** States are occupation vectors of length L with entries from zero through states_per_site minus one; a fixed-particle sector restricts their sum.

**Identities.**
- Creation on occupation n has amplitude sqrt(n+1) unless the cutoff is reached
- Annihilation on occupation n has amplitude sqrt(n) and annihilates n=0

**Stated properties** (the gate probes these on randomized inputs — they must hold, not just on the examples):
- (domain) Every enumerated local occupation is below states_per_site.
- (conservation) Every state in a fixed-particle sector has the requested total occupation.
- (algebra) Creation and annihilation use square-root occupation amplitudes and respect the cutoff.
- (bijection) Occupation states round-trip through index and state.

## `basis::SpinlessFermionBasis1D`

**Why it exists.** Represent single-flavor lattice fermions while making antisymmetric signs part of local operator action.

**Definition.** States are occupation bitstrings; applying a creation or annihilation operator multiplies by minus one to the number of occupied sites preceding that site in the declared ordering.

**Identities.**
- Without symmetry reduction, len = binomial(L,N)
- Two fermionic creation operators anticommute
- Number operators are diagonal with eigenvalues zero or one

**Stated properties** (the gate probes these on randomized inputs — they must hold, not just on the examples):
- (combinatorial) The fixed-particle unsymmetrized dimension is binomial(L,N).
- (algebra) Local creation and annihilation obey canonical fermionic anticommutation signs.
- (exclusion) Creating on an occupied site or annihilating an empty site returns no destination.
- (symmetry) Symmetry projection preserves physical matrix elements including orbit phases and fermionic signs.

## `basis::SpinfulFermionBasis1D`

**Why it exists.** Represent two fermion flavors for Hubbard-type models with an explicit, stable orbital ordering.

**Definition.** A state is a pair of up and down occupation bitstrings. The sign of every ladder operation follows the documented global ordering of the two flavor orbitals.

**Identities.**
- Without symmetry reduction, len = binomial(L,N_up) times binomial(L,N_down)
- Each flavor separately obeys canonical anticommutation
- Double occupation is allowed and counted by n_up n_down

**Stated properties** (the gate probes these on randomized inputs — they must hold, not just on the examples):
- (combinatorial) The fixed-sector dimension factorizes into the two binomial dimensions.
- (algebra) The public orbital-order convention determines all inter-flavor fermionic signs.
- (conservation) Each enumerated state has the requested particle count for each flavor.
- (bijection) Paired occupation states round-trip through index and state.

## `basis::UserBasis`

**Why it exists.** Support constrained Hilbert spaces without adding a model-specific assembler for every new physical problem.

**Definition.** A user basis combines a state enumerator or filter with named local-action callbacks; after construction it obeys the same Basis contract as built-in families.

**Identities.**
- The universal OperatorBuilder consumes UserBasis through Basis only
- Filtering changes the represented state set, not the semantics of registered local actions

**Stated properties** (the gate probes these on randomized inputs — they must hold, not just on the examples):
- (architecture) UserBasis uses the same operator assembly and matrix-free algorithms as built-in bases.
- (determinism) Repeated construction from the same filter and callbacks yields the same physical state set.
- (error) An unregistered operator name or invalid callback result returns a structured error.
- (bijection) Accepted states round-trip and rejected states are absent.

## `operator::OperatorTerm`

**Why it exists.** Parse and validate a reusable local-operator pattern once before it is applied to many basis states.

**Definition.** An operator term is an ordered local operator string plus a list of complex coefficients and zero-based site tuples of matching arity.

**Identities.**
- The total operator is linear in all coupling coefficients
- Operator-string order fixes fermionic operator order
- Parsing cost is independent of basis dimension and is not repeated per state

**Stated properties** (the gate probes these on randomized inputs — they must hold, not just on the examples):
- (validation) Every coupling site tuple has exactly the operator-string arity.
- (indexing) All public site indices are zero-based and checked against the basis length at build time.
- (linearity) Scaling all couplings by a complex scalar scales the represented linear map by that scalar.
- (complexity) A successfully constructed term stores parsed operator actions and does not reparse its string inside state enumeration.

## `operator::Coupling`

**Why it exists.** Carry one local-term coefficient and its ordered sites without coupling the API to a sparse-matrix crate.

**Definition.** A coupling is a complex scalar and an ordered zero-based site tuple; the tuple order matches the operator-string order.

**Identities.**
- Changing the coefficient scales only this contribution
- Permuting sites changes the operator unless the local actions commute

**Stated properties** (the gate probes these on randomized inputs — they must hold, not just on the examples):
- (domain) Real and complex coefficients are represented without precision loss beyond Complex64.
- (indexing) Site order is preserved exactly as supplied.
- (indexing) Public sites are zero-based.
- (validation) OperatorTerm validates the site count against the operator string.

## `operator::MatrixFormat`

**Why it exists.** Make storage and matrix-free execution an explicit backend choice rather than the definition of a Hamiltonian.

**Definition.** MatrixFormat selects materialization only; all formats represent the same linear map within numeric tolerance.

**Identities.**
- Format conversion preserves shape and action
- Sparse formats are assembled without first allocating a dense matrix

**Stated properties** (the gate probes these on randomized inputs — they must hold, not just on the examples):
- (equivalence) Every supported format produces the same action on the same vector within tolerance.
- (complexity) CSC, CSR, DIA, and MatrixFree construction do not allocate a dense rows-by-columns intermediary.
- (error) An unsupported format and operator combination returns a structured error.
- (identity) A built operator reports the format actually used.

## `operator::LinearOperator`

**Why it exists.** Give stored and matrix-free algorithms one rectangular-capable narrow waist.

**Definition.** A linear operator maps a vector of length shape.1 to one of length shape.0; apply overwrites an output buffer with that action.

**Identities.**
- apply(a x + b y) = a apply(x) + b apply(y)
- For a square Hamiltonian, shape.0 = shape.1
- A sector-changing probe may have different row and column dimensions

**Stated properties** (the gate probes these on randomized inputs — they must hold, not just on the examples):
- (linearity) Apply is complex-linear in the input vector.
- (shape) Input length equals column count and output length equals row count; mismatches return an error.
- (mutation) Apply overwrites every output entry and does not depend on its previous contents.
- (equivalence) Stored and matrix-free implementations with the same terms agree on action.

## `operator::OperatorBuilder`

**Why it exists.** Assemble square Hamiltonians and rectangular sector-changing probes through one universal path.

**Definition.** For each source basis state and term coupling, apply the local operator, map any destination into the target basis, and accumulate its coefficient into a format-specific sink.

**Identities.**
- The assembled map is the sum of all term contributions
- Duplicate row-column contributions are summed
- on(basis) is equivalent to between(basis,basis)

**Stated properties** (the gate probes these on randomized inputs — they must hold, not just on the examples):
- (architecture) The same traversal handles spin, boson, fermion, and user bases through the Basis contract.
- (shape) Between source and target bases produces shape target.len by source.len.
- (algebra) Contributions to the same matrix entry are summed exactly once in the final map.
- (validation) Hermiticity, particle conservation, and symmetry compatibility checks are individually configurable and return structured errors.
- (complexity) Sparse construction streams nonzero contributions and never materializes a dense intermediate.

## `solve::EigshOptions`

**Why it exists.** Make solver accuracy, targeting, resource limits, and reproducibility explicit at each call.

**Definition.** Options select the number and spectral region of eigenpairs, the convergence tolerance and iteration cap, an optional Krylov dimension, and a deterministic seed.

**Identities.**
- The same options and operator define the same target subspace within tolerance
- A shift target is ordered by distance to the shift

**Stated properties** (the gate probes these on randomized inputs — they must hold, not just on the examples):
- (validation) Eigenpair count, tolerance, and maximum iterations are positive.
- (validation) When supplied, Krylov dimension is greater than the requested eigenpair count.
- (validation) Shift targets are finite real values.
- (determinism) Seed controls all randomized initialization used by the solver.

## `solve::eigsh`

**Why it exists.** Compute selected eigenpairs without dense diagonalization of the full Hilbert space.

**Definition.** Return Ritz values and normalized Ritz vectors for the requested spectral target, together with residual norms and convergence metadata.

**Identities.**
- For every returned pair, residual = norm(A v - lambda v)
- Returned eigenvectors are orthonormal within tolerance
- Shift targeting orders results by distance to the shift

**Stated properties** (the gate probes these on randomized inputs — they must hold, not just on the examples):
- (residual) Every converged eigenpair satisfies its requested residual tolerance up to documented numerical slack.
- (normalization) Returned eigenvectors have unit norm and are mutually orthogonal within tolerance.
- (determinism) With the same operator, options, and seed, eigenvalues and invariant subspaces are reproducible within tolerance.
- (error) Rectangular or declared non-Hermitian operators are rejected by the Hermitian solver.
- (error) Non-convergence returns iteration count and final residual rather than a panic or silent partial success.

## `solve::EvolutionOptions`

**Why it exists.** Define the requested time samples and Krylov error controls without changing the operator abstraction.

**Definition.** Options contain an absolute time grid and numerical controls; hamiltonian=true interprets the supplied operator as H and evolves with minus i H.

**Identities.**
- Output states correspond one-to-one and in order with times
- Time zero is allowed
- Hamiltonian mode is norm preserving for Hermitian input

**Stated properties** (the gate probes these on randomized inputs — they must hold, not just on the examples):
- (validation) Times are finite and nondecreasing.
- (validation) Krylov dimension, tolerance, and maximum substeps are positive.
- (shape) The trajectory contains exactly one state for each requested time.
- (semantics) Hamiltonian mode applies the factor minus i exactly once.

## `solve::evolve`

**Why it exists.** Apply matrix-exponential time evolution without densifying a sparse or matrix-free generator.

**Definition.** For a time-independent generator A, the state at time t is exp(t A) applied to the initial state; Hermitian Hamiltonian evolution uses A = -i H.

**Identities.**
- The state at time zero equals the initial state
- For Hermitian H, norm is conserved
- Evolution over t1 followed by t2 equals evolution over t1+t2 within solver tolerance

**Stated properties** (the gate probes these on randomized inputs — they must hold, not just on the examples):
- (special_value) Time zero returns the initial vector exactly up to representation-preserving copy semantics.
- (conservation) Hermitian Hamiltonian evolution preserves vector norm within requested tolerance.
- (composition) For a time-independent generator, consecutive intervals compose to their summed interval within tolerance.
- (complexity) Evolution consumes LinearOperator action and does not require a dense materialization.
- (error) A rectangular operator or mismatched initial vector returns a structured dimension error.

## `dynamics::DriveStep`

**Why it exists.** Pair one piecewise-constant Hamiltonian with the duration for which it acts.

**Definition.** A drive step denotes evolution by exp(-i H duration) and carries shared ownership of H so a Floquet sequence can be reused.

**Identities.**
- Zero duration is the identity
- Equal Hamiltonian and duration define equal physical action

**Stated properties** (the gate probes these on randomized inputs — they must hold, not just on the examples):
- (validation) Duration is finite and nonnegative.
- (validation) The Hamiltonian is square.
- (validation) The Hamiltonian satisfies the Hermitian contract required for unitary step evolution.
- (ownership) A step can share one immutable operator across reusable Floquet sequences.

## `dynamics::Floquet`

**Why it exists.** Represent one period of a piecewise-constant drive without requiring a precomputed dense unitary.

**Definition.** For ordered steps from first to last, the one-period map is the product exp(-i H_last dt_last) through exp(-i H_first dt_first), acting on states from right to left.

**Identities.**
- Every Hermitian-step Floquet map is unitary
- Zero-duration steps are identities
- Repeated periods compose by repeated application of the same map

**Stated properties** (the gate probes these on randomized inputs — they must hold, not just on the examples):
- (composition) Drive steps are applied in the declared physical time order.
- (unitarity) A sequence of Hermitian drive steps preserves norm and has unit-modulus eigenphases within tolerance.
- (special_value) A zero-duration step has no effect.
- (error) Mismatched step dimensions or rectangular generators return structured errors.
- (complexity) Period action can remain matrix-free even when optional full-unitary materialization is requested only for small systems.

## `dynamics::SpectrumOptions`

**Why it exists.** Define a response-frequency grid and the resolvent broadening and numerical controls.

**Definition.** Options specify omega values, the source reference energy E0, positive Lorentzian broadening eta, and Krylov resolution controls.

**Identities.**
- Output order matches frequency order
- Positive broadening keeps every resolvent point off the real-axis poles

**Stated properties** (the gate probes these on randomized inputs — they must hold, not just on the examples):
- (validation) Frequencies are finite and the grid is nonempty.
- (validation) Broadening is finite and strictly positive.
- (validation) Krylov dimension and tolerance are positive.
- (shape) The returned spectrum has one value for each frequency in input order.

## `dynamics::spectral_function`

**Why it exists.** Compute response spectra in the same sector or after a sector-changing probe without full diagonalization.

**Definition.** For probe-created vector phi and reference energy E0, the broadened spectrum is minus one over pi times the imaginary part of the resolvent expectation of phi at omega plus E0 plus i eta.

**Identities.**
- The spectrum is nonnegative for positive broadening and a Hermitian target Hamiltonian
- A zero probe gives the zero spectrum
- The probe shape is target dimension by source dimension

**Stated properties** (the gate probes these on randomized inputs — they must hold, not just on the examples):
- (shape) Probe columns match the source-state length and probe rows match the target Hamiltonian dimension.
- (positivity) For Hermitian input and positive broadening, the spectral density is nonnegative up to numerical tolerance.
- (special_value) A probe that annihilates the source produces zero at every frequency.
- (numerics) Finite valid inputs produce finite outputs across the requested grid.
- (complexity) The resolvent or Krylov evaluation uses operator action and does not require full eigenvectors.

## `measure::Subspace`

**Why it exists.** Represent an eigenspace independently of a particular eigenvector gauge.

**Definition.** Subspace stores or derives an orthonormal column basis for the span of the supplied vectors and exposes ambient dimension and numerical rank.

**Identities.**
- Right multiplication by an invertible matrix preserves the subspace
- Orthonormalization preserves the projector onto the span

**Stated properties** (the gate probes these on randomized inputs — they must hold, not just on the examples):
- (invariance) Equivalent full-rank spanning sets produce the same orthogonal projector within tolerance.
- (normalization) The internal or returned basis columns are orthonormal within tolerance.
- (validation) Input storage length equals ambient dimension times declared rank.
- (validation) Numerically rank-deficient or empty spanning sets return a structured error.

## `measure::subspace_fidelity`

**Why it exists.** Compare possibly degenerate eigenspaces without depending on arbitrary eigenvector gauges or rotations.

**Definition.** After orthonormalizing bases Q_left and Q_right, fidelity is the mean squared singular value of Q_left adjoint times Q_right, equivalently the normalized trace of the product of their projectors.

**Identities.**
- Fidelity lies from zero through one
- It is invariant under unitary rotations within either subspace
- It equals one for identical subspaces and zero for orthogonal equal-rank subspaces

**Stated properties** (the gate probes these on randomized inputs — they must hold, not just on the examples):
- (bounds) The result lies in the closed interval from zero to one up to tolerance.
- (invariance) Independent unitary rotations of either spanning set leave the result unchanged.
- (symmetry) Swapping left and right leaves the result unchanged.
- (special_value) Two spanning sets for the same subspace have fidelity one.
- (error) Ambient-dimension mismatch or an empty subspace returns a structured error.

## `workflow::LindbladGenerator`

**Why it exists.** Express Markovian open-system evolution through the same matrix-free linear-operator algorithms as closed-system dynamics.

**Definition.** The generator acts as -i times the commutator of H and rho plus the sum over jumps of L rho L adjoint minus one half of the anticommutator of L adjoint L with rho.

**Identities.**
- The generator preserves trace
- Hermiticity is preserved
- With no jumps it reduces to unitary commutator evolution

**Stated properties** (the gate probes these on randomized inputs — they must hold, not just on the examples):
- (shape) A Hilbert dimension d produces a square linear operator of dimension d squared.
- (conservation) The derivative of density-matrix trace is zero within tolerance.
- (symmetry) Hermitian density matrices have Hermitian derivatives.
- (special_value) With no jumps, action equals the vectorized Hamiltonian commutator generator.
- (complexity) The generator can apply H and jumps without explicitly materializing the d squared by d squared superoperator.
