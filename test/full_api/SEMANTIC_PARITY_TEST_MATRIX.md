# Semantic parity test matrix

The verification suite separates public unit regressions from end-to-end
integration checks. A capability is complete only when its implementation,
public regression, and independent verification row all pass.

| Capability | Public regression | Verification integration | Oracle |
| --- | --- | --- | --- |
| Dynamic Hamiltonian transforms | `test/semantic_parity_regressions.jl` | `integration/semantic_parity.jl` | Evaluate the original and transformed operators at several times |
| Dynamic Hamiltonian algebra | `test/semantic_parity_regressions.jl` | `integration/semantic_parity.jl` | Dense matrix sum, product, commutator, and anticommutator |
| Batched Schrödinger evolution | `test/semantic_parity_regressions.jl` | `integration/semantic_parity.jl` | Evolve each state independently |
| Liouville-von Neumann evolution | `test/semantic_parity_regressions.jl` | `operators/hamiltonian_constructor.jl` | Unitary density-matrix propagation and trace preservation |
| Symmetry-changing operator application | `test/semantic_parity_regressions.jl` | `integration/semantic_parity.jl` | Explicit `P_target' * O_full * P_source` |
| Bra/ket operator map | `test/semantic_parity_regressions.jl` | `basis/user_basis.jl` | Direct encoded-state transition table |
| Arbitrary general-basis maps | `test/semantic_parity_regressions.jl` | `integration/general_basis_hamiltonians.jl` | Projected finite-order 2D translations |
| UserBasis callback symmetries | `test/semantic_parity_regressions.jl` | `basis/user_basis.jl` | Callback eigen-sector projection |
| Higher-spin sectors and operators | `test/semantic_parity_regressions.jl` | `basis/discrete_bases.jl` | Exact SU(2) ladder algebra |
| Multi-factor tensor operators | `test/semantic_parity_regressions.jl` | `basis/composite_bases.jl` | Explicit three-factor Kronecker product |
| Photon total-excitation sector | `test/semantic_parity_regressions.jl` | `basis/composite_bases.jl` | Fixed particle plus photon number |
| Spinful sectors and Majorana operators | `test/semantic_parity_regressions.jl` | `basis/discrete_bases.jl` | Multi-sector enumeration and Clifford algebra |
| Dynamic block Hamiltonians | `test/semantic_parity_regressions.jl` | `tools/block_tools.jl` | Block projection at multiple times and evolution |
| Expanded x/y operator forms | `test/semantic_parity_regressions.jl` | `integration/semantic_parity.jl` | Reconstruct original matrices from raising/lowering terms |
| Batched pure-state entropy | `test/semantic_parity_regressions.jl` | `tools/measurements.jl` | Independent entropy call for every state column |
| Harmonic-oscillator basis | `test/semantic_parity_regressions.jl` | `basis/photon_helpers.jl` | Truncated Fock-space ladder matrices |
| Stacked-real and scalar-time ODE evolution | `test/semantic_parity_regressions.jl` | `tools/evolution.jl` | Analytic complex phase rotation |
| Hamiltonian matrix arithmetic and integer powers | `test/semantic_parity_regressions.jl` | `integration/semantic_parity.jl` | Dense evaluation at multiple times |
| QuantumOperator default parameters and arithmetic | `test/semantic_parity_regressions.jl` | `operators/quantum_operator.jl` | Missing parameters equal one; component-wise dense oracle |
| Advanced fermion particle-hole maps | `test/semantic_parity_regressions.jl` | `basis/discrete_bases.jl` | Signed canonical transformation and complementary sectors |
| Sparse projection output controls | `test/semantic_parity_regressions.jl` | `basis/spin_basis_1d.jl`, `basis/discrete_bases.jl` | Sparse/dense type and value for vector and matrix inputs |
| Deferred general-basis construction | `test/semantic_parity_regressions.jl` | `integration/general_basis_hamiltonians.jl` | `make_basis=false` lifecycle against independently eager construction |
| Wide-integer general spin basis | `test/semantic_parity_regressions.jl` | `integration/general_basis_hamiltonians.jl` | `UInt256` states, high-site local action, and sparse Hamiltonian assembly |
| Parallel UserBasis filtering | `test/semantic_parity_regressions.jl` | `basis/user_basis.jl` | Threaded and serial constrained-state enumeration produce identical ordered states |
| General discrete and UserBasis sector shifts | `test/semantic_parity_regressions.jl` | `basis/discrete_bases.jl`, `basis/user_basis.jl` | Direct parent-sector application against full-basis projectors and callback transitions |
| Particle-conserving projection output | `test/semantic_parity_regressions.jl` | `basis/spin_basis_1d.jl` | `pcon=true` returns parent-sector coordinates instead of the full Hilbert space |
| Python-compatible operator archives | `test/runtests.jl` | `operators/archive.jl` | Julia round trip plus bidirectional interchange with Python QuSpin dense and sparse NPZ entries |

New parity work must add a failing public test first, then an independent
integration row here before the implementation is accepted.
