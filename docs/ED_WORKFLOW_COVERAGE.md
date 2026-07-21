# ED workflow coverage

This verification layer starts from 13 scientific application areas identified
in a literature search for `exact diagonalization` on 2026-07-20. LKM was used
as one paper index during the search; its internal result identifiers are not
scientific references. The paper sources for the six timed workflows are
linked in the README. Each row below states explicitly whether the current
QuSpin.jl scenario is direct, assisted, or only a proxy.

- **Direct**: the package expresses the model and the main ED observable.
- **Assisted**: the ED kernel is covered, but a paper workflow still needs
  external parameter sweeps, integral generation, state tracking, or
  post-processing.
- **Proxy**: the small test exercises relevant algebra but does not establish
  coverage of the named paper domain.

The integration sizes are intentionally small and deterministic. They are
correctness and user-scenario checks, not attempted paper reproductions.

| ID | Literature topic | ED application | Integration size | Verified quantity | Coverage | Representative paper-scale range |
|---:|---|---|---|---|---|---|
| 01 | spin-chain spectra | Heisenberg dimer | `L=2`, dim 4 | singlet/triplet energies and degeneracy | Direct | dimer is exact; clusters `N=2–16` |
| 02 | spin-chain symmetries | symmetry-resolved XXZ spectrum | `L=6`, `Nup=3`, `k=0` | sector spectrum embedding and sparse `eigsh` | Direct | chains `L≈20–36` with symmetries |
| 03 | frustrated magnetism | frustrated `J1-J2` spin gap | `L=8` | `Sz=1` minus `Sz=0` ground energy | Direct | 1D `L≈24–36`; 2D clusters `N≈20–40` |
| 04 | quantum phase transitions | transverse-field Ising scan | `L=7`, dim 128 | phase-dependent finite-size gap | Direct | `L≈20–28` |
| 05 | finite-temperature spin ED | finite-temperature susceptibility | dimer, dim 4 | thermal `β⟨Mz²⟩` against the analytic result | Assisted | full spectrum `N≈2–16`; FTLM `N≈20–32` |
| 06 | variational certification | variational-state ED certificate | `L=6`, `Nup=3` | `Evar ≥ E0` | Assisted | `L≈16–36`; ansatz optimization is external |
| 07 | Floquet spin dynamics | step-driven Floquet chain | `L=4`, dim 16 | unitary and unit-circle quasienergies | Assisted | explicit full `UF`, usually `L≈12–18` |
| 08 | Floquet effective theory | high-frequency effective Hamiltonian | two-level system | `HF` approaches the average Hamiltonian | Assisted | spin chains `L≈8–16`; branch choice matters |
| 09 | disordered localization | Anderson localization | `L=10`, one particle | disorder-enhanced IPR | Assisted | specialized single-particle solvers reach much larger `L` |
| 10 | many-body localization | MBL mid-spectrum statistics | `L=8`, `Nup=4` | shift-invert states and adjacent-gap ratios | Assisted | `L≈16–22`, plus disorder averaging |
| 11 | SSH topology | SSH edge modes | `L=8`, one particle | near-zero modes and edge weight | Direct | single-particle `L≈10²–10⁴` |
| 12 | interacting SSH systems | interacting extended SSH | `L=6`, `Nf=3` | interaction-dependent many-body gap | Direct | `L≈12–20` |
| 13 | chiral lattice models | chiral SSH spectrum | `L=7`, one particle | `E ↔ -E` pairing and zero mode | Direct | single-particle `L≈10²–10⁴` |
| 14 | topological flux response | flux-insertion diagnostic | `L=7`, one particle | `2π` spectral periodicity | Assisted | many-body `L≈12–24`; Berry tracking is external |
| 15 | few-electron quantum dots | few-electron dot shell filling | 4 orbitals, 2 electrons | addition energy and diagonal sparsity | Assisted | 10–30 orbitals, 2–8 electrons |
| 16 | few-site Hubbard systems | two-electron Hubbard dot | 2 sites, dim 4 | analytic singlet ground energy | Assisted | 2–20 orbitals for two electrons |
| 17 | multiorbital Hubbard clusters | Cu–O orbital-hole fractions | 3 orbitals, dim 9 | total and orbital-resolved occupancy | Assisted | 5–16 orbital clusters |
| 18 | correlated-cluster gaps | cluster charge gap | 5 orbitals, `N=1,2,3` | `E(N+1)+E(N-1)-2E(N)` | Assisted | Hubbard clusters around 10–16 sites |
| 19 | fractional quantum Hall spectra | interacting flux spectral flow | 6 orbitals, 2 fermions | degenerate multiplet and `2π` flow | **Proxy** | FQHE needs Landau-level and magnetic-translation structure |
| 20 | collective emission | collective bright/dark modes | 4 emitters, one excitation | non-Hermitian decay eigenvalues and dark-state count | Assisted | one-excitation `N≈10²–10⁴`; Lindblad is not covered |
| 21 | fractional Chern insulators | many-body entanglement diagnostic | 6 orbitals, 3 fermions | reduced density matrix and entropy | **Proxy** | FCI needs projected Chern bands and 2D momentum sectors |

## What the suite proves

`test/full_api/integration/ed_workflow_catalog.jl` proves that the current
candidate can execute 21 deterministic ED user scenarios and satisfy 38
physical/numerical invariants. It does **not** prove complete FQHE, FCI,
Lindblad, two-dimensional space-group, SU(2), Berry/Chern, or Wilson-loop
support.

The principal current boundaries are:

- spin models are spin-1/2 and built-in spatial symmetries are one-dimensional;
- FQHE/FCI tests are honest algebraic proxies, not paper-level implementations;
- collective emission uses an effective non-Hermitian Hamiltonian, not a
  density-matrix Lindblad solver;
- Floquet full-spectrum paths materialize a dense unitary and scale as
  `O(dim²)` memory and `O(dim³)` diagonalization time;
- paper studies normally require sweep/disorder/statistical orchestration above
  the package API.

## Timing layers

1. `julia_ed_workflow_benchmarks.jl` times all small end-to-end scenarios after
   a correctness preflight. These timings detect regressions and are not
   Python/Julia speedup claims.
2. The paired `python_paper_workflow_benchmarks.py` and
   `julia_paper_workflow_benchmarks.jl` scripts cover six common paper-shaped
   workloads:

   - MBL shift-invert, `L=14`, `Nup=7`, dim 3432;
   - XXZ Lanczos quench, `L=16`, `Nup=8`, dim 12870;
   - Floquet full-unitary heating, `L=9`, dim 512;
   - spinful Hubbard, `L=8`, `Nf=(4,4)`, dim 4900;
   - interacting SSH, `L=16`, `Nf=8`, dim 12870;
   - translation-sector XXZ, `L=18`, `Nup=9`.

Paper benchmarks run only when `benchmark_tier=paper` is selected. Each
operation validates residual, norm, or unitarity before measuring one warm-up
and five independent samples. The CSV artifacts retain raw samples; speedup is
reported as Python median divided by Julia median on the same runner.
