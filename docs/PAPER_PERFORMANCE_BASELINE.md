# Original six-workflow performance baseline

This file preserves the first six-case snapshot for historical comparison.
The live `paper` tier now contains twelve paired workflows; the six
literature-derived medium-size additions first appear in the CSV/SVG artifacts
from the next successful `paper`-tier Actions run and are not backfilled with
invented local numbers here.

This is a local, paired five-sample baseline for the six paper-shaped workflows
defined in `benchmark/python_paper_workflow_benchmarks.py` and
`benchmark/julia_paper_workflow_benchmarks.jl`. It is committed as an
interpretation aid; GitHub Actions artifacts remain the authoritative
same-runner result for a candidate commit.

## Protocol

- machine: Apple arm64 MacBook Pro, macOS 14.6.1;
- Python 3.12.8, packaged QuSpin reporting version 0.3.7;
- Julia 1.12.1, QuSpin.jl based on `bb7f150` plus pivoted sparse
  shift-invert and consolidated sparse-Hamiltonian assembly;
- one Julia, BLAS, OpenMP, and native math thread;
- complete basis + Hamiltonian + solver/observable pipeline per sample;
- physical correctness preflight, one warm-up, then five samples;
- identical model parameters, boundaries, solver tolerances, `ncv`,
  `maxiter`, and basis-order-independent uniform solver seed;
- Python measured first and Julia second, without concurrent workload.

Speedup is `Python median / Julia median`; values above one mean Julia is
faster. IQR is the 25th–75th percentile range.

| Workflow | Size | Python median [IQR] (ms) | Julia median [IQR] (ms) | Julia speedup |
|---|---|---:|---:|---:|
| MBL mid-spectrum shift-invert | `L=14`, `Nup=7`, dim 3432, `k=6` | 226.124 [224.592–230.268] | 86.318 [86.098–86.552] | **2.62×** |
| XXZ Lanczos quench | `L=16`, `Nup=8`, dim 12870, `m=80` | 168.630 [166.537–170.894] | 141.422 [138.378–141.597] | **1.19×** |
| Floquet full-unitary heating | `L=9`, dim 512, two steps | 360.041 [354.363–363.748] | 290.404 [286.947–296.803] | **1.24×** |
| Spinful Hubbard low-energy spectrum | `L=8`, `Nf=(4,4)`, dim 4900, `k=6` | 23.485 [23.455–23.656] | 35.735 [35.722–35.833] | 0.66× |
| interacting SSH low-energy spectrum | `L=16`, `Nf=8`, dim 12870, `k=6` | 43.766 [42.857–43.970] | 84.842 [84.645–84.871] | 0.52× |
| translation-sector XXZ spectrum | `L=18`, `Nup=9`, `k=0`, `k_eigs=4` | 19.190 [19.154–19.233] | 40.402 [39.270–40.892] | 0.47× |

## Interpretation

The result is mixed rather than a blanket Julia win:

- Julia has demonstrated end-to-end gains for the MBL shift-invert, Lanczos
  quench, and dense Floquet workflows.
- Spinful and spinless interacting-fermion construction plus eigensolve, and
  translation-sector construction, remain slower than Python and are the next
  optimization targets. Consolidated assembly reduced their Julia medians
  from 45.237, 161.250, and 60.468 ms to 35.735, 84.842, and 40.402 ms.
- The geometric-mean speedup across all six cases is `0.92×`, so these data do
  **not** support an overall performance-gain claim yet.
- The MBL benchmark initially exposed a real correctness failure:
  `sigma=0` selected unpivoted CHOLMOD LDLᵀ and raised a zero-pivot error on an
  invertible indefinite matrix. QuSpin.jl now forms a pivoted-LU
  matrix-free inverse operator and the paper case passes its residual check.

Hosted runners vary, so no performance threshold is used as a correctness
gate. The committed CI workflow records raw samples and reports missing
language counterparts explicitly.

## Same-runner GitHub Actions snapshot

The first complete `paper`-tier run on GitHub Actions
([run 29749615993](https://github.com/matrixlab-research/QMBED-benchmark/actions/runs/29749615993))
passed the oracle, all 670 verification assertions, all 21 small workflows,
and all six paired paper workflows. On that Ubuntu x86-64 runner the paired
medians were:

| Workflow | Python (ms) | Julia (ms) | Julia speedup |
|---|---:|---:|---:|
| MBL mid-spectrum shift-invert | 346.837 | 181.030 | **1.92×** |
| XXZ Lanczos quench | 225.234 | 165.977 | **1.36×** |
| Floquet full-unitary heating | 461.019 | 381.832 | **1.21×** |
| Spinful Hubbard low-energy spectrum | 44.101 | 45.935 | 0.96× |
| interacting SSH low-energy spectrum | 75.870 | 92.765 | 0.82× |
| translation-sector XXZ spectrum | 51.539 | 74.088 | 0.70× |

The geometric-mean speedup is `1.09×` on this runner. Together with the local
`0.92×` result, this shows that the candidate has real workload-specific gains
but does not justify a hardware-independent blanket speedup claim.
