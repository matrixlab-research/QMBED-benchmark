# Full API verification layout

The migration denominator is `../full_api_migration_plan.json`, generated from
the documented QuSpin 1.0.1 public surface.

Tests are grouped by dependency layer:

```text
basis/       state representations, sectors, projections, entropy
operators/   Hamiltonians, linear operators, exponentials, persistence
tools/       evolution, measurements, Lanczos, Floquet, block operations
```

Each documented top-level object eventually owns one Julia test file. A file is
not considered complete until its object and documented members have:

1. a public-surface test;
2. normal and boundary/error unit tests;
3. property tests for algebraic invariants;
4. an integration test when upstream examples use it;
5. a pinned Python-oracle comparison for numerical behavior.

The candidate repository never imports Python at runtime. Python is available
only in the separate oracle CI job.
