#!/usr/bin/env julia
# the gate: bring in the package the agent's MR proposes AT THE MR COMMIT via
# Pkg.add(url, rev) — julia 1.12 removed `rev` from Pkg.develop — then run the
# private suite against it. A_REPO_URL is pinned in .gitlab-ci.yml (never a
# trigger variable, so a forged trigger can't point the gate elsewhere); A_REF
# is injected by the trigger.
using Pkg
url = get(ENV, "A_REPO_URL", "")
ref = get(ENV, "A_REF", "")
isempty(url) && error("A_REPO_URL not set — pin it in .gitlab-ci.yml")
isempty(ref) && error("A_REF not set — trigger me from the task bridge")

Pkg.activate(".")
Pkg.add(Pkg.PackageSpec(url=url, rev=ref))
Pkg.instantiate()
Pkg.test()
