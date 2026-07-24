#!/usr/bin/env julia
# Install the requested Julia compatibility candidate at an explicit ref via
# Pkg.add(url, rev) — Julia 1.12 removed `rev` from Pkg.develop — then run the
# public conformance suite against it. The workflow supplies both inputs.
using Pkg
url = get(ENV, "A_REPO_URL", "")
ref = get(ENV, "A_REF", "")
isempty(url) && error("A_REPO_URL not set")
isempty(ref) && error("A_REF not set")

Pkg.activate(".")
Pkg.add(Pkg.PackageSpec(url=url, rev=ref))
Pkg.instantiate()
Pkg.test()
