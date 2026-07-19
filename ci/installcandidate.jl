#!/usr/bin/env julia

using Pkg

url = get(ENV, "A_REPO_URL", "")
ref = get(ENV, "A_REF", "")
isempty(url) && error("A_REPO_URL not set")
isempty(ref) && error("A_REF not set")

Pkg.activate(".")
Pkg.add(Pkg.PackageSpec(url=url, rev=ref))
Pkg.instantiate()
Pkg.precompile()
