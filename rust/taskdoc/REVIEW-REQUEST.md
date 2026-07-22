# Independent review — task package for quspin

You are an INDEPENDENT spec reviewer. You had no part in writing these
documents, and you must not request or consult any held-out material (the
private verification repo / its test suite, or the original implementation's
source). Your
inputs are exactly: `MOTIVATION.md`, `CONTRACT.md`, `TESTS.md`, and
`audit-report.md` (the mechanical audit; triage its warnings). The package
derives from frozen bundle `sha256:a9093788a8b2687689113102daaef231e58d33f6188b74927e9b60c097ee1721`.

The implementing agent will see ONLY these documents. Anything wrong, missing,
or leaked here goes straight into the rewrite. Work through the four
questions; be adversarial.

## 1. Deviation — do the three documents tell one consistent story?

- Does every definition in MOTIVATION.md agree with the worked examples in
  TESTS.md? Spot-check by computing several examples from the definition
  alone (high-precision reasoning). A mismatch means the doc drifted from the
  bundle even if every number is verbatim.
- Do CONTRACT.md signatures/domains agree with how MOTIVATION.md describes
  each function's arguments?

## 2. Correctness — is the content right in its own domain?

- Are the stated identities/properties actually true of the defined behavior
  (mathematical, algorithmic, or formatting — whatever the package's domain)?
- Where `num` comparisons are used: are tolerances sane for the value ranges
  shown (an atol of 1e-15 around values of 1e300 is a red flag)? Where
  `isequal`/`str`/`custom` are used: is the expected value's form unambiguous
  (encoding, ordering, whitespace)?
- Would YOU be able to pin down each function unambiguously from these
  documents alone? Name any function where two reasonable implementers could
  build observably different things.

## 3. Leakage — information that must not be here

- No implementation code, in any language, in any form — including "pseudo
  code", coefficient/lookup tables, magic constants, branch thresholds (e.g.
  "switch formula at x = 33.3"), or step-by-step algorithm recipes.
  Natural-language explanation of the intended behavior is fine; a recipe
  precise enough to transcribe into code is not.
- No numeric values beyond those in TESTS.md's table (the audit scans for
  known held-out values, but cannot catch paraphrased or rounded leaks — that
  is YOUR job).
- Nothing that reveals how the ORIGINAL package implements the function.

## 4. Scope — exactly the API, nothing else

- The set of functions the documents ask for must be exactly CONTRACT.md's
  symbol table. Flag any prose that smuggles in an extra deliverable ("you
  will also need a helper that...", "should also export...").
- Flag requirements that the gate does not measure (style, naming, internal
  structure) — the contract's promise is "green pipeline = accepted", so
  unmeasured demands must be removed, not left as vibes.

## Output

Write `REVIEW.md` next to these documents:

- First line: `VERDICT: PASS` or `VERDICT: BLOCK`.
- Then findings, one bullet each: file, location, what is wrong, and the
  minimal fix. No finding, no bullet — do not pad.
- Constraint on your own output: never write a numeric value that does not
  already appear in TESTS.md. Refer to other values by description ("the
  large-negative anchor for `log1pexp`"). Your REVIEW.md will itself be
  re-audited for leaks (`minos audit` scans every .md in the directory).

BLOCK on: any leakage finding, any incorrect mathematics, any scope
violation, or any function you judge unpinnable from the documents alone.
Deviations with an obvious mechanical fix may PASS with findings only if
they cannot mislead the implementer.
