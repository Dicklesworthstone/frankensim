# Claim-Integrity Defect Class

Bead `frankensim-extreal-program-f85xj.2.1`. This file is normative: the E02
sweep (`.2.2`), the promotion gate (`.2.3`), and the gate drills (`.2.4`)
consume the decision rules, severity rules, and label taxonomy below verbatim.

FrankenSim's stated value is not that it computes numbers; it is that a number
arrives with an honest account of what makes it trustworthy. A defect that lets
a surface assert a *stronger* epistemic state than its evidence supports is
therefore worse than an ordinary wrong answer: a wrong answer is one bad value,
while a false claim corrupts every downstream composition that trusts the claim
and, because the composition rules are supposed to prevent exactly this, it
does so while looking correct. AGENTS.md already says it: "A false certificate
is worse than an ordinary wrong answer." This document turns that sentence into
a countable, gateable defect class.

## Definition

A **claim-integrity defect** exists when any public surface can assert a
stronger epistemic state than its actual evidence establishes.

The three terms are load-bearing:

- **Public surface** — an API return type or its documented meaning, an
  evidence color, a certificate or receipt, a report/notebook line, a
  `CONTRACT.md` or `README.md` sentence, a WASM/CLI export, a ledger row, or a
  package claim. Anything a caller or reader is invited to rely on.
- **Stronger epistemic state** — a claim higher on the evidence order than
  what was proved. The order this repo uses: refusal/no-claim < Estimated <
  Validated (in its stated domain) < Verified, plus the quantitative analogues:
  a bound is stronger than an estimate, an exact count is stronger than a lower
  bound, an enclosure is stronger than a rounded interval, "the boundary was
  located" is stronger than "no transition was observed", and "the sweep
  agreed" is stronger than "both sweeps were truncated at the same endpoint".
- **Actual evidence** — what the executed code path actually establishes on
  the inputs that can reach it, not what the algorithm establishes in exact
  real arithmetic, not what the test fixtures happen to exercise, and not what
  the author intended.

The defect is in the *gap*, not in the wrongness. A value may be numerically
correct on every input ever tried and still be a claim-integrity defect if the
surface asserts an authority the code cannot deliver on some reachable input.

## Decision rules

Apply these in order. The procedure is deliberately mechanical so two auditors
reach the same verdict.

1. **State the strongest permitted inference.** Read the type, the doc comment,
   the CONTRACT line, and the report text. Write down the strongest claim a
   competent caller is entitled to infer. Use the caller's reading, not the
   author's intent: if `component_count: usize` is returned, the caller is
   entitled to read it as an exact count.
2. **State what the code establishes.** Identify the actual argument: which
   enclosure, residual, invariant, sampled set, or admitted receipt backs the
   value, and under which preconditions.
3. **Hunt for a reachable gap instance.** Look for an input or state where (2)
   is weaker than (1). Reachability is what separates severities, so record how
   it is reached (default path, feature flag, direct API call, doc only).
4. **Classify the gap.** It is a claim-integrity defect if the gap is an
   *epistemic* overstatement. It is an ordinary bug if the surface's claim is
   right and the value is merely wrong (see the exclusions below).
5. **Write the honest claim.** Every filed defect must name the claim the
   surface *should* make (`ComponentCountLowerBound(1)` rather than
   `component_count = 1`). A finding without an honest replacement is an
   incomplete finding: the fix bead needs a destination.

A gap is present, and the defect is real, whenever any of the following hold.
These are the recurring shapes, each drawn from a confirmed instance in this
repo's history:

- **Existence read as exactness.** A bounded/local argument proves *at least
  one* or *within this frame*, and the surface emits an exact global count or
  an unqualified total. (`frankensim-0547u`: an interval frame proving a
  connected negative subset exists is published as `component_count = 1`.)
- **Truncation read as agreement.** Two predicates are compared only at their
  sampled extremes, so a sweep that never observed the phenomenon reports
  agreement because both truncated maxima equal the sweep endpoint.
  (`frankensim-flutter-boundary-witness-hukmw`: `boundaries_agree` true with no
  witnessed stable-to-unstable transition.)
- **Real-arithmetic theory read as an executable bound.** A textbook bound that
  holds in exact arithmetic is emitted as `Verified` after being computed in
  round-to-nearest `f64` with no outward rounding. (`frankensim-y6yv`:
  fs-surrogate ladder residual/QoI diagnostics.)
- **Laundering.** Composition, aggregation, or admission yields a claim
  stronger than its weakest required ingredient, including the specific
  numerical-to-model-form channel where a tightly converged solve is read as
  evidence the closure model is right. (`frankensim-f5pr`: disjoint-regime
  composition recording dispersion `0.0`, i.e. perfect agreement, where there
  is no jointly supported state.)
- **Silent-pass / fail-open gates.** A checker, admission path, or CI gate that
  returns pass on malformed, empty, missing, or unparseable evidence. A gate
  that errors open is itself a claim-integrity defect, because its green is
  read as evidence. (`frankensim-2gs5h`: an empty golden `schema_fingerprint`
  accepted as matching.)
- **Forgeable certificates.** A certificate predicate that a caller-controlled
  input can satisfy without the certified property holding — non-symmetric
  input to a symmetric-form eigenvalue test, a coefficient-norm tolerance used
  as a value bound, a hash collision from ambiguous serialization.
  (`frankensim-wa8i`, `frankensim-4fbv`.)
- **Domain escape.** Validity asserted outside a model card's or correlation's
  stated domain, or an intersection of validity domains treated as implying
  in-domain validation.
- **Misused pinning authority.** A golden hash, checksum, or content address
  treated as authority for semantics it does not pin — a checksum detecting
  accidental corruption read as authenticity, or a metric-only golden read as
  proof of physical correctness.
- **Fail-open on unresolved state.** An unresolved, budget-exhausted, or
  cancelled computation reported as a clean negative result (a "miss" rather
  than "unresolved"). (`frankensim-certified-query-raycast-qduf`.)

### No-promotion doctrine: verification is not validation

Numerical verification and empirical validation are independent evidence
axes. A `Verified` enclosure may prove that the implementation solved the
declared equations to a stated numerical bound. It does not prove that those
equations, closures, material states, boundary conditions, or observation
operators represent the physical system. Consequently:

- no collection of numerical certificates may mint `Validated`;
- tightening roundoff, algebraic, or discretization terms may not rewrite the
  model-form or measurement terms of an engineering uncertainty budget;
- intersecting declared validity domains proves only a smaller declared
  domain, not agreement with data inside it; and
- a derived `Validated` package claim needs an already-`Validated` parent and
  an independently admitted dataset anchor. A derivation proof cannot replace
  either source of authority.

The canonical reviewer example is RANS. A RANS solve can have a bit-exact
deterministic execution, a tiny nonlinear residual, a converged mesh ladder,
and rigorous roundoff bounds. Those facts can make its *numerical solution of
the selected RANS equations* `Verified`. The predicted pressure drop or heat
transfer remains model-form `Estimated` until a named, independently retained
experimental dataset supports that exact QoI, regime, geometry class, boundary
condition interpretation, and comparison protocol. Smaller residuals do not
tighten the turbulence-closure discrepancy band; new admitted fidelity pairs
or experiments do.

This is the no-promotion row that the E09 certificate-regime table must retain:
proof of the discrete solve and evidence about physical adequacy compose, but
neither axis substitutes for the other.

### Not a claim-integrity defect

Precision matters as much as coverage; an inventory that absorbs every bug
gates nothing. The following are ordinary defects and must NOT carry the label:

- A numerically wrong value whose surface already claims only what it delivers
  (a wrong Estimated value is a correctness bug, not a claim-integrity bug).
- A panic, hang, leak, or unsound `unsafe` block with no epistemic
  overstatement — those are correctness/safety defects with their own labels.
- A claim that is *weaker* than the evidence supports. Under-claiming is a
  usability or opportunity issue, never a claim-integrity defect. This
  asymmetry is intentional and matches the fail-closed doctrine.
- A missing feature, an unimplemented path, or an explicit `NoClaim`/refusal —
  honest refusal is the correct behavior, not a defect.
- Prose that is vague or incomplete but not *stronger* than the code. Only
  overstatement counts.

## Severity rules

Severity is determined by **reachability**, because reachability is what
decides whether a false claim can escape into a decision. It is not a judgment
of how wrong the mathematics is.

| Severity label | Priority | Rule |
| --- | --- | --- |
| `severity:default-path` | P0 | The stronger-than-evidence claim is reachable on a default or public path: default features, a public API used as documented, a shipped report/export, or a gate that runs in `check-all`/DSR. |
| `severity:gated` | P1 | Reachable only behind a non-default feature flag (`frontier-*`, `moonshot-*`), an explicitly opt-in API, or a dev/test-only lane that no shipped path consumes. |
| `severity:doc-only` | P2 | The code is honest; a `README.md`/`CONTRACT.md`/report sentence overstates it. Fixing the sentence fully resolves it. |

Severity escalates, never averages: if any reachable route to the claim is a
default path, the defect is `severity:default-path` even when other routes are
gated. When reachability is genuinely ambiguous, assign the **stronger**
severity and say so in the bead — fail closed, consistent with the gate's
ambiguous-scope rule in `.2.3`.

`severity:default-path` is the gating severity: an open P0 claim-integrity bead
blocks capability-maturity promotion for the capabilities it scopes.

## Label taxonomy

Every claim-integrity defect bead is filed as **`--type=bug`** and carries:

1. `claim-integrity` — mandatory class membership. `br list -l claim-integrity`
   is the live inventory, and the gate reads exactly this label.
2. Exactly one severity label from the canonical set: `severity:default-path`,
   `severity:gated`, `severity:doc-only`. Zero or two or more is a taxonomy
   error that `scripts/ci/claim_integrity_inventory.sh` reports and fails on.
3. `crate:<name>` scope labels for every crate whose surface can emit the
   claim. Scope drives the gate's overlap test; a defect with no crate scope is
   treated as **global** (fail closed, blocking every promotion) rather than as
   unscoped-and-therefore-harmless.

Priority should equal the severity label's priority for open beads. The
inventory script reports mismatches. Closed historical beads are retro-tagged
with the severity that describes the defect as it existed and keep their
original priority; history is a record, not a thing to rewrite.

**Defects versus program beads.** The `claim-integrity` label also marks the
E02 program itself — this epic, its sweep, its gate, and the doctrine features
that enforce the class. Those are work, not exposure. The bead **type** is the
discriminator: a `bug` labeled `claim-integrity` is an inventory entry and is
subject to the severity and ownership rules; any other type (`epic`, `task`,
`feature`) is a program bead, listed for context and exempt. Without this
split, the gate would count its own epic as an open P0 and block every
promotion forever — a gate that can never go green teaches everyone to bypass
it, which is a worse outcome than no gate. File defects as
`br create --type=bug -l claim-integrity -l severity:<...> -l crate:<...>`.

The bead body must contain a minimal repro (the reaching input/state) and the
honest claim the surface should make instead. The sweep bead `.2.2` files
findings in exactly this shape.

## Audit method

This is the procedure the `.2.2` sweep follows, per surface:

1. Write the strongest permitted inference (decision rule 1) for each public
   claim the surface emits. Do this from the type and the contract text before
   reading the implementation, so the implementation cannot anchor the reading.
2. Construct adversarial instances aimed at the gap shapes listed above.
   **Reading the tests is not enough.** The grace-masks-the-mechanism trap —
   tests green while the claim fails — is the normal state of this defect
   class: the fixtures were chosen by the same author who chose the claim, so
   they exercise the region where claim and evidence coincide. Prefer inputs at
   truncation endpoints, empty/degenerate evidence, non-finite values,
   asymmetric or malformed structures, and states where a loop exits on budget
   rather than on convergence.
3. Record a verdict for the surface, clean or not, with the audit note. A
   surface with no recorded verdict is unaudited, which is not the same as
   clean, and the sweep report must distinguish them.
4. File one bead per finding with the labels and body fields above. Do not fix
   in the sweep (beyond trivial doc lines): fixes are separate beads so the
   inventory count stays an accurate measure of open exposure.

## Known instances

The retro-tagged inventory, applied under bead `.2.1`. These predate the label
and are tagged so history is queryable and the sweep can measure its own recall
against a known-answer set.

| Bead | Severity | Shape |
| --- | --- | --- |
| `frankensim-0547u` | `severity:default-path` | Existence read as exactness: WASM `component_count = 1` from a bounded-frame existence argument. Open. |
| `frankensim-flutter-boundary-witness-hukmw` | `severity:default-path` | Truncation read as agreement: `boundaries_agree` with no witnessed transition. Open. |
| `frankensim-2gs5h` | `severity:default-path` | Silent-pass gate: empty golden `schema_fingerprint` accepted as matching. Open. |
| `frankensim-wa8i` | `severity:default-path` | Forgeable certificates across the certificate spine (fs-sos PSD/Lyapunov on non-symmetric input; coefficient-norm read as a value bound; NaN-dropping composition). Closed. |
| `frankensim-y6yv` | `severity:gated` | Real-arithmetic theory read as an executable bound in the feature-gated fs-surrogate ladder; also Estimated coverage laundered into a Verified claim. Closed. |
| `frankensim-f5pr` | `severity:default-path` | Laundering: disjoint-regime composition recorded dispersion `0.0`; NaN comparisons left a Validated claim Validated. Closed. |
| `frankensim-4fbv` | `severity:default-path` | Forged skip certificates from malformed bounds, ambiguous record hashing, and stale plans. Closed. |
| `frankensim-certified-query-raycast-qduf` | `severity:default-path` | Fail-open on unresolved state: unresolved raycast returned as a clean miss; `\|f\|/L` treated as a Euclidean proximity bound. Closed. |

This table is the known-answer set, not the complete inventory. The complete
live inventory is whatever `br list -l claim-integrity` returns; run
`scripts/ci/claim_integrity_inventory.sh` for the checked report.

The first full pass over the claim-emitting surfaces is recorded in
[`docs/CLAIM_INTEGRITY_SWEEP_2026-07.md`](CLAIM_INTEGRITY_SWEEP_2026-07.md),
including its per-surface verdicts, the surfaces it left **unaudited**, and an
honest recall check against the known-answer set above.

## Enforcement

- `cargo run -p xtask -- check-claims` (part of `check-all`) lints that this
  file exists with its required sections and that `docs/CONVENTIONS.md`
  documents the taxonomy with the canonical severity labels. The lint keeps the
  definition present and structurally intact; it does not judge audit quality.
- `scripts/ci/claim_integrity_inventory.sh` prints the live inventory with
  severities and fails when an open `severity:default-path`/P0 defect has no
  owner, when a bead carries zero or multiple severity labels, or when the
  beads store cannot be read (fail closed — an inventory that cannot be read is
  not an empty inventory).
- Bead `.2.3` adds the promotion gate that blocks capability-maturity promotion
  while an in-scope P0 is open; `.2.4` proves that gate fails closed under
  seeded faults.
