# CONTRACT: fs-evidence

## Purpose and layer
`Evidence<T>` / `Certified<T>` (patch Rev B): numerical, statistical, and
MODEL-FORM certificates carried inside values, with a conservative
composition algebra, model cards + the registration lint, two-fidelity
discrepancy models with out-of-distribution refusal, model bracketing, and
decision-aware escalation. The reason this crate exists: without model
evidence the system can produce beautifully certified WRONG answers (mesh
error 0.7%, closure discrepancy 10%) — being able to SAY that is the
product. Layer: UTIL (usable by every layer; the bead label said L6, but
the bead scope explicitly demands a low-layer home — this is it). Depends
on fs-obs only.

## Public types and semantics
- `Evidence<T> { value, qoi, numerical, statistical, model, sensitivity,
  provenance, adjoint_ref }` — the traveling noun. `Certified<T>` is the
  opaque newtype whose constructor discipline is `Evidence::certified()`: rigorous
  numerics (Exact/Enclosure), valid statistical parameters, and in-domain
  non-negative model discrepancy; pure math
  certifies with `ModelEvidence::none()` (the explicit "no model involved"
  statement); refusals are structured `CertifyError`s.
- `NumericalCertificate { kind: Exact|Enclosure|Estimate|NoClaim, lo, hi }`
  — the plan Appendix B `Certified` fields (value + interval bound +
  provenance + adjoint hook) kept intact as the numerical slice. Severity
  is ordered; float composition never claims Exact; NoClaim absorbs.
- `StatisticalCertificate { None | EValue{e, alpha} | HalfWidth{...} }` —
  finite non-negative e-values and widths with levels/confidences strictly in
  `(0,1)`; statistical v1 composition is conservative-weakest (see no-claims).
- `ModelEvidence { cards, assumptions, validity, discrepancy_rel,
  in_domain }`; `ValidityDomain` — named-parameter boxes with
  intersection/containment; `SensitivitySummary` — d(qoi)/d(param)
  headlines, merged by magnitude.
- `Evidence::combine(op, a, b, value)` — Add/Sub/Mul/Min/Max on the QoI
  with certificates composed conservatively and provenance chained
  (`ProvenanceHash::chain`, order-sensitive, FNV until the ledger hash).
- `Evidence::assess(threshold_rel) -> DecisionStatus` and
  `UncertaintyBreakdown` — per-source relative bands, first-order-sum
  total, dominant source with the declaration-order tie law (ModelForm
  first: ties escalate the model, the band cheap refinement cannot fix).
  `escalation_advice` maps dominance to RefineNumerics /
  GatherMoreSamples / EscalateModelFidelity (the HELM governor hook).
- `ModelCard` (name, version, ambition tag, assumptions, validity, known
  failures, calibration provenance, discrepancy band) and `ModelRegistry`
  — `register_solver` REFUSES without a registered card (the lint).
- `DiscrepancyModel::fit(&[FidelityPair])` — observed parameter box +
  mean/max relative discrepancy; `query`/`evidence_at` refuse
  out-of-distribution points and any query whose key set differs from the exact
  training schema, with the first missing, unexpected, non-finite, or
  out-of-range parameter named (`OutOfDomain`).
- `ModelBracket` — N plausible models; evidence = midrange value, an
  enclosure spanning every member, spread as the model band, and a
  bracket-spread sensitivity entry (the vessel flagship's contact-line
  mitigation).
- `to_ledger_row_json` on evidence and cards — the `evidence` /
  `model_cards` table rows (canonical order, no clocks, no addresses).

- `color` module (bead qmao.1): the THREE-COLOR epistemic schema —
  `Color::{Verified{lo,hi}, Validated{regime: ValidityDomain, dataset},
  Estimated{estimator, dispersion}}` with the `ColorRank` lattice
  (verified > validated > estimated), the TOTAL conservative pairwise
  `compose` (the result never outranks the weaker operand; verified intervals
  combine per `IntervalOp` with outward-rounded arithmetic; validated regimes
  INTERSECT, and a disjoint intersection demotes further with
  infinite/no-claim dispersion; estimated absorbs everything with additive
  dispersion),
  `check_regime` (validated is a REGIONAL property: exiting, failing to
  report a regime axis, supplying a non-finite state, or declaring an
  empty/non-finite/inverted regime AUTO-DEMOTES to estimated with a
  `Demotion` flag), `regime_demotion` (the borrowed form used for bounded
  multi-parent admission preflight), `verified_from`
  (the only door to a verified color — non-enclosure certificates
  refuse with the laundering teaching error), and `color_of` (the
  honest bridge from existing Evidence receipts: model-free enclosures may be
  `Verified`, while plain `ModelEvidence` is always `Estimated`). Model cards,
  simulated fidelity pairs, discrepancy bands, and declared validity boxes do
  not authenticate experimental membership and cannot mint `Validated`;
  promotion awaits a typed, independently checkable anchored-source receipt.
  `Color::payload_json`
  escapes caller-controlled strings and represents non-finite floats as
  tagged JSON strings, never invalid bare numeric tokens. The distinct
  `Color::canonical_bytes` identity encoding is versioned (v2), structurally
  length-prefixed, deterministically ordered, and preserves every IEEE-754 bit;
  display rounding therefore never aliases color identity or authorization.
  Color algebra v2 also reserves every `derived:` leaf identity and emits
  bounded composition identities below `derived:v2:` so a caller cannot re-root
  computed evidence as an independently anchored source. Readable derived
  identities are emitted only when every component already satisfies the shared
  grammar; invalid/sentinel components use the domain-separated compact form,
  so even empty-regime demotion returns a structurally valid Estimated payload.
  `color_of` recognizes model absence only for the exact
  `ModelEvidence::none()` shape; an empty card list cannot hide discrepancy,
  assumptions, validity restrictions, or an out-of-domain model behind a
  `Verified` numerical interval.
  Write-time enforcement lives HELM-side in fs-ledger over these types.

- `falsify` module (bead qmao.4): bounded declaration catalog and diagnostic
  telemetry — `FalsifierRegistry` refuses a class with zero declarations and
  `standard()` provides seven intended checker families, separating retained
  sampled-interface replay from the much stronger continuum-watertightness
  proposal (certified oriented intersections, winding/degree, and
  coverage-complete subdivision). `catalog_gate` reports missing declarations
  only; it is explicitly not executable or release authority.
  `FalsifierHistory` ingests bounded, source-referencing `FalsifierAttempt`
  values carrying an idempotency ID, class/regime/falsifier, caller-asserted
  claim-revision and retained-artifact references, seed, positive compute
  charge, and typed outcome. It
  rejects undeclared class/falsifier pairs, treats byte-identical retries as
  no-op replays, and rejects conflicting attempt-ID reuse. `doubt` is the
  empirical discrepancy rate widened by a per-row time-uniform union-bound
  Hoeffding heuristic, with cold-start maximum and a never-zero floor.
  Discrepancy attempts emit opaque escaped, fixed-order
  `fs-evidence/falsifier-candidate` schema-version-1
  tombstone/estimator-bug *candidate* projections, correlated by attempt and
  exact seed/compute bit strings; ingestion neither authenticates the caller's
  references nor adjudicates or persists the candidates.
  `allocate_budget` validates
  and max-rescales consequence × doubt × class-review-share weights before
  normalization. `rent_review` is a preliminary class-level diagnostic over
  fixed ordered `RENT_VOLUME`-attempt windows whose zero-discrepancy windows
  decay toward a nonzero floor. Window closure is independent of review-call
  cadence; subthreshold reviews do not erase observations.

## Invariants
1. Conservativeness (G0, evd-001): composed enclosures contain every
   propagation of operand-enclosed true values (300k seeded samples);
   composed validity is exactly the per-parameter intersection;
   assumptions union sorted; discrepancy bands add; `in_domain` is a
   conjunction. Indeterminate IEEE endpoint arithmetic widens to the whole
   real line instead of discarding NaN corners.
2. Severity monotonicity: composition kind = max operand severity, floored
   at Enclosure for float ops; NoClaim absorbs to infinite bounds; an
   Estimate anywhere poisons `certified()` downstream (evd-006).
3. No card, no solver: `ModelRegistry::register_solver` refuses unknown
   cards with teaching text (evd-002).
4. Out-of-distribution discrepancy queries refuse with the violated
   parameter named — never silent extrapolation. Query keys must equal the
   training schema exactly, so an untrained physical dimension cannot be
   supplied and silently ignored; non-finite training or query coordinates are
   unusable, not a way to synthesize or enter a trained box (evd-004).
5. Dominance ties break in declaration order (ModelForm, Statistical,
   Numerical) — deterministic verdicts.
6. Ledger rows and provenance chains are deterministic (repeat-identical).
   Evidence and model-card rows stay valid JSON under hostile metadata and
   tagged non-finite no-claim values; evidence rows retain model assumptions
   and sensitivity headlines rather than dropping those semantic slices.
   Public set-like card, assumption, and known-failure vectors are sorted and
   deduplicated again at the durable rendering boundary, so caller mutation or
   insertion order cannot change row identity.
   Provenance chaining is order-sensitive.
7. Color regime checks fail closed: no empty, inverted, or non-finite regime
   and no non-finite current state can retain `Validated`; disjoint validated
   composition and regime exit both carry infinite/no-claim dispersion.
8. Certified trust boundary (gp3.2.1, evd-012): `Certified<T>` is an
   OPAQUE newtype over `Evidence<T>` — private inner, no `DerefMut`, no
   field access for writing. The ONLY constructor is
   `Evidence::certified()`, which validates the ACTUAL numbers, not the
  constructor that claimed them: scalar evidence requires bit-identical carried
  value and QoI; Exact requires a finite QoI with
   bit-identical bounds; Enclosure requires finite ordered bounds that
   CONTAIN the QoI; statistical e-values, levels, widths, and confidence
   parameters must satisfy their finite domains; model discrepancy must be
   non-negative (positive infinity is the explicit unbounded claim); empty,
   inverted, or non-finite model-validity domains refuse even when a public
   literal asserts `in_domain: true`; Estimate/NoClaim and out-of-domain models
   refuse. Decision breakdown applies the same validity check and assigns an
   infinite model band to an impossible domain.
   Reads flow through `Deref<Target = Evidence<T>>`;
   `Certified::into_evidence()` is the explicit downgrade — the mark is
   lost and any reconstruction must re-enter `certified()`
   (re-validated round trip). Escape hatches are plain `Evidence<T>` or
   `NumericalCertificate::no_claim()`, never a `Certified<T>`. Certification
   requires an owned/`'static` payload so the boundary can detect scalar `f64`
   values and bind them bit-exactly to their QoI; borrowed payloads remain
   plain evidence until promoted to owned values.
   MIGRATION: `Certified<T>` was a type alias for `Evidence<T>`;
   callers that mutated or moved fields of a certified value now call
   `.into_evidence()` first (one workspace site: fs-geom conformance),
   and callers that only read fields are unchanged via Deref. This
   crate has no serializer; persisted evidence re-enters through
   `certified()` on ingest by construction.
9. Decision assessment fails closed (evd-013): malformed or negative
   uncertainty becomes an infinite band; infinite totals and malformed
   thresholds cannot become `DecisionGrade`.
10. Color provenance identities are bounded, byte-length-framed, and
    domain-separated (evd-014). Every generated identity is in the reserved
    `derived:v2:` namespace, while all `derived:` identities are forbidden at
    source admission so computed evidence cannot be re-rooted as a leaf.
11. The shared public color validator and `color_of` bridge fail closed on
    NaN/inverted intervals, malformed identities and regimes, negative/NaN
    dispersion, and malformed evidence/model inputs (evd-015). Plain model
    evidence never self-promotes to `Validated`, including in-domain
    two-fidelity discrepancy evidence (evd-004). Ordered infinite Verified
    endpoints remain sound but vacuous enclosures.
12. `demotion_estimator_identity` is total over arbitrary strings and always
    emits a bounded identity accepted by `color_identity_reason`; invalid
    readable inputs are hash-compacted rather than interpolated. Malformed
    model-card diagnostics bind the complete sorted/deduplicated card set and
    each entry's validation reason in a v2 domain-separated streaming hash;
    distinct residual cards cannot alias merely because the first invalid card
    is the same.

## Error model
Structured teaching errors throughout: `CertifyError`, `RegistryError`,
`OutOfDomain`, `FitError`, `FalsifyError` — all `core::error::Error` with actionable
Display text. Constructors are total (enclosure bounds normalize by
swapping); no panics cross the boundary.

## Determinism class
Deterministic: pure values and mutable diagnostic state machines produce the
same results for the same ordered call sequence; renderings use sorted
(`BTreeMap`) order; there are no clocks, addresses, or hidden randomness.
Bit-stable across runs and platforms up to fs-math-class scalar-arithmetic
divergence.

## Cancellation behavior
Core certificate/color algebra is bounded small synchronous work. Falsifier
allocation and history review iterate caller data without a `Cx`; allocation
length and distinct history rows are defensively capped, but these diagnostic
APIs are not P7 hot-kernel or cancellation-authoritative paths. Callers must not
place large reviews inside latency-bounded tile loops.

## Unsafe boundary
None. `unsafe_code` denied workspace-wide.

## Feature flags
None. The mechanisms are `[S]`-grade bookkeeping; the models they DESCRIBE
carry their own ambition tags in their cards.

## Conformance tests
tests/conformance.rs, cases evd-001..evd-015 (JSON-line verdicts; seeded
cases carry seeds): the G0 conservativeness battery, the registration
lint, the worked model-discrepancy-dominates example (10% closure vs 0.7%
mesh at a 5% threshold → NotDecisionGrade{ModelForm} + escalation advice,
flipping to DecisionGrade with a 2% calibrated closure), the
out-of-distribution refusal on a synthetic two-fidelity corpus, bracketing
spread reporting with deterministic schema-valid rows, and certification
poisoning. The color cases cover disjoint-regime demotion with infinite
dispersion, outward-rounded verified arithmetic, non-finite state/regime
demotion, and deterministic escaping/tagging of hostile JSON payloads.
In-module suites cover the certificate algebra, validity laws, tie-breaking,
provenance chaining, and card rendering.
`evd-013` exercises the public `Evidence` layer against indeterminate interval
arithmetic and malformed statistical/model uncertainty.
`evd-014` locks bounded collision-resistant derived identities and source
namespace separation; `evd-015` locks the shared color structural invariant
across every public variant and bridge input, including the malformed empty-
regime demotion sentinel and complete-set malformed-card identity stability.
`evd-016` uses `fs-propcheck` seed `0xE71D_4A48_0001` for 512 shrinkable color
composition recipes: all 27 operand-kind/operation combinations, both operand
orders, overlapping and disjoint validated regimes, and structural validation
of every input and output. It locks the no-laundering rank law while preserving
the fixed cases above as regression pins.

- Falsifier registration is total: malformed, duplicate, and empty declaration
  sets refuse; the catalog lint names each distinct undeclared class once in
  canonical sorted order without making a release claim.
- Fixed and adversarial falsifier cases check that budget allocation is
  monotone in consequence AND doubt, with honest
  boundaries (cold-start max, asymptotic doubt floor, dependent-free floor,
  empty-job zero), finite max-rescaled extreme weights, and validation even
  when the requested total budget is zero.
- Every accepted discrepancy attempt produces BOTH correlated pending candidate
  payloads; neither candidate is an active tombstone or an adjudicated
  estimator bug. Exact retries are idempotent, conflicting ID reuse refuses,
  malformed/undeclared attempts leave telemetry unchanged, and two clean rent
  windows decay identically under batched versus incremental review calls.

## Admitted scientific color (bead 6pf9, stage S1)

- `Color` is a DECLARATION: publicly constructible, structurally validated,
  never authority. `AdmittedColor` is the opaque positive-evidence handle:
  private fields, single constructor `AdmittedColor::from_receipt(color,
  receipt, verifier)`.
- Local gates fire before the injected capability: malformed payloads
  (`validate_color_payload`), non-positive ranks (Estimated), and
  stale-algebra receipts (`color_algebra_version !=
  COLOR_ALGEBRA_VERSION`) refuse even under an accept-everything verifier.
- An accepting `AdmissionDecision` must name exactly the policy fingerprint
  committed by the receipt. A mismatch returns structured `PolicyMismatch`
  rather than silently detaching the admitted value from its lineage; retained
  tests cover both a matching authority and an accepting policy substitution.
- `AdmissionReceipt` is plain data (node provenance hash, row schema
  version, algebra version, policy fingerprint). Authority lives in the
  `AdmissionVerifier` capability; the default `NoAdmissionVerifier` is
  deny-all, so at this layer NOTHING admits. The authenticating verifier is
  HELM-side (`fs_ledger::LedgerColorAdmissionVerifier`), keeping this crate
  free of upward dependencies.
- No-claim: `AdmittedColor` is capability-gated, not cryptographically
  unforgeable. A lying verifier at the composition root can admit anything —
  the same trust model as `WaiverVerifier` and `SourceOriginVerifier`, and
  visible at the same audit surface. Stage S1 is additive: positive-evidence
  consumer APIs migrate to require `AdmittedColor` in later stages of bead
  6pf9.

## No-claim boundaries (colors)

- Verified-interval composition here covers outward-rounded Add/Mul and exact
  endpoint Hull; the full ledger operation algebra composes through fs-ivl
  when wired.
- Estimated dispersion combines additively (conservative); calibrated
  dispersion algebra joins the color-probes bead.
- `ModelEvidence` carries model-form diagnostics, not experimental authority.
  Even an in-domain discrepancy model trained on paired simulations remains
  `Estimated`; this crate has no authenticated anchor type that can admit a new
  `Validated` leaf.
- Estimator identity chains are human-readable up to
  `MAX_COLOR_IDENTITY_BYTES` and then collapse to a domain-separated BLAKE3
  composition identity. Human-readable pairs are byte-length-framed, so an
  identity containing `+` or `@` cannot alias a different composition tree.
  Every binary composition retains both framed operands, including repeated
  identities, and mixed `Verified`/`Validated`/`Estimated` operand classes use
  distinct domain labels. A source literally named `verified` therefore cannot
  impersonate a Verified operand. Pass-through operands receive a bounded
  `derived:v2:` identity rather than retaining a source-leaf identity. Single
  model-card identities remain source identities until an operation composes
  them.
  This bounds provenance payload growth without changing the Estimated rank or
  its accumulated dispersion. Compact identities are deterministic provenance
  labels, not signatures or source authority.

## No-claim boundaries
- Statistical composition is CONSERVATIVE-WEAKEST v1 (half-widths add,
  confidences min, mixed kinds keep the width-bearing certificate);
  proper e-value arithmetic (products under independence, e-BH) is
  fs-eproc's contract and will replace this composition rule behind the
  same API.
- Discrepancy models are honest BOOKKEEPING (observed box + mean/max
  band), not learning; trained/learned discrepancy models (FrankenTorch)
  arrive with fs-surrogate and will implement the same query/refusal
  surface.
- The adjoint hook is carried, never composed here — composed tapes are
  fs-ad's contract.
- First-order band addition is conservative for small relative bands; it
  is NOT a rigorous product-form bound for large ones — fs-ivl composition
  should be used for the numerical slice when bands are large (the
  numerical slice already does interval arithmetic; the TOTAL across
  sources is the first-order sum).
- Ledger persistence is row RENDERING only; the `evidence` /
  `model_cards` tables land with fs-ledger (rows are shaped for that
  migration).
- `ProvenanceHash` is FNV-1a until the BLAKE3-class ledger hash supersedes
  it (same upgrade path as fs-obs).

## No-claim boundaries (falsifiers)

- The registry stores falsifier IDENTITIES and stated methods; executing
  a falsifier (running rays, FD probes, full solves) is each consumer
  kernel's code. A public method string cannot prove executable binding or
  checker independence, and `catalog_gate` cannot authorize release.
- Catalog class/spec identifiers are bounded ASCII slugs; human prose belongs
  in the separately bounded method/detail fields. Class count, declarations per
  class, per-class and total declaration bytes, lint request/output size,
  allocation length, history rows/key bytes, attempt count, and retained
  attempt bytes all have explicit ceilings. These are defensive synchronous
  caps, not a substitute for a budgeted `Cx` API.
- Exact-instance authority requires a typed retained receipt bound to the claim,
  specification, policy, implementation/TCB, seed, stopping rule, budget,
  outcome, and artifacts, admitted by fs-package/fs-checker. That successor
  architecture is not implemented by this telemetry module.
- `consequence` is supplied by the caller (ledger-DAG dependent-weight
  traversal is HELM-side); the allocator's contract is what it does with
  the number, floors included.
- Candidate payloads are opaque escaped
  `fs-evidence/falsifier-candidate` schema-version-1 fixed-order JSON
  projections, with the full-width seed and compute charge represented as exact
  hexadecimal bit strings; they are not authenticated content addresses.
  History mutation is not transactional with an external
  sink; downstream adjudication and atomic persistence are mandatory before a
  candidate can invalidate a claim or count as an estimator defect.
- Stable attempt IDs prevent identical retry/replay from double-counting and
  conflicting reuse fails closed. They remain caller assertions rather than
  authenticated ledger identities, so history still cannot be trusted against
  fabricated or selectively omitted attempts. `doubt` assumes
  an admitted stationary Bernoulli process and controls time multiplicity only
  within one row, not the family of class/regime rows; here it is planning
  telemetry only. Candidate discrepancy reports conservatively raise this
  telemetry but cannot activate authority.
- Rent decay is per-class rather than per-falsifier. Fixed-volume windows are
  invariant to review-call cadence, but are not authenticated calendar/policy
  windows and have no restoration policy. Per-falsifier
  receipts, independence threat graphs, cadence enforcement, and restoration
  policy remain successor work (xpck.6).
